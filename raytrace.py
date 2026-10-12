"""Schwarzschild ray tracing for an opaque, equatorial, Keplerian thin disk.

Distances are in r_g = GM/c^2. The observer is static at infinity. We integrate
the exact null-orbit equation u'' = 3 u^2 - u, u = 1/r, with fourth-order
Runge--Kutta; the numerical lookup is reused for every observer inclination.
The first emitting disk intersection occludes all subsequent intersections.
This module contains no invented emission, image effects or fluid evolution.

The emitter frequency shift is g = sqrt(1-3/r)/(1-Omega*lambda), with
Omega = r^(-3/2), lambda = -alpha sin(i). The disk must lie outside r_ISCO=6.
Specific intensity transforms as I_nu,obs = g^3 I_nu,em(nu_obs/g); bolometric
intensity transforms as g^4. Emissivity is the renderer's responsibility.

References: Luminet (1979), A&A 75, 228; James et al. (2015),
arXiv:1502.03808 (ray-traced appearance and relativistic transfer). This is
an original Schwarzschild implementation, not DNGR, AART or a GRMHD solver.
"""
from __future__ import annotations

import argparse
import gc
import json
import math
import time

import numpy as np

CRITICAL_IMPACT = math.sqrt(27.0)


def _orbit_table(impacts, max_angle, angle_step):
    """Return u(psi,b), total escape angle and actual uniform angular step."""
    steps = int(math.ceil(max_angle / angle_step))
    h = max_angle / steps
    u = np.zeros(impacts.size, dtype=np.float64)
    v = 1.0 / np.asarray(impacts, dtype=np.float64)
    active = np.ones(impacts.size, dtype=bool)
    table = np.empty((steps + 1, impacts.size), dtype=np.float32)
    table[0] = 0.0
    escape = np.full(impacts.size, np.nan)
    for step in range(1, steps + 1):
        k1u, k1v = v, 3*u*u-u
        u2, v2 = u+h*.5*k1u, v+h*.5*k1v
        k2u, k2v = v2, 3*u2*u2-u2
        u3, v3 = u+h*.5*k2u, v+h*.5*k2v
        k3u, k3v = v3, 3*u3*u3-u3
        u4, v4 = u+h*k3u, v+h*k3v
        k4u, k4v = v4, 3*u4*u4-u4
        un = u+h/6*(k1u+2*k2u+2*k3u+k4u)
        vn = v+h/6*(k1v+2*k2v+2*k3v+k4v)
        escaped = active & (un <= 0) & (step > 1)
        captured = active & (un >= .5)
        if escaped.any():
            escape[escaped] = (step-1+u[escaped]/(u[escaped]-un[escaped]))*h
        un = np.where(active, un, u)
        un[captured], un[escaped] = .5, 0.
        active &= ~(escaped | captured)
        u, v = un, np.where(active, vn, 0.)
        table[step] = u
    return table, escape, h


class SchwarzschildTracer:
    """Cache a ray-orbit table and isotropic rectangular detector grid.

    fov is the FULL horizontal screen span in r_g. A rectangular detector has
    the same physical pixel size in both directions. Positive roll rotates the
    disk's projected horizontal axis counterclockwise in the displayed image.
    trace() returns arrays with shape (height,width), first row at image top.
    """
    def __init__(self, width, height, fov=32, roll_deg=25, rin=6, rout=40,
                 layers=4, radial_samples=8192, angle_step=.002):
        if not isinstance(width, int) or not isinstance(height, int) or min(width,height)<8:
            raise ValueError('width and height must be integers >=8')
        if not (fov>0 and 6<=rin<rout and 1<=layers<=6 and radial_samples>=256 and angle_step>0):
            raise ValueError('invalid disk, field of view or integration settings')
        started = time.perf_counter()
        self.width, self.height = width, height
        self.fov, self.roll_deg = float(fov), float(roll_deg)
        self.rin, self.rout, self.layers = float(rin), float(rout), int(layers)
        scale = fov/width
        ax = (np.arange(width,dtype=np.float64)+.5-width/2)*scale
        ay = -(np.arange(height,dtype=np.float64)+.5-height/2)*scale
        xd, yd = np.meshgrid(ax, ay)
        roll = math.radians(roll_deg)
        # Inverse detector rotation expresses every pixel in physical axes.
        x = np.cos(roll)*xd+np.sin(roll)*yd
        y = -np.sin(roll)*xd+np.cos(roll)*yd
        b = np.hypot(x,y)
        safe_b = np.maximum(b,1e-12)
        self.x = x.astype(np.float32)
        self.y = y.astype(np.float32)
        self.ex = (x/safe_b).astype(np.float32)
        self.ey = (y/safe_b).astype(np.float32)
        self.impact = b.astype(np.float32)
        self.captured = b<CRITICAL_IMPACT
        # Logarithmic refinement resolves the sensitive neighbourhood of the
        # photon orbit without requiring an enormous uniform impact grid.
        dense = np.geomspace(1e-7,.3,1024)
        impacts = np.unique(np.concatenate((
            np.linspace(.01,float(b.max())+scale,radial_samples),
            CRITICAL_IMPACT-dense, CRITICAL_IMPACT+dense,
        )))
        self.table, _, self.h = _orbit_table(impacts,layers*math.pi,angle_step)
        ri = np.clip(np.searchsorted(impacts,b)-1,0,impacts.size-2)
        self.ri = ri.astype(np.int32)
        self.rw = np.clip((b-impacts[ri])/(impacts[ri+1]-impacts[ri]),0,1).astype(np.float32)
        self.min_b_mask = b>=impacts[0]
        self.metadata = {
            'model':'Schwarzschild null geodesics; first opaque Keplerian thin-disk hit',
            'width':width,'height':height,'fov_horizontal_rg':fov,
            'fov_vertical_rg':height*scale,'roll_deg':roll_deg,
            'rin_rg':rin,'rout_rg':rout,'layers':layers,
            'radial_samples_total':int(impacts.size),'angle_step':self.h,
            'critical_impact_rg':CRITICAL_IMPACT,
            'orbit_cache_mb':self.table.nbytes/1024**2,
            'setup_seconds':time.perf_counter()-started,
        }

    def trace(self, inclination_deg):
        """Return r, phi, g, mask, hit_layer and shadow for one inclination.

        Invalid pixels have r=phi=0, g=1 and hit_layer=-1. 'shadow' denotes
        captured rays not occulted by foreground emitting material. A black
        image region can also contain escaping rays that encounter no source.
        """
        if not 0<=inclination_deg<89.9:
            raise ValueError('inclination_deg must lie in [0,89.9)')
        i = math.radians(inclination_deg)
        si, ci = math.sin(i), math.cos(i)
        first = np.arctan2(ci,-si*self.ey)
        r = np.zeros((self.height,self.width),dtype=np.float32)
        phi = np.zeros_like(r)
        g = np.ones_like(r)
        hit_layer = np.full(r.shape,-1,dtype=np.int8)
        mask = np.zeros(r.shape,dtype=bool)
        lam = -self.x*si
        for layer in range(self.layers):
            psi = first+layer*math.pi
            ai = np.minimum((psi/self.h).astype(np.int32),self.table.shape[0]-2)
            aw = psi/self.h-ai
            ri, rw = self.ri,self.rw
            u0 = self.table[ai,ri]*(1-rw)+self.table[ai,ri+1]*rw
            u1 = self.table[ai+1,ri]*(1-rw)+self.table[ai+1,ri+1]*rw
            u = u0*(1-aw)+u1*aw
            radius = np.divide(1.,u,out=np.full_like(u,np.inf),where=u>0)
            new = (~mask)&self.min_b_mask&(radius>=self.rin)&(radius<=self.rout)
            if new.any():
                rr = radius[new]
                angle = psi[new]
                qx = self.ex[new]*np.sin(angle)
                qy = -si*np.cos(angle)+ci*self.ey[new]*np.sin(angle)
                r[new] = rr
                phi[new] = np.arctan2(qy,qx)
                g[new] = np.sqrt(1.-3./rr)/(1.-lam[new]/rr**1.5)
                hit_layer[new] = layer
                mask[new] = True
        return {'r':r,'phi':phi,'g':g,'mask':mask,'hit_layer':hit_layer,
                'shadow':self.captured&~mask}


def diagnostics():
    """Check bending, symmetry, handedness and numerical convergence."""
    started = time.perf_counter()
    tracer = SchwarzschildTracer(192,256,roll_deg=0,radial_samples=4096,angle_step=.004)
    low = tracer.trace(68)
    metadata = tracer.metadata.copy()
    mask = low['mask']
    if not np.array_equal(mask,np.fliplr(mask)):
        raise AssertionError('nonspinning geometry lost left/right symmetry')
    if not np.allclose(low['r'],np.fliplr(low['r']),atol=1e-6):
        raise AssertionError('disk radii lost reflection symmetry')
    if not np.all((low['r'][mask]>=6)&(low['r'][mask]<=40)):
        raise AssertionError('emission outside prescribed stable disk')
    approaching = mask&(tracer.x<0)
    receding = mask&(tracer.x>0)
    ga, gr = float(low['g'][approaching].mean()),float(low['g'][receding].mean())
    # An arbitrary positive symmetric radial emissivity is enough to test the
    # relativistic-transfer handedness without testing an emission model here.
    flux = np.zeros(mask.shape,dtype=np.float64)
    flux[mask] = low['g'][mask]**4/low['r'][mask]**3
    flux_ratio = float(flux[approaching].sum()/flux[receding].sum())
    if not ga>gr or not flux_ratio>1:
        raise AssertionError('Doppler/bolometric-transfer handedness is wrong')
    _,escape,_ = _orbit_table(np.array([200.]),4.,.001)
    deflection = float(escape[0]-math.pi)
    if abs(deflection-.02)/.02>.03:
        raise AssertionError('weak-field deflection disagrees with 4M/b')
    del tracer
    gc.collect()
    refined = SchwarzschildTracer(192,256,roll_deg=0,radial_samples=8192,angle_step=.002)
    high = refined.trace(68)
    common = mask&high['mask']&(low['hit_layer']==high['hit_layer'])
    rel_r = np.abs(low['r'][common]-high['r'][common])/high['r'][common]
    mask_disagreement = float(np.mean(mask!=high['mask']))
    p99 = float(np.percentile(rel_r,99))
    high_flux = np.zeros(mask.shape,dtype=np.float64)
    high_flux[high['mask']] = high['g'][high['mask']]**4/high['r'][high['mask']]**3
    flux_rel = float(abs(flux.sum()-high_flux.sum())/high_flux.sum())
    if mask_disagreement>.002 or p99>.005 or flux_rel>.005:
        raise AssertionError('resolution refinement did not converge')
    return {'passed':True,'geometry_reflection_symmetric':True,
            'approaching_mean_g':ga,'receding_mean_g':gr,
            'approaching_to_receding_test_flux':flux_ratio,
            'weak_field_deflection_b200':deflection,'weak_field_leading_order':.02,
            'refinement_mask_disagreement_fraction':mask_disagreement,
            'refinement_radius_p99_relative_error':p99,
            'refinement_integrated_test_flux_relative_error':flux_rel,
            'visible_first_hit_layers':{str(k):int((low['hit_layer']==k).sum()) for k in range(4)},
            'coarse_cache':metadata,'elapsed_seconds':time.perf_counter()-started}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',help='optional JSON diagnostics output path')
    args = parser.parse_args()
    result = diagnostics()
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
