"""Schwarzschild thin-disk geodesics with relative light-travel times.

This is an original implementation, not NASA source code. The geometry,
opaque-disk selection and redshift extend the existing SchwarzschildTracer
to a continuous, full-circle observer orientation. The extra quantity is a
relative propagation time, in t_g=GM/c^3, with one
common additive constant removed at the observer at infinity.

Derivation in G=M=c=1 units, with psi increasing backwards from the observer:
    u=1/r, v=du/dpsi, v^2=1/b^2-u^2+2u^3, v'=3u^2-u;
    d(ell)/dpsi=1/[b*u^2*(1-2u)], r*=r+2 log(r/2-1).
Here ell is positive elapsed backward coordinate time. The finite quantity
q=ell+r*(r)-r*(observer) starts at zero and obeys
    q'=(1/b-v)/[u^2*(1-2u)]=1/(1/b+v).
The reciprocal form avoids cancellation as u tends to zero on the incoming
branch; the first expression avoids cancellation on the outgoing branch.
Thus delay=q-r*(r_hit)=ell-r*(observer). Evaluate variable emission at
t_em=t_obs-delay_tg, with the same arbitrary origin in all pixels. Negative
relative delays are valid; they do not mean negative physical travel times.

The regularization is derived here from the Schwarzschild null equations,
not quoted as a formula implemented by NASA. Background references:
Muller & Frauendiener (2012), https://arxiv.org/abs/1206.4259;
Gralla & Lupsasca (2020), https://arxiv.org/abs/1910.12881 (general Kerr
geodesic equations, of which this is the nonspinning specialization).
Luminet (1979), A&A 75, 228 is the thin-disk imaging reference.
"""
from __future__ import annotations

import argparse
import gc
import json
import math
import time
from pathlib import Path

import numpy as np

from raytrace import CRITICAL_IMPACT, SchwarzschildTracer


def tortoise_radius(radius, mass=1.0):
    """Schwarzschild tortoise coordinate; mass=0 is only for flat-space QA."""
    radius = np.asarray(radius, dtype=np.float64)
    if mass == 0:
        return radius
    if np.any(radius <= 2*mass):
        raise ValueError('tortoise radius is defined only outside the horizon')
    return radius + 2*mass*np.log(radius/(2*mass)-1)


def _regularized_derivative(u, v, inverse_b, mass):
    # Values at numerical RK stages beyond escape are never committed. Keep
    # their evaluation finite so discarded trial steps do not poison a table.
    outgoing = (v < 0) & (u > 0) & (1-2*mass*u > 0)
    result = np.divide(1., inverse_b+v, out=np.zeros_like(v),
                       where=(~outgoing) & (inverse_b+v > 0))
    numerator = inverse_b-v
    denominator = u*u*(1-2*mass*u)
    np.divide(numerator, denominator, out=result, where=outgoing)
    return result


def _orbit_table_with_time(impacts, max_angle, angle_step, mass=1.0):
    """Integrate u and finite q using the same RK4 stages as the baseline.

    The extra mass argument exists solely to exercise the exact flat-space
    limit in diagnostics; production always uses mass=1 in gravitational units.
    Escaped/captured rays freeze. Their endpoint table entries do not encode a
    physical disk hit and must only be used together with the radius mask.
    """
    steps = int(math.ceil(max_angle/angle_step))
    h = max_angle/steps
    inverse_b = 1/np.asarray(impacts, dtype=np.float64)
    u = np.zeros(inverse_b.size, dtype=np.float64)
    v = inverse_b.copy()
    q = np.zeros_like(u)
    active = np.ones(u.size, dtype=bool)
    radial = np.empty((steps+1, u.size), dtype=np.float32)
    regularized = np.empty_like(radial)
    radial[0] = regularized[0] = 0
    capture_u = 1/(2*mass) if mass else math.inf
    for step in range(1, steps+1):
        k1u, k1v = v, 3*mass*u*u-u
        k1q = _regularized_derivative(u, v, inverse_b, mass)
        u2, v2 = u+h*.5*k1u, v+h*.5*k1v
        k2u, k2v = v2, 3*mass*u2*u2-u2
        k2q = _regularized_derivative(u2, v2, inverse_b, mass)
        u3, v3 = u+h*.5*k2u, v+h*.5*k2v
        k3u, k3v = v3, 3*mass*u3*u3-u3
        k3q = _regularized_derivative(u3, v3, inverse_b, mass)
        u4, v4 = u+h*k3u, v+h*k3v
        k4u, k4v = v4, 3*mass*u4*u4-u4
        k4q = _regularized_derivative(u4, v4, inverse_b, mass)
        un = u+h/6*(k1u+2*k2u+2*k3u+k4u)
        vn = v+h/6*(k1v+2*k2v+2*k3v+k4v)
        qn = q+h/6*(k1q+2*k2q+2*k3q+k4q)
        escaped = active & (un <= 0) & (step > 1)
        captured = active & (un >= capture_u)
        # At escape the true q diverges. Freeze at the last finite point;
        # our finite emitting disk is far inside these large-radius cells.
        q = np.where(active & ~escaped, qn, q)
        un = np.where(active, un, u)
        un[captured], un[escaped] = capture_u, 0.
        active &= ~(escaped | captured)
        u, v = un, np.where(active, vn, 0.)
        radial[step] = u
        regularized[step] = q
    return radial, regularized, h


class SchwarzschildDelayTracer(SchwarzschildTracer):
    """Baseline tracer API plus delay_tg, zero for every invalid pixel.

    All delays share one zero point; DO NOT independently normalize each
    image order or frame. A single scene-wide time shift is harmless.
    """
    def __init__(self, width, height, fov=70, roll_deg=0, rin=6, rout=30,
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
        x = np.cos(roll)*xd+np.sin(roll)*yd
        y = -np.sin(roll)*xd+np.cos(roll)*yd
        b = np.hypot(x,y)
        safe_b = np.maximum(b,1e-12)
        self.x, self.y = x.astype(np.float32), y.astype(np.float32)
        self.ex, self.ey = (x/safe_b).astype(np.float32), (y/safe_b).astype(np.float32)
        self.impact = b.astype(np.float32)
        self.captured = b<CRITICAL_IMPACT
        dense = np.geomspace(1e-7,.3,1024)
        impacts = np.unique(np.concatenate((
            np.linspace(.01,float(b.max())+scale,radial_samples),
            CRITICAL_IMPACT-dense, CRITICAL_IMPACT+dense,
        )))
        self.table, self.time_table, self.h = _orbit_table_with_time(
            impacts,layers*math.pi,angle_step)
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
            'time_cache_mb':self.time_table.nbytes/1024**2,
            'delay_units':'t_g=GM/c^3',
            'delay_definition':'coordinate propagation time minus common observer r*; q-r*(hit)',
            'invalid_delay_value':0,
            'emission_time_rule':'t_em=t_obs-delay_tg; one common additive origin',
            'camera_angle_convention':'signed degrees about disk x axis; n=(0,-sin(theta),cos(theta))',
            'camera_screen_basis':'e_x=(1,0,0); e_y=(0,cos(theta),sin(theta)); continuous through both poles',
            'disk_plane_crossing':'positive root psi=pi for noncoplanar rays; coplanar rays first hit outer rim',
            'disk_surfaces':'same prescribed isotropic surface emission on both faces',
            'setup_seconds':time.perf_counter()-started,
        }

    def trace(self, inclination_deg):
        """Trace a signed observer angle, with a continuous camera screen basis.

        The argument is an oriented angle, not an inclination restricted to
        [0, pi/2]. A ray has z/r = cos(theta) cos(psi) +
        sin(theta) e_y sin(psi). Its disk-plane roots are separated by pi;
        selecting the first positive root generalizes both faces without
        reflecting an image or reversing the source's evolution.

        At an exactly equatorial observer, psi=0 is the observer at infinity,
        not a source hit. Noncoplanar rays first cross the disk at psi=pi.
        Exactly coplanar rays instead enter the opaque disk at its outer rim.
        The razor-thin direct image has zero angular height at this instant;
        the strongly lensed far-side image can still have finite area.
        """
        angle_deg = float(inclination_deg)
        if not math.isfinite(angle_deg):
            raise ValueError('observer angle must be finite')
        angle = math.radians(math.remainder(angle_deg,360.))
        si, ci = math.sin(angle), math.cos(angle)
        # Make the exact cardinals reproducible instead of letting the sign
        # of floating-point cos(pi/2) choose different disk-plane branches.
        if abs(si) < 2e-15:
            si = 0.
        if abs(ci) < 2e-15:
            ci = 0.
        first = np.arctan2(ci,-si*self.ey)
        if ci < 0:
            # Negating both coefficients yields the positive root directly,
            # avoiding a float32 subtraction of pi near the opposite pole.
            first = np.arctan2(-ci,si*self.ey)
        elif ci == 0:
            first = np.full_like(self.ey,math.pi)
        coplanar = (ci == 0) & (self.ey == 0)
        eligible = self.min_b_mask & ~coplanar
        r = np.zeros((self.height,self.width),dtype=np.float32)
        phi = np.zeros_like(r)
        g = np.ones_like(r)
        delay = np.zeros_like(r)
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
            new = (~mask)&eligible&(radius>=self.rin)&(radius<=self.rout)
            if not new.any():
                continue
            rr, path_angle = radius[new], psi[new]
            qx = self.ex[new]*np.sin(path_angle)
            qy = -si*np.cos(path_angle)+ci*self.ey[new]*np.sin(path_angle)
            r[new] = rr
            phi[new] = np.arctan2(qy,qx)
            g[new] = np.sqrt(1.-3./rr)/(1.-lam[new]/rr**1.5)
            hit_layer[new] = layer
            mask[new] = True
            aj, wj, rj, wr = ai[new], aw[new], ri[new], rw[new]
            q0 = self.time_table[aj,rj]*(1-wr)+self.time_table[aj,rj+1]*wr
            q1 = self.time_table[aj+1,rj]*(1-wr)+self.time_table[aj+1,rj+1]*wr
            delay[new] = q0*(1-wj)+q1*wj-tortoise_radius(r[new])

        # Production uses an even detector, so this special set is empty.
        # Supporting it explicitly also makes odd detector sizes well defined.
        rows, columns = np.nonzero(coplanar & self.min_b_mask)
        for offset in range(0,len(rows),128):
            ys, xs = rows[offset:offset+128], columns[offset:offset+128]
            ri, rw = self.ri[ys,xs], self.rw[ys,xs]
            radial = self.table[:,ri]*(1-rw)+self.table[:,ri+1]*rw
            entered = radial >= 1/self.rout
            valid = entered.any(axis=0)
            if not valid.any():
                continue
            k = entered.argmax(axis=0)[valid]
            ys, xs, ri, rw = ys[valid], xs[valid], ri[valid], rw[valid]
            j = np.flatnonzero(valid)
            fraction = (1/self.rout-radial[k-1,j])/(radial[k,j]-radial[k-1,j])
            path_angle = (k-1+fraction)*self.h
            q0 = self.time_table[k-1,ri]*(1-rw)+self.time_table[k-1,ri+1]*rw
            q1 = self.time_table[k,ri]*(1-rw)+self.time_table[k,ri+1]*rw
            r[ys,xs] = self.rout
            phi[ys,xs] = np.arctan2(-si*np.cos(path_angle),
                                      self.ex[ys,xs]*np.sin(path_angle))
            g[ys,xs] = math.sqrt(1-3/self.rout)/(1-lam[ys,xs]/self.rout**1.5)
            delay[ys,xs] = q0*(1-fraction)+q1*fraction-tortoise_radius(self.rout)
            hit_layer[ys,xs] = 0
            mask[ys,xs] = True
        radial_center = coplanar & (self.impact == 0)
        if radial_center.any():
            # b=0 has no orbital angle, but the incoming radial null ray
            # reaches the near outer rim with ell-r*(observer)=-r*(rout).
            r[radial_center] = self.rout
            phi[radial_center] = math.atan2(-si,0.)
            g[radial_center] = math.sqrt(1-3/self.rout)
            delay[radial_center] = -tortoise_radius(self.rout)
            hit_layer[radial_center] = 0
            mask[radial_center] = True
        if not np.all(np.isfinite(delay)):
            raise RuntimeError('nonfinite relative travel time in disk-hit map')
        return {'r':r,'phi':phi,'g':g,'mask':mask,'hit_layer':hit_layer,
                'shadow':self.captured&~mask,'delay_tg':delay}


def _production_view_diagnostics(full_resolution=False):
    """Refine lookup tables at the delivered scene's geometry.

    The inexpensive views sample both geodesic tables on the same detector.
    The optional 2x full-resolution check separately exercises the actual
    production array sizes; it does not claim unresolved image orders are
    individually resolved by an Instagram-sized detector.
    """
    options = dict(width=324, height=300, fov=60, roll_deg=0, rin=6,
                   rout=28, layers=4, radial_samples=8192, angle_step=.002)
    inclinations = (82., 83.5, 85.)
    tracer = SchwarzschildDelayTracer(**options)
    coarse = {inclination: tracer.trace(inclination) for inclination in inclinations}
    del tracer
    gc.collect()
    tracer = SchwarzschildDelayTracer(**{**options, 'radial_samples':16384,
                                        'angle_step':.001})
    refinement = []
    for inclination in inclinations:
        low, high = coarse[inclination], tracer.trace(inclination)
        common = low['mask'] & high['mask'] & (low['hit_layer'] == high['hit_layer'])
        difference = np.abs(low['delay_tg'][common]-high['delay_tg'][common])
        phase = difference/high['r'][common]**1.5
        radius_error = np.abs(low['r'][common]/high['r'][common]-1)
        disagreement = float(np.mean(low['mask'] != high['mask']))
        p99 = float(np.percentile(difference, 99))
        max_phase = float(phase.max())
        if disagreement > .001 or p99 > .01 or max_phase > .005:
            raise AssertionError('delay/geometry refinement failed at a production inclination')
        refinement.append({
            'inclination_deg':inclination,
            'mask_disagreement_fraction':disagreement,
            'radius_p99_relative_error':float(np.percentile(radius_error, 99)),
            'delay_p99_absolute_error_tg':p99,
            'delay_max_absolute_error_tg':float(difference.max()),
            'max_keplerian_phase_error_rad':max_phase,
            'delay_range_tg':[float(high['delay_tg'][high['mask']].min()),
                              float(high['delay_tg'][high['mask']].max())],
            'visible_first_hit_layers':{str(k):int((high['hit_layer']==k).sum())
                                       for k in range(options['layers'])},
        })
    del tracer, coarse
    gc.collect()
    full = None
    if full_resolution:
        started = time.perf_counter()
        tracer = SchwarzschildDelayTracer(**{**options, 'width':2160, 'height':2000})
        full = {'detector_width':2160, 'detector_height':2000,
                'delivery_width':1080, 'delivery_height':1000,
                'supersampling':2, 'views':[]}
        for inclination in inclinations:
            maps = tracer.trace(inclination)
            valid = maps['mask']
            finite = all(np.isfinite(maps[key]).all() for key in ('r','phi','g','delay_tg'))
            symmetric = np.allclose(maps['delay_tg'],np.fliplr(maps['delay_tg']),atol=1e-5)
            invalid_zero = bool(np.all(maps['delay_tg'][~valid] == 0))
            in_disk = bool(np.all((maps['r'][valid] >= 6) & (maps['r'][valid] <= 28)))
            if not (finite and symmetric and invalid_zero and in_disk):
                raise AssertionError('production-sized map failed finite, symmetry or disk bounds check')
            full['views'].append({
                'inclination_deg':inclination,
                'finite_maps':finite, 'delay_reflection_symmetric':bool(symmetric),
                'invalid_delay_is_zero':invalid_zero, 'hits_within_disk':in_disk,
                'visible_first_hit_layers':{str(k):int((maps['hit_layer']==k).sum())
                                           for k in range(options['layers'])},
                'delay_range_tg':[float(maps['delay_tg'][valid].min()),
                                  float(maps['delay_tg'][valid].max())],
            })
            del maps
        full['elapsed_seconds'] = time.perf_counter()-started
        del tracer
        gc.collect()
    return {'scene':options, 'lookup_refinement':refinement,
            'full_resolution_smoke_test':full,
            'interpretation':'maximum phase errors refer to local Keplerian motion, Omega=r^(-3/2); '
                             'these tests do not resolve every higher-order image or validate gas dynamics'}


def diagnostics(full_resolution=False):
    """Analytic limits, coordinate-time integral, unchanged geometry, refinement."""
    started = time.perf_counter()
    # Minkowski space has u=sin(psi)/b, q=b*tan(psi/2), delay=-b*cot(psi).
    # This tests both incoming and outgoing branches of the regularization.
    b = 12.
    radial,q,h = _orbit_table_with_time(np.array([b]),2.7,.0005,mass=0)
    psi = np.arange(radial.shape[0])*h
    select = (psi>.1)&(psi<2.6)
    flat_radius = 1/radial[select,0].astype(np.float64)
    exact_q = b*np.tan(psi[select]/2)
    flat_q_relative = float(np.max(np.abs(q[select,0]/exact_q-1)))
    flat_delay_error = float(np.max(np.abs(q[select,0]-flat_radius+b/np.tan(psi[select]))))
    if flat_q_relative>2e-6 or flat_delay_error>3e-5:
        raise AssertionError('regularized time failed its exact flat-space limit')

    # On a finite segment compare the q formula against direct integration of
    # dt/dpsi. This uses a distinct, unregularized coordinate-time integrand.
    b = 8.
    radial,q,h = _orbit_table_with_time(np.array([b]),2.7,.0001)
    psi = np.arange(radial.shape[0])*h
    select = (psi>=.2)&(psi<=2.6)
    u = radial[select,0].astype(np.float64)
    radii = 1/u
    direct_elapsed = float(np.trapezoid(1/(b*u*u*(1-2*u)),psi[select]))
    elapsed_from_q = float(q[select,0][-1]-q[select,0][0]-
                           (tortoise_radius(radii[-1])-tortoise_radius(radii[0])))
    elapsed_relative = abs(elapsed_from_q/direct_elapsed-1)
    if elapsed_relative>2e-6:
        raise AssertionError('regularized delay disagrees with direct coordinate-time integral')

    # At large impact parameters the scaled delays approach Euclidean travel.
    weak_b=1.e5
    radial,q,h = _orbit_table_with_time(np.array([weak_b]),2.4,.0005)
    psi = np.arange(radial.shape[0])*h
    select = (psi>.2)&(psi<2.3)
    radii=1/radial[select,0].astype(np.float64)
    weak_scaled_error=float(np.max(np.abs(q[select,0]-tortoise_radius(radii)+
                                        weak_b/np.tan(psi[select])))/weak_b)
    if weak_scaled_error>.001:
        raise AssertionError('relative time failed its weak-field scaling limit')

    options=dict(width=192,height=200,fov=70,roll_deg=0,rin=6,rout=30,
                 layers=4,radial_samples=4096,angle_step=.004)
    baseline=SchwarzschildTracer(**options)
    reference=baseline.trace(80)
    del baseline
    gc.collect()
    tracer=SchwarzschildDelayTracer(**options)
    coarse=tracer.trace(80)
    exact_geometry=all(np.array_equal(reference[key],coarse[key])
                       for key in ('r','phi','g','mask','hit_layer','shadow'))
    symmetric=np.allclose(coarse['delay_tg'],np.fliplr(coarse['delay_tg']),atol=1e-5)
    invalid_zero=bool(np.all(coarse['delay_tg'][~coarse['mask']]==0))
    if not (exact_geometry and symmetric and invalid_zero):
        raise AssertionError('adding time changed baseline mapping, symmetry or invalid pixels')
    del tracer
    gc.collect()
    fine=SchwarzschildDelayTracer(**{**options,'radial_samples':8192,'angle_step':.002}).trace(80)
    common=coarse['mask']&fine['mask']&(coarse['hit_layer']==fine['hit_layer'])
    error=np.abs(coarse['delay_tg'][common]-fine['delay_tg'][common])
    phase_error=error/fine['r'][common]**1.5
    p99=float(np.percentile(error,99))
    max_phase=float(phase_error.max())
    if p99>.02 or max_phase>.01:
        raise AssertionError('time lookup has not converged at the tested resolution')
    return {
        'passed':True,'flat_space_q_max_relative_error':flat_q_relative,
        'flat_space_delay_max_absolute_error_tg':flat_delay_error,
        'finite_segment_direct_elapsed_tg':direct_elapsed,
        'finite_segment_regularized_elapsed_tg':elapsed_from_q,
        'finite_segment_elapsed_relative_error':elapsed_relative,
        'weak_field_delay_error_over_impact':weak_scaled_error,
        'baseline_geometry_bitwise_identical':exact_geometry,
        'delay_reflection_symmetric':bool(symmetric),'invalid_delay_is_zero':invalid_zero,
        'refinement_delay_p99_absolute_error_tg':p99,
        'refinement_delay_max_absolute_error_tg':float(error.max()),
        'refinement_max_keplerian_phase_error_rad':max_phase,
        'delay_range_tg':[float(fine['delay_tg'][fine['mask']].min()),
                          float(fine['delay_tg'][fine['mask']].max())],
        'visible_first_hit_layers':{str(k):int((fine['hit_layer']==k).sum()) for k in range(4)},
        'production_scene':_production_view_diagnostics(full_resolution=full_resolution),
        'scope':'tests of this idealized geodesic-delay model, not validation of an accretion flow',
        'elapsed_seconds':time.perf_counter()-started,
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('nasa-raytrace-diagnostics.json'))
    parser.add_argument('--production-grid',action='store_true',
                        help='also check 2160x2000 maps at 82, 83.5 and 85 degrees')
    args=parser.parse_args()
    result=diagnostics(full_resolution=args.production_grid)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
