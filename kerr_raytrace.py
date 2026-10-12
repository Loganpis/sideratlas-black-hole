"""True Kerr null geodesics for a steady, opaque equatorial accretion disk.

This original implementation uses the exact separated Kerr potentials from
Gralla & Lupsasca (2020), PRD101 044032, equations 7, 8 and 18:
https://doi.org/10.1103/PhysRevD.101.044032
https://arxiv.org/abs/1910.12881

Prograde circular orbit and ISCO expressions: Bardeen, Press & Teukolsky
(1972), ApJ178 347, especially equations 2.16 and 2.21:
https://doi.org/10.1086/151796
https://adsabs.harvard.edu/pdf/1972ApJ...178..347B

It integrates u=1/r and z=cos(theta) from an observer at infinity, using
fourth-order Runge--Kutta and cubic-Hermite disk intersection locations.
For backward Mino parameter tau:
 u'^2 = 1 + A*u^2 + 2*K*u^3 - a^2*eta*u^4,
 z'^2 = eta + A*z^2 - a^2*z^4,
 A=a^2-lambda^2-eta, K=eta+(lambda-a)^2.
Their second-order derivatives are regular at radial and polar turning points.
Screen constants: lambda=-alpha*sin(i), eta=beta^2+(alpha^2-a^2)*cos(i)^2.
u(0)=0, u'(0)=1, z(0)=cos(i), z'(0)=beta*sin(i).

The first intersection at r_ISCO<=r<=rout is opaque. Earlier intersections in
the central hole are traversed, so higher-order disk images emerge naturally.
Emission-angle cosine in the orbiting fluid frame: mu_e=g*sqrt(eta)/r.
g=1/[u^t(1-Omega*lambda)], Omega=1/(r^(3/2)+a).

No fluid texture, camera glow or spectral model is introduced here. This
axisymmetric renderer does not integrate azimuth or coordinate travel time;
neither enters its stationary radial emission model. Distances use GM/c^2.
The C++ kernel compiles locally on first import using clang++/g++ and uses
standard threads. No downloaded code or Numba dependency is required.
"""
from __future__ import annotations
import argparse
import ctypes
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parent

def isco_radius(spin):
    """Prograde equatorial ISCO (Bardeen, Press & Teukolsky 1972, Eq.2.21)."""
    a=float(spin)
    if not 0<=a<1: raise ValueError('spin must lie in [0,1)')
    z1=1+(1-a*a)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3))
    z2=math.sqrt(3*a*a+z1*z1)
    return 3+z2-math.sqrt((3-z1)*(3+z1+2*z2))

def _load_kernel():
    source=ROOT/'kerr_raytrace.cpp'
    binary=ROOT/('_kerr_raytrace.dylib' if sys.platform=='darwin' else '_kerr_raytrace.so')
    if not binary.exists() or binary.stat().st_mtime<source.stat().st_mtime:
        compiler=shutil.which('clang++') or shutil.which('g++')
        if not compiler: raise RuntimeError('Kerr tracing requires clang++ or g++ to compile its supplied C++ kernel.')
        command=[compiler,'-std=c++17','-O3','-fPIC','-pthread',str(source),'-o',str(binary)]
        command.insert(1,'-dynamiclib' if sys.platform=='darwin' else '-shared')
        subprocess.run(command,check=True,capture_output=True,text=True)
    lib=ctypes.CDLL(str(binary))
    function=lib.trace_kerr
    f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    f32=np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS')
    i16=np.ctypeslib.ndpointer(dtype=np.int16,flags='C_CONTIGUOUS')
    u8=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
    function.argtypes=[f64,f64,ctypes.c_longlong,*([ctypes.c_double]*5),ctypes.c_int,ctypes.c_int,f32,f32,f32,i16,u8,f64]
    function.restype=None
    return lib,function

_LIBRARY=None

class KerrTracer:
    """Reusable detector; independently integrate every ray for every view.

    fov is the full horizontal impact-plane span in GM/c^2, with square
    isotropic detector pixels. roll_deg rotates the image counterclockwise.
    Defaults resolve the near-critical images without inventing a bright rim.
    """
    def __init__(self,width,height,fov=30,roll_deg=20,spin=.9,rin=None,rout=500,
                 angle_step=.006,max_angle=20*math.pi,threads=None):
        if not isinstance(width,int) or not isinstance(height,int) or min(width,height)<8:
            raise ValueError('width and height must be integers >=8')
        self.spin=float(spin)
        self.rin=isco_radius(spin) if rin is None else float(rin)
        if not (fov>0 and 0<=spin<1 and self.rin>=isco_radius(spin)-1e-10 and rout>self.rin and 0<angle_step<=.03):
            raise ValueError('invalid camera, disk, spin or integration setting')
        self.width,self.height=width,height
        self.fov,self.roll_deg=float(fov),float(roll_deg)
        self.rout,self.angle_step=float(rout),float(angle_step)
        self.max_steps=math.ceil(max_angle/angle_step)
        self.threads=int(threads or min(os.cpu_count() or 1,8))
        scale=fov/width
        ax=(np.arange(width,dtype=np.float64)+.5-width/2)*scale
        ay=-(np.arange(height,dtype=np.float64)+.5-height/2)*scale
        xd,yd=np.meshgrid(ax,ay)
        roll=math.radians(roll_deg)
        self.x=np.ascontiguousarray(np.cos(roll)*xd+np.sin(roll)*yd)
        self.y=np.ascontiguousarray(-np.sin(roll)*xd+np.cos(roll)*yd)
        self.metadata={'model':'Kerr null geodesics; first opaque prograde circular thin-disk intersection',
                       'spin':spin,'width':width,'height':height,'fov_horizontal_rg':fov,
                       'fov_vertical_rg':height*scale,'roll_deg':roll_deg,'rin_rg':self.rin,
                       'rout_rg':rout,'scaled_mino_step':angle_step,'maximum_scaled_mino_parameter':max_angle,
                       'integrator':'RK4 with cubic-Hermite equatorial intersections',
                       'observer':'static at infinity','threads':self.threads,
                       'azimuth_and_travel_time':'not required by stationary axisymmetric emission; not integrated'}

    def trace(self,inclination_deg):
        global _LIBRARY
        if not 0<inclination_deg<89.9: raise ValueError('inclination must lie in (0,89.9)')
        if _LIBRARY is None: _LIBRARY=_load_kernel()
        started=time.perf_counter()
        shape=(self.height,self.width)
        r=np.empty(shape,np.float32); g=np.empty_like(r); mu=np.empty_like(r)
        layer=np.empty(shape,np.int16); status=np.empty(shape,np.uint8);residual=np.empty(shape,np.float64)
        _LIBRARY[1](self.x,self.y,self.width*self.height,self.spin,math.radians(inclination_deg),
                    self.rin,self.rout,self.angle_step,self.max_steps,self.threads,r,g,mu,layer,status,residual)
        if np.any(status>=3):
            raise RuntimeError(f'Kerr integration left {np.count_nonzero(status>=3)} unresolved rays; refine integration budget.')
        mask=status==1
        self.metadata['last_trace_seconds']=time.perf_counter()-started
        self.metadata['last_inclination_deg']=inclination_deg
        self.metadata['last_constraint_residual_max']=float(residual.max())
        return {'r':r,'g':g,'mu':mu,'mask':mask,'hit_layer':layer,'shadow':status==2,
                'status':status,'constraint_residual':residual}

def diagnostics():
    from raytrace import SchwarzschildTracer
    started=time.perf_counter()
    commonargs=dict(width=192,height=256,fov=32,roll_deg=0,rin=6,rout=40)
    reference=SchwarzschildTracer(**commonargs,radial_samples=16384,angle_step=.001).trace(68)
    zero=KerrTracer(**commonargs,spin=0,angle_step=.006).trace(68)
    shared=reference['mask']&zero['mask']&(reference['hit_layer']==zero['hit_layer'])
    radius_error=np.abs(reference['r'][shared]-zero['r'][shared])/reference['r'][shared]
    g_error=np.abs(reference['g'][shared]-zero['g'][shared])/reference['g'][shared]
    disagreed=float(np.mean(reference['mask']!=zero['mask']))
    if disagreed>.002 or np.percentile(radius_error,99)>.001: raise AssertionError('Kerr a=0 does not reproduce Schwarzschild')
    if not np.allclose(zero['r'],np.fliplr(zero['r']),atol=1e-5): raise AssertionError('a=0 geometry not symmetric')
    coarse_tracer=KerrTracer(256,320,spin=.9,roll_deg=21,angle_step=.006)
    coarse=coarse_tracer.trace(72)
    fine=KerrTracer(256,320,spin=.9,roll_deg=21,angle_step=.003).trace(72)
    shared=coarse['mask']&fine['mask']&(coarse['hit_layer']==fine['hit_layer'])
    rr=np.abs(coarse['r'][shared]-fine['r'][shared])/fine['r'][shared]
    gg=np.abs(coarse['g'][shared]-fine['g'][shared])/fine['g'][shared]
    disagreement=float(np.mean(coarse['mask']!=fine['mask']))
    def flux(m):
        use=m['mask'];r=m['r'][use]
        return float(np.sum(m['g'][use]**4*(1-np.sqrt(isco_radius(.9)/r))/r**3))
    flux_error=abs(flux(coarse)-flux(fine))/flux(fine)
    mask=coarse['mask']
    if disagreement>.002 or np.percentile(rr,99)>.001 or flux_error>.002: raise AssertionError('Kerr refinement failed')
    if not np.all((coarse['mu'][mask]>0)&(coarse['mu'][mask]<=1.0001)): raise AssertionError('nonphysical fluid-frame emission angle')
    if not np.any(coarse['hit_layer']>=1): raise AssertionError('higher-order disk image missing')
    if coarse['constraint_residual'].max()>1e-5: raise AssertionError('Kerr potential constraint drift')
    # Independent capture test: classify the real roots of the quartic radial
    # potential instead of following an ODE. With no nearby emitting disk,
    # a ray escapes iff its radial potential has a turning point outside r+.
    capture_tracer=KerrTracer(128,160,spin=.9,fov=15,roll_deg=0,rin=1e5,rout=1e5+1)
    capture=capture_tracer.trace(72)
    lam=-capture_tracer.x*math.sin(math.radians(72))
    eta=capture_tracer.y**2+(capture_tracer.x**2-.9**2)*math.cos(math.radians(72))**2
    A=.9**2-lam**2-eta
    K=eta+(lam-.9)**2
    horizon=1+math.sqrt(1-.9**2)
    analytic_capture=np.empty(lam.size,dtype=bool)
    for j,(aa,kk,qq) in enumerate(zip(A.ravel(),K.ravel(),eta.ravel())):
        roots=np.roots([1,0,aa,2*kk,-.9**2*qq])
        turns=(np.abs(roots.imag)<1e-7)&(roots.real>horizon+1e-7)
        analytic_capture[j]=not np.any(turns)
    analytic_capture=analytic_capture.reshape(lam.shape)
    capture_disagreement=float(np.mean(analytic_capture!=capture['shadow']))
    if capture_disagreement>1e-4: raise AssertionError('Kerr silhouette disagrees with analytic radial turning points')
    high_inclination_checks={}
    for inc in [80,84]:
        hi_coarse=KerrTracer(384,512,spin=.9,fov=26,roll_deg=23,angle_step=.006).trace(inc)
        hi_fine=KerrTracer(384,512,spin=.9,fov=26,roll_deg=23,angle_step=.003).trace(inc)
        both=hi_coarse['mask']&hi_fine['mask']&(hi_coarse['hit_layer']==hi_fine['hit_layer'])
        radius_rel=np.abs(hi_coarse['r'][both]-hi_fine['r'][both])/hi_fine['r'][both]
        mask_diff=float(np.mean(hi_coarse['mask']!=hi_fine['mask']))
        flux_diff=abs(flux(hi_coarse)-flux(hi_fine))/flux(hi_fine)
        p99=float(np.percentile(radius_rel,99))
        if mask_diff>.002 or p99>.001 or flux_diff>.002: raise AssertionError(f'Kerr inclination {inc} refinement failed')
        high_inclination_checks[str(inc)]={
            'mask_disagreement':mask_diff,'radius_p99_relative_error':p99,
            'integrated_test_flux_relative_error':flux_diff,
            'maximum_potential_constraint_residual':float(hi_coarse['constraint_residual'].max()),
            'mu_max':float(hi_coarse['mu'].max()),
            'highest_visible_disk_crossing':int(hi_coarse['hit_layer'].max()),
        }
    return {'passed':True,'isco_spin0':isco_radius(0),'isco_spin09':isco_radius(.9),
            'schwarzschild_mask_disagreement':disagreed,
            'schwarzschild_radius_p99_relative_error':float(np.percentile(radius_error,99)),
            'schwarzschild_g_p99_relative_error':float(np.percentile(g_error,99)),
            'kerr_refinement_mask_disagreement':disagreement,
            'kerr_refinement_radius_p99_relative_error':float(np.percentile(rr,99)),
            'kerr_refinement_g_p99_relative_error':float(np.percentile(gg,99)),
            'kerr_refinement_integrated_test_flux_relative_error':flux_error,
            'potential_constraint_max_relative_residual':float(coarse['constraint_residual'].max()),
            'kerr_capture_mask_disagreement_with_quartic_roots':capture_disagreement,
            'high_inclination_refinement_checks':high_inclination_checks,
            'emission_mu_range':[float(coarse['mu'][mask].min()),float(coarse['mu'][mask].max())],
            'visible_layers':{str(k):int((coarse['hit_layer']==k).sum()) for k in np.unique(coarse['hit_layer'][mask])},
            'coarse_metadata':coarse_tracer.metadata,'elapsed_seconds':time.perf_counter()-started}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='kerr-raytrace-diagnostics.json')
    args=parser.parse_args()
    result=diagnostics()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
