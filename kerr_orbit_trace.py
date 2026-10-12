"""Kerr geodesics, disk azimuth, redshift and relative coordinate delays.

Original RK4 integration of Gralla & Lupsasca (2020), equations (6)-(9),(18):
https://arxiv.org/html/1910.12881v3
Prograde circular orbits and ISCO: Bardeen, Press & Teukolsky (1972),
https://doi.org/10.1086/151796, equations (2.16),(2.21).

Units G=M=c=1; positive black-hole spin and positive disk angular velocity
are aligned about +z. An opaque infinitesimally thin disk begins at the
prograde ISCO. The camera is a sequence of static observers at infinity;
its signed angle carries a continuous screen basis through both poles.

The radial equation uses u=1/r and backward Mino parameter tau:
 u'^2 = P^2-D*K*u^2; P=1+(a^2-a*lambda)u^2, D=1-2u+a^2u^2,
 K=eta+(lambda-a)^2.
The angular vector N is unit length and satisfies
 N''=-(eta+lambda^2+2*a^2*Nz^2)N+a^2*Nz*e_z.
Its angular momentum about z is -lambda. Therefore atan2(Ny,Nx) integrates
the backward angular term -lambda/sin(theta)^2 without a pole singularity.
The extra radial azimuth derivative is -a*u*(2-a*lambda*u)/D.

Elapsed backward coordinate time ell is regularized as
 q=ell+r*_Kerr(r)-r*_Kerr(observer). On the incoming branch,
 q'=(1+a^2*u^2)*K/(P+u')+a*lambda-a^2*(1-Nz^2),
 using P^2-u'^2=D*K*u^2 to avoid subtracting infinities.
On the outgoing branch use the equivalent direct difference. Return
 delay_tg=q-r*_Kerr(hit), t_emit=t_obs-delay_tg.
One common observer constant is removed for every pixel and camera view.
Negative relative delays are allowed, with invalid rays set to zero.
No finite-observer approximation, independent order time offset, fabricated
lensed ring, polar image reflection, or physical camera aberration is used.
"""
from __future__ import annotations
import ctypes
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
    a=float(spin)
    if not 0<=a<1:
        raise ValueError('aligned prograde spin must lie in [0,1)')
    z1=1+(1-a*a)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3))
    z2=math.sqrt(3*a*a+z1*z1)
    return 3+z2-math.sqrt((3-z1)*(3+z1+2*z2))

def tortoise_radius(radius,spin):
    """Kerr tortoise with the same a=0 additive convention as our old model."""
    r=np.asarray(radius,dtype=np.float64)
    rp=1+math.sqrt(1-spin*spin); rm=1-math.sqrt(1-spin*spin)
    if np.any(r<=rp):
        raise ValueError('Kerr tortoise radius is defined outside the outer horizon')
    return r+2*rp/(rp-rm)*np.log((r-rp)/2)-2*rm/(rp-rm)*np.log((r-rm)/2)

def _load_kernel():
    source=ROOT/'kerr_orbit_trace.cpp'
    binary=ROOT/('_kerr_orbit_trace.dylib' if sys.platform=='darwin' else '_kerr_orbit_trace.so')
    if not binary.exists() or binary.stat().st_mtime<source.stat().st_mtime:
        compiler=shutil.which('clang++') or shutil.which('g++')
        if not compiler:
            raise RuntimeError('Kerr tracing requires clang++ or g++ for the supplied C++ kernel')
        cmd=[compiler,'-dynamiclib' if sys.platform=='darwin' else '-shared',
             '-std=c++17','-O3','-fPIC','-pthread',str(source),'-o',str(binary)]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
    lib=ctypes.CDLL(str(binary)); function=lib.trace_kerr_orbit
    f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    f32=np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS')
    i16=np.ctypeslib.ndpointer(dtype=np.int16,flags='C_CONTIGUOUS')
    u8=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
    function.argtypes=[f64,f64,ctypes.c_longlong,*([ctypes.c_double]*6),
                       ctypes.c_int,ctypes.c_int,f32,f32,f32,f32,f32,i16,u8,f64]
    function.restype=None
    return lib,function

_LIBRARY=None

class KerrOrbitTracer:
    """Full-circle observer at infinity; trace every detector ray independently."""
    def __init__(self,width,height,fov=72,roll_deg=0,spin=.9,rin=None,rout=28,
                 angle_step=.012,max_angle=20*math.pi,threads=None):
        if not isinstance(width,int) or not isinstance(height,int) or min(width,height)<8:
            raise ValueError('width and height must be integers >=8')
        self.spin=float(spin)
        self.rin=isco_radius(spin) if rin is None else float(rin)
        if not (fov>0 and self.rin>=isco_radius(spin)-1e-10 and rout>self.rin
                and 0<angle_step<=.03 and max_angle>0):
            raise ValueError('invalid camera, disk, spin or integration setting')
        self.width,self.height=width,height
        self.fov,self.roll_deg=float(fov),float(roll_deg)
        self.rout,self.angle_step=float(rout),float(angle_step)
        self.max_steps=math.ceil(max_angle/angle_step)
        self.threads=int(threads or min(os.cpu_count() or 1,8))
        scale=fov/width
        ax=(np.arange(width,dtype=np.float64)+.5-width/2)*scale
        ay=-(np.arange(height,dtype=np.float64)+.5-height/2)*scale
        xd,yd=np.meshgrid(ax,ay); roll=math.radians(roll_deg)
        self.x=np.ascontiguousarray(np.cos(roll)*xd+np.sin(roll)*yd)
        self.y=np.ascontiguousarray(-np.sin(roll)*xd+np.cos(roll)*yd)
        self.metadata={
            'model':'Kerr null geodesics; first opaque prograde circular thin-disk intersection',
            'spin':self.spin,'spin_definition':'a*=cJ/(GM^2)',
            'spin_alignment':'black-hole angular momentum and disk angular momentum both +z; prograde',
            'width':width,'height':height,'fov_horizontal_rg':fov,
            'fov_vertical_rg':height*scale,'roll_deg':roll_deg,
            'rin_rg':self.rin,'rout_rg':rout,'outer_horizon_rg':1+math.sqrt(1-spin*spin),
            'scaled_mino_step':angle_step,'maximum_scaled_mino_parameter':self.max_steps*angle_step,
            'integrator':'double-precision RK4; Hermite equatorial root followed by partial RK4 event step',
            'observer':'static at infinity','threads':self.threads,
            'camera_angle_convention':'signed degrees about disk x axis; n=(0,-sin(theta),cos(theta))',
            'camera_screen_basis':'e_x=(1,0,0); e_y=(0,cos(theta),sin(theta)); continuous through both poles',
            'disk_plane_crossing':'observer-at-infinity crossing excluded; coplanar rays first enter outer rim',
            'delay_units':'t_g=GM/c^3',
            'delay_definition':'elapsed coordinate propagation time minus one common observer Kerr r*; q-r*(hit)',
            'invalid_delay_value':0,'emission_time_rule':'t_em=t_obs-delay_tg; one common additive origin',
            'azimuth':'pole-regular angular vector plus backward Kerr radial frame-dragging integral',
            'disk_surfaces':'same prescribed isotropic surface emission on both faces',
        }

    def trace(self,inclination_deg):
        global _LIBRARY
        angle=float(inclination_deg)
        if not math.isfinite(angle):
            raise ValueError('observer angle must be finite')
        if _LIBRARY is None:
            _LIBRARY=_load_kernel()
        inc=math.radians(math.remainder(angle,360.)); si,ci=math.sin(inc),math.cos(inc)
        if abs(si)<2e-15: si=0.
        if abs(ci)<2e-15: ci=0.
        started=time.perf_counter();shape=(self.height,self.width)
        r=np.empty(shape,np.float32); phi=np.empty_like(r); g=np.empty_like(r)
        delay=np.empty_like(r); mu=np.empty_like(r)
        layer=np.empty(shape,np.int16);status=np.empty(shape,np.uint8)
        residual=np.empty(shape,np.float64)
        _LIBRARY[1](self.x,self.y,self.width*self.height,self.spin,si,ci,self.rin,self.rout,
                    self.angle_step,self.max_steps,self.threads,r,phi,g,delay,mu,layer,status,residual)
        if np.any(status>=3):
            raise RuntimeError(f'Kerr integration left {np.count_nonzero(status>=3)} unresolved rays; refine the integration budget')
        for key,values in [('radius',r),('azimuth',phi),('shift',g),('delay',delay)]:
            if not np.isfinite(values).all():
                raise RuntimeError(f'Kerr {key} map contains nonfinite values')
        mask=status==1
        self.metadata['last_trace_seconds']=time.perf_counter()-started
        self.metadata['last_view_angle_deg']=angle
        self.metadata['last_constraint_residual_max']=float(residual.max())
        return {'r':r,'phi':phi,'g':g,'delay_tg':delay,'mu':mu,'mask':mask,
                'hit_layer':layer,'shadow':status==2,'status':status,'constraint_residual':residual}
