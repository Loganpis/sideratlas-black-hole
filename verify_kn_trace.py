#!/usr/bin/env python3
"""Independent limits, capture, symmetry, delay quadrature and RK4 refinement QA.

Scope is numerical consistency of neutral rays in the specified idealized
Kerr-Newman geometry. This is not validation of an accreting charged plasma.
"""
from __future__ import annotations
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from kn_orbit_trace import KerrNewmanOrbitTracer,tortoise_radius,antipodal_maps
from kn_orbits import isco_radius,horizon_radius,orbital_quantities
ROOT=Path(__file__).resolve().parent
RESULT={'scope':__doc__,'checks':[]}
def check(name,passed,**metrics):
    RESULT['checks'].append(dict(name=name,passed=bool(passed),**metrics))
    print(('PASS ' if passed else 'FAIL ')+name,flush=True)
def angular_error(a,b):return np.abs(np.angle(np.exp(1j*(a.astype(float)-b))))
def err(a,b,key,shared):
    e=angular_error(a[key][shared],b[key][shared]) if key=='phi' else abs(a[key][shared]-b[key][shared])
    return np.percentile(e,[50,99,100]).tolist()

def compare_limits():
    from nasa_raytrace import SchwarzschildDelayTracer
    from kerr_raytrace import KerrTracer
    kw=dict(width=144,height=120,fov=72,roll_deg=0,rin=6,rout=28)
    baseline=SchwarzschildDelayTracer(**kw,radial_samples=16384,angle_step=.001)
    zero=KerrNewmanOrbitTracer(**kw,spin=0,charge_q=0)
    for angle in (83,0,-90,-97,-180,-270):
        a,b=baseline.trace(angle),zero.trace(angle)
        common=a['mask']&b['mask']&(a['hit_layer']==b['hit_layer'])
        errors={k:err(a,b,k,common) for k in ('r','phi','g','delay_tg')}
        mismatch=float(np.mean(a['mask']!=b['mask']))
        # Legacy table interpolation has largest errors immediately adjacent
        # to the critical curve; report maxima but constrain the 99th percentile.
        check(f'Schwarzschild_limit_{angle}',mismatch<.001 and errors['r'][1]<.001 and
              errors['phi'][1]<1e-5 and errors['g'][1]<1e-4 and errors['delay_tg'][1]<.002,
              mask_disagreement=mismatch,absolute_errors_p50_p99_max=errors)
    for angle in (45,83):
        kw=dict(width=180,height=150,fov=72,roll_deg=0,spin=.9,rout=28,angle_step=.006)
        a=KerrTracer(**kw).trace(angle); b=KerrNewmanOrbitTracer(**kw,charge_q=0).trace(angle)
        common=a['mask']&b['mask']
        errors={k:err(a,b,k,common) for k in ('r','g','mu')}
        check(f'Kerr_q0_independent_scalar_latitude_{angle}',np.array_equal(a['mask'],b['mask']) and
              all(e[2]<2e-5 for e in errors.values()),absolute_errors_p50_p99_max=errors)

def capture_tests():
    # The far-away virtual annulus excludes disk hits, revealing the geodesic
    # capture set. Capture is independently determined by polynomial roots.
    for spin,charge in ((0,.4),(.9,.4)):
        t=KerrNewmanOrbitTracer(80,64,fov=16,spin=spin,charge_q=charge,rin=10000,rout=10001)
        for angle in (0,61,90,-119):
            m=t.trace(angle);si=math.sin(math.radians(angle));ci=math.cos(math.radians(angle))
            lam=-t.x*si; eta=t.y*t.y+(t.x*t.x-spin*spin)*ci*ci
            K=eta+(lam-spin)**2; A=spin*spin-lam*lam-eta
            predicted=np.empty(m['mask'].shape,bool)
            for ind in np.ndindex(predicted.shape):
                coeff=[1,0,A[ind],2*K[ind],-(spin*spin*eta[ind]+charge*charge*K[ind])]
                roots=np.roots(coeff)
                predicted[ind]=not any(abs(x.imag)<1e-7 and x.real>horizon_radius(spin,charge)+1e-7 for x in roots)
            check(f'quartic_capture_a{spin}_q{charge}_view{angle}',
                  np.array_equal(predicted,m['shadow']),disagreeing_pixels=int(np.sum(predicted!=m['shadow'])))
            if spin==0:
                rph=(3+math.sqrt(9-8*charge*charge))/2
                bcrit=rph/math.sqrt(1-2/rph+charge*charge/rph**2)
                expected=t.x*t.x+t.y*t.y<bcrit*bcrit
                check(f'RN_analytic_photon_sphere_view{angle}',np.array_equal(m['shadow'],expected),
                      photon_sphere_rg=rph,critical_impact_rg=bcrit)

def refinement_and_symmetry():
    coarse=KerrNewmanOrbitTracer(240,200,angle_step=.012)
    fine=KerrNewmanOrbitTracer(240,200,angle_step=.006)
    results={}
    for angle in (83,0,-90,-97,-180,-270,89.999,90.001):
        a,b=coarse.trace(angle),fine.trace(angle);results[angle]=a
        shared=a['mask']&b['mask']&(a['hit_layer']==b['hit_layer'])
        errors={k:err(a,b,k,shared) for k in ('r','phi','g','delay_tg')}
        residual=float(a['constraint_residual'].max())
        check(f'RK4_refinement_{angle}',np.array_equal(a['mask'],b['mask']) and
              all(e[2]<.002 for e in errors.values()) and residual<2e-7,
              errors_p50_p99_max=errors,constraint_residual_max=residual)
        check(f'physical_and_finite_maps_{angle}',all(np.isfinite(a[k]).all() for k in ('r','phi','g','delay_tg','mu'))
              and np.all(a['r'][a['mask']]>=isco_radius()-.000001)
              and np.all((a['mu'][a['mask']]>=0)&(a['mu'][a['mask']]<=1.000001))
              and np.all(a['g'][a['mask']]>0) and np.all(a['delay_tg'][~a['mask']]==0))
    for angle in (83,0,-90,89.999,90.001):
        original=results[angle];opposite=coarse.trace(angle-180); transformed=antipodal_maps(original)
        exact=all(np.array_equal(opposite[k],transformed[k]) for k in ('mask','hit_layer','status','r','g','delay_tg','mu'))
        phierror=float(angular_error(opposite['phi'],transformed['phi']).max())
        check(f'antipodal_exact_scalar_maps_{angle}',exact and phierror<5e-7,azimuth_max_radians=phierror)
    neg=KerrNewmanOrbitTracer(240,200,charge_q=-.4).trace(83)
    check('charge_sign_even_neutral_geodesics',all(np.array_equal(neg[k],results[83][k]) for k in neg))
    periodic=coarse.trace(83-360)
    check('full_orbit_returns_exactly',all(np.array_equal(periodic[k],results[83][k]) for k in periodic))
    for angle in (0,-180):
        m=results[angle]; sl=slice(20,220)
        check(f'pole_axisymmetry_{angle}',all(np.array_equal(m[k][:,sl],np.rot90(m[k][:,sl])) for k in ('mask','r','g','delay_tg')))
        expected=1/orbital_quantities(m['r'][m['mask']])['ut']
        check(f'pole_time_dilation_{angle}',np.max(abs(m['g'][m['mask']]-expected))<2e-6)
    odd=KerrNewmanOrbitTracer(81,65)
    for angle in (-90,90):
        m=odd.trace(angle); row=m['r'][32]
        # A coplanar ray misses the finite disk if its outermost radial
        # turning point is outside the outer rim. Determine this separately
        # from quartic roots instead of assuming every detector ray hits it.
        predicted=[]
        for alpha in odd.x[32]:
            lam=-alpha*math.sin(math.radians(angle));K=(lam-.9)**2
            roots=np.roots([1,0,.81-lam*lam,2*K,-.16*K])
            turning=[z.real for z in roots if abs(z.imag)<1e-8 and z.real>horizon_radius(.9,.4)]
            predicted.append(not turning or max(turning)<=28)
        predicted=np.asarray(predicted)
        check(f'exact_coplanar_rays_hit_outer_rim_{angle}',
              np.array_equal(m['mask'][32],predicted) and np.all(row[predicted]==28)
              and np.all(m['hit_layer'][32,predicted]==0),hit_pixels=int(predicted.sum()),
              miss_pixels=int((~predicted).sum()))
    rr=np.geomspace(horizon_radius(.9,.4)+.02,1e3,200);h=1e-6*rr
    derivative=(tortoise_radius(rr+h,.9,.4)-tortoise_radius(rr-h,.9,.4))/(2*h)
    analytic=(rr*rr+.81)/(rr*rr-2*rr+.81+.16)
    check('KN_tortoise_derivative',np.max(abs(derivative/analytic-1))<2e-8,max_relative_error=float(np.max(abs(derivative/analytic-1))))

# Independently integrate scalar latitude plus angular quadratures. No angular
# vector, radial ODE, kernel delay derivative or kernel azimuth derivative used.
def angular_integrals(alpha,beta,inc,crossing,a=.9,step=.0001):
    si,ci=math.sin(inc),math.cos(inc);lam=-alpha*si;eta=beta*beta+(alpha*alpha-a*a)*ci*ci
    A=a*a-lam*lam-eta; B=math.sqrt(alpha*alpha+beta*beta+1)
    y=np.array([ci,beta*si/B,-math.pi/2,0.]);s=0.;n=0
    def f(y):
        z,w,phi,t=y
        return np.array([w,(A*z-2*a*a*z**3)/(B*B),-lam/(1-z*z)/B,(a*lam-a*a+a*a*z*z)/B])
    def rk(y,h):
        k1=f(y);k2=f(y+h*k1/2);k3=f(y+h*k2/2);k4=f(y+h*k3)
        return y+h/6*(k1+2*k2+2*k3+k4)
    while s<60:
        nxt=rk(y,step)
        if y[0]*nxt[0]<0:
            if n==crossing:
                lo,hi=0.,step
                for _ in range(30):
                    mid=(lo+hi)/2; part=rk(y,mid)
                    if part[0]*y[0]>0:lo=mid
                    else:hi=mid
                hit=rk(y,(lo+hi)/2)
                return (s+(lo+hi)/2)/B,hit[2],hit[3]
            n+=1
        y=nxt;s+=step
    raise RuntimeError('angular reference never reached requested crossing')

_GX,_GW=np.polynomial.legendre.leggauss(256)
def gauss(func,left,right):
    x=(left+right)/2+(right-left)/2*_GX
    return float((right-left)/2*np.dot(_GW,func(x)))

def radial_quadratures(alpha,beta,inc,rhit,target_tau,a=.9,q=.4):
    lam=-alpha*math.sin(inc); eta=beta*beta+(alpha*alpha-a*a)*math.cos(inc)**2
    A=a*a-lam*lam-eta;K=eta+(lam-a)**2;C=a*a*eta+q*q*K;uhit=1/rhit
    def U(u):return 1+A*u*u+2*K*u**3-C*u**4
    roots=np.roots([-C,2*K,A,0,1])
    turns=[x.real for x in roots if abs(x.imag)<1e-8 and x.real>0 and x.real<1/horizon_radius(a,q)]
    ut=min(turns) if turns else None
    def base(u,kind):
        root=np.sqrt(np.maximum(U(u),1e-30));P=1+(a*a-a*lam)*u*u;D=1-2*u+(a*a+q*q)*u*u;H=1+a*a*u*u
        if kind=='mino':return 1/root
        if kind=='phi':return -a*u*(2-(q*q+a*lam)*u)/D/root
        if kind=='qin':return H*K/(P+root)/root
        if kind=='qout':return H*(P+root)/(u*u*D)/root
    if ut is None:
        integrate=lambda kind,l,r:gauss(lambda u:base(u,kind),l,r)
    else:
        def integrate(kind,l,r):
            zl=math.sqrt(max(0,1-r/ut));zh=math.sqrt(max(0,1-l/ut))
            return gauss(lambda z:base(ut*(1-z*z),kind)*2*ut*z,zl,zh)
    incoming_tau=integrate('mino',0,uhit)
    outgoing=ut is not None and abs(target_tau-incoming_tau)>1e-5
    if not outgoing:return incoming_tau,integrate('phi',0,uhit),integrate('qin',0,uhit),False
    total_tau=integrate('mino',0,ut)+integrate('mino',uhit,ut)
    phi=integrate('phi',0,ut)+integrate('phi',uhit,ut)
    regularized=integrate('qin',0,ut)+integrate('qout',uhit,ut)
    return total_tau,phi,regularized,True

def independent_delay_and_azimuth():
    tracer=KerrNewmanOrbitTracer(144,120,angle_step=.006)
    maps=tracer.trace(83);inc=math.radians(83)
    candidates=np.argwhere(maps['mask']&(abs(tracer.x)>2)&(abs(tracer.x)<10))
    chosen=[]
    for layer in (0,1):
        rows=[tuple(i) for i in candidates if maps['hit_layer'][tuple(i)]==layer]
        for ix in np.linspace(0,len(rows)-1,min(5,len(rows))).astype(int):chosen.append(rows[ix])
    diagnostics=[]
    for ind in chosen:
        alpha,beta=float(tracer.x[ind]),float(tracer.y[ind]);layer=int(maps['hit_layer'][ind]);r=float(maps['r'][ind])
        tau,phiangular,tangular=angular_integrals(alpha,beta,inc,layer)
        radialtau,phirad,qrad,outgoing=radial_quadratures(alpha,beta,inc,r,tau)
        predicted_delay=qrad+tangular-float(tortoise_radius(r,.9,.4))
        predicted_phi=phirad+phiangular
        derr=abs(predicted_delay-float(maps['delay_tg'][ind]));perr=float(angular_error(np.array(predicted_phi),maps['phi'][ind]))
        terr=abs(tau-radialtau)
        diagnostics.append(dict(pixel=[int(i) for i in ind],image_order=layer,outgoing=outgoing,delay_error_tg=derr,
                                azimuth_error_radians=perr,mino_arrival_error=terr))
        check(f'independent_separated_quadratures_pixel_{ind}',derr<.0003 and perr<3e-5 and terr<2e-6,**diagnostics[-1])
    RESULT['separated_quadrature_samples']=diagnostics
    check('quadratures_include_outgoing_and_higher_order_paths',any(d['outgoing'] for d in diagnostics) and any(d['image_order']>=1 for d in diagnostics))

def production_detector_benchmark():
    tracer=KerrNewmanOrbitTracer(2160,1800,angle_step=.012)
    maps=tracer.trace(83)
    check('production_detector_resolved_and_finite',maps['mask'].any() and
          all(np.isfinite(maps[k]).all() for k in ('r','phi','g','delay_tg','mu')) and
          float(maps['constraint_residual'].max())<2e-7,
          width=2160,height=1800,seconds=tracer.metadata['last_trace_seconds'],
          hit_pixels=int(maps['mask'].sum()),max_constraint_residual=float(maps['constraint_residual'].max()))

def main():
    started=time.perf_counter()
    for fn in (compare_limits,capture_tests,refinement_and_symmetry,independent_delay_and_azimuth,production_detector_benchmark):
        try:fn()
        except Exception as e:check(fn.__name__+'_completed',False,error=repr(e))
    RESULT['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('kn_orbit_trace.py','kn_orbit_trace.cpp','kn_orbits.py','verify_kn_trace.py')}
    RESULT['elapsed_seconds']=time.perf_counter()-started
    RESULT['passed']=all(x['passed'] for x in RESULT['checks'])
    (ROOT/'kn-orbit-trace-diagnostics.json').write_text(json.dumps(RESULT,indent=2)+'\n')
    print('PASSED' if RESULT['passed'] else 'FAILED',len(RESULT['checks']),flush=True)
    return 0 if RESULT['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
