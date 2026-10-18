#!/usr/bin/env python3
"""Numerical validation of layered Kerr–Newman null rays, independent of gas art.

Checks preserve the established opaque-disk solution, compare null capture
with independent polynomial roots and compare hit azimuth/delay with separated
quadratures. Critical-curve samples follow exact spherical null-orbit constants;
their multiple image orders emerge from integration, without ring templates.
These tests validate the specified geometric/transfer inputs, not a GRMHD disk.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from kn_accretion_trace import KerrNewmanAccretionTracer
from kn_orbit_trace import KerrNewmanOrbitTracer, tortoise_radius, antipodal_maps
from kn_orbits import horizon_radius
from verify_kn_trace import angular_integrals, radial_quadratures

ROOT=Path(__file__).resolve().parent
RESULT={'scope':__doc__, 'checks':[]}

def check(name,passed,**metrics):
    RESULT['checks'].append(dict(name=name,passed=bool(passed),**metrics))
    print(('PASS ' if passed else 'FAIL ')+name,flush=True)

def angle_error(a,b):
    return np.abs(np.angle(np.exp(1j*(a.astype(float)-b))))

def compatibility():
    old=KerrNewmanOrbitTracer(160,132,fov=28)
    new=KerrNewmanAccretionTracer(160,132,fov=28,max_hits=8)
    for view in (83,0,-90,-97,-180,-270,89.999,90.001):
        a,b=old.trace(view),new.trace(view)
        disk=b['hit_kind']==1
        found=disk.any(axis=-1)
        ix=np.argmax(disk,axis=-1)
        errors={}
        common=a['mask']&found
        for key in ('r','phi','delay_tg'):
            value=np.take_along_axis(b[key],ix[...,None],axis=-1)[...,0]
            d=angle_error(value[common],a[key][common]) if key=='phi' else abs(value[common]-a[key][common])
            errors[key]=float(d.max()) if d.size else 0.
        check(f'opaque_outer_disk_compatibility_view_{view}',np.array_equal(a['mask'],found) and
              errors['r']<.002 and errors['phi']<.001 and errors['delay_tg']<.003,
              max_errors=errors,plunge_intersections=int(np.sum(b['hit_kind']==2)))
        check(f'events_outside_black_horizon_view_{view}',
              np.all(b['r'][b['mask']]>=new.rmin-1e-7) and np.all(np.diff(b['hit_layer'],axis=-1)[b['mask'][...,1:]]>0),
              minimum_emitting_radius_rg=float(b['r'][b['mask']].min()),cutoff_rg=new.rmin)
        check(f'opaque_disk_is_final_event_view_{view}',
              np.all(np.sum(disk,axis=-1)<=1) and not np.any((np.cumsum(disk,axis=-1)>0)&(b['hit_kind']==2)))

def capture_and_null_momentum():
    for spin,charge in ((0,0),(0,.4),(.9,0),(.9,.4)):
        tracer=KerrNewmanAccretionTracer(48,40,fov=16,spin=spin,charge_q=charge,rin=10000,rout=10001,max_hits=12)
        m=tracer.trace(61)
        lam,eta=m['lambda'],m['eta'];K=eta+(lam-spin)**2;A=spin*spin-lam*lam-eta
        expected=np.empty(lam.shape,bool)
        for ix in np.ndindex(lam.shape):
            roots=np.roots([1,0,A[ix],2*K[ix],-(spin*spin*eta[ix]+charge*charge*K[ix])])
            expected[ix]=not any(abs(v.imag)<1e-7 and v.real>horizon_radius(spin,charge)+1e-7 for v in roots)
        check(f'independent_radial_quartic_capture_a{spin}_q{charge}',np.array_equal(expected,m['shadow']))
        mask=m['mask'];r=m['r'][mask].astype(float)
        lam=np.broadcast_to(lam[...,None],m['r'].shape)[mask]
        eta=np.broadcast_to(eta[...,None],m['r'].shape)[mask]
        kr=m['kr_over_E'][mask].astype(float)
        D=r*r-2*r+spin*spin+charge*charge
        tt=-((r*r+spin*spin)**2-spin*spin*D)/(r*r*D)
        cross=2*spin*(2*r-charge*charge)*lam/(r*r*D)
        pp=(D-spin*spin)*lam*lam/(r*r*D)
        rr=D/(r*r)*kr*kr;th=eta/(r*r)
        norm=np.abs(tt+cross+pp+rr+th)/(1+abs(tt)+abs(cross)+abs(pp)+abs(rr)+abs(th))
        check(f'independent_inverse_metric_photon_nullity_a{spin}_q{charge}',float(norm.max())<3e-7,
              normalized_null_residual_max=float(norm.max()))

def critical_samples(min_offset=1e-8):
    # Exact R=R'=0 spherical-null solution, derived independently from the
    # radial potential. Screen offsets probe both sides of its critical curve.
    a,q=.9,.4;inc=math.radians(83)
    r=np.linspace(1.18,4.3,160);D=r*r-2*r+a*a+q*q
    lam=(r*r+a*a-2*r*D/(r-1))/a
    eta=4*r*r*D/(r-1)**2-(lam-a)**2
    beta2=eta+a*a*math.cos(inc)**2-lam*lam/math.tan(inc)**2
    use=beta2>0
    al=-lam[use]/math.sin(inc);be=np.sqrt(beta2[use])
    e=np.r_[-np.logspace(-2,math.log10(min_offset),16),np.logspace(-2,math.log10(min_offset),16)]
    return (al[:,None]*(1+e)).ravel(),(be[:,None]*(1+e)).ravel()

def higher_orders_and_storage():
    x,y=critical_samples()
    tracer=KerrNewmanAccretionTracer(8,8,angle_step=.006,max_hits=16)
    m=tracer.trace_samples(x,y,83)
    check('spherical_critical_samples_naturally_produce_multiple_orders',int(m['hit_layer'].max())>=5 and
          int(m['mask'].sum(axis=-1).max())>=6,
          maximum_image_order=int(m['hit_layer'].max()),maximum_events=int(m['mask'].sum(axis=-1).max()),
          order_histogram={str(k):int(v) for k,v in zip(*np.unique(m['hit_layer'][m['mask']],return_counts=True))},
          samples=int(x.size))
    more=KerrNewmanAccretionTracer(8,8,angle_step=.006,max_hits=24).trace_samples(x,y,83)
    keys=('r','phi','delay_tg','kr_over_E','hit_layer','hit_kind','mask')
    check('critical_sample_storage_capacity_16_to_24_converged',
          all(np.array_equal(m[k],more[k][...,:16]) for k in keys) and not np.any(more['mask'][...,16:]) and
          np.array_equal(m['status'],more['status']))
    small=KerrNewmanAccretionTracer(8,8,angle_step=.006,max_hits=1)
    omitted=small.trace_samples(x,y,83,allow_unresolved=True)
    check('event_storage_overflow_is_explicit',np.any(omitted['overflow']) and
          np.array_equal(omitted['overflow']>0,omitted['status']==5),overflow_count=int(omitted['overflow'].sum()))
    raised=False
    try:small.trace_samples(x,y,83)
    except RuntimeError:raised=True
    check('event_storage_overflow_raises_by_default',raised)

    # Finite angular step should be checked away from mathematically divergent
    # exact critical rays; there is no finite-resolution pointwise error bound
    # at the singular limiting ray itself.
    x,y=critical_samples(min_offset=1e-4)
    coarse=tracer.trace_samples(x,y,83)
    fine=KerrNewmanAccretionTracer(8,8,angle_step=.003,max_hits=16).trace_samples(x,y,83)
    check('near_critical_step_refinement_event_topology',np.array_equal(coarse['hit_kind'],fine['hit_kind']) and
          np.array_equal(coarse['hit_layer'],fine['hit_layer']))
    common=coarse['mask']&fine['mask']
    errors={key:float((angle_error(coarse[key][common],fine[key][common]) if key=='phi' else
                       abs(coarse[key][common]-fine[key][common])).max()) for key in ('r','phi','delay_tg')}
    check('near_critical_step_refinement_positions_and_delays',errors['r']<.003 and errors['phi']<.002 and
          errors['delay_tg']<.01,max_errors=errors)


def symmetry_and_quadratures():
    tracer=KerrNewmanAccretionTracer(240,200,fov=14,angle_step=.006,max_hits=12)
    a=tracer.trace(83);b=tracer.trace(-97);reused=antipodal_maps(a)
    shared=b['mask']&reused['mask']
    check('layered_antipodal_isometry_preserves_events_and_momenta',
          np.array_equal(b['hit_kind'],reused['hit_kind']) and np.array_equal(b['hit_layer'],reused['hit_layer']) and
          np.max(abs(b['r']-reused['r']))<1e-5 and np.max(abs(b['kr_over_E']-reused['kr_over_E']))<1e-4 and
          np.max(angle_error(b['phi'][shared],reused['phi'][shared]))<1e-5 and
          np.max(abs(b['lambda']-reused['lambda']))<1e-12 and np.max(abs(b['eta']-reused['eta']))<1e-12)
    candidates=np.argwhere(a['mask']&(abs(tracer.x)[...,None]>1.2)&(abs(tracer.y)[...,None]>.15))
    chosen=[]
    for kind in (1,2):
        for order in (0,1,2):
            rows=[tuple(ix) for ix in candidates if a['hit_kind'][tuple(ix)]==kind and a['hit_layer'][tuple(ix)]==order]
            if rows:chosen.append(rows[len(rows)//2])
    diagnostics=[]
    for ind in chosen:
        pixel=ind[:-1];al=float(tracer.x[pixel]);be=float(tracer.y[pixel]);r=float(a['r'][ind]);order=int(a['hit_layer'][ind])
        tau,pa,ta=angular_integrals(al,be,math.radians(83),order)
        tr,pr,qr,outgoing=radial_quadratures(al,be,math.radians(83),r,tau)
        predicted_delay=qr+ta-float(tortoise_radius(r,.9,.4))
        derr=abs(predicted_delay-float(a['delay_tg'][ind]));perr=float(angle_error(np.array(pr+pa),a['phi'][ind]))
        record=dict(pixel=[int(i) for i in pixel],kind=int(a['hit_kind'][ind]),image_order=order,outgoing=bool(outgoing),
                    delay_error_tg=derr,azimuth_error_rad=perr,mino_arrival_error=abs(tau-tr))
        diagnostics.append(record)
        check(f'independent_separated_quadrature_kind{record["kind"]}_order{order}',
              derr<.001 and perr<.0001 and abs(tau-tr)<3e-6,**record)
    check('quadratures_include_plunge_and_multiple_images',any(d['kind']==2 for d in diagnostics) and
          any(d['image_order']>=1 for d in diagnostics))
    RESULT['separated_quadrature_samples']=diagnostics


def horizon_cutoff_and_samples():
    base=KerrNewmanAccretionTracer(256,212,fov=14,angle_step=.006,max_hits=12,horizon_offset_rg=1e-4)
    check('capture_cutoff_below_innermost_exterior_photon_turn',
          base.rmin<base.minimum_scattering_turn_rg and abs(base.minimum_scattering_turn_rg-1.2147958843207893)<1e-10,
          innermost_scattering_turn_rg=base.minimum_scattering_turn_rg,capture_cutoff_rg=base.rmin)
    rejected=False
    try:KerrNewmanAccretionTracer(8,8,horizon_offset_rg=.1)
    except ValueError:rejected=True
    check('unsafe_capture_cutoff_above_photon_turn_is_rejected',rejected)
    a=base.trace(83)
    b=KerrNewmanAccretionTracer(256,212,fov=14,angle_step=.006,max_hits=12,horizon_offset_rg=1e-5).trace(83)
    shared=a['mask']&b['mask']
    check('tenfold_horizon_cutoff_refinement_capture_and_common_events',np.array_equal(a['status'],b['status']) and
          np.all(a['mask']<=b['mask']) and np.max(abs(a['r'][shared]-b['r'][shared]))<1e-6 and
          np.max(abs(a['delay_tg'][shared]-b['delay_tg'][shared]))<1e-5,
          additional_near_horizon_events=int(np.sum(b['mask']&~a['mask'])))
    ix=(np.array([3,31,76]),np.array([17,87,193]))
    sparse=base.trace_samples(base.x[ix],base.y[ix],83)
    check('independent_sparse_pixel_sampling_matches_grid',all(np.array_equal(sparse[key],a[key][ix]) for key in a))
    check('all_diagnostic_rays_finite_resolved',np.all(a['status']<3) and
          all(np.isfinite(a[k]).all() for k in ('r','phi','delay_tg','kr_over_E')) and
          float(a['constraint_residual'].max())<2e-7,
          max_constraint_residual=float(a['constraint_residual'].max()))


def main():
    started=time.perf_counter()
    for fn in (compatibility,capture_and_null_momentum,higher_orders_and_storage,symmetry_and_quadratures,horizon_cutoff_and_samples):
        try:fn()
        except Exception as error:check(fn.__name__+'_completed',False,error=repr(error))
    RESULT['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ('kn_accretion_trace.py','kn_accretion_trace.cpp','kn_orbit_trace.py','kn_orbit_trace.cpp','kn_orbits.py','verify_kn_trace.py','verify_kn_accretion_trace.py')}
    RESULT['elapsed_seconds']=time.perf_counter()-started
    RESULT['passed']=all(record['passed'] for record in RESULT['checks'])
    (ROOT/'kn-accretion-trace-diagnostics.json').write_text(json.dumps(RESULT,indent=2)+'\n')
    print('PASSED' if RESULT['passed'] else 'FAILED',len(RESULT['checks']),flush=True)
    return 0 if RESULT['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
