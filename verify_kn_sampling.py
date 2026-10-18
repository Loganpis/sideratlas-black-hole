#!/usr/bin/env python3
"""Independent critical-curve coverage and pixel-quadrature checks for v1.2.

The critical curve is used only to allocate trajectories. Its spherical
null-orbit conditions R=R'=0 are checked directly, and a zero-emissivity
source must remain dark through the complete sampled radiance pipeline.
Actual AA2/AA3/AA4 renders and a separately traced AA12 critical-pixel
reference quantify finite-resolution differences. These checks do not
claim that every infinitely thin photon subring is individually resolved.
"""
from __future__ import annotations
import gc
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from kn_accretion_trace import KerrNewmanAccretionTracer
from kn_orbits import horizon_radius
from render_kn_accretion import (KerrNewmanAccretionRenderer, critical_sampling_mask,
                                 trace_resolved, PreparedRayEmission)

ROOT=Path(__file__).resolve().parent
RESULT={'scope':__doc__,'checks':[]}


def check(name,passed,**metrics):
    RESULT['checks'].append(dict(name=name,passed=bool(passed),**metrics))
    print(('PASS ' if passed else 'FAIL ')+name,flush=True)


def independent_critical_coordinates(view,spin=.9,charge=.4):
    """Spherical radial potential and infinity detector relation, no driver calls."""
    theta=math.radians(view);si,ci=math.sin(theta),math.cos(theta)
    if abs(spin)<1e-12:
        rr=(3+math.sqrt(9-8*charge*charge))/2
        radius=rr/math.sqrt(1-2/rr+charge*charge/(rr*rr))
        angle=np.linspace(0,2*math.pi,2001)
        alpha,beta=radius*np.cos(angle),radius*np.sin(angle)
        r=np.full(alpha.shape,rr)
        lam=-alpha*si;eta=beta*beta+alpha*alpha*ci*ci
    elif abs(si)<1e-12:
        # lambda=0 implies this cubic directly from the spherical condition.
        roots=np.roots([1,-3,spin*spin+2*charge*charge,spin*spin])
        radius_orbit=max(v.real for v in roots if abs(v.imag)<1e-10 and v.real>horizon_radius(spin,charge))
        dd=radius_orbit*radius_orbit-2*radius_orbit+spin*spin+charge*charge
        K=4*radius_orbit*radius_orbit*dd/(radius_orbit-1)**2
        angle=np.linspace(0,2*math.pi,2001)
        alpha,beta=math.sqrt(K)*np.cos(angle),math.sqrt(K)*np.sin(angle)
        r=np.full(alpha.shape,radius_orbit)
        lam=np.zeros(alpha.shape);eta=np.full(alpha.shape,K-spin*spin)
    else:
        rr=np.linspace(horizon_radius(spin,charge)+1e-7,6,50000)
        dd=rr*rr-2*rr+spin*spin+charge*charge
        P=2*rr*dd/(rr-1)
        ll=(rr*rr+spin*spin-P)/spin
        ee=4*rr*rr*dd/(rr-1)**2-(ll-spin)**2
        b2=ee+spin*spin*ci*ci-ll*ll*ci*ci/(si*si)
        use=b2>0
        alpha=np.tile(-ll[use]/si,2)
        beta=np.r_[np.sqrt(b2[use]),-np.sqrt(b2[use])]
        r=np.tile(rr[use],2);lam=np.tile(ll[use],2);eta=np.tile(ee[use],2)
    return alpha,beta,r,lam,eta


def critical_geometry_and_coverage():
    records=[]
    for spin,charge,view in ((.9,.4,83),(.9,.4,38),(.9,.4,0),(.9,.4,-97),
                             (0,.4,83),(0,.4,0),(0,0,38)):
        alpha,beta,r,lam,eta=independent_critical_coordinates(view,spin,charge)
        delta=r*r-2*r+spin*spin+charge*charge
        P=r*r+spin*spin-spin*lam;K=eta+(lam-spin)**2
        radial=P*P-delta*K
        derivative=4*r*P-(2*r-2)*K
        rerr=float(np.max(abs(radial)/(1+P*P+abs(delta*K))))
        derr=float(np.max(abs(derivative)/(1+abs(4*r*P)+abs((2*r-2)*K))))
        check(f'independent_spherical_null_conditions_a{spin}_q{charge}_view{view}',
              rerr<2e-13 and derr<2e-13,normalized_R_max=rerr,normalized_Rprime_max=derr)
        width,height,fov=405,338,72.
        mask=critical_sampling_mask(width,height,fov,view,spin,charge)
        # Detector pixels are boundary-centered exactly as traced. Round-off
        # at a critical point's integer boundary is covered by the wider
        # conservative sampling band, without changing any photon path.
        x=np.floor(width/2+alpha*width/fov).astype(int)
        y=np.floor(height/2-beta*width/fov).astype(int)
        valid=(x>=0)&(x<width)&(y>=0)&(y<height)
        fraction=float(mask[y[valid],x[valid]].mean())
        check(f'critical_curve_pixel_allocation_coverage_a{spin}_q{charge}_view{view}',
              fraction==1 and mask.dtype==bool and 0<mask.mean()<.1,
              coverage_fraction=fraction,sampling_pixel_fraction=float(mask.mean()),samples=int(valid.sum()))
        records.append(dict(spin=spin,charge=charge,view=view,coverage=fraction,
                            normalized_R_max=rerr,normalized_Rprime_max=derr))
    RESULT['critical_curve_samples']=records


def preview_quadrature_convergence():
    renderer=KerrNewmanAccretionRenderer(270,1.5,refinement=2,temporal_samples=1)
    records=[]
    for phase,label in ((0.,'side83'),(83/360,'pole0'),(.125,'view38')):
        samples={};timing={}
        for n in (2,3,4):
            renderer.refinement=n
            renderer.current_inclination=None
            started=time.perf_counter()
            values=renderer.radiance(phase)
            shown=renderer.viewport_radiance(values)
            rgb=np.asarray(renderer.compose(values,'color'))
            timing[str(n)]=time.perf_counter()-started
            samples[n]=(values,shown,rgb)
        reference=samples[4]
        metrics={}
        for n in (2,3):
            values,shown,rgb=samples[n]
            metrics[str(n)]={
                'native_total_flux_relative_difference':float(abs(values.astype(float).sum()-reference[0].astype(float).sum())/reference[0].astype(float).sum()),
                'native_L1_relative_to_total_flux':float(abs(values.astype(float)-reference[0]).sum()/reference[0].astype(float).sum()),
                'viewport_total_flux_relative_difference':float(abs(shown.astype(float).sum()-reference[1].astype(float).sum())/reference[1].astype(float).sum()),
                'viewport_L1_relative_to_total_flux':float(abs(shown.astype(float)-reference[1]).sum()/reference[1].astype(float).sum()),
                'RGB_mean_absolute_difference_0_to_255':float(np.mean(abs(rgb.astype(float)-reference[2]))),
                'RGB_max_absolute_difference_0_to_255':int(np.max(abs(rgb.astype(int)-reference[2].astype(int))))}
        a,b=metrics['2'],metrics['3']
        check(f'AA2_to_AA4_preview_global_flux_and_RGB_{label}',
              a['native_total_flux_relative_difference']<.002 and a['native_L1_relative_to_total_flux']<.005 and
              a['viewport_L1_relative_to_total_flux']<.005 and a['RGB_mean_absolute_difference_0_to_255']<.05,
              **a)
        check(f'AA3_to_AA4_reduces_preview_integrated_absolute_error_{label}',
              b['native_L1_relative_to_total_flux']<a['native_L1_relative_to_total_flux'] and
              b['viewport_L1_relative_to_total_flux']<a['viewport_L1_relative_to_total_flux'],**b)
        check(f'positive_finite_real_trajectory_pixel_quadrature_{label}',
              all(np.isfinite(v[0]).all() and np.all(v[0]>=0) for v in samples.values()),
              native_detector_shape=list(reference[0].shape),refined_parent_pixels=int(len(renderer.refined_pixels)))
        records.append(dict(view=83-360*phase,phase=phase,relative_to_AA4=metrics,seconds=timing))
        del samples
        gc.collect()
    RESULT['preview_pixel_quadrature']={
        'final_frame_size':[270,480],'native_supersample':1.5,'native_detector_shape':[338,405],
        'temporal_samples':1,'comparison':'same source and geodesic settings; only independent within-pixel trajectory quadrature changes',
        'views':records,
        'limitation':'preview evidence for integrated image stability; does not imply individual convergence of arbitrarily narrow subrings'}
    return renderer


def critical_pixel_reference(renderer):
    renderer.refinement=2;renderer.current_inclination=None;renderer.set_view(0)
    mask=critical_sampling_mask(renderer.tracer.width,renderer.tracer.height,renderer.fov,83)
    selected=np.flatnonzero(mask)
    selected=selected[np.unique(np.linspace(0,len(selected)-1,96).astype(int))]
    # This reference uses denser independently integrated trajectories than
    # the driver, h=.003 and 12x12 quadrature, at the same detector footprint.
    tracer=KerrNewmanAccretionTracer(8,8,angle_step=.003,max_hits=16)
    values={};orders={};counts={}
    for n in (2,4,8,12):
        grid=(np.arange(n)+.5)/n-.5
        dx,dy=np.meshgrid(grid,grid)
        spacing=renderer.fov/renderer.tracer.width
        alpha=renderer.tracer.x.ravel()[selected,None]+dx.ravel()*spacing
        beta=renderer.tracer.y.ravel()[selected,None]+dy.ravel()*spacing
        events,diag=trace_resolved(tracer,alpha,beta,83)
        ray=PreparedRayEmission(events,renderer.source,renderer.flow,alpha.size)
        values[n]=ray.evaluate(0).reshape(-1,n*n).mean(axis=1)
        orders[str(n)]=int(diag['maximum_image_order'])
        counts[str(n)]={'rays':int(alpha.size),'plunging_events':int(diag['plunging_events'])}
    reference=values[12];total=float(reference.sum())
    metrics={str(n):{
        'total_flux_relative_difference':float(abs(values[n].sum()-total)/total),
        'L1_relative_to_patch_flux':float(abs(values[n]-reference).sum()/total),
        'maximum_relative_pixel_difference':float(np.max(abs(values[n]-reference)/(1e-6+reference)))} for n in (2,4,8)}
    check('independent_AA12_critical_pixel_reference_contains_real_plunge_and_higher_images',
          orders['12']>=3 and counts['12']['plunging_events']>0,
          maximum_image_orders=orders,sampling=counts)
    check('critical_pixel_AA4_improves_integrated_error_over_AA2',
          metrics['4']['L1_relative_to_patch_flux']<metrics['2']['L1_relative_to_patch_flux'],
          relative_to_AA12=metrics)
    check('critical_pixel_AA8_to_AA12_integrated_L1_below_half_percent',
          metrics['8']['L1_relative_to_patch_flux']<.005,
          **metrics['8'])
    RESULT['critical_pixel_reference']={
        'view_deg':83,'detector_pixel_indices':selected.tolist(),'independent_RK4_step':.003,
        'reference_subrays_per_pixel':144,'relative_to_AA12':metrics,
        'maximum_image_orders':orders,'sampling':counts,
        'limitation':'isolated critical pixels can retain large relative differences despite very small integrated image error; no order is artificially broadened or brightened'}


def zero_source_remains_dark(renderer):
    class ZeroSource:
        def __init__(self,original):self.original=original
        def evaluate(self,time_tg,structured=True):return np.zeros_like(self.original.evaluate(time_tg,structured=structured))
    saved=[]
    for ray in (renderer.coarse,renderer.fine):
        if ray is not None:
            saved.append((ray,ray.prepared))
            ray.prepared=ZeroSource(ray.prepared)
    try:
        image=renderer._radiance_at_time(73.15)
        check('critical_curve_sampling_adds_no_light_for_zero_emissivity',
              np.array_equal(image,np.zeros_like(image)),maximum_radiance=float(image.max()))
    finally:
        for ray,prepared in saved:ray.prepared=prepared


def main():
    started=time.perf_counter()
    try:
        critical_geometry_and_coverage()
        renderer=preview_quadrature_convergence()
        critical_pixel_reference(renderer)
        zero_source_remains_dark(renderer)
    except Exception as error:
        check('sampling_diagnostics_completed',False,error=repr(error))
    RESULT['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
        ('verify_kn_sampling.py','render_kn_accretion.py','render_kn_orbit.py','kn_accretion_trace.py','kn_accretion_trace.cpp',
         'kn_flow.py','kn_dynamic_emission.py','kn_orbit_emission.py','kn_orbits.py')}
    RESULT['elapsed_seconds']=time.perf_counter()-started
    RESULT['passed']=all(item['passed'] for item in RESULT['checks'])
    (ROOT/'kn-sampling-diagnostics.json').write_text(json.dumps(RESULT,indent=2)+'\n')
    print('PASSED' if RESULT['passed'] else 'FAILED',len(RESULT['checks']),flush=True)
    return 0 if RESULT['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
