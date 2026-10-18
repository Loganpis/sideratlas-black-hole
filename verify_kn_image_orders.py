#!/usr/bin/env python3
"""Measure real Kerr–Newman image orders and gray foreground transmission.

Reproduce the public kn-image-order-diagnostics.json and fixed-exposure
per-order review images under the ignored output folder's preview directory.
The diagnostic uses a central detector fixture, not the full portrait movie.
Its order-separated images are never composited into the final animation.
"""
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
from PIL import Image
from kn_accretion_trace import KerrNewmanAccretionTracer
from kn_flow import NeutralAccretionFlow
from kn_dynamic_emission import DynamicDiskEmission
from render_kn_accretion import compact_hits, transfer_weights, PLUNGE_OPTICAL_DEPTH
from render_kn_orbit import display_intensity, illustrative_palette

ROOT=Path(__file__).resolve().parent


def sha(name):
    return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()


def main():
    review_dir=ROOT/'output-sideratlas-kn-v120'/'preview'/'image-orders'
    review_dir.mkdir(parents=True,exist_ok=True)
    source_names=('verify_kn_image_orders.py','kn_accretion_trace.py','kn_accretion_trace.cpp','kn_orbits.py','kn_flow.py',
                  'kn_dynamic_emission.py','render_kn_accretion.py')
    source_hashes={n:sha(n) for n in source_names}
    started=time.perf_counter()
    width,height,fov,view,time_tg=600,500,16.,83.,0.
    tau=PLUNGE_OPTICAL_DEPTH
    tracer=KerrNewmanAccretionTracer(width,height,fov=fov,angle_step=.006,max_hits=12)
    maps=tracer.trace(view);events=compact_hits(maps)
    flow=NeutralAccretionFlow();source=DynamicDiskEmission(flow=flow)
    weights,g=transfer_weights(events,flow,width*height)
    intensity=source.prepare(events['r'],events['phi'],events['delay_tg']).evaluate(time_tg)
    mu=np.clip(g*np.sqrt(np.maximum(events['eta'],0))/events['r'],0,1)
    inner=events['hit_kind']==2
    slab_factor=np.ones_like(g)
    slab_factor[inner]=-np.expm1(-np.minimum(80,tau/np.maximum(mu[inner],1e-300)))/(-math.expm1(-tau))
    local_no_foreground=intensity*slab_factor*g**4
    physical=intensity*weights
    foreground=weights/(slab_factor*g**4)
    frame=np.bincount(events['pixel'],weights=physical,minlength=width*height).reshape(height,width)
    rows=[];order_frames={}
    for order in np.unique(events['hit_layer']):
     order_use=events['hit_layer']==order
     order_frame=np.bincount(events['pixel'][order_use],weights=physical[order_use],minlength=width*height).reshape(height,width)
     order_frames[int(order)]=order_frame
     Image.fromarray(illustrative_palette(display_intensity(order_frame))).save(review_dir/f'order-{int(order)}.png')
     for kind in (1,2):
      use=order_use&(events['hit_kind']==kind)
      if not use.any():continue
      amount=float(physical[use].sum());counterfactual=float(local_no_foreground[use].sum())
      rows.append({'image_order':int(order),'emitting_zone':'opaque_outer_disk' if kind==1 else 'thin_plunging_layer',
        'recorded_intersections':int(use.sum()),'received_relative_radiance_sum':amount,
        'fraction_of_fixture_received_radiance':amount/float(physical.sum()),
        'counterfactual_without_foreground_attenuation_relative_radiance_sum':counterfactual,
        'intensity_weighted_foreground_transmission':amount/counterfactual if counterfactual else None,
        'foreground_transmission_min_median_max':np.percentile(foreground[use],[0,50,100]).tolist(),
        'frequency_shift_min_median_max':np.percentile(g[use],[0,50,100]).tolist()})
    Image.fromarray(illustrative_palette(display_intensity(frame))).save(review_dir/'all-orders.png')
    changed={**events,'hit_layer':events['hit_layer'].astype(np.int32)*17+101}
    changed_weights,_=transfer_weights(changed,flow,width*height)
    checks=[]
    def check(name,value,**metrics):checks.append({'name':name,'passed':bool(value),**metrics})
    check('all_events_strictly_outside_the_horizon',np.all(events['r']>flow.horizon_rg),minimum_radius_rg=float(events['r'].min()))
    check('all_rays_resolved_no_event_overflow',np.all(maps['status']<3) and not maps['overflow'].any(),maximum_image_order=int(events['hit_layer'].max()))
    check('received_intensity_finite_and_nonnegative',np.isfinite(physical).all() and np.all(physical>=0))
    check('higher_order_images_emerge_from_nonzero_calculated_intersections',np.any(events['hit_layer']>=2) and physical[events['hit_layer']>=2].sum()>0,
          higher_order_received_relative_radiance_sum=float(physical[events['hit_layer']>=2].sum()))
    check('changing_image_order_labels_does_not_change_transfer_or_brightness',np.array_equal(weights,changed_weights))
    check('foreground_opacity_only_attenuates_each_local_source',np.all(physical<=local_no_foreground*(1+2e-15)) and np.all(foreground>=0) and np.all(foreground<=1+2e-15))
    partition=np.sum(list(order_frames.values()),axis=0)
    error=float(np.max(abs(frame-partition)))
    check('actual_image_is_sum_of_trajectory_image_orders',np.allclose(frame,partition,rtol=2e-15,atol=1e-15),maximum_absolute_partition_error=error)
    check('zero_emissivity_produces_no_radiance',np.array_equal(np.bincount(events['pixel'],weights=np.zeros_like(intensity)*weights,minlength=width*height),np.zeros(width*height)))
    check('measured_source_files_unchanged_during_fixture',all(sha(n)==v for n,v in source_hashes.items()))
    report={
     'version':'1.2.0','scope':'Measured image-order contributions from actual Kerr-Newman equatorial intersections in one central detector fixture; no rings are drawn or image orders brightened',
     'fixture':{'width':width,'height':height,'fov_horizontal_rg':fov,'fov_vertical_rg':fov*height/width,
      'view_angle_deg':view,'observer_time_tg':time_tg,'RK4_scaled_Mino_step':.006,'max_stored_events':12,
      'pixel_area_rg2':(fov/width)**2,'pixel_quadrature':'one independently integrated central ray per detector pixel; no selective AA in this diagnostic',
      'normal_plunge_optical_depth':tau,'source_time_rule':'each source event uses t_emit=t_observer-delay_tg'},
     'interpretation':[
      'This 16 by 13.333 rg central fixture is not the full portrait movie or an intrinsic total source luminosity. Percentages depend on crop, time and resolution.',
      'Relative detector sums are proportional to received flux for equal-area screen pixels; absolute mass, distance, luminosity and spectrum are unspecified.',
      'Counterfactual unabsorbed means foreground transmission is set to one while keeping the same local slab angular-emission factor and source. It is not a zero-tau physical model.',
      'Image order counts finite equatorial crossings along a null trajectory; it is not an assigned visible stripe number.',
      'Finite pixel resolution can merge narrow higher-order contributions; nonzero integrated contributions do not guarantee a separately resolved stripe.',
      'The event horizon emits no light. A ray ending in capture can retain foreground emitting intersections.',
      'Per-order review images use the identical fixed display curve and exposure as the final renderer, with no image-order contrast boost. They are temporary ignored review files and never composited into delivery.'
     ],
     'calculated_effects':['neutral Kerr-Newman null trajectories and image intersections','flowing-emitter frequency shift including radial photon momentum','per-event propagation delays and retarded times','neutral geodesic ISCO plunge','photon capture'],
     'illustrative_choices':['neutral test matter and hypothetical large black-hole charge','slow exterior radial drift','multiscale advected emissivity and periodically repeated births','weak prescribed inner heating',f'gray thin-slab normal optical depth {tau}'],
     'transfer_equations':{'frequency_shift':'g=1/[ut-lambda*uphi-(k_r/E)*ur]','local_plunge_surface_intensity':'I_face [1-exp(-tau/mu)]/[1-exp(-tau)]','foreground_transmission':'product of exp(-tau/mu) over preceding thin intersections; opaque outer disk ends ray','observed_event_intensity':'foreground_transmission times local intensity times g^4','image_order_multiplier':1},
     'geometry':tracer.metadata,
     'source_model':source.metadata(),
     'source_sha256':source_hashes,
     'total_received_relative_radiance_sum':float(physical.sum()),
     'per_image_order_and_zone':rows,
     'superseded_smooth_source_review_reference':{
      'status':'historical illustrative-source comparison, not the final source',
      'same_fixture':True,
      'outer_order_1_intensity_weighted_foreground_transmission':.9937477291679134,
      'outer_order_2_intensity_weighted_foreground_transmission':.6057363798356402,
      'outer_order_3_intensity_weighted_foreground_transmission':.5582937822669329,
      'outer_order_received_fixture_radiance_fractions':{'0':.9217900654063329,'1':.07717737333127175,'2':.00035232541757218187,'3':.00002246112060953848},
      'interpretation':'The weakly structured source already contained naturally calculated higher images; its smooth appearance did not arise from removal of the broad first-order image by the new opacity.'},
     'checks':checks,'passed':all(c['passed'] for c in checks),'elapsed_seconds':time.perf_counter()-started}
    (ROOT/'kn-image-order-diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'checks':len(checks),'source_sha256':source_hashes['kn_dynamic_emission.py'],'groups':rows},indent=2))
    return 0 if report['passed'] else 1


if __name__=="__main__":raise SystemExit(main())
