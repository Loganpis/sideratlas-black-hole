#!/usr/bin/env python3
"""Verify the current Kerr renderer, saved numerical diagnostics and final files.

These checks measure numerical consistency of the specified idealized model,
not physical validation of a turbulent flow or agreement with NASA's movie.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from PIL import Image
from render_kerr_orbit import (KerrOrbitRenderer,ROOT,OUTPUT_DIR,VIDEO_NAMES,source_hash,
    SPIN,RIN,FOV,OUTER_RADIUS,INCLINATION,LOOP_TG,view_angle_at_phase)
from verify import check,verify_video
from verify_nasa import display_checks
from instagram_export import metadata as export_metadata


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finish(record):
    record['passed']=all(c['passed'] for c in record.get('checks',[]))
    return record


def actual_renderer_checks():
    renderer=KerrOrbitRenderer(width=270,supersample=2)
    result={'scope':'Fresh checks with current rendering code at 270x480 and 2x sampling',
            'geometry':renderer.geometry,'emission':renderer.source.metadata()}
    check(result,'spin_and_inner_boundary',renderer.spin==.9 and abs(renderer.source.rin-RIN)<1e-12 and abs(RIN-2.320883041761887)<1e-12)
    maps={k:v.copy() for k,v in renderer.maps.items()}
    initial=renderer.prepared.evaluate(0)
    end=renderer.prepared.evaluate(LOOP_TG)
    check(result,'source_periodic_at_unwrapped_retarded_times',np.allclose(initial,end,atol=1e-12,rtol=0),maximum_error=float(np.max(abs(initial-end))))
    nondelayed=renderer.source.prepare(maps['r'],maps['phi'],0).evaluate(0)
    check(result,'light_arrival_delays_affect_evolving_source',not np.allclose(initial,nondelayed),rms_difference=float(np.sqrt(np.mean((initial-nondelayed)**2))))
    frames0=renderer.frames(0)
    frames1=renderer.frames(1)
    for mode in VIDEO_NAMES:
        check(result,mode+'_uncompressed_loop_endpoints_match',np.array_equal(np.asarray(frames0[mode]),np.asarray(frames1[mode])))
    phases=[0,.125,.25,.375,.5,.625,.75,.875,INCLINATION/360,(INCLINATION+180)/360]
    for phase in phases:
        radiance=renderer.radiance(phase)
        m=renderer.maps; hit=m['mask']
        check(result,f'phase_{phase:g}_finite_positive_physics',
              all(np.isfinite(m[k]).all() for k in ('r','phi','g','delay_tg')) and
              np.isfinite(radiance).all() and np.all(radiance>=0) and hit.any() and
              np.all(m['r'][hit]>=RIN-2e-6) and np.all(m['r'][hit]<=OUTER_RADIUS+2e-6) and
              np.all(m['g'][hit]>0) and np.all(m['delay_tg'][~hit]==0),
              angle_deg=view_angle_at_phase(phase),hit_count=int(hit.sum()),
              largest_constraint_residual=float(m['constraint_residual'].max()))
    smooth=[]
    for pole,phase in [('north',INCLINATION/360),('south',(INCLINATION+180)/360)]:
        smooth.append(renderer.radiance(phase,structured=False))
        m=renderer.maps;hit=m['mask'];r=m['r'][hit].astype(float)
        # Independent circular-orbit time dilation at lambda=0.
        expected=np.sqrt(r**3-3*r**2+2*SPIN*r**1.5)/(r**1.5+SPIN)
        error=np.max(abs(m['g'][hit]-expected))
        check(result,pole+'_pole_shift_matches_circular_orbit_time_dilation',error<2e-6,max_absolute_error=float(error))
        h,w=hit.shape; sl=slice((w-h)//2,(w+h)//2)
        check(result,pole+'_axisymmetric_geometry',np.array_equal(hit[:,sl],np.rot90(hit[:,sl])) and all(np.allclose(m[k][:,sl],np.rot90(m[k][:,sl]),rtol=1e-5,atol=1e-5) for k in ('r','g','delay_tg')))
    check(result,'identical_disk_faces_have_same_smooth_pole_image',np.allclose(smooth[0],smooth[1],rtol=1e-5,atol=1e-5))
    smooth0=renderer.radiance(0,structured=False)
    left,right=float(smooth0[:,:135].sum()),float(smooth0[:,135:].sum())
    check(result,'approaching_side_brighter',left>right,left_right_ratio=left/right)
    opposite=renderer.radiance(.5,structured=False)
    check(result,'opposite_view_reverses_side_asymmetry',np.allclose(smooth0,np.fliplr(opposite),atol=3e-5,rtol=3e-5))
    check(result,'higher_order_disk_images_present',np.any(maps['hit_layer']>=1),orders={str(k):int(np.sum(maps['hit_layer']==k)) for k in np.unique(maps['hit_layer'][maps['mask']])})
    check(result,'spin_breaks_side_view_left_right_geometric_symmetry',not np.allclose(maps['r'],np.fliplr(maps['r'])))
    for angle in (-90,-270):
        m=renderer.tracer.trace(angle)
        check(result,f'exact_plane_{angle}_finite',all(np.isfinite(m[k]).all() for k in ('r','phi','g','delay_tg')) and m['mask'].any())
    m=renderer.tracer.trace(INCLINATION-360)
    check(result,'full_orbit_maps_return',np.array_equal(m['mask'],maps['mask']) and all(np.allclose(m[k],maps[k],rtol=0,atol=1e-5) for k in ('r','phi','g','delay_tg')))
    return finish(result)


def metadata_checks(output):
    d=json.loads((output/'render-metadata.json').read_text());result={}
    check(result,'recorded_source_matches_current_code',d['source_sha256']==source_hash())
    check(result,'delivery_dimensions_and_timing',(d['width'],d['height'],d['fps'],d['frames'],d['duration_seconds'],d['supersample'])==(1080,1920,30,240,8,2))
    g=d['geometry'];spin=d['spin_configuration']
    check(result,'kerr_spin_disk_and_detector_parameters',spin['a_star']==.9 and spin['alignment']=='aligned, prograde equatorial disk' and g['spin']==.9 and abs(g['rin_rg']-RIN)<1e-12 and g['rout_rg']==28 and g['fov_horizontal_rg']==72 and (g['width'],g['height'])==(2160,1800))
    check(result,'camera_completes_continuous_revolution',d['view_angle_degrees']==[-277,83] and d['camera_orbit']['angular_span_deg']==-360)
    check(result,'static_observer_and_disk_face_assumptions_disclosed','static observers' in d['camera_orbit']['observer'] and 'both faces' in d['disk_faces'])
    check(result,'illustrative_color_disclosed','not thermal true color' in d['display']['color'])
    check(result,'chosen_spin_disclosed','not a fit' in spin['choice'])
    check(result,'source_parameters_match',d['emission']==__import__('kerr_orbit_emission').ShearingDiskEmission().metadata())
    check(result,'export_settings_match',d['encoding']==export_metadata())
    check(result,'brand_and_ui_clearance',d['layout']['brand']=='@sideratlas' and d['layout']['header_top']>=250 and d['layout']['last_credit_top']<1600)
    check(result,'files_fully_rendered',d.get('elapsed_seconds',0)>0 and not list(output.glob('*.partial.mp4')))
    return finish(result)


def output_assets(output):
    result={'files':{}}
    names=list(VIDEO_NAMES.values())+['render-metadata.json','sideratlas-kerr-black-hole-color-cover.png','sideratlas-kerr-black-hole-mono-cover.png','steady-reference.png','contact-sheet.jpg']
    for name in names:
        path=output/name
        check(result,name+'_exists',path.is_file() and path.stat().st_size>0)
        if path.is_file():
            result['files'][name]={'sha256':sha256(path),'size_bytes':path.stat().st_size}
            if path.suffix in ('.png','.jpg'):
                with Image.open(path) as im:
                    im.load()
                    check(result,name+'_image_decodes',im.size==((1080,1280) if name=='contact-sheet.jpg' else (1080,1920)))
    return finish(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR)
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args();args.output.mkdir(exist_ok=True,parents=True)
    start=time.perf_counter()
    result={'scope':__doc__,'preflight_only':args.preflight,'verifier_sha256':sha256(__file__),'videos':{}}
    functions=[('display',display_checks),('actual_renderer',actual_renderer_checks)]
    if not args.preflight:
        functions += [('metadata',lambda:metadata_checks(args.output)),('assets',lambda:output_assets(args.output))]
    for name,fn in functions:
        print('Checking '+name,flush=True)
        try:result[name]=fn()
        except Exception as error:result[name]={'passed':False,'error':str(error)}
    if not args.preflight:
        for mode,name in VIDEO_NAMES.items():
            print('Decoding '+name,flush=True)
            try:result['videos'][mode]=verify_video(args.output/name,shutil.which('ffmpeg'),shutil.which('ffprobe'),hero_region=(425,1145),expected_transfer='bt709',expected_range='tv')
            except Exception as error:result['videos'][mode]={'passed':False,'error':str(error)}
    result['passed']=all(result[n]['passed'] for n,_ in functions) and all(r['passed'] for r in result['videos'].values())
    result['elapsed_seconds']=time.perf_counter()-start
    destination=args.output/('verification-preflight.json' if args.preflight else 'verification.json')
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(('PASS: ' if result['passed'] else 'FAIL: ')+str(destination),flush=True)
    if not result['passed']:
        for n,r in list(result.items())+list(result['videos'].items()):
            if isinstance(r,dict):
                if 'error' in r:print(n,r['error'])
                for c in r.get('checks',[]):
                    if not c['passed']:print(n,c)
    return 0 if result['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
