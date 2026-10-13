#!/usr/bin/env python3
"""Verify the current Kerr–Newman renderer, saved numerical diagnostics and final files.

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
import subprocess
from fractions import Fraction
import numpy as np
from PIL import Image
from render_kn_orbit import (KerrNewmanRenderer,ROOT,OUTPUT_DIR,VIDEO_NAMES,source_hash,
    SPIN,CHARGE,RIN,FOV,OUTER_RADIUS,INCLINATION,LOOP_TG,view_angle_at_phase,LOOP_VIDEO_NAMES,DURATION,LOOP_DURATION,FPS,REPEATS)
from verify import check,probe_video,mp4_atoms,decode_statistics
from render_kn_orbit import display_intensity,illustrative_palette
from instagram_export import metadata as export_metadata


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def portable_path(path):
    """Keep published reports free of local machine paths."""
    try:
        return str(Path(path).resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return Path(path).name


def finish(record):
    record['passed']=all(c['passed'] for c in record.get('checks',[]))
    return record


def display_checks():
    values=np.linspace(0,1,65537)
    palette=illustrative_palette(values).astype(float)/255
    linear=np.where(palette<=.04045,palette/12.92,((palette+.055)/1.055)**2.4)
    luminance=linear@np.array([.2126,.7152,.0722])
    curve=display_intensity(np.geomspace(1e-12,1e12,10000))
    result={'scope':'Current KN renderer display functions; color is illustrative brightness'}
    check(result,'black_exact',np.array_equal(palette[0],np.zeros(3)))
    check(result,'palette_luminance_monotonic',np.all(np.diff(luminance)>=-1e-15))
    check(result,'fixed_display_curve_monotonic_bounded',np.isfinite(curve).all() and np.all(curve>=0) and np.all(curve<=1) and np.all(np.diff(curve)>=0))
    return finish(result)


def numerical_reports():
    result={'scope':'Saved numerical diagnostics with every recorded source hash checked against the current files','reports':{}}
    for name in ('kn-orbit-trace-diagnostics.json','kn-orbit-emission-diagnostics.json'):
        path=ROOT/name
        report=json.loads(path.read_text())
        passed=report.get('passed',report.get('status')=='passed')
        hashes=report['source_sha256']
        matches={k:sha256(ROOT/k)==v for k,v in hashes.items()}
        result['reports'][name]={'sha256':sha256(path),'checks':len(report['checks']),'source_matches':matches}
        check(result,name+'_passed_with_current_sources',passed and all(matches.values()) and all(c['passed'] for c in report['checks']))
    return finish(result)


def actual_renderer_checks():
    renderer=KerrNewmanRenderer(width=270,supersample=2)
    result={'scope':'Fresh checks with current rendering code at 270x480 and 2x sampling',
            'geometry':renderer.geometry,'emission':renderer.source.metadata()}
    check(result,'spin_and_inner_boundary',renderer.spin==.9 and abs(renderer.source.rin-RIN)<1e-12 and abs(RIN-1.577504566512094)<1e-12 and renderer.charge==.4)
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
              radiance.shape==m['r'].shape and np.isfinite(radiance).all() and np.all(radiance>=0) and hit.any() and
              np.all(m['r'][hit]>=RIN-2e-6) and np.all(m['r'][hit]<=OUTER_RADIUS+2e-6) and
              np.all(m['g'][hit]>0) and np.all(m['delay_tg'][~hit]==0),
              angle_deg=view_angle_at_phase(phase),hit_count=int(hit.sum()),
              largest_constraint_residual=float(m['constraint_residual'].max()))
    smooth=[]
    for pole,phase in [('north',INCLINATION/360),('south',(INCLINATION+180)/360)]:
        smooth.append(renderer.radiance(phase,structured=False))
        m=renderer.maps;hit=m['mask'];r=m['r'][hit].astype(float)
        # Independent circular-orbit time dilation at lambda=0.
        omega=np.sqrt(r-CHARGE**2)/(r*r+SPIN*np.sqrt(r-CHARGE**2))
        gtt=-(1-2/r+CHARGE**2/r**2)
        gtphi=-SPIN*(2/r-CHARGE**2/r**2)
        gphiphi=r*r+SPIN**2+SPIN**2*(2/r-CHARGE**2/r**2)
        expected=np.sqrt(-gtt-2*omega*gtphi-omega**2*gphiphi)
        error=np.max(abs(m['g'][hit]-expected))
        check(result,pole+'_pole_shift_matches_circular_orbit_time_dilation',error<2e-6,max_absolute_error=float(error))
        h,w=hit.shape; sl=slice((w-h)//2,(w+h)//2)
        check(result,pole+'_axisymmetric_geometry',np.array_equal(hit[:,sl],np.rot90(hit[:,sl])) and all(np.allclose(m[k][:,sl],np.rot90(m[k][:,sl]),rtol=1e-5,atol=1e-5) for k in ('r','g','delay_tg')))
    check(result,'identical_disk_faces_have_same_smooth_pole_image',np.allclose(smooth[0],smooth[1],rtol=1e-5,atol=1e-5))
    for phase in (0,.125,INCLINATION/360):
        reused=renderer.antipodal_radiance(phase)
        traced=renderer.radiance(phase+.5)
        error=float(np.max(abs(reused-traced)))
        check(result,f'antipodal_retarded_radiance_matches_fresh_trace_{phase:g}',np.allclose(reused,traced,rtol=5e-5,atol=2e-6),maximum_absolute_error=error)
    smooth0=renderer.radiance(0,structured=False)
    midpoint=smooth0.shape[1]//2
    left,right=float(smooth0[:,:midpoint].sum()),float(smooth0[:,midpoint:].sum())
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


def composition_checks():
    """Independent synthetic fields expose distortion, framing and footer regressions."""
    renderer=KerrNewmanRenderer.__new__(KerrNewmanRenderer)
    renderer.width,renderer.height,renderer.scale=270,480,.25
    renderer.spin,renderer.charge,renderer.fov=SPIN,CHARGE,FOV
    renderer._configure_layout()
    renderer.overlays={mode:renderer._overlay(mode) for mode in VIDEO_NAMES}
    result={'scope':'Actual compositor applied to synthetic native float-radiance fields; no ray tracing needed'}
    for mode,overlay in renderer.overlays.items():
        alpha=np.asarray(overlay)[...,3]
        box=overlay.getbbox()
        check(result,mode+'_only_header_text_remains',box is not None and
              box[1]>=round(250*renderer.scale) and not alpha[round(500*renderer.scale):].any(),
              overlay_bbox=list(box) if box else None)
    check(result,'color_and_intensity_have_same_header',np.array_equal(
        np.asarray(renderer.overlays['color']),np.asarray(renderer.overlays['mono'])))
    # Removing the overlay here isolates spatial sampling from text/shadow pixels.
    renderer.overlays={mode:Image.new('RGBA',(renderer.width,renderer.height)) for mode in VIDEO_NAMES}
    shape=(450,540)
    constant=np.full(shape,.02,np.float32)
    image=np.asarray(renderer.compose(constant,'color'))
    expected=illustrative_palette(display_intensity(np.array(.02)))
    check(result,'uniform_radiance_fills_entire_portrait',np.all(image==expected),
          maximum_channel_error=int(np.abs(image.astype(int)-expected.astype(int)).max()))
    spacing=FOV/shape[1]
    x=(np.arange(shape[1])+.5-shape[1]/2)*spacing
    y=-(np.arange(shape[0])+.5-shape[0]/2)*spacing
    xx,yy=np.meshgrid(x,y)
    px,py=np.meshgrid(np.arange(renderer.width)+.5,np.arange(renderer.height)+.5)
    def image_moments(field):
        weights=np.asarray(renderer.compose(field.astype(np.float32),'mono'))[...,0].astype(float)
        total=weights.sum()
        cx,cy=float((weights*px).sum()/total),float((weights*py).sum()/total)
        dx,dy=px-cx,py-cy
        covariance=np.array([[(weights*dx*dx).sum(),(weights*dx*dy).sum()],
                             [(weights*dx*dy).sum(),(weights*dy*dy).sum()]])/total
        return cx,cy,covariance
    cx,cy,cov=image_moments(.02*np.exp(-(xx**2+yy**2)/(2*1.2**2)))
    check(result,'centered_circular_source_is_not_stretched',
          abs(cx-renderer.width/2)<.1 and abs(cy-renderer.height/2)<.1 and
          abs(cov[0,0]/cov[1,1]-1)<.01 and abs(cov[0,1])/np.trace(cov)<.005,
          centroid_pixels=[cx,cy],variance_ratio=float(cov[0,0]/cov[1,1]))
    _,_,cov=image_moments(.02*np.exp(-yy**2/(2*.18**2))*(abs(xx)<25))
    _,vectors=np.linalg.eigh(cov)
    axis=vectors[:,-1]
    angle=math.degrees(math.atan2(axis[1],axis[0]))%180
    expected_angle=math.degrees(math.atan2(renderer.height,renderer.width))
    check(result,'edge_on_emission_spans_corner_diagonal',abs(angle-expected_angle)<1,
          measured_angle_degrees=angle,expected_angle_degrees=expected_angle)
    return finish(result)


def metadata_checks(output):
    d=json.loads((output/'render-metadata.json').read_text());result={}
    check(result,'recorded_source_matches_current_code',d['source_sha256']==source_hash())
    check(result,'delivery_dimensions_and_timing',(d['width'],d['height'],d['fps'],d['frames'],d['duration_seconds'],d['supersample'])==(1080,1920,30,960,32,2.5))
    g=d['geometry'];spin=d['spin_configuration']
    check(result,'kerr_newman_spin_charge_disk_and_detector_parameters',spin['a_star']==.9 and spin['charge_q']==.4 and g['charge_q']==.4 and spin['alignment']=='aligned, prograde equatorial disk' and g['spin']==.9 and abs(g['rin_rg']-RIN)<1e-12 and g['rout_rg']==28 and g['fov_horizontal_rg']==72 and g['roll_deg']==0 and (g['width'],g['height'])==(2700,2250))
    check(result,'camera_completes_continuous_revolution',d['view_angle_degrees']==[-277,83] and d['camera_orbit']['angular_span_deg']==-360)
    check(result,'static_observer_and_disk_face_assumptions_disclosed','static observers' in d['camera_orbit']['observer'] and 'both faces' in d['disk_faces'])
    check(result,'illustrative_color_disclosed','not thermal true color' in d['display']['color'])
    check(result,'chosen_spin_disclosed','not a fit' in spin['choice'])
    check(result,'source_parameters_match',d['emission']==__import__('kn_orbit_emission').ShearingDiskEmission().metadata())
    check(result,'export_settings_match',d['encoding']==export_metadata())
    layout=d['layout']
    check(result,'brand_header_and_caption_only_explanation',layout['brand']=='@sideratlas' and layout['header_top']>=250 and layout['lower_text_in_caption'] and 'last_credit_top' not in layout)
    check(result,'simulation_uses_entire_portrait',layout['full_frame'] and
          (layout['hero_x'],layout['hero_y'],layout['hero_width'],layout['hero_height'])==(0,0,1080,1920))
    roll=math.radians(layout['display_roll_degrees'])
    fov=layout['viewport_fov_horizontal_rg']
    half_w=fov/2;half_h=half_w*d['height']/d['width']
    extent_x=abs(math.cos(roll))*half_w+abs(math.sin(roll))*half_h
    extent_y=abs(math.sin(roll))*half_w+abs(math.cos(roll))*half_h
    check(result,'diagonal_viewport_stays_inside_traced_domain',
          fov==28 and abs(roll+math.atan2(d['height'],d['width']))<1e-12 and
          extent_x<g['fov_horizontal_rg']/2 and extent_y<g['fov_vertical_rg']/2,
          viewport_half_extent_rg=[extent_x,extent_y],
          traced_half_extent_rg=[g['fov_horizontal_rg']/2,g['fov_vertical_rg']/2])
    native_density=g['width']/g['fov_horizontal_rg']
    display_density=d['width']/fov
    check(result,'native_detector_density_near_delivery_pixel_density',native_density/display_density>=.97,
          native_samples_per_final_pixel=native_density/display_density)
    check(result,'files_fully_rendered',d.get('elapsed_seconds',0)>0 and not list(output.glob('*.partial.mp4')))
    return finish(result)


def output_assets(output):
    result={'files':{}}
    names=list(VIDEO_NAMES.values())+list(LOOP_VIDEO_NAMES.values())+['render-metadata.json','sideratlas-kerr-newman-color-cover.png','sideratlas-kerr-newman-mono-cover.png','steady-reference.png','contact-sheet.jpg']
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


def video_checks(path):
    ffmpeg,ffprobe=shutil.which('ffmpeg'),shutil.which('ffprobe')
    result={'file':portable_path(path),'sha256':sha256(path),'size_bytes':path.stat().st_size}
    probe,s=probe_video(path,ffprobe)
    result['stream']={k:s.get(k) for k in ['codec_name','profile','width','height','pix_fmt',
      'r_frame_rate','avg_frame_rate','nb_frames','duration','start_time','color_space',
      'color_transfer','color_primaries','color_range','level','bit_rate','chroma_location']}
    check(result,'1080x1920_h264_high_420',s['width']==1080 and s['height']==1920 and s['codec_name']=='h264' and s['profile']=='High' and s['pix_fmt']=='yuv420p')
    check(result,'thirty_fps_960_frames_32_seconds',Fraction(s['r_frame_rate'])==30 and Fraction(s['avg_frame_rate'])==30 and int(s['nb_frames'])==960 and abs(float(s['duration'])-32)<1e-6)
    check(result,'zero_start_silent',abs(float(s['start_time']))<1e-9 and not any(x['codec_type']=='audio' for x in probe['streams']))
    check(result,'rec709_limited_range',all(s[k]=='bt709' for k in ('color_space','color_transfer','color_primaries')) and s['color_range']=='tv')
    atoms=mp4_atoms(path);result['atoms']=atoms
    check(result,'faststart',atoms.get('moov',1e20)<atoms.get('mdat',-1))
    result['decode']=decode_statistics(path,ffmpeg,1080,1920,hero_region=(0,1920))
    check(result,'all_960_frames_decode_cleanly',result['decode']['decoded_frames']==960 and not result['decode']['decoder_errors'])
    check(result,'seam_consistent_with_adjacent_motion',result['decode']['seam_consistent_with_nearby_transitions'])
    # Check decoded YUV hashes, including the IDR at each eight-second seam.
    raw=subprocess.run([ffmpeg,'-v','error','-i',str(path),'-map','0:v:0','-f','framemd5','-'],capture_output=True,text=True,check=True).stdout
    hashes=[line.split(',')[-1].strip() for line in raw.splitlines() if line and not line.startswith('#')]
    check(result,'four_pixel_identical_decoded_cycles',len(hashes)==960 and all(hashes[j*240:(j+1)*240]==hashes[:240] for j in range(1,4)))
    result['first_cycle_framemd5']=hashes[:240]
    return finish(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR)
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args();args.output.mkdir(exist_ok=True,parents=True)
    start=time.perf_counter()
    result={'scope':__doc__,'preflight_only':args.preflight,'verifier_sha256':sha256(__file__),'videos':{}}
    functions=[('numerical_reports',numerical_reports),('display',display_checks),('composition',composition_checks),('actual_renderer',actual_renderer_checks)]
    if not args.preflight:
        functions += [('metadata',lambda:metadata_checks(args.output)),('assets',lambda:output_assets(args.output))]
    for name,fn in functions:
        print('Checking '+name,flush=True)
        try:result[name]=fn()
        except Exception as error:result[name]={'passed':False,'error':str(error)}
    if not args.preflight:
        for mode,name in VIDEO_NAMES.items():
            print('Decoding '+name,flush=True)
            try:result['videos'][mode]=video_checks(args.output/name)
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
