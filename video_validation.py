"""Presentation geometry and export integrity checks for the current Reel."""
from pathlib import Path
from fractions import Fraction
import hashlib
import math
import shutil
import subprocess
import numpy as np
from PIL import Image
from render_kn_orbit import (KerrNewmanRenderer, ROOT, VIDEO_NAMES, LOOP_VIDEO_NAMES, SPIN, CHARGE, FOV, display_intensity, illustrative_palette)
from render import encode_srgb
from video_io import check, probe_video, mp4_atoms, decode_statistics


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
