#!/usr/bin/env python3
"""Original Kerr–Newman Reel with a neutral prograde test disk.

Four identical eight-second loops form the 32-second delivery.

No NASA frames or source code are used. Evolving emission is an explicitly
specified kinematic model; the orange palette encodes bolometric brightness.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import inspect
import json
import math
from pathlib import Path
import time
import shutil
import subprocess

import numpy as np
from PIL import Image, ImageDraw

from kn_orbit_trace import KerrNewmanOrbitTracer
from kn_orbits import isco_radius
from kn_orbit_emission import ShearingDiskEmission
from render import Renderer as FontProvider, ROOT, FPS, encode_srgb
from instagram_export import encoder, metadata as encoding_metadata

FOV = 72.0
VIEWPORT_FOV = 28.0
DISPLAY_ROLL = -math.degrees(math.atan2(16, 9))
DEFAULT_SUPERSAMPLE = 2.5
OUTER_RADIUS = 28.0
INCLINATION = 83.0
LOOP_TG = 800.0
EXPOSURE = 0.02
SPIN = 0.9
CHARGE = 0.4
RIN = isco_radius(SPIN, CHARGE)
LOOP_DURATION = 8
REPEATS = 4
DURATION = LOOP_DURATION * REPEATS
OUTPUT_DIR = ROOT/'output-sideratlas-kn'
VIDEO_NAMES = {'color': 'sideratlas-kerr-newman-reel.mp4',
               'mono': 'sideratlas-kerr-newman-intensity.mp4'}
LOOP_VIDEO_NAMES = {mode:Path(name).stem+'-loop.mp4' for mode,name in VIDEO_NAMES.items()}


def view_angle_at_phase(phase, inclination=INCLINATION):
    """One continuous meridional revolution, with forward source time.

    This is an oriented camera angle, not a polar inclination restricted to
    [0, 180]. The tracer carries the screen basis continuously through both
    poles and both disk-plane crossings. Modulo one turn gives an exact loop.
    """
    return inclination - 360*(phase % 1)


# Retain the previous public helper name for scripts using this renderer.
inclination_at_phase = view_angle_at_phase


def display_intensity(radiance):
    """One fixed monotone brightness curve, shared by every frame."""
    v = np.maximum(np.asarray(radiance, dtype=np.float64), 0) / EXPOSURE
    return v / (1 + v)


def illustrative_palette(values):
    """Original, monotone-luminance orange display palette, not a spectrum."""
    stops = np.array([0, .06, .18, .35, .55, .76, .92, 1.])
    colors = np.array([[0,0,0], [17,0,0], [73,3,0], [151,15,0],
                       [225,48,0], [255,121,7], [255,200,46], [255,246,177]], float)
    return np.uint8(np.rint(np.stack([np.interp(values, stops, colors[:,c])
                                     for c in range(3)], axis=-1)))


class KerrNewmanRenderer:
    font = FontProvider.font

    def __init__(self, width=1080, supersample=DEFAULT_SUPERSAMPLE, inclination=INCLINATION,
                 fov=FOV, outer_radius=OUTER_RADIUS, spin=SPIN, charge=CHARGE, angle_step=.012):
        if width < 270 or width % 2 or supersample < 1:
            raise ValueError('even width >=270 and supersample >=1 required')
        self.width, self.height = width, round(width*16/9)//2*2
        self.scale = width/1080
        self.fov = float(fov)
        self._configure_layout()
        self.inclination = inclination
        self.spin, self.charge = float(spin), float(charge)
        self.source = ShearingDiskEmission(spin_a_star=self.spin, charge_q=self.charge, outer_radius_rg=outer_radius, loop_period_tg=LOOP_TG)
        self.tracer = KerrNewmanOrbitTracer(round(width*supersample), round(self.hero_h*supersample),
                                      fov=fov, roll_deg=0, spin=self.spin, charge_q=self.charge, rout=outer_radius,
                                      angle_step=angle_step)
        self.geometry = self.tracer.metadata.copy()
        self.current_inclination = None
        self.set_view(0)
        self.overlays = {mode: self._overlay(mode) for mode in VIDEO_NAMES}

    def _configure_layout(self):
        # Trace a reflection-symmetric, unrolled domain so the antipodal
        # geodesic reuse remains exact. A fixed rotated portrait viewport is
        # sampled from the native float radiance, before the display curve.
        self.hero_h = round(900*self.scale)
        self.hero_x, self.hero_y = 0, 0
        self.hero_display_size = (self.width, self.height)
        self.viewport_fov = VIEWPORT_FOV
        self.display_roll = DISPLAY_ROLL

    def set_view(self, phase):
        angle=view_angle_at_phase(phase,self.inclination)
        if self.current_inclination is not None and abs(angle-self.current_inclination)<1e-12:
            return
        self.maps=self.tracer.trace(angle)
        self.prepared=self.source.prepare(self.maps['r'],self.maps['phi'],self.maps['delay_tg'])
        self.transfer=self.maps['g'].astype(np.float64)**4
        self.current_inclination=angle

    def _overlay(self, mode):
        overlay = Image.new('RGBA', (self.width, self.height))
        d, s = ImageDraw.Draw(overlay), self.scale
        white, soft = (245,245,245), (208,208,208)
        def text(x,y,value,size=28,fill=white,role='regular',anchor='lt'):
            d.text((x*s,y*s),value,font=self.font(size,role),fill=fill,anchor=anchor,
                   stroke_width=max(1,round(2*s)),stroke_fill=(0,0,0,220))
        text(72,270,'@sideratlas',34)
        text(965,277,'NO.001',25,soft,anchor='rt')
        text(64,333,'THE BLACK HOLE',82,role='bold')
        return overlay

    def radiance(self, phase, structured=True):
        self.set_view(phase)
        t = (phase % 1)*LOOP_TG
        values = self.prepared.evaluate(t, structured=structured)*self.transfer
        if not np.isfinite(values).all() or np.any(values < 0):
            raise AssertionError('nonfinite or negative radiance')
        return values.astype('float32')

    def antipodal_radiance(self, phase):
        """Exact isometry maps the observer to phase+1/2; no image warping.

        Camera roll is zero. Neutral KN and the equatorial disk are invariant
        under a pi azimuth rotation plus reflection through the equator.
        Source time still advances, and the shifted disk azimuth is evaluated
        afresh: this does not mirror a previously rendered image.
        """
        from kn_orbit_trace import antipodal_maps
        self.set_view(phase)
        if self.tracer.roll_deg != 0:
            raise ValueError('antipodal screen mapping requires zero roll')
        maps=antipodal_maps(self.maps)
        prepared=self.source.prepare(maps['r'],maps['phi'],maps['delay_tg'])
        t=((phase+.5)%1)*LOOP_TG
        values=prepared.evaluate(t)*maps['g'].astype(np.float64)**4
        if not np.isfinite(values).all() or np.any(values<0):
            raise AssertionError('invalid antipodal radiance')
        return values.astype('float32')

    def viewport_affine(self, radiance_shape, output_size):
        """Inverse rigid image-plane rotation plus one isotropic pixel scale.

        PIL transforms pixel centers in boundary coordinates. The source
        detector has square pixels of size fov/source_width in r_g.
        """
        source_h, source_w = radiance_shape
        output_w, output_h = output_size
        pixel_scale = self.viewport_fov/output_w * source_w/self.fov
        angle = math.radians(self.display_roll)
        a, b = math.cos(angle)*pixel_scale, -math.sin(angle)*pixel_scale
        d, e = math.sin(angle)*pixel_scale, math.cos(angle)*pixel_scale
        return (a,b,source_w/2-a*output_w/2-b*output_h/2,
                d,e,source_h/2-d*output_w/2-e*output_h/2)

    def viewport_radiance(self, radiance):
        # Bilinear interpolation is nonnegative and cannot introduce bright
        # ringing around a photon image. Average a 2x portrait grid in radiance
        # before the nonlinear display transform; never rotate encoded RGB.
        size = (self.width*2,self.height*2)
        sampled = Image.fromarray(np.asarray(radiance,dtype='float32')).transform(
            size,Image.Transform.AFFINE,self.viewport_affine(radiance.shape,size),
            resample=Image.Resampling.BILINEAR,fillcolor=0)
        return np.asarray(sampled.resize((self.width,self.height),Image.Resampling.BOX))

    def compose(self, radiance, mode):
        presented = self.viewport_radiance(radiance)
        level = display_intensity(presented)
        rgb = illustrative_palette(level) if mode=='color' else np.repeat(encode_srgb(level)[...,None],3,axis=2)
        canvas = Image.fromarray(rgb).convert('RGBA')
        return Image.alpha_composite(canvas,self.overlays[mode]).convert('RGB')

    def frames(self, phase, modes=tuple(VIDEO_NAMES), structured=True):
        values = self.radiance(phase,structured=structured)
        return {mode:self.compose(values,mode) for mode in modes}

    def metadata(self):
        self.set_view(0)
        hit = self.maps['mask']
        return {'width':self.width,'height':self.height,'fps':FPS,'frames':FPS*DURATION,
                'duration_seconds':DURATION,'loop_duration_seconds':LOOP_DURATION,
                'repetitions':REPEATS,'unique_frames':FPS*LOOP_DURATION,'inclination_degrees':[0,180],
                'view_angle_degrees':[self.inclination-360,self.inclination],
                'camera_orbit':{'law':'theta = theta_start - 360° phase (modulo one turn)',
                                'start_deg':self.inclination,'angular_span_deg':-360,
                                'end_deg':self.inclination,
                                'north_pole_phase':self.inclination/360,
                                'south_pole_phase':(self.inclination+180)/360,
                                'observer':'sequence of static observers at infinity; no camera-velocity aberration',
                                'screen_basis':'continuous rotation in a plane perpendicular to the disk',
                                'gas_time_direction':'forward throughout the full revolution'},
                'geometry':self.geometry,'emission':self.source.metadata(),
                'geodesic_reuse':'120 independently traced views; other 120 use exact antipodal isometry with new retarded source evaluation, no spatial interpolation',
                'spin_configuration':{'a_star':self.spin, 'charge_q':self.charge,
                                      'charge_definition':'q=Q/[M sqrt(4 pi epsilon_0 G)]',
                                      'horizon_condition':self.spin**2+self.charge**2,
                                      'matter':'neutral test disk; electromagnetic forces on plasma are not modeled', 'alignment':'aligned, prograde equatorial disk',
                                      'choice':'user-chosen hypothetical spin and large net charge, not a fit to an observed black hole',
                                      'angular_momentum_direction':'+z for both black hole and disk'},
                'layout':{'brand':'@sideratlas','header_top':270,'hero_y':0,'hero_height':1920,
                          'hero_width':1080,'hero_x':0,'full_frame':True,
                          'radiance_shape':[self.tracer.height,self.tracer.width],
                          'display_roll_degrees':self.display_roll,
                          'viewport_fov_horizontal_rg':self.viewport_fov,
                          'text_left':72,'text_right':965,
                          'lower_text_in_caption':True,'caption_file':'caption-kerr-newman.text',
                          'viewport':'fixed rigid image-plane rotation and isotropic crop; outer disk can leave the frame near pole views',
                          'safe_area_note':'Title over the full-frame simulation; interface overlays vary by device and caption state'},
                'encoding':encoding_metadata(),
                'display':{'exposure':EXPOSURE,'curve':'x/(1+x), x=relative radiance/exposure',
                           'color':'original illustrative orange brightness palette, not thermal true color',
                           'averaging':'native float detector radiance; bilinear rotated viewport at 2x output size, BOX averaged before display transform'},
                'scene':'original calculation; no NASA images, frames or source code used',
                'time':'800 t_g shown in each 8-second cycle; same encoded cycle repeated four times for 32 seconds; periodic source forcing, full 360-degree viewpoint revolution; no physical mass assigned',
                'disk_faces':'same isotropic modeled emission field on both faces of a geometrically thin opaque disk',
                'ray_maps':{'inclination_degrees':self.current_inclination,
                            'delay_tg_range':[float(self.maps['delay_tg'][hit].min()),float(self.maps['delay_tg'][hit].max())],
                            'frequency_shift_range':[float(self.maps['g'][hit].min()),float(self.maps['g'][hit].max())],
                            'image_layers':{str(k):int(np.sum(self.maps['hit_layer']==k)) for k in np.unique(self.maps['hit_layer'][hit])}},
                'outputs':VIDEO_NAMES,'loop_outputs':LOOP_VIDEO_NAMES}


def source_hash():
    digest=hashlib.sha256()
    for name in ('render_kn_orbit.py','instagram_export.py','kn_orbit_emission.py','kn_orbit_trace.py','kn_orbit_trace.cpp','nasa_emission.py','kn_orbits.py','render.py','emission.py','raytrace.py'):
        digest.update((ROOT/name).read_bytes())
    for path in sorted((ROOT/'fonts').glob('*.otf')):
        digest.update(path.read_bytes())
    return digest.hexdigest()


def radiance_signature(metadata):
    """Identify raw physical frames independently of typography and encoding."""
    digest=hashlib.sha256()
    for name in ('kn_orbit_emission.py','kn_orbit_trace.py','kn_orbit_trace.cpp','nasa_emission.py','kn_orbits.py','emission.py'):
        digest.update((ROOT/name).read_bytes())
    for function in (view_angle_at_phase,KerrNewmanRenderer.set_view,KerrNewmanRenderer.radiance,KerrNewmanRenderer.antipodal_radiance):
        digest.update(inspect.getsource(function).encode())
    physical={key:metadata[key] for key in
              ('width','height','fps','frames','supersample','inclination_degrees','camera_orbit','emission','spin_configuration')}
    physical['geometry']={k:v for k,v in metadata['geometry'].items() if not k.startswith('last_') and k!='setup_seconds'}
    digest.update(json.dumps(physical,sort_keys=True).encode())
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--supersample',type=float,default=DEFAULT_SUPERSAMPLE)
    parser.add_argument('--preview',action='store_true')
    parser.add_argument('--spin',type=float,default=SPIN)
    parser.add_argument('--charge',type=float,default=CHARGE)
    parser.add_argument('--angle-step',type=float,default=.012)
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR)
    parser.add_argument('--radiance-cache',type=Path,
                        help='reuse a physical-frame cache with a matching provenance manifest')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    print('Tracing Kerr–Newman null geodesics, azimuth and arrival times for the full viewpoint revolution.',flush=True)
    renderer=KerrNewmanRenderer(args.width,args.supersample,spin=args.spin,charge=args.charge,angle_step=args.angle_step)
    first=renderer.frames(0)
    for mode,frame in first.items(): frame.save(args.output/f'sideratlas-kerr-newman-{mode}-cover.png')
    renderer.frames(0,structured=False)['color'].save(args.output/'steady-reference.png')
    sheet=Image.new('RGB',(360*3,640*2))
    for col,phase in enumerate((0,.25,.75)):
        frames=first if phase==0 else renderer.frames(phase)
        for row,mode in enumerate(VIDEO_NAMES):
            sheet.paste(frames[mode].resize((360,640),Image.Resampling.LANCZOS),(360*col,640*row))
    sheet.save(args.output/'contact-sheet.jpg',quality=95)
    metadata=renderer.metadata()
    metadata['source_sha256']=source_hash()
    metadata['supersample']=args.supersample
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    if args.preview:
        print(f'Preview complete in {time.perf_counter()-started:.1f}s.',flush=True)
        return
    signature=radiance_signature(metadata)
    cache=args.radiance_cache or args.output/'.render-cache'/signature[:20]
    cache.mkdir(parents=True,exist_ok=True)
    cache_manifest=None
    if args.radiance_cache:
        manifest_file=cache/'radiance-cache-manifest.json'
        cache_manifest=json.loads(manifest_file.read_text())
        if cache_manifest['radiance_signature'] != signature:
            raise ValueError('cached radiance was calculated with different physical code or settings')
        metadata['radiance_cache_provenance']={
            'radiance_signature':signature,
            'manifest_sha256':hashlib.sha256(manifest_file.read_bytes()).hexdigest(),
            'source_render_sha256':cache_manifest.get('source_render_sha256'),
            'note':'Reuse of original float radiance; no frames decoded from a previous video'}
    metadata['radiance_signature']=signature
    total=FPS*LOOP_DURATION
    encoders={mode:encoder(args.output/(Path(name).stem+'.partial.mp4'),renderer.width,renderer.height,FPS)
              for mode,name in LOOP_VIDEO_NAMES.items()}
    success=False
    try:
        for n in range(total):
            file=cache/f'{n:03d}.npy'
            if file.exists():
                if cache_manifest and hashlib.sha256(file.read_bytes()).hexdigest()!=cache_manifest['frames_sha256'][file.name]:
                    raise ValueError('cached frame does not match its provenance manifest')
                values=np.load(file,allow_pickle=False)
                if values.shape!=(renderer.tracer.height,renderer.tracer.width) or not np.isfinite(values).all() or np.any(values<0):
                    raise ValueError('invalid cached radiance')
            else:
                values=renderer.radiance(n/total)
                temp=file.with_suffix('.partial.npy')
                np.save(temp,values,allow_pickle=False)
                temp.replace(file)
                if n < total//2:
                    opposite_file=cache/f'{n+total//2:03d}.npy'
                    if not opposite_file.exists():
                        opposite_values=renderer.antipodal_radiance(n/total)
                        opposite_temp=opposite_file.with_suffix('.partial.npy')
                        np.save(opposite_temp,opposite_values,allow_pickle=False)
                        opposite_temp.replace(opposite_file)
            for mode,process in encoders.items():
                process.stdin.write(renderer.compose(values,mode).tobytes())
            if n%15==0: print(f'Frame {n+1}/{total}; elapsed {time.perf_counter()-started:.0f}s',flush=True)
        success=True
    finally:
        for process in encoders.values():
            process.stdin.close()
        codes=[process.wait() for process in encoders.values()]
    if not success or any(codes): raise RuntimeError(f'Incomplete video export: {codes}')
    for name in LOOP_VIDEO_NAMES.values():
        (args.output/(Path(name).stem+'.partial.mp4')).replace(args.output/name)
    # Repeat the identical encoded 8-second cycle without another lossy encode.
    # Each cycle starts with an IDR frame and has no B-frames/edit lists.
    for mode,name in VIDEO_NAMES.items():
        subprocess.run([shutil.which('ffmpeg') or 'ffmpeg','-hide_banner','-loglevel','error','-y',
                        '-stream_loop',str(REPEATS-1),'-i',str(args.output/LOOP_VIDEO_NAMES[mode]),
                        '-map','0:v:0','-an','-c:v','copy','-movflags','+faststart','-use_editlist','0',
                        str(args.output/(Path(name).stem+'.partial.mp4'))],check=True)
        (args.output/(Path(name).stem+'.partial.mp4')).replace(args.output/name)
    manifest={'radiance_signature':signature,'source_render_sha256':metadata['source_sha256'],
              'frames_sha256':{f'{n:03d}.npy':hashlib.sha256((cache/f'{n:03d}.npy').read_bytes()).hexdigest()
                               for n in range(total)}}
    (cache/'radiance-cache-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    metadata['repetition_encoding']='Four copies of one 8-second H.264 stream; lossless packet copy, no extra encode'
    metadata['elapsed_seconds']=time.perf_counter()-started
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Finished in {time.perf_counter()-started:.1f}s.',flush=True)


if __name__=='__main__': main()
