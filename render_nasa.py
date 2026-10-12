#!/usr/bin/env python3
"""Original Schwarzschild Reel using NASA's documented geometric features.

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

import numpy as np
from PIL import Image, ImageDraw

from nasa_raytrace import SchwarzschildDelayTracer
from nasa_emission import ShearingDiskEmission
from render import Renderer as FontProvider, ROOT, FPS, DURATION, encode_srgb
from instagram_export import encoder, metadata as encoding_metadata

FOV = 72.0
OUTER_RADIUS = 28.0
INCLINATION = 83.0
LOOP_TG = 800.0
EXPOSURE = 0.65
OUTPUT_DIR = ROOT/'output-sideratlas'
VIDEO_NAMES = {'color': 'sideratlas-black-hole-reel.mp4',
               'mono': 'sideratlas-black-hole-intensity.mp4'}


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


class NasaInspiredRenderer:
    font = FontProvider.font

    def __init__(self, width=1080, supersample=2, inclination=INCLINATION,
                 fov=FOV, outer_radius=OUTER_RADIUS):
        if width < 270 or width % 2 or supersample < 1:
            raise ValueError('even width >=270 and supersample >=1 required')
        self.width, self.height = width, round(width*16/9)//2*2
        self.scale = width/1080
        # The pole view is circular: leave enough black space around the full
        # disk to keep it clear of the title and explanatory text at all phases.
        self.hero_h, self.hero_y = round(900*self.scale), round(425*self.scale)
        # Preserve the physical detector and full field of view. A smaller,
        # isotropic presentation leaves the lower explanation clear of UI.
        self.hero_display_size = (round(864*self.scale), round(720*self.scale))
        self.hero_x = (self.width-self.hero_display_size[0])//2
        self.inclination = inclination
        self.source = ShearingDiskEmission(outer_radius_rg=outer_radius, loop_period_tg=LOOP_TG)
        self.tracer = SchwarzschildDelayTracer(round(width*supersample), round(self.hero_h*supersample),
                                          fov=fov, roll_deg=0, rin=6, rout=outer_radius,
                                          layers=4, radial_samples=8192, angle_step=.002)
        self.geometry = self.tracer.metadata.copy()
        self.current_inclination = None
        self.set_view(0)
        self.overlays = {mode: self._overlay(mode) for mode in VIDEO_NAMES}

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
        white, soft, muted = (245,245,245), (208,208,208), (165,165,165)
        def text(x,y,value,size=28,fill=white,role='regular',anchor='lt'):
            d.text((x*s,y*s),value,font=self.font(size,role),fill=fill,anchor=anchor)
        def line(y):
            d.line((72*s,y*s,965*s,y*s),fill=(55,55,55),width=max(1,round(s)))
        def equation(y,runs,size=31):
            x=72*s
            for value,script in runs:
                face=self.font(size if script is None else size*.64,'math')
                dy=0 if script is None else (size*.45 if script=='sub' else -size*.16)
                d.text((x,(y+dy)*s),value,font=face,fill=soft,anchor='lt')
                x+=d.textlength(value,font=face)
        text(72,270,'@sideratlas',34)
        text(965,277,'NO.001',25,soft,anchor='rt')
        text(64,333,'THE BLACK HOLE',82,role='bold')
        line(411)
        label = 'Illustrative color' if mode=='color' else 'Total-light intensity'
        text(72,1168,'Original relativistic model  –  Schwarzschild / a* = 0  –  '+label,
             24,soft)
        line(1200)
        text(72,1220,'Bends light into multiple images of its disk',31)
        text(72,1256,'Approaching gas appears brighter',31)
        text(72,1292,'Gravity redshifts light and delays its arrival',31)
        text(72,1328,'Inner orbits shear bright structures faster',31)
        equation(1380,[('I',None),('obs','sub'),(' = g',None),('4','sup'),
                       (' I',None),('em','sub'),('     Ω = r',None),('−3/2','sup')])
        text(72,1423,'G = M = c = 1   ·   rISCO = 6 rg   ·   360° viewpoint orbit',25,soft,'math')
        text(72,1462,'Page & Thorne 1974 · Luminet 1979 · Bozza 2002',23,muted)
        text(72,1487,'Müller & Frauendiener 2012 · James et al. 2015',23,muted)
        text(72,1512,'Schnittman & Bertschinger 2004 · Gralla & Lupsasca 2020',23,muted)
        text(72,1537,'Cárdenas-Avendaño et al. 2023 · ICC sRGB',23,muted)
        text(72,1562,'Visual reference: NASA / Goddard / J. Schnittman 2019',23,muted)
        return overlay

    def radiance(self, phase, structured=True):
        self.set_view(phase)
        t = (phase % 1)*LOOP_TG
        values = self.prepared.evaluate(t, structured=structured)*self.transfer
        if not np.isfinite(values).all() or np.any(values < 0):
            raise AssertionError('nonfinite or negative radiance')
        small = np.asarray(Image.fromarray(values.astype('float32')).resize(
            (self.width,self.hero_h),Image.Resampling.BOX))
        return small

    def compose(self, radiance, mode):
        presented = np.asarray(Image.fromarray(radiance.astype('float32')).resize(
            self.hero_display_size,Image.Resampling.BOX))
        level = display_intensity(presented)
        rgb = illustrative_palette(level) if mode=='color' else np.repeat(encode_srgb(level)[...,None],3,axis=2)
        canvas = Image.new('RGB',(self.width,self.height))
        canvas.paste(Image.fromarray(rgb),(self.hero_x,self.hero_y))
        canvas.paste(self.overlays[mode],(0,0),self.overlays[mode])
        return canvas

    def frames(self, phase, modes=tuple(VIDEO_NAMES), structured=True):
        values = self.radiance(phase,structured=structured)
        return {mode:self.compose(values,mode) for mode in modes}

    def metadata(self):
        self.set_view(0)
        hit = self.maps['mask']
        return {'width':self.width,'height':self.height,'fps':FPS,'frames':FPS*DURATION,
                'duration_seconds':DURATION,'inclination_degrees':[0,180],
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
                'layout':{'brand':'@sideratlas','header_top':270,'hero_y':425,'hero_height':720,
                          'hero_width':864,'hero_x':108,'radiance_height':900,
                          'text_left':72,'text_right':965,
                          'last_credit_top':1562,
                          'safe_area_note':'Extra top clearance; interface overlays vary by device and caption state'},
                'encoding':encoding_metadata(),
                'display':{'exposure':EXPOSURE,'curve':'x/(1+x), x=relative radiance/exposure',
                           'color':'original illustrative orange brightness palette, not thermal true color',
                           'averaging':'physical bolometric radiance averaged before display transform'},
                'scene':'original calculation; no NASA images, frames or source code used',
                'time':'800 t_g shown in 8 seconds; periodic source forcing, full 360-degree viewpoint revolution; no physical mass assigned',
                'disk_faces':'same isotropic modeled emission field on both faces of a geometrically thin opaque disk',
                'ray_maps':{'inclination_degrees':self.current_inclination,
                            'delay_tg_range':[float(self.maps['delay_tg'][hit].min()),float(self.maps['delay_tg'][hit].max())],
                            'frequency_shift_range':[float(self.maps['g'][hit].min()),float(self.maps['g'][hit].max())],
                            'image_layers':{str(k):int(np.sum(self.maps['hit_layer']==k)) for k in range(4)}},
                'outputs':VIDEO_NAMES}


def source_hash():
    digest=hashlib.sha256()
    for name in ('render_nasa.py','instagram_export.py','nasa_emission.py','nasa_raytrace.py','raytrace.py','emission.py','render.py'):
        digest.update((ROOT/name).read_bytes())
    for path in sorted((ROOT/'fonts').glob('*.otf')):
        digest.update(path.read_bytes())
    return digest.hexdigest()


def radiance_signature(metadata):
    """Identify raw physical frames independently of typography and encoding."""
    digest=hashlib.sha256()
    for name in ('nasa_emission.py','nasa_raytrace.py','raytrace.py','emission.py'):
        digest.update((ROOT/name).read_bytes())
    for function in (view_angle_at_phase,NasaInspiredRenderer.set_view,NasaInspiredRenderer.radiance):
        digest.update(inspect.getsource(function).encode())
    physical={key:metadata[key] for key in
              ('width','height','fps','frames','supersample','inclination_degrees','camera_orbit','emission')}
    physical['geometry']={k:v for k,v in metadata['geometry'].items() if k!='setup_seconds'}
    digest.update(json.dumps(physical,sort_keys=True).encode())
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--supersample',type=float,default=2)
    parser.add_argument('--preview',action='store_true')
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR)
    parser.add_argument('--radiance-cache',type=Path,
                        help='reuse a physical-frame cache with a matching provenance manifest')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    print('Preparing Schwarzschild orbit tables for the full 360-degree viewpoint revolution.',flush=True)
    renderer=NasaInspiredRenderer(args.width,args.supersample)
    first=renderer.frames(0)
    for mode,frame in first.items(): frame.save(args.output/f'sideratlas-black-hole-{mode}-cover.png')
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
    total=FPS*DURATION
    encoders={mode:encoder(args.output/(Path(name).stem+'.partial.mp4'),renderer.width,renderer.height,FPS)
              for mode,name in VIDEO_NAMES.items()}
    success=False
    try:
        for n in range(total):
            file=cache/f'{n:03d}.npy'
            if file.exists():
                if cache_manifest and hashlib.sha256(file.read_bytes()).hexdigest()!=cache_manifest['frames_sha256'][file.name]:
                    raise ValueError('cached frame does not match its provenance manifest')
                values=np.load(file,allow_pickle=False)
                if values.shape!=(renderer.hero_h,renderer.width) or not np.isfinite(values).all() or np.any(values<0):
                    raise ValueError('invalid cached radiance')
            else:
                values=renderer.radiance(n/total)
                temp=file.with_suffix('.partial.npy')
                np.save(temp,values,allow_pickle=False)
                temp.replace(file)
            for mode,process in encoders.items():
                process.stdin.write(renderer.compose(values,mode).tobytes())
            if n%15==0: print(f'Frame {n+1}/{total}; elapsed {time.perf_counter()-started:.0f}s',flush=True)
        success=True
    finally:
        for process in encoders.values():
            process.stdin.close()
        codes=[process.wait() for process in encoders.values()]
    if not success or any(codes): raise RuntimeError(f'Incomplete video export: {codes}')
    for name in VIDEO_NAMES.values():
        (args.output/(Path(name).stem+'.partial.mp4')).replace(args.output/name)
    metadata['elapsed_seconds']=time.perf_counter()-started
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Finished in {time.perf_counter()-started:.1f}s.',flush=True)


if __name__=='__main__': main()
