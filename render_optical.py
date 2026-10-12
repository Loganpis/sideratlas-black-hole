#!/usr/bin/env python3
"""Render the wavelength-resolved @atlastro Kerr black-hole Reel.

This version preserves the Kerr disk's physical parameters, computes camera
diffraction separately at every visible wavelength, and traces an overscan
margin. Its fixed display contrast curve preserves in-gamut RGB ratios.
See README-OPTICAL.md for the model, limitations, and every source used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from PIL import Image, ImageDraw

from kerr_emission import KerrDiskParameters, KerrEmission
from kerr_raytrace import KerrTracer
from optics import color_correction, limb_factor
from spectral_optics import SpectralCamera, FIRST_NULL_FACTOR
from render import Renderer as BaseRenderer, ROOT, FPS, DURATION, encoder, encode_srgb

EXPOSURE = 1/1.8e8
FOV = 24.
APERTURE = 10.
DISTANCE = 1e6
SUPPORT = 24.
RAY_STEP = .018
MODES = ('camera', 'ideal')
VIDEO_NAMES = {'camera':'atlastro-truecolor-optics-reel.mp4',
               'ideal':'atlastro-truecolor-ideal-reel.mp4'}


def inclination(phase):
    return 86.5 + 1.5*math.cos(2*math.pi*(phase % 1))


def display_rgb(radiance, exposure=EXPOSURE):
    """A fixed display contrast curve, not an astrophysical emission law.

    v=E*max(RGB,0), m=max(v), output=v*m/(1+m²). The common multiplier
    preserves in-gamut chromaticity. It expands dark contrast and compresses
    highlights without adding light, changing exposure, or painting a halo.
    """
    values=np.maximum(np.asarray(radiance,dtype=np.float64),0)*exposure
    peak=values.max(axis=-1,keepdims=True)
    return encode_srgb(values*peak/(1+peak*peak))


class OpticalRenderer:
    font = BaseRenderer.font

    def __init__(self,width=1080,supersample=1.5,wavelength_step_nm=10,fast_spectrum=True):
        if width<270 or width%2 or supersample<1:
            raise ValueError('width must be even and >=270; supersample must be >=1')
        self.width=width
        self.height=round(width*16/9)//2*2
        self.scale=width/1080
        self.fast_spectrum=fast_spectrum
        self.hero_h=round(1620*self.scale)
        self.hero_y=round(30*self.scale)
        self.parameters=KerrDiskParameters(mass_solar=1e10,eddington_ratio=.05,spin_a_star=.9)
        self.disk=KerrEmission(self.parameters)
        # This extra field supplies every photon that can reach the visible
        # viewport through the finite-support camera PSF. It is removed only
        # AFTER diffraction, never black-padded at the final image boundary.
        self.padding=max(2,math.ceil(FIRST_NULL_FACTOR*830e-9/APERTURE /
                                    (FOV/(DISTANCE*width))*SUPPORT))
        camera_w=width+2*self.padding
        camera_h=self.hero_h+2*self.padding
        camera_fov=FOV*camera_w/width
        self.camera=SpectralCamera(camera_w,camera_h,camera_fov,
                                  aperture_m=APERTURE,observer_distance_rg=DISTANCE,
                                  wavelength_step_nm=wavelength_step_nm,support_first_nulls=SUPPORT)
        if self.camera.padding_pixels>self.padding:
            raise AssertionError('traced overscan is smaller than the diffraction support')
        self.tracer=KerrTracer(round(camera_w*supersample),round(camera_h*supersample),
                               fov=camera_fov,roll_deg=14,spin=.9,rout=500,angle_step=RAY_STEP)
        self.overlays={mode:self._overlay(mode) for mode in MODES}
        self.last_stats={}

    def _overlay(self,mode):
        image=Image.new('RGBA',(self.width,self.height))
        d,s=ImageDraw.Draw(image),self.scale
        white,soft,muted=(245,245,245),(210,210,210),(170,170,170)
        def text(x,y,value,size=28,fill=white,role='regular',anchor='lt'):
            d.text((x*s,y*s),value,font=self.font(size,role),fill=fill,anchor=anchor)
        def line(y):
            d.line((72*s,y*s,965*s,y*s),fill=(60,60,60),width=max(1,round(s)))
        def equation(y,runs,size=30):
            x=72*s
            for value,script in runs:
                face=self.font(size if script is None else size*.64,'math')
                dy=0 if script is None else (size*.45 if script=='sub' else -size*.16)
                d.text((x,(y+dy)*s),value,font=face,fill=soft,anchor='lt')
                x+=d.textlength(value,font=face)
        text(72,79,'@atlastro',34)
        text(965,85,'NO.001',25,soft,anchor='rt')
        text(64,156,'THE BLACK HOLE',84,role='bold')
        text(72,268,'GRAVITY TYPE',25,soft)
        text(965,269,'KERR  /  a* = 0.90',25,soft,'math',anchor='rt')
        line(309)
        text(72,1269,'CALCULATED VISIBLE COLOR',21,soft)
        text(965,1269,'STEADY DISK · VIEW SCAN',21,soft,anchor='rt')
        line(1310)
        text(72,1337,'Signature Move',28,soft,'italic')
        text(72,1382,'It bends light and drags spacetime around it.',31)
        text(72,1423,'Its orbiting gas boosts and shifts the disk light.',31)
        text(72,1464,'Its atmosphere shapes the spectrum; its horizon traps light.',29)
        camera_line='The camera diffracts its light, wavelength by wavelength.'
        if mode=='ideal': camera_line='This reference shows its light without camera diffraction.'
        text(72,1505,camera_line,29)
        equation(1560,[('T',None),('obs','sub'),(' = g f',None),('col','sub'),(' T',None),('eff','sub'),
                       ('     r',None),('+','sub'),(' = 1.44r',None),('g','sub')],30)
        text(72,1623,'Page & Thorne 1974 · Bardeen et al. 1972 · Li et al. 2005',23,muted)
        text(72,1655,'Gralla & Lupsasca 2020 · Done et al. 2012',23,muted)
        text(72,1687,'Luminet 1979 · Müller & Frauendiener 2012',23,muted)
        text(72,1719,'James et al. 2015 · Bozza 2002 · Airy 1835 · CIE 2019 / ICC',23,muted)
        text(72,1770,'M = 10¹⁰ M⊙   ·   visible spectrum: 360–830 nm',24,soft,'math')
        footer='10 m circular aperture · wavelength-dependent diffraction'
        if mode=='ideal': footer='Ideal angular resolution · same spectrum and display curve'
        text(72,1812,footer,22,muted,'italic')
        return image

    def physical_images(self,angle,modes=MODES):
        step=RAY_STEP
        self.tracer.angle_step=step
        self.tracer.max_steps=math.ceil(20*math.pi/step)
        try:
            m=self.tracer.trace(angle)
        except RuntimeError:
            m=None
        if m is None or m['constraint_residual'].max()>1e-5:
            step=.006
            self.tracer.angle_step=step
            self.tracer.max_steps=math.ceil(20*math.pi/step)
            m=self.tracer.trace(angle)
        if m['constraint_residual'].max()>1e-5:
            raise AssertionError('Kerr potential constraint did not meet the tolerance after refinement')
        temperature=self.disk.temperature(m['r'])
        fcol=color_correction(temperature)
        observed=m['g']*fcol*temperature
        dilution=limb_factor(m['mu'])/fcol**4
        if 'camera' in modes:
            evaluate=self.camera.apply_fast if self.fast_spectrum else self.camera.apply
            images=evaluate(observed,dilution,return_ideal=True)
        else:
            images={'ideal':self.camera.apply(observed,dilution,diffraction=False)}
        pad=self.padding
        crop=(slice(pad,pad+self.hero_h),slice(pad,pad+self.width))
        use=m['mask']
        self.last_stats={'angle_deg':angle,
                         'ray_integration_step_used':step,
                         'constraint_residual_max':float(m['constraint_residual'].max()),
                         'temperature_max_k':float(temperature.max()),
                         'observed_color_temperature_max_k':float(observed.max()),
                         'frequency_shift_min':float(m['g'][use].min()),
                         'frequency_shift_max':float(m['g'][use].max()),
                         'visible_crossing_layers':{str(k):int(np.sum(m['hit_layer']==k))
                                                   for k in np.unique(m['hit_layer'][use])},
                         'emission_outside_isco':bool(np.all(m['r'][use]>=self.parameters.isco_rg-1e-6)),
                         'spectrum_camera_seconds':self.camera.metadata.get('last_apply_seconds')}
        result={}
        for mode in modes:
            linear=images[mode][crop]
            if not np.isfinite(linear).all(): raise AssertionError('nonfinite output radiance')
            self.last_stats[f'{mode}_radiance_max']=float(linear.max())
            result[mode]=Image.fromarray(display_rgb(linear))
        return result

    def frames(self,phase,modes=MODES):
        angle=inclination(phase)
        frames={}
        for mode,hero in self.physical_images(angle,modes).items():
            frame=Image.new('RGB',(self.width,self.height),(0,0,0))
            frame.paste(hero,(0,self.hero_y))
            frame.paste(self.overlays[mode],(0,0),self.overlays[mode])
            ImageDraw.Draw(frame).text((72*self.scale,1206*self.scale),f'i = {angle:.1f}°',
                                      font=self.font(24,'math'),fill=(210,210,210),anchor='lt')
            frames[mode]=frame
        return frames


def _atomic_json(path,data):
    temporary=path.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(path)


def _cache_signature(args,modes):
    names=['render_optical.py','spectral_optics.py','kerr_raytrace.py','kerr_raytrace.cpp',
           'kerr_emission.py','emission.py','optics.py','render.py']
    files=[ROOT/name for name in names]+sorted((ROOT/'fonts').glob('*.otf'))+[
        ROOT/'data'/'CIE_xyz_1931_2deg.csv']
    parameters={'width':args.width,'supersample':args.supersample,'wavelength_step':args.wavelength_step,
                'modes':list(modes),'fast_spectrum':not args.direct_spectrum,
                'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    digest=hashlib.sha256(json.dumps(parameters,sort_keys=True).encode()).hexdigest()
    return digest,parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--supersample',type=float,default=1.5)
    parser.add_argument('--wavelength-step',type=float,default=10)
    parser.add_argument('--preview',action='store_true')
    parser.add_argument('--direct-spectrum',action='store_true',help='use the slower direct wavelength sum for reference')
    parser.add_argument('--mode',choices=['camera','ideal','both'],default='both')
    parser.add_argument('--output',type=Path,default=ROOT/'output-optical')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    modes=MODES if args.mode=='both' else (args.mode,)
    start=time.perf_counter()
    print('Preparing Kerr rays and wavelength-resolved camera optics.',flush=True)
    renderer=OpticalRenderer(args.width,args.supersample,args.wavelength_step,not args.direct_spectrum)
    count=FPS*DURATION
    signature,parameters=_cache_signature(args,modes)
    cache=args.output/'.render-cache'/signature[:16]
    cache.mkdir(parents=True,exist_ok=True)
    _atomic_json(cache/'inputs.json',parameters)
    targets=(0,count//4,count//2) if args.preview else range(count//2+1)
    stats=[]
    reused=0
    for n in targets:
        record_path=cache/f'{n:03d}.json'
        record=None
        if record_path.exists():
            candidate=json.loads(record_path.read_text())
            intact=all((cache/item['name']).exists() and
                       hashlib.sha256((cache/item['name']).read_bytes()).hexdigest()==item['sha256']
                       for item in candidate['images'].values())
            if intact and set(candidate['images'])==set(modes): record=candidate
        if record is None:
            frames=renderer.frames(n/count,modes)
            images={}
            for mode,frame in frames.items():
                path=cache/f'{mode}-{n:03d}.png'
                temporary=path.with_suffix('.tmp.png')
                frame.save(temporary,compress_level=1)
                temporary.replace(path)
                images[mode]={'name':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
            record={'stats':renderer.last_stats.copy(),'images':images}
            _atomic_json(record_path,record)
            _atomic_json(cache/'engine.json',{'camera':renderer.camera.metadata,'raytracer':renderer.tracer.metadata})
        else:
            reused+=1
        stats.append(record['stats'])
        if n==0:
            for mode in modes:
                with Image.open(cache/record['images'][mode]['name']) as frame:
                    frame.save(args.output/f'atlastro-optical-{mode}-cover.png')
        if n%5==0 or args.preview:
            print(f'View {n+1}/{count//2+1}; {reused} restored; elapsed {time.perf_counter()-start:.0f}s',flush=True)
    sheet=Image.new('RGB',(360*3,640*len(modes)))
    for column,n in enumerate((0,count//4,count//2)):
        for row,mode in enumerate(modes):
            with Image.open(cache/f'{mode}-{n:03d}.png') as frame:
                sheet.paste(frame.resize((360,640),Image.Resampling.LANCZOS),(column*360,row*640))
    sheet.save(args.output/'atlastro-optical-contact-sheet.jpg',quality=95)
    print(f'Previews complete after {time.perf_counter()-start:.1f}s.',flush=True)
    if args.preview: return
    temporary_paths={mode:args.output/VIDEO_NAMES[mode].replace('.mp4','.rendering.mp4') for mode in modes}
    outputs={mode:encoder(temporary_paths[mode],renderer.width,renderer.height,FPS) for mode in modes}
    print('Encoding the two complete loops from saved views.',flush=True)
    try:
        for n in range(count):
            key=min(n,count-n)
            for mode in modes:
                with Image.open(cache/f'{mode}-{key:03d}.png') as frame:
                    outputs[mode].stdin.write(frame.convert('RGB').tobytes())
    finally:
        for process in outputs.values(): process.stdin.close()
        codes=[p.wait() for p in outputs.values()]
    if any(codes): raise RuntimeError(f'Encoder failure: {codes}')
    for mode in modes: temporary_paths[mode].replace(args.output/VIDEO_NAMES[mode])
    engine=json.loads((cache/'engine.json').read_text())
    metadata={'width':renderer.width,'height':renderer.height,'fps':FPS,'frames':count,'duration_seconds':DURATION,
              'supersample':args.supersample,'disk':renderer.disk.metadata(),'raytracer':engine['raytracer'],
              'camera':engine['camera'],'checkpoint_signature':signature,'restored_views':reused,
              'composition':{'visible_fov_horizontal_rg':FOV,'hero_height':renderer.hero_h,'hero_y':renderer.hero_y,
                             'overscan_each_edge_pixels':renderer.padding},
              'animation':'85–88 degree cosine view study, i=86.5+1.5cos(2pi phase); independent static observers',
              'display':{'exposure':EXPOSURE,'curve':'v=max(RGB,0)*E; m=max(v); out=v*m/(1+m^2); sRGB encoding',
                         'color_role':'calculated visible color for the specified model; fixed display contrast, no automatic white balance',
                         'source_radiance':'diluted shifted Planck spectrum, not bolometric false color'},
              'frame_stats':stats,'elapsed_seconds':time.perf_counter()-start,
              'outputs':{mode:VIDEO_NAMES[mode] for mode in modes}}
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Finished in {time.perf_counter()-start:.1f}s.',flush=True)


if __name__=='__main__':
    main()
