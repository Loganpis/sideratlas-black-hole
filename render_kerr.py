#!/usr/bin/env python3
"""Render the second @atlastro Reel: Kerr geodesics and visible spectral color.

No gas texture or luminous ring is painted onto the solution. The optional
camera adds only a normalized, explicitly specified diffraction PSF. See
README-KERR.md for the changed physical parameters and all approximations.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import tempfile
import time

import numpy as np
from PIL import Image, ImageDraw

from emission import ThermalSpectrum
from kerr_emission import KerrDiskParameters, KerrEmission
from kerr_raytrace import KerrTracer
from optics import atmosphere_rgb, chromaticity_preserving_tonemap, AiryCamera, CameraParameters
from render import Renderer as OriginalRenderer, ROOT, FPS, DURATION, encoder, encode_srgb

MODES = ('camera', 'ideal')
VIDEO_NAMES = {'camera':'atlastro-kerr-truecolor-reel.mp4',
               'ideal':'atlastro-kerr-ideal-optics-reel.mp4'}
EXPOSURE = 1 / 1.8e8


def inclination(phase):
    return 83.0 + 3.0*math.cos(2*math.pi*(phase % 1))


class KerrRenderer:
    font = OriginalRenderer.font

    def __init__(self, width=1080, supersample=1.5, camera=True):
        if width < 270 or width % 2 or supersample < 1:
            raise ValueError('width must be even and >=270; supersample must be >=1')
        self.width = width
        self.height = round(width*16/9)//2*2
        self.scale = width/1080
        self.hero_h = round(940*self.scale)
        self.hero_y = round(322*self.scale)
        self.parameters = KerrDiskParameters(mass_solar=1e10,eddington_ratio=.05,spin_a_star=.9)
        self.disk = KerrEmission(self.parameters)
        self.spectrum = ThermalSpectrum()
        self.tracer = KerrTracer(round(width*supersample),round(self.hero_h*supersample),
                                 fov=26,roll_deg=14,spin=.9,rout=500,angle_step=.006)
        self.camera = AiryCamera(width,self.hero_h,fov_rg=26,parameters=CameraParameters()) if camera else None
        self.overlays = {mode:self._overlay(mode) for mode in MODES}
        self.last_stats = {}

    def _overlay(self,mode):
        image = Image.new('RGBA',(self.width,self.height))
        d,s = ImageDraw.Draw(image),self.scale
        white,soft,muted = (245,245,245),(214,214,214),(170,170,170)
        def text(x,y,value,size=28,fill=white,role='regular',anchor='lt'):
            d.text((x*s,y*s),value,font=self.font(size,role),fill=fill,anchor=anchor)
        def line(y):
            d.line((72*s,y*s,965*s,y*s),fill=(62,62,62),width=max(1,round(s)))
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
        text(965,269,'KERR  /  a* = 0.90',25,soft,role='math',anchor='rt')
        line(309)
        # A clearly bounded scientific viewport: the extended disk may leave
        # the field of view, but is never faded or cut off at a false radius.
        text(72,1287,'CALCULATED VISIBLE COLOR',21,soft)
        text(965,1287,'STEADY DISK · VIEW SCAN',21,soft,anchor='rt')
        line(1328)
        text(72,1350,'Signature Move',28,soft,'italic')
        text(72,1394,'It drags spacetime and bends light around it.',31)
        text(72,1436,'Its orbiting disk is Doppler-brightened and shifted.',31)
        text(72,1478,'Its disk atmosphere changes color and brightness.',31)
        text(72,1520,'It traps light inside its horizon.',31)
        equation(1571,[('T',None),('obs','sub'),(' = g f',None),('col','sub'),(' T',None),('eff','sub'),
                       ('     r',None),('+','sub'),(' = 1.44r',None),('g','sub')],30)
        text(72,1629,'Page & Thorne 1974 · Bardeen et al. 1972 · Li et al. 2005',23,muted)
        text(72,1661,'Gralla & Lupsasca 2020 · Done et al. 2012',23,muted)
        text(72,1693,'Luminet 1979 · Müller & Frauendiener 2012',23,muted)
        text(72,1725,'James et al. 2015 · Bozza 2002 · Airy 1835 · CIE 2019 / ICC',23,muted)
        text(72,1774,'M = 10¹⁰ M⊙   ·   visible spectrum: 360–830 nm',24,soft,'math')
        label='20 m ideal camera · 550 nm diffraction approximation' if mode=='camera' else 'Ideal angular resolution · no camera diffraction'
        text(72,1816,label,22,muted,'italic')
        return image

    def physical_images(self,angle,modes=MODES):
        m=self.tracer.trace(angle)
        temp=self.disk.temperature(m['r'])
        rgb=atmosphere_rgb(temp,m['g'],m['mu'],self.spectrum)
        self.last_stats={'angle_deg':angle,'maximum_local_temperature_k':float(temp.max()),
                         'maximum_spectral_rgb_radiance':float(rgb.max()),
                         'emitting_fraction':float(m['mask'].mean()),
                         'constraint_residual_max':float(m['constraint_residual'].max()),
                         'visible_crossing_layers':{str(k):int(np.sum(m['hit_layer']==k))
                                                   for k in np.unique(m['hit_layer'][m['mask']])}}
        # Detector subpixels average energy, before the display transfer.
        small=np.stack([np.asarray(Image.fromarray(rgb[...,c].astype('float32')).resize(
            (self.width,self.hero_h),Image.Resampling.BOX)) for c in range(3)],axis=-1)
        output={}
        for mode in modes:
            linear=self.camera.apply(small) if mode=='camera' and self.camera else small
            output[mode]=Image.fromarray(encode_srgb(chromaticity_preserving_tonemap(linear,EXPOSURE)))
        return output

    def frames(self,phase,modes=MODES):
        angle=inclination(phase)
        images={}
        for mode,hero in self.physical_images(angle,modes).items():
            frame=Image.new('RGB',(self.width,self.height),(0,0,0))
            frame.paste(hero,(0,self.hero_y))
            frame.paste(self.overlays[mode],(0,0),self.overlays[mode])
            ImageDraw.Draw(frame).text((72*self.scale,1207*self.scale),f'i = {angle:.1f}°',
                font=self.font(24,'math'),fill=(214,214,214),anchor='lt')
            images[mode]=frame
        return images


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--supersample',type=float,default=1.5)
    parser.add_argument('--preview',action='store_true')
    parser.add_argument('--mode',choices=['camera','ideal','both'],default='both')
    parser.add_argument('--output',type=Path,default=ROOT/'output-kerr')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    modes=MODES if args.mode=='both' else (args.mode,)
    started=time.perf_counter()
    print('Preparing Kerr, atmosphere, spectrum and diffraction calculations.',flush=True)
    renderer=KerrRenderer(args.width,args.supersample,camera='camera' in modes)
    start_frames=renderer.frames(0,modes)
    start_stats=renderer.last_stats.copy()
    for mode,frame in start_frames.items():
        frame.save(args.output/f'atlastro-kerr-{mode}-cover.png')
    sheet=Image.new('RGB',(360*3,640*len(modes)))
    for column,phase in enumerate((0,.25,.5)):
        frames=start_frames if phase==0 else renderer.frames(phase,modes)
        for row,mode in enumerate(modes):
            sheet.paste(frames[mode].resize((360,640),Image.Resampling.LANCZOS),(column*360,row*640))
    sheet.save(args.output/'atlastro-kerr-contact-sheet.jpg',quality=95)
    if args.preview:
        print('Preview complete.',flush=True)
        return
    total=FPS*DURATION
    encoders={mode:encoder(args.output/VIDEO_NAMES[mode],renderer.width,renderer.height,FPS) for mode in modes}
    all_stats=[]
    with tempfile.TemporaryDirectory(prefix='atlastro-kerr-') as temporary:
        cache=Path(temporary)
        try:
            for n in range(total):
                key=min(n,total-n)
                if n<=total//2:
                    frames=start_frames if n==0 else renderer.frames(n/total,modes)
                    all_stats.append(start_stats.copy() if n==0 else renderer.last_stats.copy())
                    for mode,frame in frames.items():
                        frame.save(cache/f'{mode}-{key:03d}.png',compress_level=1)
                else:
                    frames={mode:Image.open(cache/f'{mode}-{key:03d}.png').convert('RGB') for mode in modes}
                for mode,frame in frames.items():
                    encoders[mode].stdin.write(frame.tobytes())
                if n%10==0:
                    print(f'Frame {n+1}/{total}; elapsed {time.perf_counter()-started:.0f}s',flush=True)
        finally:
            for process in encoders.values(): process.stdin.close()
            codes=[process.wait() for process in encoders.values()]
        if any(codes): raise RuntimeError(f'Encoding failed: {codes}')
    metadata={'width':renderer.width,'height':renderer.height,'frames':total,'fps':FPS,
              'duration_seconds':DURATION,'supersample':args.supersample,
              'animation':'Independent static observers, i=83+3cos(2pi phase); no evolving fluid',
              'disk':renderer.disk.metadata(),'raytracer':renderer.tracer.metadata,
              'spectrum':renderer.spectrum.metadata,'camera':renderer.camera.metadata if renderer.camera else None,
              'atmosphere':'Done2012 color correction, flux-normalized Li2005 limb law; approximate diluted Planck spectrum',
              'display':{'exposure':EXPOSURE,'curve':'v/(1+max(v)) with common RGB divisor, then sRGB',
                         'subpixel_average':'linear physical spectral radiance',
                         'color_tags':'BT709 primaries/matrix, sRGB transfer'},
              'frame_stats':all_stats,'elapsed_seconds':time.perf_counter()-started}
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Finished in {time.perf_counter()-started:.1f}s.',flush=True)


if __name__=='__main__':
    main()
