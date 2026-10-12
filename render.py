#!/usr/bin/env python3
"""Render @atlastro's reproducible, relativistic black-hole view study.

The main Reel shows bolometric radiance in monochrome. A companion shows
the CIE-integrated visible thermal spectrum. Both use the same geodesics,
Page--Thorne disk and one fixed exposure, with no invented texture or glow.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from emission import DEFAULT_PARAMETERS, ThermalSpectrum, bolometric_intensity, observed_temperature
from raytrace import SchwarzschildTracer

ROOT = Path(__file__).resolve().parent
FPS, DURATION = 30, 8
BOLO_WHITE = 1.30e12  # W m^-2 sr^-1 mapped to display white, constant for every frame
THERMAL_WHITE = 1.70e9  # CIE-weighted radiance mapped to linear-sRGB 1
WHITE, SOFT, MUTED = (243, 243, 243), (204, 204, 204), (156, 156, 156)


def inclination(phase):
    """A periodic viewing study of a steady solution, not fluid evolution."""
    return 82.5 + 4.5 * math.cos(2 * math.pi * (phase % 1))


def encode_srgb(linear):
    linear = np.clip(linear, 0, 1)
    encoded = np.where(linear <= .0031308, 12.92*linear,
                       1.055*np.power(linear, 1/2.4)-.055)
    return np.uint8(np.rint(np.clip(encoded, 0, 1)*255))


class Renderer:
    def __init__(self, width=1080, supersample=1.5):
        if width < 270 or width % 2:
            raise ValueError('width must be an even integer >=270')
        if supersample < 1:
            raise ValueError('supersample must be >=1')
        self.width = width
        self.height = int(round(width*16/9)) // 2 * 2
        self.scale = width/1080
        self.hero_h = int(round(1200*self.scale))
        self.hero_y = int(round(230*self.scale))
        self.tracer = SchwarzschildTracer(
            int(round(width*supersample)), int(round(self.hero_h*supersample)),
            fov=32, roll_deg=15, rin=6, rout=40,
            radial_samples=8192, angle_step=.002, layers=4)
        self.spectrum = ThermalSpectrum()
        self.overlays = {mode: self._overlay(mode) for mode in ('bolometric', 'thermal')}
        self.last_stats = {}

    def font(self, size, role='regular'):
        names = {'regular':'lmroman10-regular.otf', 'bold':'lmroman12-bold.otf',
                 'italic':'lmroman10-italic.otf', 'math':'STIXTwoMath-Regular.otf'}
        return ImageFont.truetype(str(ROOT/'fonts'/names[role]), max(8, round(size*self.scale)))

    def _overlay(self, mode):
        overlay = Image.new('RGBA', (self.width, self.height))
        d = ImageDraw.Draw(overlay)
        s = self.scale

        def text(x, y, value, size, color=WHITE, role='regular', anchor='lt'):
            face = self.font(size, role)
            d.text((round(x*s), round(y*s)), value, font=face, fill=color, anchor=anchor)

        def line(y, color=(65,65,65)):
            d.line((72*s,y*s,965*s,y*s), fill=color, width=max(1,round(s)))

        def mathline(y, runs, size=30):
            x = 72*s
            for value, script in runs:
                face = self.font(size if script is None else size*.66, 'math')
                offset = 0 if script is None else (size*.46 if script=='sub' else -size*.15)
                d.text((x,(y+offset)*s),value,font=face,fill=SOFT,anchor='lt')
                x += d.textlength(value,font=face)

        text(72, 79, '@atlastro', 34)
        text(965, 85, 'NO.001', 25, SOFT, anchor='rt')
        text(64, 156, 'THE BLACK HOLE', 84, role='bold')
        text(72, 268, 'GRAVITY TYPE', 25, SOFT)
        line(315)

        # These labels sit in empty sky; they never paint a halo onto the disk.
        text(72, 345, 'Schwarzschild', 24, MUTED, 'italic')
        text(965, 345, 'a* = 0', 26, MUTED, 'math', anchor='rt')

        text(72, 1265, 'STEADY THIN DISK', 21, SOFT)
        text(965, 1265, 'VIEWING-ANGLE SCAN', 21, SOFT, anchor='rt')
        line(1312)
        text(72, 1337, 'Signature Move', 29, SOFT, 'italic')
        text(72, 1387, 'It bends light into multiple images of its hot disk.', 32)
        text(72, 1431, 'Its orbiting gas is brighter on the approaching side.', 32)
        text(72, 1475, 'Its gravity redshifts light; its horizon traps it.', 32)
        mathline(1531,[('r',None),('h','sub'),(' = 2r',None),('g','sub'),
                       ('     r',None),('ISCO','sub'),(' = 6r',None),('g','sub'),
                       ('     r',None),('g','sub'),(' = GM/c',None),('2','sup')])
        if mode == 'bolometric':
            mathline(1580,[('I',None),('obs','sub'),(' = g',None),('4','sup'),
                           (' I',None),('em','sub'),('   ·   total radiated light, fixed exposure',None)],27)
        else:
            mathline(1580,[('T',None),('obs','sub'),(' = g T',None),('eff','sub'),
                           ('   ·   visible thermal light, fixed exposure',None)],27)
        text(72, 1640, 'Page & Thorne (1974) · Luminet (1979) · Bozza (2002)', 25, MUTED)
        text(72, 1675, 'Müller & Frauendiener (2012) · James et al. (2015)', 25, MUTED)
        text(72, 1710, 'Thermal color: CIE (2019) · ICC sRGB', 25, MUTED)
        text(72, 1771, 'M = 10⁸ M⊙   ·   idealized relativistic model', 25, SOFT, 'math')
        text(72, 1810, 'Studies & model assumptions in the caption.', 23, MUTED, 'italic')
        return overlay

    def physical_images(self, angle, modes):
        m = self.tracer.trace(angle)
        output = {}
        for mode in modes:
            if mode == 'bolometric':
                values = bolometric_intensity(m['r'], m['g']) / BOLO_WHITE
            else:
                values = self.spectrum.lookup_rgb_linear(observed_temperature(m['r'], m['g'])) / THERMAL_WHITE
            self.last_stats[mode] = {'max_linear_display_value':float(values.max()),
                                     'negative_fraction':float(np.mean(values<0)),
                                     'over_white_fraction':float(np.mean(values>1))}
            # Average physical radiance across each detector pixel BEFORE
            # the nonlinear sRGB transfer. No bloom, rim, vignette or stars.
            if values.ndim == 2:
                small = np.asarray(Image.fromarray(values.astype('float32')).resize(
                    (self.width,self.hero_h), Image.Resampling.BOX))
                rgb = np.repeat(encode_srgb(small)[...,None], 3, axis=2)
            else:
                small = np.stack([np.asarray(Image.fromarray(values[...,c].astype('float32')).resize(
                    (self.width,self.hero_h), Image.Resampling.BOX)) for c in range(3)], axis=-1)
                rgb = encode_srgb(small)
            output[mode] = Image.fromarray(rgb)
        return output

    def compose(self, hero, mode, angle):
        canvas = Image.new('RGB', (self.width,self.height), (0,0,0))
        canvas.paste(hero, (0,self.hero_y))
        canvas.paste(self.overlays[mode], (0,0), self.overlays[mode])
        d = ImageDraw.Draw(canvas)
        # Numerical readout makes the source of the animation unambiguous.
        d.text((72*self.scale,1218*self.scale), f'i = {angle:.1f}°',
               font=self.font(25,'math'), fill=MUTED, anchor='lt')
        return canvas

    def frames(self, phase, modes):
        angle = inclination(phase)
        return {mode:self.compose(hero,mode,angle)
                for mode,hero in self.physical_images(angle,modes).items()}


def encoder(path, width, height, fps):
    ffmpeg = shutil.which('ffmpeg')
    if not ffmpeg:
        raise RuntimeError('FFmpeg is required on PATH to export the video.')
    command = [ffmpeg,'-hide_banner','-loglevel','error','-y',
               '-f','rawvideo','-pix_fmt','rgb24','-s',f'{width}x{height}',
               '-r',str(fps),'-i','-','-an','-c:v','libx264','-preset','medium',
               '-crf','17','-maxrate','18M','-bufsize','36M','-pix_fmt','yuv420p',
               '-profile:v','high','-g',str(fps*2),'-keyint_min',str(fps*2),
               '-bf','0','-sc_threshold','0','-flags','+cgop',
               '-color_primaries','bt709','-color_trc','iec61966-2-1',
               '-colorspace','bt709','-sws_flags','accurate_rnd+full_chroma_int',
               '-vf','scale=out_color_matrix=bt709:out_range=tv',
               '-bsf:v','h264_metadata=colour_primaries=1:transfer_characteristics=13:matrix_coefficients=1',
               '-movflags','+faststart','-use_editlist','0',str(path)]
    return subprocess.Popen(command,stdin=subprocess.PIPE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--supersample',type=float,default=1.5)
    parser.add_argument('--preview',action='store_true',help='save a poster and contact sheet only')
    parser.add_argument('--mode',choices=['bolometric','thermal','both'],default='both')
    parser.add_argument('--output',type=Path,default=ROOT/'output')
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    modes = ['bolometric','thermal'] if args.mode=='both' else [args.mode]
    started=time.perf_counter()
    print('Preparing the geodesic and thermal-spectrum tables.',flush=True)
    renderer=Renderer(args.width,args.supersample)
    frames0=renderer.frames(0,modes)
    for mode,frame in frames0.items():
        frame.save(args.output/f'atlastro-{mode}-cover.png')
    # A poster view at each extremum and at the midpoint.
    sheet=Image.new('RGB',(360*3,640*len(modes)))
    for col,phase in enumerate((0,.25,.5)):
        frames=frames0 if phase==0 else renderer.frames(phase,modes)
        for row,mode in enumerate(modes):
            sheet.paste(frames[mode].resize((360,640),Image.Resampling.LANCZOS),(col*360,row*640))
    sheet.save(args.output/'atlastro-contact-sheet.jpg',quality=95)
    if args.preview:
        print('Preview complete.',flush=True)
        return
    total=FPS*DURATION
    names={'bolometric':'atlastro-black-hole-reel.mp4','thermal':'atlastro-visible-thermal-reel.mp4'}
    encoders={mode:encoder(args.output/names[mode],renderer.width,renderer.height,FPS) for mode in modes}
    # cos(theta)=cos(2pi-theta): lossless caching avoids tracing the same
    # viewpoint twice, without quantizing the camera or holding frames.
    with tempfile.TemporaryDirectory(prefix='atlastro-frames-') as temp:
        cache=Path(temp)
        try:
            for n in range(total):
                key=min(n,total-n)
                if n<=total//2:
                    frames=frames0 if n==0 else renderer.frames(n/total,modes)
                    for mode,frame in frames.items():
                        frame.save(cache/f'{mode}-{key:03d}.png',compress_level=1)
                else:
                    frames={mode:Image.open(cache/f'{mode}-{key:03d}.png').convert('RGB') for mode in modes}
                for mode,frame in frames.items():
                    encoders[mode].stdin.write(frame.tobytes())
                if n%15==0:
                    print(f'Frame {n+1}/{total}; elapsed {time.perf_counter()-started:.0f}s',flush=True)
        finally:
            for process in encoders.values():
                if process.stdin:
                    process.stdin.close()
            codes=[process.wait() for process in encoders.values()]
        if any(codes):
            raise RuntimeError(f'Video encoding failed: {codes}')
    metadata={'width':renderer.width,'height':renderer.height,'fps':FPS,'frames':total,
              'duration_seconds':DURATION,'supersample':args.supersample,
              'inclination_degrees':[78,87],'animation':'periodic static-viewpoint scan; no fluid evolution',
              'disk':DEFAULT_PARAMETERS.metadata(),'raytracer':renderer.tracer.metadata,
              'spectrum':renderer.spectrum.metadata,
              'display':{'bolometric_white_w_m2_sr':BOLO_WHITE,'thermal_white_cie_weighted':THERMAL_WHITE,
                         'encoding':'sRGB transfer, BT.709 YCbCr matrix; fixed global exposure; no additional effects',
                         'subpixel_average':'linear physical radiance, box filter'},
              'last_frame_stats':renderer.last_stats,'elapsed_seconds':time.perf_counter()-started,
              'outputs':{m:names[m] for m in modes}}
    (args.output/'render-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Finished in {time.perf_counter()-started:.1f}s.',flush=True)


if __name__=='__main__':
    main()
