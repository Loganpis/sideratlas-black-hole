"""Disposable composition study; does not change the production renderer."""
from pathlib import Path
import time
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from kerr_raytrace import KerrTracer
from kerr_emission import KerrDiskParameters, KerrEmission
from emission import ThermalSpectrum
from optics import atmosphere_rgb, chromaticity_preserving_tonemap
from render import encode_srgb

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'output-optical-study'
OUT.mkdir(exist_ok=True)
W,H=480,600
disk=KerrEmission(KerrDiskParameters(mass_solar=1e10,eddington_ratio=.05,spin_a_star=.9))
spectrum=ThermalSpectrum()
font=ImageFont.truetype(str(ROOT/'fonts/lmroman10-regular.otf'),20)
sheets={k:Image.new('RGB',(3*W,3*(H+44))) for k in ('current','linear','contrast')}
started=time.perf_counter()
for col,fov in enumerate((28,36,44)):
    tracer=KerrTracer(W,H,fov=fov,roll_deg=14,spin=.9,rout=500,angle_step=.006)
    for row,inc in enumerate((70,80,87)):
        m=tracer.trace(inc)
        rgb=atmosphere_rgb(disk.temperature(m['r']),m['g'],m['mu'],spectrum)
        np.savez_compressed(OUT/f'physical-fov{fov}-i{inc}.npz',rgb=rgb.astype('float32'))
        values=np.maximum(rgb,0)/3e8
        maxval=values.max(-1,keepdims=True)
        # Shared-RGB scaling keeps chromaticity. This is solely a display curve.
        linear=values/np.maximum(maxval,1)
        contrast=values*np.divide(maxval,1+maxval**2,out=np.zeros_like(maxval),where=True)
        ims={'current':chromaticity_preserving_tonemap(rgb,1/1.8e8),
             'linear':linear,'contrast':contrast}
        for key,values in ims.items():
            img=Image.fromarray(encode_srgb(values))
            img.save(OUT/f'{key}-fov{fov}-i{inc}.png')
            sheet=sheets[key]
            sheet.paste(img,(col*W,row*(H+44)+44))
            ImageDraw.Draw(sheet).text((col*W+15,row*(H+44)+10),f'{key}: fov={fov}, i={inc}',font=font,fill='white')
        print(f'fov {fov} i {inc}: {time.perf_counter()-started:.1f}s; RGB max {rgb.max():.3g}',flush=True)
for key,sheet in sheets.items():
    sheet.save(OUT/f'{key}-contact-sheet.jpg',quality=94)
print(str(OUT),flush=True)
