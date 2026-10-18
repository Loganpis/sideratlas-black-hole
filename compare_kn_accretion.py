#!/usr/bin/env python3
"""Compare emissivity prescriptions on identical v1.2 photon paths.

This is a source-isolation diagnostic, not a comparison of complete historical
models: all columns share the current flow, opacity, sampling and compositor.
The smooth experiment is reconstructed by default; an optional saved original
source file can replace that reference.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

from PIL import Image, ImageDraw

from kn_orbit_emission import ShearingDiskEmission
from compare_kn_emission import SupersededSmoothEmission
from render_kn_accretion import (KerrNewmanAccretionRenderer, PreparedRayEmission,
                                ROOT, source_hash)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=540)
    parser.add_argument('--supersample',type=float,default=2.5)
    parser.add_argument('--smooth-reference',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'output-sideratlas-kn-v120'/'preview'/'texture-comparison')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    renderer=KerrNewmanAccretionRenderer(args.width,args.supersample)
    sources=[('v1.1 emissivity',ShearingDiskEmission())]
    reference_sha=None
    smooth_source=SupersededSmoothEmission(flow=renderer.flow)
    if args.smooth_reference:
        spec=importlib.util.spec_from_file_location('_superseded_source',args.smooth_reference)
        module=importlib.util.module_from_spec(spec)
        sys.modules[spec.name]=module
        spec.loader.exec_module(module)
        smooth_source=module.DynamicDiskEmission(flow=renderer.flow)
        reference_sha=hashlib.sha256(args.smooth_reference.read_bytes()).hexdigest()
    sources.append(('Superseded smooth',smooth_source))
    sources.append(('Restored multiscale',renderer.source))
    phases=(0.,.2,.4)
    cell_w,cell_h,label_h=360,640,44
    sheet=Image.new('RGB',(cell_w*len(sources),(cell_h+label_h)*len(phases)))
    draw=ImageDraw.Draw(sheet)
    original_source=renderer.source
    records=[]
    for row,phase in enumerate(phases):
        renderer.set_view(phase)
        coarse_events=renderer.coarse.events
        fine_events=renderer.fine.events if renderer.fine is not None else None
        for col,(label,source) in enumerate(sources):
            renderer.coarse=PreparedRayEmission(coarse_events,source,renderer.flow,
                                               renderer.tracer.width*renderer.tracer.height)
            renderer.fine=(PreparedRayEmission(fine_events,source,renderer.flow,
                len(renderer.refined_pixels)*renderer.refinement**2) if fine_events is not None else None)
            radiance=renderer._exposed_radiance(phase)
            frame=renderer.compose(radiance,'color')
            name=f'phase-{phase:.2f}-source-{col}.png'
            frame.save(args.output/name)
            sheet.paste(frame.resize((cell_w,cell_h),Image.Resampling.LANCZOS),
                        (col*cell_w,row*(cell_h+label_h)+label_h))
            draw.text((col*cell_w+10,row*(cell_h+label_h)+12),f'{label} / phase {phase:.1f}',
                      font=renderer.font(28),fill=(225,225,225))
            records.append({'source':label,'phase':phase,'file':name,
                            'total_native_radiance':float(radiance.sum())})
        renderer.source=original_source
    sheet.save(args.output/'emission-comparison.jpg',quality=96)
    report={'scope':__doc__,'source_sha256':source_hash(),
            'comparison_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'original_emissivity_sha256':hashlib.sha256((ROOT/'kn_orbit_emission.py').read_bytes()).hexdigest(),
            'optional_smooth_reference_sha256':reference_sha,
            'smooth_reference_reconstruction_sha256':hashlib.sha256((ROOT/'compare_kn_emission.py').read_bytes()).hexdigest(),
            'width':renderer.width,'height':renderer.height,'supersample':args.supersample,
            'fixed_across_columns':['null paths','emitter four-velocity','per-event delays','gray opacity',
                'independent selected subpixel rays','source shutter','exposure','tone curve','portrait composition'],
            'caveat':'v1.1 emissivity has no inner source; this column is not the complete historical v1.1 renderer',
            'frames':records}
    (args.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'contact_sheet':str(args.output/'emission-comparison.jpg'),'frames':len(records)}))


if __name__=='__main__':main()
