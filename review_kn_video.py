#!/usr/bin/env python3
"""Color-managed decoded-frame review of the delivered Kerr-Newman video.

Compare full-portrait decoded samples to their native float-radiance
compositions; this measures the local export, not Instagram transcoding.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess

import numpy as np
from PIL import Image,ImageDraw

from instagram_export import BT709_TO_SRGB_FILTER
from render_kn_orbit import KerrNewmanRenderer,OUTPUT_DIR,VIDEO_NAMES,source_hash


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR)
    parser.add_argument('--radiance-cache',type=Path,
                        help='native float-frame cache; defaults to this render signature')
    args=parser.parse_args()
    out=args.output
    meta=json.loads((out/'render-metadata.json').read_text())
    cache=args.radiance_cache or out/'.render-cache'/meta['radiance_signature'][:20]
    manifest=json.loads((cache/'radiance-cache-manifest.json').read_text())
    if manifest['radiance_signature']!=meta['radiance_signature']:
        raise ValueError('review cache signature does not match the delivered render')
    if meta['source_sha256']!=source_hash():
        raise ValueError('delivered video was rendered with different current source files')
    # Build only the compositor; every physical radiance frame already exists.
    renderer=KerrNewmanRenderer.__new__(KerrNewmanRenderer)
    renderer.width,renderer.height=meta['width'],meta['height']
    renderer.scale=renderer.width/1080
    renderer.spin=meta['spin_configuration']['a_star']
    renderer.charge=meta['spin_configuration']['charge_q']
    renderer.fov=meta['geometry']['fov_horizontal_rg']
    renderer._configure_layout()
    renderer.overlays={mode:renderer._overlay(mode) for mode in VIDEO_NAMES}
    overlay=np.asarray(renderer.overlays['color'])[...,3]>0
    # Exclude the entire title band, including any codec bleed around glyphs.
    overlay_box=renderer.overlays['color'].getbbox()
    header_free=np.ones((renderer.height,renderer.width),bool)
    header_free[:overlay_box[3]+8]=False
    indices=list(range(0,240,20))
    thumbnail_w,thumbnail_h,cell_h=270,480,520
    sheet=Image.new('RGB',(thumbnail_w*4,cell_h*3))
    draw=ImageDraw.Draw(sheet)
    records=[]
    for j,n in enumerate(indices):
        raw=subprocess.run([shutil.which('ffmpeg') or 'ffmpeg','-v','error',
            '-i',str(out/VIDEO_NAMES['color']),
            '-vf',f'select=eq(n\\,{n}),'+BT709_TO_SRGB_FILTER,
            '-frames:v','1','-f','rawvideo','-pix_fmt','rgb24','-'],
            capture_output=True,check=True).stdout
        decoded=np.frombuffer(raw,np.uint8).reshape(renderer.height,renderer.width,3)
        image=Image.fromarray(decoded)
        if n==60:
            image.save(out/'decoded-review.png')
        sheet.paste(image.resize((thumbnail_w,thumbnail_h),Image.Resampling.LANCZOS),
                    ((j%4)*thumbnail_w,(j//4)*cell_h))
        draw.text(((j%4)*thumbnail_w+14,(j//4)*cell_h+490),f'{n/30:.2f} s',
                  font=renderer.font(22),fill=(215,215,215))
        path=cache/f'{n:03d}.npy'
        frame_sha=hashlib.sha256(path.read_bytes()).hexdigest()
        if frame_sha!=manifest['frames_sha256'][path.name]:
            raise ValueError(f'{path.name} does not match the native radiance manifest')
        radiance=np.load(path,allow_pickle=False)
        if radiance.shape!=(meta['geometry']['height'],meta['geometry']['width']):
            raise ValueError('review requires native detector samples, not an older downsampled cache')
        source=np.asarray(renderer.compose(radiance,'color'))
        difference=abs(source.astype(float)-decoded.astype(float))
        lit=(np.max(source,axis=2)>10)&header_free&~overlay
        if not lit.any():
            raise ValueError(f'frame {n} has no measurable lit simulation pixels')
        records.append({'frame':n,'seconds':n/30,
            'native_radiance_sha256':frame_sha,
            'full_rgb_mae_255':float(difference.mean()),
            'lit_disk_rgb_mae_255':float(difference[lit].mean()),
            'lit_disk_pixels_excluding_header':int(lit.sum())})
    sheet.save(out/'orbit-review.jpg',quality=95)
    report={
        'video_file':str(Path(out.name)/VIDEO_NAMES['color']),
        'video_sha256':hashlib.sha256((out/VIDEO_NAMES['color']).read_bytes()).hexdigest(),
        'source_sha256':source_hash(),
        'radiance_signature':meta['radiance_signature'],
        'scope':__doc__,
        'frames':records,
        'overlay_bbox':list(overlay_box),
        'bottom_clearance_pixels':renderer.height-overlay_box[3],
        'full_portrait_review':True,
        'footer_overlay_absent':not overlay[round(500*renderer.scale):].any().item(),
        'lit_error_region':'Pixels above brightness 10 below the title band; no assumption of a black header',
    }
    (out/'visual-review.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
