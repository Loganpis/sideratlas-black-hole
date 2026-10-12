#!/usr/bin/env python3
"""Color-managed decoded-frame review of the delivered Kerr-Newman video."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
from PIL import Image,ImageDraw
from instagram_export import BT709_TO_SRGB_FILTER
from render_kn_orbit import KerrNewmanRenderer,OUTPUT_DIR,VIDEO_NAMES,source_hash,display_intensity,illustrative_palette
root=Path(__file__).resolve().parent
out=OUTPUT_DIR
meta=json.loads((out/'render-metadata.json').read_text())
cache=out/'.render-cache'/meta['radiance_signature'][:20]
# Build only the compositor; the audited raw radiance is already cached.
r=KerrNewmanRenderer.__new__(KerrNewmanRenderer)
r.width,r.height,r.scale,r.spin,r.charge=1080,1920,1,.9,.4
r.hero_h,r.hero_y,r.hero_x=900,425,108
r.hero_display_size=(864,720)
r.overlays={mode:r._overlay(mode) for mode in VIDEO_NAMES}
indices=list(range(0,240,20))
sheet=Image.new('RGB',(432*4,430*3));draw=ImageDraw.Draw(sheet)
records=[]
for j,n in enumerate(indices):
 raw=subprocess.run(['ffmpeg','-v','error','-i',str(out/VIDEO_NAMES['color']),'-vf',f'select=eq(n\\,{n}),'+BT709_TO_SRGB_FILTER,'-frames:v','1','-f','rawvideo','-pix_fmt','rgb24','-'],capture_output=True,check=True).stdout
 decoded=np.frombuffer(raw,np.uint8).reshape(1920,1080,3)
 image=Image.fromarray(decoded)
 if n==60:image.save(out/'decoded-review.png')
 crop=image.crop((108,425,972,1145)).resize((432,360),Image.Resampling.LANCZOS)
 sheet.paste(crop,((j%4)*432,(j//4)*430))
 draw.text(((j%4)*432+20,(j//4)*430+373),f'{n/30:.2f} s',font=r.font(22),fill=(215,215,215))
 source=np.asarray(r.compose(np.load(cache/f'{n:03d}.npy',allow_pickle=False),'color'))
 diff=abs(source.astype(float)-decoded.astype(float)); hero=diff[425:1145,108:972]
 lit=np.max(source[425:1145,108:972],axis=2)>10
 records.append({'frame':n,'seconds':n/30,'full_rgb_mae_255':float(diff.mean()),'lit_disk_rgb_mae_255':float(hero[lit].mean()),'black_top_max':int(decoded[:200].max())})
sheet.save(out/'orbit-review.jpg',quality=95)
report={'video_sha256':hashlib.sha256((out/VIDEO_NAMES['color']).read_bytes()).hexdigest(),
 'source_sha256':source_hash(),'scope':'Color-managed decoded samples versus the corresponding float-radiance compositions; not an Instagram-upload measurement',
 'frames':records,'overlay_bbox':list(r.overlays['color'].getbbox()),'bottom_clearance_pixels':1920-r.overlays['color'].getbbox()[3],
 'exact_black_header':all(d['black_top_max']==0 for d in records)}
(out/'visual-review.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
