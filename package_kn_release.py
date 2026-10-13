#!/usr/bin/env python3
"""Package the verified public project, excluding large raw-frame caches."""
from pathlib import Path
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'output-sideratlas-kn'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
report=json.loads((OUT/'verification.json').read_text())
if report.get('passed') is not True:raise RuntimeError('Final verification must pass before packaging')
output_names='''sideratlas-kerr-newman-reel.mp4 sideratlas-kerr-newman-intensity.mp4 sideratlas-kerr-newman-reel-loop.mp4 sideratlas-kerr-newman-intensity-loop.mp4 sideratlas-kerr-newman-color-cover.png sideratlas-kerr-newman-mono-cover.png steady-reference.png contact-sheet.jpg render-metadata.json verification.json decoded-review.png orbit-review.jpg visual-review.json'''.split()
for name in output_names:
    if not (OUT/name).is_file():raise FileNotFoundError(OUT/name)
files=sorted(p for p in ROOT.iterdir() if p.is_file() and p.suffix in ('.py','.cpp','.md','.txt','.text','.json'))
files.append(ROOT/'.gitignore')
excluded={'render-process.json','checks-process.json','delivery-manifest.json','release-files.json'}
for directory in ('output','output-kerr','output-optical','output-schwarzschild','output-sideratlas','output-sideratlas-kn'):
    folder=ROOT/directory
    files+=sorted(p for p in folder.iterdir() if p.is_file() and p.suffix in ('.mp4','.png','.jpg','.json') and p.name not in excluded and '.partial' not in p.name)
    provenance=folder/'provenance'
    if provenance.is_dir():files+=sorted(p for p in provenance.iterdir() if p.is_file() and p.suffix in ('.json','.md','.txt'))
files+=sorted(p for p in (ROOT/'fonts').iterdir() if p.is_file())
files+=sorted(p for p in (ROOT/'data').iterdir() if p.is_file())
for path in files:
    if not path.is_file():raise FileNotFoundError(path)
manifest={'release':'@sideratlas Kerr-Newman a*=.9 q=.4; four identical eight-second cycles',
          'scope':'Complete published source and completed media, including historical versions; generated caches excluded',
          'files':{str(p.relative_to(ROOT)):{'sha256':sha(p),'size_bytes':p.stat().st_size} for p in files}}
(OUT/'release-files.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(OUT/'release-files.json')
archive=ROOT/'sideratlas-kerr-newman-source.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path in files:z.write(path,'sideratlas-kerr-newman/'+str(path.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
    bad=z.testzip()
    if bad:raise RuntimeError('Archive CRC failure: '+bad)
    members=len(z.namelist())
delivery={'archive':{'path':archive.name,'sha256':sha(archive),'size_bytes':archive.stat().st_size,'members':members},
          'main_video':{'path':output_names[0],'sha256':sha(OUT/output_names[0])},
          'verification_sha256':sha(OUT/'verification.json'),
          'metadata_sha256':sha(OUT/'render-metadata.json')}
(OUT/'delivery-manifest.json').write_text(json.dumps(delivery,indent=2)+'\n')
print(json.dumps(delivery,indent=2))
