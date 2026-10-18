#!/usr/bin/env python3
"""Bundle the verified v1.2 release and historical sources without caches."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'output-sideratlas-kn-v120'
RELEASE_INPUTS=('caption-kerr-newman.text','caption-kerr-newman.txt',
    'kn-flow-diagnostics.json','kn-accretion-trace-diagnostics.json',
    'kn-dynamic-emission-diagnostics.json','kn-sampling-diagnostics.json',
    'kn-image-order-diagnostics.json','kn-emission-comparison-diagnostics.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    report=json.loads((OUT/'verification.json').read_text())
    if report.get('passed') is not True:raise RuntimeError('Current delivery verification must pass')
    if report.get('verifier_sha256')!=sha(ROOT/'verify_kn_accretion.py'):
        raise RuntimeError('Delivery verification was produced by a different verifier')
    meta=json.loads((OUT/'render-metadata.json').read_text())
    from render_kn_accretion import source_hash
    if meta['source_sha256']!=source_hash():raise RuntimeError('Current source must match the delivered render')
    assembly=meta.get('assembly',{})
    if (assembly.get('script')!='assemble_kn_accretion.py' or
        assembly.get('script_sha256')!=sha(ROOT/'assemble_kn_accretion.py') or
        assembly.get('radiance_unchanged') is not True):
        raise RuntimeError('Cached assembly provenance must match the published utility')
    expected_release_inputs={name:sha(ROOT/name) for name in RELEASE_INPUTS}
    if report.get('release_input_sha256')!=expected_release_inputs:
        raise RuntimeError('Caption or numerical diagnostics changed after verification')
    if list(OUT.glob('*.partial.mp4')):raise RuntimeError('Partial videos cannot be packaged')
    assembly_required=('sideratlas-kerr-newman-reel.mp4','sideratlas-kerr-newman-intensity.mp4',
        'sideratlas-kerr-newman-reel-loop.mp4','sideratlas-kerr-newman-intensity-loop.mp4',
        'sideratlas-kerr-newman-color-cover.png','sideratlas-kerr-newman-mono-cover.png',
        'contact-sheet.jpg','steady-reference.png','steady-reference-manifest.json',
        'native-radiance-manifest.json',
        'view-diagnostics.json','render-metadata.json')
    required=assembly_required+('verification.json','visual-review.json','orbit-review.jpg',
        'decoded-review.png')
    for name in required:
        if not (OUT/name).is_file():raise FileNotFoundError(OUT/name)
    verified_assets=report.get('asset_sha256',{})
    for name in assembly_required:
        if verified_assets.get(name)!=sha(OUT/name):
            raise RuntimeError('Delivery verification is stale for '+name)
    for mode,name in (('color',assembly_required[0]),('mono',assembly_required[1])):
        if report['videos'][mode]['sha256']!=sha(OUT/name):
            raise RuntimeError('Delivery verification is stale for '+name)
    steady=json.loads((OUT/'steady-reference-manifest.json').read_text())
    recorded_steady=meta.get('steady_reference',{})
    if (steady.get('file')!='steady-reference.png' or
        steady.get('sha256')!=sha(OUT/'steady-reference.png') or
        steady.get('source_render_sha256')!=source_hash() or
        steady.get('generator')!='generate_kn_steady_reference.py' or
        steady.get('generator_sha256')!=sha(ROOT/'generate_kn_steady_reference.py') or
        steady.get('dimensions')!=[meta.get('width'),meta.get('height')] or
        recorded_steady.get('manifest')!='steady-reference-manifest.json' or
        recorded_steady.get('manifest_sha256')!=sha(OUT/'steady-reference-manifest.json') or
        any(recorded_steady.get(key)!=value for key,value in steady.items())):
        raise RuntimeError('Steady reference must match its current source and generator')
    review=json.loads((OUT/'visual-review.json').read_text())
    if (review['source_sha256']!=meta['source_sha256'] or
        review['radiance_signature']!=meta['radiance_signature'] or
        review['video_sha256']!=sha(OUT/assembly_required[0]) or
        review.get('reviewer')!='review_kn_video.py' or
        review.get('reviewer_sha256')!=sha(ROOT/'review_kn_video.py') or
        review.get('entrypoint')!='review_kn_accretion_video.py' or
        review.get('entrypoint_sha256')!=sha(ROOT/'review_kn_accretion_video.py') or
        any(review.get('review_artifact_sha256',{}).get(name)!=sha(OUT/name)
            for name in ('orbit-review.jpg','decoded-review.png'))):
        raise RuntimeError('Decoded-frame review must match the current render and video')
    files=sorted(p for p in ROOT.iterdir() if p.is_file() and p.suffix in ('.py','.cpp','.md','.txt','.text','.json'))
    files.append(ROOT/'.gitignore')
    excluded={'render-process.json','checks-process.json','verification-preflight.json',
              'delivery-manifest.json','release-files.json'}
    for folder in sorted(ROOT.glob('output*')):
        if not folder.is_dir():continue
        files+=sorted(p for p in folder.iterdir() if p.is_file() and
            p.suffix in ('.mp4','.png','.jpg','.json') and p.name not in excluded and '.partial' not in p.name)
        provenance=folder/'provenance'
        if provenance.is_dir():files+=sorted(p for p in provenance.iterdir() if p.is_file() and p.suffix in ('.json','.md','.txt'))
    for folder in ('fonts','data'):
        files+=sorted(p for p in (ROOT/folder).iterdir() if p.is_file())
    manifest={'version':'1.2.0','release':'Sideratlas dynamic neutral Kerr-Newman accretion',
        'scope':'All published source and completed media; illustrative source choices disclosed; caches excluded',
        'files':{str(p.relative_to(ROOT)):{'sha256':sha(p),'size_bytes':p.stat().st_size} for p in files}}
    (OUT/'release-files.json').write_text(json.dumps(manifest,indent=2)+'\n')
    files.append(OUT/'release-files.json')
    archive=ROOT/'sideratlas-kerr-newman-v1.2.0-source.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as bundle:
        for path in files:bundle.write(path,'sideratlas-black-hole/'+str(path.relative_to(ROOT)))
    with zipfile.ZipFile(archive) as bundle:
        bad=bundle.testzip()
        if bad:raise RuntimeError('Archive CRC failure:'+bad)
        members=len(bundle.namelist())
    delivery={'version':'1.2.0','archive':{'path':archive.name,'sha256':sha(archive),
        'size_bytes':archive.stat().st_size,'members':members},
        'videos':{name:{'path':name,'sha256':sha(OUT/name),'size_bytes':(OUT/name).stat().st_size} for name in assembly_required[:4]},
        'description':{'path':'caption-kerr-newman.text','sha256':sha(ROOT/'caption-kerr-newman.text')},
        'verification_sha256':sha(OUT/'verification.json'),'metadata_sha256':sha(OUT/'render-metadata.json')}
    (OUT/'delivery-manifest.json').write_text(json.dumps(delivery,indent=2)+'\n')
    print(json.dumps(delivery,indent=2))


if __name__=='__main__':main()
