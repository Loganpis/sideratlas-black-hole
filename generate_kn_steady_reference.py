#!/usr/bin/env python3
"""Generate the v1.2 steady-emission comparison with complete provenance."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from render_kn_accretion import (KerrNewmanAccretionRenderer, OUTPUT_DIR,
                                  source_hash)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=OUTPUT_DIR / 'steady-reference.png')
    parser.add_argument('--manifest', type=Path,
                        help='defaults to steady-reference-manifest.json beside the image')
    args = parser.parse_args()
    destination = args.output.resolve()
    manifest_path = (args.manifest.resolve() if args.manifest else
                     destination.with_name('steady-reference-manifest.json'))
    destination.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    script = Path(__file__).resolve()
    source_render_sha = source_hash()
    generator_sha = sha256(script)
    renderer = KerrNewmanAccretionRenderer(
        1080, 2.5, refinement=2, temporal_samples=3)
    image = renderer.frames(0, structured=False)['color']
    image.save(destination)
    if source_hash() != source_render_sha or sha256(script) != generator_sha:
        raise RuntimeError('source changed while the steady reference was being generated')
    manifest = {
        'file': destination.name,
        'sha256': sha256(destination),
        'source_render_sha256': source_render_sha,
        'generator': script.name,
        'generator_sha256': generator_sha,
        'method': ('phase 0 production Kerr-Newman renderer with '
                   'structured=False; unchanged geodesics, transfer and compositor'),
        'dimensions': [image.width, image.height],
        'renderer': {
            'width': 1080,
            'supersample': 2.5,
            'refinement': 2,
            'temporal_samples': 3,
        },
    }
    temporary = manifest_path.with_suffix(manifest_path.suffix + '.partial')
    temporary.write_text(json.dumps(manifest, indent=2) + '\n')
    temporary.replace(manifest_path)
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    main()
