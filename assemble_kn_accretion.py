#!/usr/bin/env python3
"""Assemble the v1.2 delivery from verified native radiance frames.

This utility performs no ray tracing and changes no physical radiance. It
reuses the production compositor and Instagram encoder, records the exact
native-frame bytes, and builds the four-repeat 32-second delivery.
"""
from __future__ import annotations

import argparse
from io import BytesIO
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import numpy as np
from PIL import Image

from instagram_export import encoder
from render_kn_accretion import OUTPUT_DIR, source_hash
from render_kn_orbit import (KerrNewmanRenderer, FPS, REPEATS, VIDEO_NAMES,
                             LOOP_VIDEO_NAMES)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def compositor(metadata: dict) -> KerrNewmanRenderer:
    """Construct only the unchanged v1.1/v1.2 display compositor."""
    renderer = KerrNewmanRenderer.__new__(KerrNewmanRenderer)
    renderer.width = int(metadata['width'])
    renderer.height = int(metadata['height'])
    renderer.scale = renderer.width / 1080
    renderer.spin = float(metadata['spin_configuration']['a_star'])
    renderer.charge = float(metadata['spin_configuration']['charge_q'])
    renderer.fov = float(metadata['geometry']['fov_horizontal_rg'])
    renderer._configure_layout()
    renderer.overlays = {mode: renderer._overlay(mode) for mode in VIDEO_NAMES}
    return renderer


def load_native(path: Path, expected_shape: tuple[int, int]) -> tuple[np.ndarray, str]:
    """Read each frame once, hashing the exact bytes used for composition."""
    data = path.read_bytes()
    digest = sha256_bytes(data)
    values = np.load(BytesIO(data), allow_pickle=False)
    if values.dtype != np.dtype('float32') or values.shape != expected_shape:
        raise ValueError(f'{path.name}: expected float32 {expected_shape}, got {values.dtype} {values.shape}')
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f'{path.name}: radiance must be finite and nonnegative')
    return values, digest


def repeat_loop(loop_path: Path, delivery_path: Path) -> None:
    ffmpeg = shutil.which('ffmpeg') or 'ffmpeg'
    temporary = delivery_path.with_name(delivery_path.stem + '.partial.mp4')
    temporary.unlink(missing_ok=True)
    subprocess.run([
        ffmpeg, '-v', 'error', '-y', '-stream_loop', str(REPEATS - 1),
        '-i', str(loop_path), '-map', '0:v:0', '-an', '-c:v', 'copy',
        '-movflags', '+faststart', '-use_editlist', '0', str(temporary),
    ], check=True)
    temporary.replace(delivery_path)


def transactional_commit(staging: Path, output: Path, names: list[str]) -> None:
    """Replace the complete artifact set, restoring every old file on error."""
    backup = Path(tempfile.mkdtemp(prefix='.assembly-backup-', dir=output))
    backed_up: list[tuple[Path, Path]] = []
    installed: list[Path] = []
    try:
        for name in names:
            source = staging / name
            target = output / name
            if not source.is_file():
                raise FileNotFoundError(source)
            if target.exists():
                saved = backup / name
                target.replace(saved)
                backed_up.append((saved, target))
            source.replace(target)
            installed.append(target)
    except BaseException:
        for target in reversed(installed):
            target.unlink(missing_ok=True)
        for saved, target in reversed(backed_up):
            saved.replace(target)
        raise
    finally:
        shutil.rmtree(backup, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT_DIR)
    parser.add_argument('--steady-reference', type=Path,
                        default=OUTPUT_DIR / 'steady-reference.png')
    parser.add_argument('--steady-reference-manifest', type=Path,
                        help='producer manifest; defaults beside the steady reference')
    parser.add_argument('--verify-cover', action='store_true',
                        help='require existing covers to match cached frame 000 before replacement')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    metadata_path = output / 'render-metadata.json'
    metadata = json.loads(metadata_path.read_text())
    current_source = source_hash()
    if metadata.get('source_sha256') != current_source:
        raise RuntimeError('render metadata does not match the current frozen physical source')
    signature = metadata.get('radiance_signature')
    if not isinstance(signature, str) or len(signature) != 64:
        raise RuntimeError('render metadata has no full radiance signature')
    if (metadata.get('unique_frames'), metadata.get('fps'), metadata.get('repetitions')) != (240, FPS, REPEATS):
        raise RuntimeError('assembly requires the 240-frame, 30 fps, four-repeat production contract')

    cache = output / '.render-cache' / signature[:20]
    diagnostic_path = cache / 'view-diagnostics.json'
    diagnostics = json.loads(diagnostic_path.read_text())
    if (not isinstance(diagnostics, list) or len(diagnostics) != 120 or
            [entry.get('frame') for entry in diagnostics] != list(range(120))):
        raise RuntimeError('all 120 independent viewpoint diagnostics are required')
    expected_shape = (int(metadata['geometry']['height']), int(metadata['geometry']['width']))
    native_paths = [cache / f'{n:03d}.npy' for n in range(240)]
    missing = [path.name for path in native_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f'missing {len(missing)} native radiance frames: {missing[:8]}')

    trusted_manifest_path = cache / 'radiance-cache-manifest.json'
    trusted_manifest_bytes = trusted_manifest_path.read_bytes()
    trusted_manifest = json.loads(trusted_manifest_bytes)
    trusted_hashes = trusted_manifest.get('frames_sha256', {})
    expected_names = {path.name for path in native_paths}
    if (trusted_manifest.get('radiance_signature') != signature or
            trusted_manifest.get('source_render_sha256') != current_source or
            set(trusted_hashes) != expected_names):
        raise RuntimeError('a complete producer radiance manifest matching this render is required')

    renderer = compositor(metadata)
    retained: dict[int, np.ndarray] = {}
    frame_hashes: dict[str, str] = {}
    presentation_indices = {0, 60, 120}
    for n in sorted(presentation_indices):
        values, digest = load_native(native_paths[n], expected_shape)
        if digest != trusted_hashes[native_paths[n].name]:
            raise RuntimeError(f'{native_paths[n].name} differs from the producer manifest')
        retained[n] = values
        frame_hashes[native_paths[n].name] = digest
    composed = {n: {mode: renderer.compose(retained[n], mode) for mode in VIDEO_NAMES}
                for n in sorted(retained)}
    if args.verify_cover:
        for mode, label in (('color', 'color'), ('mono', 'mono')):
            existing = output / f'sideratlas-kerr-newman-{label}-cover.png'
            if not existing.is_file():
                raise FileNotFoundError(f'cover equivalence requested but {existing.name} is absent')
            if not np.array_equal(np.asarray(Image.open(existing).convert('RGB')),
                                  np.asarray(composed[0][mode])):
                raise RuntimeError(f'{existing.name} does not match cached frame 000')

    reference_source = args.steady_reference.resolve()
    if not reference_source.is_file():
        raise FileNotFoundError('a freshly generated steady-reference.png is required')
    if reference_source.name != 'steady-reference.png':
        raise ValueError('the steady reference must retain its published filename')
    with Image.open(reference_source) as reference:
        reference.load()
        if reference.size != (renderer.width, renderer.height):
            raise ValueError('steady-reference.png has the wrong delivery dimensions')
    reference_sha = sha256_file(reference_source)
    reference_manifest_source = (args.steady_reference_manifest.resolve()
        if args.steady_reference_manifest else
        reference_source.with_name('steady-reference-manifest.json'))
    if not reference_manifest_source.is_file():
        raise FileNotFoundError('the steady-reference producer manifest is required')
    reference_manifest_bytes = reference_manifest_source.read_bytes()
    reference_manifest = json.loads(reference_manifest_bytes)
    reference_generator = Path(__file__).resolve().parent / 'generate_kn_steady_reference.py'
    expected_reference = {
        'file': reference_source.name,
        'sha256': reference_sha,
        'source_render_sha256': current_source,
        'generator': reference_generator.name,
        'generator_sha256': sha256_file(reference_generator),
        'dimensions': [renderer.width, renderer.height],
    }
    if any(reference_manifest.get(key) != value
           for key, value in expected_reference.items()):
        raise RuntimeError('steady-reference image or manifest does not match its current producer')

    staging = Path(tempfile.mkdtemp(prefix='.assembly-staging-', dir=output))
    try:
        for mode, label in (('color', 'color'), ('mono', 'mono')):
            composed[0][mode].save(staging / f'sideratlas-kerr-newman-{label}-cover.png')
        sheet = Image.new('RGB', (1080, 1280))
        for column, n in enumerate((0, 60, 120)):
            for row, mode in enumerate(VIDEO_NAMES):
                sheet.paste(composed[n][mode].resize((360, 640), Image.Resampling.LANCZOS),
                            (column * 360, row * 640))
        sheet.save(staging / 'contact-sheet.jpg', quality=95)
        shutil.copy2(reference_source, staging / 'steady-reference.png')
        (staging / 'steady-reference-manifest.json').write_bytes(reference_manifest_bytes)
        (staging / 'native-radiance-manifest.json').write_bytes(trusted_manifest_bytes)

        processes = {}
        loop_partials = {}
        error = None
        try:
            for mode, name in LOOP_VIDEO_NAMES.items():
                partial = staging / (Path(name).stem + '.partial.mp4')
                loop_partials[mode] = partial
                processes[mode] = encoder(partial, renderer.width, renderer.height, FPS)
            for n, path in enumerate(native_paths):
                if n in retained:
                    values = retained[n]
                    digest = frame_hashes[path.name]
                else:
                    values, digest = load_native(path, expected_shape)
                    frame_hashes[path.name] = digest
                if digest != trusted_hashes[path.name]:
                    raise RuntimeError(f'{path.name} differs from the producer manifest')
                for mode, process in processes.items():
                    process.stdin.write(renderer.compose(values, mode).tobytes())
                if n % 20 == 0 or n == 239:
                    print(f'Composed {n + 1}/240 native frames', flush=True)
        except BaseException as exc:
            error = exc
        finally:
            for process in processes.values():
                if process.stdin is not None and not process.stdin.closed:
                    try:
                        process.stdin.close()
                    except BrokenPipeError:
                        pass
            codes = {mode: process.wait() for mode, process in processes.items()}
        if error is not None:
            raise error
        if any(code != 0 for code in codes.values()):
            raise RuntimeError(f'incomplete encodes: {codes}')
        for mode, name in LOOP_VIDEO_NAMES.items():
            loop_partials[mode].replace(staging / name)
        for mode, name in VIDEO_NAMES.items():
            repeat_loop(staging / LOOP_VIDEO_NAMES[mode], staging / name)

        (staging / 'view-diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
        script_sha = sha256_file(Path(__file__))
        metadata['assembly'] = {
            'script': Path(__file__).name,
            'script_sha256': script_sha,
            'method': 'producer-manifest-verified native radiance, unchanged compositor, single-generation loop encoding',
            'radiance_unchanged': True,
            'native_manifest': 'native-radiance-manifest.json',
            'native_manifest_sha256': sha256_bytes(trusted_manifest_bytes),
        }
        metadata['steady_reference'] = dict(reference_manifest)
        metadata['steady_reference'].update({
            'file': 'steady-reference.png',
            'manifest': 'steady-reference-manifest.json',
            'manifest_sha256': sha256_bytes(reference_manifest_bytes),
        })
        metadata['repetition_encoding'] = {
            'method': 'four bitstream-identical repetitions of the encoded 8-second loop',
            'stream_loop_count': REPEATS - 1,
            'video_reencoded': False,
            'audio': 'none',
        }
        metadata['elapsed_seconds'] = time.perf_counter() - started
        metadata['elapsed_scope'] = ('native-cache validation, composition and encoding only; '
                                     'excludes prior geodesic tracing and radiance synthesis')
        (staging / 'render-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
        names = list(LOOP_VIDEO_NAMES.values()) + list(VIDEO_NAMES.values()) + [
            'sideratlas-kerr-newman-color-cover.png',
            'sideratlas-kerr-newman-mono-cover.png', 'contact-sheet.jpg',
            'steady-reference.png', 'steady-reference-manifest.json',
            'native-radiance-manifest.json',
            'view-diagnostics.json', 'render-metadata.json']
        transactional_commit(staging, output, names)
        for stale in ('verification.json', 'visual-review.json', 'orbit-review.jpg',
                      'decoded-review.png'):
            (output / stale).unlink(missing_ok=True)
        for partial in output.glob('*.partial.mp4'):
            partial.unlink()
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    print(f'Assembly complete in {metadata["elapsed_seconds"]:.1f}s', flush=True)


if __name__ == '__main__':
    main()
