#!/usr/bin/env python3
"""Verify the physical model, seamless loop, and both default Reel exports.

Run after render.py. --skip-physics repeats the movie checks without rerunning
the independent emission and ray-tracing diagnostics. No network is used.
"""
from __future__ import annotations

import argparse
import gc
import json
import math
import shutil
import subprocess
import tempfile
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

from render import DURATION, FPS, ROOT, Renderer, inclination


VIDEO_NAMES = {
    'bolometric': 'atlastro-black-hole-reel.mp4',
    'thermal': 'atlastro-visible-thermal-reel.mp4',
}


def check(record, name, passed, **details):
    record.setdefault('checks', []).append({'name': name, 'passed': bool(passed), **details})


def mp4_atoms(path):
    """Read top-level MP4 atom offsets without loading the movie into memory."""
    atoms = {}
    length = path.stat().st_size
    with path.open('rb') as stream:
        while stream.tell() + 8 <= length:
            position = stream.tell()
            header = stream.read(8)
            size = int.from_bytes(header[:4], 'big')
            kind = header[4:8].decode('ascii', errors='replace')
            header_size = 8
            if size == 1:
                size = int.from_bytes(stream.read(8), 'big')
                header_size = 16
            elif size == 0:
                size = length - position
            if size < header_size or position + size > length:
                raise ValueError(f'Invalid MP4 atom {kind!r} at byte {position}')
            atoms.setdefault(kind, position)
            stream.seek(position + size)
    return atoms


def probe_video(path, ffprobe):
    result = subprocess.run(
        [ffprobe, '-v', 'error', '-show_streams', '-show_format',
         '-of', 'json', str(path)], capture_output=True, text=True, check=True)
    probe = json.loads(result.stdout)
    streams = [stream for stream in probe['streams'] if stream['codec_type'] == 'video']
    if len(streams) != 1:
        raise ValueError(f'Expected one video stream, found {len(streams)}')
    return probe, streams[0]


def mean_difference(first, second):
    return float(np.abs(np.subtract(first, second, dtype=np.int16)).mean())


def decode_statistics(path, ffmpeg, width, height, hero_region=(380, 1210)):
    """Decode every frame; retain only endpoints and transition measurements.

    Differences are mean absolute errors in decoded 8-bit RGB code values.
    The hero region excludes the fixed title and most explanatory typography.
    """
    frame_bytes = width * height * 3
    first = previous = last = None
    transitions, hero_transitions = [], []
    count = 0
    hero = slice(*(round(edge * height / 1920) for edge in hero_region))
    command = [ffmpeg, '-hide_banner', '-loglevel', 'error', '-xerror',
               '-i', str(path), '-map', '0:v:0', '-an', '-sn',
               '-fps_mode', 'passthrough', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-']
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        try:
            while True:
                data = process.stdout.read(frame_bytes)
                if not data:
                    break
                if len(data) != frame_bytes:
                    raise ValueError(f'Decoder returned an incomplete frame: {len(data)} bytes')
                frame = np.frombuffer(data, dtype=np.uint8).reshape(height, width, 3)
                if first is None:
                    first = frame
                if previous is not None:
                    transitions.append(mean_difference(previous, frame))
                    hero_transitions.append(mean_difference(previous[hero], frame[hero]))
                previous = last = frame
                count += 1
            code = process.wait()
        except BaseException:
            process.kill()
            process.wait()
            raise
        finally:
            process.stdout.close()
        errors.seek(0)
        diagnostics = errors.read().decode('utf-8', errors='replace').strip()
    if code != 0:
        raise RuntimeError(f'FFmpeg decode failed ({code}): {diagnostics}')
    if count < 3:
        raise ValueError(f'Only {count} movie frames decoded')

    # Compare the wrap with its adjacent transitions: a moving viewpoint need
    # not have the same apparent image speed throughout the loop.
    endpoint_steps = transitions[:3] + transitions[-3:]
    endpoint_hero_steps = hero_transitions[:3] + hero_transitions[-3:]
    seam = mean_difference(last, first)
    hero_seam = mean_difference(last[hero], first[hero])
    # Lossy interframe encoding gives the first I-frame and final P-frame
    # different errors. A two-code-value MAE allowance avoids requiring exact
    # decoded pixels. The measured values remain in the report for review.
    tolerance = 2.0
    seam_limit = max(endpoint_steps) + tolerance
    hero_limit = max(endpoint_hero_steps) + tolerance
    return {
        'decoded_frames': count,
        'decoder_errors': diagnostics,
        'difference_units': 'mean absolute difference, RGB code values on a 0–255 scale',
        'all_transition_mae': {
            'minimum': float(min(transitions)),
            'median': float(np.median(transitions)),
            'maximum': float(max(transitions)),
        },
        'first_three_transition_mae': transitions[:3],
        'last_three_transition_mae': transitions[-3:],
        'seam_mae': seam,
        'hero_region_y': [hero.start, hero.stop],
        'hero_first_three_transition_mae': hero_transitions[:3],
        'hero_last_three_transition_mae': hero_transitions[-3:],
        'hero_seam_mae': hero_seam,
        'codec_mae_allowance': tolerance,
        'seam_limit_from_endpoint_transitions': seam_limit,
        'hero_seam_limit_from_endpoint_transitions': hero_limit,
        'seam_consistent_with_nearby_transitions': seam <= seam_limit and hero_seam <= hero_limit,
        'seam_check_scope': 'A perceptual discontinuity screen with codec tolerance; exact periodicity is checked separately on uncompressed renderer frames.',
    }


def verify_video(path, ffmpeg, ffprobe, hero_region=(380, 1210),
                 expected_transfer='iec61966-2-1', expected_range=None):
    result = {'file': str(path), 'size_bytes': path.stat().st_size}
    probe, stream = probe_video(path, ffprobe)
    fields = ('codec_name', 'profile', 'width', 'height', 'pix_fmt',
              'r_frame_rate', 'avg_frame_rate', 'nb_frames', 'duration',
              'start_time', 'color_space', 'color_transfer', 'color_primaries',
              'color_range', 'level', 'bit_rate')
    result['stream'] = {key: stream.get(key) for key in fields}
    check(result, '1080x1920_reel', (stream['width'], stream['height']) == (1080, 1920))
    check(result, 'h264_yuv420p', stream['codec_name'] == 'h264' and stream['pix_fmt'] == 'yuv420p')
    check(result, '30_fps', Fraction(stream['r_frame_rate']) == FPS and Fraction(stream['avg_frame_rate']) == FPS)
    check(result, '240_declared_frames', int(stream.get('nb_frames', '-1')) == FPS * DURATION)
    check(result, 'eight_second_duration', math.isclose(float(stream.get('duration', -1)), DURATION, abs_tol=1e-6))
    check(result, 'zero_start_time', math.isclose(float(stream.get('start_time', -1)), 0, abs_tol=1e-9))
    check(result, 'no_audio_stream', not any(s['codec_type'] == 'audio' for s in probe['streams']))
    check(result, 'specified_color_tags',
          stream.get('color_space') == 'bt709' and stream.get('color_primaries') == 'bt709'
          and stream.get('color_transfer') == expected_transfer,
          expected_transfer=expected_transfer)
    if expected_range is not None:
        check(result, 'specified_color_range', stream.get('color_range') == expected_range,
              expected_range=expected_range)
    atoms = mp4_atoms(path)
    result['mp4_atom_offsets'] = atoms
    check(result, 'fast_start', 'moov' in atoms and 'mdat' in atoms and atoms['moov'] < atoms['mdat'])
    result['decode'] = decode_statistics(path, ffmpeg, stream['width'], stream['height'], hero_region)
    check(result, 'all_240_frames_decode', result['decode']['decoded_frames'] == FPS * DURATION)
    check(result, 'loop_seam_with_codec_tolerance', result['decode']['seam_consistent_with_nearby_transitions'])
    result['passed'] = all(item['passed'] for item in result['checks'])
    return result


def verify_periodicity():
    result = {}
    total = FPS * DURATION
    angles = [inclination(frame / total) for frame in range(total)]
    mirrored_error = max(abs(inclination(frame / total) - inclination((total-frame) / total))
                         for frame in range(total))
    result['default_animation'] = {
        'fps': FPS, 'frames': total, 'duration_seconds': DURATION,
        'inclination_min_deg': min(angles), 'inclination_max_deg': max(angles),
        'phase_zero_deg': inclination(0), 'phase_one_deg': inclination(1),
        'maximum_mirrored_phase_angle_error_deg': mirrored_error,
    }
    check(result, 'default_phase_math_periodic', inclination(0) == inclination(1) and mirrored_error < 1e-12)
    renderer = Renderer(width=270, supersample=1)
    modes = list(VIDEO_NAMES)
    beginning = renderer.frames(0, modes)
    end = renderer.frames(1, modes)
    for mode in modes:
        equal = np.array_equal(np.asarray(beginning[mode]), np.asarray(end[mode]))
        check(result, f'{mode}_uncompressed_phase_zero_equals_one', equal,
              dimensions=[renderer.width, renderer.height])
    result['scope'] = 'Actual phase0/1 frame equality at 270x480; all default 240-frame camera phases checked at their full-precision values. Movie encoding is checked separately.'
    result['passed'] = all(item['passed'] for item in result['checks'])
    del renderer
    gc.collect()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'output',
                        help='directory containing both default exports; also receives verification.json')
    parser.add_argument('--skip-physics', action='store_true',
                        help='skip the independent emission and ray-tracing diagnostics')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {'physics_skipped': args.skip_physics, 'videos': {}, 'errors': []}
    ffmpeg, ffprobe = shutil.which('ffmpeg'), shutil.which('ffprobe')
    if not ffmpeg or not ffprobe:
        report['errors'].append('Both ffmpeg and ffprobe must be available on PATH.')
    missing = [name for name in VIDEO_NAMES.values() if not (args.output / name).is_file()]
    if missing:
        report['errors'].append('Render both default movies first. Missing: ' + ', '.join(missing))
    if not report['errors']:
        if not args.skip_physics:
            print('Checking emission and light propagation.', flush=True)
            try:
                from emission import diagnostics as emission_diagnostics
                from raytrace import diagnostics as raytrace_diagnostics
                report['physics'] = {'emission': emission_diagnostics(), 'raytrace': raytrace_diagnostics(), 'passed': True}
            except Exception as error:
                report['physics'] = {'passed': False, 'error': str(error)}
        print('Checking uncompressed periodicity.', flush=True)
        try:
            report['periodicity'] = verify_periodicity()
        except Exception as error:
            report['periodicity'] = {'passed': False, 'error': str(error)}
        for mode, name in VIDEO_NAMES.items():
            print(f'Decoding and checking {name}.', flush=True)
            try:
                report['videos'][mode] = verify_video(args.output / name, ffmpeg, ffprobe)
            except Exception as error:
                report['videos'][mode] = {'passed': False, 'error': str(error)}
    report['passed'] = (
        not report['errors']
        and (args.skip_physics or report.get('physics', {}).get('passed', False))
        and report.get('periodicity', {}).get('passed', False)
        and len(report['videos']) == len(VIDEO_NAMES)
        and all(video['passed'] for video in report['videos'].values())
    )
    report['elapsed_seconds'] = time.perf_counter() - started
    destination = args.output / 'verification.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(f'{"PASS" if report["passed"] else "FAIL"}: {destination}', flush=True)
    if not report['passed']:
        for error in report['errors']:
            print(error)
        for mode, result in report['videos'].items():
            if 'error' in result:
                print(f'{mode}: {result["error"]}')
            for item in result.get('checks', []):
                if not item['passed']:
                    print(f'{mode}: failed {item["name"]}')
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
