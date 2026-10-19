"""Streaming video inspection and loop-boundary checks for Sideratlas."""
from pathlib import Path
import json
import subprocess
import tempfile
import numpy as np


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
