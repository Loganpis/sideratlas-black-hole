"""Audit local SDR conversion and H.264 export with an independent patch chart."""
from pathlib import Path
import hashlib
import json
import subprocess

import numpy as np
from PIL import Image

import instagram_export as export


def main():
    root = Path(__file__).resolve().parent
    out = root / 'output-schwarzschild' / 'export-audit'
    out.mkdir(parents=True, exist_ok=True)
    width, height = 1080, 1920
    chart = np.zeros((height, width, 3), np.uint8)
    colors = [(0,0,0), (1,1,1), (4,4,4), (8,8,8), (16,16,16), (32,32,32),
              (64,64,64), (128,128,128), (192,192,192), (255,255,255),
              (255,128,0), (240,80,12), (180,20,0), (255,215,140),
              (240,240,234), (128,80,20)]
    for index, color in enumerate(colors):
        row, col = divmod(index, 4)
        chart[row*480:(row+1)*480, col*270:(col+1)*270] = color
    Image.fromarray(chart).save(out / 'test-chart.png')
    video = out / 'test-chart.mp4'
    process = export.encoder(video, width, height, 30)
    for _ in range(30):
        process.stdin.write(chart.tobytes())
    process.stdin.close()
    if process.wait() != 0:
        raise RuntimeError('Chart encoding failed')
    stream = json.loads(subprocess.check_output([
        'ffprobe', '-v', 'error', '-show_streams', '-of', 'json', str(video)
    ]))['streams'][0]
    decoded = np.frombuffer(subprocess.check_output([
        'ffmpeg', '-v', 'error', '-i', str(video), '-frames:v', '1',
        '-vf', export.BT709_TO_SRGB_FILTER, '-pix_fmt', 'rgb24', '-f', 'rawvideo', '-'
    ]), np.uint8).reshape(height, width, 3)
    Image.fromarray(decoded).save(out / 'decoded-chart-srgb.png')
    yuv = np.frombuffer(subprocess.check_output([
        'ffmpeg', '-v', 'error', '-i', str(video), '-frames:v', '1',
        '-pix_fmt', 'yuv420p', '-f', 'rawvideo', '-'
    ]), np.uint8)
    size = width * height
    y = yuv[:size].reshape(height, width)
    u = yuv[size:size + size//4].reshape(height//2, width//2)
    v = yuv[size + size//4:].reshape(height//2, width//2)
    patches = []
    gray_errors = []
    for index, color in enumerate(colors):
        row, col = divmod(index, 4)
        interior = (slice(row*480+32, (row+1)*480-32),
                    slice(col*270+32, (col+1)*270-32))
        error = np.abs(decoded[interior].astype(float) - chart[interior])
        patch = {'input_srgb': color,
                 'decoded_srgb_median': np.median(decoded[interior], axis=(0,1)).tolist(),
                 'interior_mae_codes': float(error.mean()),
                 'interior_max_error_codes': float(error.max())}
        if color[0] == color[1] == color[2]:
            # Independent display-referred prediction, not another zscale call.
            srgb = color[0] / 255.0
            linear = srgb / 12.92 if srgb <= .04045 else ((srgb + .055)/1.055)**2.4
            expected_y = 16 + 219 * linear**(1/2.4)
            measured_y = float(np.median(y[interior]))
            patch.update(expected_luma_unrounded=expected_y, measured_luma=measured_y)
            gray_errors.append(abs(measured_y-expected_y))
        patches.append(patch)
    checks = {
        'h264_high_yuv420p': stream['codec_name'] == 'h264' and stream['profile'] == 'High'
                              and stream['pix_fmt'] == 'yuv420p',
        'native_dimensions_and_rate': stream['width'] == width and stream['height'] == height
                                     and stream['avg_frame_rate'] == '30/1'
                                     and stream['nb_frames'] == '30',
        'declared_color_tags': all(stream.get(k) == val for k, val in
                                  export.metadata()['output_color_tags'].items()),
        'black_srgb_exact': bool(np.all(decoded[32:448, 32:238] == 0)),
        'black_yuv_exact': bool(np.all(y[32:448,32:238] == 16)
                               and np.all(u[16:224,16:119] == 128)
                               and np.all(v[16:224,16:119] == 128)),
        'patch_mean_error_below_3_codes': max(p['interior_mae_codes'] for p in patches) < 3,
        'independent_gray_transfer_within_1_code': max(gray_errors) <= 1,
        'actual_transfer_conversion': abs(int(y[720,945]) - (16 + 219*128/255)) > 3,
        'no_b_frames': stream['has_b_frames'] == 0,
    }
    report = {
        'passed': all(checks.values()), 'checks': checks,
        'scope': 'Local synthetic uniform-patch encode/decode test; not an Instagram upload test',
        'pipeline': export.metadata(),
        'pipeline_source_sha256': hashlib.sha256(Path(export.__file__).read_bytes()).hexdigest(),
        'audit_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'test_file': str(video.relative_to(root)),
        'stream': {k: stream.get(k) for k in ['codec_name', 'profile', 'pix_fmt', 'width',
                   'height', 'avg_frame_rate', 'nb_frames', 'color_space', 'color_transfer',
                   'color_primaries', 'color_range', 'chroma_location', 'bit_rate']},
        'patches': patches,
        'maximum_independent_gray_luma_error_codes': max(gray_errors),
        'black_yuv_codes': [int(y[240,135]), int(u[120,67]), int(v[120,67])],
        'midgray_128_encoded_luma': int(y[720,945]),
    }
    (root / 'instagram-export-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': report['passed'], 'checks': checks,
                      'maximum_independent_gray_luma_error_codes': max(gray_errors)}, indent=2))
    if not report['passed']:
        raise SystemExit('Local export audit failed')


if __name__ == '__main__':
    main()
