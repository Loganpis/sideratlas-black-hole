"""Single-generation SDR H.264 export of full-range sRGB animation frames.

The color transform changes sample values as well as labels. Feed RGB24 frames
to the returned process's stdin, then close stdin and check its return code.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import shutil
import subprocess


CRF = 14
MAXRATE = "20M"
BUFSIZE = "40M"

# Expand RGB to 16-bit planar storage before the actual zimg color transform.
# No resize, grain, sharpening, or stochastic dither is introduced here.
SRGB_TO_BT709_FILTER = (
    "format=gbrp16le,"
    "zscale=primariesin=709:transferin=iec61966-2-1:matrixin=gbr:rangein=full:"
    "primaries=709:transfer=709:matrix=709:range=limited:"
    "chromal=left:filter=spline36:agamma=0:dither=none,"
    "format=yuv420p"
)

# Use this for comparisons with the original sRGB PNG/RGB frames. A matrix-only
# RGB decode leaves the Rec.709 transfer encoded and is not a color-managed QA.
BT709_TO_SRGB_FILTER = (
    "zscale=primariesin=709:transferin=709:matrixin=709:rangein=limited:chromalin=left:"
    "primaries=709:transfer=iec61966-2-1:matrix=gbr:range=full:"
    "filter=spline36:agamma=0:dither=none,"
    "format=gbrp16le,format=rgb24"
)


@lru_cache(maxsize=1)
def _ffmpeg() -> str:
    executable = shutil.which("ffmpeg")
    if executable is None:
        raise RuntimeError("FFmpeg with the zscale filter is required for SDR export")
    probe = subprocess.run(
        [executable, "-hide_banner", "-h", "filter=zscale"],
        check=True, capture_output=True, text=True,
    )
    if "Filter zscale" not in probe.stdout + probe.stderr:
        raise RuntimeError("This FFmpeg build needs libzimg/zscale for the color conversion")
    return executable


def encoder(path: str | Path, width: int, height: int, fps: int):
    """Start an encoder accepting full-range, sRGB-transfer RGB24 frames."""
    if not all(isinstance(value, int) and value > 0 for value in (width, height, fps)):
        raise ValueError("Width, height and frame rate must be positive integers")
    if width % 2 or height % 2:
        raise ValueError("4:2:0 output requires even frame dimensions")
    command = [
        _ffmpeg(), "-hide_banner", "-loglevel", "error", "-y",
        "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{width}x{height}",
        "-framerate", str(fps), "-i", "-", "-an",
        "-vf", SRGB_TO_BT709_FILTER,
        "-c:v", "libx264", "-preset", "slow", "-crf", str(CRF),
        "-maxrate", MAXRATE, "-bufsize", BUFSIZE,
        "-profile:v", "high", "-pix_fmt", "yuv420p",
        "-r", str(fps), "-fps_mode", "cfr",
        "-g", str(2 * fps), "-keyint_min", str(2 * fps),
        "-bf", "0", "-sc_threshold", "0", "-flags", "+cgop",
        "-color_primaries", "bt709", "-color_trc", "bt709",
        "-colorspace", "bt709", "-color_range", "tv",
        "-chroma_sample_location", "left",
        "-bsf:v", "h264_metadata=video_full_range_flag=0:colour_primaries=1:"
        "transfer_characteristics=1:matrix_coefficients=1",
        "-movflags", "+faststart", "-use_editlist", "0", str(path),
    ]
    return subprocess.Popen(command, stdin=subprocess.PIPE)


def metadata() -> dict:
    """Describe the export contract; reported settings are not measured results."""
    return {
        "pipeline_version": 1,
        "container": "MP4", "codec": "H.264", "encoder": "libx264",
        "profile": "High", "preset": "slow", "crf": CRF,
        "maximum_rate_bps": 20_000_000, "vbv_buffer_bits": 40_000_000,
        "rate_control": "quality-targeted VBR with a VBV cap; CRF is not a platform requirement",
        "pixel_format": "yuv420p", "progressive": True,
        "frame_rate_mode": "constant", "gop_seconds": 2,
        "closed_gop": True, "b_frames": 0, "scene_cut_keyframes": False,
        "faststart": True, "edit_lists": False, "audio": "none",
        "input": {"pixel_format": "rgb24", "primaries": "bt709",
                  "transfer": "iec61966-2-1", "range": "full"},
        "output_color_tags": {"color_space": "bt709", "color_transfer": "bt709",
                              "color_primaries": "bt709", "color_range": "tv",
                              "chroma_location": "left"},
        "color_conversion": "sRGB RGB to Rec.709 limited-range YCbCr via zscale; pixels transformed before tagging",
        "conversion_reference": "display-referred: sRGB EOTF to linear light, then inverse BT.1886 (gamma 2.4) for Rec.709 display encoding",
        "filter": SRGB_TO_BT709_FILTER,
        "resize": "none", "added_sharpening": False, "added_grain": False,
        "instagram_processing": "Instagram generates its own delivery encodings; playback quality is not guaranteed by local export",
    }
