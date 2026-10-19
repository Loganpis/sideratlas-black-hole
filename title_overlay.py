"""Presentation-only title layout for the final Sideratlas delivery.

The simulation's frozen physical source and native radiance are independent
of this module. Assembly applies these titles before the single video encode.
"""
from pathlib import Path
import hashlib

from PIL import Image, ImageDraw

RELEASE_VERSION = '1.2.1'


def centered_overlay(renderer):
    overlay = Image.new('RGBA', (renderer.width, renderer.height))
    draw = ImageDraw.Draw(overlay)
    scale = renderer.scale
    stroke = max(1, round(2*scale))
    options = {'anchor':'lt', 'stroke_width':stroke,
               'stroke_fill':(0, 0, 0, 220)}
    draw.text((72*scale, 270*scale), '@sideratlas',
              font=renderer.font(34), fill=(245,245,245), **options)
    draw.text((965*scale, 277*scale), 'NO.001',
              font=renderer.font(25), fill=(208,208,208),
              **{**options, 'anchor':'rt'})
    title = 'THE BLACK HOLE'
    font = renderer.font(82, 'bold')
    left, _, right, _ = draw.textbbox((0, 0), title, font=font,
                                      anchor='lt', stroke_width=stroke)
    x = (renderer.width-left-right)/2
    draw.text((x, 333*scale), title, font=font,
              fill=(245,245,245), **options)
    return overlay


def apply_title_layout(renderer):
    overlay = centered_overlay(renderer)
    renderer.overlays = {mode:overlay.copy() for mode in ('color','mono')}
    return renderer


def provenance():
    path = Path(__file__).resolve()
    return {'version':RELEASE_VERSION, 'module':path.name,
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'title':'THE BLACK HOLE', 'alignment':'horizontal center',
            'title_y_at_1080':333, 'account_and_number_positions_unchanged':True}
