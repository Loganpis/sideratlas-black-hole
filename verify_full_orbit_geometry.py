"""Regression and symmetry checks for the full-circle Schwarzschild camera.

This checks ray geometry and relative propagation times, not a fluid model or
the undocumented implementation of NASA's visualization. Run directly to
refresh full-orbit-geometry-diagnostics.json beside this file.
"""
from __future__ import annotations

import gc
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from nasa_raytrace import SchwarzschildDelayTracer, tortoise_radius
from raytrace import SchwarzschildTracer


def diagnostics():
    started = time.perf_counter()
    options = dict(width=432,height=360,fov=72,roll_deg=0,rin=6,rout=28,
                   layers=4,radial_samples=8192,angle_step=.002)
    tracer = SchwarzschildDelayTracer(**options)
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise AssertionError(name)

    # The old tracer can use this same orbit table. These comparisons isolate
    # the camera generalization from an unrelated integration refinement.
    for angle in (0.,30.,68.,83.,89.):
        old, new = SchwarzschildTracer.trace(tracer,angle), tracer.trace(angle)
        check(f'baseline_above_disk_bitwise_{angle:g}',
              all(np.array_equal(old[key],new[key]) for key in old))

    views = []
    angles = (83.,0.,-7.,-89.999999,-90.,-90.000001,-97.,-180.,
              -269.999999,-270.,-270.000001,-277.,90.,180.,270.,443.)
    for angle in angles:
        maps = tracer.trace(angle)
        valid = maps['mask']
        label = f'{angle:.9g}'
        check(f'finite_and_positive_transfer_{label}',
              all(np.isfinite(maps[key]).all() for key in ('r','phi','g','delay_tg'))
              and np.all(maps['g'][valid] > 0))
        check(f'disk_bounds_and_invalid_values_{label}',
              np.all((maps['r'][valid] >= 6) & (maps['r'][valid] <= 28))
              and np.all(maps['delay_tg'][~valid] == 0)
              and np.all(maps['hit_layer'][~valid] == -1))
        check(f'left_right_geometric_symmetry_{label}',
              all(np.array_equal(maps[key],np.fliplr(maps[key]))
                  for key in ('r','mask','delay_tg')))
        views.append({'angle_deg':angle,'disk_pixels':int(valid.sum()),
                      'first_hit_layers':{str(k):int((maps['hit_layer']==k).sum())
                                          for k in range(4)},
                      'delay_range_tg':[float(maps['delay_tg'][valid].min()),
                                        float(maps['delay_tg'][valid].max())]})

    side, opposite = tracer.trace(83.), tracer.trace(-97.)
    check('antipodal_radius_delay_mask_bitwise',
          all(np.array_equal(side[key],opposite[key])
              for key in ('r','mask','delay_tg','hit_layer')))
    check('antipodal_doppler_reflection_bitwise',
          np.array_equal(side['g'],np.fliplr(opposite['g'])))
    check('antipodal_emission_azimuth_reverses',
          np.allclose(np.sin(side['phi']),-np.sin(opposite['phi']),atol=2e-7)
          and np.allclose(np.cos(side['phi']),np.cos(opposite['phi']),atol=2e-7))
    for angle in (-277.,443.):
        repeat = tracer.trace(angle)
        check(f'full_revolution_identity_{angle:g}',
              all(np.array_equal(side[key],repeat[key]) for key in side))

    # At the poles, orbital motion is transverse: no line-of-sight Doppler
    # term. The geometrical image and propagation time are circular.
    for angle in (0.,-180.):
        maps = tracer.trace(angle)
        valid = maps['mask']
        check(f'polar_frequency_shift_{angle:g}',
              np.allclose(maps['g'][valid],np.sqrt(1-3/maps['r'][valid]),atol=2e-7))
        left = (options['width']-options['height'])//2
        square = (slice(None),slice(left,left+options['height']))
        check(f'polar_quarter_turn_geometry_{angle:g}',
              all(np.array_equal(maps[key][square],np.rot90(maps[key][square]))
                  for key in ('r','mask','delay_tg')))

    # For an even detector no ray is exactly coplanar. The exact-edge maps
    # approach the same physical image from above and below, even though the
    # root-count label changes when the observer-at-infinity root is skipped.
    edge_limits = []
    for edge in (-90.,-270.):
        exact = tracer.trace(edge)
        for offset in (-1e-6,1e-6):
            near = tracer.trace(edge+offset)
            common = exact['mask'] & near['mask']
            radius_error = float(np.max(np.abs(exact['r'][common]-near['r'][common])))
            delay_error = float(np.max(np.abs(exact['delay_tg'][common]-near['delay_tg'][common])))
            disagreement = float(np.mean(exact['mask'] != near['mask']))
            check(f'edge_limit_{edge:g}_{offset:g}',
                  disagreement < 1e-4 and radius_error < .002 and delay_error < .002)
            edge_limits.append({'edge_angle_deg':edge,'offset_deg':offset,
                                'mask_disagreement_fraction':disagreement,
                                'max_radius_difference_rg':radius_error,
                                'max_delay_difference_tg':delay_error})

    del tracer
    gc.collect()
    odd = SchwarzschildDelayTracer(**{**options,'width':217,'height':181})
    rim_impact = options['rout']/math.sqrt(1-2/options['rout'])
    for angle in (-90.,90.):
        maps = odd.trace(angle)
        row = 90
        expected = odd.impact[row] <= rim_impact
        check(f'coplanar_first_outer_rim_{angle:g}',
              np.array_equal(maps['mask'][row],expected)
              and np.all(maps['r'][row,expected] == options['rout'])
              and np.all(maps['hit_layer'][row,expected] == 0))
        check(f'central_radial_delay_{angle:g}',
              abs(float(maps['delay_tg'][90,108])+float(tortoise_radius(28))) < 2e-6)

    root = Path(__file__).resolve().parent
    source_sha256 = {name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                     for name in ('nasa_raytrace.py','raytrace.py',Path(__file__).name)}
    return {'passed':True,'scene':options,'checks':checks,'views':views,
            'source_sha256':source_sha256,
            'disk_plane_limits':edge_limits,
            'edge_treatment':'Infinitesimally thin opaque disk. At the disk plane, '
                             'noncoplanar rays first intersect at psi=pi; lensed rings persist. '
                             'Exactly coplanar rays first enter the outer rim. '
                             'This is a mathematical thin-surface limit, not a finite-height flow.',
            'scope':'Geometry and travel-time regression/symmetry checks; not a fluid-model '
                    'validation or a reproduction of unpublished NASA parameters.',
            'elapsed_seconds':time.perf_counter()-started}


if __name__ == '__main__':
    result = diagnostics()
    output = Path(__file__).with_name('full-orbit-geometry-diagnostics.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':result['passed'],'checks':len(result['checks']),
                      'elapsed_seconds':result['elapsed_seconds'],'output':str(output)},indent=2))
