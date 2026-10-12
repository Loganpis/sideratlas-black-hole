#!/usr/bin/env python3
"""Check the original Schwarzschild Reel and its disclosed numerical model.

--preflight runs current-code renderer checks without requiring completed MP4s.
The saved ray/source diagnostics are recorded as prior evidence, not silently
rerun or relabeled as validation of the delivered scene. These checks establish
numerical consistency and export integrity, not GRMHD or a NASA reconstruction.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
from pathlib import Path
import shutil
import time

import numpy as np
from PIL import Image

from instagram_export import metadata as export_metadata
from render_nasa import (
    DURATION, FPS, ROOT, FOV, OUTER_RADIUS, INCLINATION, LOOP_TG,
    NasaInspiredRenderer, VIDEO_NAMES, OUTPUT_DIR, display_intensity,
    illustrative_palette, source_hash, view_angle_at_phase,
)
from verify import check, verify_video


def finish(record):
    record['passed'] = all(item['passed'] for item in record.get('checks', []))
    return record


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1024*1024), b''):
            digest.update(data)
    return digest.hexdigest()


def saved_diagnostics():
    result = {'scope': (
        'Prior saved diagnostics, not rerun by this verifier. Their files have no '
        'original source hashes, so these hashes identify the reports, not the code '
        'that produced them. Saved production-view tests use 82/83.5/85 degrees '
        'and a different hero height; current full-orbit camera checks are separate.'),
        'reports': {}}
    for name in ('nasa-raytrace-diagnostics.json', 'nasa-emission-diagnostics.json'):
        path = ROOT/name
        data = json.loads(path.read_text())
        result['reports'][name] = {'sha256': sha256(path), 'data': data}
    rays = result['reports']['nasa-raytrace-diagnostics.json']['data']
    source = result['reports']['nasa-emission-diagnostics.json']['data']
    check(result, 'saved_ray_diagnostics_report_pass', rays.get('passed') is True)
    check(result, 'saved_source_diagnostics_report_consistent_invariants',
          source['azimuthal_mean_max_relative_error'] < 1e-12 and
          source['periodicity_max_absolute_error'] < 1e-12 and
          source['independent_scalar_equation_max_absolute_error'] < 1e-12 and
          source['prepared_vs_direct_max_absolute_error'] < 1e-12 and
          source['positive_emission_min'] >= 0 and
          source['smooth_mode_max_absolute_error'] == 0)
    return finish(result)


def display_checks():
    values = np.linspace(0, 1, 65537)
    palette = illustrative_palette(values).astype(float)/255
    linear = np.where(palette <= .04045, palette/12.92,
                      ((palette+.055)/1.055)**2.4)
    luminance = linear @ np.array([.2126, .7152, .0722])
    transformed = display_intensity(np.geomspace(1e-12, 1e12, 10000))
    result = {'palette_scope': 'Display colors encode brightness; no spectrum is inferred.'}
    check(result, 'black_is_exactly_black', np.array_equal(palette[0], np.zeros(3)))
    check(result, 'palette_linear_srgb_luminance_is_nondecreasing',
          np.all(np.diff(luminance) >= -1e-15),
          minimum_luminance_step=float(np.diff(luminance).min()))
    check(result, 'fixed_tone_curve_is_finite_bounded_and_monotonic',
          np.isfinite(transformed).all() and np.all(transformed >= 0) and
          np.all(transformed <= 1) and np.all(np.diff(transformed) >= 0))
    return finish(result)


def actual_renderer_checks():
    renderer = NasaInspiredRenderer(width=270, supersample=2)
    result = {'scope': (
        'Fresh checks using the actual render path at 270x480, 2x spatial '
        'supersampling, including a full 360-degree observer revolution. Equal endpoint '
        'frames alone do not imply an evolving source is physically periodic.'),
        'geometry': renderer.geometry, 'emission': renderer.source.metadata()}
    phases = (0, .125, .25, .375, .5, .625, .75, .875, 1)
    angles = [view_angle_at_phase(phase) for phase in phases]
    north_phase, south_phase = INCLINATION/360, (INCLINATION+180)/360
    check(result, 'camera_completes_one_revolution_with_both_poles',
          np.allclose(angles, [83,38,-7,-52,-97,-142,-187,-232,83], rtol=0, atol=1e-12) and
          math.isclose(view_angle_at_phase(north_phase), 0, abs_tol=1e-12) and
          math.isclose(view_angle_at_phase(south_phase), -180, abs_tol=1e-12),
          sample_phases=list(phases), sample_view_angles_degrees=angles,
          north_pole_phase=north_phase, south_pole_phase=south_phase)
    def camera_basis(phase):
        angle = math.radians(view_angle_at_phase(phase))
        si, ci = math.sin(angle), math.cos(angle)
        return np.array([[0,-si,ci], [1,0,0], [0,ci,si]])
    epsilon = 1e-5
    velocity_before = (camera_basis(0)-camera_basis(-epsilon))/epsilon
    velocity_after = (camera_basis(epsilon)-camera_basis(0))/epsilon
    check(result, 'camera_orientation_and_angular_velocity_join_at_loop_boundary',
          np.allclose(camera_basis(0), camera_basis(1), rtol=0, atol=1e-12) and
          np.allclose(velocity_before, velocity_after, rtol=0, atol=4e-4),
          maximum_basis_velocity_difference=float(np.max(np.abs(velocity_before-velocity_after))),
          scope='Continuous detector basis, not continuity of a wrapped angle label.')

    # Save the side view before changing the observer. Source periodicity is
    # checked in physical time, independently of the renderer's phase modulo.
    maps, mask = renderer.maps, renderer.maps['mask']
    start = renderer.prepared.evaluate(0.0)
    finish_source = renderer.prepared.evaluate(LOOP_TG)
    error = float(np.max(np.abs(start-finish_source)))
    check(result, 'unwrapped_retarded_source_is_periodic', error < 1e-12,
          maximum_absolute_relative_intensity_error=error)
    steady0 = renderer.prepared.evaluate(0, structured=False)
    steady1 = renderer.prepared.evaluate(.375*LOOP_TG, structured=False)
    check(result, 'smooth_source_is_time_independent_at_fixed_view',
          np.array_equal(steady0, steady1))
    check(result, 'missed_rays_are_physically_black_before_resampling', np.all(start[~mask] == 0))
    no_delay = renderer.source.prepare(maps['r'], maps['phi'], 0).evaluate(0)
    delay_difference = float(np.sqrt(np.mean((start[mask]-no_delay[mask])**2)))
    check(result, 'retarded_emission_differs_from_instantaneous_emission', delay_difference > 0,
          rms_source_difference_on_emitting_rays=delay_difference,
          relative_delay_range_tg=[float(maps['delay_tg'][mask].min()), float(maps['delay_tg'][mask].max())])
    layers = {str(k): int(np.sum(maps['hit_layer'] == k)) for k in range(4)}
    check(result, 'side_view_primary_and_secondary_disk_images_present',
          layers['0'] > 0 and layers['1'] > 0, first_visible_hit_counts=layers,
          scope='First opaque hit order; this does not assert that every photon subring is resolved.')

    first = renderer.frames(0)
    last = renderer.frames(1)
    for mode in VIDEO_NAMES:
        check(result, mode+'_uncompressed_phase_zero_equals_one',
              np.array_equal(np.asarray(first[mode]), np.asarray(last[mode])))
    radiances = []
    for phase, angle in zip(phases, angles):
        radiances.append(renderer.radiance(phase))
        view = renderer.maps
        hit = view['mask']
        detail = {'phase': phase, 'view_angle_degrees': angle}
        check(result, f'phase_{phase:g}_maps_are_finite',
              all(np.isfinite(view[k]).all() for k in ('r', 'phi', 'g', 'delay_tg')), **detail)
        check(result, f'phase_{phase:g}_emitting_hits_inside_disk_with_positive_shift',
              hit.any() and np.all(view['r'][hit] >= 6) and
              np.all(view['r'][hit] <= OUTER_RADIUS) and np.all(view['g'][hit] > 0), **detail)
        check(result, f'phase_{phase:g}_missed_rays_have_zero_delay',
              np.all(view['delay_tg'][~hit] == 0), **detail)
        check(result, f'phase_{phase:g}_geometric_mask_and_delay_reflect_left_right',
              np.array_equal(hit, np.fliplr(hit)) and
              np.allclose(view['delay_tg'], np.fliplr(view['delay_tg']), rtol=0, atol=1e-5), **detail)
    smooth_poles = []
    for pole_name, pole_phase in (('north', north_phase), ('south', south_phase)):
        smooth_pole = renderer.radiance(pole_phase, structured=False)
        smooth_poles.append(smooth_pole)
        view, hit = renderer.maps, renderer.maps['mask']
        # Face-on gas has no line-of-sight orbital velocity. Gravitational
        # and transverse-Doppler redshift remain on either disk face.
        expected_shift = np.sqrt(1-3/view['r'][hit].astype(np.float64))
        shift_error = float(np.max(np.abs(view['g'][hit]-expected_shift)))
        check(result, pole_name+'_pole_has_correct_gravitational_and_transverse_doppler_shift',
              shift_error < 1e-6, maximum_absolute_frequency_shift_error=shift_error)
        h, w = hit.shape
        square = slice((w-h)//2, (w+h)//2)
        check(result, pole_name+'_pole_disk_geometry_has_quarter_turn_symmetry',
              np.array_equal(hit[:, square], np.rot90(hit[:, square])) and
              all(np.allclose(view[k][:, square], np.rot90(view[k][:, square]),
                              rtol=0, atol=1e-5) for k in ('r','g','delay_tg')))
        pole_left = float(smooth_pole[:, :135].sum())
        pole_right = float(smooth_pole[:, 135:].sum())
        check(result, pole_name+'_pole_smooth_reference_has_no_azimuthal_doppler_asymmetry',
              math.isclose(pole_left, pole_right, rel_tol=1e-6),
              left_integrated_intensity=pole_left, right_integrated_intensity=pole_right)
    check(result, 'identical_disk_faces_produce_equal_axisymmetric_pole_images',
          np.allclose(smooth_poles[0], smooth_poles[1], rtol=1e-6, atol=1e-6))
    changes = [float(np.sqrt(np.mean((frame-radiances[0])**2))) for frame in radiances[1:]]
    check(result, 'sampled_radiance_is_finite_nonnegative_and_dynamic',
          all(np.isfinite(frame).all() and np.all(frame >= 0) for frame in radiances)
          and max(changes) > 0,
          sample_phases=list(phases), rms_changes_from_phase_zero=changes)
    smooth_side = renderer.radiance(0, structured=False)
    check(result, 'smooth_reference_changes_between_side_and_pole',
          not np.array_equal(smooth_side, smooth_poles[0]))
    left, right = float(smooth_side[:, :135].sum()), float(smooth_side[:, 135:].sum())
    check(result, 'side_smooth_reference_has_approaching_side_doppler_asymmetry', left > right,
          left_integrated_intensity=left, right_integrated_intensity=right,
          left_right_ratio=left/right)
    smooth_opposite = renderer.radiance(.5, structured=False)
    opposite = renderer.maps
    check(result, 'opposite_observer_keeps_radius_delay_and_reverses_doppler_side',
          np.array_equal(maps['mask'], opposite['mask']) and
          all(np.allclose(maps[k], opposite[k], rtol=0, atol=1e-5)
              for k in ('r', 'delay_tg')) and
          np.allclose(maps['g'], np.fliplr(opposite['g']), rtol=0, atol=1e-6) and
          np.allclose(smooth_side, np.fliplr(smooth_opposite), rtol=1e-5, atol=1e-6),
          opposite_left_integrated_intensity=float(smooth_opposite[:, :135].sum()),
          opposite_right_integrated_intensity=float(smooth_opposite[:, 135:].sum()),
          scope='Axisymmetric source and identical disk faces with a continuous camera basis.')
    check(result, 'structured_source_evolves_at_a_fixed_view_during_full_orbit',
          not np.array_equal(start, renderer.source.prepare(maps['r'], maps['phi'], maps['delay_tg']).evaluate(.5*LOOP_TG)))

    # Exact plane crossings are not removed or faded. For an infinitesimally
    # thin opaque disk, noncoplanar rays first meet its plane again at psi=pi.
    # A finite spatial detector may not resolve all the resulting thin images.
    for plane_angle in (-90., -270.):
        crossing = []
        for offset in (-1e-6, 0., 1e-6):
            view = renderer.tracer.trace(plane_angle+offset)
            hit = view['mask']
            valid = (all(np.isfinite(view[k]).all() for k in ('r','phi','g','delay_tg')) and
                     np.all(view['r'][hit] >= 6) and np.all(view['r'][hit] <= OUTER_RADIUS) and
                     np.all(view['g'][hit] > 0) and np.all(view['delay_tg'][~hit] == 0))
            crossing.append({'offset_degrees': offset, 'emitting_pixels': int(hit.sum()), 'valid': bool(valid)})
        check(result, f'plane_{plane_angle:g}_exact_and_two_sided_limits_are_finite',
              all(item['valid'] for item in crossing) and crossing[1]['emitting_pixels'] > 0,
              samples=crossing,
              scope='Tests finite physical maps, not a guarantee that narrow image orders are spatially resolved.')
    full_turn = renderer.tracer.trace(INCLINATION-360)
    check(result, 'geodesic_maps_return_after_full_camera_revolution',
          np.array_equal(maps['mask'], full_turn['mask']) and
          all(np.allclose(maps[k], full_turn[k], rtol=0, atol=1e-5)
              for k in ('r','phi','g','delay_tg')))

    # Detector refinement is reported, not turned into a made-up accuracy bound.
    # Both images cover identical physical fields. Only detector sampling changes.
    reference_low = radiances[2].copy()
    del renderer, maps, start, finish_source, radiances, no_delay, view
    gc.collect()
    high = NasaInspiredRenderer(width=540, supersample=2)
    reference_high = high.radiance(.25)
    reduced = np.asarray(Image.fromarray(reference_high).resize((270,225), Image.Resampling.BOX))
    error = reference_low-reduced
    peak = float(reduced.max())
    result['detector_refinement'] = {
        'phase': .25, 'delivered_widths': [270,540], 'rays_per_pixel_per_axis': 2,
        'comparison': '540-pixel radiance area-averaged to 270 pixels before either display transform',
        'rms_difference_relative_to_high_resolution_peak': float(np.sqrt(np.mean(error**2))/peak),
        'maximum_absolute_difference_relative_to_high_resolution_peak': float(np.max(np.abs(error))/peak),
        'integrated_radiance_relative_difference': float(abs(reference_low.sum()/reduced.sum()-1)),
        'scope': 'Sampling sensitivity at one phase, not an error bound for the full-resolution movie; thin higher-order images may remain unresolved.',
    }
    check(result, 'detector_refinement_produces_finite_nonnegative_images',
          np.isfinite(reduced).all() and np.all(reduced >= 0) and peak > 0)
    del high
    gc.collect()
    return finish(result)


def metadata_checks(output):
    path = output/'render-metadata.json'
    data = json.loads(path.read_text())
    result = {'file': str(path), 'sha256': sha256(path)}
    check(result, 'recorded_source_matches_current_renderer_dependencies', data['source_sha256'] == source_hash())
    check(result, 'delivery_dimensions_timing_and_supersampling',
          (data['width'],data['height'],data['fps'],data['frames'],data['duration_seconds'],data['supersample']) ==
          (1080,1920,30,240,8,2))
    geom = data['geometry']
    check(result, 'selected_schwarzschild_geometry',
          data['inclination_degrees'] == [0, 180] and INCLINATION == 83 and
          geom['width'] == 2160 and geom['height'] == 1800 and geom['rin_rg'] == 6 and
          geom['rout_rg'] == OUTER_RADIUS == 28 and geom['fov_horizontal_rg'] == FOV == 72 and
          geom['roll_deg'] == 0 and geom['layers'] == 4 and
          math.isclose(geom['critical_impact_rg'], math.sqrt(27), abs_tol=1e-12))
    orbit = data['camera_orbit']
    check(result, 'recorded_camera_orbit_matches_full_revolution',
          data['view_angle_degrees'] == [-277,83] and
          [orbit['start_deg'], orbit['angular_span_deg'], orbit['end_deg']] == [83,-360,83] and
          math.isclose(orbit['north_pole_phase'], 83/360, abs_tol=1e-12) and
          math.isclose(orbit['south_pole_phase'], 263/360, abs_tol=1e-12) and
          '360' in orbit['law'] and data['ray_maps']['inclination_degrees'] == 83)
    check(result, 'both_disk_faces_and_static_observer_assumptions_disclosed',
          'both faces' in data['disk_faces'] and
          'static observers' in orbit['observer'] and 'forward' in orbit['gas_time_direction'])
    check(result, 'recorded_source_matches_current_parameters',
          data['emission'] == __import__('nasa_emission').ShearingDiskEmission().metadata())
    check(result, 'all_four_image_layers_hit_in_production_map',
          all(data['ray_maps']['image_layers'].get(str(k),0) > 0 for k in range(4)))
    check(result, 'production_frequency_shifts_and_delays_finite',
          np.isfinite(data['ray_maps']['frequency_shift_range']).all() and
          np.isfinite(data['ray_maps']['delay_tg_range']).all() and
          0 < data['ray_maps']['frequency_shift_range'][0] < 1 < data['ray_maps']['frequency_shift_range'][1])
    check(result, 'both_output_names_match', data['outputs'] == VIDEO_NAMES)
    check(result, 'recorded_encoding_matches_current_export_pipeline',
          data['encoding'] == export_metadata())
    result['layout'] = data['layout']
    check(result, 'brand_and_header_have_requested_top_clearance',
          data['layout']['brand'] == '@sideratlas' and
          data['layout']['header_top'] >= 250)
    check(result, 'recorded_hero_placement_matches_video_review_region',
          data['layout']['hero_y'] == 425 and data['layout']['hero_height'] == 720 and
          data['layout']['hero_width'] == 864 and data['layout']['hero_x'] == 108)
    check(result, 'hero_display_preserves_radiance_aspect_ratio',
          data['layout']['radiance_height'] == 900 and
          math.isclose(data['layout']['hero_width']/data['layout']['hero_height'],
                       data['width']/data['layout']['radiance_height'], abs_tol=1e-12))
    check(result, 'render_has_finished', data.get('elapsed_seconds',0) > 0)
    return finish(result)


def output_assets(output):
    result = {'scope': 'Hashes identify the reviewed output files; no preexisting checksum manifest is assumed.', 'files': {}}
    names = list(VIDEO_NAMES.values()) + ['render-metadata.json', 'sideratlas-black-hole-color-cover.png',
        'sideratlas-black-hole-mono-cover.png', 'steady-reference.png', 'contact-sheet.jpg']
    for name in names:
        path = output/name
        check(result, name+'_exists', path.is_file() and path.stat().st_size > 0)
        if path.is_file():
            result['files'][name] = {'sha256': sha256(path), 'size_bytes': path.stat().st_size}
            if path.suffix in ('.png','.jpg'):
                with Image.open(path) as picture:
                    picture.load()
                    check(result, name+'_image_decodes', picture.size == ((1080,1280) if name == 'contact-sheet.jpg' else (1080,1920)))
    check(result, 'no_incomplete_mp4_remains', not list(output.glob('*.partial.mp4')))
    return finish(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT_DIR)
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {'preflight_only': args.preflight, 'videos': {}, 'errors': [],
              'scope': 'Numerical consistency and export integrity of a disclosed idealized model; not physical validation of an accretion flow or an exact NASA reconstruction.',
              'verifier_sha256': sha256(Path(__file__))}
    functions = [('saved_diagnostics',saved_diagnostics),('display',display_checks),
                 ('actual_renderer',actual_renderer_checks)]
    if not args.preflight:
        functions += [('render_metadata', lambda: metadata_checks(args.output)),
                      ('output_assets', lambda: output_assets(args.output))]
    for name, function in functions:
        print(f'Checking {name.replace("_"," ")}.', flush=True)
        try:
            report[name] = function()
        except Exception as error:
            report[name] = {'passed': False, 'error': str(error)}
    if not args.preflight:
        ffmpeg, ffprobe = shutil.which('ffmpeg'), shutil.which('ffprobe')
        if not ffmpeg or not ffprobe:
            report['errors'].append('FFmpeg and ffprobe must be available on PATH.')
        else:
            for mode, filename in VIDEO_NAMES.items():
                print(f'Decoding {filename}.', flush=True)
                try:
                    report['videos'][mode] = verify_video(args.output/filename,ffmpeg,ffprobe,
                                                         hero_region=(425,1145),
                                                         expected_transfer='bt709', expected_range='tv')
                except Exception as error:
                    report['videos'][mode] = {'passed': False, 'error': str(error)}
    report['passed'] = not report['errors'] and all(report[name]['passed'] for name,_ in functions)
    if not args.preflight:
        report['passed'] = report['passed'] and len(report['videos']) == len(VIDEO_NAMES) and all(
            item['passed'] for item in report['videos'].values())
    report['elapsed_seconds'] = time.perf_counter()-started
    destination = args.output/('verification-preflight.json' if args.preflight else 'verification.json')
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(f'{"PASS" if report["passed"] else "FAIL"}: {destination}',flush=True)
    if not report['passed']:
        for name, section in list(report.items())+list(report['videos'].items()):
            if isinstance(section,dict):
                if 'error' in section:
                    print(f'{name}: {section["error"]}')
                for item in section.get('checks',[]):
                    if not item['passed']:
                        print(f'{name}: failed {item["name"]}')
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
