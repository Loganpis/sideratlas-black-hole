#!/usr/bin/env python3
"""Verify the wavelength-resolved Reel, its physical model and optical pipeline.

--preflight runs fresh spectral-optics diagnostics, selected-model ray checks,
uncompressed periodicity, display chromaticity and 10 versus 5 nm convergence.
The default also audits all 121 rendered views and fully decodes both movies.
Saved earlier disk/ray diagnostics are identified as baseline provenance, not
misrepresented as new checks of every rendered ray. No network is used.
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

import render_optical as renderer_module
from render_optical import (OpticalRenderer, ROOT, FPS, DURATION, VIDEO_NAMES,
                            inclination, EXPOSURE, FOV, APERTURE, DISTANCE, SUPPORT, RAY_STEP)
from emission import XYZ_TO_LINEAR_SRGB
from kerr_emission import KerrDiskParameters, KerrEmission
from kerr_raytrace import KerrTracer
from optics import color_correction, limb_factor
from spectral_optics import SpectralCamera, FIRST_NULL_FACTOR
from verify import check, verify_video


def _finished(result):
    result['passed'] = all(item['passed'] for item in result.get('checks', []))
    return result


def verify_selected_model():
    """Exercise the chosen 10^10 solar mass, a*=0.9 model through 88 degrees."""
    parameters = KerrDiskParameters(mass_solar=1e10, eddington_ratio=.05, spin_a_star=.9)
    disk = KerrEmission(parameters)
    result = {'parameters': parameters.metadata(), 'inclinations_deg': [85, 86.5, 88]}
    radii = np.geomspace(parameters.isco_rg+1e-4, 500, 10001)
    reference = KerrEmission(KerrDiskParameters())
    measured = disk.flux(radii)/reference.flux(radii)
    expected = (.05/1e10)/(.1/1e8)
    scaling_error = float(np.max(np.abs(measured/expected-1)))
    check(result, 'mass_and_accretion_flux_scaling', scaling_error < 1e-12,
          maximum_relative_error=scaling_error)
    dense = np.linspace(parameters.isco_rg, 20, 100001)
    maximum_temperature = float(disk.temperature(dense).max())
    result['peak_effective_temperature_k'] = maximum_temperature
    result['peak_color_hardening_factor'] = float(color_correction(maximum_temperature))
    tracer = KerrTracer(320, 480, fov=FOV, roll_deg=14, spin=.9, rout=500, angle_step=RAY_STEP)
    maps = {}
    results = {}
    for angle in result['inclinations_deg']:
        mapping = tracer.trace(angle)
        use = mapping['mask']
        temp = disk.temperature(mapping['r'])
        observed = mapping['g']*color_correction(temp)*temp
        check(result, f'{angle:g}_degree_rays_resolved', not np.any(mapping['status'] >= 3))
        check(result, f'{angle:g}_degree_emission_angles_physical',
              np.all((mapping['mu'][use] > 0) & (mapping['mu'][use] <= 1)))
        check(result, f'{angle:g}_degree_isco_boundary',
              np.all(mapping['r'][use] >= parameters.isco_rg-1e-6))
        check(result, f'{angle:g}_degree_potential_constraints',
              mapping['constraint_residual'].max() < 1e-5,
              maximum_relative_residual=float(mapping['constraint_residual'].max()))
        check(result, f'{angle:g}_degree_temperature_and_shift_finite',
              np.all(np.isfinite(observed)) and np.all(observed >= 0) and
              np.all(mapping['g'][use] > 0) and temp.max() <= maximum_temperature*1.0001)
        results[str(angle)] = {
            'frequency_shift_range': [float(mapping['g'][use].min()), float(mapping['g'][use].max())],
            'maximum_effective_temperature_k': float(temp.max()),
            'maximum_observed_color_temperature_k': float(observed.max()),
            'maximum_emission_angle_cosine': float(mapping['mu'][use].max()),
            'highest_visible_disk_crossing': int(mapping['hit_layer'][use].max()),
        }
        if angle == 88:
            maps['coarse'] = mapping
    fine = KerrTracer(320, 480, fov=FOV, roll_deg=14, spin=.9, rout=500, angle_step=.006).trace(88)
    coarse = maps['coarse']
    shared = coarse['mask'] & fine['mask'] & (coarse['hit_layer'] == fine['hit_layer'])
    relative_r = np.abs(coarse['r'][shared]-fine['r'][shared])/fine['r'][shared]
    relative_g = np.abs(coarse['g'][shared]-fine['g'][shared])/fine['g'][shared]
    disagree = float(np.mean(coarse['mask'] != fine['mask']))
    check(result, '88_degree_integration_step_refinement',
          disagree < .002 and np.percentile(relative_r, 99) < .001 and np.percentile(relative_g, 99) < .001,
          mask_disagreement=disagree,
          radius_p99_relative_error=float(np.percentile(relative_r, 99)),
          frequency_shift_p99_relative_error=float(np.percentile(relative_g, 99)),
          scope='Map radii and frequency shifts are stored as float32; zero differences do not imply exact trajectories.')
    result['sampled_views'] = results
    return _finished(result)


def verify_display():
    """Inspect the actual display curve before 8-bit encoding, not a duplicate."""
    result = {'declared_exposure': EXPOSURE}
    rng = np.random.default_rng(12)
    radiance = 10.**rng.uniform(-1, 14, (4000, 3))
    original_encode = renderer_module.encode_srgb
    try:
        renderer_module.encode_srgb = lambda x: x
        linear = renderer_module.display_rgb(radiance)
        more_exposed = renderer_module.display_rgb(radiance, EXPOSURE*2)
        zero = renderer_module.display_rgb(np.zeros((4, 3)))
    finally:
        renderer_module.encode_srgb = original_encode
    before = radiance/radiance.sum(axis=-1, keepdims=True)
    after = linear/linear.sum(axis=-1, keepdims=True)
    ratio_error = float(np.max(np.abs(before-after)))
    check(result, 'actual_curve_preserves_in_gamut_chromaticity', ratio_error < 1e-12,
          maximum_RGB_ratio_error=ratio_error)
    check(result, 'actual_curve_stays_in_display_range',
          np.all(np.isfinite(linear)) and np.all((linear >= 0) & (linear <= 1)))
    check(result, 'actual_curve_exposure_is_monotone', np.all(more_exposed >= linear-1e-14))
    check(result, 'zero_radiance_remains_black', np.all(zero == 0))
    check(result, 'sRGB_encoding_follows_the_checked_curve',
          np.array_equal(renderer_module.display_rgb(radiance), original_encode(linear)))
    result['scope'] = ('Chromaticity check precedes quantization. Clipping negative out-of-sRGB primaries '
                       'is the documented display-gamut limitation; the curve is not a physical emission law.')
    return _finished(result)


def verify_periodicity_and_optics():
    result = {}
    total = FPS*DURATION
    angles = np.array([inclination(n/total) for n in range(total)])
    mirrored_error = max(abs(inclination(n/total)-inclination((total-n)/total)) for n in range(total))
    check(result, 'phase_function_periodic', inclination(0) == inclination(1) and mirrored_error < 1e-12,
          maximum_mirrored_phase_error_deg=mirrored_error)
    check(result, 'inclination_scan_is_85_to_88_degrees', angles.min() == 85 and angles.max() == 88)
    renderer = OpticalRenderer(width=270, supersample=1.5, wavelength_step_nm=10)
    beginning = renderer.frames(0)
    end = renderer.frames(1)
    for mode in VIDEO_NAMES:
        check(result, f'{mode}_uncompressed_phase_zero_equals_phase_one',
              np.array_equal(np.asarray(beginning[mode]), np.asarray(end[mode])),
              frame_dimensions=[renderer.width, renderer.height])
    camera = renderer.camera
    check(result, 'overscan_covers_every_nonzero_PSF_offset',
          renderer.padding >= max(k.shape[0]//2 for k in camera.kernels),
          traced_margin_pixels=renderer.padding, largest_kernel_radius_pixels=camera.padding_pixels)
    scale_error = abs(camera.fov_rg/camera.width-FOV/renderer.width)
    check(result, 'overscan_preserves_angular_pixel_scale', scale_error < 1e-14,
          field_per_pixel_absolute_error=scale_error)

    # Trace one common set of supersampled rays. Refining the wavelength grid
    # must not change geometry or subpixel locations in this comparison.
    mapping = renderer.tracer.trace(86.5)
    temp = renderer.disk.temperature(mapping['r'])
    fcol = color_correction(temp)
    observed = mapping['g']*fcol*temp
    dilution = limb_factor(mapping['mu'])/fcol**4
    coarse = camera.apply(observed, dilution, return_ideal=True)
    accelerated = camera.apply_fast(observed, dilution, return_ideal=True)
    refined_camera = SpectralCamera(camera.width, camera.height, camera.fov_rg,
                                   aperture_m=APERTURE, observer_distance_rg=DISTANCE,
                                   wavelength_step_nm=5, support_first_nulls=SUPPORT)
    refined = refined_camera.apply(observed, dilution, return_ideal=True)
    pad = renderer.padding
    crop = (slice(pad, pad+renderer.hero_h), slice(pad, pad+renderer.width))
    inverse_matrix = np.linalg.inv(XYZ_TO_LINEAR_SRGB)
    comparisons = {}
    for mode in VIDEO_NAMES:
        a, b = coarse[mode][crop], refined[mode][crop]
        fast = accelerated[mode][crop]
        fast_error=float(np.max(np.abs(fast-a))/np.max(np.abs(a)))
        check(result,f'{mode}_accelerated_spectrum_matches_direct',fast_error<2e-5,
              maximum_RGB_error_relative_to_frame_peak=fast_error)
        ya = (a @ inverse_matrix.T)[..., 1]
        yb = (b @ inverse_matrix.T)[..., 1]
        total_y_error = float(abs(ya.sum()/yb.sum()-1))
        peak_y_error = float(np.max(np.abs(ya-yb))/yb.max())
        peak_rgb_error = float(np.max(np.abs(a-b))/np.max(np.abs(b)))
        bright = yb > yb.max()*1e-3
        local_y_p99 = float(np.percentile(np.abs(ya[bright]/yb[bright]-1), 99))
        check(result, f'{mode}_10_to_5nm_spectral_camera_convergence',
              total_y_error < .001 and peak_y_error < .001 and peak_rgb_error < .003,
              total_CIE_Y_relative_error=total_y_error,
              maximum_CIE_Y_error_relative_to_peak=peak_y_error,
              maximum_RGB_error_relative_to_peak=peak_rgb_error,
              CIE_Y_p99_relative_error_where_Y_exceeds_0_001_of_peak=local_y_p99)
        comparisons[mode] = {
            'total_CIE_Y_relative_error': total_y_error,
            'maximum_CIE_Y_error_relative_to_peak': peak_y_error,
            'maximum_RGB_error_relative_to_peak': peak_rgb_error,
            'CIE_Y_p99_relative_error_where_Y_exceeds_0_001_of_peak': local_y_p99,
        }
    result['spectral_comparison'] = comparisons
    result['scope'] = ('Exact frame periodicity and whole-image 10-versus-5-nm comparison at 270x480, '
                       'using 1.5x supersampled rays. The same finite-support pupil is used at each wavelength. '
                       'Full-resolution encoded video is checked separately.')
    del renderer, camera, refined_camera, coarse, refined
    gc.collect()
    return _finished(result)


def verify_metadata(path, selected):
    data = json.loads(path.read_text())
    result = {'file': str(path)}
    disk, ray, camera = data['disk']['parameters'], data['raytracer'], data['camera']
    composition = data['composition']
    check(result, 'rendered_physical_parameters', disk['mass_solar'] == 1e10 and
          disk['eddington_ratio'] == .05 and disk['spin_a_star'] == .9)
    check(result, 'raytracer_disk_and_spin_match_emission', ray['spin'] == disk['spin_a_star'] and
          abs(ray['rin_rg']-disk['isco_rg']) < 1e-10 and ray['rout_rg'] == 500)
    check(result, 'reel_delivery_format',
          (data['width'], data['height'], data['fps'], data['frames'], data['duration_seconds']) ==
          (1080, 1920, 30, 240, 8))
    check(result, 'scientific_viewport_and_subpixel_setting',
          composition['visible_fov_horizontal_rg'] == FOV and composition['hero_height'] == 1620 and
          composition['hero_y'] == 30 and data['supersample'] == 1.5)
    pad = composition['overscan_each_edge_pixels']
    required_pad = math.ceil(FIRST_NULL_FACTOR*830e-9/APERTURE/(FOV/(DISTANCE*data['width']))*SUPPORT)
    detector_w, detector_h = data['width']+2*pad, composition['hero_height']+2*pad
    check(result, 'overscan_contains_declared_PSF_support', pad >= required_pad and
          pad >= camera['recommended_overscan_pixels_each_edge'],
          supplied_margin_pixels=pad, independently_calculated_required_margin_pixels=required_pad)
    check(result, 'camera_detector_and_overscan_dimensions_match',
          camera['detector_width_pixels'] == detector_w and camera['detector_height_pixels'] == detector_h)
    check(result, 'traced_subpixels_match_camera',
          ray['width'] == round(detector_w*data['supersample']) and
          ray['height'] == round(detector_h*data['supersample']))
    check(result, 'overscan_field_matches_visible_pixel_angle',
          abs(camera['field_horizontal_rg']/detector_w-FOV/data['width']) < 1e-14 and
          abs(ray['fov_horizontal_rg']-camera['field_horizontal_rg']) < 1e-12)
    check(result, 'wavelength_resolved_camera_is_explicit',
          camera['aperture_m'] == APERTURE and camera['observer_distance_rg'] == DISTANCE and
          camera['wavelength_nm'] == [360, 830] and camera['wavelength_step_nm'] == 10 and
          camera['wavelength_samples'] == 48 and camera['kernel_support_first_nulls'] == SUPPORT and
          camera['astrophysical_effect'] is False)
    check(result, 'fixed_display_exposure', data['display']['exposure'] == EXPOSURE)
    acceleration=camera.get('acceleration',{})
    check(result,'spectral_acceleration_keeps_same_optical_model',
          acceleration.get('optical_model_changed') is False and
          (acceleration.get('method')=='direct fallback' or
           (acceleration.get('rank') == 18 and
            acceleration.get('temperature_domain_k', [0, 0])[1] == 100000 and
            acceleration.get('max_monochromatic_error_relative_to_continuous_Planck_peak',1)<2e-4)),
          scope='Shape error includes the direct Planck numerical-zero discontinuity near 25 K; fresh spectral_acceleration checks separately require rendered RGB and Y errors below 2e-6 of frame peak.',
          acceleration=acceleration)
    stats = data['frame_stats']
    check(result,'all_views_use_checked_ray_steps',all(
          item.get('ray_integration_step_used') in (RAY_STEP,.006) for item in stats))
    expected = np.array([inclination(n/(FPS*DURATION)) for n in range(FPS*DURATION//2+1)])
    actual = np.array([item['angle_deg'] for item in stats])
    angle_error = float(np.max(np.abs(actual-expected))) if len(actual) == len(expected) else math.inf
    check(result, 'all_121_independently_traced_views_recorded', len(stats) == 121 and angle_error < 1e-10,
          count=len(stats), maximum_expected_angle_error_deg=angle_error)
    residuals = np.array([item['constraint_residual_max'] for item in stats])
    check(result, 'every_view_has_small_geodesic_constraint_residual',
          len(residuals) == 121 and np.isfinite(residuals).all() and residuals.max() < 1e-5,
          maximum_relative_residual=float(residuals.max()))
    temperatures = np.array([item['temperature_max_k'] for item in stats])
    peak = selected['peak_effective_temperature_k']
    check(result, 'all_views_use_the_selected_thermal_solution',
          np.isfinite(temperatures).all() and np.all(temperatures <= peak*1.0001) and
          np.all(temperatures > peak*.98),
          independent_analytic_peak_k=peak,
          rendered_peak_range_k=[float(temperatures.min()), float(temperatures.max())])
    check(result, 'no_view_emits_inside_ISCO', all(item['emission_outside_isco'] for item in stats))
    check(result, 'all_recorded_frequency_shifts_physical', all(
          math.isfinite(item['frequency_shift_min']) and math.isfinite(item['frequency_shift_max']) and
          0 < item['frequency_shift_min'] <= item['frequency_shift_max'] for item in stats))
    check(result, 'all_views_contain_higher_order_disk_images', all(
          any(int(layer) > 0 and count > 0 for layer, count in item['visible_crossing_layers'].items())
          for item in stats))
    check(result, 'both_outputs_have_finite_positive_radiance_every_view', all(
          math.isfinite(item[f'{mode}_radiance_max']) and item[f'{mode}_radiance_max'] > 0
          for item in stats for mode in VIDEO_NAMES))
    result['unresolved_ray_guard'] = ('KerrTracer.trace raises for any unresolved or nonfinite ray. '
                                     'Saved per-view statistics imply that guard passed; full maps are not archived.')
    return _finished(result)


def baseline_provenance():
    result = {'scope': 'Previously saved diagnostics for unchanged modules; not rerun by this verifier.', 'reports': {}}
    for filename in ('kerr-emission-diagnostics.json', 'kerr-raytrace-diagnostics.json', 'optics-diagnostics.json'):
        path = ROOT/filename
        raw = path.read_bytes()
        data = json.loads(raw)
        result['reports'][filename] = {'passed': data.get('passed', False),
                                       'sha256': hashlib.sha256(raw).hexdigest()}
    result['passed'] = all(item['passed'] for item in result['reports'].values())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'output-optical')
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {'preflight_only': args.preflight, 'videos': {}, 'errors': []}
    from spectral_optics import diagnostics as spectral_diagnostics, fast_diagnostics
    for name, function in [('baseline_provenance', baseline_provenance),
                           ('spectral_optics', spectral_diagnostics),
                           ('spectral_acceleration', fast_diagnostics),
                           ('selected_model', verify_selected_model),
                           ('display', verify_display),
                           ('periodicity_and_optics', verify_periodicity_and_optics)]:
        print(f'Checking {name.replace("_", " ")}.', flush=True)
        try:
            report[name] = function()
        except Exception as error:
            report[name] = {'passed': False, 'error': str(error)}
    required = ['baseline_provenance', 'spectral_optics', 'spectral_acceleration', 'selected_model', 'display', 'periodicity_and_optics']
    if not args.preflight:
        required.append('render_metadata')
        try:
            report['render_metadata'] = verify_metadata(args.output/'render-metadata.json', report['selected_model'])
        except Exception as error:
            report['render_metadata'] = {'passed': False, 'error': str(error)}
        ffmpeg, ffprobe = shutil.which('ffmpeg'), shutil.which('ffprobe')
        if not ffmpeg or not ffprobe:
            report['errors'].append('FFmpeg and ffprobe must be available on PATH.')
        else:
            for mode, filename in VIDEO_NAMES.items():
                print(f'Decoding {filename}.', flush=True)
                try:
                    report['videos'][mode] = verify_video(args.output/filename, ffmpeg, ffprobe)
                except Exception as error:
                    report['videos'][mode] = {'passed': False, 'error': str(error)}
    report['passed'] = not report['errors'] and all(report[name].get('passed', False) for name in required)
    if not args.preflight:
        report['passed'] = report['passed'] and len(report['videos']) == len(VIDEO_NAMES) and all(
            item['passed'] for item in report['videos'].values())
    report['elapsed_seconds'] = time.perf_counter()-started
    destination = args.output/('verification-preflight.json' if args.preflight else 'verification.json')
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print(f'{"PASS" if report["passed"] else "FAIL"}: {destination}', flush=True)
    if not report['passed']:
        for name, section in list(report.items())+list(report['videos'].items()):
            if isinstance(section, dict):
                if 'error' in section:
                    print(f'{name}: {section["error"]}')
                for item in section.get('checks', []):
                    if not item['passed']:
                        print(f'{name}: failed {item["name"]}')
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
