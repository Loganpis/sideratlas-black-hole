#!/usr/bin/env python3
"""Verify the delivered Kerr Reels, their physical parameters and periodicity.

Run after render_kerr.py. --preflight checks numerical physics and exact
uncompressed periodicity before the full movies are ready, and writes a
separate verification-preflight.json. The default run also audits all 121
rendered views and decodes all 240 frames of each MP4. No network is used.
"""
from __future__ import annotations
import argparse
import gc
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np

from render_kerr import KerrRenderer, ROOT, FPS, DURATION, VIDEO_NAMES, inclination, EXPOSURE
from verify import check, verify_video


def verify_periodicity():
    result={}
    total=FPS*DURATION
    angles=np.array([inclination(n/total) for n in range(total)])
    mirrored=max(abs(inclination(n/total)-inclination((total-n)/total)) for n in range(total))
    result['animation']={'frames':total,'fps':FPS,'duration_seconds':DURATION,
                         'inclination_min_deg':float(angles.min()),'inclination_max_deg':float(angles.max()),
                         'maximum_mirrored_phase_angle_error_deg':mirrored,
                         'phase_zero_deg':inclination(0),'phase_one_deg':inclination(1)}
    check(result,'camera_phase_is_periodic',inclination(0)==inclination(1) and mirrored<1e-12)
    check(result,'inclination_scan_is_80_to_86_degrees',angles.min()==80 and angles.max()==86)
    renderer=KerrRenderer(width=270,supersample=1)
    beginning=renderer.frames(0)
    end=renderer.frames(1)
    for mode in VIDEO_NAMES:
        check(result,f'{mode}_uncompressed_phase_zero_equals_phase_one',
              np.array_equal(np.asarray(beginning[mode]),np.asarray(end[mode])),
              dimensions=[renderer.width,renderer.height])
    result['scope']='Exact rendered frame equality at 270x480; full-precision camera phases checked for all 240 frames. MP4 decoding is checked separately.'
    result['passed']=all(c['passed'] for c in result['checks'])
    del renderer
    gc.collect()
    return result


def verify_selected_model():
    """Check delivery parameters independently of the modules' demo defaults."""
    from kerr_emission import KerrDiskParameters,KerrEmission
    from kerr_raytrace import KerrTracer
    from optics import color_correction
    from emission import ThermalSpectrum
    parameters=KerrDiskParameters(mass_solar=1e10,eddington_ratio=.05,spin_a_star=.9)
    disk=KerrEmission(parameters)
    result={'parameters':parameters.metadata()}
    rr=np.geomspace(parameters.isco_rg+1e-4,500,10001)
    default_disk=KerrEmission(KerrDiskParameters())
    measured=disk.flux(rr)/default_disk.flux(rr)
    expected=(.05/1e10)/(.1/1e8)
    scale_error=float(np.max(np.abs(measured/expected-1)))
    check(result,'actual_mass_and_accretion_normalization',scale_error<1e-12,
          expected_flux_ratio_to_module_default=expected,maximum_relative_error=scale_error)
    dense_r=np.linspace(parameters.isco_rg,20,100001)
    dense_t=disk.temperature(dense_r)
    peak=float(dense_t.max())
    result['selected_model_peak_effective_temperature_k']=peak
    result['selected_model_peak_color_correction']=float(color_correction(peak))
    # The delivered upper inclination is tested independently, including its
    # fluid-frame angles and exact spectral transfer with the chosen mass.
    tracer=KerrTracer(384,384,spin=.9,fov=26,roll_deg=14,rout=500)
    mapping=tracer.trace(86)
    use=mapping['mask']
    temperature=disk.temperature(mapping['r'])
    shifted=mapping['g']*color_correction(temperature)*temperature
    spectrum=ThermalSpectrum()
    rgb=spectrum.lookup_rgb_linear(shifted[use])
    check(result,'upper_inclination_has_no_unresolved_rays',not np.any(mapping['status']>=3))
    check(result,'fluid_emission_angles_are_physical',np.all((mapping['mu'][use]>0)&(mapping['mu'][use]<=1)))
    check(result,'no_emission_inside_isco',np.all(mapping['r'][use]>=parameters.isco_rg-1e-6))
    check(result,'spectral_color_transform_is_finite',np.all(np.isfinite(rgb)))
    # A sufficiently cool Planck spectrum lies outside the sRGB primary
    # triangle. Its negative matrix components are a display-gamut issue,
    # not negative physical spectral intensity; the documented display curve
    # clips only those components before its common-scalar compression.
    result['sampled_srgb_out_of_gamut_fraction']=float(np.mean(np.any(rgb<0,axis=-1)))
    result['gamut_scope']='CIE-to-sRGB can yield negative primaries for cool out-of-gamut colors. The renderer clips those primaries in its documented display transform.'
    check(result,'upper_inclination_constraint_residual',mapping['constraint_residual'].max()<1e-5,
          maximum_relative_residual=float(mapping['constraint_residual'].max()))
    result['frequency_shift_range']=[float(mapping['g'][use].min()),float(mapping['g'][use].max())]
    result['observed_color_temperature_range_k']=[float(shifted[use].min()),float(shifted[use].max())]
    result['scope']='Chosen M=10^10 solar masses, normalization eta*Mdot*c^2/L_Edd=0.05, a*=0.9; demo defaults are not substituted.'
    result['passed']=all(c['passed'] for c in result['checks'])
    return result


def verify_metadata(path,selected_model):
    data=json.loads(path.read_text())
    result={'file':str(path)}
    disk=data['disk']['parameters'];ray=data['raytracer'];camera=data['camera']
    check(result,'rendered_model_parameters',disk['mass_solar']==1e10 and disk['eddington_ratio']==.05 and disk['spin_a_star']==.9,
          mass_solar=disk['mass_solar'],eddington_ratio=disk['eddington_ratio'],spin_a_star=disk['spin_a_star'])
    check(result,'raytracer_spin_and_disk_boundaries_match',
          ray['spin']==disk['spin_a_star'] and abs(ray['rin_rg']-disk['isco_rg'])<1e-10 and ray['rout_rg']==500,
          rin_rg=ray['rin_rg'],rout_rg=ray['rout_rg'])
    check(result,'reel_and_sampling_configuration',
          (data['width'],data['height'],data['fps'],data['frames'],data['duration_seconds'])==(1080,1920,30,240,8)
          and data['supersample']==1.5)
    check(result,'raytracer_matches_detector_subpixels',ray['width']==1620 and ray['height']==1410)
    check(result,'fixed_declared_display_exposure',data['display']['exposure']==EXPOSURE)
    check(result,'camera_configuration_is_explicit',camera['aperture_m']==20 and camera['wavelength_nm']==550
          and camera['observer_distance_rg']==1e6 and camera['astrophysical_effect'] is False)
    stats=data['frame_stats']
    expected=np.array([inclination(n/(FPS*DURATION)) for n in range(FPS*DURATION//2+1)])
    actual=np.array([v['angle_deg'] for v in stats])
    angle_error=float(np.max(np.abs(actual-expected))) if len(actual)==len(expected) else math.inf
    check(result,'all_121_independently_traced_views_recorded',len(stats)==121 and angle_error<1e-10,
          count=len(stats),maximum_expected_angle_error_deg=angle_error)
    residuals=np.array([v['constraint_residual_max'] for v in stats])
    check(result,'every_view_has_small_potential_constraint_residual',
          np.all(np.isfinite(residuals)) and np.max(residuals)<1e-5,
          maximum_relative_residual=float(residuals.max()))
    temps=np.array([v['maximum_local_temperature_k'] for v in stats])
    analytic_peak=selected_model['selected_model_peak_effective_temperature_k']
    check(result,'recorded_temperatures_use_selected_mass_and_accretion_rate',
          np.all(np.isfinite(temps)) and np.all(temps<=analytic_peak*1.0001) and np.all(temps>=analytic_peak*.98),
          independently_calculated_peak_k=analytic_peak,
          recorded_peak_temperature_range_k=[float(temps.min()),float(temps.max())])
    check(result,'all_views_contain_higher_order_disk_images',
          all(any(int(k)>0 and count>0 for k,count in v['visible_crossing_layers'].items()) for v in stats))
    check(result,'radiance_is_finite_and_positive',
          all(math.isfinite(v['maximum_spectral_rgb_radiance']) and v['maximum_spectral_rgb_radiance']>0 for v in stats))
    result['unresolved_ray_guard']='KerrTracer.trace raises before returning any map with an unresolved or nonfinite trajectory. Completed frame records therefore passed that guard; full per-ray maps are not saved.'
    result['atmosphere_declared']=data['atmosphere']
    result['passed']=all(c['passed'] for c in result['checks'])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'output-kerr')
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    report={'preflight_only':args.preflight,'videos':{},'errors':[]}
    print('Checking Kerr propagation, emission, atmosphere and camera optics.',flush=True)
    try:
        from kerr_raytrace import diagnostics as geodesic_diagnostics
        from kerr_emission import diagnostics as disk_diagnostics
        from emission import diagnostics as spectrum_diagnostics
        from optics import diagnostics as optics_diagnostics
        report['physics']={'raytrace':geodesic_diagnostics(),'kerr_emission':disk_diagnostics(),
                           'planck_and_colorimetry':spectrum_diagnostics(),'atmosphere_and_optics':optics_diagnostics(),
                           'passed':True}
    except Exception as error:
        report['physics']={'passed':False,'error':str(error)}
    print('Checking chosen physical parameters and exact loop periodicity.',flush=True)
    for name,function in [('selected_model',verify_selected_model),('periodicity',verify_periodicity)]:
        try: report[name]=function()
        except Exception as error: report[name]={'passed':False,'error':str(error)}
    if not args.preflight:
        try: report['render_metadata']=verify_metadata(args.output/'render-metadata.json',report['selected_model'])
        except Exception as error: report['render_metadata']={'passed':False,'error':str(error)}
        ffmpeg,ffprobe=shutil.which('ffmpeg'),shutil.which('ffprobe')
        if not ffmpeg or not ffprobe:
            report['errors'].append('FFmpeg and ffprobe must be available on PATH.')
        else:
            for mode,name in VIDEO_NAMES.items():
                print(f'Decoding and checking {name}.',flush=True)
                try: report['videos'][mode]=verify_video(args.output/name,ffmpeg,ffprobe)
                except Exception as error: report['videos'][mode]={'passed':False,'error':str(error)}
    required=['physics','selected_model','periodicity']+([] if args.preflight else ['render_metadata'])
    report['passed']=not report['errors'] and all(report.get(k,{}).get('passed',False) for k in required)
    if not args.preflight:
        report['passed']=report['passed'] and len(report['videos'])==len(VIDEO_NAMES) and all(v['passed'] for v in report['videos'].values())
    report['elapsed_seconds']=time.perf_counter()-started
    destination=args.output/('verification-preflight.json' if args.preflight else 'verification.json')
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(f'{"PASS" if report["passed"] else "FAIL"}: {destination}',flush=True)
    if not report['passed']:
        for name,section in list(report.items())+list(report['videos'].items()):
            if isinstance(section,dict):
                if 'error' in section: print(f'{name}: {section["error"]}')
                for item in section.get('checks',[]):
                    if not item['passed']: print(f'{name}: failed {item["name"]}')
    return 0 if report['passed'] else 1

if __name__=='__main__':
    raise SystemExit(main())
