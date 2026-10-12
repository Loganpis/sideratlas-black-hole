"""Integrator speed/refinement study; leaves the production settings alone."""
from pathlib import Path
import json,time
import numpy as np
from kerr_raytrace import KerrTracer
from kerr_emission import KerrDiskParameters,KerrEmission
from emission import ThermalSpectrum
from optics import atmosphere_rgb
from render_optical import display_rgb

root=Path(__file__).resolve().parent
disk=KerrEmission(KerrDiskParameters(mass_solar=1e10,eddington_ratio=.05,spin_a_star=.9))
spec=ThermalSpectrum()
settings=dict(width=540,height=764,fov=24*1300/1080,roll_deg=14,spin=.9,rout=500)
tracers={step:KerrTracer(**settings,angle_step=step) for step in (.006,.012,.018)}
results={}
start=time.perf_counter()
for angle in (85.,86.5,88.):
    then=time.perf_counter();reference=tracers[.006].trace(angle);refsecs=time.perf_counter()-then
    ref_rgb=atmosphere_rgb(disk.temperature(reference['r']),reference['g'],reference['mu'],spec)
    ref_image=display_rgb(ref_rgb).astype('int16')
    flux_ref=ref_rgb.sum((0,1))
    details={}
    for step in (.012,.018):
        then=time.perf_counter();m=tracers[step].trace(angle);secs=time.perf_counter()-then
        same=m['mask']&reference['mask']&(m['hit_layer']==reference['hit_layer'])
        rr=np.abs(m['r'][same].astype('float64')-reference['r'][same])/reference['r'][same]
        gg=np.abs(m['g'][same].astype('float64')-reference['g'][same])/reference['g'][same]
        rgb=atmosphere_rgb(disk.temperature(m['r']),m['g'],m['mu'],spec)
        diff=np.abs(rgb-ref_rgb)
        display_diff=np.abs(display_rgb(rgb).astype('int16')-ref_image)
        flux=np.abs(rgb.sum((0,1))-flux_ref)/flux_ref
        # Normalize mean pixel radiance error by mean physical radiance;
        # avoid singular relative errors on physically black pixels.
        image_l1=diff.sum((0,1))/flux_ref
        residual=float(m['constraint_residual'].max())
        mismatch=float(np.mean(m['mask']!=reference['mask']))
        layer_mismatch=float(np.mean(m['hit_layer']!=reference['hit_layer']))
        guard=residual<1e-5 and mismatch<.002 and float(np.max(image_l1))<.002 and np.percentile(rr,99)<.001
        details[str(step)]={'passed_standard_geometry_guard':bool(guard),
            'trace_seconds':secs,'speedup_vs_0_006':refsecs/secs,
            'max_potential_residual':residual,'mask_disagreement':mismatch,'layer_disagreement':layer_mismatch,
            'radius_relative_error_p99':float(np.percentile(rr,99)),
            'radius_relative_error_max':float(rr.max()),
            'g_relative_error_p99':float(np.percentile(gg,99)),
            'g_relative_error_max':float(gg.max()),
            'integrated_rgb_flux_relative_error':flux.tolist(),
            'image_rgb_L1_relative_error':image_l1.tolist(),
            'display_rgb8_absolute_error_max':int(display_diff.max()),
            'display_rgb8_absolute_error_p99_9':float(np.percentile(display_diff,99.9)),
            'fraction_pixels_with_any_rgb8_error_over_one_code':float(np.mean(display_diff.max(-1)>1)),
            'mu_range':[float(m['mu'][m['mask']].min()),float(m['mu'][m['mask']].max())],
            'higher_order_pixels':int(np.sum(m['hit_layer']>=1)),
            'highest_disk_crossing':int(m['hit_layer'].max()),
        }
        print(f'i {angle} step {step}: guard={guard}, {secs:.3f}s, residual {residual:.3g}, mask {mismatch:.3g}, L1 {np.max(image_l1):.3g}, RGB8 max {display_diff.max()}',flush=True)
    results[str(angle)]={'reference_trace_seconds':refsecs,'steps':details}
report={'camera':settings,'reference_step':.006,'angles':results,
        'scope':'Integrator-step comparison at identical 540x764 detector pixels, including physical visible RGB and the production display curve. No diffraction applied: PSF averaging would usually reduce individual-pixel differences. Not a proof for every ray of the full detector.',
        'precision_note':'r/g/mu outputs are float32. Tiny differences below this quantization cannot be resolved.',
        'elapsed_seconds':time.perf_counter()-start}
path=root/'output-optical'/'step-speed-study.json'
path.write_text(json.dumps(report,indent=2)+'\n')
print(str(path),flush=True)
