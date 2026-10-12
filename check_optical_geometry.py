"""Independent small-detector refinement checks for the final optical scan."""
from pathlib import Path
import json
import time
import numpy as np
from kerr_raytrace import KerrTracer, isco_radius
from kerr_emission import KerrDiskParameters, KerrEmission
from emission import ThermalSpectrum
from optics import atmosphere_rgb, color_correction

root=Path(__file__).resolve().parent
output=root/'output-optical'/'geometry-checks.json'
output.parent.mkdir(exist_ok=True)
start=time.perf_counter()
parameters=KerrDiskParameters(mass_solar=1e10,eddington_ratio=.05,spin_a_star=.9)
disk=KerrEmission(parameters)
spectrum=ThermalSpectrum()
base=dict(width=384,height=576,fov=24,roll_deg=14,spin=.9,rout=500)
coarse=KerrTracer(**base,angle_step=.006)
fine=KerrTracer(**base,angle_step=.003)
results={}
all_pass=True
for angle in (85.,86.5,88.):
    c=coarse.trace(angle);f=fine.trace(angle)
    same=c['mask']&f['mask']&(c['hit_layer']==f['hit_layer'])
    mask=c['mask']
    rr=np.abs(c['r'][same].astype('float64')-f['r'][same])/f['r'][same]
    gg=np.abs(c['g'][same].astype('float64')-f['g'][same])/f['g'][same]
    mm=np.abs(c['mu'][same].astype('float64')-f['mu'][same])
    temp=disk.temperature(c['r']);tf=disk.temperature(f['r'])
    rgb=atmosphere_rgb(temp,c['g'],c['mu'],spectrum)
    rgbfine=atmosphere_rgb(tf,f['g'],f['mu'],spectrum)
    total=rgb.sum(axis=(0,1));totalfine=rgbfine.sum(axis=(0,1))
    rel_flux=np.abs(total-totalfine)/np.maximum(np.abs(totalfine),1e-30)
    observed_temp=temp*c['g']*color_correction(temp)
    mask_diff=float(np.mean(c['mask']!=f['mask']))
    crossing_diff=float(np.mean(c['hit_layer']!=f['hit_layer']))
    finite=all(np.isfinite(c[k]).all() and np.isfinite(f[k]).all() for k in ('r','g','mu','constraint_residual'))
    criteria={
        'finite_geometry_and_spectrum':bool(finite and np.isfinite(rgb).all() and np.isfinite(rgbfine).all()),
        'all_rays_resolved':bool(np.max(c['status'])<3 and np.max(f['status'])<3),
        'emitting_radius_inside_disk':bool(np.all(c['r'][mask]>=isco_radius(.9)-1e-6) and np.all(c['r'][mask]<=500+1e-4)),
        'emission_angle_in_physical_range':bool(np.all(c['mu'][mask]>0) and np.all(c['mu'][mask]<=1.000001)),
        'positive_frequency_shift':bool(np.all(c['g'][mask]>0)),
        'mask_refinement_below_0_2_percent':mask_diff<.002,
        'radius_refinement_p99_below_0_1_percent':bool(np.percentile(rr,99)<.001),
        'frequency_refinement_p99_below_0_1_percent':bool(np.percentile(gg,99)<.001),
        'integrated_spectral_rgb_refinement_below_0_2_percent':bool(np.max(rel_flux)<.002),
        'potential_constraint_below_1e_5':bool(max(c['constraint_residual'].max(),f['constraint_residual'].max())<1e-5),
        'at_least_one_higher_order_disk_image':bool(np.any(c['hit_layer'][mask]>=1)),
        'zero_radiance_for_non_emitting_pixels':bool(np.all(rgb[~mask]==0)),
    }
    passed=all(criteria.values());all_pass &=passed
    result={
        'passed':passed,'criteria':criteria,
        'emitting_fraction':float(mask.mean()),'shadow_fraction':float(c['shadow'].mean()),
        'mask_disagreement_fraction':mask_diff,'crossing_layer_disagreement_fraction':crossing_diff,
        'radius_relative_error':{'p99':float(np.percentile(rr,99)),'max':float(rr.max())},
        'frequency_shift_relative_error':{'p99':float(np.percentile(gg,99)),'max':float(gg.max())},
        'mu_absolute_error':{'p99':float(np.percentile(mm,99)),'max':float(mm.max())},
        'integrated_visible_rgb_relative_error':rel_flux.tolist(),
        'potential_constraint_max':{'step_0_006':float(c['constraint_residual'].max()),'step_0_003':float(f['constraint_residual'].max())},
        'radius_rg_range':[float(c['r'][mask].min()),float(c['r'][mask].max())],
        'frequency_shift_range':[float(c['g'][mask].min()),float(c['g'][mask].max())],
        'emission_mu_range':[float(c['mu'][mask].min()),float(c['mu'][mask].max())],
        'effective_temperature_k_range':[float(temp[mask].min()),float(temp[mask].max())],
        'observed_color_temperature_k_range':[float(observed_temp[mask].min()),float(observed_temp[mask].max())],
        'visible_crossing_layers':{str(k):int((c['hit_layer'][mask]==k).sum()) for k in np.unique(c['hit_layer'][mask])},
    }
    results[str(angle)]=result
    print(f'i={angle}: passed={passed}, mask disagreement={mask_diff:.6g}, RGB relative error={np.max(rel_flux):.6g}, residual={c["constraint_residual"].max():.6g}',flush=True)
report={'passed':all_pass,'camera':base,'steps':[.006,.003],'angles':results,
        'precision_note':'The kernel integrates in float64 but returns r, g and mu in float32. Zero refinement differences mean agreement at float32 output precision, not exact continuum solutions.',
        'scope':'Integrator refinement at fixed 384x576 pixels; not a detector-resolution or physical-model convergence test.',
        'elapsed_seconds':time.perf_counter()-start}
output.write_text(json.dumps(report,indent=2)+'\n')
print(str(output),flush=True)
if not all_pass: raise SystemExit(1)
