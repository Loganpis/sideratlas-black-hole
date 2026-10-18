"""Source invariants and independent scalar checks for v1.2 emission.

These test the declared phenomenological source and its characteristic
advection, not GRMHD, plasma heating, or agreement with an observed disk.
Emitter normalization and geodesic dynamics are checked in verify_kn_flow.py.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from kn_dynamic_emission import DynamicDiskEmission, DynamicParameters
from kn_flow import NeutralAccretionFlow
from kn_orbit_emission import ShearingDiskEmission as PreviousSource

ROOT = Path(__file__).resolve().parent


def diagnostics():
    checks = []

    def check(name, value, passed):
        checks.append({"name": name, "value": value, "passed": bool(passed)})
        if not passed:
            raise AssertionError(f"{name}: {value}")

    source = DynamicDiskEmission()
    old = PreviousSource()
    p = source.parameters
    rng = np.random.default_rng(120013326)
    r = source.rin+np.geomspace(1e-7, 28-source.rin, 2048)
    error = float(np.max(np.abs(source.radial_profile(r)-old.radial_profile(r))))
    check("unchanged_Page_Thorne_exterior_radial_background", error, error < 1e-14)
    check("unchanged_Page_Thorne_peak_radius", source.peak_radius_rg,
          source.peak_radius_rg == old.peak_radius_rg)
    # Source-time azimuthal means are distinct from observer-time images,
    # whose phi-dependent propagation delays and transfer need not average
    # to the radial background.
    rmean = np.concatenate((np.linspace(source.flow.inner_radius_rg+1e-5,
        source.rin-.0001, 12), np.geomspace(source.rin+.0001, 28, 24)))[:, None]
    phi = (2*math.pi*np.arange(1024)/1024)[None, :]
    mean_error = 0.
    for t in (-103.4, 0., 11.7, 24.99999, 399.2, 799.999):
        values = source.evaluate(rmean, phi, t)
        mean_error = max(mean_error, float(np.max(np.abs(
            values.mean(axis=1)/source.radial_profile(rmean[:, 0])-1))))
    check("exact_emitted_azimuthal_mean_at_each_radius_and_time", mean_error, mean_error < 2e-13)
    rr = rng.uniform(source.flow.inner_radius_rg, 28, 12000)
    ph = rng.uniform(-10, 10, rr.size)
    tt = rng.uniform(-2400, 2400, rr.size)
    emission = source.evaluate(rr, ph, tt)
    background = source.radial_profile(rr)
    modulation = emission/background
    check("positive_finite_intensity", [float(emission.min()), float(emission.max())],
          np.all(np.isfinite(emission)) and np.all(emission >= 0))
    check("analytically_bounded_contrast", [float(modulation.min()), float(modulation.max())],
          modulation.min() >= source.minimum_modulation-1e-14
          and modulation.max() <= source.kernel_maximum+1e-14)
    delay = rng.uniform(-80, 500, rr.size)
    prepared = source.prepare(rr, ph, delay)
    t_obs = 117.15
    difference = float(np.max(np.abs(prepared.evaluate(t_obs)-source.evaluate(rr, ph, t_obs-delay))))
    check("per_ray_retarded_time_sign_and_broadcasting", difference, difference == 0)
    difference = float(np.max(np.abs(emission-source.evaluate(rr, ph, tt+p.loop_period_tg))))
    check("exact_constructed_birth_periodicity", difference, difference < 3e-12)
    difference = float(np.max(np.abs(prepared.evaluate(127, structured=False)-background)))
    check("unstructured_mode_is_radial_background", difference, difference == 0)
    zero = DynamicDiskEmission(contrast=0)
    check("zero_contrast_disables_perturbations", None,
          np.array_equal(zero.evaluate(rr, ph, tt), zero.radial_profile(rr)))
    check("multiscale_packet_widths_and_irregular_birth_centers", None,
          source.packet_widths.min() < .06 and source.packet_widths.max() > .35
          and not np.array_equal(source.packet_centers[0], source.packet_centers[1]))

    # Independent scalar transcription of the documented source equation.
    # It does not call prepare/evaluate/_birth_field or copy its field arrays.
    def scalar(radius, azimuth, emitted_time):
        interval = p.loop_period_tg/p.births_per_loop
        t = emitted_time % p.loop_period_tg
        index = math.floor(t/interval)
        age = t-index*interval
        current_weight = math.sin(.5*math.pi*age/interval)**2
        material = source.flow.streamline(radius)
        field = 0.
        def smooth(x):
            x = min(max(x, 0.), 1.)
            return x*x*(3-2*x)
        def kernel(theta, bandwidth, cap):
            ks = {n: (1+math.cos(theta))**n/(math.comb(2*n, n)/2**n)
                  for n in (1, 2, 4, 8, 12)}
            w2 = smooth((bandwidth-2)/.25)
            w4 = smooth((bandwidth-4)/.6)
            if cap <= 4:
                middle = (1-w4)*ks[2]+w4*ks[4]
            else:
                w8, w12 = smooth((bandwidth-8)/1.2), smooth((bandwidth-12)/2)
                high = (1-w12)*ks[8]+w12*ks[12]
                middle = (1-w4)*ks[2]+w4*((1-w8)*ks[4]+w8*high)
            return (1-w2)*ks[1]+w2*middle
        for birth_age, birth_index, weight in [(age, index, current_weight),
                (age+interval, (index-1) % p.births_per_loop, 1-current_weight)]:
            upstream_age = float(material["age_tg"])-birth_age
            r0 = float(source.flow.radius_at_age(upstream_age))
            a0 = float(source.flow.azimuth_at_age(upstream_age))
            phi0 = (azimuth-float(material["azimuth_rad"])+a0) % (2*math.pi)
            band = float(np.interp(upstream_age, source._bandwidth_age, source._bandwidth))
            spiral = 1+source.spiral_contrast*smooth((band-2)/.25)*math.cos(
                2*phi0-3*math.log(r0/source.rin)+float(source.spiral_phases[birth_index]))
            numerator = 1+p.contrast*source.spiral_weight*spiral
            denominator = 1+p.contrast*source.spiral_weight
            for j in range(len(source.packet_caps)):
                center, width = source.packet_centers[birth_index, j], source.packet_widths[birth_index, j]
                x = (math.log(r0)-center)/width
                if abs(x) >= 1:
                    continue
                radial = math.exp(1-1/(1-x*x))
                w = p.contrast*source.packet_amplitudes[birth_index, j]*radial
                theta = phi0+source.packet_slopes[birth_index, j]*(math.log(r0)-center)+source.packet_phases[birth_index, j]
                numerator += w*kernel(theta, band, int(source.packet_caps[j]))
                denominator += w
            field += weight*numerator/denominator
        return float(source.radial_profile(radius))*field
    indices = np.arange(80)
    scalar_reference = np.array([scalar(float(rr[i]), float(ph[i]), float(t_obs-delay[i])) for i in indices])
    error = float(np.max(np.abs(prepared.evaluate(t_obs)[indices]-scalar_reference)))
    check("independent_scalar_source_equation", error, error < 3e-13)

    # Follow the same material element forward, including across the ISCO;
    # its backward-traced birth position and azimuth must remain unchanged.
    radii = np.array([source.rin+.003, source.rin+.0001, source.rin-.003,
                      (source.rin+source.horizon_rg)/2, 2., 5., 20.])
    c0 = source.flow.streamline(radii)
    dt = 3.7
    next_radii = source.flow.radius_at_age(c0["age_tg"]+dt)
    c1 = source.flow.streamline(next_radii)
    phi0 = np.linspace(-.5, 1., len(radii))
    phi1 = phi0+c1["azimuth_rad"]-c0["azimuth_rad"]
    birth_age = 17.3
    rb0 = source.flow.radius_at_age(c0["age_tg"]-birth_age)
    rb1 = source.flow.radius_at_age(c1["age_tg"]-(birth_age+dt))
    pb0 = phi0-c0["azimuth_rad"]+source.flow.azimuth_at_age(c0["age_tg"]-birth_age)
    pb1 = phi1-c1["azimuth_rad"]+source.flow.azimuth_at_age(c1["age_tg"]-birth_age-dt)
    radial_error, angle_error = float(np.max(np.abs(rb0-rb1))), float(np.max(np.abs(pb0-pb1)))
    check("exact_backward_material_birth_radius_including_ISCO_crossing", radial_error, radial_error < 2e-13)
    check("exact_backward_material_birth_azimuth_including_ISCO_crossing", angle_error, angle_error < 3e-11)
    check("radial_advection_is_inward", next_radii.tolist(), np.all(next_radii < radii))
    # The entire newborn modulation, including its radial weights and
    # angular-resolution choice, is material, not only its phase. Following
    # the same parcel leaves it unchanged before the separate birth blend.
    upstream0 = source.prepare(radii, phi0)
    upstream1 = source.prepare(next_radii, phi1)
    same_birth = np.full(len(radii), 7, dtype=np.int64)
    m0 = upstream0._birth_modulation(np.full(len(radii), birth_age), same_birth)
    m1 = upstream1._birth_modulation(np.full(len(radii), birth_age+dt), same_birth)
    material_error = float(np.max(abs(m0-m1)))
    check("whole_packet_modulation_is_materially_advected", material_error, material_error < 3e-10)
    permutation = rng.permutation(len(rr))
    shuffled = source.prepare(rr[permutation], ph[permutation], delay[permutation]).evaluate(t_obs)
    permutation_error = float(np.max(abs(shuffled-prepared.evaluate(t_obs)[permutation])))
    check("radius_sorted_preparation_preserves_original_event_permutation",
          permutation_error, permutation_error == 0)
    check("large_prepared_arrays_use_monotone_radius_order", None,
          np.all(np.diff(prepared.radius) >= 0))

    angles = 2*math.pi*np.arange(4096)/4096
    kernel_prepared = source.prepare([source.peak_radius_rg], [0.])
    kernel_mean_error, kernel_minimum, leakage = 0., 1., 0.
    for band in [1.8, 2., 2.125, 2.25, 3.8, 4.1, 4.6, 8.4, 9.2, 12.5, 20.]:
        kernel_values = kernel_prepared._kernel(angles, np.full(angles.shape, band), 12)
        kernel_mean_error = max(kernel_mean_error, abs(float(kernel_values.mean())-1))
        kernel_minimum = min(kernel_minimum, float(kernel_values.min()))
        spectrum = abs(np.fft.rfft(kernel_values))/len(angles)
        highest_allowed = max(1, min(12, math.floor(band)))
        leakage = max(leakage, float(spectrum[highest_allowed+1:].max()))
    check("positive_kernel_blending_preserves_analytic_unit_angular_mean",
          kernel_mean_error, kernel_mean_error < 2e-15)
    check("positive_kernel_blending_never_introduces_negative_emission",
          kernel_minimum, kernel_minimum >= 0)
    check("finite_angular_bandwidth_respects_local_advection_guide",
          leakage, leakage < 2e-14)
    # These are explicit presentation metrics at one common emissivity
    # normalization; they are not claims about an observed disk's variability.
    times = [0., 12.5, 137.3, 366.7]
    textures = np.stack([source.evaluate(source.peak_radius_rg, angles, t)
        /source.radial_profile(source.peak_radius_rg) for t in times])
    emitted_rms = textures.std(axis=1)
    levels = (textures*.7**4/.02)/(1+textures*.7**4/.02)
    level_rms = 255*levels.std(axis=1)
    check("restored_resolved_inner_emission_richness", emitted_rms.tolist(),
          float(emitted_rms.mean()) > .30 and float(emitted_rms.max()) < 1.)
    check("texture_survives_fixed_monotone_tone_mapping", level_rms.tolist(),
          float(level_rms.mean()) > 7.)
    shifted = np.stack([source.evaluate(source.peak_radius_rg*math.exp(.04), angles, t)
        /source.radial_profile(source.peak_radius_rg*math.exp(.04)) for t in times])
    radial_structure = float(np.sqrt(np.mean((shifted-textures)**2)))
    check("resolved_nonaxisymmetric_radial_fine_structure", radial_structure, radial_structure > .08)
    turns = float(source.angular_velocity(10)*800/(2*math.pi))
    check("gas_frequency_not_rounded_to_close_loop", turns, abs(turns-round(turns)) > .01)

    # Every sin^2/cos^2 birth boundary, including the 800tg seam, is C1.
    small_r = rng.uniform(source.flow.inner_radius_rg+.001, 12., 120)
    small_phi = rng.uniform(-math.pi, math.pi, len(small_r))
    derivative_errors = []
    for dt in [1e-3, 5e-4]:
        worst = 0.
        for seam in [0., 25., 50., 375., 775., 800.]:
            center = source.evaluate(small_r, small_phi, seam)
            left = (center-source.evaluate(small_r, small_phi, seam-dt))/dt
            right = (source.evaluate(small_r, small_phi, seam+dt)-center)/dt
            worst = max(worst, float(np.max(np.abs(left-right))))
        derivative_errors.append(worst)
    ratio = derivative_errors[1]/derivative_errors[0]
    check("C1_birth_and_loop_seam_refinement_ratio", ratio, .43 < ratio < .57)

    rh, ri = source.horizon_rg, source.rin
    inner_r = rh+(ri-rh)*np.array([.01, .1, .5, .9, .99])
    expected = p.plunge_peak_fraction*(1-(inner_r-rh)/(ri-rh))**2
    error = float(np.max(np.abs(source.radial_profile(inner_r)-expected)))
    check("declared_weak_inner_heating_equation", error, error < 1e-16)
    check("finite_nonzero_intrinsic_exterior_horizon_limit", float(source.radial_profile(rh+1e-7)),
          abs(float(source.radial_profile(rh+1e-7))/p.plunge_peak_fraction-1) < 1e-6)
    check("horizon_and_interior_source_black", source.radial_profile([rh-1e-3, rh, ri]).tolist(),
          np.array_equal(source.radial_profile([rh-1e-3, rh, ri]), np.zeros(3)))
    edge_offsets = np.array([1e-5, 5e-6])
    edge_values = source.radial_profile(ri-edge_offsets)
    check("inner_zero_value_and_slope_at_ISCO", (edge_values[1]/edge_values[0]).item(),
          abs(edge_values[1]/edge_values[0]-.25) < 1e-10)
    dark = DynamicDiskEmission(plunge_peak_fraction=0)
    check("inner_heating_can_be_disabled", None,
          np.all(dark.evaluate(inner_r, .1, 0) == 0))
    bad = [0., rh, source.flow.inner_radius_rg-1e-8, 28.1, np.nan, np.inf]
    check("invalid_and_nonemitting_source_radii_zero", source.evaluate(bad, 0., 0.).tolist(),
          np.array_equal(source.evaluate(bad, 0., 0.), np.zeros(len(bad))))
    reject = 0
    for kwargs in ({"contrast": 21.}, {"contrast": -.1}, {"plunge_peak_fraction": -.01},
                   {"births_per_loop": 3}, {"loop_period_tg": 0.}):
        try:
            DynamicParameters(**kwargs)
        except ValueError:
            reject += 1
    check("invalid_source_parameters_rejected", reject, reject == 5)
    try:
        DynamicDiskEmission(flow=NeutralAccretionFlow(outer_radius_rg=28))
        upstream_rejected = False
    except ValueError:
        upstream_rejected = True
    check("upstream_domain_clipping_not_silently_allowed", upstream_rejected, upstream_rejected)
    try:
        prepared.evaluate(np.nan)
        invalid_time_rejected = False
    except ValueError:
        invalid_time_rejected = True
    check("nonfinite_observer_time_rejected", invalid_time_rejected, invalid_time_rejected)
    check("phenomenological_physics_disclosed", None,
          not source.metadata()["plunge_emissivity_is_calculated_from_thermodynamics"]
          and "GRMHD turbulence or MRI" in source.metadata()["not_modeled"])
    # Timing is diagnostic only: host scheduling must not determine physics
    # test outcomes. It also documents the hit-only source cost for rendering.
    start = time.perf_counter()
    for t in (0., 15., 71.):
        prepared.evaluate(t)
    seconds_per_12000_hits = (time.perf_counter()-start)/3
    filenames = ["kn_dynamic_emission.py", "verify_kn_dynamic_emission.py",
                 "kn_flow.py", "kn_orbit_emission.py", "kn_orbits.py"]
    return {"status": "passed", "checks_passed": len(checks), "checks": checks,
        "metadata": source.metadata(), "seconds_per_12000_hits": seconds_per_12000_hits,
        "source_sha256": {n: hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in filenames},
        "scope": "Declared bounded emission, exact exterior azimuthal mean, characteristic advection and retarded times; not GRMHD or plasma-heating validation"}


if __name__ == "__main__":
    result = diagnostics()
    path = ROOT/"kn-dynamic-emission-diagnostics.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "checks_passed": result["checks_passed"],
                      "report": str(path), "seconds_per_12000_hits": result["seconds_per_12000_hits"]}, indent=2))
