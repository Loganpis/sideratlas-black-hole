"""Independent numerical limits and source invariants for the Kerr animation."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from emission import nt_flux_shape
from kerr_emission import isco_radius, orbital_quantities
from kerr_orbit_emission import DimensionlessPageThorne, ShearingDiskEmission
from nasa_emission import ShearingDiskEmission as SchwarzschildSource

ROOT = Path(__file__).resolve().parent


def diagnostics():
    checks = []

    def check(name, value, passed):
        checks.append({"name": name, "value": value, "passed": bool(passed)})
        if not passed:
            raise AssertionError(f"{name}: {value}")

    source = ShearingDiskEmission()
    zero = ShearingDiskEmission(spin_a_star=0)
    old = SchwarzschildSource()
    rng = np.random.default_rng(314159)
    check("prograde_ISCO_spin_0.9", source.rin, abs(source.rin-2.320883041761887) < 2e-14)
    check("Schwarzschild_ISCO", zero.rin, zero.rin == 6.0)
    check("retrograde_ISCO_spin_minus_0.9", isco_radius(-0.9), abs(isco_radius(-0.9)-8.717352279606489) < 2e-14)
    radii = np.geomspace(6.001, 28, 4001)
    # The old radial factor extracted 3/(8pi), whereas the Kerr integral
    # extracts 1/(4pi): its dimensionless value must be exactly 3/2 times old.
    zero_error = float(np.max(np.abs(zero.flux_model.flux_shape(radii)/(1.5*nt_flux_shape(radii))-1)))
    check("raw_flux_Schwarzschild_analytic_limit_max_relative_error", zero_error, zero_error < 2e-6)
    profile_error = float(np.max(np.abs(zero.radial_profile(radii)/old.radial_profile(radii)-1)))
    check("normalized_profile_Schwarzschild_limit_max_relative_error", profile_error, profile_error < 2e-6)
    r = rng.uniform(6.01, 28, 1000)
    phi = rng.uniform(0, 2*math.pi, 1000)
    t = rng.uniform(-1000, 1000, 1000)
    structured_limit = float(np.max(np.abs(zero.evaluate(r, phi, t)/old.evaluate(r, phi, t)-1)))
    check("complete_structured_source_Schwarzschild_limit_max_relative_error", structured_limit, structured_limit < 2e-6)
    check("peak_normalization", float(source.radial_profile(source.peak_radius_rg)), abs(float(source.radial_profile(source.peak_radius_rg))-1) < 1e-15)
    check("Kerr_radial_peak_moves_inside_old_ISCO", source.peak_radius_rg, source.rin < source.peak_radius_rg < 6.0)

    refined = DimensionlessPageThorne(0.9, 28, samples=131073)
    r = source.rin+np.geomspace(1e-5, 28-source.rin, 5001)
    refinement = float(np.max(np.abs(source.flux_model.flux_shape(r)/refined.flux_shape(r)-1)))
    check("flux_doubled_grid_max_relative_error", refinement, refinement < 2e-6)
    integral_errors = []
    for radius in [source.rin+0.001, 3, 5, 10, 28]:
        rr = np.linspace(source.rin, radius, 100001)
        q = orbital_quantities(rr, 0.9)
        integral = np.trapezoid(q["energy_minus_omega_l"]*q["dangular_momentum_dr"], rr)
        reference = -q["domega_dr"][-1]/(radius*q["energy_minus_omega_l"][-1]**2)*integral
        integral_errors.append(abs(float(source.flux_model.flux_shape(radius))/reference-1))
    check("independent_linear_grid_integral_max_relative_error", max(integral_errors), max(integral_errors) < 2e-6)

    # These algebraically different circular-orbit expressions are independent
    # of the production analytic derivative implementation.
    rr = np.geomspace(source.rin+0.01, 28, 1001)
    def energy_angmom(x):
        denominator = x**0.75*np.sqrt(x**1.5-3*np.sqrt(x)+2*0.9)
        return ((x**1.5-2*np.sqrt(x)+0.9)/denominator,
                (x*x-2*0.9*np.sqrt(x)+0.9**2)/denominator)
    q = orbital_quantities(rr, 0.9)
    energies, momenta = energy_angmom(rr)
    check("independent_specific_energy_formula", float(np.max(np.abs(energies-q["energy"]))), np.allclose(energies, q["energy"], rtol=0, atol=2e-15))
    check("independent_specific_angular_momentum_formula", float(np.max(np.abs(momenta-q["angular_momentum"]))), np.allclose(momenta, q["angular_momentum"], rtol=0, atol=4e-15))
    h = 1e-5*rr
    ep, lp = energy_angmom(rr+h)
    em, lm = energy_angmom(rr-h)
    derivative_error = float(np.max(np.abs((lp-lm)/(2*h)-q["dangular_momentum_dr"])))
    check("analytic_angular_momentum_derivative_vs_finite_difference", derivative_error, derivative_error < 5e-9)
    conservation = float(np.max(np.abs(q["denergy_dr"]-q["omega"]*q["dangular_momentum_dr"])))
    check("circular_orbit_dE_equals_Omega_dL", conservation, conservation < 5e-15)
    exact_omega = 1/(rr**1.5+0.9)
    check("exact_Kerr_orbital_frequency", float(np.max(np.abs(source.angular_velocity(rr)-exact_omega))), np.array_equal(source.angular_velocity(rr), exact_omega))
    check("lifetime_uses_Kerr_orbital_period", float(np.max(np.abs(source.lifetimes-np.minimum(1.5*2*math.pi*(source.radii**1.5+0.9), 640)))), np.array_equal(source.lifetimes, np.minimum(1.5*2*math.pi*(source.radii**1.5+0.9), 640)))

    radii = np.geomspace(source.rin+0.001, 28, 80)[:, None]
    phi = np.arange(4096)[None, :]*(2*math.pi/4096)
    mean_errors = []
    for time in [0, 79.31, 311.71, 799.999]:
        values = source.evaluate(radii, phi, time)
        mean_errors.append(float(np.max(np.abs(values.mean(axis=1)/source.radial_profile(radii[:, 0])-1))))
    check("structured_azimuthal_mean_equals_Page_Thorne", max(mean_errors), max(mean_errors) < 1e-12)
    r = rng.uniform(source.rin+0.001, 28, 8000)
    phi = rng.uniform(-3*math.pi, 3*math.pi, 8000)
    t = rng.uniform(-1600, 1600, 8000)
    values = source.evaluate(r, phi, t)
    periodic = float(np.max(np.abs(values-source.evaluate(r, phi, t+800))))
    check("source_birth_periodicity", periodic, periodic < 2e-12)
    check("finite_positive_emission", [float(np.min(values)), float(np.max(values))], np.all(np.isfinite(values)) and np.all(values >= 0))
    delay = rng.uniform(-80, 250, 8000)
    prepared = source.prepare(r, phi, delay)
    observer_time = 137.5
    retarded = float(np.max(np.abs(prepared.evaluate(observer_time)-source.evaluate(r, phi, observer_time-delay))))
    check("retarded_time_sign_and_broadcasting", retarded, retarded == 0)
    check("prepared_Kerr_omega_for_each_radius", float(np.max(np.abs(prepared.omega-1/(r**1.5+0.9)))), np.array_equal(prepared.omega, 1/(r**1.5+0.9)))

    # Independent scalar transcription of the declared source equation tests
    # the actual retarded phase, not two aliases of the same evaluator.
    def scalar(radius, azimuth, emission_time):
        numerator = denominator = 1.0
        p = source.parameters
        for j, center in enumerate(source.log_centers):
            x = (math.log(radius)-center)/p.radial_log_half_width
            if abs(x) >= 1:
                continue
            radial_weight = math.exp(-2*x*x/(1-x*x))
            for birth, phase in zip(source.births[j], source.phases[j]):
                age = (emission_time-birth) % p.loop_period_tg
                if age >= source.lifetimes[j]:
                    continue
                weight = p.contrast*radial_weight*math.sin(math.pi*age/source.lifetimes[j])**4
                theta = azimuth-phase-age/(radius**1.5+p.spin_a_star)
                numerator += weight*math.exp(p.azimuthal_kappa*math.cos(theta))/np.i0(p.azimuthal_kappa)
                denominator += weight
        return float(source.radial_profile(radius))*numerator/denominator
    errors = []
    for time in [-11.3, 0.0, 799.999, 823.1]:
        reference = np.array([scalar(float(r[i]), float(phi[i]), time-float(delay[i])) for i in range(50)])
        errors.append(float(np.max(np.abs(prepared.evaluate(time)[:50]-reference))))
    check("independent_scalar_retarded_source_equation", max(errors), max(errors) < 2e-12)
    masks = source.evaluate([0, source.rin, 28.1, np.nan, np.inf], 0, 0).tolist()
    check("nonemitting_invalid_radii_zero", masks, masks == [0.0]*5)
    smooth_error = float(np.max(np.abs(prepared.evaluate(137.5, structured=False)-source.radial_profile(r))))
    check("smooth_mode_exact_radial_profile", smooth_error, smooth_error == 0)
    deltas = []
    for dt in [1e-3, 5e-4]:
        f0 = source.evaluate(r, phi, 0)
        left = (f0-source.evaluate(r, phi, -dt))/dt
        right = (source.evaluate(r, phi, dt)-f0)/dt
        deltas.append(float(np.max(np.abs(left-right))))
    ratio = deltas[1]/deltas[0]
    check("seam_derivative_continuity_halving_ratio", ratio, 0.45 < ratio < 0.55)
    turns = float(source.angular_velocity(10)*800/(2*math.pi))
    check("orbital_frequency_not_quantized_to_loop", turns, abs(turns-round(turns)) > 0.01)
    no_structure = ShearingDiskEmission(contrast=0)
    check("zero_contrast_gives_radial_profile", float(np.max(np.abs(no_structure.evaluate(r, phi, t)-no_structure.radial_profile(r)))), np.array_equal(no_structure.evaluate(r, phi, t), no_structure.radial_profile(r)))

    filenames = ["kerr_orbit_emission.py", "verify_kerr_orbit_emission.py", "kerr_emission.py", "nasa_emission.py", "emission.py"]
    return {"status": "passed", "checks_passed": len(checks), "checks": checks,
            "metadata": source.metadata(),
            "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in filenames},
            "scope": "Kerr emission quadrature and kinematic source invariants; not GRMHD validation or a fit to NASA"}


if __name__ == "__main__":
    result = diagnostics()
    path = ROOT/"kerr-orbit-emission-diagnostics.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "checks_passed": result["checks_passed"], "report": str(path), "peak_radius_rg": result["metadata"]["radial_peak_rg"], "inner_radius_rg": result["metadata"]["inner_radius_rg"]}, indent=2))
