"""Independent numerical limits and source invariants for the neutral Kerr--Newman animation."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from emission import nt_flux_shape
from kn_orbits import isco_radius, orbital_quantities, horizon_radius, radial_potential_derivatives
from kerr_emission import orbital_quantities as kerr_quantities, isco_radius as kerr_isco
from kerr_orbit_emission import ShearingDiskEmission as KerrSource
from kn_orbit_emission import DimensionlessPageThorne, ShearingDiskEmission
from nasa_emission import ShearingDiskEmission as SchwarzschildSource

ROOT = Path(__file__).resolve().parent


def diagnostics():
    checks = []

    def check(name, value, passed):
        checks.append({"name": name, "value": value, "passed": bool(passed)})
        if not passed:
            raise AssertionError(f"{name}: {value}")

    source = ShearingDiskEmission()
    zero = ShearingDiskEmission(spin_a_star=0, charge_q=0)
    old = SchwarzschildSource()
    rng = np.random.default_rng(314159)
    check("Schwarzschild_ISCO", zero.rin, abs(zero.rin-6.0) < 3e-14)
    check("subextremal_chosen_parameters", .9**2+.4**2, .9**2+.4**2 < 1)
    check("horizon_rg", horizon_radius(), abs(horizon_radius()-(1+math.sqrt(.03))) < 1e-15)
    # Independent published ISCO equation, Liu, Lee & Lin2017 Eq22. The
    # production root finder instead uses a fixed-E,L radial second derivative.
    def published_isco_equation(r, a, q):
        return r*(6*r-r*r-9*q*q+3*a*a)+4*q*q*(q*q-a*a)-8*a*(r-q*q)**1.5
    residuals = []
    for a, charge in [(0,0),(.9,0),(-.9,0),(0,.4),(.9,.4),(.4,.8),(-.4,.8)]:
        ri = isco_radius(a, charge)
        orbit = orbital_quantities(ri,a,charge)
        residuals.append(abs(published_isco_equation(ri,a,charge))/ri**3)
        R, Rp, Rpp = radial_potential_derivatives(ri,orbit["energy"],orbit["angular_momentum"],a,charge)
        check(f"ISCO_R_Rprime_Rsecond_a{a}_q{charge}", [float(R),float(Rp),float(Rpp)],
              abs(R)/ri**4 < 2e-14 and abs(Rp)/ri**3 < 2e-14 and abs(Rpp)/ri**2 < 2e-14)
        for delta, stable in [(-.001,False),(.001,True)]:
            o = orbital_quantities(ri+delta,a,charge)
            curvature = float(radial_potential_derivatives(ri+delta,o["energy"],o["angular_momentum"],a,charge)[2])
            check(f"ISCO_stability_side_a{a}_q{charge}_d{delta}",curvature,(curvature<0)==stable)
        if charge == 0:
            check(f"independent_Kerr_ISCO_a{a}",ri,abs(ri-kerr_isco(a))<6e-14)
        elif a == 0:
            roots = np.roots([1,-6,9*charge**2,-4*charge**4])
            outer = max(float(x.real) for x in roots if abs(x.imag)<1e-12)
            check("Reissner_Nordstrom_outer_ISCO_cubic", ri,abs(ri-outer)<2e-14)
    check("published_ISCO_equation_scaled_residual",max(residuals),max(residuals)<3e-14)
    check("neutral_charge_sign_symmetry_ISCO",isco_radius(.9,-.4),isco_radius(.9,-.4)==source.rin)
    check("charge_moves_chosen_ISCO_inward",source.rin,source.rin<isco_radius(.9,0))
    kerr = KerrSource()
    qzero = ShearingDiskEmission(charge_q=0)
    rr = np.geomspace(kerr.rin+.001,28,2001)
    limit = float(np.max(np.abs(qzero.radial_profile(rr)/kerr.radial_profile(rr)-1)))
    check("normalized_flux_qzero_Kerr_limit",limit,limit<2e-10)
    for key, values in orbital_quantities(rr,.9,0).items():
        if key in kerr_quantities(rr,.9):
            error = float(np.max(np.abs(values-kerr_quantities(rr,.9)[key])))
            check(f"orbital_qzero_Kerr_limit_{key}",error,error<2e-13)
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
    check("KN_radial_peak_moves_inside_Schwarzschild_ISCO", source.peak_radius_rg, source.rin < source.peak_radius_rg < 6.0)

    refined = DimensionlessPageThorne(0.9, 0.4, 28, samples=131073)
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

    # Independent metric construction: radial geodesic equation solves
    # g_tt,r + 2 Omega g_tphi,r + Omega^2 g_phiphi,r = 0. Then normalize u.u=-1
    # and contract Killing vectors to obtain E,L; no closed orbit formula used.
    rr = np.geomspace(source.rin+0.01,28,1001)
    def metric_orbit(x):
        a, charge = .9,.4
        H = 2/x-charge**2/x**2
        Hp = -2/x**2+2*charge**2/x**3
        gtt, gtp, gpp = -1+H, -a*H, x*x+a*a+a*a*H
        gttp, gtpp, gppp = Hp,-a*Hp,2*x+a*a*Hp
        om = (-gtpp+np.sqrt(gtpp*gtpp-gttp*gppp))/gppp
        ut = 1/np.sqrt(-(gtt+2*om*gtp+om*om*gpp))
        return -(gtt+om*gtp)*ut,(gtp+om*gpp)*ut,om,ut
    quantities = orbital_quantities(rr,.9,.4)
    en, angmom, om, ut = metric_orbit(rr)
    for key, vals in [("energy",en),("angular_momentum",angmom),("omega",om),("ut",ut)]:
        error=float(np.max(np.abs(vals-quantities[key])))
        check("independent_metric_"+key,error,error<5e-13)
    ce,cl,co,cu = metric_orbit(rr.astype(complex)+1e-20j)
    for key, vals in [("denergy_dr",ce.imag/1e-20),("dangular_momentum_dr",cl.imag/1e-20),("domega_dr",co.imag/1e-20)]:
        error=float(np.max(np.abs(vals-quantities[key])))
        check("metric_complex_step_"+key,error,error<2e-12)
    identity=float(np.max(np.abs(quantities["denergy_dr"]-om*quantities["dangular_momentum_dr"])))
    check("circular_orbit_dE_equals_Omega_dL",identity,identity<3e-14)
    inverseut=float(np.max(np.abs(en-om*angmom-1/ut)))
    check("energy_minus_OmegaL_equals_inverse_ut",inverseut,inverseut<2e-15)
    signflip=orbital_quantities(rr,.9,-.4)
    check("neutral_charge_sign_symmetry_all_orbit_quantities",None,all(np.array_equal(signflip[k],v) for k,v in quantities.items()))
    exact_omega = np.sqrt(rr-.4**2)/(rr**2+.9*np.sqrt(rr-.4**2))
    check("exact_KN_orbital_frequency",float(np.max(np.abs(source.angular_velocity(rr)-exact_omega))),np.array_equal(source.angular_velocity(rr),exact_omega))
    lifetimes=np.minimum(1.5*2*math.pi/source.angular_velocity(source.radii),640)
    check("lifetime_uses_KN_orbital_period",float(np.max(np.abs(source.lifetimes-lifetimes))),np.array_equal(source.lifetimes,lifetimes))
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
    check("prepared_KN_omega_for_each_radius",float(np.max(np.abs(prepared.omega-source.angular_velocity(r)))),np.array_equal(prepared.omega,source.angular_velocity(r)))

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
                theta = azimuth-phase-age*math.sqrt(radius-p.charge_q**2)/(radius**2+p.spin_a_star*math.sqrt(radius-p.charge_q**2))
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

    filenames = ["kn_orbits.py", "kn_orbit_emission.py", "verify_kn_orbit_emission.py", "kerr_orbit_emission.py", "kerr_emission.py", "nasa_emission.py", "emission.py"]
    return {"status": "passed", "checks_passed": len(checks), "checks": checks,
            "metadata": source.metadata(),
            "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in filenames},
            "scope": "Neutral Kerr-Newman orbits, emission quadrature and kinematic source invariants; not plasma/GRMHD validation or a fit to NASA"}


if __name__ == "__main__":
    result = diagnostics()
    path = ROOT/"kn-orbit-emission-diagnostics.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "checks_passed": result["checks_passed"], "report": str(path), "peak_radius_rg": result["metadata"]["radial_peak_rg"], "inner_radius_rg": result["metadata"]["inner_radius_rg"]}, indent=2))
