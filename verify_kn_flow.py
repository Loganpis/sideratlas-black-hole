"""Independent metric, limiting-case and transport checks for kn_flow.py.

Run ``python3 verify_kn_flow.py``.  These are numerical consistency checks
of a neutral test-flow model, not validation against a fluid simulation.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from kn_flow import FlowParameters, NeutralAccretionFlow, equatorial_metric
from kn_orbits import orbital_quantities

ROOT = Path(__file__).resolve().parent


def diagnostics():
    checks = []

    def check(name, value, passed):
        checks.append({"name": name, "value": value, "passed": bool(passed)})
        if not passed:
            raise AssertionError(f"{name}: {value}")

    # The general metric is constructed independently in factored ds^2 form.
    def independent_metric(r, a, q):
        delta = r*r-2*r+a*a+q*q
        g = np.zeros(np.shape(r)+(3, 3))
        g[..., 0, 0] = (a*a-delta)/(r*r)
        g[..., 0, 2] = g[..., 2, 0] = a*(delta-r*r-a*a)/(r*r)
        g[..., 1, 1] = r*r/delta
        g[..., 2, 2] = ((r*r+a*a)**2-a*a*delta)/(r*r)
        return g

    flow = NeutralAccretionFlow()
    p = flow.parameters
    check("selected_model_subextremal", p.spin_a_star**2+p.charge_q**2, p.spin_a_star**2+p.charge_q**2 < 1)
    check("nonzero_inward_launch_speed", flow.launch_radial_speed, -.01 < flow.launch_radial_speed < 0)
    check("finite_positive_ISCO_crossing_time", float(flow.streamline(flow.inner_radius_rg)["age_tg"]),
          0 < float(flow.streamline(flow.inner_radius_rg)["age_tg"]) < 1000)

    for a, q in [(.9,.4), (.9,0), (0,0), (0,.4), (-.5,.2)]:
        f = NeutralAccretionFlow(spin_a_star=a, charge_q=q, radial_samples=8193)
        label = f"a{a}_q{q}"
        rr = f.horizon_rg+np.geomspace(.001, 28-f.horizon_rg, 2001)
        v = f.four_velocity(rr)
        g = independent_metric(rr, a, q)
        u = np.stack((v["ut"], v["ur"], v["uphi"]), axis=-1)
        norm = np.einsum("...i,...ij,...j->...", u, g, u)
        error = float(np.max(np.abs(norm+1)))
        check(f"independent_metric_timelike_normalization_{label}", error, error < 3e-8)
        check(f"future_directed_inward_flow_{label}", [float(v["ut"].min()), float(v["ur"].max())],
              np.all(v["ut"] > 0) and np.all(v["ur"] < 0))
        E = -(g[...,0,0]*v["ut"]+g[...,0,2]*v["uphi"])
        L = g[...,2,0]*v["ut"]+g[...,2,2]*v["uphi"]
        error = float(max(np.max(abs(E-v["energy"])), np.max(abs(L-v["angular_momentum"]))))
        check(f"independent_Killing_vector_energy_momentum_{label}", error, error < 2e-11)
        inner = rr < f.isco_rg
        check(f"conserved_geodesic_E_and_L_in_plunge_{label}",
              [float(np.ptp(E[inner])), float(np.ptp(L[inner]))],
              np.ptp(E[inner]) < 2e-11 and np.ptp(L[inner]) < 3e-11)
        # Matrix inversion is independent of the separated first integrals.
        covariant = np.stack((-v["energy"], v["angular_momentum"]), axis=-1)
        temporal_block = g[..., (0,0,2,2), (0,2,0,2)].reshape((-1,2,2))
        solved = np.linalg.solve(temporal_block, covariant[...,None])[...,0]
        # u^phi changes sign for a retrograde hole; a unit-added scale avoids
        # dividing by the real zero where an emitter momentarily has Omega=0.
        err = float(max(np.max(abs(solved[...,0]-v["ut"])/(1+abs(v["ut"]))),
                        np.max(abs(solved[...,1]-v["uphi"])/(1+abs(v["uphi"])))))
        check(f"independent_metric_inverse_velocity_{label}", err, err < 1e-11)
        # Direct quartic potential versus the stable triple-root construction.
        r = np.linspace(f.horizon_rg+.0001, f.isco_rg, 2001)
        Ee, Ll = f.plunge_energy, f.plunge_angular_momentum
        delta = r*r-2*r+a*a+q*q
        raw = (Ee*(r*r+a*a)-a*Ll)**2-delta*(r*r+(Ll-a*Ee)**2)
        diff = float(np.max(abs(raw-f.plunge_potential(r)))/(1+np.max(abs(raw))))
        check(f"independent_timelike_radial_quartic_{label}", diff, diff < 3e-13)
        check(f"plunge_has_no_radial_turning_point_{label}", float(f.plunge_potential(r).min()),
              np.all(f.plunge_potential(r) > 0))
        # Compare both sides to the common launch, without a zero-speed stall.
        vi = f.four_velocity(f.isco_rg)
        vl = f.four_velocity(f.isco_rg-1e-9)
        err = float(max(abs(float(vl[k]/vi[k])-1) for k in ["ut","ur","uphi"]))
        check(f"continuous_ISCO_four_velocity_{label}", err, err < 3e-8)
        # Finite radial perturbations must tend to the hole's horizon angular
        # velocity. This limit follows from horizon regularity, independently
        # of the launch energy or angular momentum.
        rh = f.horizon_rg
        near = f.four_velocity(rh+1e-7)["omega"]
        omega_h = a/(rh*rh+a*a)
        check(f"frame_dragging_horizon_angular_velocity_{label}", float(abs(near-omega_h)), abs(near-omega_h) < 1e-7)

    # Schwarzschild zero-drift ISCO plunge has a closed radial expression:
    # u^r=-sqrt(2/(3*rI))*(rI/r-1)^(3/2), rI=6.  A tiny finite launch
    # recovers it away from the zero-speed ISCO endpoint.
    zero = NeutralAccretionFlow(spin_a_star=0, charge_q=0,
                                alpha_viscosity=1e-7, radial_samples=8193)
    r = np.linspace(2.01, 5.9, 500)
    analytic = -(1/3)*(6/r-1)**1.5
    err = float(np.max(abs(zero.four_velocity(r)["ur"]/analytic-1)))
    check("Schwarzschild_analytic_ISCO_plunge_radial_velocity_limit", err, err < 2e-12)
    check("Schwarzschild_analytic_ISCO_energy_limit", zero.plunge_energy, abs(zero.plunge_energy-math.sqrt(8/9)) < 2e-15)
    check("Schwarzschild_analytic_ISCO_momentum_limit", zero.plunge_angular_momentum, abs(zero.plunge_angular_momentum-math.sqrt(12)) < 2e-14)

    # Independent radial geodesic equation from metric derivatives and finite
    # differences of u^r; no separated-potential derivatives are used here.
    r = np.linspace(flow.horizon_rg+.015, flow.isco_rg-.005, 1001)
    v = flow.four_velocity(r)
    a, q = p.spin_a_star, p.charge_q
    de = r*r-2*r+a*a+q*q
    fp = -2/(r*r)+2*q*q/(r*r*r)
    dgtt, dgtp, dgpp = fp, -a*fp, 2*r+a*a*fp
    dgrr = (2*r*de-r*r*(2*r-2))/(de*de)
    h = 2e-6
    dur = (flow.four_velocity(r-2*h)["ur"]-8*flow.four_velocity(r-h)["ur"]
           +8*flow.four_velocity(r+h)["ur"]-flow.four_velocity(r+2*h)["ur"])/(12*h)
    cov_accel = v["ur"]*dur+.5*de/(r*r)*(dgrr*v["ur"]**2
        -dgtt*v["ut"]**2-2*dgtp*v["ut"]*v["uphi"]-dgpp*v["uphi"]**2)
    err = float(np.max(abs(cov_accel)))
    check("independent_radial_Christoffel_geodesic_equation", err, err < 1e-9)

    # Independent full -p_mu*u^mu contraction and circular-limit redshift.
    r = np.linspace(flow.isco_rg+.001, 28, 1001)
    lam = np.linspace(-2,2,r.size)
    eta = .5
    de = r*r-2*r+a*a+q*q
    R = (r*r+a*a-a*lam)**2-de*((lam-a)**2+eta)
    kr = np.sqrt(R)/de
    for radial_sign in [-1,1]:
        v = flow.four_velocity(r)
        covp = np.stack((-np.ones_like(r),radial_sign*kr,lam),axis=-1)
        contrau = np.stack((v["ut"],v["ur"],v["uphi"]),axis=-1)
        expected = -1/np.einsum("...i,...i->...",covp,contrau)
        got = flow.redshift(r,lam,radial_sign*kr)
        err = float(np.max(abs(got-expected)))
        check(f"independent_photon_contraction_radial_sign_{radial_sign}",err,err < 1e-15)
    small = NeutralAccretionFlow(alpha_viscosity=1e-7, radial_samples=8193)
    circ = orbital_quantities(r,.9,.4)
    g_circ = 1/(circ["ut"]*(1-lam*circ["omega"]))
    err = float(np.max(abs(small.redshift(r,lam,kr)/g_circ-1)))
    check("zero_drift_limit_recovers_existing_circular_redshift",err,err < 2e-9)
    # Outgoing escaping equatorial photons redshift away at the horizon;
    # this is calculated by their four-momentum, not by a screen-space fade.
    gap = np.array([.01,.001,.0001,.00001])
    r = flow.horizon_rg+gap
    de = (r-flow.horizon_rg)*(r-(2-flow.horizon_rg))
    kr = np.sqrt((r*r+a*a)**2-de*a*a)/de
    shifts = flow.redshift(r,0,kr)
    check("outgoing_near_horizon_photons_redshift_to_zero", shifts.tolist(),
          np.all(np.diff(shifts) < 0) and shifts[-1] < 1e-5)
    check("near_horizon_outgoing_redshift_linear_gap_limit",float(np.ptp(shifts[-2:]/gap[-2:])),
          np.ptp(shifts[-2:]/gap[-2:]) < 4e-5)
    # A constant-time surface at the mathematical horizon has no emitter.
    rejected = 0
    for radius in [flow.horizon_rg,flow.horizon_rg-.001,float("nan")]:
        try:
            flow.four_velocity(radius)
        except ValueError:
            rejected += 1
    check("no_emitter_on_inside_or_invalid_horizon",rejected,rejected == 3)

    # Independent Gauss--Legendre radial quadratures in log(r-r+) resolve
    # the coordinate-time divergence and the slow but finite ISCO crossing.
    nodes, weights = np.polynomial.legendre.leggauss(256)
    errors = {k:[] for k in ["age_tg","azimuth_rad","proper_age_tg"]}
    for lower in [flow.inner_radius_rg,.5*(flow.horizon_rg+flow.isco_rg),flow.isco_rg-.01]:
        z0,z1 = np.log(lower-flow.horizon_rg),np.log(flow.isco_rg-flow.horizon_rg)
        z = .5*(z1-z0)*nodes+.5*(z0+z1)
        r = flow.horizon_rg+np.exp(z)
        v = flow.four_velocity(r)
        for name,numerator in [("age_tg",v["ut"]),("azimuth_rad",v["uphi"]),("proper_age_tg",np.ones_like(r))]:
            reference = .5*(z1-z0)*np.dot(weights,numerator/(-v["ur"])*np.exp(z))
            errors[name].append(abs(float(flow.streamline(lower)[name])-reference))
    for name,errs in errors.items():
        check("independent_gauss_legendre_streamline_"+name,max(errs),max(errs) < .0003)
    r = flow.horizon_rg+np.geomspace(p.horizon_offset_rg,28-flow.horizon_rg,2001)
    refined = NeutralAccretionFlow(radial_samples=131073)
    for name in ["age_tg","azimuth_rad","proper_age_tg"]:
        coarse,fine = flow.streamline(r)[name],refined.streamline(r)[name]
        error = float(np.max(abs(coarse-fine)))
        check("doubled_grid_streamline_"+name,error,error < .0003)
    # A finite material displacement follows dr/dt and dphi/dt calculated
    # from the metric, independently of the invariant-table labels.
    r = np.linspace(flow.horizon_rg+.05,28,2000)
    dr = 2e-5
    plus,minus = flow.streamline(r+dr),flow.streamline(r-dr)
    dT = (plus["age_tg"]-minus["age_tg"])/(2*dr)
    dA = (plus["azimuth_rad"]-minus["azimuth_rad"])/(2*dr)
    v = flow.four_velocity(r)
    eT = float(np.max(abs(dT*v["dr_dt"]-1)))
    eA = float(np.max(abs(dA*v["dr_dt"]-v["omega"])))
    check("streamline_time_derivative_matches_physical_radial_velocity",eT,eT < 5e-4)
    check("streamline_azimuth_derivative_matches_physical_angular_velocity",eA,eA < 2e-4)
    r = np.concatenate((np.linspace(flow.inner_radius_rg,flow.isco_rg,501),np.linspace(flow.isco_rg,28,501)))
    ages = flow.streamline(r)["age_tg"]
    inv = float(np.max(abs(flow.radius_at_age(ages)-r)))
    check("streamline_age_inverse_radius",inv,inv < 2e-13)
    angle = float(np.max(abs(flow.azimuth_at_age(ages)-flow.streamline(r)["azimuth_rad"])))
    check("streamline_age_inverse_azimuth",angle,angle < 2e-11)
    check("ISCO_flow_age_and_azimuth_origin",[float(v) for v in flow.streamline(flow.isco_rg).values()],
          all(float(v) == 0 for v in flow.streamline(flow.isco_rg).values()))
    check("source_discloses_prescribed_outer_flow_and_unsolved_emissivity",None,
          flow.metadata()["outer_drift_is_prescribed"] and "plunge emissivity" in flow.metadata()["not_solved"])

    return {"passed": all(v["passed"] for v in checks),"check_count":len(checks),"checks":checks,
        "flow":flow.metadata(),
        "scope":"numerical consistency of idealized neutral test-matter kinematics; not a viscous, MHD or thermodynamic accretion solution",
        "source_sha256":{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                         for name in ["kn_flow.py","kn_orbits.py","verify_kn_flow.py"]}}


if __name__ == "__main__":
    report = diagnostics()
    (ROOT/"kn-flow-diagnostics.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"passed":report["passed"],"checks":report["check_count"],
                      "report":"kn-flow-diagnostics.json"},indent=2))
