"""Dimensionless, structured emission for an ORIGINAL Schwarzschild disk.

This is not NASA software, an inferred NASA texture, or GRMHD. The NASA-like
feature is a thin disk containing evolving, sheared bright and dark lanes.
The azimuthal mean is the zero-torque Page--Thorne flux (from emission.py).
There is no assigned mass, accretion rate, temperature, spectrum, or RGB color.
The caller ray traces the surface and applies the bolometric g**4 transfer.

The optional structure is a deliberately small kinematic source model:

  J(r,phi,t) = J_PT(r) [1+c sum_j w_j(r,t) A_j(r,phi,t)] /
                       [1+c sum_j w_j(r,t)]
  A_j = exp[kappa cos(phi-phi_j-Omega(r) age_j)] / I0(kappa)
  Omega(r) = r**(-3/2),  age_j = (t-birth_j) mod P.

Each nonnegative weight has compact radial support and a sin**4 lifetime
envelope. A_j has azimuthal mean one, so <J>_phi=J_PT at every r and t.
Different radii use their exact Schwarzschild circular-orbit frequencies;
there is no rounding to loop harmonics. The source births repeat after P,
and die smoothly before reappearing. Thus the SOURCE is periodically forced
for the video; gas is not claimed to complete a closed trajectory in P.

Sources informing the choices (not copied implementations):
* Page & Thorne (1974), doi:10.1086/152990: mean radiative flux.
* Schnittman & Bertschinger (2004), doi:10.1086/383180,
  https://arxiv.org/abs/astro-ph/0309458: ray-traced orbiting/sheared hot spots.
* Cardenas-Avendano, Lupsasca & Zhu (2023), doi:10.1103/PhysRevD.107.043030,
  https://arxiv.org/abs/2211.07469: separate fluctuating equatorial emissivity,
  retarded emission times, orbital-time correlation, and Eq.181 for Omega.
  This module DOES NOT implement their Matern/inoisy stochastic field.

Widths, contrast, lifetimes, seed, radial cutoff, and loop period are exposed
assumptions, not quantities derived from those papers or fitted to data.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import time

import numpy as np

from emission import nt_flux_shape


@dataclass(frozen=True)
class StructureParameters:
    outer_radius_rg: float = 28.0
    loop_period_tg: float = 800.0
    contrast: float = 1.6
    seed: int = 13326
    radial_bands: int = 64
    radial_log_half_width: float = 0.04
    azimuthal_kappa: float = 2.5
    lifetime_orbits: float = 1.5
    maximum_lifetime_loop_fraction: float = 0.8
    mean_overlapping_births: float = 2.0

    def __post_init__(self):
        numeric = [self.outer_radius_rg, self.loop_period_tg, self.contrast,
                   self.radial_log_half_width, self.azimuthal_kappa,
                   self.lifetime_orbits, self.maximum_lifetime_loop_fraction,
                   self.mean_overlapping_births]
        if not all(np.isfinite(numeric)):
            raise ValueError("structure parameters must be finite")
        if self.outer_radius_rg <= 6 or self.loop_period_tg <= 0 or self.contrast < 0:
            raise ValueError("outer radius >6, positive period, nonnegative contrast required")
        if (self.radial_bands < 2 or self.radial_log_half_width <= 0 or
                not 0 < self.azimuthal_kappa <= 50 or self.lifetime_orbits <= 0 or
                not 0 < self.maximum_lifetime_loop_fraction < 1 or
                self.mean_overlapping_births <= 0):
            raise ValueError("invalid structure scales")


class ShearingDiskEmission:
    """Relative bolometric surface intensity; smooth radial maximum is one.

    ``evaluate(r, phi, emission_time_tg)`` accepts broadcast-compatible arrays.
    Radius is in GM/c**2, phi in radians, time in GM/c**3. The azimuth convention
    must match the ray tracer's positive orbital direction. Invalid/missed rays
    (r<=6 or r>outer radius) emit zero. For a static ray map, ``prepare`` avoids
    recalculating radial neighborhoods on every movie frame.
    """

    def __init__(self, outer_radius_rg=28.0, loop_period_tg=800.0,
                 contrast=1.6, seed=13326, **structure_options):
        self.parameters = StructureParameters(
            outer_radius_rg=outer_radius_rg, loop_period_tg=loop_period_tg,
            contrast=contrast, seed=seed, **structure_options)
        p = self.parameters
        rng = np.random.default_rng(p.seed)
        spacing = math.log(p.outer_radius_rg / 6.0) / (p.radial_bands - 1)
        centers = np.linspace(math.log(6.0), math.log(p.outer_radius_rg), p.radial_bands)
        # Interior jitter avoids an exactly regular stack of rings. Keep the
        # endpoint bands so the entire emitting interval remains covered.
        centers[1:-1] += rng.uniform(-0.22, 0.22, p.radial_bands - 2) * spacing
        self.log_centers = centers
        self.radii = np.exp(centers)
        self.lifetimes = np.minimum(
            p.lifetime_orbits * 2 * math.pi * self.radii**1.5,
            p.maximum_lifetime_loop_fraction * p.loop_period_tg)
        self.births = []
        self.phases = []
        for lifetime in self.lifetimes:
            count = max(2, int(math.ceil(p.mean_overlapping_births * p.loop_period_tg / lifetime)))
            # Stratified birth times avoid accidental long intervals with no
            # structure. They are a reproducible illustrative forcing pattern.
            births = ((np.arange(count) + rng.uniform(0.15, 0.85, count)) /
                      count * p.loop_period_tg)
            self.births.append(births)
            self.phases.append(rng.uniform(0, 2 * math.pi, count))
        self.angular_normalization = float(np.i0(p.azimuthal_kappa))
        # A global unit conversion only, independent of frame and modulation.
        # One-dimensional golden-section maximization keeps this normalization
        # independent of ray resolution or camera framing.
        a, b = 6.0, p.outer_radius_rg
        q = (math.sqrt(5) - 1) / 2
        for _ in range(90):
            x, y = b - q * (b - a), a + q * (b - a)
            if float(nt_flux_shape(x)) < float(nt_flux_shape(y)):
                a = x
            else:
                b = y
        self.peak_radius_rg = (a + b) / 2
        self.radial_peak = float(nt_flux_shape(self.peak_radius_rg))

    def radial_profile(self, radius_rg):
        r = np.asarray(radius_rg, dtype=np.float64)
        good = np.isfinite(r) & (r > 6) & (r <= self.parameters.outer_radius_rg)
        safe = np.where(good, r, 6.0)
        return np.where(good, nt_flux_shape(safe) / self.radial_peak, 0.0)

    def prepare(self, radius_rg, phi_rad, delay_tg=0.0):
        """Cache static ray coordinates; emission time is observer time-delay.

        Any COMMON constant may have been subtracted from every ray delay.
        The sign convention is explicit: longer delays sample earlier source
        times. No photon flight time is calculated or approximated here.
        """
        return PreparedEmission(self, radius_rg, phi_rad, delay_tg)

    def evaluate(self, radius_rg, phi_rad, emission_time_tg, structured=True):
        r, phi, t = np.broadcast_arrays(
            np.asarray(radius_rg, dtype=np.float64),
            np.asarray(phi_rad, dtype=np.float64),
            np.asarray(emission_time_tg, dtype=np.float64))
        # Reuse exactly the prepared evaluator so scalar and retarded-time
        # paths cannot silently implement different emission models.
        return self.prepare(r, phi, -t).evaluate(0.0, structured=structured)

    def metadata(self):
        return {
            **asdict(self.parameters),
            "model": "original kinematic source, not NASA software or GRMHD",
            "intensity": "relative bolometric comoving surface intensity",
            "radial_profile": "zero-torque Schwarzschild Page-Thorne flux; peak=1",
            "radial_peak_rg": self.peak_radius_rg,
            "unscaled_radial_peak": self.radial_peak,
            "angular_velocity_tg_inverse": "Omega(r)=r^(-3/2), never quantized",
            "radial_weight": "exp[-2 x^2/(1-x^2)], |x|<1; x=ln(r/r_j)/half_width",
            "radial_e_folding_fraction_approx": self.parameters.radial_log_half_width / math.sqrt(3),
            "angular_weight": "exp[kappa*cos(phi-phi_j-Omega(r)*age)]/I0(kappa)",
            "lifetime_envelope": "sin(pi*age/lifetime)^4 on 0<age<lifetime; zero outside",
            "lifetime_tg_range": [float(np.min(self.lifetimes)), float(np.max(self.lifetimes))],
            "birth_count": int(sum(len(x) for x in self.births)),
            "normalization": "exact azimuthal mean equals Page-Thorne at each emission time and radius",
            "loop_assumption": "periodically repeated finite-lived source births; not repeating gas orbits",
            "retarded_time": "evaluate at t_em=t_obs-delay for each ray",
            "not_modeled": ["magnetic fields or MRI", "fluid dynamics", "disk thickness",
                            "radial inflow", "plunging-region emission", "thermal spectra",
                            "polarization", "corona", "returning radiation"],
            "sources": [
                {"doi": "10.1086/152990", "role": "mean thin-disk flux"},
                {"doi": "10.1086/383180", "role": "orbiting/sheared inhomogeneities as a kinematic model"},
                {"doi": "10.1103/PhysRevD.107.043030", "role": "variable equatorial source, retarded time, orbital scales; not its stochastic solver"},
            ],
        }


class PreparedEmission:
    """Sparse radial neighborhoods for a fixed array of traced disk hits."""

    def __init__(self, source, radius_rg, phi_rad, delay_tg):
        self.source = source
        r, phi, delay = np.broadcast_arrays(
            np.asarray(radius_rg, dtype=np.float64),
            np.asarray(phi_rad, dtype=np.float64),
            np.asarray(delay_tg, dtype=np.float64))
        self.shape = r.shape
        r, phi, delay = r.ravel(), phi.ravel(), delay.ravel()
        self.base = source.radial_profile(r)
        good = self.base > 0
        if np.any(~np.isfinite(phi[good])) or np.any(~np.isfinite(delay[good])):
            raise ValueError("emitting rays require finite azimuth and emission delay")
        self.valid_indices = np.flatnonzero(good)
        rr = r[good]
        self.phi = phi[good].copy()
        self.delay = delay[good].copy()
        self.omega = rr**(-1.5)
        logr = np.log(rr)
        self.groups = []
        for j, center in enumerate(source.log_centers):
            x = (logr - center) / source.parameters.radial_log_half_width
            idx = np.flatnonzero(np.abs(x) < 1)
            xx = x[idx]
            weight = np.exp(-2 * xx * xx / (1 - xx * xx))
            self.groups.append((j, idx, weight))

    def evaluate(self, observer_time_tg, structured=True):
        if not np.ndim(observer_time_tg) == 0 or not np.isfinite(observer_time_tg):
            raise ValueError("prepared evaluation requires a finite scalar observer time")
        source, p = self.source, self.source.parameters
        out = self.base.copy()
        if not structured or p.contrast == 0 or not len(self.valid_indices):
            return out.reshape(self.shape)
        numerator = np.ones(len(self.valid_indices), dtype=np.float64)
        denominator = np.ones_like(numerator)
        emission_time = float(observer_time_tg) - self.delay
        for j, idx, radial_weight in self.groups:
            if not len(idx):
                continue
            subtime, subphi, omega = emission_time[idx], self.phi[idx], self.omega[idx]
            life = source.lifetimes[j]
            radial_weight = p.contrast * radial_weight
            for birth, phase in zip(source.births[j], source.phases[j]):
                age = np.remainder(subtime - birth, p.loop_period_tg)
                active = age < life
                if not np.any(active):
                    continue
                ia = idx[active]
                age = age[active]
                envelope = np.sin(math.pi * age / life)**4
                weight = radial_weight[active] * envelope
                angle = subphi[active] - phase - omega[active] * age
                angular = np.exp(p.azimuthal_kappa * np.cos(angle)) / source.angular_normalization
                numerator[ia] += weight * angular
                denominator[ia] += weight
        out[self.valid_indices] *= numerator / denominator
        return out.reshape(self.shape)


def diagnostics():
    """Numerical source checks, not validation against NASA or GRMHD."""
    source = ShearingDiskEmission()
    rng = np.random.default_rng(271828)
    radii = np.geomspace(6.001, 28, 80)[:, None]
    phi = np.arange(4096)[None, :] * (2 * math.pi / 4096)
    mean_errors = []
    for t in [0, 79.31, 311.71, 799.99]:
        sampled = source.evaluate(radii, phi, t)
        mean_errors.append(float(np.max(np.abs(sampled.mean(axis=1) /
                                              source.radial_profile(radii[:, 0]) - 1))))
    r = rng.uniform(6.01, 28, 8000)
    ph = rng.uniform(-3 * math.pi, 3 * math.pi, r.size)
    t = rng.uniform(-1600, 1600, r.size)
    a = source.evaluate(r, ph, t)
    b = source.evaluate(r, ph, t + source.parameters.loop_period_tg)
    delay = rng.uniform(0, 180, r.size)
    prep = source.prepare(r, ph, delay)
    reference = source.evaluate(r, ph, 137.5 - delay)
    cached = prep.evaluate(137.5)
    # Independent scalar implementation of the stated equation. Comparing
    # evaluate() with prepare() alone would not catch a shared implementation
    # error because they deliberately use the same production evaluator.
    def scalar_reference(radius, azimuth, emission_time):
        numerator, denominator = 1.0, 1.0
        p = source.parameters
        for j, center in enumerate(source.log_centers):
            x = (math.log(radius) - center) / p.radial_log_half_width
            if abs(x) >= 1:
                continue
            radial = math.exp(-2 * x * x / (1 - x * x))
            for birth, phase in zip(source.births[j], source.phases[j]):
                age = (emission_time - float(birth)) % p.loop_period_tg
                if age >= source.lifetimes[j]:
                    continue
                weight = p.contrast * radial * math.sin(math.pi * age / source.lifetimes[j])**4
                angle = azimuth - float(phase) - age / radius**1.5
                numerator += weight * math.exp(p.azimuthal_kappa * math.cos(angle)) / source.angular_normalization
                denominator += weight
        return float(source.radial_profile(radius)) * numerator / denominator
    scalar_times = [0.0, -11.3, 799.999, 823.1]
    scalar_errors = []
    for observer_time in scalar_times:
        reference_samples = np.array([scalar_reference(float(r[i]), float(ph[i]),
                                                      observer_time - float(delay[i]))
                                      for i in range(50)])
        scalar_errors.append(float(np.max(np.abs(prep.evaluate(observer_time)[:50] - reference_samples))))
    eps = 1e-3
    seam_left = (source.evaluate(r, ph, 0) - source.evaluate(r, ph, -eps)) / eps
    seam_right = (source.evaluate(r, ph, eps) - source.evaluate(r, ph, 0)) / eps
    seam_difference = float(np.max(np.abs(seam_left - seam_right)))
    half_eps = eps / 2
    seam_left_half = (source.evaluate(r, ph, 0) - source.evaluate(r, ph, -half_eps)) / half_eps
    seam_right_half = (source.evaluate(r, ph, half_eps) - source.evaluate(r, ph, 0)) / half_eps
    seam_difference_half = float(np.max(np.abs(seam_left_half - seam_right_half)))
    # A genuine derivative jump would not tend to zero when halving dt.
    # These finite-lived sin^4 patches are C^3 across each birth/death.
    orbital_period_at_10_rg = 2 * math.pi * 10**1.5
    dphi_over_loop = source.parameters.loop_period_tg / 10**1.5
    no_structure = prep.evaluate(137.5, structured=False)
    changed = prep.evaluate(137.5 + source.parameters.loop_period_tg / 240)
    # A realistic large static map includes radius/phi/delay varying per ray.
    br = rng.uniform(6.01, 28, 300000)
    bp = rng.uniform(0, 2 * math.pi, br.size)
    bd = rng.uniform(0, 180, br.size)
    started = time.perf_counter()
    prepared = source.prepare(br, bp, bd)
    prepare_seconds = time.perf_counter() - started
    started = time.perf_counter()
    evaluated = prepared.evaluate(123.45)
    frame_seconds = time.perf_counter() - started
    result = {
        "metadata": source.metadata(),
        "azimuthal_mean_max_relative_error": max(mean_errors),
        "periodicity_max_absolute_error": float(np.max(np.abs(a - b))),
        "prepared_vs_direct_max_absolute_error": float(np.max(np.abs(reference - cached))),
        "independent_scalar_equation_max_absolute_error": max(scalar_errors),
        "seam_left_right_derivative_max_difference_dt_1e_minus3": seam_difference,
        "seam_left_right_derivative_max_difference_dt_5e_minus4": seam_difference_half,
        "seam_derivative_difference_halving_ratio": seam_difference_half / seam_difference,
        "smooth_mode_max_absolute_error": float(np.max(np.abs(no_structure - source.radial_profile(r)))),
        "rms_source_change_over_one_video_frame": float(np.sqrt(np.mean((changed - cached)**2))),
        "orbital_kinematics_example": {
            "radius_rg": 10.0,
            "period_tg": orbital_period_at_10_rg,
            "loop_period_tg": source.parameters.loop_period_tg,
            "orbital_turns_per_loop": dphi_over_loop / (2 * math.pi),
            "interpretation": "Noninteger orbital turns: the source births repeat; the orbit frequency is not forced to loop."
        },
        "positive_emission_min": float(np.min(a)),
        "structured_to_mean_factor_range": [float(np.min(a / source.radial_profile(r))), float(np.max(a / source.radial_profile(r)))],
        "masked_radius_output": source.evaluate([0, 2, 6, 28.1, np.nan], 0, 0).tolist(),
        "benchmark": {"emitting_rays": int(br.size), "prepare_seconds": prepare_seconds,
                      "frame_seconds": frame_seconds, "all_finite": bool(np.all(np.isfinite(evaluated)))},
        "scope": "source-model invariants only; not a NASA reconstruction or GRMHD validation",
    }
    assert result["azimuthal_mean_max_relative_error"] < 1e-11
    assert result["periodicity_max_absolute_error"] < 1e-11
    assert result["prepared_vs_direct_max_absolute_error"] < 1e-12
    assert result["independent_scalar_equation_max_absolute_error"] < 1e-12
    assert result["smooth_mode_max_absolute_error"] == 0
    assert 0.45 < result["seam_derivative_difference_halving_ratio"] < 0.55
    assert result["rms_source_change_over_one_video_frame"] > 0
    assert result["positive_emission_min"] >= 0
    assert result["masked_radius_output"] == [0.0] * 5
    assert result["benchmark"]["all_finite"]
    return result


if __name__ == "__main__":
    report = diagnostics()
    Path(__file__).with_name("nasa-emission-diagnostics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
