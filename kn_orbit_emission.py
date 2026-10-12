"""Relative bolometric thin-disk emission in a Kerr--Newman geometry.

This is an explicitly idealized, NEUTRAL test-matter disk. The hole has charge
q=Q/(M sqrt(4 pi epsilon_0 G)), while the orbiting emitters have zero specific
charge. Charge changes the metric, circular orbits, ISCO and radial flux.
No charged-particle Lorentz forces, plasma currents or charge evolution are
calculated. A large q=.4 is a theoretical example, not an inferred realistic
astrophysical charge. The hole and disk align for the shipped a*=+.9.

Page & Thorne's general stationary-axisymmetric thin-disk conservation law
(1974, Eqs. 11b,12) is reevaluated with neutral Kerr--Newman E,L,Omega and the
near-equatorial determinant factor sqrt(-g)=r:

 F(r) ~ -Omega'(r)/[r (E-Omega L)^2] integral_ISCO^r (E-Omega L)L' dr
 Omega(r)=sqrt(r-q^2)/(r^2+a*sqrt(r-q^2)), in G=M=c=1

The stress vanishes at the computed neutral ISCO; no emission is assigned
inside it. Only the relative flux shape is used, normalized to maximum one.
Mass, accretion rate, temperature and spectrum remain unspecified. Both faces
are opaque and emit isotropically. No returning-radiation heating, inflow,
finite thickness or electromagnetic disk dynamics is included.

Finite-lived patches reuse the earlier source prescription, now advected at
the exact Kerr--Newman frequency. Their azimuthal mean equals the radial flux.
Births repeat every 800 GM/c^3 to construct a loop; orbital frequencies are
never rounded. Traced per-ray delays enter as t_em=t_obs-delay.

Sources:
* Page & Thorne1974, doi:10.1086/152990, Eqs.11b,12:
  https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P
* Pugliese, Quevedo & Ruffini2013, neutral Kerr--Newman orbits:
  https://arxiv.org/abs/1303.6250
* Ivanov & Prodanov2005, neutral E,L Eqs.19--20:
  https://arxiv.org/abs/gr-qc/0504025
* Schnittman & Bertschinger2004, orbiting/sheared hot spots:
  https://arxiv.org/abs/astro-ph/0309458
* Cardenas-Avendano, Lupsasca & Zhu2023, retarded source methodology
  (not their stochastic solver): https://arxiv.org/abs/2211.07469
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import numpy as np

from kn_orbits import isco_radius, orbital_quantities
from nasa_emission import ShearingDiskEmission as SchwarzschildShearingDiskEmission
from nasa_emission import PreparedEmission as SchwarzschildPreparedEmission


@dataclass(frozen=True)
class StructureParameters:
    spin_a_star: float = 0.9
    charge_q: float = 0.4
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

    @property
    def inner_radius_rg(self):
        return isco_radius(self.spin_a_star, self.charge_q)

    @property
    def isco_rg(self):
        return self.inner_radius_rg

    def __post_init__(self):
        rin = self.inner_radius_rg
        numeric = [self.outer_radius_rg, self.loop_period_tg, self.contrast,
                   self.radial_log_half_width, self.azimuthal_kappa,
                   self.lifetime_orbits, self.maximum_lifetime_loop_fraction,
                   self.mean_overlapping_births]
        if not all(np.isfinite(numeric)):
            raise ValueError("structure parameters must be finite")
        if self.outer_radius_rg <= rin or self.loop_period_tg <= 0 or self.contrast < 0:
            raise ValueError("outer radius must exceed ISCO; positive period and nonnegative contrast required")
        if (not isinstance(self.radial_bands, int) or self.radial_bands < 2 or
                self.radial_log_half_width <= 0 or not 0 < self.azimuthal_kappa <= 50 or
                self.lifetime_orbits <= 0 or not 0 < self.maximum_lifetime_loop_fraction < 1 or
                self.mean_overlapping_births <= 0):
            raise ValueError("invalid structure scales")


class DimensionlessPageThorne:
    """Neutral Kerr--Newman flux shape only; no dimensional disk parameters are constructed."""

    def __init__(self, spin_a_star=0.9, charge_q=0.4, outer_radius_rg=28.0, samples=65537):
        self.spin_a_star = float(spin_a_star)
        self.charge_q = float(charge_q)
        self.rin = isco_radius(spin_a_star, charge_q)
        self.rout = float(outer_radius_rg)
        if not np.isfinite(self.rout) or self.rout <= self.rin or samples < 1025:
            raise ValueError("invalid integration interval or sample count")
        self.samples = int(samples)
        # Resolve the exact quadratic stress-free edge behavior with a grid
        # logarithmic in r-rISCO, rather than wasting samples at large radii.
        delta = np.geomspace(1e-7*self.rin, self.rout-self.rin, self.samples-1)
        rr = np.concatenate(([self.rin], self.rin+delta))
        q = orbital_quantities(rr, self.spin_a_star, self.charge_q)
        integrand = q["energy_minus_omega_l"]*q["dangular_momentum_dr"]
        integrand[0] = 0.0  # Analytically L'(ISCO)=0; remove roundoff only here.
        if np.any(integrand[1:] <= 0):
            raise ValueError("nonpositive angular-momentum derivative outside ISCO")
        integral = np.concatenate(([0.0], np.cumsum(
            0.5*(integrand[1:]+integrand[:-1])*np.diff(rr))))
        shape = -q["domega_dr"]/(rr*q["energy_minus_omega_l"]**2)*integral
        if not np.all(np.isfinite(shape)) or np.any(shape[1:] <= 0):
            raise ValueError("invalid dimensionless flux")
        self._log_delta = np.log(delta)
        self._log_shape = np.log(shape[1:])

    def flux_shape(self, radius_rg):
        r = np.asarray(radius_rg, dtype=np.float64)
        if np.any(~np.isfinite(r)) or np.any(r > self.rout):
            raise ValueError("flux radii must be finite and within integration grid")
        delta = r-self.rin
        log_delta = np.log(np.maximum(delta, np.finfo(float).tiny))
        log_shape = np.interp(log_delta, self._log_delta, self._log_shape)
        # Below the first tiny interval, use the analytic quadratic edge
        # scaling. Negative/zero delta is masked after safe evaluation.
        log_shape = np.where(log_delta < self._log_delta[0],
                             self._log_shape[0]+2*(log_delta-self._log_delta[0]), log_shape)
        return np.where(delta > 0, np.exp(log_shape), 0.0)


class ShearingDiskEmission(SchwarzschildShearingDiskEmission):
    """Drop-in relative bolometric source for a neutral Kerr--Newman disk.

    evaluate(r,phi,t_em), prepare(r,phi,delay).evaluate(t_obs), radial_profile(r),
    and metadata() retain the Schwarzschild renderer's source interface.
    """

    def __init__(self, outer_radius_rg=28.0, loop_period_tg=800.0,
                 contrast=1.6, seed=13326, spin_a_star=0.9, charge_q=0.4, **structure_options):
        self.parameters = StructureParameters(
            spin_a_star=spin_a_star, charge_q=charge_q, outer_radius_rg=outer_radius_rg,
            loop_period_tg=loop_period_tg, contrast=contrast, seed=seed,
            **structure_options)
        p = self.parameters
        self.rin = p.inner_radius_rg
        self.flux_model = DimensionlessPageThorne(p.spin_a_star, p.charge_q, p.outer_radius_rg)
        rng = np.random.default_rng(p.seed)
        spacing = math.log(p.outer_radius_rg/self.rin)/(p.radial_bands-1)
        centers = np.linspace(math.log(self.rin), math.log(p.outer_radius_rg), p.radial_bands)
        centers[1:-1] += rng.uniform(-0.22, 0.22, p.radial_bands-2)*spacing
        self.log_centers = centers
        self.radii = np.exp(centers)
        self.lifetimes = np.minimum(
            p.lifetime_orbits*2*math.pi/self.angular_velocity(self.radii),
            p.maximum_lifetime_loop_fraction*p.loop_period_tg)
        self.births, self.phases = [], []
        for lifetime in self.lifetimes:
            count = max(2, int(math.ceil(p.mean_overlapping_births*p.loop_period_tg/lifetime)))
            births = ((np.arange(count)+rng.uniform(0.15, 0.85, count))/count*p.loop_period_tg)
            self.births.append(births)
            self.phases.append(rng.uniform(0, 2*math.pi, count))
        self.angular_normalization = float(np.i0(p.azimuthal_kappa))
        a, b = self.rin, p.outer_radius_rg
        q = (math.sqrt(5)-1)/2
        for _ in range(90):
            x, y = b-q*(b-a), a+q*(b-a)
            if float(self.flux_model.flux_shape(x)) < float(self.flux_model.flux_shape(y)):
                a = x
            else:
                b = y
        self.peak_radius_rg = (a+b)/2
        self.radial_peak = float(self.flux_model.flux_shape(self.peak_radius_rg))

    def angular_velocity(self, radius_rg):
        r = np.asarray(radius_rg, dtype=np.float64)
        s = np.sqrt(r-self.parameters.charge_q**2)
        return s/(r*r+self.parameters.spin_a_star*s)

    def radial_profile(self, radius_rg):
        r = np.asarray(radius_rg, dtype=np.float64)
        good = np.isfinite(r) & (r > self.rin) & (r <= self.parameters.outer_radius_rg)
        safe = np.where(good, r, self.rin)
        return np.where(good, self.flux_model.flux_shape(safe)/self.radial_peak, 0.0)

    def prepare(self, radius_rg, phi_rad, delay_tg=0.0):
        return PreparedEmission(self, radius_rg, phi_rad, delay_tg)

    def metadata(self):
        data = super().metadata()
        data.update({
            **asdict(self.parameters),
            "inner_radius_rg": self.rin,
            "spin_alignment": ("aligned / prograde" if self.parameters.spin_a_star > 0 else
                               "antialigned / retrograde" if self.parameters.spin_a_star < 0 else
                               "nonrotating hole; positive disk angular velocity"),
            "spin_definition": "a*=Jc/(G M^2); chosen model parameter, not a measurement",
            "orbital_direction": "positive Boyer-Lindquist azimuth",
            "charge_definition": "q=Q/(M*sqrt(4*pi*epsilon_0*G)); chosen theoretical parameter",
            "emitter_specific_charge": 0.0,
            "matter_assumption": "neutral test-matter disk in a fixed charged geometry; no Lorentz force or plasma dynamics",
            "radial_profile": "general Page-Thorne conservation integral with neutral Kerr-Newman circular orbits; zero ISCO stress; peak=1",
            "radial_quadrature": "cumulative trapezoid with analytic circular-orbit derivatives",
            "radial_quadrature_samples": self.flux_model.samples,
            "radial_lookup": "log flux versus log(r-rISCO); quadratic unresolved edge limit",
            "angular_velocity_tg_inverse": "Omega(r)=sqrt(r-q^2)/(r^2+a*sqrt(r-q^2)), never quantized",
            "dimensional_normalization": "none; mass, accretion rate, temperature and spectrum unspecified",
            "surface_emission": "two identical opaque, isotropically emitting faces",
            "flux_equation": "F proportional to -Omega'/(r*(E-Omega*L)^2) integral_ISCO^r (E-Omega*L)*L' dr",
        })
        data["not_modeled"] = ["disk magnetohydrodynamics or MRI" if v == "magnetic fields or MRI" else v for v in data["not_modeled"]]
        data["not_modeled"] += ["charged-particle Lorentz forces", "plasma currents",
                                "charge neutralization or evolution", "disk backreaction"]
        data["sources"] += [
            {"arxiv": "1303.6250", "role": "neutral Kerr-Newman circular orbits and metric"},
            {"arxiv": "gr-qc/0504025", "role": "neutral specific energy and angular momentum, Eqs.19-20"},
            {"arxiv": "1706.05466", "role": "independent published neutral ISCO equation benchmark, not the root solver"},
        ]
        return data


class PreparedEmission(SchwarzschildPreparedEmission):
    """Reuse identical patch weights and retarded times, with neutral Kerr--Newman advection."""

    def __init__(self, source, radius_rg, phi_rad, delay_tg):
        super().__init__(source, radius_rg, phi_rad, delay_tg)
        r = np.broadcast_arrays(np.asarray(radius_rg, dtype=np.float64),
                                np.asarray(phi_rad, dtype=np.float64),
                                np.asarray(delay_tg, dtype=np.float64))[0].ravel()
        self.omega = source.angular_velocity(r[self.valid_indices])
