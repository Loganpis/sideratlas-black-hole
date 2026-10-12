"""Dimensionless Page--Thorne emission and finite-lived Kerr disk structure.

The production Kerr orbit animation uses this relative bolometric source, with
no assigned mass, accretion rate, temperature, spectrum, or visible RGB color.
The ray tracer supplies Kerr disk coordinates, g, and photon travel delays; the
renderer applies I_obs = g**4 I_em before its illustrative intensity palette.

Positive disk angular velocity is the convention. spin_a_star > 0 therefore
means an equatorial prograde disk aligned with the hole's angular momentum;
spin_a_star < 0 is the antialigned branch. The shipped choice is +0.9, a model
parameter, not a measurement of NASA's visualization or a real black hole.

Page & Thorne (1974), Eqs. (11b), (12), and (15), give the flux shape in
G=M=c=1, apart from a radius-independent normalization:

  F(r) ~ -Omega'(r) / [r (E-Omega L)**2]
         * integral_rISCO^r (E-Omega L) L'(x) dx
  Omega(r) = 1 / (r**(3/2) + a*)

We normalize its radial maximum to one. This removes an undetermined overall
scale without imposing a physical mass or accretion rate. The inner stress is
zero and there is no emission inside the ISCO. Both faces emit identically and
isotropically; returning-radiation heating, plunging emission, radial inflow,
finite thickness, and magnetic/plasma evolution are omitted.

The existing finite-lived patch evaluator is reused, with the exact Kerr
Omega substituted in both lifetimes and phase advection. Its azimuthal mean is
one at every radius and emission time. Patch births repeat every 800 GM/c**3;
orbital frequencies are never rounded to force a loop. Each pixel evaluates
its source at t_em=t_obs-delay, with delays supplied by the Kerr ray tracer.

Primary sources:
* Page & Thorne (1974), doi:10.1086/152990, Eqs. 11b, 12, 15:
  https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P
* Schnittman & Bertschinger (2004), doi:10.1086/383180, kinematic hot spots:
  https://arxiv.org/abs/astro-ph/0309458
* Cardenas-Avendano, Lupsasca & Zhu (2023), doi:10.1103/PhysRevD.107.043030,
  retarded equatorial emission and circular orbital frequency (not their
  stochastic field implementation): https://arxiv.org/abs/2211.07469
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import numpy as np

from kerr_emission import isco_radius, orbital_quantities
from nasa_emission import ShearingDiskEmission as SchwarzschildShearingDiskEmission
from nasa_emission import PreparedEmission as SchwarzschildPreparedEmission


@dataclass(frozen=True)
class StructureParameters:
    spin_a_star: float = 0.9
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
        return isco_radius(self.spin_a_star)

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
    """Kerr flux shape only; no dimensional disk parameters are constructed."""

    def __init__(self, spin_a_star=0.9, outer_radius_rg=28.0, samples=65537):
        self.spin_a_star = float(spin_a_star)
        self.rin = isco_radius(spin_a_star)
        self.rout = float(outer_radius_rg)
        if not np.isfinite(self.rout) or self.rout <= self.rin or samples < 1025:
            raise ValueError("invalid integration interval or sample count")
        self.samples = int(samples)
        # Resolve the exact quadratic stress-free edge behavior with a grid
        # logarithmic in r-rISCO, rather than wasting samples at large radii.
        delta = np.geomspace(1e-7*self.rin, self.rout-self.rin, self.samples-1)
        rr = np.concatenate(([self.rin], self.rin+delta))
        q = orbital_quantities(rr, self.spin_a_star)
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
    """Drop-in relative bolometric source for a Kerr disk with a chosen spin.

    evaluate(r,phi,t_em), prepare(r,phi,delay).evaluate(t_obs), radial_profile(r),
    and metadata() retain the Schwarzschild renderer's source interface.
    """

    def __init__(self, outer_radius_rg=28.0, loop_period_tg=800.0,
                 contrast=1.6, seed=13326, spin_a_star=0.9, **structure_options):
        self.parameters = StructureParameters(
            spin_a_star=spin_a_star, outer_radius_rg=outer_radius_rg,
            loop_period_tg=loop_period_tg, contrast=contrast, seed=seed,
            **structure_options)
        p = self.parameters
        self.rin = p.inner_radius_rg
        self.flux_model = DimensionlessPageThorne(p.spin_a_star, p.outer_radius_rg)
        rng = np.random.default_rng(p.seed)
        spacing = math.log(p.outer_radius_rg/self.rin)/(p.radial_bands-1)
        centers = np.linspace(math.log(self.rin), math.log(p.outer_radius_rg), p.radial_bands)
        centers[1:-1] += rng.uniform(-0.22, 0.22, p.radial_bands-2)*spacing
        self.log_centers = centers
        self.radii = np.exp(centers)
        self.lifetimes = np.minimum(
            p.lifetime_orbits*2*math.pi*(self.radii**1.5+p.spin_a_star),
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
        return 1.0/(r**1.5+self.parameters.spin_a_star)

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
            "radial_profile": "zero-torque Kerr Page-Thorne flux; radial peak=1",
            "radial_quadrature": "cumulative trapezoid with analytic circular-orbit derivatives",
            "radial_quadrature_samples": self.flux_model.samples,
            "radial_lookup": "log flux versus log(r-rISCO); quadratic unresolved edge limit",
            "angular_velocity_tg_inverse": "Omega(r)=1/(r^(3/2)+a*), never quantized",
            "dimensional_normalization": "none; mass, accretion rate, temperature and spectrum unspecified",
            "surface_emission": "two identical opaque, isotropically emitting faces",
            "flux_equation": "F proportional to -Omega'/(r*(E-Omega*L)^2) integral_ISCO^r (E-Omega*L)*L' dr",
        })
        return data


class PreparedEmission(SchwarzschildPreparedEmission):
    """Reuse identical patch weights and retarded times, with Kerr advection."""

    def __init__(self, source, radius_rg, phi_rad, delay_tg):
        super().__init__(source, radius_rg, phi_rad, delay_tg)
        r = np.broadcast_arrays(np.asarray(radius_rg, dtype=np.float64),
                                np.asarray(phi_rad, dtype=np.float64),
                                np.asarray(delay_tg, dtype=np.float64))[0].ravel()
        self.omega = source.angular_velocity(r[self.valid_indices])
