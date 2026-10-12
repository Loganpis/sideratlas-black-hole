"""Zero-torque Page--Thorne emission for a circular, equatorial Kerr disk.

This module changes the physical disk solution with spin; it does not tint a
Schwarzschild image. Radius is in GM/c^2, spin is a*=Jc/(GM^2), and orbits
have positive angular velocity. Positive spin therefore denotes prograde
orbits. In geometrized units G=M=c=1, Page & Thorne (1974), Eqs. (11b),
(12), and (15), give the local flux from ONE face of the disk:

 F(r) = Mdot*c^6/(4*pi*G^2*M^2) * [-Omega'(r)/(r*(E-Omega*L)^2)]
        * integral_ISCO^r (E-Omega*L)*L'(r') dr'.

Analytic orbital functions and their derivatives are used; only the final
radial integral is numerical. A dense logarithmic distance-from-ISCO grid
resolves its quadratic zero at the stress-free inner edge. The cumulative
trapezoid integral and log-log interpolation are checked against the exact
Schwarzschild solution, direct quadrature, and doubled grid resolution.

Each surface element has T_eff=(F/sigma_SB)^(1/4). For the f_col=1 thermal
assumption, the existing emission.ThermalSpectrum maps T_obs=g*T_eff to
CIE-visible color. Any angular emissivity law must be applied separately
and normalized to conserve this local flux. Do not multiply spectral RGB
by g^3 or g^4 a second time.

This is a steady, razor-thin model with zero ISCO stress. It omits returning
radiation and photon-capture feedback on heating, an atmosphere calculation,
a corona, and magnetic/plasma evolution. It is not the full KERRBB model.

Primary source actually used:
 Page & Thorne (1974), ApJ 191, 499--506, doi:10.1086/152990,
 https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P
The existing emission module documents the Planck and colorimetry sources.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import argparse
import json
import math
from pathlib import Path

import numpy as np

from emission import C, G, GM_SUN, PROTON_MASS, SIGMA_SB, THOMSON_CROSS_SECTION


def isco_radius(spin_a_star: float) -> float:
    """ISCO for the positive-angular-velocity orbit branch, Page--Thorne 15l."""
    a = float(spin_a_star)
    if not math.isfinite(a) or abs(a) >= 1.0:
        raise ValueError("spin_a_star must be finite and strictly between -1 and 1")
    z1 = 1.0 + np.cbrt(1.0-a*a)*(np.cbrt(1.0+a)+np.cbrt(1.0-a))
    z2 = math.sqrt(3.0*a*a + z1*z1)
    return float(3.0+z2-np.sign(a)*math.sqrt((3.0-z1)*(3.0+z1+2.0*z2)))


def orbital_quantities(radius_rg, spin_a_star: float) -> dict:
    """Specific E, L, Omega and analytic radial derivatives, G=M=c=1.

    Inputs must lie on the physical timelike circular-orbit branch. The
    disk solution only evaluates this at and outside its ISCO. L is in
    GM/c; Omega is in c^3/(GM); primes use dimensionless r/rg.
    """
    r = np.asarray(radius_rg, dtype=np.float64)
    a = float(spin_a_star)
    if not math.isfinite(a) or abs(a) >= 1.0:
        raise ValueError("invalid dimensionless Kerr spin")
    if np.any(~np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("orbital radii must be finite and positive")
    root_r = np.sqrt(r)
    metric_c = 1.0-3.0/r+2.0*a/r**1.5
    if np.any(metric_c <= 0):
        raise ValueError("radius is outside the timelike circular-orbit branch")
    sqrt_c = np.sqrt(metric_c)
    energy_numerator = 1.0-2.0/r+a/r**1.5
    angular_numerator = 1.0-2.0*a/r**1.5+a*a/r**2
    dc = 3.0/r**2-3.0*a/r**2.5
    de_num = 2.0/r**2-1.5*a/r**2.5
    dl_num = 3.0*a/r**2.5-2.0*a*a/r**3
    energy = energy_numerator/sqrt_c
    angular_momentum = root_r*angular_numerator/sqrt_c
    omega = 1.0/(r**1.5+a)
    denergy = (de_num-0.5*energy_numerator*dc/metric_c)/sqrt_c
    dangular_momentum = root_r/sqrt_c * (
        angular_numerator/(2.0*r)+dl_num-0.5*angular_numerator*dc/metric_c)
    domega = -1.5*root_r/(r**1.5+a)**2
    # This equivalent expression avoids subtraction of E and Omega*L.
    energy_minus_omega_l = sqrt_c/(1.0+a/r**1.5)
    return {"energy": energy, "angular_momentum": angular_momentum,
            "omega": omega, "denergy_dr": denergy,
            "dangular_momentum_dr": dangular_momentum, "domega_dr": domega,
            "energy_minus_omega_l": energy_minus_omega_l}


@dataclass(frozen=True)
class KerrDiskParameters:
    """lambda=eta*Mdot*c^2/L_Edd, eta=1-E_ISCO (capture uncorrected)."""

    mass_solar: float = 1.0e8
    eddington_ratio: float = 0.1
    spin_a_star: float = 0.9

    def __post_init__(self):
        if not math.isfinite(self.mass_solar) or self.mass_solar <= 0:
            raise ValueError("mass_solar must be finite and positive")
        if not math.isfinite(self.eddington_ratio) or self.eddington_ratio <= 0:
            raise ValueError("eddington_ratio must be finite and positive")
        isco_radius(self.spin_a_star)

    @property
    def gm_si(self):
        return GM_SUN*self.mass_solar

    @property
    def mass_kg(self):
        return self.gm_si/G

    @property
    def rg_m(self):
        return self.gm_si/C**2

    @property
    def tg_s(self):
        return self.gm_si/C**3

    @property
    def isco_rg(self):
        return isco_radius(self.spin_a_star)

    @property
    def horizon_rg(self):
        return 1.0+math.sqrt(1.0-self.spin_a_star**2)

    @property
    def efficiency(self):
        return 1.0-float(orbital_quantities(self.isco_rg, self.spin_a_star)["energy"])

    @property
    def eddington_luminosity_w(self):
        return 4.0*math.pi*self.gm_si*PROTON_MASS*C/THOMSON_CROSS_SECTION

    @property
    def accretion_rate_kg_s(self):
        return self.eddington_ratio*self.eddington_luminosity_w/(self.efficiency*C**2)

    @property
    def flux_normalization_w_m2(self):
        return self.accretion_rate_kg_s*C**6/(4.0*math.pi*self.gm_si**2)

    def metadata(self):
        return {**asdict(self), "mass_kg": self.mass_kg, "rg_m": self.rg_m,
                "tg_s": self.tg_s, "isco_rg": self.isco_rg,
                "horizon_rg": self.horizon_rg,
                "efficiency_binding_energy": self.efficiency,
                "eddington_luminosity_w": self.eddington_luminosity_w,
                "accretion_rate_kg_s": self.accretion_rate_kg_s,
                "flux_normalization_w_m2": self.flux_normalization_w_m2,
                "eddington_ratio_definition": "eta*Mdot*c^2/L_Edd; no photon capture correction",
                "spectrum": "local Planck base; renderer may apply a separately documented atmosphere",
                "inner_boundary": "stress-free ISCO; no emission inside the ISCO"}


class KerrEmission:
    """Physical emitted flux and temperature; construction builds one tiny LUT."""

    def __init__(self, parameters: KerrDiskParameters | None = None,
                 samples: int = 65537, maximum_radius_rg: float = 1.0e8):
        self.parameters = parameters or KerrDiskParameters()
        self.rin = self.parameters.isco_rg
        if samples < 1025 or maximum_radius_rg <= self.rin:
            raise ValueError("invalid radial integration grid")
        self.samples = int(samples)
        self.maximum_radius_rg = float(maximum_radius_rg)
        # The log offset handles flux proportional to (r-r_ISCO)^2. Starting
        # at 1e-7*r_ISCO keeps analytic derivative cancellation negligible.
        delta = np.geomspace(1.0e-7*self.rin, maximum_radius_rg-self.rin,
                             self.samples-1)
        self._log_delta = np.log(delta)
        rr = np.concatenate(([self.rin], self.rin+delta))
        orbital = orbital_quantities(rr, self.parameters.spin_a_star)
        integrand = orbital["energy_minus_omega_l"]*orbital["dangular_momentum_dr"]
        integrand[0] = 0.0  # L'(ISCO)=0, removing roundoff at the endpoint.
        if np.min(integrand) < -1e-10:
            raise ValueError("negative angular-momentum derivative outside the ISCO")
        integral = np.concatenate(([0.0], np.cumsum(
            0.5*(integrand[:-1]+integrand[1:])*np.diff(rr))))
        shape = (-orbital["domega_dr"] /
                 (rr*orbital["energy_minus_omega_l"]**2))*integral
        if np.any(shape[1:] <= 0):
            raise ValueError("nonpositive Page--Thorne flux beyond the inner edge")
        self._log_shape = np.log(shape[1:])

    def flux_shape(self, radius_rg):
        """Dimensionless flux factor; zero inside ISCO or for zero sentinels.

        Radii beyond maximum_radius_rg raise instead of silently clamping.
        The unresolved <1e-7*r_ISCO interval uses the exact quadratic edge
        scaling; it is far smaller than a rendered pixel's radial footprint.
        """
        r = np.asarray(radius_rg, dtype=np.float64)
        if np.any(~np.isfinite(r)):
            raise ValueError("radii must be finite")
        if np.any(r > self.maximum_radius_rg):
            raise ValueError("requested radius exceeds emission integration grid")
        delta = r-self.rin
        log_delta = np.log(np.maximum(delta, np.finfo(float).tiny))
        log_shape = np.interp(log_delta, self._log_delta, self._log_shape)
        log_shape = np.where(log_delta < self._log_delta[0],
                             self._log_shape[0]+2.0*(log_delta-self._log_delta[0]),
                             log_shape)
        return np.where(delta > 0, np.exp(log_shape), 0.0)

    def flux(self, radius_rg):
        """Local bolometric energy flux from ONE disk face, W m^-2."""
        return self.parameters.flux_normalization_w_m2*self.flux_shape(radius_rg)

    def temperature(self, radius_rg):
        """Effective surface temperature, K, before atmosphere corrections."""
        return (self.flux(radius_rg)/SIGMA_SB)**0.25

    temp = temperature
    effective_temperature = temperature

    def observed_temperature(self, radius_rg, frequency_shift):
        g = np.asarray(frequency_shift, dtype=np.float64)
        if np.any(~np.isfinite(g)) or np.any(g < 0):
            raise ValueError("frequency shifts must be finite and nonnegative")
        return g*self.temperature(radius_rg)

    def bolometric_intensity(self, radius_rg, frequency_shift):
        """Isotropic-equivalent observed radiance; angular law not included."""
        g = np.asarray(frequency_shift, dtype=np.float64)
        if np.any(~np.isfinite(g)) or np.any(g < 0):
            raise ValueError("frequency shifts must be finite and nonnegative")
        return g**4*self.flux(radius_rg)/math.pi

    def metadata(self):
        return {"parameters": self.parameters.metadata(),
                "solution": "Page--Thorne 1974, Eqs. 11b, 12, 15",
                "integral_grid_samples": self.samples,
                "integral_maximum_radius_rg": self.maximum_radius_rg,
                "quadrature": "cumulative trapezoid; analytic orbital derivatives",
                "lookup": "log flux versus log distance from ISCO",
                "omissions": ["returning-radiation heating", "ISCO stress",
                              "finite thickness", "disk-atmosphere solution",
                              "plunging emission", "corona", "GRMHD evolution"]}


def diagnostics() -> dict:
    """Meaningful numerical checks against independent limits/quadratures."""
    from emission import DiskParameters, emitted_flux

    schwarzschild = KerrEmission(KerrDiskParameters(spin_a_star=0.0))
    test_r = np.geomspace(6.001, 1.0e6, 10001)
    reference = emitted_flux(test_r, DiskParameters())
    schwarz_error = float(np.max(np.abs(schwarzschild.flux(test_r)/reference-1.0)))
    if schwarz_error > 1.0e-5:
        raise AssertionError("Kerr a=0 flux fails its exact Schwarzschild limit")

    disk = KerrEmission()
    refined = KerrEmission(samples=2*disk.samples-1)
    radii = disk.rin+np.geomspace(1.0e-5, 1000.0, 10001)
    convergence_error = float(np.max(np.abs(disk.flux(radii)/refined.flux(radii)-1.0)))
    if convergence_error > 1.0e-5:
        raise AssertionError("Kerr flux quadrature has not converged")

    # Integrate at individual test radii on a different linear grid, without
    # the cumulative lookup or its interpolation.
    quadrature_errors = []
    for radius in (disk.rin+0.001, 3.0, 5.0, 10.0, 40.0, 100.0):
        rr = np.linspace(disk.rin, radius, 100001)
        q = orbital_quantities(rr, disk.parameters.spin_a_star)
        integral = np.trapezoid(q["energy_minus_omega_l"]*q["dangular_momentum_dr"], rr)
        reference_shape = (-q["domega_dr"][-1]/
                           (radius*q["energy_minus_omega_l"][-1]**2))*integral
        quadrature_errors.append(abs(float(disk.flux_shape(radius))/reference_shape-1.0))
    if max(quadrature_errors) > 1.0e-5:
        raise AssertionError("Kerr flux lookup disagrees with direct quadrature")

    orbital = orbital_quantities(np.geomspace(disk.rin+0.01, 1000.0, 1001), 0.9)
    identity_error = float(np.max(np.abs(orbital["denergy_dr"]/
                              (orbital["omega"]*orbital["dangular_momentum_dr"])-1.0)))
    if identity_error > 1.0e-10:
        raise AssertionError("analytic orbit derivatives violate dE=Omega*dL")
    if not np.all(disk.flux([0.0, disk.parameters.horizon_rg, disk.rin]) == 0.0):
        raise AssertionError("nonzero flux at or inside the stress-free ISCO")
    far_radius = 1.0e8
    newtonian_leading_shape = 1.5/far_radius**3
    far_ratio = float(disk.flux_shape(far_radius)/newtonian_leading_shape)
    if abs(far_ratio-1.0) > 0.001:
        raise AssertionError("flux fails its Newtonian large-radius limit")
    peak_grid = np.linspace(disk.rin, 20.0, 100001)
    peak_radius = float(peak_grid[np.argmax(disk.flux(peak_grid))])
    return {"passed": True,
            "schwarzschild_closed_form_max_relative_error": schwarz_error,
            "doubling_grid_max_relative_error": convergence_error,
            "independent_linear_quadrature_max_relative_error": max(quadrature_errors),
            "dE_equals_Omega_dL_max_relative_error": identity_error,
            "flux_to_leading_newtonian_ratio_at_1e8_rg": far_ratio,
            "flux_peak_radius_rg": peak_radius,
            "effective_temperature_peak_k": float(disk.temperature(peak_radius)),
            "emission": disk.metadata()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(diagnostics(), indent=2)
    if args.output:
        args.output.write_text(report+"\n")
    print(report)
