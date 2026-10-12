"""Physical emission for a Schwarzschild, zero-torque, optically thick thin disk.

The Page--Thorne radiative-flux solution is evaluated in closed form for a*=0.
The disk has its stress-free inner edge at r=6 GM/c^2. Each surface element is
an isotropic blackbody with T_eff=(F/sigma_SB)^1/4 (color correction = 1).
For an observer at infinity, Liouville transfer gives

    I_nu,obs(nu) = g^3 B_nu(nu/g, T_eff) = B_nu(nu, g*T_eff)
    I_bol,obs    = g^4 F/pi.

These are idealized thin-disk assumptions, not GRMHD, a disk-atmosphere model,
or a prediction for any particular observed black hole. Returning radiation,
limb darkening, scattering/coronae, magnetic fields, and photon-capture feedback
on the disk temperature are omitted. Nothing in this module invents texture.

Sources actually used:
* Page & Thorne (1974), ApJ 191, 499, doi:10.1086/152990.
* Luminet (1979), A&A 75, 228, section 4 (flux and frequency transfer).
* Mueller & Frauendiener (2012), arXiv:1206.4259, section 2.3, Eq. (11).
* CIE (2019), 1931 2-degree observer, doi:10.25039/CIE.DS.xvudnb9b.
* ICC, Interpretation of sRGB, section 7 (XYZ-to-linear-sRGB matrix),
  https://registry.color.org/rgb-registry/files/sRGB.pdf.

Only NumPy and the Python standard library are required. The accompanying
CIE table is loaded locally; running this module never accesses the network.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
import json
import math

import numpy as np


# SI constants. h, k_B and c are exact SI definitions; G, m_p and sigma_T
# use CODATA 2018 values. GM_sun is the IAU nominal solar mass parameter.
C = 299_792_458.0
G = 6.67430e-11
H = 6.62607015e-34
KB = 1.380649e-23
SIGMA_SB = 2.0 * math.pi**5 * KB**4 / (15.0 * H**3 * C**2)
PROTON_MASS = 1.67262192369e-27
THOMSON_CROSS_SECTION = 6.6524587321e-29
GM_SUN = 1.3271244e20
ETA_SCHWARZSCHILD = 1.0 - math.sqrt(8.0 / 9.0)

CIE_MD5 = "17cca777db64b17170f06f67ce9d3ab7"
CIE_FILENAME = Path(__file__).resolve().parent / "data/CIE_xyz_1931_2deg.csv"
XYZ_TO_LINEAR_SRGB = np.array([
    [3.2406255, -1.5372080, -0.4986286],
    [-0.9689307, 1.8757561, 0.0415175],
    [0.0557101, -0.2040211, 1.0569959],
], dtype=np.float64)


@dataclass(frozen=True)
class DiskParameters:
    """Accretion normalization; ``eddington_ratio`` is eta*Mdot*c^2/L_Edd.

    eta=1-sqrt(8/9) is the binding-energy efficiency at the Schwarzschild
    ISCO. This is a normalization convention, not the ray-traced luminosity
    reaching infinity after photon capture. L_Edd assumes ionized hydrogen.
    """

    mass_solar: float = 1.0e8
    eddington_ratio: float = 0.1

    def __post_init__(self):
        if not (np.isfinite(self.mass_solar) and self.mass_solar > 0):
            raise ValueError("mass_solar must be finite and positive")
        if not (np.isfinite(self.eddington_ratio) and self.eddington_ratio > 0):
            raise ValueError("eddington_ratio must be finite and positive")

    @property
    def gm_si(self) -> float:
        return GM_SUN * self.mass_solar

    @property
    def mass_kg(self) -> float:
        return self.gm_si / G

    @property
    def rg_m(self) -> float:
        return self.gm_si / C**2

    @property
    def tg_s(self) -> float:
        return self.gm_si / C**3

    @property
    def efficiency(self) -> float:
        return ETA_SCHWARZSCHILD

    @property
    def eddington_luminosity_w(self) -> float:
        return 4.0 * math.pi * self.gm_si * PROTON_MASS * C / THOMSON_CROSS_SECTION

    @property
    def accretion_rate_kg_s(self) -> float:
        return self.eddington_ratio * self.eddington_luminosity_w / (self.efficiency * C**2)

    @property
    def flux_normalization_w_m2(self) -> float:
        return 3.0 * C**6 * self.accretion_rate_kg_s / (8.0 * math.pi * self.gm_si**2)

    def metadata(self) -> dict:
        return {
            **asdict(self),
            "spin_a_star": 0.0,
            "mass_kg": self.mass_kg,
            "rg_m": self.rg_m,
            "tg_s": self.tg_s,
            "isco_rg": 6.0,
            "horizon_rg": 2.0,
            "efficiency_binding_energy": self.efficiency,
            "eddington_luminosity_w": self.eddington_luminosity_w,
            "accretion_rate_kg_s": self.accretion_rate_kg_s,
            "flux_normalization_w_m2": self.flux_normalization_w_m2,
            "spectrum": "isotropic Planck blackbody; color correction f_col=1",
            "eddington_ratio_definition": "eta*Mdot*c^2/L_Edd; does not correct for photon capture",
        }


DEFAULT_PARAMETERS = DiskParameters()


def nt_flux_shape(radius_rg):
    """Dimensionless Page--Thorne flux shape, zero at/below the ISCO.

    For x=r/rg >= 6, F=(3 c^6 Mdot / 8 pi G^2 M^2)*f(x), where

      f = J/[x^(5/2)(x-3)]
      J = sqrt(x)-sqrt(6) - sqrt(3)/2 *
          ln[((sqrt(x)-sqrt(3))/(sqrt(x)+sqrt(3))) /
              ((sqrt(6)-sqrt(3))/(sqrt(6)+sqrt(3)))] .

    The equivalent Page--Thorne integral is
      J=int_6^x (E-Omega*L) dL/dr dr
       =int_6^x (r-6)/[2 sqrt(r)(r-3)] dr.

    ``radius_rg`` may contain zero sentinels for rays that miss the disk.
    The physical emitting outer radius is determined by the ray tracer.
    """
    r = np.asarray(radius_rg, dtype=np.float64)
    rr = np.maximum(r, 6.0)
    s = np.sqrt(rr)
    a, b = math.sqrt(3.0), math.sqrt(6.0)
    delta = s - b
    # log1p evaluates the logarithmic ratio without dividing nearly equal
    # numbers. The remaining cancellation at the exact ISCO is safely zero.
    log_ratio = np.log1p(delta / (b-a)) - np.log1p(delta / (b+a))
    integral = delta - 0.5 * a * log_ratio
    shape = np.maximum(integral, 0.0) / (rr**2.5 * (rr-3.0))
    return np.where((r > 6.0) & np.isfinite(r), shape, 0.0)


def emitted_flux(radius_rg, parameters: DiskParameters = DEFAULT_PARAMETERS):
    """Flux emerging from one disk face, in local W m^-2."""
    return parameters.flux_normalization_w_m2 * nt_flux_shape(radius_rg)


def effective_temperature(radius_rg, parameters: DiskParameters = DEFAULT_PARAMETERS):
    """Local effective temperature in kelvin, from F=sigma_SB T^4."""
    return np.power(emitted_flux(radius_rg, parameters) / SIGMA_SB, 0.25)


def observed_temperature(radius_rg, g, parameters: DiskParameters = DEFAULT_PARAMETERS):
    """Equivalent observed blackbody temperature T_obs=g*T_eff, kelvin."""
    shift = np.asarray(g, dtype=np.float64)
    if np.any(~np.isfinite(shift)) or np.any(shift < 0):
        raise ValueError("frequency-shift factors must be finite and nonnegative")
    return shift * effective_temperature(radius_rg, parameters)


def bolometric_intensity(radius_rg, g, parameters: DiskParameters = DEFAULT_PARAMETERS):
    """Observed bolometric radiance g^4 F/pi, W m^-2 sr^-1."""
    return np.asarray(g, dtype=np.float64)**4 * emitted_flux(radius_rg, parameters) / math.pi


def planck_nu(frequency_hz, temperature_k):
    """B_nu in W m^-2 sr^-1 Hz^-1, with NumPy broadcasting."""
    nu, temp = np.broadcast_arrays(np.asarray(frequency_hz, dtype=np.float64),
                                   np.asarray(temperature_k, dtype=np.float64))
    if np.any(nu <= 0) or np.any(~np.isfinite(nu)):
        raise ValueError("frequencies must be finite and positive")
    if np.any(temp < 0) or np.any(~np.isfinite(temp)):
        raise ValueError("temperatures must be finite and nonnegative")
    exponent = H * nu / (KB * np.maximum(temp, np.finfo(float).tiny))
    denominator = np.expm1(np.minimum(exponent, 700.0))
    out = (2.0 * H / C**2) * nu**3 / denominator
    return np.where((temp > 0) & (exponent < 700.0), out, 0.0)


def planck_lambda(wavelength_m, temperature_k):
    """B_lambda in W m^-2 sr^-1 m^-1, with NumPy broadcasting."""
    wavelength, temp = np.broadcast_arrays(np.asarray(wavelength_m, dtype=np.float64),
                                           np.asarray(temperature_k, dtype=np.float64))
    if np.any(wavelength <= 0) or np.any(~np.isfinite(wavelength)):
        raise ValueError("wavelengths must be finite and positive")
    if np.any(temp < 0) or np.any(~np.isfinite(temp)):
        raise ValueError("temperatures must be finite and nonnegative")
    exponent = H * C / (wavelength * KB * np.maximum(temp, np.finfo(float).tiny))
    denominator = np.expm1(np.minimum(exponent, 700.0))
    out = (2.0 * H * C**2) / (wavelength**5 * denominator)
    return np.where((temp > 0) & (exponent < 700.0), out, 0.0)


class ThermalSpectrum:
    """Absolute CIE-weighted visible radiance mapped to linear sRGB primaries.

    Full 360--830 nm CIE table, 1 nm trapezoidal spectral quadrature.
    No per-pixel luminance, color, or maximum normalization is performed.
    XYZ = integral B_lambda(lambda,T) [xbar,ybar,zbar] d lambda.
    RGB = XYZ_TO_LINEAR_SRGB @ XYZ, in weighted W m^-2 sr^-1.
    Multiply by one fixed global exposure before tone mapping and encoding.

    These XYZ values omit the common 683 lm/W photometric scale since a fixed
    display exposure absorbs it. The Y component times 683 is luminance cd/m^2.
    This is spectral colorimetry, not a bolometric pseudo-color map. Because
    gravity/Doppler shift changes a Planck spectrum to T_obs=g*T_eff, do NOT
    multiply the resulting RGB by g^3 or g^4 again.

    A log-T lookup accelerates rendering; input values outside the table are
    calculated directly (rather than silently being clipped). Tiny negative
    out-of-sRGB-gamut channel values are returned as such for the display
    pipeline to handle. Most of the hot disk is in gamut; cooler pixels very
    close to the stress-free inner edge can require display-gamut clipping.
    """

    def __init__(self, cie_path=None, bins=4096, min_temperature=100.0,
                 max_temperature=1.0e6):
        if bins < 32 or not 0 < min_temperature < max_temperature:
            raise ValueError("invalid temperature lookup grid")
        path = Path(cie_path) if cie_path is not None else CIE_FILENAME
        data = path.read_bytes()
        if hashlib.md5(data).hexdigest() != CIE_MD5:
            raise ValueError("CIE table checksum disagrees with official dataset")
        table = np.loadtxt(path, delimiter=",")
        if table.shape != (471, 4) or not np.array_equal(table[:, 0], np.arange(360,831)):
            raise ValueError("unexpected CIE table wavelength grid")
        self.wavelength_m = table[:, 0] * 1.0e-9
        self.cmf = table[:, 1:4]
        self.log_temperatures = np.linspace(math.log(min_temperature),
                                           math.log(max_temperature), bins)
        self.temperatures = np.exp(self.log_temperatures)
        self.rgb_table = self._direct(self.temperatures)
        self.metadata = {
            "cie_source": "https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer",
            "cie_doi": "10.25039/CIE.DS.xvudnb9b",
            "cie_md5": CIE_MD5,
            "cie_wavelength_nm": [360, 830],
            "spectral_step_nm": 1,
            "lookup_bins": bins,
            "lookup_temperature_k": [min_temperature, max_temperature],
            "radiance_units": "CIE-weighted W m^-2 sr^-1; no per-pixel normalization",
            "display_matrix": "ICC sRGB interpretation, section 7",
        }

    def _direct(self, temperature_k):
        temperatures = np.asarray(temperature_k, dtype=np.float64)
        flat = temperatures.ravel()
        result = np.empty((flat.size, 3), dtype=np.float64)
        # Keep memory bounded for out-of-grid values and diagnostics.
        for start in range(0, flat.size, 2048):
            t = flat[start:start+2048]
            spectra = planck_lambda(self.wavelength_m[None, :], t[:, None])
            xyz = np.trapezoid(spectra[:, :, None] * self.cmf[None, :, :],
                               x=self.wavelength_m, axis=1)
            result[start:start+t.size] = xyz @ XYZ_TO_LINEAR_SRGB.T
        return result.reshape(temperatures.shape + (3,))

    def lookup_rgb_linear(self, temperature_k):
        temperatures = np.asarray(temperature_k, dtype=np.float64)
        if np.any(~np.isfinite(temperatures)) or np.any(temperatures < 0):
            raise ValueError("temperatures must be finite and nonnegative")
        flat = temperatures.ravel()
        log_t = np.log(np.maximum(flat, self.temperatures[0]))
        rgb = np.column_stack([np.interp(log_t, self.log_temperatures, self.rgb_table[:, c])
                               for c in range(3)])
        outside = (flat > 0) & ((flat < self.temperatures[0]) | (flat > self.temperatures[-1]))
        if np.any(outside):
            rgb[outside] = self._direct(flat[outside])
        rgb[flat == 0] = 0.0
        return rgb.reshape(temperatures.shape + (3,))

    # Compatibility spelling requested by the renderer coordinator.
    lookupRGBlinear = lookup_rgb_linear


def diagnostics(parameters: DiskParameters = DEFAULT_PARAMETERS) -> dict:
    """Independent integral, Planck-transfer, energy, and color checks."""
    radii = np.array([6.01, 6.1, 7.0, 9.55, 20.0, 40.0, 100.0])
    integral_shapes = []
    for radius in radii:
        rr = np.linspace(6.0, radius, 20001)
        e = (rr - 2.0) / np.sqrt(rr * (rr - 3.0))
        ell = rr / np.sqrt(rr - 3.0)
        omega = rr**-1.5
        dell_dr = (rr - 6.0) / (2.0 * (rr - 3.0)**1.5)
        integral = np.trapezoid((e - omega*ell) * dell_dr, rr)
        integral_shapes.append(integral / (radius**2.5 * (radius-3.0)))
    closed = nt_flux_shape(radii)
    integral_error = float(np.max(np.abs(closed / integral_shapes - 1)))
    if integral_error > 1.0e-7:
        raise AssertionError("closed-form flux disagrees with the Page--Thorne integral")
    peak_grid = np.linspace(6.0, 15.0, 100001)
    peak_index = int(np.argmax(nt_flux_shape(peak_grid)))
    peak_r = float(peak_grid[peak_index])
    if abs(peak_r-9.55) > 0.01:
        raise AssertionError("flux maximum disagrees with Luminet / Mueller")
    if not np.all(nt_flux_shape([0.0, 2.0, 5.99, 6.0]) == 0):
        raise AssertionError("nonzero stress-free flux inside ISCO")
    # All-frequency transfer and Stefan--Boltzmann are independent of the
    # CIE-visible-band calculation. Integrate Planck's x=h*nu/kT form.
    frequencies = np.geomspace(1.0e7, 1.0e19, 100001)
    test_temperature = 50000.0
    integrated = float(math.pi*np.trapezoid(planck_nu(frequencies, test_temperature), frequencies))
    sb_expected = SIGMA_SB * test_temperature**4
    sb_error = abs(integrated/sb_expected - 1)
    if sb_error > 2e-7:
        raise AssertionError("Planck integral disagrees with Stefan--Boltzmann")
    g = 0.68
    nu = np.geomspace(1e13, 1e17, 100)
    transfer_error = float(np.max(np.abs(
        g**3 * planck_nu(nu/g, test_temperature) /
        planck_nu(nu, g*test_temperature) - 1)))
    if transfer_error > 2e-12:
        raise AssertionError("Planck temperature shift disagrees with Liouville transfer")
    spectrum = ThermalSpectrum()
    test_t = np.geomspace(2000, 5e5, 151)
    exact_rgb = spectrum._direct(test_t)
    lookup_rgb = spectrum.lookup_rgb_linear(test_t)
    color_error = float(np.max(np.abs(lookup_rgb/exact_rgb-1)))
    if color_error > 0.001:
        raise AssertionError("color lookup exceeds 0.1 percent interpolation error")
    return {
        "parameters": parameters.metadata(),
        "page_thorne_integral_max_relative_error": integral_error,
        "flux_peak_radius_rg": peak_r,
        "flux_peak_shape": float(nt_flux_shape(peak_r)),
        "flux_peak_w_m2": float(emitted_flux(peak_r, parameters)),
        "temperature_peak_k": float(effective_temperature(peak_r, parameters)),
        "stefan_boltzmann_integral_relative_error": sb_error,
        "planck_liouville_relative_error": transfer_error,
        "cie_lookup_max_relative_error_2000_to_500000k": color_error,
        "rgb_at_peak_temperature": spectrum.lookup_rgb_linear(effective_temperature(peak_r, parameters)).tolist(),
        "spectrum": spectrum.metadata,
    }


if __name__ == "__main__":
    print(json.dumps(diagnostics(), indent=2))
