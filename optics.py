"""Atmosphere and display optics for the relativistic thermal disk.

The atmosphere is an explicitly approximate, color-corrected Planck model;
this module does not solve vertical radiative transfer. The optional Airy PSF
belongs to a hypothetical camera, not to the black hole. It uses a single
reference wavelength for every RGB channel (an achromatic approximation).

Primary sources:
* Li, Zimmerman, Narayan & McClintock (2005), ApJS 157, 335,
  https://arxiv.org/abs/astro-ph/0411583, Eq. D20 and E5--E7.
* Done, Davis, Jin, Blaes & Ward (2012), MNRAS 420, 1848,
  https://arxiv.org/abs/1107.5429, Eq. 1--2 and their radial application.
* Airy (1835), Transactions of the Cambridge Philosophical Society 5, 283:
  the clear circular aperture intensity [2 J1(pi D theta/lambda)/(...) ]^2.
  Formula independently documented by Astropy's AiryDisk2D:
  https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.AiryDisk2D.html

The tone curve is an original display choice, not an astrophysical equation.
A fixed exposure and common scalar compression preserve in-gamut RGB ratios.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np

# Temperature in kelvin to k_B T in keV, from exact SI constants.
KELVIN_TO_KEV = 1.380649e-23 / 1.602176634e-16


def limb_factor(mu_em):
    """Flux-normalized Li et al. (2005) D20 scattering limb approximation.

    mu_em is the cosine of the emission angle in the orbiting gas rest frame.
    Upsilon=1/2+3/4 mu; 2 integral_0^1 Upsilon(mu) mu dmu=1.
    This is the KERRBB angular law, not a fully computed AGN atmosphere.
    """
    mu = np.asarray(mu_em, dtype=np.float64)
    if np.any(~np.isfinite(mu)) or np.any(mu < -1e-8) or np.any(mu > 1+1e-8):
        raise ValueError('emission-angle cosines must be in [0,1]')
    return .5 + .75 * np.clip(mu, 0, 1)


def color_correction(temperature_k):
    """Continuous form of Done et al. (2012) approximate atmosphere fit.

    f_col=1 at T<=30000 K. Above that, take the smaller of the rising
    (T/30000)^0.82 branch and the saturated (72/T_keV)^(1/9) branch.
    The published rounded fits differ slightly at their quoted 100000 K
    boundary; joining at their actual intersection (~101749 K) prevents
    an artificial intensity step without changing either fitted branch.
    At unrealistically high temperatures cap the minimum at unity.
    """
    temp = np.asarray(temperature_k, dtype=np.float64)
    if np.any(~np.isfinite(temp)) or np.any(temp < 0):
        raise ValueError('temperatures must be finite and nonnegative')
    safe = np.maximum(temp, 1e-30)
    rising = (safe / 3e4)**.82
    saturated = (72 / (safe * KELVIN_TO_KEV))**(1/9)
    return np.maximum(1., np.minimum(rising, saturated))


def atmosphere_rgb(temperature_k, g, mu_em, spectrum, *, hardening=True,
                   limb_darkening=True):
    """Observed CIE-integrated linear RGB radiance, before camera/display.

    I_nu,obs = Upsilon(mu) f_col^-4 B_nu(nu, g f_col T_eff).
    The f_col^-4 dilution preserves the bolometric local flux. Relativistic
    intensity transfer is already in the shifted Planck spectrum: do not
    multiply this RGB by g^3 or g^4 a second time.
    """
    temp, shift, mu = np.broadcast_arrays(
        np.asarray(temperature_k, dtype=np.float64),
        np.asarray(g, dtype=np.float64), np.asarray(mu_em, dtype=np.float64))
    if np.any(~np.isfinite(shift)) or np.any(shift < 0):
        raise ValueError('frequency shifts must be finite and nonnegative')
    fcol = color_correction(temp) if hardening else np.ones_like(temp)
    angular = limb_factor(mu) if limb_darkening else np.ones_like(temp)
    return spectrum.lookup_rgb_linear(shift*fcol*temp) * (angular/fcol**4)[...,None]


def chromaticity_preserving_tonemap(linear_rgb, exposure=1.):
    """Fixed exposure, bounded output, and no independent channel clipping.

    For v=exposure*max(RGB,0), output=v/[1+max(v)]. Every channel has the
    same denominator. RGB ratios, hence XYZ chromaticity, are preserved for
    in-gamut colors. Negative out-of-gamut primaries alone are clipped before
    the curve; this rare display-gamut limitation is not a physical effect.
    Exposure must stay identical in every frame. Output is linear display RGB.
    """
    if not np.isfinite(exposure) or exposure <= 0:
        raise ValueError('exposure must be finite and positive')
    values = np.asarray(linear_rgb, dtype=np.float64)
    if values.shape[-1] != 3 or np.any(~np.isfinite(values)):
        raise ValueError('expected finite linear RGB')
    values = np.maximum(values, 0) * exposure
    return values / (1 + values.max(axis=-1, keepdims=True))


@lru_cache(maxsize=4)
def _quadrature(nodes):
    x, w = np.polynomial.legendre.leggauss(nodes)
    return (x+1)*math.pi/2, w/2


def bessel_j1(values):
    """J1 from its defining integral, using Gauss--Legendre quadrature.

    J1(x)=1/pi integral_0^pi cos(t-x sin(t)) dt. Kernel arguments are small
    enough for this direct quadrature; no copied approximation or SciPy needed.
    """
    a = np.asarray(values, dtype=np.float64)
    flat = a.ravel()
    # Nodes rise with the highest oscillation frequency in the integrand.
    nodes = max(128, int(np.ceil(float(np.max(np.abs(flat), initial=0))))+48)
    angles, weights = _quadrature(nodes)
    out = np.empty_like(flat)
    for start in range(0, len(flat), 512):
        end = min(start+512, len(flat))
        out[start:end] = np.cos(angles[None,:]-flat[start:end,None]*np.sin(angles)[None,:]) @ weights
    return out.reshape(a.shape)


@dataclass(frozen=True)
class CameraParameters:
    """Hypothetical ideal camera used only for the optional optical preview.

    Distance is in rg, making pixel_angle=fov_rg/(distance_rg*width).
    The default 10^6 rg distance is near the infinity-observer limit of the
    ray tracer. At 1080 pixels across 32 rg, D=20m gives first-null radius
    1.13 pixels at 550nm. The film is an angle scan, not observer travel.
    """
    aperture_m: float = 20.0
    wavelength_nm: float = 550.0
    observer_distance_rg: float = 1e6
    support_first_nulls: float = 40.0
    kernel_subsamples: int = 3

    def __post_init__(self):
        if min(self.aperture_m, self.wavelength_nm, self.observer_distance_rg,
               self.support_first_nulls) <= 0 or self.kernel_subsamples < 1:
            raise ValueError('camera parameters must be positive')

    def first_null_pixels(self, fov_rg, width):
        pixel_radian = fov_rg/(self.observer_distance_rg*width)
        return 1.2196698912665045*self.wavelength_nm*1e-9/self.aperture_m/pixel_radian

    def metadata(self, fov_rg, width):
        return {**asdict(self), 'first_null_pixels':self.first_null_pixels(fov_rg,width),
                'model':'unobstructed circular pupil; monochromatic Airy PSF applied to linear RGB',
                'limitation':'550nm reference wavelength applied equally to RGB; not full polychromatic camera',
                'boundary':'zero-padded linear convolution; flux beyond cropped image is discarded',
                'kernel':'pixel-integrated finite support, renormalized to unit sum',
                'astrophysical_effect':False}


def airy_kernel(first_null_pixels, support_first_nulls=40., subsamples=3):
    """Pixel-integrated positive unit-sum Airy kernel, with a finite cutoff.

    The finite disk retains >99% of the analytical encircled energy for the
    default support. We renormalize the retained energy, so convolution is
    flux-conserving before output cropping. No extra Gaussian/bloom is added.
    """
    if first_null_pixels <= 0 or support_first_nulls < 2 or subsamples < 1:
        raise ValueError('invalid Airy kernel parameters')
    radius = max(2, int(math.ceil(first_null_pixels*support_first_nulls)))
    y,x = np.mgrid[-radius:radius+1,-radius:radius+1]
    kernel = np.zeros(x.shape, dtype=np.float64)
    offsets = (np.arange(subsamples)+.5)/subsamples-.5
    for dy in offsets:
        for dx in offsets:
            r = np.hypot(x+dx,y+dy)
            arg = 3.8317059702075125*r/first_null_pixels
            j1 = bessel_j1(arg)
            ratio = np.ones_like(arg)
            np.divide(2*j1,arg,out=ratio,where=arg>1e-10)
            intensity = ratio**2
            intensity[r > radius+.5] = 0
            kernel += intensity
    kernel /= kernel.sum()
    return kernel


def _next_fast_length(n):
    while True:
        factor=n
        for prime in (2,3,5):
            while factor%prime==0:
                factor//=prime
        if factor==1:
            return n
        n+=1


class AiryCamera:
    """Reusable linear-light camera convolution, independent of GR tracing."""
    def __init__(self, width, height, fov_rg=32., parameters=None):
        self.width,self.height = int(width),int(height)
        self.parameters = parameters or CameraParameters()
        first = self.parameters.first_null_pixels(fov_rg,width)
        self.kernel = airy_kernel(first,self.parameters.support_first_nulls,
                                  self.parameters.kernel_subsamples)
        self.metadata = self.parameters.metadata(fov_rg,width)
        self.metadata['kernel_size_pixels'] = list(self.kernel.shape)
        extent=self.kernel.shape[0]-1
        self.fft_shape=(_next_fast_length(height+extent),_next_fast_length(width+extent))
        self.kernel_fft=np.fft.rfft2(self.kernel,s=self.fft_shape)

    def apply(self, linear_rgb):
        values=np.asarray(linear_rgb,dtype=np.float64)
        if values.shape != (self.height,self.width,3) or np.any(~np.isfinite(values)):
            raise ValueError('RGB dimensions do not match camera')
        output=np.empty_like(values)
        origin=self.kernel.shape[0]//2
        for channel in range(3):
            full=np.fft.irfft2(np.fft.rfft2(values[...,channel],s=self.fft_shape)*self.kernel_fft,
                              s=self.fft_shape)
            output[...,channel]=full[origin:origin+self.height,origin:origin+self.width]
        # A positive input and kernel cannot create negative light. Permit
        # legitimate negative out-of-gamut RGB from the spectral transform.
        if np.min(values)>=0:
            np.maximum(output,0,out=output)
        return output


def diagnostics():
    """Numerical checks of physical normalization and display invariants."""
    mu=np.linspace(0,1,10001)
    angular_flux=float(2*np.trapezoid(limb_factor(mu)*mu,mu))
    temps=np.geomspace(1000,1e7,10000)
    f=color_correction(temps)
    # Analytical thermal integral scales T^4: dilution must cancel f_col^4.
    conservation=float(np.max(np.abs((f*temps)**4/f**4/temps**4-1)))
    rng=np.random.default_rng(9)
    rgb=10**rng.uniform(-8,7,(2000,3))
    mapped=chromaticity_preserving_tonemap(rgb,2.0)
    old_ratio=rgb/rgb.sum(axis=1,keepdims=True)
    new_ratio=mapped/mapped.sum(axis=1,keepdims=True)
    hue_error=float(np.max(np.abs(old_ratio-new_ratio)))
    j1_zero=float(abs(bessel_j1(np.array([3.8317059702075125]))[0]))
    # Enough space for the compact kernel around a central impulse.
    camera=AiryCamera(192,192,32,CameraParameters(aperture_m=4.,support_first_nulls=24))
    impulse=np.zeros((192,192,3));impulse[96,96]=[1.,2.,3.]
    convolved=camera.apply(impulse)
    flux_error=float(np.max(np.abs(convolved.sum(axis=(0,1))-[1,2,3])))
    # Linearity prevents intensity-threshold artifacts in faint/bright rings.
    linearity=float(np.max(np.abs(camera.apply(impulse*3)-convolved*3)))
    checks={
        'limb_angular_flux_normalization':angular_flux,
        'color_correction_bolometric_relative_error':conservation,
        'tone_curve_chromaticity_max_error':hue_error,
        'Airy_first_null_J1_absolute_error':j1_zero,
        'PSF_kernel_sum':float(camera.kernel.sum()),
        'PSF_impulse_flux_absolute_error':flux_error,
        'PSF_linearity_absolute_error':linearity,
        'camera':camera.metadata,
    }
    checks['passed']=bool(abs(angular_flux-1)<1e-7 and conservation<1e-12 and
                          hue_error<1e-12 and j1_zero<1e-10 and flux_error<1e-10 and
                          linearity<1e-12)
    if not checks['passed']:
        raise AssertionError(checks)
    return checks


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    report=diagnostics()
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
