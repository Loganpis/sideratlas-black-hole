# Color, atmosphere and optional camera optics

The thermal disk has a calculated visible spectrum. There is no universal “true color” of every black hole: the black hole emits no visible surface light here, while the disk spectrum depends on its mass, accretion rate, spin and atmosphere. This model predicts its visible color conditional on those choices, using the supplied CIE 1931 observer table and an sRGB display transform.

## Implemented atmosphere approximation

Li et al. (2005), *Multi-Temperature Blackbody Spectrum of a Thin Accretion Disk Around a Kerr Black Hole: Model Computations and Comparison with Observations*, ApJS 157, 335, [arXiv:astro-ph/0411583](https://arxiv.org/abs/astro-ph/0411583), gives the gray limb law in Eq. D20 and the diluted color-corrected thermal spectrum in Eqs. E5–E7. We use

- `Upsilon(mu) = 1/2 + 3 mu/4`, normalized so `2 integral_0^1 Upsilon(mu) mu dmu = 1`.
- `I_nu,obs = Upsilon(mu) f_col^-4 B_nu(nu, g f_col T_eff)`.

The emission cosine must be evaluated in the circularly orbiting gas frame. At an equatorial crossing of a Kerr null geodesic, `mu = g sqrt(Q)/r`, for Carter constant `Q/E_infinity^2` in geometric units. This follows from contracting the photon momentum with the disk unit normal and dividing by its locally measured energy. It is not the observer's inclination cosine.

Done et al. (2012), *Intrinsic disc emission and the soft X-ray excess in active galactic nuclei*, MNRAS 420, 1848, [arXiv:1107.5429](https://arxiv.org/abs/1107.5429), Section 2, Eqs. 1–2, supplies an approximate temperature-dependent color correction. Their discussion explicitly applies the formulae at each disk radius. We use `max(1, min((T/30000)^0.82, (72/T_keV)^(1/9)))`. The two published rounded branches are joined at their actual intersection, 101749 K, rather than making a small discontinuity at the rounded 100000 K boundary. The `f_col^-4` factor is essential: it preserves bolometric thermal flux while redistributing spectral energy. Applying `g^3` or `g^4` to the resulting visible RGB again would double-count relativistic transfer.

These formulae are an approximate atmosphere prescription, not a solution for its vertical structure, frequency-dependent absorption edges, polarization, or returning radiation. The gray scattering limb law is a standard thin-disk approximation; real AGN atmospheres can have more complex angular emission, especially in the cooler outer disk.

## Optional camera

`AiryCamera` models the diffraction of an ideal unobstructed circular pupil. Its intensity point-spread function is

`P(theta) proportional to [2 J1(pi D theta/lambda)/(pi D theta/lambda)]^2`.

This is the circular-aperture diffraction result of Airy (1835), *On the Diffraction of an Object-glass with Circular Aperture*, Transactions of the Cambridge Philosophical Society 5, 283. The formula and first-null constant are also documented by the [Astropy AiryDisk2D reference](https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.AiryDisk2D.html). The implementation independently evaluates the defining Bessel integral with Gauss–Legendre quadrature; it does not depend on Astropy or SciPy.

The hypothetical camera uses a 20 m clear aperture at distance 1 million gravitational radii. At 1080 pixels across 32 gravitational radii, the 550 nm first dark Airy ring is 1.13 pixels from the point-source center. The infinity-observer ray-tracing approximation is appropriate at this distance. The same 550 nm reference kernel is applied to all three linear RGB channels: **this is an explicitly achromatic approximation, not full wavelength-resolved diffraction**. It is an optional camera view, not an extra black-hole light source or gravitational halo. No spikes, glare threshold, arbitrary bloom, or lens flare are added.

The kernel is integrated over its pixels, cut off at 40 first-null radii (more than 99% analytical encircled energy), and normalized to unit sum. Zero-padded convolution avoids wrap-around light. Cropping naturally loses light outside the frame. Camera convolution occurs in linear light, before the display tone curve and sRGB encoding.

## Display

One fixed exposure is applied throughout the loop. The original tone curve uses `v = exposure max(RGB,0)`, then `v/(1+max(v))` per pixel. The denominator is common to all channels, preserving in-gamut RGB ratios and hence CIE chromaticity; it does not normalize all pixels to the same brightness. Tiny out-of-gamut negative sRGB primaries are clipped before this display curve. This compression is a display choice and does not alter the stored physical flux calculations.

## Validation

Run `python optics.py --output optics-diagnostics.json`. Checks cover hemisphere-integrated limb normalization, flux conservation by the diluted blackbody spectrum, tone-curve chromaticity preservation, the first zero of the Bessel function, PSF normalization, impulse flux, and convolution linearity. The saved diagnostic report records their numerical residuals.
