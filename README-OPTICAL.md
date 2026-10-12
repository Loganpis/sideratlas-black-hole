# @atlastro — NO.001 — THE BLACK HOLE: spectral optics

This Reel calculates the visible spectrum of a relativistic accretion disk around a spinning black hole, then passes that spectrum through a specified camera. It retains the Kerr light paths and thermal disk model of the previous version. The optical upgrade is **wavelength-dependent diffraction before color integration**, with a closer, taller view of the lensed disk.

“True color” means the color predicted for **this particular physical model**, converted through standard human color-matching functions and adapted to screen brightness. The black hole does not emit the visible disk light. Its disk's color depends on its mass, accretion rate, atmosphere and viewing geometry; these are not universal black-hole colors or a photograph of a measured object.

## Outputs

- `output-optical/atlastro-truecolor-optics-reel.mp4`: the main Reel, with wavelength-resolved camera diffraction.
- `output-optical/atlastro-truecolor-ideal-reel.mp4`: the same physical model, spectral integration and display settings with the camera diffraction removed.
- `caption-optical.txt`: a ready-to-use caption crediting every study, dataset and display reference used.
- `render_optical.py`: rendering entry point.
- `spectral_optics.py`: wavelength integration and the camera response.
- `verify_optical.py`: numerical and delivery checks.

Both films use a black background, Latin Modern and STIX Two Math typography, @atlastro, NO.001, THE BLACK HOLE and GRAVITY TYPE. Their format is 1080 × 1920 pixels, eight seconds, 240 frames, 30 fps, H.264/yuv420p, without audio. The loop scans the viewing angle of a steady model; it is not eight seconds of evolving accretion gas.

## Physical model and framing

| Quantity | Value or assumption |
| --- | --- |
| Spacetime | Kerr; dimensionless spin a\* = 0.9 |
| Mass | 10¹⁰ solar masses |
| Accretion normalization | ηṀc²/L_Edd = 0.05 |
| Efficiency | η = 1 − E_ISCO = 0.155752992 |
| Radius unit | r_g = GM/c² |
| Horizon coordinate radius | r_+ = 1.435889894 r_g |
| Emitting disk | Prograde, steady, axisymmetric, opaque, razor-thin |
| Disk edges | Stress-free ISCO at 2.320883042 r_g; outer edge 500 r_g |
| Ray-tracing observer | Static at infinity |
| Viewing-angle scan | i = 86.5° + 1.5° cos(2πt/8 s), hence 85°–88° |
| Image roll | 14° |
| Horizontal field | 24 r_g; a close crop of the larger disk |
| Hero region | 1080 × 1620 pixels, starting at y = 30 |

These mass, accretion and spin values are unchanged from the previous Kerr version. The local effective temperature peaks near 32,438 K at 3.4419 r_g, and the approximate atmosphere reaches a color-hardening factor of about 1.0662. Their displayed appearance follows the shifted spectrum and the specified exposure, not an assigned orange or blue palette.

The extreme inclination is a view into the geometry of this thin-disk model. A disk with finite thickness could obscure more of its inner emission at these angles; the film does not model that additional structure.

## Light paths and emission

`kerr_raytrace.py` and its original C++ kernel integrate the separated Kerr null-geodesic potentials with fourth-order Runge–Kutta integration and interpolated disk crossings. The source equations are collected by Gralla & Lupsasca (2020); the circular orbits and ISCO follow Bardeen, Press & Teukolsky (1972).

A ray stops at its first crossing of opaque emitting material. A crossing inside the empty inner disk edge can continue to another image of the disk. Foreground disk material therefore hides some farther images. Frame dragging, lensing and orbital motion affect the calculated image; no bright ring is drawn over it. The dark image is not the horizon's direct outline.

In geometric units G = M = c = 1:

```text
Ω = 1 / (r^(3/2) + a)
g = ν_obs/ν_em = 1 / [u^t (1 − Ωλ)]
μ_em = g sqrt(q) / r
```

Here λ and q are the ray's conserved angular-momentum and Carter parameters, and μ_em is its emission-angle cosine in the orbiting gas frame. The model is stationary and axisymmetric, so azimuthal phase and travel time are not needed for its intensity calculation. There is no hidden gas-orbit or time-delay animation.

`kerr_emission.py` computes the zero-torque Page–Thorne (1974) disk flux from one disk face:

```text
F(r) = [Ṁc⁶/(4πG²M²)] × [−Ω,r / {r(E − ΩL)²}]
       × ∫[r_ISCO to r] (E − ΩL)L,r dr
T_eff = (F/σ_SB)^(1/4)
```

The orbit quantities and their derivatives are analytic; a dense radial quadrature evaluates the integral. Accretion normalization uses the hydrogen-electron-scattering Eddington luminosity and η = 1 − E_ISCO. It does not include a further correction for emitted light later captured by the hole.

The atmosphere is explicitly approximate: the temperature-dependent color-hardening branches follow Done et al. (2012), and the flux-normalized limb law follows Li et al. (2005): Υ(μ) = 1/2 + 3μ/4. The observed spectrum is

```text
I_ν,obs = Υ(μ_em) f_col⁻⁴ B_ν(ν, g f_col T_eff).
```

The dilution f_col⁻⁴ preserves the local bolometric thermal flux. Frequency transfer is already contained in the shifted Planck spectrum. Neither another g³/g⁴ multiplier nor bolometric intensity is applied to the resulting visible RGB. See `OPTICS-NOTES.md` for the atmosphere prescription; its older camera and tone-curve sections describe the previous version, superseded here.

## Wavelength-resolved optical effects

The new camera samples **48 wavelengths from 360 to 830 nm at 10 nm intervals**, using trapezoidal integration with half weights at both endpoints. The official checksum-verified CIE 1931 2° observer data have 1 nm spacing; these matching functions are sampled at the actual integration wavelengths.

At each wavelength λ, the pipeline is:

1. Calculate the observed Planck radiance at T_obs = g f_col T_eff and multiply by Υ/f_col⁴.
2. Average the supersampled **radiance**, not its temperature, into detector pixels.
3. Convolve that monochromatic image with its own wavelength-dependent Airy intensity point-spread function.
4. Integrate the convolved spectrum against the CIE color-matching functions, then convert XYZ to linear sRGB.

The hypothetical camera has a **10 m unobstructed circular aperture** at a distance of **10⁶ r_g**. That distance sets the angular image scale; the astronomical light-path calculation uses the infinity-observer approximation. Its circular-pupil response is

```text
P_λ(θ) ∝ [2 J₁(πDθ/λ)/(πDθ/λ)]²
θ_first-null = 1.21967 λ/D.
```

At the delivered image scale, the first null is approximately 1.98 pixels at 360 nm, 3.02 pixels at 550 nm and 4.56 pixels at 830 nm. Red light is therefore diffracted more broadly than blue light. This is a spectral calculation, not an RGB-channel blur with an arbitrarily assigned wavelength.

Each pixel-integrated kernel is cut off at 24 first-null radii and normalized to unit sum; truncating the remote wings is a numerical approximation. A traced overscan margin extends beyond the delivered crop by the largest kernel radius, retaining the external disk light that can diffract into the image. Zero-padded convolution avoids wrap-around.

The default export accelerates this same linear spectral calculation with a rank-18 basis of normalized Planck spectra over approximately 24.764–100,000 K. The lower boundary is the direct Planck implementation's existing numerical-zero threshold, not a new cutoff for faint light. It restores the absolute Planck amplitude at each source subpixel before averaging; the internal normalization never changes the disk brightness. The wavelength-specific Airy kernels and CIE weights are precombined with the basis. Against the direct wavelength sum across the tested scenes and Kerr views, the largest RGB difference is about 1.23 × 10⁻⁷ of the image peak. The normalized spectral-shape error reaches 1.16 × 10⁻⁴ near the numerical-zero discontinuity at about 25 K; absolute radiance there is negligible, and the separate rendered-radiance checks enforce a 2 × 10⁻⁶ error bound. Temperatures outside the basis domain, and entirely cold scenes, use the direct method. `--direct-spectrum` keeps that reference calculation available. `spectral-optics-fast-diagnostics.json` records the tested range and numerical errors.

Diffraction is an **instrument response**, distinct from the gravitational lensing around the black hole. The ideal-optics companion removes only this response. Neither film adds a glare threshold, atmospheric seeing, lens flare, decorative photon rim, invented surface texture or an emitting halo in the central dark region.

## Screen brightness and color

Both companions share one fixed exposure and one display transform. With linear sRGB vector RGB, the default exposure E = 1/(1.8 × 10⁸) and m = max(v):

```text
v = E × max(RGB, 0)
display_linear = v × m / (1 + m²)
```

Standard sRGB encoding follows. This common scalar curve preserves in-gamut linear RGB ratios, and therefore chromaticity, while compressing brightness and reducing faint emission on screen. It is a display-contrast choice, not a change in emitted energy. Negative out-of-gamut primaries are clipped before the curve. The camera's wavelength-dependent mixing can physically change the local spectrum; the subsequent display curve does not separately recolor its channels.

The Reel cannot reproduce the actual luminous intensity or dynamic range of an accretion disk. It does not simulate an observer's adaptation, a measured telescope, or a calibrated photograph of a named black hole.

## Checks and reproduction

Install Python 3 with NumPy and Pillow, a C++17 compiler (`clang++` or `g++`), and FFmpeg with libx264. The included original C++ source compiles locally on first use. No remote simulation or image needs to be downloaded.

```sh
python3 -m pip install -r requirements.txt
python3 render_optical.py --help
python3 render_optical.py
python3 verify_optical.py
```

Rendering saves each independently traced view, its checks and image hashes in `output-optical/.render-cache/`. Running the same command after an interruption restores completed views. Checkpoints are keyed by the rendering parameters and hashes of the source, fonts and color data; different calculations cannot silently reuse them. The final MP4 files are replaced only after encoding completes.

The new renderer uses a scaled integration step of 0.018, compared against 0.006 through the selected viewing range. Every production view must satisfy a geodesic-potential residual below 10⁻⁵; a failed view is retraced at 0.006 and rejected if it still fails. `output-optical/step-speed-study.json` documents the comparison, including unchanged masks and displayed pixels in its sampled views. This is a checked numerical speed improvement within the same Kerr equations.

Keep the source modules, `fonts/` and `data/` directories together. Font licenses and dataset provenance accompany the assets. The original Schwarzschild and Kerr renderers remain separate; this version does not overwrite their documented outputs.

`spectral-optics-diagnostics.json` records the wavelength-grid, chromatic diffraction and numerical checks. Across the tested 2,000–100,000 K blackbody range, the 10 nm quadrature differs from the 1 nm reference by at most 0.1092% of the peak RGB channel, 0.001141% in CIE Y and 0.000139 in CIE x or y. These are the tested color-integration errors, not a claim about total astrophysical-model accuracy. Other checks cover wavelength scaling of the Airy first null, kernel quadrature, impulse flux, convolution linearity, finite-precision Fourier kernels and supersampled radiance averaging.

The existing Kerr diagnostics separately test disk-flux quadrature and grid refinement, reduction to the independent Schwarzschild model, light-path refinement, geodesic-potential constraints and capture classification. `verify_optical.py` checks the new outputs and loop. Consult its saved report for the actual delivery checks and measured results.

This is an idealized relativistic thin-disk calculation. It omits finite disk thickness, turbulence, magnetic fields, jets, a corona, emission inside the ISCO, returning-radiation heating, polarization and a full atmosphere solution. It is not a GRMHD evolution. Narrow higher-order images may remain unresolved at Reel resolution.

## Every study, dataset and display reference used

The first five entries provide the current astrophysical model. Entries 6–9 supply explicitly credited appearance references and independent-baseline validation; they are not additional simulation engines. Entries 10–12 provide the camera and color calculations.

1. **Gralla, S. E. & Lupsasca, A. (2020).** *The Null Geodesics of the Kerr Exterior.* Physical Review D **101**, 044032. Kerr potentials and geodesic equations; integrated numerically here rather than through the paper's elliptic-function solution. [DOI](https://doi.org/10.1103/PhysRevD.101.044032) · [arXiv](https://arxiv.org/abs/1910.12881).
2. **Page, D. N. & Thorne, K. S. (1974).** *Disk-Accretion onto a Black Hole. I. Time-Averaged Structure of Accretion Disk.* The Astrophysical Journal **191**, 499–506. Relativistic disk flux, equations 11b, 12 and 15. [DOI](https://doi.org/10.1086/152990) · [Paper](https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P).
3. **Bardeen, J. M., Press, W. H. & Teukolsky, S. A. (1972).** *Rotating Black Holes: Locally Nonrotating Frames, Energy Extraction, and Scalar Synchrotron Radiation.* The Astrophysical Journal **178**, 347–370. Circular orbit and ISCO expressions. [DOI](https://doi.org/10.1086/151796) · [Paper](https://articles.adsabs.harvard.edu/pdf/1972ApJ...178..347B).
4. **Li, L.-X., Zimmerman, E. R., Narayan, R. & McClintock, J. E. (2005).** *Multi-Temperature Blackbody Spectrum of a Thin Accretion Disk around a Kerr Black Hole: Model Computations and Comparison with Observations.* The Astrophysical Journal Supplement Series **157**, 335–370. Limb darkening and diluted thermal transfer; the full KERRBB code and returning-radiation solution are not used. [DOI](https://doi.org/10.1086/428089) · [arXiv](https://arxiv.org/abs/astro-ph/0411583).
5. **Done, C., Davis, S. W., Jin, C., Blaes, O. & Ward, M. (2012).** *Intrinsic disc emission and the soft X-ray excess in active galactic nuclei.* Monthly Notices of the Royal Astronomical Society **420**, 1848–1860. Approximate spectral-hardening branches, equations 1–2; the corona/soft-excess model is not implemented. [DOI](https://doi.org/10.1111/j.1365-2966.2011.19779.x) · [arXiv](https://arxiv.org/abs/1107.5429).
6. **Luminet, J.-P. (1979).** *Image of a spherical black hole with thin accretion disk.* Astronomy and Astrophysics **75**, 228–235. Thin-disk appearance and transfer reference for the Schwarzschild baseline. [Paper](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L).
7. **Müller, T. & Frauendiener, J. (2012).** *Interactive visualization of a thin disc around a Schwarzschild black hole.* European Journal of Physics **33**, 955–963. Baseline rendering equations and Schwarzschild comparison. [DOI](https://doi.org/10.1088/0143-0807/33/4/955) · [arXiv](https://arxiv.org/abs/1206.4259).
8. **James, O., von Tunzelmann, E., Franklin, P. & Thorne, K. S. (2015).** *Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar.* Classical and Quantum Gravity **32**, 065001. Appearance and transfer reference; neither DNGR code nor imagery is used. [DOI](https://doi.org/10.1088/0264-9381/32/6/065001) · [arXiv](https://arxiv.org/abs/1502.03808).
9. **Bozza, V. (2002).** *Gravitational lensing in the strong field limit.* Physical Review D **66**, 103001. Strong-deflection validation of the independent Schwarzschild baseline, equation 50. [DOI](https://doi.org/10.1103/PhysRevD.66.103001) · [arXiv](https://arxiv.org/abs/gr-qc/0208075).
10. **Airy, G. B. (1835).** *On the diffraction of an object-glass with circular aperture.* Transactions of the Cambridge Philosophical Society **5**, 283–291. Circular-pupil diffraction. The standard formula is also documented in the [Astropy AiryDisk2D reference](https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.AiryDisk2D.html); this code does not depend on Astropy.
11. **CIE (2019).** *Colour-matching functions of CIE 1931 standard colorimetric observer.* International Commission on Illumination. Visible spectral color integration. [Dataset DOI](https://doi.org/10.25039/CIE.DS.xvudnb9b) · [Official dataset](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer).
12. **International Color Consortium.** *How to interpret the sRGB color space (specified in IEC 61966-2-1) for ICC profiles.* XYZ-to-linear-sRGB matrix and display encoding. [Official specification](https://registry.color.org/rgb-registry/files/sRGB.pdf).

The contrast curve and layout are original display choices. No NASA, DNGR, AART or GRMHD frames are incorporated. Other papers in the earlier reading guide remain background reading rather than unreported calculations in this film.
