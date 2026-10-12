# @atlastro — NO.001 — THE BLACK HOLE: Kerr and visible color

This version computes visible color from the emitted spectrum of an idealized accretion disk around a spinning black hole. Its light paths follow the Kerr metric. Relativistic lensing, frequency shifts, fluid-frame viewing angles, and occultation determine the image; no photon rim, halo, jets, stars, face, or surface texture is painted onto it.

The black hole itself does not emit the visible disk light. “True color” here means the color predicted for **this specified disk and viewing model**, integrated through standard human color-matching functions and adapted to a finite-brightness screen. It does not mean every black hole has these colors, or that this is a photograph of a particular object.

## Deliverables

- `output-kerr/atlastro-kerr-truecolor-reel.mp4`: main Reel, including the subtle diffraction response of the specified hypothetical camera.
- `output-kerr/atlastro-kerr-ideal-optics-reel.mp4`: companion with ideal angular resolution and no camera point-spread function. The same astrophysical model and spectral color pipeline are used.
- `caption-kerr.txt`: Instagram caption including all studies, data, and display references used.
- `render_kerr.py`: composition and rendering entry point; it reuses font and encoding helpers from the original `render.py`.

The films are eight seconds, 240 frames, 30 fps, 1080 × 1920 pixels, H.264 with yuv420p and no audio. The periodic camera-angle scan makes a loop; it is not eight seconds of simulated gas evolution.

## Defined physical system

| Quantity | Value or assumption |
| --- | --- |
| Spacetime | Kerr; positive-spin, prograde equatorial disk |
| Mass | **10¹⁰ solar masses** |
| Dimensionless spin | **a\* = 0.9** |
| Accretion normalization | **ηṀc²/L_Edd = 0.05** |
| Binding-energy efficiency | η = 1 − E_ISCO = 0.155752992 |
| Length unit | r_g = GM/c² |
| Horizon coordinate radius | r_+ = 1.435889894 r_g |
| Disk inner edge | Stress-free ISCO at 2.320883042 r_g |
| Disk outer edge | 500 r_g; opaque at every emitting intersection |
| Geometry | Steady, axisymmetric, razor-thin disk |
| View inclination | Smooth cosine scan, 80°–86° from the spin axis |
| Image roll | 14° |
| Horizontal field | 26 r_g; a close crop of the larger disk |
| Astrophysical observer | Static at infinity |
| Display | One fixed exposure and RGB-ratio-preserving compression |

**The physical parameters changed from the earlier film:** that example used a non-spinning, 10⁸-solar-mass hole, ηṀc²/L_Edd = 0.1, and an outer disk radius of 40 r_g. This example uses a larger, cooler system and a much more extended disk. The color change follows those explicit assumptions, together with the improved transfer and display pipeline; it is not a claim that the former mass intrinsically has the new spectrum.

For these new parameters, the local effective temperature peaks near **32,438 K** at **3.4419 r_g**. The atmosphere approximation gives a maximum color-hardening factor of approximately **1.0662**. Gravitational and orbital frequency shifts then change the spectrum seen from each image location. The thermal spectrum can look blue-white, white, or warmer in different regions without assigning arbitrary colors to “approaching” and “receding” material.

## Light paths and the dark center

`kerr_raytrace.py` and its original `kerr_raytrace.cpp` kernel integrate the separated Kerr null-geodesic potentials. They use fourth-order Runge–Kutta integration, regular second-order turning-point equations, and cubic-Hermite disk-intersection positions. The source equations are those collected by Gralla & Lupsasca (2020). Screen coordinates determine the conserved angular momentum λ and Carter parameter q; no spherical shadow mask is pasted over the image.

The ray stops at its **first intersection with opaque disk material**. A crossing through the empty inner hole continues onward and can reveal a higher-order disk image. The extended foreground disk blocks light from some more distant images, even when those images would be visually attractive. Consequently an unobstructed, uniformly bright ring is not guaranteed. The dark image is also not a direct geometric outline of the event horizon.

With G = M = c = 1, the orbital frequency and frequency transfer are

```text
Ω = 1 / (r^(3/2) + a)
g = ν_obs/ν_em = 1 / [u^t (1 − Ω λ)]
μ_em = g sqrt(q) / r
```

The spin-dependent circular orbits and ISCO follow Bardeen, Press & Teukolsky (1972). Frame dragging affects the computed trajectories and image asymmetry. The stationary axisymmetric emission does not require azimuth or travel time, so those coordinates are not integrated; there is no time-delay animation or moving hotspot hidden in the calculation.

## Disk energy and visible spectrum

`kerr_emission.py` evaluates the stress-free Page–Thorne disk flux from one disk face. With r dimensionless in r_g and E and L the circular-orbit integrals,

```text
F(r) = [Ṁ c⁶/(4πG²M²)] × [−Ω,r / {r(E − ΩL)²}]
       × ∫[r_ISCO to r] (E − ΩL) L,r dr

T_eff = (F/σ_SB)^(1/4)
```

The orbital quantities and derivatives are analytic; a dense radial grid evaluates the final integral. The accretion-rate convention uses η = 1 − E_ISCO and the hydrogen-electron-scattering Eddington luminosity. It is not a correction for the fraction of emitted light later captured by the hole.

`optics.py` adds explicitly approximate atmosphere effects. It uses the temperature-dependent color-hardening branches of Done et al. (2012), joined at their actual intersection near 101,749 K to avoid an artificial step. It uses the flux-normalized limb law from Li et al. (2005), Υ(μ) = 1/2 + 3μ/4. This is an approximate angular/spectral prescription, not a solved AGN atmosphere or the full KERRBB implementation.

The observed specific intensity is

```text
I_ν,obs = Υ(μ_em) f_col⁻⁴ B_ν(ν, g f_col T_eff)
```

The f_col⁻⁴ dilution preserves local bolometric flux. The relativistic transfer is already included in the shifted Planck spectrum: the renderer does **not** multiply the resulting RGB by another g³ or g⁴. Bolometric boosting is not substituted for visible-band brightness.

`emission.ThermalSpectrum` integrates the spectrum over the official **CIE 1931 2° observer, 360–830 nm, 1 nm spacing**, and converts XYZ to linear sRGB using the ICC specification. The bundled dataset is checksum-verified. A common scalar tone curve preserves in-gamut RGB ratios, while reducing radiance to screen range. Exposure remains fixed through the loop. These display choices do not reconstruct the scene's full physical dynamic range or a particular observer's state of visual adaptation.

## Camera diffraction is separate from black-hole lensing

The main Reel applies the intensity point-spread function of a hypothetical clear circular **20 m aperture**, with **550 nm** as the reference wavelength and angular scale corresponding to a distance of **10⁶ r_g**. Its first diffraction null is approximately **1.39 pixels** at a 1080-pixel horizontal field of 26 r_g. The Airy kernel is pixel-integrated, normalized, and convolved in linear light before the display transform.

The same 550 nm kernel is applied to all RGB channels: this is an **achromatic reference-wavelength approximation**, not a full polychromatic telescope model. There is no glare threshold, Gaussian bloom, or hand-painted glow. Faint diffraction wings are a camera response, not radiation emerging from the horizon or an extra general-relativistic halo. The ideal-optics companion omits this response so its contribution can be inspected directly.

## Scope and validation

This is a more detailed relativistic thin-disk calculation, **not a GRMHD evolution or a universal prediction of an accreting black hole's appearance**. It omits finite disk thickness, turbulence, magnetic fields, corona/Comptonization, jets, emission from plunging material, returning-radiation heating, and a full atmosphere calculation. Finite thickness could matter especially for the chosen high inclinations. An exact mathematical photon subring can be narrower than a rendered pixel; no guarantee is made that every order is resolved at Instagram resolution.

The included diagnostics exercise numerical and physical invariants rather than just checking that an image exists:

- `kerr-emission-diagnostics.json`: Kerr a=0 agrees with the independent Schwarzschild closed form to 4.73 × 10⁻⁸ in the sampled radii; independent Kerr quadrature agrees to 3.98 × 10⁻⁸; doubling the grid changes sampled flux by at most 3.72 × 10⁻⁸. The large-radius limit and zero ISCO flux also pass. Those diagnostics use their documented baseline mass; the dimensionless flux checks carry over to the film's mass scaling.
- `kerr-raytrace-diagnostics.json`: comparison with the independent Schwarzschild tracer, step refinement, geodesic-potential constraints, fluid-frame angle bounds, and capture classification from independent radial-polynomial roots. Consult the report for exact resolutions and sampled viewing angles.
- `optics.py --output optics-diagnostics.json`: limb-law normalization, color-hardening energy conservation, tone-curve chromaticity, Airy null, kernel normalization, and convolution linearity.
- Original diagnostics and `PHYSICS-AUDIT.md` remain as provenance for the Schwarzschild baseline and its strong-deflection checks; their results alone are not presented as validation of the Kerr calculation.

## Reproduce

Use Python 3 with NumPy and Pillow, a local C++17 compiler (`clang++` or `g++`), and FFmpeg with libx264. The supplied C++ kernel compiles locally on first use. The calculation does not fetch a remote model, simulation, font, or spectrum.

```sh
python3 -m pip install -r requirements.txt
python3 render_kerr.py --help
python3 render_kerr.py
python3 verify_kerr.py
```

Useful independent checks:

```sh
python3 kerr_emission.py --output kerr-emission-diagnostics.json
python3 kerr_raytrace.py --output kerr-raytrace-diagnostics.json
python3 optics.py --output optics-diagnostics.json
```

Keep `render.py`, `raytrace.py`, `emission.py`, the new Kerr modules, `optics.py`, and the bundled `fonts/` and `data/` directories together. Latin Modern and STIX Two Math provide the scientific typography. Font and data licenses/provenance accompany their assets.

`verify_kerr.py` checks the chosen physical parameters, both full video streams,
their color tags and loop joins, all 121 independently traced viewing angles,
and exact uncompressed periodicity at preview resolution. Its final report is
`output-kerr/verification.json`. The source archive includes the original C++
source; its small platform-specific compiled library is rebuilt locally.

## Every study, dataset, and display reference used

The first five entries provide the new model equations. The next four are retained as explicitly identified appearance references and independent-baseline validation. They are credited rather than silently represented as additional simulation engines.

1. **Gralla, S. E. & Lupsasca, A. (2020).** *The Null Geodesics of the Kerr Exterior.* Physical Review D **101**, 044032. Separated Kerr potentials and geodesic equations; the local implementation integrates them numerically rather than using the paper's elliptic-function solution. [DOI](https://doi.org/10.1103/PhysRevD.101.044032) · [arXiv](https://arxiv.org/abs/1910.12881).
2. **Page, D. N. & Thorne, K. S. (1974).** *Disk-Accretion onto a Black Hole. I. Time-Averaged Structure of Accretion Disk.* The Astrophysical Journal **191**, 499–506. Local relativistic disk flux, equations 11b, 12 and 15. [DOI](https://doi.org/10.1086/152990) · [Paper](https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P).
3. **Bardeen, J. M., Press, W. H. & Teukolsky, S. A. (1972).** *Rotating Black Holes: Locally Nonrotating Frames, Energy Extraction, and Scalar Synchrotron Radiation.* The Astrophysical Journal **178**, 347–370. Circular orbit and ISCO expressions. [DOI](https://doi.org/10.1086/151796) · [Paper](https://articles.adsabs.harvard.edu/pdf/1972ApJ...178..347B).
4. **Li, L.-X., Zimmerman, E. R., Narayan, R. & McClintock, J. E. (2005).** *Multi-Temperature Blackbody Spectrum of a Thin Accretion Disk around a Kerr Black Hole: Model Computations and Comparison with Observations.* The Astrophysical Journal Supplement Series **157**, 335–370. Angular limb-darkening and diluted thermal transfer conventions. The full KERRBB code and its returning-radiation solution are not used. [DOI](https://doi.org/10.1086/428089) · [arXiv](https://arxiv.org/abs/astro-ph/0411583).
5. **Done, C., Davis, S. W., Jin, C., Blaes, O. & Ward, M. (2012).** *Intrinsic disc emission and the soft X-ray excess in active galactic nuclei.* Monthly Notices of the Royal Astronomical Society **420**, 1848–1860. Approximate temperature-dependent color-hardening branches, equations 1–2; their corona/soft-excess model is not implemented. [DOI](https://doi.org/10.1111/j.1365-2966.2011.19779.x) · [arXiv](https://arxiv.org/abs/1107.5429).
6. **Luminet, J.-P. (1979).** *Image of a spherical black hole with thin accretion disk.* Astronomy and Astrophysics **75**, 228–235. Original thin-disk appearance and transfer reference for the Schwarzschild baseline. [Paper](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L).
7. **Müller, T. & Frauendiener, J. (2012).** *Interactive visualization of a thin disc around a Schwarzschild black hole.* European Journal of Physics **33**, 955–963. Baseline disk-rendering equations and Schwarzschild comparison. [DOI](https://doi.org/10.1088/0143-0807/33/4/955) · [arXiv](https://arxiv.org/abs/1206.4259).
8. **James, O., von Tunzelmann, E., Franklin, P. & Thorne, K. S. (2015).** *Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar.* Classical and Quantum Gravity **32**, 065001. Appearance and relativistic-transfer reference; neither DNGR code nor its imagery is used. [DOI](https://doi.org/10.1088/0264-9381/32/6/065001) · [arXiv](https://arxiv.org/abs/1502.03808).
9. **Bozza, V. (2002).** *Gravitational lensing in the strong field limit.* Physical Review D **66**, 103001. Strong-deflection validation of the independent Schwarzschild baseline, equation 50. [DOI](https://doi.org/10.1103/PhysRevD.66.103001) · [arXiv](https://arxiv.org/abs/gr-qc/0208075).
10. **Airy, G. B. (1835).** *On the diffraction of an object-glass with circular aperture.* Transactions of the Cambridge Philosophical Society **5**, 283–291. Circular-pupil intensity response. The standard formula is also documented by the [Astropy AiryDisk2D reference](https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.AiryDisk2D.html); this implementation does not depend on Astropy.
11. **CIE (2019).** *Colour-matching functions of CIE 1931 standard colorimetric observer.* International Commission on Illumination. Spectral color integration. [Dataset DOI](https://doi.org/10.25039/CIE.DS.xvudnb9b) · [Official dataset](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer).
12. **International Color Consortium.** *How to interpret the sRGB color space (specified in IEC 61966-2-1) for ICC profiles.* XYZ-to-linear-sRGB matrix and display transfer. [Official specification](https://registry.color.org/rgb-registry/files/sRGB.pdf).

No NASA, DNGR, AART, or GRMHD frames are incorporated. Additional papers in the earlier reading guide remain background reading, not unreported calculations in this film.
