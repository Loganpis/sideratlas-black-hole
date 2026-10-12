# @atlastro · NO.001 · THE BLACK HOLE

A silent, eight-second Instagram Reel generated from numerical light propagation around a non-spinning black hole. The large image, black background, Latin Modern lettering and STIX mathematical symbols are produced entirely by code. There is no character, star field, painted glow or invented disk texture.

The main film is a **bolometric intensity visualization**: its grayscale brightness represents radiation summed over all frequencies. The optional thermal-color film instead integrates the model's visible spectrum. These are two different observables of the same idealized disk, not interchangeable color treatments.

## Render it

Use Python 3.10 or later with NumPy 2.x and Pillow, plus an FFmpeg installation with the `libx264` encoder. FFmpeg and FFprobe must be available on your executable search path. The data and fonts are bundled; rendering needs no network connection.

```sh
python -m pip install -r requirements.txt
python render.py --preview
python render.py --mode both
python verify.py
```

The default canvas is 1080 × 1920 pixels, with 240 frames at 30 frames per second. The MP4 uses H.264, 4:2:0 color sampling and a fast-start container, without an audio track. Rendering and verification results are saved in `output/`:

- `atlastro-black-hole-reel.mp4` — main bolometric Reel.
- `atlastro-visible-thermal-reel.mp4` — visible thermal-color companion.
- `atlastro-bolometric-cover.png` — full-resolution cover for the main Reel.
- `atlastro-thermal-cover.png` and `atlastro-contact-sheet.jpg` — companion cover and comparison views.
- `render-metadata.json` — model, sampling, exposure and output settings.

Use `python render.py --help` for the actual defaults and options. A smaller preview can be generated with `--preview --width 540`; `--supersample` controls spatial sampling. `--mode bolometric` selects the grayscale image, `--mode thermal` selects visible thermal color, and `--mode both` renders both. `--output` selects another output directory. Keep the bundled `fonts/` and `data/` directories beside the Python modules.

## What the animation shows

The luminous material is one flat, opaque accretion disk. Curved light paths make it appear above and below the central darkness: those arcs are additional images of the disk, not material rising vertically from it. Relativistic orbital motion makes the approaching side brighter. Gravity and motion together shift the radiation's frequency and intensity.

The event horizon itself emits no light in this model. Its radius is **2 r_g**; the boundary between captured and escaping rays for a distant observer is at screen impact parameter **√27 r_g**. These are different quantities. A dark image pixel can also be a ray that escapes but misses the finite luminous disk. Foreground disk material can cover part of the captured-ray region, so the visible dark shape should not be identified uncritically with either the horizon or a perfectly circular shadow.

The disk is steady and axisymmetric. Its gas has the circular-orbit velocity required by the frequency-transfer calculation, but the film contains **no simulated time-dependent fluid motion**. Animation comes from a repeating viewing-angle scan:

```text
i(p) = 82.5° + 4.5° cos(2πp),  p = frame / 240.
```

Inclination is measured from the disk normal, giving 78°–87°. Each frame represents a separate static observer at infinity. Eight seconds is the presentation duration, not a physical journey or eight seconds of black-hole evolution. The endpoint at phase 1 equals phase 0 mathematically; the encoded film omits that duplicate endpoint so playback advances by one ordinary frame step through the join.

## Physical parameters and assumptions

| Quantity | Value / convention |
| --- | --- |
| Spacetime | Schwarzschild; dimensionless spin a* = 0 |
| Black-hole mass | 10⁸ solar masses |
| Distance unit | r_g = GM/c² |
| Horizon radius | 2 r_g |
| Disk inner edge | 6 r_g, the innermost stable circular orbit |
| Disk outer edge | 40 r_g |
| Accretion normalization | η Ṁ c² / L_Edd = 0.1 |
| Efficiency used in that normalization | η = 1 − √(8/9) ≈ 0.0572 |
| Eddington luminosity convention | Fully ionized hydrogen; L_Edd = 4πGM m_p c / σ_T |
| Surface | Infinitesimally thin, opaque, isotropic emission from either face |
| Disk heating | Zero-torque Page–Thorne solution for a* = 0 |
| Spectrum for thermal version | Planck blackbody, color correction f_col = 1 |
| Observer | Static at infinity for each frame |
| Image roll | 15° |
| Horizontal field of view | 32 r_g across the ray-traced image |

The luminosity ratio defines Ṁ through the binding-energy efficiency; it is **not** a claim that the finite, ray-traced disk delivers exactly 0.1 L_Edd to infinity after photon capture. The close composition crops the outer disk at the image sides.

## Calculations

`raytrace.py` integrates the exact Schwarzschild null-orbit equation in dimensionless units using fourth-order Runge–Kutta:

```text
u″(ψ) = 3u² − u,    u = r_g/r,
u(0) = 0,          u′(0) = 1/b.
```

Here b is the screen impact parameter in units of r_g, and ψ is the angle in each ray's orbital plane. A precomputed radial lookup is interpolated across detector pixels and inclinations. For each ray, the first intersection with emitting material between 6 and 40 r_g occludes subsequent intersections. Up to four disk-plane crossings are evaluated. The implementation refines sampling around the critical impact parameter √27; finite pixel resolution still limits the narrowest higher-order images. This is original code, not DNGR or AART.

With x = r/r_g, the frequency ratio of received to emitted light is

```text
g = ν_obs / ν_em = √(1 − 3/x) / (1 − λ x^(−3/2)),
λ = −α sin(i).
```

α is the detector's horizontal coordinate in the disk-aligned screen, in r_g. This expression includes gravitational and orbital frequency shifts for circular Keplerian emission and a static observer at infinity.

`emission.py` evaluates the Schwarzschild specialization of the Page–Thorne surface flux:

```text
F(x) = [3 c⁶ Ṁ / (8π G² M²)] J(x) / [x^(5/2)(x − 3)],

J(x) = √x − √6 − (√3/2) ln{
         [(√x − √3)/(√x + √3)] /
         [(√6 − √3)/(√6 + √3)] }.
```

F is the local emitted power per unit area of one face; it is zero at and below the stress-free inner edge. The local temperature satisfies F = σ_SB T_eff⁴. No luminosity pattern or hot spots are painted onto this solution.

For the **main bolometric film**, isotropic emission and relativistic intensity transfer give

```text
I_bol,obs = g⁴ F / π.
```

One fixed exposure maps **1.30 × 10¹² W m⁻² sr⁻¹** to linear display white, followed by the sRGB display transfer function. The same exposure applies to every pixel and frame; the grayscale is a presentation of all-frequency energy, not a prediction of naked-eye color. Display clipping and quantization are distinct from the calculated radiance.

For the **visible thermal-color film**, Liouville transfer gives

```text
I_ν,obs(ν) = g³ B_ν(ν/g, T_eff) = B_ν(ν, g T_eff).
```

The code integrates the Planck spectrum at T_obs = g T_eff against the CIE 1931 two-degree color-matching functions from 360–830 nm in 1 nm steps, then converts XYZ to linear sRGB. It divides all three channels by the fixed reference **1.70 × 10⁹** in CIE-weighted radiance units before display encoding. There is no second multiplication of RGB values by g³ or g⁴, and no per-pixel color normalization. The hot disk appears blue-white in this visible-band model; much of its total emission lies outside visible wavelengths. CIE data and ICC conversion details are included in `data/SOURCES.md`.

The default render samples the hero image at 1.5 times its output dimensions, then averages linear physical radiance before the nonlinear display curve. The encoded MP4 is tagged with sRGB transfer, BT.709 primaries and the BT.709 RGB-to-YCbCr matrix. Compression and social-platform re-encoding can alter the displayed pixels.

## Scientific scope

This is a numerically ray-traced, idealized relativistic thin disk. It is not a full general-relativistic magnetohydrodynamic (GRMHD) simulation, an observed black-hole image or a prediction for a named astronomical object. The equations are physically motivated, while the mass, accretion rate, disk truncation, camera scan, crop and display exposure are explicit choices for this visualization.

It omits black-hole spin, disk thickness and warps, turbulence, magnetic fields, jets, a corona, scattering, disk-atmosphere corrections, limb darkening, returning radiation, emission inside the inner edge, and feedback of photon capture on the disk temperature. In particular, an infinitely thin disk viewed nearly edge-on is a useful idealization; a disk with appreciable thickness can obscure different regions. Adding more plausible phenomena would require additional physical inputs rather than decorative effects.

## Verification

`python verify.py` checks the delivered media and the scientific modules. `python raytrace.py` and `python emission.py` also run their own diagnostics independently.

The scientific checks cover left/right symmetry of the non-spinning geometry, emission-radius bounds, Doppler handedness, weak-field bending against 4/b, and refinement of the numerical ray lookup. The emission checks compare the closed-form flux against an independently integrated Page–Thorne expression; verify its maximum near 9.55 r_g and its zero inner-edge flux; integrate Planck's law against Stefan–Boltzmann; verify spectral frequency transfer; and bound the color-lookup interpolation error. These tests support this implementation's stated model, not omitted astrophysics or every unresolved image detail.

An additional physics audit compares strong-field bending against the Schwarzschild asymptotic result in Bozza (2002), equation (50). `physics-audit.json` retains the audit results. This is a validation reference, not an alternative rendering approximation.

The video checks inspect dimensions, duration, frame count, codec and pixel format, decode the complete streams, and assess the loop join. Exact settings and measured results are retained with the outputs rather than inferred from the filename.

## Every study and standard used

The five papers below are credited on the Reel. Their roles differ: the first three inform the physical model and equations; the fourth provides a ray-traced appearance and transfer reference; the fifth supports numerical validation. The color bonus additionally uses the CIE dataset and ICC specification. Full source names and links are also in `caption.txt`.

1. **Don N. Page & Kip S. Thorne (1974)**, *Disk-Accretion onto a Black Hole. I. Time-Averaged Structure of Accretion Disk*, The Astrophysical Journal **191**, 499–506. Zero-torque relativistic disk flux, specialized here to Schwarzschild. [DOI: 10.1086/152990](https://doi.org/10.1086/152990) · [Paper](https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P).
2. **Jean-Pierre Luminet (1979)**, *Image of a spherical black hole with thin accretion disk*, Astronomy and Astrophysics **75**, 228–235. Thin-disk appearance, light bending and relativistic brightness transfer. [Paper](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L). No DOI is asserted for this archival paper.
3. **Thomas Müller & Jörg Frauendiener (2012)**, *Interactive visualization of a thin disc around a Schwarzschild black hole*, European Journal of Physics **33**, 955–963. Schwarzschild disk rendering and the radial flux expression, including §2.3, equation (11). [DOI: 10.1088/0143-0807/33/4/955](https://doi.org/10.1088/0143-0807/33/4/955) · [Paper](https://arxiv.org/abs/1206.4259). This renderer uses numerical integration rather than their analytic GPU implementation.
4. **Oliver James, Eugénie von Tunzelmann, Paul Franklin & Kip S. Thorne (2015)**, *Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar*, Classical and Quantum Gravity **32**, 065001. Visual and relativistic-transfer reference. [DOI: 10.1088/0264-9381/32/6/065001](https://doi.org/10.1088/0264-9381/32/6/065001) · [Paper](https://arxiv.org/abs/1502.03808). Its Kerr spacetime and DNGR ray-bundle implementation are not implemented here.
5. **Valerio Bozza (2002)**, *Gravitational lensing in the strong field limit*, Physical Review D **66**, 103001. Equation (50) provides a strong-deflection reference for the numerical audit. [DOI: 10.1103/PhysRevD.66.103001](https://doi.org/10.1103/PhysRevD.66.103001) · [Paper](https://arxiv.org/abs/gr-qc/0208075).
6. **CIE (2019)**, *Colour-matching functions of CIE 1931 standard colorimetric observer*, International Commission on Illumination, Vienna. Spectral color integration for the thermal version. [DOI: 10.25039/CIE.DS.xvudnb9b](https://doi.org/10.25039/CIE.DS.xvudnb9b) · [Official dataset](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer). The original data and metadata are bundled; the loader verifies the published MD5 checksum.
7. **International Color Consortium**, *How to interpret the sRGB color space (specified in IEC 61966-2-1) for ICC profiles*. Section A.7 supplies XYZ-to-linear-sRGB conversion; A.8 supplies the display transfer function. [Official specification](https://registry.color.org/rgb-registry/files/sRGB.pdf). A display standard, not an astrophysics study.

Other papers in the earlier research reading guide are background reading, not calculations silently implemented in this film. No NASA, DNGR, AART or GRMHD imagery is incorporated.

## Included assets

`render.py` composes and encodes the Reel, `raytrace.py` calculates null geodesics, `emission.py` calculates surface emission and spectral color, and `verify.py` verifies the result. `caption.txt` is ready-to-copy post text. Font attribution and the original GUST/LPPL and SIL Open Font License files are in `fonts/`; data provenance is in `data/`. Latin Modern is based on Computer Modern, the familiar LaTeX family; equations use STIX Two Math.
