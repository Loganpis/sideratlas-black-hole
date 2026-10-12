# @sideratlas — original Schwarzschild disk Reel

This release replaces the previous smooth Kerr/thermal-camera scene with an **original, NASA-inspired Schwarzschild visualization**: a full viewpoint revolution above and below a thin disk, approaching-side brightening, gravitationally bent disk images, narrow higher-order images, and moving sheared brightness structures. No NASA frames, texture, simulation dataset or software are incorporated.

NASA's [SVS 13326 visualization](https://svs.gsfc.nasa.gov/13326/) documents a non-rotating black hole, evolving disk structure, and a 360° viewpoint revolution in a plane perpendicular to the disk. Its public description does not specify the complete emissivity, color, or fluid model. This recreation follows those documented features; its timing and chosen parameters do not reproduce an unpublished NASA setup. It is not a full magnetized-fluid simulation and has not been calibrated against a measured black hole.

## Deliverables

- `output-sideratlas/sideratlas-black-hole-reel.mp4`: orange illustrative brightness palette.
- `output-sideratlas/sideratlas-black-hole-intensity.mp4`: monochrome view of the same calculated bolometric radiance.
- `output-sideratlas/steady-reference.png`: same geometry and mean emission with the brightness structures switched off.
- `output-sideratlas/sideratlas-black-hole-color-cover.png`: cover image.
- `output-sideratlas/verification.json`: checks of the finished files and numerical model.
- `caption-schwarzschild.txt`: Instagram caption, assumptions and all source credits.
- `INSTAGRAM-EXPORT.md`: export choices and Instagram upload guidance.
- `sideratlas-black-hole-source.zip`: current videos, code, documentation and scientific credits.

The videos are 1080 × 1920, 30 fps, eight seconds, H.264/yuv420p, no audio. The design uses @sideratlas, NO.001, THE BLACK HOLE and Latin Modern/STIX mathematics fonts. The header sits farther below the top edge to leave more room for Instagram's interface; interface placement can vary by device and viewing mode. A single information line combines “Original relativistic model”, “Schwarzschild / a* = 0” and the color description. The four explanatory lines have no introductory label or sentence-ending periods:

```text
Bends light into multiple images of its disk
Approaching gas appears brighter
Gravity redshifts light and delays its arrival
Inner orbits shear bright structures faster
```

The former type badge and bottom model-summary footer have been removed from the artwork. Their scientific context remains documented here and in the caption. Earlier releases remain available in their own folders.

## What is calculated, and what is chosen

| Component | Calculation or explicit choice |
| --- | --- |
| Spacetime | Schwarzschild: spin zero, as documented by NASA |
| Observer | Sequence of static observers at infinity; directed viewing angle θ = 83° − 360° phase, completing one revolution perpendicular to the disk |
| Disk | Opaque, infinitesimally thin, equatorial; circular Keplerian motion; identical modeled emission on both faces |
| Emitting edges | Inner edge 6 r_g, stress-free; outer edge 28 r_g, a chosen finite cutoff |
| Image | Fixed 72 r_g horizontal span; continuous screen basis throughout the revolution, with no added roll; same physical scale on both axes; the full face-on disk fits the view |
| Mean emitted brightness | Relativistic zero-torque Page–Thorne radial flux; normalized to peak one |
| Transfer | Bolometric intensity multiplied by g⁴; first opaque disk intersection wins |
| Photon propagation | Numerically integrated Schwarzschild null geodesics and relative coordinate travel times |
| Variable brightness | Original positive, sheared finite-lived patches; explicitly kinematic, not simulated magnetic turbulence |
| Color | Original orange intensity palette; **not a temperature, measured spectrum, or true visible color** |
| Time | 800 t_g represented in 8 screen seconds; periodic source births and one constant-speed viewpoint revolution make the loop |
| Mass and accretion rate | No numerical values assigned; their overall flux scale cancels in relative brightness |
| Camera response | Ideal spatial sampling, without the previous arbitrary aperture/diffraction model |

Here r_g = GM/c² and t_g = GM/c³. The displayed eight seconds are not a claim about a real black hole's elapsed physical time. A chosen mass would be necessary to assign seconds.

The finite outer edge, viewpoint path, exposure and brightness-pattern statistics remain assumptions. Using a normalized bolometric image avoids having to choose a physical temperature, atmosphere or visible spectrum merely to obtain an orange appearance. It does not eliminate all astrophysical uncertainty.

## Relativity and photon arrival times

In G = M = c = 1 units the orbit solver integrates

```text
u = 1/r,     u'' = 3u² − u
Ω(r) = r^(−3/2)
g = sqrt(1 − 3/r) / [1 − Ω λ],     λ = −α sin(θ)
I_obs = g⁴ I_em
```

The observer-plane impact parameter at the critical photon orbit is sqrt(27) r_g; the horizon radius is 2 r_g. A disk image or black region is not simply a drawing of the horizon. Foreground disk material occludes subsequent intersections, including some higher-order images. The implementation follows up to four candidate equatorial intersections; still narrower images can remain unresolved.

`nasa_raytrace.py` extends the previously checked `raytrace.py` mapping with a finite relative travel-time integral. Let ψ increase backward from the observer, b be impact parameter, v = du/dψ, and r* = r + 2 ln(r/2 − 1). If ℓ is elapsed backward coordinate travel time, define q = ℓ + r*(r) − r*(observer). Then

```text
dq/dψ = (1/b − v) / [u²(1 − 2u)] = 1 / (1/b + v)
delay = q − r*(r_hit)
t_em = t_obs − delay
```

The incoming branch uses the cancellation-safe reciprocal; the outgoing branch uses the first expression. One common infinite observer offset is removed. Negative relative delays are legitimate. No per-image or per-frame shift is introduced. Consequently the direct and bent images can show the emitting surface at different source times.

## Full viewpoint revolution

The eight-second viewpoint path completes a full revolution in a plane perpendicular to the disk:

```text
phase = t / 8 s
θ(t) = 83° − 360° phase
```

Here θ is a directed angle, not an inclination restricted to 0°–180°. The conventional inclination is arccos(cos θ). The observer direction and screen basis rotate continuously; they do not reset at either pole. It starts near the disk plane at θ = 83°, passes the north pole at about 1.84 seconds, crosses the opposite side of the disk plane at 3.84 seconds, passes the south pole at 5.84 seconds, and crosses the plane again at 7.84 seconds. At eight seconds it returns to the initial view. The endpoint is not duplicated in the 240-frame export.

This follows the type of full revolution described for NASA's movie, with our own constant angular speed and eight-second duration. It does not copy NASA's pauses or timing. The field of view remains fixed, so the apparent shape change comes from viewing geometry rather than a changing zoom. For this layout, the 1080 × 900 radiance image is reduced uniformly to 864 × 720 and centered above the explanatory text. The whole disk remains in view at the poles. This changes only the display framing; the ray tracing, disk dynamics and viewpoint motion are unchanged.

This is a scan through virtual static observers, not the trajectory of a physically moving camera; camera-velocity aberration is not modeled. Each viewpoint gets its own recomputed disk intersections, relativistic transfer and relative photon travel times, with the same observer-time convention. Changing lensed arcs and foreground occultation follow those intersections; no visibility fade or mirrored footage is added. Apparent disk rotation can reverse when the disk is viewed from below. The source continues evolving forward throughout the revolution, and the black hole remains non-rotating.

Both faces use the same modeled brightness at a given disk position and emission time. This is an explicit idealization, not a calculated disk atmosphere. Exactly edge-on views are especially sensitive to the infinitesimally thin disk assumption and finite pixel sampling; a disk with thickness would obscure light differently.

At the exact disk-plane crossing, the direct image of this infinitesimally thin disk has zero height, while lensed far-side rings can persist; the ray calculation therefore does not force the whole disk image to disappear.

## Moving structure: the extra assumption

`nasa_emission.py` keeps the exact azimuthal mean at the Page–Thorne value and adds a specified kinematic pattern:

```text
J(r,φ,t) = J_PT(r) [1 + c Σ w_j(r,t) A_j(r,φ,t)] / [1 + c Σ w_j(r,t)]
A_j = exp[κ cos(φ − φ_j − Ω(r) age_j)] / I₀(κ)
```

The angular factor has mean one. Nonnegative compact radial weights and smooth sin⁴ lifetime envelopes produce evolving bright and dark lanes without negative intensity. Different radii move at their own unrounded Keplerian frequency, so structures shear. They are evaluated at each ray's retarded emission time.

Current choices are 64 radial bands, logarithmic half-width 0.04, κ = 2.5, contrast c = 1.6 and reproducible seed 13326. Lifetimes are 1.5 local orbital periods, capped at 0.8 of the 800 t_g loop. The birth schedule targets two overlapping events per band. Full values and generated event counts are in the render metadata.

These parameters are illustrative choices, not quantities measured from NASA's movie or derived from magnetohydrodynamics. They preserve the prescribed mean radiative profile but do not solve local energy transport. There is no simulated MRI, magnetic field, radial accretion flow or gas heating.

Births repeat periodically and patches vanish smoothly before repeating. **It is the source pattern that is forced to repeat; gas orbital frequencies are not quantized, and the film does not claim that a turbulent disk physically repeats every eight seconds.**

## Display

Physical radiance is averaged over subpixels before any color or contrast transformation. The production render uses two samples per axis. The final layout reduction also averages radiance before tone mapping, preserving the image's proportions. A single fixed curve is applied throughout:

```text
x = relative observed radiance / 0.65
display level = x / (1 + x)
```

The color version maps this level through a monotone orange intensity palette. The monochrome companion uses the same level with the standard sRGB transfer function. Neither version adds bloom, an artificial photon rim, lens flare, star fields or a luminous halo. Relative brightness, not absolute emitted power, is shown. Contrast and palette are presentation choices; NASA's undocumented palette is not asserted to have the same definition.

## Verification and practical limits

`nasa-raytrace-diagnostics.json` records the prior fixed-view baseline checks: the exact flat-space limit of the time integral, a separate finite-segment coordinate-time quadrature, recovery of the earlier geometry, symmetry and numerical refinement. It also records full-size maps and higher-order disk intersections. `nasa-emission-diagnostics.json` checks positivity, mean-flux preservation, an independent direct expression, unrounded orbital frequencies, retarded sampling, periodicity and smoothness at the loop boundary.

`verify_full_orbit_geometry.py` generates `full-orbit-geometry-diagnostics.json`, which records 70 geometry and travel-time checks across the full revolution, both poles and exact disk-plane crossings, together with SHA-256 hashes of the tracer sources and verification script.

`verify_nasa.py` checks the delivered videos and the combined rendering model, including the full viewpoint revolution and loop endpoints. The baseline ray-trace report is retained for provenance; `output-sideratlas/verification.json` records the checks for the files and source hashes it identifies. A prior report does not validate a newly edited source or an unfinished render. Passing checks establish numerical consistency within the stated assumptions. They are not evidence that the toy emission pattern reproduces a real accretion flow. Fine image structure remains limited by pixel sampling and video compression.

Omitted physics includes finite disk thickness, magnetic turbulence and fluid evolution, radial inflow, plunging-region emission, returning-radiation heating, a corona, jets, scattering, polarization and a full atmosphere. At high inclination, finite thickness in particular can change which inner regions are visible. No claim of a uniquely realistic or exact NASA reconstruction is made.

## Reproduction

Python needs NumPy and Pillow. FFmpeg with libzimg/zscale and FFprobe must be on PATH; no C++ compiler is needed for this release.

```sh
python3 render_nasa.py --preview --width 540 --supersample 1.5
python3 render_nasa.py
python3 verify_nasa.py
```

The export saves per-frame radiance in the output directory's `.render-cache/`, keyed by the physical calculation, model and detector resolution. Typography and encoding are applied afterward, so changing them does not require tracing identical light paths again. An interrupted run can resume matching frames. The optional `--radiance-cache` argument accepts a prior cache only with a matching physical signature and a manifest of verified frame hashes. This release reused the original floating-point radiance; its provenance is in `output-sideratlas/provenance/`. No pixels were taken from the older compressed movie. The cache itself is omitted from the download; a normal `python3 render_nasa.py` reconstructs it. Final MP4s replace temporary files only after encoding completes. Preserve the Python modules and `fonts/` directory together. See `INSTAGRAM-EXPORT.md` for the delivery encoding and upload settings; Instagram performs its own processing after upload.

## Sources and their actual roles

1. **Page & Thorne (1974), _Disk-Accretion onto a Black Hole. I._** Mean radial flux with zero inner stress. [DOI](https://doi.org/10.1086/152990).
2. **Luminet (1979), _Image of a spherical black hole with thin accretion disk_.** Schwarzschild thin-disk imaging and transfer reference. [Paper](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L).
3. **Müller & Frauendiener (2012), _Interactive visualization of a thin disc around a Schwarzschild black hole_.** Schwarzschild rendering and propagation reference. [Paper](https://arxiv.org/abs/1206.4259).
4. **James et al. (2015), _Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar_.** Appearance/transfer background for the inherited ray tracer; DNGR is not used. [Paper](https://arxiv.org/abs/1502.03808).
5. **Bozza (2002), _Gravitational lensing in the strong field limit_.** Strong-deflection validation reference for the inherited Schwarzschild baseline. [Paper](https://arxiv.org/abs/gr-qc/0208075).
6. **Gralla & Lupsasca (2020), _The Null Geodesics of the Kerr Exterior_.** General geodesic background for the new time-integral documentation; the regularization above is derived directly, and their elliptic-function solver is not used. [Paper](https://arxiv.org/abs/1910.12881).
7. **Schnittman & Bertschinger (2004), _The Harmonic Structure of High Frequency Quasi-periodic Oscillations in Accreting Black Holes_.** Orbiting and sheared hot spots as a kinematic modeling precedent, not a copied source realization. [Paper](https://arxiv.org/abs/astro-ph/0309458), [DOI](https://doi.org/10.1086/383180).
8. **Cárdenas-Avendaño, Lupsasca & Zhu (2023), _Adaptive Analytical Ray Tracing of Black Hole Photon Rings_.** Variable equatorial emission, retarded source times and orbital scaling. Neither AART nor its inoisy/Matérn stochastic solver is used. [Paper](https://arxiv.org/abs/2211.07469).
9. **NASA Goddard / Jeremy Schnittman (2019), _Black Hole Accretion Disk Visualization_, SVS 13326.** Visual and educational reference: non-rotating hole, distorted thin disk, sheared structures and full viewpoint revolution above and below the disk. NASA did not produce or endorse this independent code. [Source](https://svs.gsfc.nasa.gov/13326/).
10. **International Color Consortium, sRGB interpretation.** Display encoding only. [Specification](https://registry.color.org/rgb-registry/files/sRGB.pdf).

The earlier Kerr/thermal release additionally used Bardeen et al., Li et al., Done et al., Airy and CIE data. Those models are not applied in this new bolometric release; their credit remains with the earlier outputs. Font sources and licenses are in `fonts/`.
