# @sideratlas — original Kerr disk Reel

This version changes the spacetime to a rotating Kerr black hole, recalculating the light paths, source azimuths, propagation delays, orbital motion and mean radial emission. The chosen dimensionless spin is **a* = Jc/(GM²) = +0.9**. The disk is equatorial and **aligned, prograde**: its angular momentum points along the black hole's spin. These are stated model choices, not measured properties of a particular object.

NASA Goddard / Jeremy Schnittman's [2019 visualization](https://svs.gsfc.nasa.gov/13326/) remains the visual reference for distorted disk images, sheared structures and a full viewpoint revolution. NASA's reference models a non-rotating black hole; the present Kerr configuration is an independent extension. No NASA footage, textures, data or code are incorporated. This is a relativistic ray-traced visualization with prescribed emission, not a GRMHD simulation or a reconstruction of NASA's unpublished setup.

## Deliverables

- `output-sideratlas-kerr/sideratlas-kerr-black-hole-reel.mp4`: orange illustrative brightness palette
- `output-sideratlas-kerr/sideratlas-kerr-black-hole-intensity.mp4`: monochrome rendering of the same calculated bolometric radiance
- `output-sideratlas-kerr/render-metadata.json`: configuration and source hashes
- `output-sideratlas-kerr/verification.json`: validation of the finished release identified in that report
- `caption-kerr-orbit.txt`: caption, assumptions and source credits

The planned export remains 1080 × 1920, 30 fps and eight seconds, with @sideratlas, NO.001 and THE BLACK HOLE. The added explanatory line is **“Spin drags spacetime around the black hole”**. The information strip identifies Kerr / a* = 0.9; the constants line states the aligned, prograde disk and r_ISCO ≈ 2.32 r_g. Color continues to be labeled illustrative. See `INSTAGRAM-EXPORT.md` for encoding and upload settings.

## Physical model and explicit choices

| Component | Calculation or choice |
| --- | --- |
| Spacetime | Kerr, a* = +0.9; fixed vacuum geometry |
| Disk alignment | Equatorial, zero tilt; angular momentum parallel to the hole's spin |
| Disk surfaces | Opaque and infinitesimally thin; identical isotropic prescribed emission on both faces |
| Emitting range | Prograde ISCO ≈ 2.32088 r_g to chosen outer cutoff 28 r_g |
| Inner boundary | Zero stress at the ISCO; no emission inside it |
| Mean intensity | Spin-dependent Page–Thorne flux, divided by its radial peak |
| Gas motion | Circular prograde Kerr geodesics; Ω = 1/(r^(3/2) + a) in G = M = c = 1 |
| Transfer | I_obs = g⁴ I_em for bolometric intensity; first opaque disk intersection wins |
| Variable structure | Positive, finite-lived, shearing brightness patches with an imposed periodic birth schedule |
| Observer | Sequence of static observers at infinity; one 360° viewpoint scan perpendicular to the disk |
| Display | Fixed exposure and monotone illustrative orange palette, with a monochrome companion |
| Physical scale | r_g = GM/c², t_g = GM/c³; no numerical mass or accretion rate assigned |

For this spin the horizon's Boyer–Lindquist radius is r_+ = 1 + sqrt(1 − a²) ≈ 1.43589 r_g. That coordinate radius is not the observed size of the dark region. Lensing, the emitting inner edge and foreground occultation determine the image. Compared with the former a* = 0 disk, emission can extend inward from 6 r_g to about 2.32 r_g. The image is not simply a smaller disk pasted onto the previous shadow.

Aligning the disk keeps the standard circular thin-disk model applicable. A tilted or precessing disk would require further dynamical assumptions. The spin value is substantial but below the extremal limit. It is not a claim that every astrophysical black hole has this spin.

## Original Kerr integration

The ray calculation uses the separated Kerr null-geodesic equations, including radial and polar turning points. The screen constants for oriented viewing angle θ are

```text
λ = −α sin θ
η = β² + (α² − a²) cos² θ
Δ = r² − 2r + a²
R(r) = [r² + a² − aλ]² − Δ[η + (λ − a)²]
```

Here θ labels the signed camera orientation, while the spacetime polar coordinate along a ray is separate. The camera uses n = (0, −sin θ, cos θ), screen-right (1, 0, 0) and screen-up (0, cos θ, sin θ). This preserves orientation through both poles and below the disk; images are not mirrored to fabricate the return motion.

The implementation integrates backward from infinity. Azimuth therefore advances with the negative of the future-directed geodesic derivative; elapsed photon travel time increases positively. The angular part is represented by a unit vector, avoiding the Boyer–Lindquist azimuth singularity at either pole. The residual radial contribution to azimuth includes frame dragging. These are original numerical implementations of the equations, not the elliptic-function solver supplied by Gralla & Lupsasca or AART.

The circular emitter's redshift factor is

```text
uᵗ = (r^(3/2) + a) / sqrt(r³ − 3r² + 2a r^(3/2))
g = 1 / [uᵗ (1 − Ωλ)]
```

The source is sampled at **t_em = t_obs − delay** for every disk hit. One common infinite-observer time offset is removed. Delays are not separately normalized for individual frames or image orders, and negative relative delays are valid.

For reproducibility, the regularized elapsed-time variable is q = ℓ + r*_Kerr(r) − r*_Kerr(observer). In units M = 1, define u = 1/r, D = 1 − 2u + a²u², P = 1 + (a² − aλ)u², K = η + (λ − a)² and v = du/dτ for backward Mino parameter τ. Its incoming derivative is evaluated in cancellation-safe form:

```text
dq/dτ = (1 + a²u²) K/(P + v) + aλ − a²(1 − z²)
delay = q − r*_Kerr(r_hit)
```

This regularization is derived directly from the geodesic equations for this implementation. The Kerr tortoise coordinate uses logarithms with denominator 2, so its a → 0 normalization matches the prior Schwarzschild implementation. Outgoing rays use the algebraically equivalent direct expression to avoid a small denominator. Numerical integration accuracy and finite pixel sampling still limit fine higher-order structure.

## Emission and motion

The mean radial profile is recomputed for the chosen spin using the zero-inner-stress Page–Thorne flux integral and Kerr circular-orbit quantities. Its absolute scale is removed by one global peak normalization. Consequently no physical temperature, thermal spectrum, mass or accretion rate is assigned. The orange mapping represents relative total-light brightness, **not true visible color**. The optional physical-temperature methods in older library code are not part of this release's image formation.

The retained original patch model evaluates each radius at its Kerr angular velocity, including in patch lifetimes and retarded-time sampling. Its angular factor has mean one; normalization preserves the prescribed radial mean at every source time. Widths, contrast, seed and lifetimes are illustrative choices, not quantities inferred from NASA's movie or simulated plasma dynamics.

The 800 t_g source-birth schedule repeats in the eight-second film, allowing a loop without quantizing orbital frequencies. Each patch dies smoothly before the schedule repeats. Real turbulent gas is not claimed to repeat. Inner accessible orbits have shorter periods, but one should not interpret spin as increasing Ω at every fixed radius: for the chosen positive prograde branch, the denominator at fixed r is r^(3/2) + a.

The viewpoint angle follows θ = 83° − 360° phase. The source keeps evolving forward while the view passes above and below the disk. This is a virtual scan through static observers, not a physical camera trajectory; observer-velocity aberration is omitted. A razor-thin disk's direct image has zero height exactly edge-on while bent images may persist. The same idealization becomes particularly consequential when viewing the innermost Kerr disk at high inclination.

## Limits and verification

This model omits finite disk thickness, disk self-gravity, magnetic fields and MRI, hydrodynamic evolution, radial inflow, plunging-region emission, a corona, jets, scattering, polarization, disk-atmosphere physics and returning-radiation heating. Returning radiation and finite thickness can matter for rapidly rotating holes. The spin upgrade adds physical effects within a specified idealized model; it does not establish a universally most realistic image.

**Release verification is pending during implementation.** The final report must identify the new Kerr source files and finished output hashes. Prior Schwarzschild or steady Kerr reports are useful benchmarks, not validation of an unfinished or changed movie. Required checks include the a = 0 limit for radius, azimuth and delay; Kerr numerical refinement; both poles and both disk-plane crossings; conserved geodesic quantities; periodic emission; and inspection/decoding of the final exports. Replace this paragraph with the actual final outcome after those checks run.

## Reproduction

The intended entry points are `render_kerr_orbit.py` and `verify_kerr_orbit.py`. Python needs NumPy and Pillow; the supplied C++ ray integrator needs clang++ or g++. FFmpeg needs libx264 and libzimg/zscale, and FFprobe must be available. Source code, fonts and their licenses must stay together. Kerr radiance is newly calculated; the old Schwarzschild frame cache is not valid for it.

## Sources and their roles

1. **Bardeen, Press & Teukolsky (1972), _Rotating Black Holes: Locally Nonrotating Frames, Energy Extraction, and Scalar Synchrotron Radiation_** — circular prograde motion and ISCO. [Paper](https://adsabs.harvard.edu/pdf/1972ApJ...178..347B), [DOI](https://doi.org/10.1086/151796)
2. **Page & Thorne (1974), _Disk-Accretion onto a Black Hole. I._** — spin-dependent mean radial flux with zero inner stress. [Paper](https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P), [DOI](https://doi.org/10.1086/152990)
3. **Gralla & Lupsasca (2020), _Null Geodesics of the Kerr Exterior_** — Kerr geodesic potentials and coordinate equations; numerical solver and regularization here are original. [Paper](https://arxiv.org/abs/1910.12881)
4. **Cárdenas-Avendaño, Lupsasca & Zhu (2023), _Adaptive Analytical Ray Tracing of Black Hole Photon Rings_** — variable equatorial emission, retarded times, and circular-orbit/transfer formulas; neither AART nor its stochastic source generator is used. [Paper](https://arxiv.org/abs/2211.07469)
5. **Schnittman & Bertschinger (2004), _The Harmonic Structure of High Frequency Quasi-periodic Oscillations in Accreting Black Holes_** — orbiting/sheared hot spots as kinematic precedent. [Paper](https://arxiv.org/abs/astro-ph/0309458)
6. **Luminet (1979), _Image of a spherical black hole with thin accretion disk_** — Schwarzschild thin-disk imaging background and inherited benchmark. [Paper](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L)
7. **Müller & Frauendiener (2012), _Interactive visualization of a thin disc around a Schwarzschild black hole_** — inherited Schwarzschild propagation benchmark. [Paper](https://arxiv.org/abs/1206.4259)
8. **Bozza (2002), _Gravitational lensing in the strong field limit_** — inherited Schwarzschild strong-deflection validation reference. [Paper](https://arxiv.org/abs/gr-qc/0208075)
9. **James et al. (2015), _Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar_** — Kerr image and transfer background; DNGR is not used. [Paper](https://arxiv.org/abs/1502.03808)
10. **NASA Goddard / Jeremy Schnittman (2019), SVS 13326** — visual reference only; its documented black hole is non-rotating. NASA did not produce or endorse this code. [Source](https://svs.gsfc.nasa.gov/13326/)
11. **International Color Consortium, sRGB interpretation** — display encoding; [specification](https://registry.color.org/rgb-registry/files/sRGB.pdf). Rec.709 delivery conversion is documented in `INSTAGRAM-EXPORT.md`

Older thermal, atmosphere and diffraction studies remain credited with the earlier outputs that used them. They are not claimed as applied models in this bolometric Kerr release.
