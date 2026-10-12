# @sideratlas — Kerr–Newman disk Reel

This release depicts the requested **Kerr–Newman black hole with a* = +0.9 and q = +0.4**. The disk is equatorial and prograde: its angular momentum is parallel to the black hole's spin. Light paths, emitter positions, propagation delays, circular motion and the mean radial emission profile are recalculated for that geometry.

**This is a hypothetical large-charge configuration with an idealized neutral test disk.** It is a visualization of a specified relativistic model, not a prediction for an ordinary astrophysical accretion plasma. The charge is deliberately large enough to affect spacetime appreciably. Plausible charges in astrophysical environments are generally much smaller; even a small net charge can nevertheless affect charged particles. The Galactic-centre study by [Zajaček et al. (2018)](https://arxiv.org/abs/1808.07327) discusses this distinction. No stability or long-term survival of the chosen charge is claimed here.

NASA Goddard / Jeremy Schnittman's [2019 visualization](https://svs.gsfc.nasa.gov/13326/) remains the visual reference for the full viewpoint revolution and distorted disk images. That reference depicts a non-rotating black hole. This original code does not use NASA footage, textures, simulation data or code, and NASA has not endorsed it.

## Deliverables

- `output-sideratlas-kn/sideratlas-kerr-newman-reel.mp4`: orange illustrative brightness, 32 seconds
- `output-sideratlas-kn/sideratlas-kerr-newman-intensity.mp4`: the same calculated radiance in monochrome, 32 seconds
- `output-sideratlas-kn/render-metadata.json`: numerical configuration, model assumptions and provenance
- `output-sideratlas-kn/verification.json`: checks of the identified finished release
- `caption-kerr-newman.txt`: ready-to-copy caption with assumptions and source credits

The delivery is **four identical repetitions of one eight-second loop**, not 32 seconds of newly simulated source evolution. Each loop contains one 360° virtual viewpoint scan, with source time advancing by 800 GM/c³. The movie remains 1080 × 1920 at 30 fps. The @sideratlas identity, LaTeX-like fonts, black background and Instagram text margins are retained. See `INSTAGRAM-EXPORT.md` for encoding and upload guidance.

The two added explanatory lines are:

> Spin drags spacetime around the black hole
>
> Electric charge changes spacetime and light paths

The information strip identifies **Kerr–Newman: a* = 0.9, q = 0.4** and illustrative color. A nearby line discloses the hypothetical large charge and neutral test disk. Photons have no electric charge: the additional deflection here follows the charged black hole's spacetime geometry, not electrical attraction of light.

## Parameter convention and disk model

Here a* = Jc/(GM²), q = Q/[M√(4πε₀G)], r_g = GM/c² and t_g = GM/c³. In geometrized units G = M = c = 1, the numerical parameters are a = 0.9 and Q = q = 0.4. The dimensionless q is **not a charge in coulombs** and is not the fraction of the maximum charge at fixed spin.

The black-hole condition is a*² + q² ≤ 1. This choice gives 0.97, so the model has an outer horizon at Boyer–Lindquist radius r₊ = 1 + √(1 − a*² − q²) ≈ 1.173205 r_g. The maximum |q| at this spin is √0.19 ≈ 0.435890. These coordinate radii are not the observed radius of the central dark region. The calculated neutral prograde ISCO is **1.577505 r_g**, compared with 2.320883 r_g for uncharged Kerr at the same spin. [Pugliese, Quevedo & Ruffini (2013)](https://arxiv.org/abs/1303.6250) supplies the neutral-orbit reference.

| Component | Calculation or explicit choice |
| --- | --- |
| Background | Fixed Kerr–Newman electrovac geometry, a* = +0.9, q = +0.4 |
| Disk alignment | Equatorial, zero tilt, prograde; disk and hole angular momenta both along +z |
| Emitting matter | Neutral test matter on circular geodesics; no Lorentz force, currents or plasma response |
| Surface | Opaque, infinitesimally thin, identical isotropic prescribed emission on both faces |
| Inner edge | Calculated neutral prograde ISCO; zero inner stress and no emission inside it |
| Outer edge | Chosen radius 28 r_g |
| Mean brightness | General Page–Thorne thin-disk flux integral evaluated with Kerr–Newman circular-orbit quantities, normalized by one radial peak |
| Transfer | Bolometric I_obs = g⁴ I_em; first opaque emitting intersection wins |
| Variable structure | Positive finite-lived shearing brightness patches with a repeating birth schedule |
| Observer | A sequence of static observers at infinity, scanned through a complete revolution |
| Display | Fixed exposure 0.02 in peak-normalized radiance units; monotone orange illustrative palette or monochrome intensity |
| Physical scale | No numerical mass, accretion rate, thermal spectrum or temperature assigned |

The disk prescription does not solve the dynamics of charged accretion gas. Treating matter as neutral test emitters avoids silently choosing a particle charge-to-mass ratio or a magnetosphere, but remains a consequential idealization at q = 0.4.

## Recalculated paths and delays

The tracer uses Kerr–Newman separated null-geodesic equations from [Wang, Lee & Lin (2022)](https://arxiv.org/abs/2208.11906), equations 9–15 and 21. Its radial potential is

```text
Δ = r² − 2r + a² + q²
R = [r² + a² − aλ]² − Δ[η + (λ − a)²]
λ = Lz/Eγ,   η = Carter_constant/Eγ²
```

The polar potential retains the Kerr form. Charge also changes the radial azimuth and travel-time terms through Δ. The implementation is an original numerical integrator, not the paper's elliptic-function solver. A pole-regular angular vector maintains a continuous camera basis through both poles. Numerical integration follows turning points and equatorial crossings; separate rings are not pasted into the image.

The circular source has g = 1/[uᵗ(1 − Ωλ)] and is sampled at **t_em = t_obs − delay**. One common infinite-observer time offset is removed for all pixels, image orders and viewpoints. Negative relative delays are valid. The original regularized delay integral is recomputed with the Kerr–Newman tortoise coordinate, including its charge-dependent horizon coefficients. The internal variable named `q` in that integration is a time regularization variable; `charge_q` denotes the black-hole charge parameter.

The signed viewing angle is θ = 83° − 360° phase. The screen basis remains continuous as the view passes above and below the disk. Source time continues forward throughout. This is a virtual scan among static observers, not a physical camera moving at that speed; observer-motion aberration is not included.

The production render traces 120 independent detector maps for the first half-revolution. The other 120 views use an exact antipodal symmetry of this aligned, two-sided model: reversing the detector's horizontal coordinate preserves the ray constants, radius, redshift and delay, while the emitting azimuth gains π. This reuses calculated geometry, not rendered images. Emission is evaluated anew at each view's own azimuth and retarded time, producing **240 distinct time samples** per eight-second cycle. No angular or temporal image interpolation is used. The completed encoded cycle is then repeated four times without a second lossy encode.

## Recalculated circular motion and emission

For the selected neutral prograde branch, in G = M = c = 1,

```text
s = √(r − q²)
Ω = s / (r² + a s)
uᵗ = (r² + a s) / [r √(r² − 3r + 2q² + 2a s)]
```

The neutral circular orbit's specific energy E and angular momentum L are recalculated from the metric. Marginal radial stability determines the ISCO. The implementation checks the neutral radial potential and the q = 0 Kerr and a = q = 0 Schwarzschild limits. Neutral geometry is unchanged by reversing the sign of q; charged particles would not share this symmetry.

The [Page–Thorne (1974)](https://doi.org/10.1086/152990) conservation-law derivation begins with a general stationary, axisymmetric, equatorially symmetric geometry before specializing to Kerr. Here its zero-inner-stress flux integral is evaluated using the neutral Kerr–Newman quantities:

```text
F(r) ∝ −Ω,r / [r(E − ΩL)²] × ∫(r_ISCO to r) (E − ΩL)L,r dr
```

The equatorial determinant factor is r in the near-disk coordinates used by that integral. The old Kerr closed-form radial profile is not reused. Peak normalization removes the unassigned luminosity scale. Consequently orange remains **illustrative total-light brightness, not true visible color**.

The display first averages the calculated radiance over the detector samples and presentation pixels. A single curve, x/(1+x) with x = radiance/0.02, then sets the brightness for every viewpoint and both display modes. There is no frame-by-frame auto-exposure. The value 0.02 is a presentation choice, not a physical luminosity, and this peak-normalized animation cannot establish which black-hole model is intrinsically brighter.

The retained patch prescription uses the new Ω(r) in shear, lifetimes and retarded-time sampling. Patches have positive brightness and a normalization that preserves the prescribed azimuthal radial mean. They die smoothly, and their birth schedule repeats every 800 t_g to close the loop without quantizing orbital frequencies. This is prescribed kinematic structure, not turbulence generated by fluid or magnetic-field equations.

## Scope and reproducibility

This model includes the electromagnetic field's effect on the background metric. It does not simulate disk charge separation, Lorentz forces, conductivity, currents, an external magnetosphere, plasma refraction, polarization or discharge of the black hole. It also omits fluid evolution, finite disk thickness, radial inflow, emission inside the ISCO, disk self-gravity, a corona, jets, scattering and returning-radiation heating. The electromagnetic field of a rotating charge exists in the Kerr–Newman solution; omitting the disk's electromagnetic response does not mean that field is absent from the spacetime.

The calculation and export entry point is `render_kn_orbit.py`. The new physical modules are `kn_orbit_trace.py`, `kn_orbit_trace.cpp`, `kn_orbits.py` and `kn_orbit_emission.py`. Python requires NumPy and Pillow; the supplied C++ kernel requires clang++ or g++. FFmpeg must support libx264 and libzimg/zscale, and FFprobe must be available. Retain bundled fonts and licenses when reproducing the output. Neither prior Schwarzschild radiance nor an uncharged Kerr cache validates or supplies this charged render.

To reproduce the release after installing `requirements.txt` and the external tools above:

```sh
python3 verify_kn_trace.py
python3 verify_kn_orbit_emission.py
python3 render_kn_orbit.py
python3 verify_kn_orbit.py
```

`--preview --width 540 --supersample 1 --output output-sideratlas-kn/preview` makes a small still-image preview. The default render independently traces a 2160 × 1800 detector for 120 views, evaluates 240 radiance frames, averages the radiance before the display transform and exports both 1080 × 1920 versions. A parameter-and-source-keyed floating-point cache allows an interrupted render to resume; no earlier Schwarzschild/Kerr frames are reused.

The finished release passed **62 ray-tracing checks**, **69 orbital/emission checks**, and **75 renderer, provenance and delivery checks**. Reports identify their exact source hashes. The independent checks include uncharged Kerr and Schwarzschild limits, charge-sign symmetry, analytic capture boundaries, radial stability, numerical refinement, both poles and disk-plane crossings, and independently integrated azimuth and arrival times. These are numerical consistency tests of the declared idealized model, not validation of charged accretion plasma.

Both final videos decode all **960 frames at 1080 × 1920, 30 fps, exactly 32 seconds**. Decoded frame hashes confirm four pixel-identical copies of the original eight-second cycle. The loop boundary is consistent with nearby motion. The delivery uses H.264 High, limited-range Rec.709 and fast-start MP4. Twelve color-managed decoded color frames were also compared with their original radiance compositions: full-frame mean absolute RGB differences range from 0.115 to 0.199 on a 0–255 scale; the largest lit-disk mean difference is 2.116. This measures the local encode, not Instagram's subsequent processing.

The sampled orbit contact sheet and full-frame review confirm that the changing disk views and text fit the intended layout. The overlay spans y = 270 through 1618, leaving 301 pixels clear below. The exact checks, video hashes and sample comparisons are in `output-sideratlas-kn/verification.json` and `visual-review.json`.

## Sources and their roles

- **Newman et al. (1965)** — historical charged, rotating metric underlying the modern geodesic formulation. [DOI](https://doi.org/10.1063/1.1704351)
- **Kerr (1963)** — historical uncharged, rotating metric used as the q = 0 limit. [DOI](https://doi.org/10.1103/PhysRevLett.11.237)
- **Wang, Lee & Lin (2022)** — Kerr–Newman null geodesic equations; original numerical implementation here. [arXiv:2208.11906](https://arxiv.org/abs/2208.11906), [DOI](https://doi.org/10.1103/PhysRevD.106.084048)
- **Pugliese, Quevedo & Ruffini (2013)** — neutral equatorial circular orbits and stability in Kerr–Newman. [arXiv:1303.6250](https://arxiv.org/abs/1303.6250)
- **Ivanov & Prodanov (2005)** — exact neutral circular-orbit energy and angular momentum, equations 19–20; the paper's charged-particle approximation and pseudo-Newtonian potential are not used. [arXiv:gr-qc/0504025](https://arxiv.org/abs/gr-qc/0504025)
- **Liu, Lee & Lin (2017)** — independent neutral ISCO equation, equation 22. Stability is checked directly against the radial potential rather than relying on the paper's prose sign convention. [arXiv:1706.05466](https://arxiv.org/abs/1706.05466)
- **Page & Thorne (1974)** — general mean thin-disk flux integral and zero-inner-stress assumption. [DOI](https://doi.org/10.1086/152990)
- **Zajaček et al. (2018)** — context for astrophysical electric charge and the distinction between geometric and charged-particle effects; not a q = 0.4 disk model. [arXiv:1808.07327](https://arxiv.org/abs/1808.07327)
- **Bardeen, Press & Teukolsky (1972)** — uncharged Kerr circular-orbit/ISCO limit. [DOI](https://doi.org/10.1086/151796)
- **Gralla & Lupsasca (2020)** — inherited Kerr null-geodesic background and q = 0 comparison; their solver is not used. [arXiv:1910.12881](https://arxiv.org/abs/1910.12881)
- **Cárdenas-Avendaño, Lupsasca & Zhu (2023)** — inherited retarded equatorial-source and radiative-transfer approach; AART and its source generator are not used. [arXiv:2211.07469](https://arxiv.org/abs/2211.07469)
- **Schnittman & Bertschinger (2004)** — kinematic hot-spot precedent, not a simulated plasma texture. [arXiv:astro-ph/0309458](https://arxiv.org/abs/astro-ph/0309458)
- **James et al. (2015)** — Kerr lensing and transfer background; DNGR is not used. [arXiv:1502.03808](https://arxiv.org/abs/1502.03808)
- **Luminet (1979)** — inherited Schwarzschild disk-image benchmark. [A&A 75, 228](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L)
- **Müller & Frauendiener (2012)** — inherited Schwarzschild propagation benchmark. [arXiv:1206.4259](https://arxiv.org/abs/1206.4259)
- **Bozza (2002)** — inherited Schwarzschild strong-deflection benchmark. [arXiv:gr-qc/0208075](https://arxiv.org/abs/gr-qc/0208075)
- **NASA/Goddard/J. Schnittman (2019)** — visual reference only. [SVS 13326](https://svs.gsfc.nasa.gov/13326/)
- **ICC sRGB** — display encoding. [Specification](https://registry.color.org/rgb-registry/files/sRGB.pdf). Rec.709 delivery conversion is documented separately in `INSTAGRAM-EXPORT.md`

Older thermal, disk-atmosphere and diffraction studies remain credited with earlier outputs that used them. They are not applied models in this bolometric Kerr–Newman release.
