# Sideratlas v1.2.0 — accretion, plunge and multiple images

Version 1.2.0 extends the original **Kerr–Newman visualization with $a_*=+0.9$ and $q=+0.4$**. It adds restrained evolving disk structure, slow radial inflow, luminous neutral gas plunging inside the ISCO, and radiative transfer through multiple images of that inner gas. The full-frame diagonal composition, typography, fixed exposure, black sky and smooth viewpoint revolution are retained.

**The spacetime, photon paths, timelike plunge, frequency shifts and delays are calculated. Brightness fluctuations, outer drift, inner heating and optical depth are declared prescriptions.** This remains an idealized neutral test-matter visualization in a hypothetical large-charge geometry. It is not a general relativistic magnetohydrodynamic (GRMHD) or thermodynamic accretion simulation.

## Delivery and version history

The [v1.2.0 release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.2.0) supplies the 32-second orange Reel, its monochrome companion, eight-second originals, description and complete project archive. Current media and reports are in `output-sideratlas-kn-v120/`. The video names remain `sideratlas-kerr-newman-reel.mp4`, `sideratlas-kerr-newman-intensity.mp4`, and their corresponding `-loop.mp4` files.

Each 32-second movie repeats the same encoded eight-second cycle four times. A cycle contains a full virtual viewpoint revolution and 800 $t_g$ of periodically driven source evolution. It is not 32 seconds of distinct physical evolution. The only on-screen words are **@sideratlas**, **NO.001**, and **THE BLACK HOLE**, over the simulation. Explanations, caveats, source credits and the repository link are supplied in [caption-kerr-newman.text](caption-kerr-newman.text); post that description with the video.

The [v1.1.0 release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.1.0) preserves the earlier circular-emitter model and diagonal layout in `output-sideratlas-kn/`. Its entry point remains `render_kn_orbit.py`, documented in [README-KERR-NEWMAN.md](README-KERR-NEWMAN.md). The [v1.0.0 release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.0.0) retains the original smaller horizontal composition and lower text. Their reports apply to those historical versions.

## Geometry and declared choices

The conventions are $a_*=Jc/(GM^2)$, $q=Q/[M\sqrt{4\pi\epsilon_0G}]$, $r_g=GM/c^2$, and $t_g=GM/c^3$. Calculations use $G=M=c=1$, metric signature $-+++$. The dimensionless charge is not measured in coulombs or given as a fraction of the spin-dependent maximum charge.

The chosen parameters give $a_*^2+q^2=0.97<1$. The outer horizon is

$$r_+=1+\sqrt{1-a_*^2-q^2}=1.1732050808\,r_g$$

The neutral prograde ISCO is $r_I=1.5775045665\,r_g$. These are Boyer–Lindquist coordinate radii, not apparent image radii. All gas angular momentum is aligned with the hole's positive spin. Reversing the sign of the background charge leaves neutral trajectories unchanged; this would not hold for charged gas.

| Component | Calculation or choice |
| --- | --- |
| Background | Stationary Kerr–Newman electrovac spacetime, fixed $a_*=0.9$, $q=0.4$ |
| Emitters | Neutral test matter, no electromagnetic force on the gas |
| Geometry of emission | Infinitesimally thin equatorial sheet, aligned and prograde, outer emitting edge $28r_g$ |
| Outer motion | Circular Kerr–Newman angular frequency plus a small normalized radial boost |
| Outer mean brightness | Original zero-ISCO-stress Page–Thorne radial profile, normalized to peak one |
| Irregularities | Positive multiscale localized packets and a coherent spiral, exact azimuthal normalization, smooth finite-lived births |
| Inner motion | Neutral timelike geodesic launched with the finite outer-flow velocity at the ISCO |
| Inner source brightness | Separately prescribed weak heating, limiting face intensity 0.08 of the outer radial peak |
| Inner absorption | Homogeneous gray slab with prescribed comoving normal optical depth 0.15 |
| Outer absorption | First opaque outer-disk intersection terminates the backward ray |
| Transfer | Per-event retardation, full flowing-emitter frequency shift, $g^4$ bolometric intensity and foreground attenuation |
| Observer | Sequence of static observers at infinity; no aberration from a physically moving camera |
| Display | Fixed exposure, orange illustrative brightness or monochrome; no assigned temperature or visible spectrum |

The choice $q=0.4$ is a theoretical experiment. The code does not establish that such a charge can persist in an accreting astrophysical plasma. [Zajaček et al. (2018)](https://arxiv.org/abs/1808.07327) provides context on astrophysical charge scales and charged-particle effects; it does not supply this large-charge disk model.

## Outer accretion and preserved Page–Thorne profile

For neutral prograde circular motion, with $s=\sqrt{r-q^2}$,

$$\Omega_c=\frac{s}{r^2+as}$$

The original mean exterior flux is retained:

$$F_{\rm PT}(r)\propto\frac{-\Omega_{c,r}}{r(E_c-\Omega_cL_c)^2}\int_{r_I}^{r}(E_c-\Omega_cL_c)L_{c,r}\,dr$$

The general stationary-axisymmetric conservation integral of [Page & Thorne (1974)](https://doi.org/10.1086/152990) is evaluated with the neutral Kerr–Newman circular quantities. Its zero-stress edge and normalized radial shape are preserved outside the ISCO. No mass, accretion rate or temperature is assigned.

A small radial velocity is prescribed using the familiar thin-disk scaling

$$v_r=-\alpha(H/r)^2v_{\phi,\mathrm{ZAMO}},\qquad\alpha=0.1,\quad H/r=0.15$$

Here $v_{\phi,\mathrm{ZAMO}}$ is the circular speed relative to a local zero-angular-momentum observer. The scaling motivates a radial boost in the circular emitter's orthonormal radial direction, not a complete viscous solution. Its four-velocity is

$$u^\mu=\gamma_r\left(u_c^\mu+v_r e_{(r)}^\mu\right),\qquad e_{(r)}^r=\sqrt{\Delta}/r,\qquad\gamma_r=(1-v_r^2)^{-1/2}$$

This preserves $u^\mu u_\mu=-1$ and the circular angular frequency while making $u^r<0$. The launch speed at the ISCO is approximately $-0.00141298c$. The radial speed law, its parameters and the neutral-matter approximation are choices; the metric conversion and normalization are calculated. [Penna et al. (2013)](https://doi.org/10.1093/mnras/sts185) discusses the drift scaling. The assigned aspect ratio only sets that speed scale; finite vertical thickness is not ray traced.

Retaining a circular Page–Thorne flux as the brightness background of a weakly drifting flow is a thin-disk approximation. The calculation does not solve a coupled stress, surface-density, heating and radial-momentum system. In particular, the zero-stress flux boundary and the prescribed finite launch are not claimed as a self-consistent viscous disk solution.

## Multiscale structure and material advection

The brightness field combines eight broad and eighteen finer angular packets with a coherent $m=2$ logarithmic spiral. Packet birth radii, angular phases, slopes and radial widths vary, rather than placing brightness on fixed concentric annuli. Their compact radial supports have logarithmic half-widths 0.27–0.41 for broad features and 0.055–0.125 for fine features. Differential rotation winds these localized patches into irregular arcs and filaments. They are emissivity structures, not photon rings.

Positive angular kernels have a known azimuthal mean:

$$K_n(\theta)=\frac{(1+\cos\theta)^n}{\binom{2n}{n}/2^n},\qquad\langle K_n\rangle_\phi=1$$

Their finite degrees $n=1,2,4,8,12$ are blended with nonnegative weights. Each birth field is normalized at every radius:

$$M=\frac{1+C\sum_j w_jK_j}{1+C\sum_j w_j},\qquad C=3,\quad w_j\geq0$$

The coherent spiral is another positive, unit-mean component. Thus intensity remains positive and bounded while its azimuthal mean is exactly the unchanged exterior Page–Thorne background. $C$ is a relative kernel weight, not a fractional contrast. The normalization prevents a radial packet envelope from adding an axisymmetric brightness ring.

An angular-advection cadence guide broadens fast inner features by blending positive kernels, instead of suppressing their contrast wholesale. It estimates the largest angular speed along nine future samples of each 50-$t_g$ parcel lifetime, with a small safety margin. Finer angular features are allowed farther out; radial detail remains near the emission peak. This is a sampling prescription, not a predicted turbulent cutoff or a proof of complete temporal convergence. Radial transport, field births, retardation and changing viewpoints can also vary in time. The renderer keeps its three-point source shutter, checked against nine-point quadrature at selected views.

These are illustrative emissivity fluctuations, not calculated turbulent velocities, MRI or plasma eddies. Entire packets, including their radial envelopes and selected angular widths, are transported along the same radial and azimuthal flow used in the frequency shift. Let $T(r)$ and $A(r)$ be coordinate travel time and accumulated azimuth relative to the ISCO. Along a streamline,

$$\frac{dT}{dr}=\frac{u^t}{u^r},\qquad\frac{dA}{dr}=\frac{u^\phi}{u^r}$$

Both increase inward, since the radius decreases. For source-field age $s_b$, the birth radius satisfies $T(r_b)=T(r)-s_b$, and the birth azimuth is $\phi_b=\phi-A(r)+A(r_b)$. Thus structure spirals inward and shears through differential rotation, rather than moving a rendered texture across the image.

Seeded births occur every $25t_g$ and crossfade over a $50t_g$ lifetime with complementary smooth windows. The 32-birth schedule repeats after $800t_g$ to make a constructed loop. Orbital frequencies are not rounded to the movie period. The analytically normalized angular kernels and convex birth blending ensure that the emitted azimuthal mean at fixed radius and emitted time equals the original exterior radial profile exactly. Observed images need not preserve that mean because lensing, Doppler shifts and retarded times weight the source differently.

This method follows the idea of transported brightness structures and retarded sources used in [Schnittman & Bertschinger (2004)](https://arxiv.org/abs/astro-ph/0309458) and [Cárdenas-Avendaño et al. (2023)](https://arxiv.org/abs/2211.07469); their source generators and simulation data are not used.

### Why the initial v1.2 texture became too smooth

The original v1.1 source used 64 narrow radial patch supports with strong angular inhomogeneity. The first v1.2 experiment replaced them with only three broad modes and suppressed the finer modes near the brightest gas. At the Page–Thorne peak, matched emitted-time samples had azimuthal RMS variations of 0.069–0.102 of the mean, compared with 0.284–0.707 in v1.1. The unchanged fixed-exposure tone curve compressed those already weak variations further. Independent subpixel sampling and the source shutter were not the principal causes.

The restored source gives 0.266–0.686 on the same samples, a 5.97-fold increase in mean emitted RMS over the smooth experiment. Its radial structure remains less fine-grained than v1.1. The [source comparison report](kn-emission-comparison-diagnostics.json) defines the samples, Fourier bands, radial structure functions and local tone sensitivities. These are source-plane diagnostics, not measured MRI or claims about NASA's simulation.

The [comparison image](output-sideratlas-kn-v120/emission-comparison.jpg) evaluates the v1.1 brightness prescription, the smooth experiment and the restored source on identical current photon paths, flow, opacity, retarded times, sampling and display settings. The first column is therefore an isolated old-emissivity comparison, not a rerender of the complete historical model. Run `python3 compare_kn_emission.py` to reproduce the source statistics and `python3 compare_kn_accretion.py` for the three-column images. The smooth reference is reconstructed from its recorded parameters and checked against captured statistics; it is never used by the production renderer.

Disk texture and higher-order images are checked separately. The [image-order report](kn-image-order-diagnostics.json) isolates actual trajectory crossings through order 4 in one central detector fixture. Renaming their image-order labels leaves the brightness unchanged. No image-order boost, added ring or painted critical curve enters the animation. Physically wound disk features can resemble arcs or bands in projection; they are not automatically higher-order images.

## Exact neutral plunge and the black horizon

At the ISCO the finite radial boost sets $E_p=\gamma_rE_I$ and $L_p=\gamma_rL_I$. These match the outer four-velocity continuously. The resulting specific constants are approximately $E_p=0.7652802623$ and $L_p=1.6681549435$.

Inside the ISCO, the gas follows the neutral equatorial timelike first integrals, with zero Carter constant:

$$\Delta=r^2-2r+a^2+q^2,\quad P=E_p(r^2+a^2)-aL_p$$

$$R_p=P^2-\Delta[r^2+(L_p-aE_p)^2],\qquad u^r=-\sqrt{R_p}/r^2$$

$$u^t=\frac{(r^2+a^2)P/\Delta+a(L_p-aE_p)}{r^2},\qquad u^\phi=\frac{aP/\Delta+(L_p-aE_p)}{r^2}$$

These are the Kerr–Newman geodesic equations of [Wang, Lee & Lin (2022)](https://arxiv.org/abs/2208.11906) and [Ko, Lee & Lin (2024)](https://arxiv.org/abs/2309.13832), numerically evaluated here. The radial potential is evaluated with its ISCO triple-root structure to avoid subtraction errors near the slow launch. An exactly circular ISCO launch with zero radial velocity has no finite inward crossing; the small finite launch avoids that degeneracy without adding a radial force to the inner geodesic.

Frame dragging follows from the metric: the inward flow's coordinate angular velocity approaches $\Omega_H=a/(r_+^2+a^2)$ close to the horizon. Coordinate time grows strongly there, while proper infall time remains finite. The simulation samples the emitted pattern at each photon's own delayed source time; it does not slow an image by hand.

Zero-stress Page–Thorne theory does not predict the inner luminosity. A separate, explicitly chosen face-intensity profile is therefore used only for $r_+<r<r_I$:

$$I_{\rm face,plunge}=0.08(1-x)^2,\qquad x=\frac{r-r_+}{r_I-r_+}$$

It has zero value and first derivative at the ISCO and a finite exterior limit near the horizon. Fluctuations multiply this background in the same bounded way. The weak source demonstrates how luminous plunging gas would be transferred to the observer; it is not a temperature, density, dissipation or cooling prediction. [Mummery & Stone (2024)](https://arxiv.org/abs/2407.02164) supplies motivation for studying spiral transport inside the ISCO, but its thick-disk GRMHD calculation is not reproduced or calibrated here.

No source is assigned on or inside the event horizon. Photon capture and calculated frequency shifts produce the observed fading; there is no luminous horizon outline. An exterior numerical cutoff at $r_++10^{-4}r_g$ terminates tracing before the Boyer–Lindquist coordinate singularity. The innermost excluded sliver is a numerical boundary, not an emitting surface. A central dark image region depends on the surrounding emission and need not coincide with the projected mathematical horizon.

## Multiple images and transfer

The new tracer reuses the original pole-regular Kerr–Newman null equations and regularized propagation-time integral. For photon constants $\lambda=L_z/E_\gamma$ and $\eta=\mathcal Q/E_\gamma^2$,

$$R_\gamma=[r^2+a^2-a\lambda]^2-\Delta[\eta+(\lambda-a)^2]$$

At every emitting intersection, the forward **covariant** radial momentum retains its physical sign:

$$k_r/E_\gamma=\pm\sqrt{R_\gamma}/\Delta$$

The flowing emitter's frequency ratio is calculated from its photon-energy contraction:

$$g=\frac{E_{\rm obs}}{E_{\rm em}}=\frac{1}{u^t-\lambda u^\phi-(k_r/E_\gamma)u^r}$$

The radial term matters for plunging gas. The old circular expression is recovered as the drift tends to zero. Bolometric specific intensity transforms as $I_{\rm obs}=g^4I_{\rm em}$. Each intersection uses $t_{\rm em}=t_{\rm obs}-\mathrm{delay}$, with one common observer-time offset across rays, images and views.

An opaque outer disk hides later intersections. The inner plunging layer is assigned a homogeneous gray-slab normal optical depth $\tau=0.15$ so that a ray can carry contributions from several crossings. This optical depth is an illustrative choice, not inferred from an accretion density. At local surface-normal cosine $\mu$, the declared face intensity is converted to source function and oblique emission by

$$S=\frac{I_{\rm face}}{1-e^{-\tau}},\qquad I_{\rm em}(\mu)=S(1-e^{-\tau/\mu}),\qquad\mathcal T(\mu)=e^{-\tau/\mu}$$

The normal-view intensity therefore remains $I_{\rm face}$. The grazing limit is finite and opaque. For ordered events from the camera outward, the observed contribution is $g_i^4I_i$ multiplied by every foreground transmission. The outer surface is isotropic and opaque. Both sides share the same prescribed emission.

Rays naturally accumulate extra equatorial crossings when they linger near unstable photon trajectories. The stored crossing order determines the separate image contributions. Their apparent position, stretching, delay and brightness follow the traced paths and transfer; four stripes or prescribed ring widths are never inserted. Event storage is finite, but exhausted rays are retried with a larger capacity rather than silently discarding images.

Higher orders become progressively thinner near the metric's critical curve. Their visibility depends on emitting material, opacity, inclination and pixel scale. Several orders can be present numerically while merging into one band at 1080-pixel delivery width. The renderer does not amplify them to guarantee a particular band count. The distinction between broad lensed images, extremely narrow higher orders and source-dependent central darkness is discussed by [Gralla, Holz & Wald (2019)](https://arxiv.org/abs/1906.00873). [Vincent et al. (2022)](https://doi.org/10.1051/0004-6361/202244339) illustrates how absorption can suppress higher-order images; its synchrotron spectrum and thick-disk models are not used here.

## Sampling and cinematic presentation

The physical viewport fills the 9:16 frame on the retained fixed diagonal. The close crop can exclude outer-disk emission near the poles. The viewpoint is the original continuous signed-angle sweep through both poles, with a continuous detector basis. It represents a sequence of stationary observers at infinity, not a camera traveling along that orbit.

Independent $2\times2$ subpixel photon trajectories refine pixels near the analytically calculated critical curve and pixels containing higher-order or plunging-region contributions. Their floating-point radiances are averaged before the display transform. This resolves thin structures where possible without drawing rings, sharpening their brightness or interpolating across image orders. The critical curve guides extra numerical sampling only; it is not an added visible component.

Three temporal Gauss samples average the evolving source over a half-frame exposure, corresponding to a 180-degree shutter. The ray geometry is held at the frame's central view while source time and its retarded evaluations advance through the exposure. This is a source-motion blur approximation, not an exact blur integral over the moving viewpoint. It reduces temporal aliasing of the fast inner flow without changing orbital frequencies or interpolating encoded frames.

Both display modes use the same calculated radiance, fixed exposure and monotone tone curve. Orange denotes illustrative bolometric brightness; it is not the predicted color of thermal accretion light. No blooming horizon, background stars, painted lens flares or automatic exposure are added. The delivery remains 1080 × 1920 at 30 fps, silent H.264 with Rec.709 tagging and fast-start playback. See [INSTAGRAM-EXPORT.md](INSTAGRAM-EXPORT.md).

## Reproduction and numerical evidence

Install Python dependencies from `requirements.txt`, retain the supplied fonts and licenses, and provide a C++17 compiler plus FFmpeg/FFprobe with libx264 and zscale. The C++ tracing kernel compiles locally. The v1.2.0 entry point is `render_kn_accretion.py`; prior entry points continue to describe prior models.

```sh
python3 verify_kn_trace.py
python3 verify_kn_orbit_emission.py
python3 verify_kn_flow.py
python3 verify_kn_dynamic_emission.py
python3 verify_kn_accretion_trace.py
python3 verify_kn_sampling.py
python3 verify_kn_image_orders.py
python3 compare_kn_emission.py
python3 render_kn_accretion.py
# Optional after a complete native cache: rebuild media without retracing rays
python3 assemble_kn_accretion.py
python3 verify_kn_accretion.py
python3 review_kn_accretion_video.py
python3 package_kn_accretion_release.py
```

The original orbit/flux checks remain relevant to the shared physical equations. New checks independently examine timelike normalization, Killing-vector energy and angular momentum, radial geodesic acceleration, the Schwarzschild plunge limit, continuous ISCO launch, horizon angular velocity, photon momentum signs, radial redshift, streamline quadrature, periodicity, bounded intensity, preserved exterior mean, multi-intersection paths and event-storage convergence.

Current numerical configuration, source hashes, transfer choices and generated media provenance are recorded in `output-sideratlas-kn-v120/render-metadata.json`; finished-file checks are in its `verification.json`. `assemble_kn_accretion.py` verifies and hashes every native radiance frame, then reuses the unchanged compositor and encoder without changing calculated radiance. Module reports identify the exact checked source files. Consult those measured reports for current results rather than applying historical check counts to this release. Cache directories contain resumable calculations and are excluded from Git. The complete download includes source, documentation and completed media; compiled libraries are rebuilt on the user's machine.

The independent sampling report `kn-sampling-diagnostics.json` compares the selected 2×2 quadrature with 3×3, 4×4 and a separately traced 12×12 critical-pixel reference. With the restored multiscale source at preview resolution, the largest full-frame mean RGB difference between 2×2 and 4×4 was 0.0359 on a 0–255 scale. Individual nearly unresolved subring pixels can have much larger relative intensity differences. This bounds the integrated presentation error for the tested cases; it does not claim that every arbitrarily thin photon subring is resolved at Instagram resolution.

These checks assess consistency, limiting cases, numerical convergence and delivery integrity within the declared model. They do not validate a charged plasma, the illustrative inner luminosity or a particular real black hole. Unresolved physical ingredients include MHD/MRI, pressure forces in the plunge, electromagnetic gas dynamics and discharge, finite thickness, a corona, scattering, polarization, returning-radiation heating and self-gravity.

## Sources and their roles

- **Newman et al. (1965)** — charged rotating metric underlying the modern geodesic equations. [DOI](https://doi.org/10.1063/1.1704351)
- **Kerr (1963)** — uncharged rotating metric and $q=0$ limit. [DOI](https://doi.org/10.1103/PhysRevLett.11.237)
- **Wang, Lee & Lin (2022)** — Kerr–Newman null and neutral timelike first integrals; an original numerical integrator implements these here. [arXiv:2208.11906](https://arxiv.org/abs/2208.11906)
- **Ko, Lee & Lin (2024)** — neutral Kerr–Newman inspiral/plunge geodesics; the elliptic-function solver is not imported. [arXiv:2309.13832](https://arxiv.org/abs/2309.13832)
- **Pugliese, Quevedo & Ruffini (2013)** — neutral equatorial orbit and stability reference. [arXiv:1303.6250](https://arxiv.org/abs/1303.6250)
- **Ivanov & Prodanov (2005)** — exact neutral circular energy and angular momentum, Eqs. 19–20; no pseudo-Newtonian potential is used. [arXiv:gr-qc/0504025](https://arxiv.org/abs/gr-qc/0504025)
- **Liu, Lee & Lin (2017)** — independent neutral ISCO-equation check. [arXiv:1706.05466](https://arxiv.org/abs/1706.05466)
- **Page & Thorne (1974)** — general exterior thin-disk flux integral and zero-stress brightness background. [DOI](https://doi.org/10.1086/152990)
- **Penna et al. (2013)** — motivation for the prescribed $\alpha(H/r)^2$ slow radial-drift scale; no viscous solution is imported. [DOI](https://doi.org/10.1093/mnras/sts185)
- **Mummery & Stone (2024)** — physical motivation for spiral transport inside the ISCO; their GRMHD, density and heating solutions are not implemented. [arXiv:2407.02164](https://arxiv.org/abs/2407.02164)
- **Schnittman & Bertschinger (2004)** — precedent for orbiting/shearing brightness structures. [arXiv:astro-ph/0309458](https://arxiv.org/abs/astro-ph/0309458)
- **Cárdenas-Avendaño, Lupsasca & Zhu (2023)** — time-dependent equatorial-source transfer and per-image retardation methodology; AART and its stochastic source solver are not used. [arXiv:2211.07469](https://arxiv.org/abs/2211.07469)
- **Gralla, Holz & Wald (2019)** — image-order interpretation, source dependence and higher-order visibility. [arXiv:1906.00873](https://arxiv.org/abs/1906.00873)
- **Vincent et al. (2022)** — source/opacity dependence of higher-order images; no spectral or thick-disk fit is applied. [DOI](https://doi.org/10.1051/0004-6361/202244339)
- **Zajaček et al. (2018)** — astrophysical charge context, not support for an ordinary $q=0.4$ accretion flow. [arXiv:1808.07327](https://arxiv.org/abs/1808.07327)
- **Bardeen, Press & Teukolsky (1972)** — inherited uncharged Kerr orbit/ISCO checks. [DOI](https://doi.org/10.1086/151796)
- **Gralla & Lupsasca (2020)** — inherited Kerr null-geodesic background and comparison limit. [arXiv:1910.12881](https://arxiv.org/abs/1910.12881)
- **James et al. (2015)** — inherited lensing/transfer background; DNGR is not used. [arXiv:1502.03808](https://arxiv.org/abs/1502.03808)
- **Luminet (1979)** — inherited Schwarzschild disk-image benchmark. [A&A 75, 228](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L)
- **Müller & Frauendiener (2012)** — inherited Schwarzschild propagation benchmark. [arXiv:1206.4259](https://arxiv.org/abs/1206.4259)
- **Bozza (2002)** — inherited strong-deflection benchmark. [arXiv:gr-qc/0208075](https://arxiv.org/abs/gr-qc/0208075)
- **NASA/Goddard/J. Schnittman (2019)** — viewpoint and visual reference only; no footage, textures, simulation data or source code is included. [SVS 13326](https://svs.gsfc.nasa.gov/13326/)
- **ICC sRGB** — display encoding reference, with Rec.709 video conversion documented separately. [Specification](https://registry.color.org/rgb-registry/files/sRGB.pdf)

Earlier thermal, atmosphere and diffraction references remain with the historical outputs that used them. They do not supply a thermal or spectral model for this bolometric release.
