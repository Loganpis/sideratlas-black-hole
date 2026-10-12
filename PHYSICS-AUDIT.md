# Physical and numerical audit

The implementation is physically consistent with the stated **Schwarzschild, steady, opaque, geometrically thin, zero-torque disk** model. This is not a GRMHD calculation or a prediction for a particular observed object. Only the display orientation changes with movie phase. The spectrum is a local isotropic blackbody, with no atmosphere hardening, corona, returning radiation or photon-capture feedback on the disk temperature.

## Checks performed

- The exact Schwarzschild null-orbit equation is integrated with RK4. The first disk hit at 6–40 gravitational radii is the only emitting surface visible along each ray. Both surfaces can emit.
- The frequency shift combines gravity and circular Keplerian motion: `g = sqrt(1-3/r)/(1-lambda/r^(3/2))`, in units `G=M=c=1`. Specific intensity obeys Liouville transfer. For a blackbody the observed spectrum is exactly the Planck spectrum at `g*T_eff`; applying another `g^3` or `g^4` factor to the computed RGB would be incorrect.
- The Page–Thorne closed-form radial flux agrees with independent numerical quadrature to about `2e-8` relative error. Its maximum is near `9.55095 r_g`. The normalization in Müller & Frauendiener (2012), Eq. (11), agrees after converting their Schwarzschild-radius units to `r_g=GM/c^2`.
- Independent Planck integration reproduces Stefan–Boltzmann to about `1.3e-8`; the spectral Liouville identity agrees to floating-point precision.
- The spectral lookup was increased from 2048 to 4096 temperature samples. Its largest tested relative RGB-channel error over 2000–500000 K is now `2.43e-4`, below the stated 0.1% tolerance. The earlier 2048-point table narrowly failed that criterion at 2000 K.
- Existing ray-tracer diagnostics pass reflection symmetry, emission-radius limits, Doppler handedness, weak-field deflection and resolution refinement.
- A targeted strong-field audit compared the orbit lookup to directly integrated rays around `b_c=sqrt(27) r_g`, at 50°, 66° and 82° inclination. The largest tested 99th-percentile radius discrepancy was `6.96e-5`. The estimated absolute bolometric error integrated over this critical annulus is less than `2.7e-8` of the full sampled image flux.
- Additional fifth and sixth disk-plane crossings occupy less than `0.0032` output pixel radially at 1080 pixels across a `32 r_g` field. Their estimated omitted bolometric flux is approximately `3–4e-7` of the full image. Thus four crossings are sufficient at this delivery resolution; this does not guarantee adequacy for an arbitrarily magnified photon-ring study.
- The escape angle on supercritical rays approaches the known logarithmic strong-deflection expression in Bozza (2002), Eq. (50). The subcritical comparison ray is captured. This is an independent check, not a method used to render the image.

Detailed numerical results are in `raytrace-diagnostics.json`, `emission-diagnostics.json`, and `physics-audit.json`. The targeted tests are spot checks and quadrature estimates, not universal error bounds. Pixel supersampling and the final video codec require separate rendering checks.

## Attribution by role

1. **Page & Thorne (1974)** — zero-torque relativistic accretion flux. [DOI](https://doi.org/10.1086/152990).
2. **Luminet (1979)** — thin-disk black-hole imaging and relativistic frequency/flux transfer. [Paper](https://articles.adsabs.harvard.edu/pdf/1979A%26A....75..228L).
3. **Müller & Frauendiener (2012)** — thin-disk geometry and the Schwarzschild closed-form flux used for cross-checking. [Paper](https://arxiv.org/abs/1206.4259).
4. **James et al. (2015)** — credited relativistic-appearance reference. No DNGR implementation is included. [Paper](https://arxiv.org/abs/1502.03808).
5. **Bozza (2002)** — independent strong-deflection validation only, Eq. (50). [Paper](https://arxiv.org/abs/gr-qc/0208075), [DOI](https://doi.org/10.1103/PhysRevD.66.103001).
6. **CIE (2019)** — official 1931 2° color-matching functions, 360–830 nm at 1 nm intervals. [Dataset](https://doi.org/10.25039/CIE.DS.xvudnb9b).
7. **ICC sRGB interpretation** — XYZ to linear-sRGB transform. [Specification](https://registry.color.org/rgb-registry/files/sRGB.pdf).

The implementation does not use AART, DNGR, a Kerr solver, a GRMHD snapshot, fabricated turbulence or star imagery. The event-horizon radius is `2 r_g`; it should not be confused with the apparent dark region or the critical impact radius `sqrt(27) r_g`.
