# @sideratlas — THE BLACK HOLE — v1.2.0

The current release depicts **Kerr–Newman: a* = +0.9, q = +0.4**, with an aligned, prograde neutral test flow. It adds restrained accreting structure, calculated gas trajectories inside the ISCO and multiple lensed images through a thin inner emitting layer.

- [Main 32-second Reel](output-sideratlas-kn-v120/sideratlas-kerr-newman-reel.mp4)
- [Monochrome companion](output-sideratlas-kn-v120/sideratlas-kerr-newman-intensity.mp4)
- Reusable eight-second cycles: the corresponding `-loop.mp4` files in `output-sideratlas-kn-v120/`
- [Instagram description and all credits (.text)](caption-kerr-newman.text), with a `.txt` compatibility copy
- [Equations, choices, limitations and study roles](README-V1.2.0.md)
- [Encoding and posting guidance](INSTAGRAM-EXPORT.md)
- [Complete v1.2.0 GitHub download](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.2.0)

The final movies are 1080 × 1920 at 30 fps and silent. Each repeats the same encoded eight-second cycle exactly four times. The viewpoint completes a full revolution per cycle, passing above and below the disk, while source time advances through 800 GM/c³. This is a repeating source/viewpoint presentation, not 32 seconds of distinct simulation evolution.

The disk fills the portrait frame on a fixed diagonal behind **@sideratlas**, **NO.001** and **THE BLACK HOLE**. These are the only on-screen words. There is no lower text block. The black sky, mathematical typography, fixed exposure and restrained orange palette are retained. The close crop can exclude part of the outer disk near the poles, and Instagram controls can cover parts of the image.

Post the supplied description with the video. It contains the explanations, parameter choices, scientific qualifications and all study credits, plus “The code for this simulation is to be found at github.com/loganpis/sideratlas-black-hole”. Orange is illustrative bolometric brightness, not a predicted thermal or visible spectrum.

**What is physical:** the Kerr–Newman geometry, photon trajectories and image orders, neutral timelike plunge, relativistic frequency shifts including radial gas motion, capture and per-image propagation delays. The event horizon itself remains non-emitting.

**What is prescribed:** the small outer radial drift, smooth finite-lived brightness irregularities and spiral modes, weak inner heating, gray optical depth and repeating source births. The existing Page–Thorne mean profile is preserved outside the ISCO as a thin-disk background; the model does not solve a self-consistent viscous flow. Nested images emerge from the paths and transfer, with no artificially drawn stripes. Their visible width and brightness depend on opacity, emission and final pixel scale.

**The large charge is hypothetical.** Emitters are neutral test matter in a fixed charged spacetime. The code does not simulate electromagnetic gas forces, MHD/MRI, charge discharge, disk thickness or thermodynamics. The documents clearly distinguish these approximations from calculated relativistic effects.

Run `render_kn_accretion.py` for the current model. After a native cache is complete, `assemble_kn_accretion.py` can verify its frames and rebuild the delivery without changing the radiance. Full setup and independent numerical checks are in [README.md](README.md). Current configuration and provenance are in `output-sideratlas-kn-v120/render-metadata.json`, with finished-video results in its `verification.json`. Numerical consistency checks do not validate an accreting charged plasma or the illustrative inner luminosity.

The previous circular-emitter model remains in `output-sideratlas-kn/`, with `render_kn_orbit.py` and [v1.1.0](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.1.0). The [v1.0.0 release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.0.0) retains the horizontal disk and lower text. Other Schwarzschild, Kerr and optical experiments remain in historical folders. No NASA footage, textures, simulation data or source code is used.
