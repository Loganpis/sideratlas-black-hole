# @sideratlas — THE BLACK HOLE

The current release is the original **Kerr–Newman model with a* = +0.9 and q = +0.4**, with an aligned, prograde neutral test disk.

- Main 32-second Reel: `output-sideratlas-kn/sideratlas-kerr-newman-reel.mp4`
- Monochrome companion: `output-sideratlas-kn/sideratlas-kerr-newman-intensity.mp4`
- Reusable eight-second cycles: the corresponding `-loop.mp4` files in the same folder
- Instagram description and all credits: `caption-kerr-newman.text` (`.txt` also supplied for compatibility)
- Equations, choices, limitations and reproduction: `README-KERR-NEWMAN.md`
- Encoding and posting guidance: `INSTAGRAM-EXPORT.md`
- Complete source and delivery bundle: `sideratlas-kerr-newman-source.zip`, available from the [v1.1.0 GitHub release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.1.0) or generated with `python3 package_kn_release.py`

The final movies are 1080 × 1920, 30 fps and silent. Each repeats exactly the same eight-second cycle four times. The viewpoint completes a full revolution per cycle, passing above and below the disk. The source keeps evolving forward, with the traced photon delays applied separately at every disk intersection. The repetitions are copied from the encoded cycle without another lossy encode.

This is a **hypothetical large-charge configuration**, not a model of the usual nearly neutral astrophysical black hole. Charge changes the spacetime, photon paths, circular orbits, ISCO and mean flux. The disk emitters are neutral test matter: the code does not simulate the electromagnetic response, discharge or fluid dynamics of accretion plasma. Orange is an illustrative brightness palette, not a visible thermal spectrum. Those distinctions are explained in the description file; post it with the video.

The enlarged disk now fills the portrait frame on a fixed diagonal, behind **@sideratlas**, **NO.001**, and **THE BLACK HOLE**. These are the only on-screen text. The account name begins 270 pixels from the top, the number at 277 pixels, and the title at 333 pixels. Latin Modern text has a subtle glyph shadow. There is no lower text block or reserved lower panel. The close crop can exclude the outer disk during some viewpoints, especially near the poles. Instagram overlays vary with viewing mode, caption expansion, and device.

The description retains the feature explanations, including “Spin drags spacetime around the black hole” and “Electric charge changes spacetime and light paths”, plus the model choices and all source credits. It adds “The code for this simulation is to be found at github.com/loganpis/sideratlas-black-hole”. The earlier layout with the horizontal disk and lower text remains available in [v1.0.0](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.0.0).

The physical model, assumptions and raw-radiance provenance are in `output-sideratlas-kn/render-metadata.json`. The final video checks are in `verification.json`; independent ray and emission reports are `kn-orbit-trace-diagnostics.json` and `kn-orbit-emission-diagnostics.json` in the source folder. `visual-review.json` records decoded-frame comparisons. These establish numerical consistency and export integrity for the stated model, not validation against a full charged-plasma simulation.

Run `render_kn_orbit.py` to reproduce the current release and `verify_kn_orbit.py` to check the finished files. Python dependencies are listed in `requirements.txt`; the new supplied C++ kernel requires clang++ or g++, and the export requires FFmpeg with libx264 and libzimg/zscale plus FFprobe. Fonts and licenses are bundled. See the README for complete commands and independent physical checks.

Earlier Schwarzschild, Kerr and optical experiments remain in the repository's historical output folders. Their source and completed media are included in the current release archive; duplicate historical archives are omitted. They are preserved historical work, not the current charged model. No NASA footage, textures, simulation data or source code is used.
