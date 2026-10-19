# Sideratlas · NO.001 — The Black Hole

A 32-second vertical black-hole animation for **@sideratlas**, with the code and scientific documentation needed to reproduce it.

The model is **Kerr–Newman: a* = 0.9, q = 0.4**, with an aligned neutral accretion flow and calculated photon paths. The disk fills the frame diagonally, behind the centered title. Each movie contains four identical eight-second loops.

**[Download v1.2.1](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.2.1)** · **[Instagram caption](caption-kerr-newman.text)** · **[Scientific model and papers](MODEL.md)**

<a href="https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.2.1"><img src="output-sideratlas-kn-v120/sideratlas-kerr-newman-color-cover.png" width="360" alt="Sideratlas black hole with a centered title over a diagonal accretion disk"></a>

## Watch and post

- [Color Reel](output-sideratlas-kn-v120/sideratlas-kerr-newman-reel.mp4) — 1080 × 1920, 30 fps, 32 seconds
- [Monochrome Reel](output-sideratlas-kn-v120/sideratlas-kerr-newman-intensity.mp4) — the same calculated radiance
- [Instagram description](caption-kerr-newman.text) — explanations, assumptions, every study credit, and the code link
- [Export and upload notes](INSTAGRAM-EXPORT.md)

The large electric charge is hypothetical. Gas motion and light propagation are calculated within the declared model; brightness fluctuations, slow outer drift, inner heating, and optical depth are prescribed. Orange represents bolometric brightness rather than a predicted visible spectrum. [MODEL.md](MODEL.md) explains the equations, limitations, and contribution of every cited study.

## Reproduce the animation

Use Python 3.12, NumPy, Pillow, a C++17 compiler (`clang++` or `g++`), and FFmpeg/FFprobe with `libx264` and `zscale`.

```sh
python3 -m pip install -r requirements.txt
python3 render_kn_accretion.py
python3 generate_kn_steady_reference.py
python3 assemble_kn_accretion.py
python3 verify_kn_accretion.py
python3 review_kn_accretion_video.py
python3 package_kn_accretion_release.py
```

The pipeline separates the physical calculation from the final titles and video export. If the native frame cache already exists, skip `render_kn_accretion.py`. Assembly verifies every saved frame and reuses its exact radiance. `--radiance-cache /path/to/cache` selects an existing cache for assembly or decoded-frame review. Raw caches and compiled libraries are excluded from Git.

## Find the code

| Files | Purpose |
| --- | --- |
| `render_kn_accretion.py` | Physical calculation and native frames |
| `kn_accretion_trace.py`, `kn_accretion_trace.cpp` | Photon trajectories, capture, multiple intersections, and delays |
| `kn_flow.py`, `kn_orbits.py` | Gas motion, circular orbits, and the ISCO plunge |
| `kn_dynamic_emission.py`, `kn_orbit_emission.py` | Emission and mean Page–Thorne brightness |
| `assemble_kn_accretion.py`, `title_overlay.py`, `instagram_export.py` | Saved-frame assembly, centered title, and final encoding |
| `verify_kn_*.py`, `video_io.py`, `video_validation.py` | Numerical and delivery checks |
| `review_kn_accretion_video.py` | Decoded-frame comparison with saved radiance |
| `compare_kn_emission.py`, `compare_kn_accretion.py` | Controlled tests of the current emission prescription |
| `output-sideratlas-kn-v120/` | Current finished movies, covers, metadata, and reports |
| `fonts/` | Bundled fonts and licenses |

Shared modules with earlier names (`render_kn_orbit.py`, `render.py`, `emission.py`, `nasa_emission.py`, `raytrace.py`, and the Kerr/Schwarzschild benchmark helpers) remain because the current calculation or its independent limit checks import them. Their frozen bytes preserve the provenance of the saved physical frames. Older animations, captions, standalone render drivers, and experimental documentation have been removed from the current tree.

See [MODEL.md](MODEL.md) for the science, [PUBLICATION.md](PUBLICATION.md) for reproducibility and checksums, and [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md) for font attribution.
