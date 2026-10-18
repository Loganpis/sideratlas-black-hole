# Publication provenance

Prepared for `Loganpis/sideratlas-black-hole` on 2026-10-12; updated for version 1.1.0 on 2026-10-13 and version 1.2.0 on 2026-10-18.

The public project contains the current Kerr–Newman source, finished videos, stills, caption, scientific documentation, and validation reports, together with the earlier completed visualization versions. Fonts and color-matching data retain their upstream attribution and licenses.

Temporary previews, compiled libraries, Python bytecode, process logs, duplicate historical archives, and several gigabytes of regenerable numerical caches are excluded. The current release archive bundles the published source and completed media. GitHub also provides a source snapshot for each release tag.

## Version 1.2.0: accretion and geodesic plunge

The new entry point `render_kn_accretion.py` and its new flow, emission and multiple-event tracing modules extend the model while preserving the previous renderer and completed media. Current outputs are in `output-sideratlas-kn-v120/`. The previous circular-emitter model remains in `output-sideratlas-kn/` and its tagged releases.

New calculated quantities include the neutral timelike plunge, full radial emitter contribution to frequency shift, ordered plunging-layer intersections and selective independent subpixel photon paths. The outer radial drift, finite-lived spiral/irregular emission, weak inner heating, gray optical depth and periodic source forcing are explicitly prescribed. The exterior Page–Thorne mean profile is retained. No GRMHD, plasma thermodynamics or charge-discharge calculation is claimed.

The [model document](README-V1.2.0.md), module reports, sampling convergence report and current delivery verification distinguish these choices from calculated effects. `assemble_kn_accretion.py` records every native-frame checksum and reuses the production compositor and encoder without altering radiance. Final media and public archive checksums are recorded in the current output folder's manifests. Runtime logs, compiled libraries, numerical caches and temporary previews remain excluded.

The new flow, emission, multiple-event tracing, sampling, image-order and source-comparison reports pass 207 checks against their recorded source files. A separate `shared-physics-regression.json` report records 131 passing checks of the unchanged shared null tracing, neutral circular orbits and Page–Thorne emission. These counts describe consistency tests within the stated model; they do not validate the prescribed emissivity or a charged plasma.

The first v1.2 emission experiment was excessively smooth because it replaced the original narrow patch supports with only three weak broad modes. Before publication it was replaced with positive, analytically normalized multiscale packets and a coherent spiral. The Page–Thorne azimuthal mean, flow, geodesics, opacity, sampling, shutter and display curve were retained. Source and image-order comparison reports document the distinction between emissivity texture and calculated multiple images; the comparison contact sheet uses identical current paths for all source prescriptions. Temporary superseded media and their radiance caches remain excluded from publication.

## Version 1.1.0: full-frame diagonal composition

The requested presentation update enlarges the simulation to fill the 1080 × 1920 frame and uses a fixed −arctan(16/9) ≈ −60.642° screen rotation. The portrait viewport has a horizontal field of view of 28 r_g. An unrolled 2700 × 2250 detector covers 72 × 60 r_g; floating-point radiance is bilinearly sampled into a 2× portrait viewport and box-averaged before tone mapping. The close framing can crop outer-disk emission, particularly near the poles. The disk remains equatorial and aligned with the black-hole spin; screen rotation does not introduce physical disk tilt.

Only @sideratlas, NO.001, and THE BLACK HOLE remain in the video, with the simulation behind the lettering. The entire lower text block, model disclosures, and bibliography move to `caption-kerr-newman.text`; the `.txt` copy is retained for compatibility. The description includes the repository link requested for Instagram.

The Kerr–Newman parameters, geodesic equations, orbital motion, emissivity prescription, photon delays, eight-second cycle, and four-repeat delivery remain the same. The larger source detector and new viewport are presentation and sampling changes. Current generated media, metadata, and delivery reports identify this composition; their hashes are recorded in the release manifests. Numerical and delivery reports should be interpreted with the exact source revision they identify.

The completed v1.1.0 delivery passed 84 renderer, composition, provenance, and video checks. Both 32-second files decode all 960 frames and contain four pixel-identical cycles. In twelve color-managed decoded samples, the largest full-frame mean absolute RGB error was 0.699 on a 0–255 scale; the largest mean error over lit disk pixels outside the title band was 0.992. The only overlay occupies the title region, with bounds (65, 268, 965, 393); no lower overlay remains. These are local export checks, not measurements after Instagram processing.

The [v1.1.0 release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.1.0) contains the diagonal composition. The previous tagged [v1.0.0 release](https://github.com/Loganpis/sideratlas-black-hole/releases/tag/v1.0.0) retains the smaller horizontal disk and lower text block, together with that delivery's source, media, and reports.

## Version 1.0.0: portable report copies

Local absolute project paths were replaced by repository-relative paths in these public report copies:

- `output/verification.json`
- `output-kerr/verification.json`
- `output-optical/verification.json`
- `output-schwarzschild/verification.json`
- `output-sideratlas/verification.json`
- `output-sideratlas/provenance/prior-verification.json`
- `output-sideratlas-kn/verification.json`

For the initial publication, measured results, numerical configurations, source hashes, and media bytes were preserved. Report file hashes changed because path strings changed. The release file manifest and archive checksum were rebuilt from the public copies. The numerical tests were performed for the original delivered render; that publication step did not represent a new physical render or rerun of every historical test. The current Kerr–Newman delivery reports are updated separately for version 1.1.0; the list above records the historical path normalization.

The **v1.0.0 renderer's historical source SHA-256** is `3103914c35fa2e799ee5eaba5c72b11d052c6ee0c557064cb7302be413c196b9`. Its input source and font bytes were unchanged by the initial publication. The **v1.0.0 orange 32-second video's historical SHA-256** is `fc5ae0ced5d396bf94fff9ab0d3882bbb550a05a9691176cc98926523ef352f1`. These hashes do not identify the version 1.1.0 composition.

The root GitHub guide, these publication notes, third-party notices, ignore rules, and release packaging were added or adapted for the initial publication. Those documentation and packaging changes did not alter the physical model or then-completed animation.

Earlier reports are retained as historical records and may identify source revisions preceding later shared-module changes. Use `output-sideratlas-kn-v120/` reports for version 1.2.0, and the earlier output folders or tagged releases for historical models.
