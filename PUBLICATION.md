# Reproducibility and publication

The current release is **v1.2.1**, prepared on 2026-10-19 for `Loganpis/sideratlas-black-hole`.

The title is horizontally centered over the same v1.2 physical simulation. The account name and number retain their positions. Caption headings are **ASSUMPTIONS** and **SOURCES**.

## Physical frames and presentation

The physical source SHA-256 remains `f638dad36408167811a34e5929286e0051222b54879a928646d14a61f158b66f`. The native radiance signature remains `af43803e4e58364866f863a7e5d0026016ed0037338905b26d4ea37d26539558`.

`assemble_kn_accretion.py` checks all 240 native frames against their original producer manifest before adding the centered title and encoding the eight-second cycle. The 32-second movie repeats that encoded cycle four times by packet copy. No photon trajectories or movie radiance are recalculated for this title change. A separate steady-emission validation still is generated with the current title.

`title_overlay.py` has its own checksum in the presentation metadata. The delivery verifier binds the caption, numerical reports, assembly utility, title module, steady-reference producer, and completed media. The decoded-frame reviewer compares the exported movie with the exact cached radiance and current compositor.

## Repository contents

The current tree contains this simulation, its required shared modules and numerical benchmarks, the current output set, scientific documentation, and font licenses. Earlier standalone animations and experimental documents are removed from the current tree. Git history and existing release tags preserve earlier publications.

Some shared modules retain their earlier names and original bytes because the current calculation and its independent limit checks use them. The active pipeline and each file's role are listed in [README.md](README.md).

## Verify a download

`output-sideratlas-kn-v120/release-files.json` lists file sizes and SHA-256 checksums for archive contents. `delivery-manifest.json` records the archive and video hashes. Raw caches, compiled libraries, process records, and temporary previews are excluded.

Tests establish numerical consistency and export integrity for the assumptions in [MODEL.md](MODEL.md). The reports identify their exact source hashes and measured results.
