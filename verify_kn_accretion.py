#!/usr/bin/env python3
"""Verify the v1.2 accretion renderer, numerical provenance and delivery.

These checks validate numerical consistency of the declared relativistic
geometry, kinematic source and thin-slab transfer. They do not establish a
GRMHD solution, plasma thermodynamics, or agreement with an observed disk.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from render_kn_accretion import (ROOT, OUTPUT_DIR, VIDEO_NAMES, FPS, DURATION,
    LOOP_DURATION, REPEATS, LOOP_TG, SPIN, CHARGE, FOV, INCLINATION,
    KerrNewmanAccretionRenderer, PreparedRayEmission, transfer_weights,
    trace_resolved, source_hash, PLUNGE_OPTICAL_DEPTH)
from kn_accretion_trace import KerrNewmanAccretionTracer
from kn_flow import NeutralAccretionFlow
from kn_dynamic_emission import DynamicDiskEmission
from instagram_export import metadata as export_metadata
from title_overlay import provenance as title_provenance
from video_io import check
from video_validation import (sha256, finish, display_checks, composition_checks,
                             output_assets, video_checks)

ASSEMBLY_ASSETS = (list(VIDEO_NAMES.values()) +
    ['sideratlas-kerr-newman-reel-loop.mp4',
     'sideratlas-kerr-newman-intensity-loop.mp4',
     'sideratlas-kerr-newman-color-cover.png',
     'sideratlas-kerr-newman-mono-cover.png', 'contact-sheet.jpg',
     'steady-reference.png', 'steady-reference-manifest.json',
     'native-radiance-manifest.json',
     'view-diagnostics.json', 'render-metadata.json'])

NUMERICAL_REPORTS = {
    "kn-flow-diagnostics.json": 70,
    "kn-accretion-trace-diagnostics.json": 51,
    "kn-dynamic-emission-diagnostics.json": 35,
    "kn-sampling-diagnostics.json": 27,
    "kn-image-order-diagnostics.json": 9,
    "kn-emission-comparison-diagnostics.json": 15,
}
RELEASE_INPUTS = ('caption-kerr-newman.text', 'caption-kerr-newman.txt',
                  *NUMERICAL_REPORTS)


def numerical_reports():
    """Require saved numerical measurements to match every current input."""
    result = {"scope": "Saved numerical diagnostics; every recorded source checksum is checked",
              "reports": {}}
    for name, minimum_checks in NUMERICAL_REPORTS.items():
        path = ROOT/name
        report = json.loads(path.read_text())
        matches = {n: sha256(ROOT/n) == digest for n, digest in report["source_sha256"].items()}
        count = len(report["checks"])
        result["reports"][name] = {"sha256": sha256(path), "checks": count,
                                    "source_matches": matches}
        passed = report.get("passed", report.get("status") == "passed")
        check(result, name+"_passed_with_current_sources", passed and count >= minimum_checks
              and all(matches.values()) and all(c["passed"] for c in report["checks"]))
    return finish(result)


def transfer_checks():
    """Independent closed-form transfer and near-horizon emitter tests."""
    result = {"scope": "Closed-form foreground transfer and outgoing-photon horizon limit"}

    class KnownShifts:
        def __init__(self, shifts):
            self.shifts = np.asarray(shifts)
        def redshift(self, *args):
            return self.shifts.copy()

    tau = PLUNGE_OPTICAL_DEPTH
    radius = np.array([1.4, 3.])
    g = np.array([.5, .8])
    mu = np.array([.4, .7])
    events = {"pixel": np.array([0, 0]), "slot": np.array([0, 1]),
        "r": radius, "lambda": np.zeros(2), "kr_over_E": np.zeros(2),
        "eta": (mu*radius/g)**2, "hit_kind": np.array([2, 1]),
        "hit_layer": np.array([0, 2])}
    weights, measured_g = transfer_weights(events, KnownShifts(g), 1, tau)
    trans = math.exp(-tau/mu[0])
    expected = np.array([(1-trans)/(1-math.exp(-tau))*g[0]**4,
                         trans*g[1]**4])
    error = float(np.max(abs(weights-expected)))
    check(result, "two_layer_slab_then_opaque_closed_form", error < 1e-15,
          maximum_absolute_error=error)
    check(result, "each_emitter_has_own_g_to_fourth", np.array_equal(g, measured_g)
          and not np.isclose(weights[1], trans*g[0]**4))
    opaque = {**events, "hit_kind": np.array([1, 2])}
    opaque_weights, _ = transfer_weights(opaque, KnownShifts(g), 1, tau)
    check(result, "opaque_foreground_blocks_all_later_events",
          opaque_weights[0] == g[0]**4 and opaque_weights[1] == 0)
    alternate_orders = {**events, "hit_layer": np.array([7, 19])}
    alternate_weights, _ = transfer_weights(alternate_orders, KnownShifts(g), 1, tau)
    check(result, "no_image_order_brightness_boost", np.array_equal(weights, alternate_weights))
    grazing = {**events, "eta": np.zeros(2)}
    grazing_weights, _ = transfer_weights(grazing, KnownShifts(g), 1, tau)
    expected_grazing = g[0]**4/(1-math.exp(-tau))
    check(result, "grazing_thin_slab_saturates_finitely",
          np.isfinite(grazing_weights).all() and abs(grazing_weights[0]-expected_grazing) < 1e-15
          and grazing_weights[1] < 1e-30)

    flow = NeutralAccretionFlow()
    source = DynamicDiskEmission(flow=flow)
    radii = flow.horizon_rg+np.array([.01, .001, .0001])
    delta = (radii-flow.horizon_rg)*(radii-(2-flow.horizon_rg))
    lam, eta = np.zeros(3), np.full(3, .5)
    potential = (radii*radii+SPIN*SPIN)**2-delta*(eta+SPIN*SPIN)
    kr = np.sqrt(potential)/delta
    shifts = flow.redshift(radii, lam, kr)
    intrinsic = source.radial_profile(radii)
    observed = intrinsic*shifts**4
    check(result, "near_horizon_intrinsic_emission_is_finite_not_faded_by_hand",
          np.isfinite(intrinsic).all() and np.all(np.diff(intrinsic) > 0)
          and intrinsic[-1] > .99*source.parameters.plunge_peak_fraction,
          relative_emission=intrinsic.tolist())
    check(result, "outgoing_plunge_g4_drives_near_horizon_fade",
          np.all(np.diff(shifts) < 0) and np.all(np.diff(observed) < 0)
          and observed[-1]/observed[0] < 1e-6,
          frequency_shifts=shifts.tolist(), observed_relative_intensity=observed.tolist())
    velocity = flow.four_velocity(radii)
    no_radial_term = 1/(velocity["ut"]-lam*velocity["uphi"])
    check(result, "radial_momentum_term_materially_changes_plunging_g",
          float(np.max(abs(shifts/no_radial_term-1))) > .1,
          largest_relative_difference=float(np.max(abs(shifts/no_radial_term-1))))
    return finish(result)


def actual_renderer_checks():
    renderer = KerrNewmanAccretionRenderer(270, 1.5, refinement=2, temporal_samples=3)
    result = {"scope": "Fresh 270x480 renderer with 1.5x detector, 2x2 selected subpixels and three source-shutter samples",
              "geometry": renderer.geometry, "view_diagnostics": {}}
    check(result, "chosen_spin_charge_and_continuous_inner_flow",
          renderer.spin == .9 and renderer.charge == .4
          and renderer.flow.horizon_rg < renderer.tracer.rmin < renderer.source.rin
          and renderer.flow.inner_radius_rg <= renderer.tracer.rmin)
    first = renderer.radiance(0)
    endpoint = renderer.radiance(1)
    check(result, "raw_radiance_loop_endpoints_equal", np.array_equal(first, endpoint))
    frames0, frames1 = renderer.frames(0), renderer.frames(1)
    for mode in VIDEO_NAMES:
        check(result, mode+"_uncompressed_loop_endpoints_equal",
              np.array_equal(np.asarray(frames0[mode]), np.asarray(frames1[mode])))
    initial = renderer.coarse.prepared.evaluate(0)
    periodic = renderer.coarse.prepared.evaluate(LOOP_TG)
    error = float(np.max(abs(initial-periodic)))
    check(result, "per_ray_retarded_source_is_periodic", error < 3e-12,
          maximum_error=error)
    events = renderer.coarse.events
    no_delays = renderer.source.prepare(events["r"], events["phi"], 0).evaluate(0)
    delay_effect = float(np.sqrt(np.mean((initial-no_delays)**2)))
    check(result, "per_image_delays_change_variable_emission", delay_effect > 1e-5,
          rms_emitted_intensity_difference=delay_effect)
    source_times = []
    original_evaluator = renderer.coarse.prepared.evaluate
    def record_evaluate(t, structured=True):
        source_times.append(float(t))
        return original_evaluator(t, structured=structured)
    renderer.coarse.prepared.evaluate = record_evaluate
    exposed = renderer.radiance(0)
    renderer.coarse.prepared.evaluate = original_evaluator
    half_width = LOOP_TG/(FPS*LOOP_DURATION)*.25
    offset = math.sqrt(3/5)*half_width
    check(result, "three_point_finite_source_exposure_uses_correct_times",
          np.allclose(source_times, [-offset, 0, offset], rtol=0, atol=1e-14),
          sampled_times_tg=source_times)
    explicit = ((5/18)*renderer._radiance_at_time(-offset)+(4/9)*renderer._radiance_at_time(0)
                +(5/18)*renderer._radiance_at_time(offset)).astype(np.float32)
    check(result, "shutter_averages_linear_radiance_before_display",
          np.array_equal(exposed, explicit))
    phases = [0., .1, .25, .375, .5, .75,
              INCLINATION/360, (INCLINATION+180)/360, (INCLINATION+90)/360]
    maximum_order, plunge_count = -1, 0
    for phase in phases:
        values = renderer.radiance(phase)
        diag = renderer.view_diagnostics
        result["view_diagnostics"][str(phase)] = diag
        maximum_order = max(maximum_order, diag["maximum_image_order"],
                            diag.get("refined", {}).get("maximum_image_order", -1))
        plunge_count += diag["plunging_events"]
        ev = renderer.coarse.events
        check(result, f"phase_{phase:g}_finite_radiance_and_exterior_emission",
              values.shape == (renderer.tracer.height, renderer.tracer.width)
              and np.isfinite(values).all() and np.all(values >= 0) and values.sum() > 0
              and np.all(ev["r"] > renderer.flow.horizon_rg)
              and np.all(ev["r"] <= 28+2e-6)
              and np.all(renderer.coarse.g > 0))
        check(result, f"phase_{phase:g}_resolved_independent_subpixel_rays",
              diag["maximum_constraint_residual"] < 2e-7
              and diag["refined_parent_pixels"] > 0
              and diag["refined"]["rays"] == 4*diag["refined_parent_pixels"]
              and diag["refined"]["maximum_constraint_residual"] < 2e-7)
        if phase in (0., .1, .25):
            nodes, weights = np.polynomial.legendre.leggauss(9)
            reference = np.zeros_like(values, dtype=np.float64)
            for node, weight in zip(nodes, weights):
                reference += .5*weight*renderer._radiance_at_time(
                    phase*LOOP_TG+float(node)*half_width)
            relative_l1 = float(np.abs(values-reference).sum()/reference.sum())
            rgb_error = float(np.abs(np.asarray(renderer.compose(values, 'color'), dtype=float)
                -np.asarray(renderer.compose(reference, 'color'), dtype=float)).mean())
            check(result, f"phase_{phase:g}_three_to_nine_point_source_shutter_convergence",
                  relative_l1 < .001 and rgb_error < .05,
                  radiance_l1_relative_to_total=relative_l1, rgb_mae_255=rgb_error,
                  scope="Source-only exposure at fixed central photon geometry; excludes camera-motion blur")
    check(result, "natural_image_orders_at_least_second_order", maximum_order >= 2,
          maximum_order=maximum_order)
    check(result, "geodesic_plunge_source_is_actually_intersected", plunge_count > 0,
          sampled_plunging_intersections=plunge_count)
    for phase in (0., .1):
        reused = renderer.antipodal_radiance(phase)
        fresh = renderer.radiance(phase+.5)
        error = float(np.max(abs(reused-fresh)))
        check(result, f"antipodal_new_retarded_emission_matches_fresh_{phase:g}",
              np.allclose(reused, fresh, atol=2e-6, rtol=5e-5), maximum_absolute_error=error)

    # Capture is a ray termination, not a command to erase every foreground
    # intersection. Conversely the exactly radial polar ray sees no gas and
    # has zero radiance: the event horizon itself is never a source surface.
    maps = renderer.tracer.trace(INCLINATION)
    captured_foreground = maps["shadow"] & maps["mask"].any(axis=-1)
    renderer.set_view(0)
    raw = renderer.coarse.evaluate(0).reshape(renderer.tracer.height, renderer.tracer.width)
    check(result, "captured_rays_retain_luminous_foreground_gas",
          captured_foreground.any() and np.any(raw[captured_foreground] > 0),
          captured_foreground_pixels=int(captured_foreground.sum()),
          luminous_captured_pixels=int(np.count_nonzero(raw[captured_foreground] > 0)))
    center = renderer.tracer.trace_samples(np.array([0.]), np.array([0.]), 0.)
    from render_kn_accretion import compact_hits
    center_emission = PreparedRayEmission(compact_hits(center), renderer.source, renderer.flow, 1)
    check(result, "central_polar_horizon_ray_black_without_emissive_surface",
          center["shadow"][0] and not center["mask"].any()
          and center_emission.evaluate(0)[0] == 0)

    # Force the actual production retry branch on a ray with multiple thin
    # intersections. Compare to an independently requested high-capacity ray.
    multiple = np.flatnonzero(maps["mask"].sum(axis=-1).ravel() > 1)
    if multiple.size:
        index = multiple[0]
        alpha, beta = renderer.tracer.x.ravel()[index:index+1], renderer.tracer.y.ravel()[index:index+1]
        limited = KerrNewmanAccretionTracer(8, 8, max_hits=1)
        recovered, stats = trace_resolved(limited, alpha, beta, INCLINATION)
        reference = KerrNewmanAccretionTracer(8, 8, max_hits=16, angle_step=.006,
                                             max_angle=40*math.pi).trace_samples(alpha, beta, INCLINATION)
        truth = compact_hits(reference)
        check(result, "overflow_retraces_without_dropping_emitting_events",
              stats["retried_rays"] == 1 and len(recovered["r"]) > 1
              and all(np.array_equal(recovered[k], truth[k]) for k in
                      ("r", "phi", "delay_tg", "hit_layer", "hit_kind")),
              recovered_events=len(recovered["r"]), retry_stats=stats)
    else:
        check(result, "overflow_fixture_has_multiple_intersections", False)
    return finish(result)


def metadata_checks(output):
    d = json.loads((output/"render-metadata.json").read_text())
    result = {"scope": "Current renderer metadata and public scientific disclosures"}
    check(result, "recorded_source_matches_current_code", d["source_sha256"] == source_hash())
    check(result, "v120_delivery_dimensions_and_four_cycles",
          (d["version"], d["width"], d["height"], d["fps"], d["frames"], d["duration_seconds"],
           d["loop_duration_seconds"], d["repetitions"], d["unique_frames"])
          == ("1.2.0", 1080, 1920, 30, 960, 32, 8, 4, 240))
    g, spin, sample = d["geometry"], d["spin_configuration"], d["sampling"]
    check(result, "kerr_newman_spin_charge_and_exterior_plunge_boundary",
          spin["a_star"] == .9 and spin["charge_q"] == .4 and g["spin"] == .9
          and g["charge_q"] == .4 and g["rout_rg"] == 28
          and g["outer_horizon_rg"] < g["plunge_inner_cutoff_rg"] < g["rin_rg"]
          and g["plunge_inner_cutoff_rg"] < g["minimum_exterior_photon_turn_rg"])
    check(result, "full_orbit_static_observer_assumption",
          d["view_angle_degrees"] == [-277, 83] and d["camera_orbit"]["angular_span_deg"] == -360
          and "static observers" in d["camera_orbit"]["observer"])
    expected_source = DynamicDiskEmission().metadata()
    check(result, "saved_emission_matches_current_model", d["emission"] == expected_source)
    check(result, "saved_flow_matches_current_model", d["flow"] == NeutralAccretionFlow().metadata())
    check(result, "illustrative_source_opacity_charge_and_palette_disclosed",
          not d["emission"]["plunge_emissivity_is_calculated_from_thermodynamics"]
          and d["radiative_transfer"]["opacity_is_prescribed"]
          and "hypothetical" in spin["choice"] and "not a visible spectrum" in d["display"]["color"])
    check(result, "independent_subpixels_and_source_only_exposure_recorded",
          sample["refinement_grid"] == 2 and sample["temporal_samples"] == 3
          and "independently traced" in sample["refinement"]
          and "fixed central-view" in sample["shutter"])
    check(result, "export_contract_unchanged", d["encoding"] == export_metadata())
    layout = d["layout"]
    check(result, "title_only_full_frame_portrait_identity",
          layout["brand"] == "@sideratlas" and layout["lower_text_in_caption"]
          and layout["full_frame"] and (layout["hero_x"], layout["hero_y"],
              layout["hero_width"], layout["hero_height"]) == (0, 0, 1080, 1920))
    from assemble_kn_accretion import compositor
    renderer = compositor(d)
    title_box = renderer.overlays['color'].crop((0,320,1080,450)).getbbox()
    check(result, "title_is_centered_on_the_portrait",
          title_box is not None and abs((title_box[0]+title_box[2])/2-540) <= 1
          and layout.get('title_alignment') == 'horizontal center'
          and layout.get('title_center_x') == 540
          and d.get('presentation') == title_provenance(),
          title_bounds_in_band=list(title_box) if title_box else None)
    roll = math.radians(layout["display_roll_degrees"])
    fov = layout["viewport_fov_horizontal_rg"]
    half_w, half_h = fov/2, fov/2*d["height"]/d["width"]
    ex = abs(math.cos(roll))*half_w+abs(math.sin(roll))*half_h
    ey = abs(math.sin(roll))*half_w+abs(math.cos(roll))*half_h
    check(result, "rigid_diagonal_viewport_covered_by_traced_domain",
          fov == 28 and abs(roll+math.atan2(d["height"], d["width"])) < 1e-12
          and ex < g["fov_horizontal_rg"]/2 and ey < g["fov_vertical_rg"]/2)
    density = (g["width"]/g["fov_horizontal_rg"])/(d["width"]/fov)
    check(result, "native_radiance_density_near_output_and_extra_critical_sampling",
          density >= .97 and sample["refinement_grid"] >= 2,
          native_samples_per_output_pixel=density)
    check(result, "complete_render_and_no_partial_video_files",
          d.get("elapsed_seconds", 0) > 0 and not list(output.glob("*.partial.mp4")))
    assembly = d.get("assembly", {})
    check(result, "cached_assembly_provenance_matches_current_utility",
          assembly.get("script") == "assemble_kn_accretion.py"
          and assembly.get("script_sha256") == sha256(ROOT/"assemble_kn_accretion.py")
          and assembly.get("radiance_unchanged") is True
          and "excludes prior geodesic tracing" in d.get("elapsed_scope", ""))
    caption = (ROOT/"caption-kerr-newman.text").read_text()
    check(result, "caption_within_instagram_limit_and_matching_txt",
          len(caption.encode("utf-16-le"))//2 <= 2200
          and caption == (ROOT/"caption-kerr-newman.txt").read_text(), characters=len(caption))
    check(result, "caption_includes_requested_code_sentence_and_model_disclosure",
          "The code for this simulation is to be found at https://github.com/loganpis/sideratlas-black-hole" in caption
          and "not GRMHD" in caption and "hypothetical" in caption)
    check(result, "caption_uses_requested_headings",
          '\nASSUMPTIONS\n' in caption and '\nSOURCES\n' in caption
          and 'ILLUSTRATIVE CHOICES' not in caption and 'SOURCES AND CONTEXT' not in caption)
    return finish(result)


def view_diagnostic_checks(output):
    records = json.loads((output/"view-diagnostics.json").read_text())
    result = {"scope": "Recorded independent production central views and critical-zone subrays"}
    check(result, "all_120_independent_views_recorded_once", isinstance(records, list)
          and len(records) == 120 and [r["frame"] for r in records] == list(range(120)))
    if not isinstance(records, list):
        return finish(result)
    maximum_order = max(r["maximum_image_order"] for r in records)
    maximum_residual = max(max(r["maximum_constraint_residual"],
        r.get("refined", {}).get("maximum_constraint_residual", 0)) for r in records)
    check(result, "production_constraints_within_verified_kernel_tolerance",
          maximum_residual < 2e-7, maximum_residual=maximum_residual)
    check(result, "multiple_image_orders_recorded_from_actual_events", maximum_order >= 2,
          maximum_image_order=maximum_order)
    check(result, "critical_and_plunging_pixels_have_independent_four_ray_refinement",
          all(r["refinement_grid"] == 2 and r["refined_parent_pixels"] > 0
              and r["refined"]["rays"] == 4*r["refined_parent_pixels"] for r in records),
          total_refined_parent_pixels=sum(r["refined_parent_pixels"] for r in records))
    check(result, "every_view_intersects_luminous_plunging_region",
          all(r["plunging_events"] > 0 for r in records),
          total_plunging_events=sum(r["plunging_events"]+r["refined"]["plunging_events"] for r in records))
    result["retried_rays"] = sum(r["retried_rays"]+r["refined"]["retried_rays"] for r in records)
    return finish(result)


def assembly_artifact_checks(output):
    """Bind every cache-to-delivery artifact to this verification report."""
    result = {"scope": "Complete transactional assembly artifact set", "files": {}}
    missing = [name for name in ASSEMBLY_ASSETS if not (output/name).is_file()]
    check(result, "all_assembly_artifacts_present", not missing, missing=missing)
    if missing:
        return finish(result)
    result["files"] = {name: {"sha256": sha256(output/name),
                              "size_bytes": (output/name).stat().st_size}
                       for name in ASSEMBLY_ASSETS}
    meta = json.loads((output/"render-metadata.json").read_text())
    assembly = meta.get("assembly", {})
    manifest = json.loads((output/"native-radiance-manifest.json").read_text())
    check(result, "published_native_manifest_matches_assembly_metadata",
          assembly.get("native_manifest") == "native-radiance-manifest.json"
          and assembly.get("native_manifest_sha256") ==
              result["files"]["native-radiance-manifest.json"]["sha256"]
          and manifest.get("radiance_signature") == meta.get("radiance_signature")
          and manifest.get("source_render_sha256") == meta.get("source_sha256")
          and set(manifest.get("frames_sha256", {})) ==
              {f"{n:03d}.npy" for n in range(240)})
    steady = meta.get("steady_reference", {})
    steady_manifest = json.loads((output/"steady-reference-manifest.json").read_text())
    generator = ROOT/"generate_kn_steady_reference.py"
    check(result, "steady_reference_bound_to_current_source_file_and_generator",
          steady.get("file") == "steady-reference.png"
          and steady.get("manifest") == "steady-reference-manifest.json"
          and steady.get("manifest_sha256") ==
              result["files"]["steady-reference-manifest.json"]["sha256"]
          and steady_manifest.get("file") == "steady-reference.png"
          and steady_manifest.get("sha256") ==
              result["files"]["steady-reference.png"]["sha256"]
          and steady_manifest.get("source_render_sha256") == source_hash()
          and steady_manifest.get("generator") == generator.name
          and steady_manifest.get("generator_sha256") == sha256(generator)
          and steady_manifest.get("dimensions") == [meta.get("width"), meta.get("height")]
          and all(steady.get(key) == value for key, value in steady_manifest.items()))
    check(result, "no_partial_delivery_files", not list(output.glob("*.partial.mp4")))
    return finish(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    result = {"scope": __doc__, "preflight_only": args.preflight,
              "verifier_sha256": sha256(__file__), "videos": {}}
    functions = [("numerical_reports", numerical_reports), ("transfer", transfer_checks),
                 ("display", display_checks), ("composition", composition_checks),
                 ("actual_renderer", actual_renderer_checks)]
    if not args.preflight:
        functions += [("metadata", lambda: metadata_checks(args.output)),
                      ("view_diagnostics", lambda: view_diagnostic_checks(args.output)),
                      ("assets", lambda: output_assets(args.output)),
                      ("assembly_artifacts", lambda: assembly_artifact_checks(args.output))]
    for name, function in functions:
        print("Checking "+name, flush=True)
        try:
            result[name] = function()
        except Exception as error:
            result[name] = {"passed": False, "error": str(error)}
    if not args.preflight:
        for mode, name in VIDEO_NAMES.items():
            print("Decoding "+name, flush=True)
            try:
                result["videos"][mode] = video_checks(args.output/name)
            except Exception as error:
                result["videos"][mode] = {"passed": False, "error": str(error)}
    result["passed"] = all(result[n]["passed"] for n, _ in functions) and all(
        r["passed"] for r in result["videos"].values())
    result["checks_passed"] = sum(sum(c["passed"] for c in r.get("checks", []))
        for r in [result[n] for n, _ in functions]+list(result["videos"].values()))
    if not args.preflight and result.get("assembly_artifacts", {}).get("files"):
        result["asset_sha256"] = {name: record["sha256"] for name, record in
                                  result["assembly_artifacts"]["files"].items()}
    if not args.preflight:
        result["release_input_sha256"] = {name: sha256(ROOT/name)
                                            for name in RELEASE_INPUTS}
    result["elapsed_seconds"] = time.perf_counter()-start
    destination = args.output/("verification-preflight.json" if args.preflight else "verification.json")
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print(("PASS: " if result["passed"] else "FAIL: ")+str(destination), flush=True)
    if not result["passed"]:
        for name, record in [(n, result[n]) for n, _ in functions]+list(result["videos"].items()):
            if "error" in record:
                print(name, record["error"], flush=True)
            for c in record.get("checks", []):
                if not c["passed"]:
                    print(name, c, flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
