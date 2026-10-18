"""Kerr–Newman rays through a thin plunging layer and an opaque outer disk.

The unchanged null equations and pole-regular backward RK4 method from
``kn_orbit_trace`` are reused in a separate kernel. Each equatorial crossing
inside the ISCO is recorded in camera-to-source order; the first opaque
outer-disk intersection terminates the ray. Additional images follow from
these intersections, with no prescribed number, width or brightness of rings.

Geometry and momentum are physical calculations in the stated stationary
metric. The emissivity and optical thickness assigned to the plunging layer
are separate illustrative model choices. All gas remains outside an explicit
near-horizon cutoff; neither the horizon nor its interior emits.

Reference: Wang, Lee & Lin (2022), PRD 106, 084048, equations (9)–(15), (21)
https://arxiv.org/abs/2208.11906
"""
from __future__ import annotations
import ctypes
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
from kn_orbit_trace import KerrNewmanOrbitTracer
from kn_orbits import horizon_radius

ROOT = Path(__file__).resolve().parent
_LIBRARY = None


def _load_kernel():
    source = ROOT/'kn_accretion_trace.cpp'
    binary = ROOT/('_kn_accretion_trace.dylib' if sys.platform == 'darwin' else '_kn_accretion_trace.so')
    if not binary.exists() or binary.stat().st_mtime < source.stat().st_mtime:
        compiler = shutil.which('clang++') or shutil.which('g++')
        if not compiler:
            raise RuntimeError('The supplied Kerr–Newman kernel requires clang++ or g++')
        subprocess.run([compiler, '-dynamiclib' if sys.platform == 'darwin' else '-shared',
                        '-std=c++17', '-O3', '-fPIC', '-pthread', str(source), '-o', str(binary)],
                       check=True, capture_output=True, text=True)
    lib = ctypes.CDLL(str(binary))
    function = lib.trace_kn_accretion
    f64 = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    f32 = np.ctypeslib.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS')
    i16 = np.ctypeslib.ndpointer(dtype=np.int16, flags='C_CONTIGUOUS')
    u8 = np.ctypeslib.ndpointer(dtype=np.uint8, flags='C_CONTIGUOUS')
    function.argtypes = [f64, f64, ctypes.c_longlong, *([ctypes.c_double]*9),
                         *([ctypes.c_int]*3), *([f32]*4), i16, u8, u8, u8, f64]
    function.restype = None
    return lib, function


class KerrNewmanAccretionTracer(KerrNewmanOrbitTracer):
    """Neutral null paths with ordered thin-plunge and opaque-disk events.

    For every ray, fields ``r, phi, delay_tg, kr_over_E, hit_layer, hit_kind,
    mask`` have shape ``(..., max_hits)``. ``lambda, eta, status, shadow,
    overflow, constraint_residual`` have the detector/sample shape ``(...)``.
    ``hit_kind`` is 0 for padding, 1 for opaque outer disk, 2 for thin plunge.
    ``hit_layer`` counts all finite equatorial crossings, including empty ones.
    ``status`` is 0 escape, 1 opaque disk, 2 capture, 3 exhausted integration,
    4 nonfinite integration, 5 exceeded event storage. Status 3–5 raises by
    default rather than silently dropping unresolved images.

    The future-directed covariant photon momentum is ``k_r/E=kr_over_E``;
    for a neutral emitter with u^theta=0, g^-1 = ut-lambda*uphi-(k_r/E)*ur.
    Its local surface-normal cosine is mu=g*sqrt(eta)/r at a plane crossing.
    """
    def __init__(self, width, height, fov=72, roll_deg=0, spin=.9, charge_q=.4,
                 rin=None, rout=28, angle_step=.012, max_angle=20*math.pi,
                 threads=None, max_hits=6, horizon_offset_rg=1e-4):
        super().__init__(width, height, fov=fov, roll_deg=roll_deg, spin=spin,
                         charge_q=charge_q, rin=rin, rout=rout, angle_step=angle_step,
                         max_angle=max_angle, threads=threads)
        if not isinstance(max_hits, int) or not 1 <= max_hits <= 64:
            raise ValueError('max_hits must be an integer in [1,64]')
        self.max_hits = max_hits
        self.horizon_offset_rg = float(horizon_offset_rg)
        self.rmin = horizon_radius(spin, charge_q) + self.horizon_offset_rg
        if not 0 < self.horizon_offset_rg < self.rin-horizon_radius(spin, charge_q):
            raise ValueError('horizon cutoff must lie strictly between the horizon and the ISCO')
        # An exterior-scattered null ray can turn no deeper than the
        # prograde equatorial photon orbit. A capture cutoff below that orbit
        # cannot terminate a ray that would have escaped after a radial turn.
        # This bound is also needed for user-supplied nondefault cutoffs.
        lo, hi = horizon_radius(spin, charge_q), 4.
        for _ in range(80):
            rr = (lo+hi)/2
            f = rr*rr-3*rr+2*charge_q*charge_q+2*abs(spin)*math.sqrt(rr-charge_q*charge_q)
            if f < 0: lo = rr
            else: hi = rr
        self.minimum_scattering_turn_rg = (lo+hi)/2
        if self.rmin >= self.minimum_scattering_turn_rg:
            raise ValueError('capture cutoff must lie below the innermost exterior photon orbit')
        self.max_parameter = float(max_angle)
        # Extra approach steps regularize the BL coordinate singularity,
        # without extending the configured total scaled Mino time budget.
        self.max_steps += 512
        self.metadata.update({
            'model': 'Kerr-Newman null geodesics; all thin-plunge intersections then first opaque outer disk',
            'kernel': 'kn_accretion_trace.cpp; original kn_orbit_trace equations preserved',
            'plunge_inner_cutoff_rg': self.rmin,
            'horizon_cutoff_offset_rg': self.horizon_offset_rg,
            'minimum_exterior_photon_turn_rg': self.minimum_scattering_turn_rg,
            'capture_cutoff_safety': 'cutoff strictly below prograde equatorial photon orbit; no exterior escaping radial turn is discarded',
            'near_horizon_step_limit': 'h <= 0.15 (1/r_plus-u)/(du/ds) on incoming approach',
            'max_emitting_intersections_stored': max_hits,
            'storage_limit_policy': 'unresolved overflow raises; test convergence by increasing capacity',
            'image_order': 'number of finite equatorial crossings before the emitting event; no artificial bands',
            'event_order': 'near to far along a backward ray',
            'disk_opacity': 'outer disk opaque; inner plunging sheet transfer prescribed separately',
            'photon_radial_momentum': 'future-directed covariant k_r/E = sign(backward du/ds) sqrt(R)/Delta',
            'emitter_redshift_rule': 'g^-1=ut-lambda*uphi-(k_r/E)*ur; utheta=0',
            'surface_normal_cosine': 'mu=g*sqrt(eta)/r',
        })

    def trace(self, inclination_deg, *, allow_unresolved=False):
        return self.trace_samples(self.x, self.y, inclination_deg, allow_unresolved=allow_unresolved)

    def trace_samples(self, alpha, beta, inclination_deg, *, allow_unresolved=False):
        """Trace arbitrary detector coordinates, including sparse refined pixels.

        ``alpha`` and ``beta`` are unrolled/rotated physical screen coordinates
        in r_g, as in ``self.x, self.y``. Broadcastable arrays preserve their
        common shape. This permits independent subpixel trajectories near the
        critical curve without interpolating image orders or curved paths.
        """
        global _LIBRARY
        angle = float(inclination_deg)
        if not math.isfinite(angle):
            raise ValueError('observer angle must be finite')
        alpha, beta = np.broadcast_arrays(np.asarray(alpha, dtype=np.float64), np.asarray(beta, dtype=np.float64))
        if not alpha.size or not np.isfinite(alpha).all() or not np.isfinite(beta).all():
            raise ValueError('sample coordinates must be nonempty and finite')
        alpha, beta = np.ascontiguousarray(alpha), np.ascontiguousarray(beta)
        inc = math.radians(math.remainder(angle, 360.))
        si, ci = math.sin(inc), math.cos(inc)
        if abs(si) < 2e-15: si = 0.
        if abs(ci) < 2e-15: ci = 0.
        if _LIBRARY is None:
            _LIBRARY = _load_kernel()
        started = time.perf_counter()
        shape = alpha.shape
        event_shape = (*shape, self.max_hits)
        r = np.empty(event_shape, np.float32)
        phi, delay, kr = (np.empty_like(r) for _ in range(3))
        layer = np.empty(event_shape, np.int16)
        kind = np.empty(event_shape, np.uint8)
        status, overflow = (np.empty(shape, np.uint8) for _ in range(2))
        residual = np.empty(shape, np.float64)
        _LIBRARY[1](alpha, beta, alpha.size, self.spin, self.charge_q, si, ci, self.rin, self.rout,
                    self.rmin, self.angle_step, self.max_parameter, self.max_steps, self.threads,
                    self.max_hits, r, phi, delay, kr, layer, kind, status, overflow, residual)
        if not allow_unresolved and np.any(status >= 3):
            raise RuntimeError(f'Kerr–Newman trace left {np.count_nonzero(status >= 3)} unresolved rays '
                               f'({np.count_nonzero(overflow)} event overflows); increase integration/event budget')
        for key, value in (('radius', r), ('azimuth', phi), ('delay', delay), ('radial momentum', kr)):
            if not np.isfinite(value).all():
                raise RuntimeError(f'Kerr–Newman {key} map contains nonfinite values')
        mask = kind > 0
        lam = np.asarray(-alpha*si, dtype=np.float64)
        eta = np.asarray(beta*beta+(alpha*alpha-self.spin*self.spin)*ci*ci, dtype=np.float64)
        # Keep the compact event positions, but reconstruct the on-shell
        # covariant momentum in double precision AT THOSE STORED POSITIONS.
        # Independently rounding r and k_r to float32 would be inconsistent
        # near Delta=0, corrupting both reconstructed nullity and redshift.
        # Work layer by layer on emitting samples to bound temporary memory.
        negative = np.signbit(kr)
        del kr
        kr = np.zeros(event_shape, dtype=np.float64)
        rp = horizon_radius(self.spin, self.charge_q)
        rm = 2-rp
        for slot in range(self.max_hits):
            use = mask[...,slot]
            rr = r[...,slot][use].astype(np.float64)
            ll, ee = lam[use], eta[use]
            u = 1/rr
            kk = ee+(ll-self.spin)**2
            potential = 1+(self.spin*self.spin-ll*ll-ee)*u*u+2*kk*u**3-(self.spin*self.spin*ee+self.charge_q*self.charge_q*kk)*u**4
            delta_scaled = (rr-rp)*(rr-rm)/(rr*rr)
            value = np.sqrt(np.maximum(0,potential))/delta_scaled
            kr[...,slot][use] = np.where(negative[...,slot][use],-value,value)
        self.metadata['radial_momentum_storage'] = 'float64 on-shell at stored float32 event radius; branch sign from integrated ray'
        self.metadata['last_trace_seconds'] = time.perf_counter()-started
        self.metadata['last_view_angle_deg'] = angle
        self.metadata['last_constraint_residual_max'] = float(residual.max())
        self.metadata['last_max_hits_observed'] = int(mask.sum(axis=-1).max())
        self.metadata['last_max_image_order_observed'] = int(layer.max())
        return {'r': r, 'phi': phi, 'delay_tg': delay, 'kr_over_E': kr, 'lambda': lam, 'eta': eta,
                'mask': mask, 'hit_layer': layer, 'hit_kind': kind, 'shadow': status == 2,
                'status': status, 'overflow': overflow, 'constraint_residual': residual}

    # Alias spelling for callers describing sampling rather than sample arrays.
    trace_sampled = trace_samples
