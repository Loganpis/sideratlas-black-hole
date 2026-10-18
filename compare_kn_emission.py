"""Reproduce matched source-plane texture comparisons for Sideratlas v1.2.

Run ``python3 compare_kn_emission.py``.  This creates the public
``kn-emission-comparison-diagnostics.json`` without rendering or reading
any video, preview, radiance cache or pre-existing comparison report.

The superseded .32-contrast source is reconstructed from its declared
equation and parameters. A small embedded historical measurement checks
that this reconstruction matches the intermediate implementation captured
before replacement. It is comparison evidence, not a physical source used
by the renderer. Current and v1.1 source models are imported unchanged.

All measurements concern the emitted source or an illustrative fixed-g
scalar tone transform. They do not establish fluid realism, image quality,
agreement with NASA, temporal convergence, or Instagram delivery quality.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np

from kn_dynamic_emission import DynamicDiskEmission
from kn_flow import NeutralAccretionFlow
from kn_orbit_emission import ShearingDiskEmission

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "kn-emission-comparison-diagnostics.json"
SOURCE_NAMES = (
    "compare_kn_emission.py", "kn_dynamic_emission.py", "kn_orbit_emission.py",
    "nasa_emission.py", "kn_flow.py", "kn_orbits.py", "emission.py",
    "render_kn_accretion.py", "render_kn_orbit.py",
)
TIMES = (0., 12.5, 137.3, 366.7)
MODE_BANDS = (("m1_2", 1, 2), ("m3_4", 3, 4), ("m5_8", 5, 8),
              ("m9_16", 9, 16), ("m17_32", 17, 32), ("m33_plus", 33, 1024))

# The original source no longer supplies the production emission. This
# immutable observation record permits an independent reconstruction check.
# It contains parameters, source hashes and source-plane statistics only.
CAPTURED_SMOOTH_REFERENCE = json.loads(r'''
{
  "source_sha256_at_reference_measurement": {
    "kn_dynamic_emission.py": "eab11eb0f8f803c4c86daa0198d66f6c00b00a3dcdf46ab0408f61dd7542dd59",
    "kn_orbit_emission.py": "c1a8587a2b82380ab9c5abc8b87ffab08513ea8dea759569a723cf8925ec68cf",
    "nasa_emission.py": "98acbba4cca4311e6d3edd0620ea0f60e32f9191929163113153b9db8cb010b4",
    "kn_flow.py": "6c925008ef8eff7c570c4de793307a945b345202a24c88a3ab60a308e0eb2c40",
    "kn_orbits.py": "84c097f23ce6427c10466a57db53395105a47047dd10d53da4569158183d08cf",
    "render_kn_accretion.py": "8a4d530517aa6fbae9a62d966feb344c045a1630ebc87fb536b7c3de3440e645",
    "render_kn_orbit.py": "94b393f626a837edf1b91bd6ce00b7db89aa753a2fa572a79ffbb472b63f4fc4"
  },
  "source_metadata": {
    "spin_a_star": 0.9,
    "charge_q": 0.4,
    "outer_radius_rg": 28.0,
    "loop_period_tg": 800.0,
    "contrast": 0.32,
    "plunge_peak_fraction": 0.08,
    "births_per_loop": 32,
    "seed": 13326,
    "inner_radius_rg": 1.1732150807568875,
    "isco_rg": 1.577504566512094,
    "horizon_rg": 1.1732050807568875,
    "radial_peak_rg": 2.2008687199235943,
    "radial_profile": "zero-stress neutral Kerr-Newman Page-Thorne relative flux outside ISCO; exact azimuthal mean, peak=1",
    "radial_quadrature_samples": 65537,
    "flow": {
      "spin_a_star": 0.9,
      "charge_q": 0.4,
      "alpha_viscosity": 0.1,
      "aspect_ratio": 0.15,
      "outer_radius_rg": 40.0,
      "horizon_offset_rg": 1e-05,
      "radial_samples": 65537,
      "isco_rg": 1.577504566512094,
      "horizon_rg": 1.1732050807568875,
      "inner_radius_rg": 1.1732150807568875,
      "neutral_emitter_specific_charge": 0.0,
      "outer_motion": "normalized circular four-velocity plus comoving radial tetrad boost; prescribed alpha(H/r)^2 ZAMO circular-speed scaling",
      "outer_drift_is_prescribed": true,
      "launch_radial_speed_comoving_c": -0.0014129787739372655,
      "launch_radial_lorentz_factor": 1.0000009982560025,
      "plunge_energy_per_rest_mass": 0.7652802622604552,
      "plunge_angular_momentum_per_rest_mass": 1.668154943505371,
      "plunge_motion": "neutral equatorial timelike geodesic with E,L fixed by the finite-drift ISCO launch; Carter constant zero",
      "plunge_radial_potential": "R=[E(r^2+a^2)-aL]^2-Delta[r^2+(L-aE)^2], u^r=-sqrt(R)/r^2",
      "frequency_shift": "g=1/[u^t-lambda*u^phi-(k_r/E)*u^r]; forward covariant photon k_r retained with its radial sign",
      "streamline_coordinate": "age and accumulated azimuth relative to ISCO; t-age(r) and phi-azimuth(r) advect along flow",
      "streamline_quadrature": "cumulative trapezoid on separate exterior geometric and inner log(r-r+) grids; float64 linear lookup",
      "event_horizon": "no emission or timelike emitter is defined on or inside r+",
      "not_solved": [
        "viscous hydrodynamics",
        "magnetohydrodynamics",
        "charged-particle forces",
        "disk thickness",
        "thermodynamics",
        "plunge emissivity"
      ],
      "sources": [
        {
          "arxiv": "2208.11906",
          "role": "Kerr-Newman metric and neutral timelike first integrals"
        },
        {
          "arxiv": "2309.13832",
          "role": "neutral Kerr-Newman inspiral/plunge equations and interpretation"
        },
        {
          "doi": "10.1093/mnras/sts185",
          "role": "motivation only for the prescribed alpha(H/r)^2 radial drift scaling"
        }
      ]
    },
    "perturbation": "three broad Fourier/logarithmic-spiral modes m=1,2,4; bounded contrast; no radial bands",
    "mode_numbers": [
      1,
      2,
      4
    ],
    "mode_weights": [
      0.3,
      0.5,
      0.2
    ],
    "radial_wavenumbers": [
      0.7,
      -3.0,
      2.1
    ],
    "mode_contrast_envelopes": "m1 full; m2 .35 to 1 smoothstep over rISCO to1.8*rISCO; m4 zero through3rg, smoothstep to full at5rg; suppresses unresolved fast fine features near hole",
    "mode_phase_rule": "m*(phi-A(r)+A(r_birth))+k*log(r_birth/rISCO)+seeded_birth_phase",
    "advection": "r_birth=age^-1(age(r)-birth_age); backward along the same normalized outer drift and neutral geodesic plunge used in photon frequency shifts",
    "coherence_lifetime_tg": 50.0,
    "coherence_choice": "phenomenological finite lifetime, about a few bright inner orbits; not a calculated MRI correlation time",
    "birth_blending": "current and previous births, complementary sin^2/cos^2 windows; continuous first time derivative",
    "periodicity": "seeded field births repeat every loop_period_tg; gas orbital frequencies are not quantized",
    "emission_time_rule": "t_em=t_obs-delay_tg separately at each photon/source event",
    "azimuthal_mean": "exactly background profile at fixed radius and emitted time, since every perturbation has integer m>=1",
    "emitted_intensity_bounds_relative_to_background": [
      0.6799999999999999,
      1.32
    ],
    "plunge_emissivity": "separate prescribed weak heating: plunge_peak_fraction*(1-x)^2, x=(r-rplus)/(rISCO-rplus); zero value/slope at ISCO, finite exterior limit at horizon; no source on/inside horizon",
    "plunge_emissivity_is_calculated_from_thermodynamics": false,
    "horizon": "no source on or inside rplus; no boundary glow; g^4, photon delays and capture applied by renderer",
    "surface_emission": "identical isotropic thin equatorial source on both sides; optical-depth handling specified by tracer",
    "dimensional_normalization": "none; mass, accretion rate, temperature and spectrum unspecified",
    "illustrative_prescriptions": [
      "alpha(H/r)^2 slow exterior radial drift",
      "finite-lived smooth emissivity fluctuations and logarithmic spiral modes",
      "weak separate plunging-region heating",
      "periodic births used to construct an eight-second loop"
    ],
    "not_modeled": [
      "GRMHD turbulence or MRI",
      "charged-particle Lorentz forces",
      "plasma currents or charge evolution",
      "thermodynamic spectrum",
      "backreaction",
      "finite disk thickness",
      "returning-radiation heating"
    ],
    "sources": [
      {
        "doi": "10.1086/152990",
        "role": "Page-Thorne stationary thin-disk conservation integral outside ISCO"
      },
      {
        "arxiv": "astro-ph/0309458",
        "role": "motivation for advected and shearing brightness structures; not that paper's specific source model"
      },
      {
        "arxiv": "2211.07469",
        "role": "methodology of time-dependent nonaxisymmetric equatorial sources with per-image retarded times; not the AART stochastic solver"
      }
    ]
  },
  "samples": [
    {
      "radius_rg": 1.7,
      "emitted_time_tg": 0.0,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.07958479262854107,
      "relative_min": 0.9161358008631433,
      "relative_max": 1.1546115203428513
    },
    {
      "radius_rg": 2.2008687199235943,
      "emitted_time_tg": 0.0,
      "relative_mean": 0.9999999999999999,
      "relative_rms_about_azimuthal_mean": 0.10167834654103645,
      "relative_min": 0.838932327221608,
      "relative_max": 1.1931064145340247
    },
    {
      "radius_rg": 3.0,
      "emitted_time_tg": 0.0,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.13193938001976516,
      "relative_min": 0.7500283323457073,
      "relative_max": 1.2016213393501416
    },
    {
      "radius_rg": 5.0,
      "emitted_time_tg": 0.0,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.13948476619330152,
      "relative_min": 0.7083777033039504,
      "relative_max": 1.2570596499370401
    },
    {
      "radius_rg": 10.0,
      "emitted_time_tg": 0.0,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.1394847661933018,
      "relative_min": 0.7191136250107396,
      "relative_max": 1.1849361490585064
    },
    {
      "radius_rg": 20.0,
      "emitted_time_tg": 0.0,
      "relative_mean": 0.9999999999999999,
      "relative_rms_about_azimuthal_mean": 0.13948476619330188,
      "relative_min": 0.701986086550845,
      "relative_max": 1.1826048151201622
    },
    {
      "radius_rg": 1.7,
      "emitted_time_tg": 12.5,
      "relative_mean": 0.9999999999999999,
      "relative_rms_about_azimuthal_mean": 0.06375039539211583,
      "relative_min": 0.9065051533541034,
      "relative_max": 1.0871088049861402
    },
    {
      "radius_rg": 2.2008687199235943,
      "emitted_time_tg": 12.5,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.07095736904814726,
      "relative_min": 0.8886090707796739,
      "relative_max": 1.1177628262948938
    },
    {
      "radius_rg": 3.0,
      "emitted_time_tg": 12.5,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.08183715691877898,
      "relative_min": 0.8419219409600442,
      "relative_max": 1.1094777260448536
    },
    {
      "radius_rg": 5.0,
      "emitted_time_tg": 12.5,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.0982917845721109,
      "relative_min": 0.8428093337097414,
      "relative_max": 1.2319461961177163
    },
    {
      "radius_rg": 10.0,
      "emitted_time_tg": 12.5,
      "relative_mean": 0.9999999999999999,
      "relative_rms_about_azimuthal_mean": 0.10164661442730623,
      "relative_min": 0.8158566264823569,
      "relative_max": 1.1651787967498828
    },
    {
      "radius_rg": 20.0,
      "emitted_time_tg": 12.5,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.12185114115481151,
      "relative_min": 0.7555698941599552,
      "relative_max": 1.2129693302677291
    },
    {
      "radius_rg": 1.7,
      "emitted_time_tg": 137.3,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.051724488603210016,
      "relative_min": 0.9217733121916104,
      "relative_max": 1.0680573236966076
    },
    {
      "radius_rg": 2.2008687199235943,
      "emitted_time_tg": 137.3,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.06869337276141695,
      "relative_min": 0.8944218188524303,
      "relative_max": 1.1006635856805391
    },
    {
      "radius_rg": 3.0,
      "emitted_time_tg": 137.3,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.09214995103956522,
      "relative_min": 0.8185155336912704,
      "relative_max": 1.127058278345638
    },
    {
      "radius_rg": 5.0,
      "emitted_time_tg": 137.3,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.09118284546616666,
      "relative_min": 0.8279560142151382,
      "relative_max": 1.1403812536620934
    },
    {
      "radius_rg": 10.0,
      "emitted_time_tg": 137.3,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.0974579415503286,
      "relative_min": 0.8118138976157654,
      "relative_max": 1.2002012900909276
    },
    {
      "radius_rg": 20.0,
      "emitted_time_tg": 137.3,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.11348128538138447,
      "relative_min": 0.7971986914706619,
      "relative_max": 1.163195161511879
    },
    {
      "radius_rg": 1.7,
      "emitted_time_tg": 366.7,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.057049381417969426,
      "relative_min": 0.9108535167697859,
      "relative_max": 1.0934153947508023
    },
    {
      "radius_rg": 2.2008687199235943,
      "emitted_time_tg": 366.7,
      "relative_mean": 0.9999999999999999,
      "relative_rms_about_azimuthal_mean": 0.09607029515592184,
      "relative_min": 0.8872706180519886,
      "relative_max": 1.1920556401960998
    },
    {
      "radius_rg": 3.0,
      "emitted_time_tg": 366.7,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.09897965312901677,
      "relative_min": 0.8144048169750396,
      "relative_max": 1.1473925815055477
    },
    {
      "radius_rg": 5.0,
      "emitted_time_tg": 366.7,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.08594798347888549,
      "relative_min": 0.8725942994324264,
      "relative_max": 1.145648587657133
    },
    {
      "radius_rg": 10.0,
      "emitted_time_tg": 366.7,
      "relative_mean": 0.9999999999999998,
      "relative_rms_about_azimuthal_mean": 0.1145694882747403,
      "relative_min": 0.8262808157915819,
      "relative_max": 1.24603778873345
    },
    {
      "radius_rg": 20.0,
      "emitted_time_tg": 366.7,
      "relative_mean": 1.0,
      "relative_rms_about_azimuthal_mean": 0.13354867120307162,
      "relative_min": 0.7662120385653721,
      "relative_max": 1.2274467578876127
    }
  ]
}
''')


def smoothstep(x):
    x = np.clip(x, 0., 1.)
    return x*x*(3-2*x)


class SupersededSmoothEmission:
    """Superseded m=1,2,4 equation for comparison, including retarded times.

    The source interface matches ``PreparedRayEmission`` so a controlled
    visual comparison can use the current ray events, flow and transfer.
    This class is never selected by the production renderer.
    """

    def __init__(self, flow=None):
        self.background = ShearingDiskEmission()
        self.flow = flow or NeutralAccretionFlow()
        if self.flow.spin_a_star != .9 or self.flow.charge_q != .4:
            raise ValueError("the historical comparison uses a*=.9, q=.4")
        self.rin = self.background.rin
        self.horizon_rg = self.flow.horizon_rg
        self.peak_radius_rg = self.background.peak_radius_rg
        rng = np.random.default_rng(13326)
        self.phases = rng.uniform(-math.pi, math.pi, (32, 3))
        self.phases[:, 1] = .45*np.sin(2*math.pi*np.arange(32)/32)

    def radial_profile(self, radius):
        r = np.asarray(radius, dtype=float)
        outer = self.background.radial_profile(r)
        inside = np.isfinite(r) & (r > self.horizon_rg) & (r < self.rin)
        x = np.clip((r-self.horizon_rg)/(self.rin-self.horizon_rg),0,1)
        return outer+np.where(inside,.08*(1-x)**2,0.)

    def _evaluate_valid(self, r, p, emitted_time, structured=True):
        base = self.radial_profile(r)
        if not structured:
            return base
        coords = self.flow.streamline(r)
        age, angle = coords["age_tg"], coords["azimuth_rad"]
        t = np.remainder(np.asarray(emitted_time,dtype=float),800)
        t = np.where(t >= 800,0.,t)
        index = np.floor(t/25).astype(np.int64)
        offset = t-index*25
        w = np.sin(.5*math.pi*offset/25)**2
        envelopes = (np.ones_like(r),
            .35+.65*smoothstep((r-self.rin)/(.8*self.rin)),
            smoothstep((r-3)/2))

        def birth(s, j):
            birth_age = age-s
            r0 = self.flow.radius_at_age(birth_age)
            p0 = np.remainder(p-angle+self.flow.azimuth_at_age(birth_age), 2*math.pi)
            lr = np.log(r0/self.rin)
            field = np.zeros_like(r)
            for k, (m, weight, kr) in enumerate(zip((1,2,4), (.3,.5,.2), (.7,-3,2.1))):
                field += weight*envelopes[k]*np.cos(m*p0+kr*lr+self.phases[j,k])
            return field

        modulation = 1+.32*(w*birth(offset,index)+(1-w)*birth(offset+25,(index-1)%32))
        return base*modulation

    def prepare(self, radius, phi, delay_tg=0., zone=None):
        return PreparedSupersededSmoothEmission(self,radius,phi,delay_tg)

    def evaluate(self, radius, phi, emitted_time, zone=None, structured=True):
        return self.prepare(radius,phi,-np.asarray(emitted_time)).evaluate(
            0.,structured=structured)

    def metadata(self):
        return CAPTURED_SMOOTH_REFERENCE["source_metadata"]


class PreparedSupersededSmoothEmission:
    """Retarded source-time adapter for the historical comparison equation."""

    def __init__(self, source, radius, phi, delay):
        self.source = source
        r,p,d = np.broadcast_arrays(np.asarray(radius,dtype=float),
                                    np.asarray(phi,dtype=float),np.asarray(delay,dtype=float))
        self.shape = r.shape
        valid = (np.isfinite(r) & np.isfinite(p) & np.isfinite(d)
                 & (r >= source.flow.inner_radius_rg) & (r <= 28))
        self.valid_indices = np.flatnonzero(valid)
        self.r = r.ravel()[self.valid_indices]
        self.phi = p.ravel()[self.valid_indices]
        self.delay = d.ravel()[self.valid_indices]

    def evaluate(self, observer_time_tg, structured=True):
        time = float(observer_time_tg)
        if not math.isfinite(time):
            raise ValueError("finite observer time required")
        result = np.zeros(math.prod(self.shape),dtype=float)
        if self.valid_indices.size:
            result[self.valid_indices] = self.source._evaluate_valid(
                self.r,self.phi,time-self.delay,structured)
        return result.reshape(self.shape)


# Retain the explicit reference spelling used by the measurement functions.
SupersededSmoothReference = SupersededSmoothEmission


def source_metrics(source, radii):
    """Fixed-r,t angular statistics and matched radial structure functions."""
    phi = np.arange(2048)*2*math.pi/2048
    samples = []
    for t in TIMES:
        J = source.evaluate(radii[:,None], phi[None,:], t)
        B = J/source.radial_profile(radii)[:,None]
        for radius, relative in zip(radii, B):
            ft = np.fft.rfft(relative-1)/relative.size
            power = 2*abs(ft)**2
            power[0], power[-1] = abs(ft[0])**2, abs(ft[-1])**2
            samples.append({
                "radius_rg": float(radius), "emitted_time_tg": t,
                "relative_mean": float(relative.mean()),
                "relative_rms_about_azimuthal_mean": float(relative.std()),
                "relative_min": float(relative.min()), "relative_max": float(relative.max()),
                "relative_percentiles_5_50_95": np.percentile(relative,[5,50,95]).tolist(),
                "angular_variance_by_mode_band": {
                    name:float(power[lo:hi+1].sum()) for name,lo,hi in MODE_BANDS}})
    pp = np.arange(256)*2*math.pi/256
    radial = []
    for delta in (.02,.04,.08,.16):
        rr = np.geomspace(1.7*math.exp(delta/2),20*math.exp(-delta/2),48)
        values = []
        for t in TIMES:
            ra,rb = rr*math.exp(-delta/2),rr*math.exp(delta/2)
            A = source.evaluate(ra[:,None],pp[None,:],t)/source.radial_profile(ra)[:,None]
            B = source.evaluate(rb[:,None],pp[None,:],t)/source.radial_profile(rb)[:,None]
            values.append(float(np.mean((A-B)**2)))
        radial.append({"delta_log_radius":delta,
            "mean_squared_relative_difference":float(np.mean(values)),
            "by_emitted_time":values})
    return {"source_metadata":source.metadata(), "samples":samples,
            "radial_structure_function":radial}


def tone_metrics(source, radii):
    """Scalar tone sensitivity at hypothetical fixed g, before RGB encoding."""
    phi = np.arange(2048)*2*math.pi/2048
    samples = []
    for t in TIMES:
        intensities = source.evaluate(radii[:,None],phi[None,:],t)
        for radius,J in zip(radii,intensities):
            by_g = []
            for g in (.3,.5,.7,1.):
                x = g**4*J/.02
                y = x/(1+x)
                by_g.append({"g":g, "tone_mean":float(y.mean()),
                    "tone_relative_rms":float(y.std()/y.mean()),
                    "tone_min":float(y.min()), "tone_max":float(y.max()),
                    "tone_percentiles_5_50_95":np.percentile(y,[5,50,95]).tolist()})
            samples.append({"radius_rg":float(radius),"emitted_time_tg":t,
                            "by_frequency_shift":by_g})
    return samples


def diagnostics():
    initial_hashes = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                      for name in SOURCE_NAMES}
    old, smooth, current = ShearingDiskEmission(), SupersededSmoothReference(), DynamicDiskEmission()
    radii = np.array([1.7,float(current.peak_radius_rg),3.,5.,10.,20.])
    models = {"v1_1_reference":old,"superseded_smooth_v1_2_reference":smooth,
              "restored_multiscale_v1_2":current}
    metrics = {key:source_metrics(value,radii) for key,value in models.items()}
    tone = {key:tone_metrics(value,radii) for key,value in models.items()}
    checks = []

    def check(name, value, passed):
        checks.append({"name":name,"value":value,"passed":bool(passed)})
        if not passed:
            raise AssertionError(f"{name}: {value}")

    ref = CAPTURED_SMOOTH_REFERENCE["samples"]
    got = metrics["superseded_smooth_v1_2_reference"]["samples"]
    error = 0.
    for recorded, measured in zip(ref,got):
        check_keys = ("relative_mean","relative_rms_about_azimuthal_mean","relative_min","relative_max")
        if recorded["radius_rg"] != measured["radius_rg"] or recorded["emitted_time_tg"] != measured["emitted_time_tg"]:
            raise AssertionError("historical source comparison uses unmatched samples")
        error = max(error,max(abs(recorded[key]-measured[key]) for key in check_keys))
    check("superseded_equation_matches_captured_statistics",error,error < 2e-12 and len(ref)==len(got)==24)
    test_r = np.array([[smooth.horizon_rg+.001,1.3,1.7],[2.2,5.,20.]])
    test_phi = np.array([[.1,1.7,-2.3]])
    test_delay = np.array([[12.7],[-33.4]])
    delayed = smooth.prepare(test_r,test_phi,test_delay).evaluate(137.9)
    scalar = np.empty_like(delayed)
    phi_grid,delay_grid = np.broadcast_arrays(test_phi+np.zeros_like(test_r),
                                              test_delay+np.zeros_like(test_r))
    for index in np.ndindex(test_r.shape):
        scalar[index] = smooth._evaluate_valid(np.asarray(test_r[index]),
            np.asarray(phi_grid[index]),np.asarray(137.9-delay_grid[index]))
    adapter_error = float(np.max(abs(delayed-scalar)))
    check("superseded_prepared_retarded_interface_and_broadcasting",adapter_error,adapter_error < 2e-14)
    smooth_mode = smooth.prepare(test_r,test_phi,test_delay).evaluate(137.9,structured=False)
    check("superseded_unstructured_mode_keeps_PT_and_inner_profile",None,
          np.array_equal(smooth_mode,smooth.radial_profile(test_r)))
    invalid = smooth.evaluate([smooth.horizon_rg-.001,smooth.horizon_rg,float("nan"),28.1],0,0)
    check("superseded_horizon_interior_and_invalid_sources_black",invalid.tolist(),
          np.array_equal(invalid,np.zeros(4)))
    inner_r = np.linspace(smooth.flow.inner_radius_rg,smooth.rin,1001)
    profile_error = float(np.max(abs(smooth.radial_profile(inner_r)-current.radial_profile(inner_r))))
    check("superseded_and_restored_use_identical_plunge_radial_background",profile_error,profile_error < 2e-15)
    for key,data in metrics.items():
        mean_error = max(abs(s["relative_mean"]-1) for s in data["samples"])
        check(key+"_Page_Thorne_azimuthal_mean",mean_error,mean_error < 5e-13)
        minimum = min(s["relative_min"] for s in data["samples"])
        check(key+"_finite_positive_sampled_intensity",minimum,minimum > 0 and
              all(math.isfinite(s["relative_max"]) for s in data["samples"]))
        parseval = max(abs(sum(s["angular_variance_by_mode_band"].values())
                           -s["relative_rms_about_azimuthal_mean"]**2) for s in data["samples"])
        check(key+"_angular_power_Parseval_consistency",parseval,parseval < 2e-13)
    final_hashes = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_NAMES}
    check("source_files_unchanged_during_measurement",None,initial_hashes==final_hashes)

    peak = float(radii[1])
    summary = {key:[s["relative_rms_about_azimuthal_mean"] for s in value["samples"]
                    if s["radius_rg"]==peak] for key,value in metrics.items()}
    average = {key:float(np.mean(value)) for key,value in summary.items()}
    ratio = average["restored_multiscale_v1_2"]/average["superseded_smooth_v1_2_reference"]
    frame_dt = 800/240
    exposure = frame_dt/2
    cadence = []
    for radius,omega in zip(radii,current.flow.four_velocity(radii)["omega"]):
        cadence.append({"radius_rg":float(radius),"omega_rad_per_tg":float(omega),
            "gas_angular_advance_per_frame_rad":float(omega*frame_dt),
            "largest_integer_m_below_local_temporal_Nyquist":int(np.floor(math.pi/(omega*frame_dt))),
            "ideal_box_shutter_amplitude_retention":{
                str(m):float(np.sinc(m*omega*exposure/(2*math.pi))) for m in (1,2,4,8)}})
    return {
        "schema_version":2,"passed":True,"check_count":len(checks),"checks":checks,
        "comparison":"v1.1 shearing source, superseded smooth v1.2 source, and restored multiscale v1.2 source at matched source-plane samples",
        "scope":"Emitted source-plane intensity and analytical/local scalar presentation sensitivity only. These measurements do not compare NASA imagery, full rendered detectors, retarded-image mixing, lensing visibility, encoded video or Instagram recompression. All fluctuation scales are illustrative, not measured turbulence or a GRMHD validation.",
        "reproduction":{"command":"python3 compare_kn_emission.py",
            "requires":"Python, NumPy and the identified source files; no videos, render cache or previous comparison report",
            "python":platform.python_version(),"numpy":np.__version__,
            "helper_sha256":final_hashes["compare_kn_emission.py"]},
        "sampling":{"radii_rg":radii.tolist(),"emitted_times_tg":list(TIMES),
            "azimuths":2048,"azimuth_grid":"phi_j=2*pi*j/2048, j=0..2047",
            "matched_source_geometry":{"spin_a_star":.9,"charge_q":.4,"outer_emitting_radius_rg":28},
            "radial_structure_function_grid":{"inner_radius_rg":1.7,"outer_radius_rg":20.,
                "geometric_midpoint_count":48,"azimuth_count":256,"times_tg":list(TIMES),
                "delta_log_radii":[.02,.04,.08,.16]}},
        "metric_definitions":{
            "relative_intensity":"B(r,phi,t)=J(r,phi,t)/J_background(r), using the common peak-normalized exterior Page-Thorne background",
            "relative_rms_about_azimuthal_mean":"sqrt(mean_phi((B-mean_phi(B))^2)); dimensionless",
            "angular_variance_by_mode_band":"rFFT of B-1 divided by N; twice squared amplitude at nonzero, non-Nyquist modes, with DC and Nyquist not doubled; band sums have Parseval variance units",
            "radial_structure_function":"mean over matched emitted times, uniform azimuth and geometric radial midpoints of [B(r*exp(delta/2),phi,t)-B(r*exp(-delta/2),phi,t)]^2; endpoints remain inside the common domain",
            "tone_relative_sensitivity":"For fixed g, x=g^4*J/.02 and y=x/(1+x), dln(y)/dln(J)=1/(1+x). Listed sensitivity examples set J=1, not a measured detector radiance",
            "tone_compressed_source_plane_metrics":"At each actual sampled J(r,phi,t), set illustrative g in {.3,.5,.7,1}, compute x=g^4*J/.02 and y=x/(1+x); record mean(y), std(y)/mean(y), range and percentiles before palette/RGB encoding. No detector g is inferred",
            "shutter_retention":"For cos[m(phi-Omega*t)] with constant local Omega and ideal box exposure D, amplitude retention=sinc(m*Omega*D/2), sinc(z)=sin(z)/z. The three-point Gaussian source shutter approximates this integral",
            "cadence_limit":"For constant-Omega circular advection at interval dt, m*Omega*dt<=pi is the local angular Nyquist guide; it is not a complete temporal convergence guarantee"},
        "source_sha256":final_hashes,
        "source_sha256_at_reference_measurement":CAPTURED_SMOOTH_REFERENCE["source_sha256_at_reference_measurement"],
        "source_sha256_at_restored_measurement":final_hashes,
        **metrics,
        "tone_compressed_source_plane_metrics":tone,
        "temporal_sampling_guide":{"frame_count_per_cycle":240,"cycle_tg":800,
            "frame_interval_tg":frame_dt,"exposure_tg":exposure,"samples":cadence},
        "tone_sensitivity_examples":[{"g":g,"J":1.,"display_exposure":.02,
            "observed_radiance":g**4,"x":g**4/.02,"retained_relative_contrast":1/(1+g**4/.02)}
            for g in (.3,.5,.7,1.,1.2)],
        "superseded_smooth_equation_reconstruction":{
            "historical_source_sha256":CAPTURED_SMOOTH_REFERENCE["source_sha256_at_reference_measurement"]["kn_dynamic_emission.py"],
            "implementation":"SupersededSmoothEmission in compare_kn_emission.py, with prepare/evaluate retarded-time interface for controlled visual comparisons",
            "equation":"B=1+.32*[w*F(s,j)+(1-w)*F(s+25,j-1)], w=sin(pi*s/50)^2. F=sum_m weight_m*envelope_m(r)*cos[m*(phi-A(r)+A(r_birth))+kr_m*ln(r_birth/rISCO)+phase_jm], r_birth=age_inverse(age(r)-s)",
            "parameters":"m=(1,2,4), weights=(.3,.5,.2), kr=(.7,-3,2.1). Envelopes: m1=1; m2=.35+.65*smoothstep[(r-rISCO)/(.8rISCO)]; m4=smoothstep[(r-3)/2]. Seed default_rng(13326), uniform(-pi,pi), shape(32,3); m2 phases replaced by .45*sin(2*pi*j/32). smoothstep=3x^2-2x^3 with x clipped[0,1]",
            "captured_statistics_checked":"All24 matched means, RMS, minima and maxima measured before replacement; numerical constants embedded in the helper",
            "maximum_absolute_statistic_difference":error},
        "matched_bright_peak_summary":{"radius_rg":peak,"emitted_times_tg":list(TIMES),
            "azimuthal_relative_RMS_by_source":summary,"mean_RMS_by_source":average,
            "restored_to_superseded_mean_RMS_ratio":ratio},
        "causal_findings":[
            "The superseded source had only three broad angular modes; m4 vanished below3rg, near the bright inner region. Its measured bright-peak RMS was much smaller than v1.1",
            "The fixed display tone curve further compresses relative source contrast in bright regions; modest fluctuations can become nearly invisible",
            "The half-frame source shutter retains about99.4% of m1 and97.5% of m2 amplitude at the mean-flux peak, so it is not the main cause of missing low-order texture",
            "The restored positive angular-normalized packet mixture adds overlapping irregular radial supports and a coherent spiral while preserving the exterior Page-Thorne mean",
            f"At the matched bright-peak samples, the restored/superseded mean emitted RMS ratio is {ratio:.6f}; this measures source contrast, not rendered-image fidelity",
            "Angular packet widths follow a declared local cadence guide; the gas frequencies and photon paths are unchanged"],
        "limitations":[
            "Emitted-time snapshots omit per-ray retardation. A rendered view samples different source times and directions at different pixels",
            "Fixed-g tone metrics are controlled source-plane examples, not a detector prediction or RGB/compression assessment",
            "Additional spatial-frequency power is not automatically more physical; source contrast, widths, coherence and periodic forcing remain prescriptions",
            "Radial packet advection, radial phase gradients, birth crossfades, retarded-time gradients and viewpoint motion also affect temporal behavior; angular Nyquist limits alone do not establish full convergence",
            "The captured smooth-source hash identifies a superseded intermediate file. Its explicit reconstruction is a comparison reference, never a production renderer source",
            "This report does not change or validate thermodynamics, charged-plasma physics, photon paths, Page-Thorne flux, the flow or the chosen display exposure"]}


if __name__ == "__main__":
    result = diagnostics()
    OUTPUT.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"passed":result["passed"],"checks":result["check_count"],
        "report":OUTPUT.name,
        "restored_to_smooth_mean_RMS_ratio":result["matched_bright_peak_summary"]["restored_to_superseded_mean_RMS_ratio"]},indent=2))
