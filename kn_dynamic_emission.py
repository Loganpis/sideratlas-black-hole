"""Restrained, advected emission for the v1.2 neutral Kerr--Newman disk.

Calculated quantities are the pre-existing Page--Thorne radial flux outside
the ISCO and the separately calculated, normalized flow characteristics in
``kn_flow``.  Brightness irregularities and inner heating are explicitly
phenomenological.  This is not a GRMHD or thermal plasma calculation.

Localized angular wave packets on several radial scales define each new
field. A coherent m=2 logarithmic spiral supplies its broad organization.
Every angular kernel is positive, has a finite Fourier bandwidth, and has
analytically unit azimuthal mean. Irregular radial weights therefore create
localized arcs and filaments without introducing axisymmetric rings. A
field is transported backwards along the actual radial/azimuthal flow to
its birth, including differential rotation and the timelike plunge. Two
finite-lived fields crossfade with complementary sin^2 weights. Births
repeat after 800 t_g without rounding any gas orbital frequency.

The perturbation has exactly zero azimuthal mean and bounded amplitude.
Positive kernels of angular degree 1,2,4,8,12 are blended convexly, using a
local angular-advection cadence guide. Fine radial structure remains in the
bright inner disk even where rapid rotation limits angular bandwidth. This
guide does not bound every temporal frequency: radial transport, births,
retardation and the changing view also vary. Temporal integration still
requires convergence checks. No turbulence dissipation scale is predicted.
Consequently its instantaneous emitted azimuthal mean outside the ISCO is
the original zero-stress Page--Thorne flux, even though its flowing emitter
velocity is a prescribed small departure from circular motion.  The thin-
disk flux is therefore used as a background prescription, not asserted as
an exact accreting-fluid solution for the boosted flow.

Inside the ISCO, a weak smooth ramp of separate prescribed emissivity is
used. It vanishes with its first derivative at the ISCO and approaches a
finite exterior limit at the horizon. It is not a continuation or prediction
of zero-stress Page--Thorne theory. Its observed near-horizon fading is
therefore calculated through the transfer, not painted into the source.
The renderer supplies the independently calculated g^4, propagation delay,
and photon capture.  No source, luminous boundary, or glow is assigned on
or inside the horizon.  Retarded times are evaluated per photon image.

Sources for methodology rather than a claim to reproduce their solutions:
Page & Thorne (1974), Eqs.11b,12, doi:10.1086/152990
Schnittman & Bertschinger (2004), https://arxiv.org/abs/astro-ph/0309458
Cardenas-Avendano, Lupsasca & Zhu (2023), https://arxiv.org/abs/2211.07469
See kn_flow.py for the separately calculated emitter motion and its sources.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from kn_flow import NeutralAccretionFlow
from kn_orbit_emission import DimensionlessPageThorne
from kn_orbits import horizon_radius, isco_radius


@dataclass(frozen=True)
class DynamicParameters:
    spin_a_star: float = .9
    charge_q: float = .4
    outer_radius_rg: float = 28.
    loop_period_tg: float = 800.
    contrast: float = 3.
    plunge_peak_fraction: float = .08
    births_per_loop: int = 32
    seed: int = 13326

    def __post_init__(self):
        ri = isco_radius(self.spin_a_star, self.charge_q)
        vals = [self.outer_radius_rg, self.loop_period_tg, self.contrast,
                self.plunge_peak_fraction]
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("source parameters must be finite")
        if not (self.outer_radius_rg > ri and self.loop_period_tg > 0
                and 0 <= self.contrast <= 20 and 0 <= self.plunge_peak_fraction <= 1
                and isinstance(self.births_per_loop, int) and self.births_per_loop >= 4):
            raise ValueError("invalid source extent, contrast, heating or birth count")


class DynamicDiskEmission:
    """Relative bolometric emissivity with true characteristic advection.

    ``prepare(r,phi,delay).evaluate(t_obs)`` evaluates the retarded source.
    ``evaluate(r,phi,t_emit)`` is its direct emitted-time equivalent.
    Optional ``zone`` is accepted for the tracer's interface; radii determine
    whether an event belongs to the outer disk or prescribed inner source.
    """

    def __init__(self, outer_radius_rg=28., loop_period_tg=800., contrast=3.,
                 seed=13326, spin_a_star=.9, charge_q=.4, flow=None,
                 plunge_peak_fraction=.08, births_per_loop=32):
        self.parameters = DynamicParameters(spin_a_star, charge_q, outer_radius_rg,
            loop_period_tg, contrast, plunge_peak_fraction, births_per_loop, seed)
        p = self.parameters
        self.rin = isco_radius(p.spin_a_star, p.charge_q)
        self.horizon_rg = horizon_radius(p.spin_a_star, p.charge_q)
        # The upstream margin permits exact backward characteristics at the
        # outer emitting rim; no source-radius clipping is used for advection.
        self.flow = flow or NeutralAccretionFlow(spin_a_star=p.spin_a_star,
            charge_q=p.charge_q, outer_radius_rg=max(40., p.outer_radius_rg+2.))
        if (self.flow.spin_a_star != p.spin_a_star or self.flow.charge_q != p.charge_q
                or self.flow.outer_radius_rg <= p.outer_radius_rg):
            raise ValueError("source and flow must share geometry, with an upstream flow margin")
        self.flux_model = DimensionlessPageThorne(p.spin_a_star, p.charge_q,
                                                 p.outer_radius_rg)
        a, b = self.rin, p.outer_radius_rg
        q = (math.sqrt(5)-1)/2
        for _ in range(90):
            x, y = b-q*(b-a), a+q*(b-a)
            if float(self.flux_model.flux_shape(x)) < float(self.flux_model.flux_shape(y)):
                a = x
            else:
                b = y
        self.peak_radius_rg = (a+b)/2
        self.radial_peak = float(self.flux_model.flux_shape(self.peak_radius_rg))
        rng = np.random.default_rng(p.seed)
        j = np.arange(p.births_per_loop)
        self.spiral_phases = .45*np.sin(2*math.pi*j/p.births_per_loop)
        self.spiral_weight = .30
        self.spiral_contrast = .85
        # Broad structures overlap throughout the disk. Finer structures
        # occupy irregular sectors and birth radii rather than fixed annuli.
        coarse, fine = 8, 18
        extent = math.log(p.outer_radius_rg/self.rin)
        centers = np.r_[np.linspace(-.08, extent+.08, coarse),
                        np.linspace(-.02, extent+.02, fine)]
        spacing = np.r_[np.full(coarse, (extent+.16)/(coarse-1)),
                        np.full(fine, (extent+.04)/(fine-1))]
        self.packet_centers = (math.log(self.rin)+centers[None, :]
            +rng.uniform(-.26, .26, (p.births_per_loop, coarse+fine))*spacing)
        self.packet_widths = np.concatenate((
            rng.uniform(.27, .41, (p.births_per_loop, coarse)),
            rng.uniform(.055, .125, (p.births_per_loop, fine))), axis=1)
        self.packet_amplitudes = np.concatenate((
            rng.uniform(.55, .95, (p.births_per_loop, coarse)),
            rng.uniform(.85, 1.40, (p.births_per_loop, fine))), axis=1)
        self.packet_phases = rng.uniform(-math.pi, math.pi,
                                         (p.births_per_loop, coarse+fine))
        self.packet_slopes = rng.uniform(-3., 3., (p.births_per_loop, coarse+fine))
        self.packet_caps = np.r_[np.full(coarse, 4), np.full(fine, 12)]
        self.kernel_degrees = (1, 2, 4, 8, 12)
        self.kernel_means = {n: math.comb(2*n, n)/2**n for n in self.kernel_degrees}
        self.kernel_maximum = 2**12/self.kernel_means[12]
        max_weight = self.spiral_weight+float(np.max(self.packet_amplitudes.sum(axis=1)))
        self.minimum_modulation = 1/(1+p.contrast*max_weight)
        self.birth_interval_tg = p.loop_period_tg/p.births_per_loop
        lifetime = 2*self.birth_interval_tg
        outer_age = float(self.flow.streamline(p.outer_radius_rg)["age_tg"])
        # Fail rather than silently clipping a backward-traced fluid element.
        self.flow.radius_at_age(outer_age-lifetime)
        self._prepare_sampling_bandwidth(lifetime)

    def _prepare_sampling_bandwidth(self, lifetime):
        """Birth-coordinate bandwidth, constant for each transported parcel.

        The same finite-lived parcel can rotate faster later in its life.
        Bound that speed on nine future characteristic samples, with a two
        percent safety margin. Clipping here is used solely for a conservative
        angular-speed bound beyond the exterior computational cutoff; actual
        source advection always uses strict, unclipped backward trajectories.
        """
        r = self.horizon_rg+np.geomspace(self.flow.inner_radius_rg-self.horizon_rg,
            self.flow.outer_radius_rg-self.horizon_rg, 32769)
        r[0], r[-1] = self.flow.inner_radius_rg, self.flow.outer_radius_rg
        age = self.flow.streamline(r)["age_tg"]
        maximum_omega = np.zeros_like(r)
        for dt in np.linspace(0, lifetime, 9):
            downstream = self.flow.radius_at_age(age+dt, clip=True)
            maximum_omega = np.maximum(maximum_omega,
                np.abs(self.flow.four_velocity(downstream)["omega"]))
        self._bandwidth_age = age[::-1].copy()
        # 240 unique source frames per delivered cycle. This changes only
        # feature width, never the gas frequency, source clock or ray paths.
        frame_interval = self.parameters.loop_period_tg/240
        self._bandwidth = (.99*math.pi/(1.02*np.maximum(maximum_omega, 1e-12)
                                      *frame_interval))[::-1].copy()

    def angular_velocity(self, radius_rg):
        return self.flow.four_velocity(radius_rg)["omega"]

    def radial_profile(self, radius_rg):
        """Exact exterior PT background plus the disclosed weak inner ramp."""
        r = np.asarray(radius_rg, dtype=np.float64)
        outer = np.isfinite(r) & (r > self.rin) & (r <= self.parameters.outer_radius_rg)
        safe = np.where(outer, r, self.rin)
        profile = np.where(outer, self.flux_model.flux_shape(safe)/self.radial_peak, 0.)
        inside = np.isfinite(r) & (r > self.horizon_rg) & (r < self.rin)
        x = np.clip((r-self.horizon_rg)/(self.rin-self.horizon_rg), 0., 1.)
        # Finite exterior emissivity at r+ lets the independently calculated
        # g^4, delay and capture produce the observed fading. Zero value and
        # slope at ISCO preserve the exterior stress-free background edge.
        ramp = self.parameters.plunge_peak_fraction*(1-x)*(1-x)
        return profile+np.where(inside, ramp, 0.)

    def prepare(self, radius_rg, phi_rad, delay_tg=0., zone=None):
        return PreparedDynamicEmission(self, radius_rg, phi_rad, delay_tg, zone)

    def evaluate(self, radius_rg, phi_rad, emission_time_tg, zone=None, structured=True):
        # t_obs=0 and delay=-t_emit exactly reproduces a directly supplied
        # source time, with all broadcasting handled in one implementation.
        return self.prepare(radius_rg, phi_rad, -np.asarray(emission_time_tg), zone).evaluate(
            0., structured=structured)

    def metadata(self):
        p = self.parameters
        return {**asdict(p), "inner_radius_rg": self.flow.inner_radius_rg,
            "isco_rg": self.rin, "horizon_rg": self.horizon_rg,
            "radial_peak_rg": self.peak_radius_rg,
            "radial_profile": "zero-stress neutral Kerr-Newman Page-Thorne relative flux outside ISCO; exact azimuthal mean, peak=1",
            "radial_quadrature_samples": self.flux_model.samples,
            "flow": self.flow.metadata(),
            "perturbation": "positive, analytically mean-normalized localized angular packets on several radial scales, plus a coherent m2 logarithmic spiral; no axisymmetric rings",
            "packet_count": len(self.packet_caps),
            "coarse_packets": 8, "fine_packets": 18,
            "radial_support_log_half_widths": {"coarse": [.27, .41], "fine": [.055, .125]},
            "packet_radial_weights": "C-infinity compact bumps exp(1-1/(1-x^2)) for |x|<1, x=(log(r_birth)-jittered_birth_center)/width",
            "angular_kernel": "K_n(theta)=(1+cos(theta))^n/[binomial(2n,n)/2^n]; positive, unit mean, finite degree n",
            "kernel_degrees": list(self.kernel_degrees),
            "kernel_blending": "convex mixtures of degrees1,2,4,8,12; broad packets capped at4; highest active angular degree follows a birth-parcel angular-advection cadence guide",
            "sampling_bandwidth": "0.99*pi/[1.02*max_future_abs(Omega)*dt_frame]; nine characteristic times per50t_g parcel lifetime; 240 unique frames/cycle; smooth positive-kernel blending, gas frequencies unchanged",
            "sampling_limit": "angular guide only, not a bound on every temporal frequency; radial transport, birth crossfades, retardation and changing viewpoint require separate temporal convergence checks",
            "packet_phase_rule": "phi_birth+k_j*(log(r_birth)-log(r_center,j))+seeded_birth_phase; coherent spiral 2*phi_birth-3*log(r_birth/rISCO)",
            "normalization": "M=(1+C*sum_j w_j*K_j)/(1+C*sum_j w_j), with nonnegative radial weights and unit-mean positive kernels; current/previous fields crossfade convexly",
            "advection": "r_birth=age^-1(age(r)-birth_age); backward along the same normalized outer drift and neutral geodesic plunge used in photon frequency shifts",
            "coherence_lifetime_tg": 2*self.birth_interval_tg,
            "coherence_choice": "phenomenological finite lifetime, about a few bright inner orbits; not a calculated MRI correlation time",
            "birth_blending": "current and previous births, complementary sin^2/cos^2 windows; continuous first time derivative",
            "periodicity": "seeded field births repeat every loop_period_tg; gas orbital frequencies are not quantized",
            "emission_time_rule": "t_em=t_obs-delay_tg separately at each photon/source event",
            "azimuthal_mean": "exactly background profile at fixed radius and emitted time, since every perturbation has integer m>=1",
            "emitted_intensity_bounds_relative_to_background": [self.minimum_modulation, self.kernel_maximum],
            "plunge_emissivity": "separate prescribed weak heating: plunge_peak_fraction*(1-x)^2, x=(r-rplus)/(rISCO-rplus); zero value/slope at ISCO, finite exterior limit at horizon; no source on/inside horizon",
            "plunge_emissivity_is_calculated_from_thermodynamics": False,
            "horizon": "no source on or inside rplus; no boundary glow; g^4, photon delays and capture applied by renderer",
            "surface_emission": "identical isotropic thin equatorial source on both sides; optical-depth handling specified by tracer",
            "dimensional_normalization": "none; mass, accretion rate, temperature and spectrum unspecified",
            "illustrative_prescriptions": ["alpha(H/r)^2 slow exterior radial drift",
                "finite-lived multiscale emissivity packets and a logarithmic spiral",
                "weak separate plunging-region heating", "periodic births used to construct an eight-second loop"],
            "not_modeled": ["GRMHD turbulence or MRI", "charged-particle Lorentz forces",
                "plasma currents or charge evolution", "thermodynamic spectrum", "backreaction",
                "finite disk thickness", "returning-radiation heating"],
            "sources": [
                {"doi": "10.1086/152990", "role": "Page-Thorne stationary thin-disk conservation integral outside ISCO"},
                {"arxiv": "astro-ph/0309458", "role": "motivation for advected and shearing brightness structures; not that paper's specific source model"},
                {"arxiv": "2211.07469", "role": "methodology of time-dependent nonaxisymmetric equatorial sources with per-image retarded times; not the AART stochastic solver"}]}


class PreparedDynamicEmission:
    """Compact hit-only source arrays; two smooth characteristic lookups/frame."""

    def __init__(self, source, radius_rg, phi_rad, delay_tg, zone=None):
        self.source = source
        r, phi, delay = np.broadcast_arrays(np.asarray(radius_rg, dtype=np.float64),
            np.asarray(phi_rad, dtype=np.float64), np.asarray(delay_tg, dtype=np.float64))
        self.shape = r.shape
        valid = (np.isfinite(r) & np.isfinite(phi) & np.isfinite(delay)
                 & (r >= source.flow.inner_radius_rg) & (r <= source.parameters.outer_radius_rg))
        self.valid_indices = np.flatnonzero(valid)
        # A stable radius ordering greatly improves the locality of repeated
        # monotone characteristic lookups. Scatter through these same indices
        # restores the caller's exact pixel/event ordering on every evaluation.
        if self.valid_indices.size >= 64:
            order = np.argsort(r.ravel()[self.valid_indices], kind="stable")
            self.valid_indices = self.valid_indices[order]
        self.radius = r.ravel()[self.valid_indices]
        self.phi = phi.ravel()[self.valid_indices]
        self.delay = delay.ravel()[self.valid_indices]
        self.background = source.radial_profile(self.radius)
        coords = source.flow.streamline(self.radius)
        self.age = coords["age_tg"]
        self.azimuth = coords["azimuth_rad"]

    @staticmethod
    def _smoothstep(x):
        x = np.clip(x, 0., 1.)
        return x*x*(3-2*x)

    def _kernel(self, theta, bandwidth, degree_cap):
        source = self.source
        c = 1+np.cos(theta)
        square = c*c
        k2 = square/source.kernel_means[2]
        w2 = self._smoothstep((bandwidth-2)/.25)
        if degree_cap <= 2:
            return (1-w2)*c+w2*k2
        fourth = square*square
        k4 = fourth/source.kernel_means[4]
        w4 = self._smoothstep((bandwidth-4)/.6)
        if degree_cap <= 4:
            middle = (1-w4)*k2+w4*k4
        else:
            eighth = fourth*fourth
            k8 = eighth/source.kernel_means[8]
            twelfth = eighth*fourth
            k12 = twelfth/source.kernel_means[12]
            w8 = self._smoothstep((bandwidth-8)/1.2)
            w12 = self._smoothstep((bandwidth-12)/2.)
            high = (1-w12)*k8+w12*k12
            middle = (1-w4)*k2+w4*((1-w8)*k4+w8*high)
        return (1-w2)*c+w2*middle

    def _birth_modulation(self, age_since_birth, birth_index):
        source = self.source
        birth_age = self.age-age_since_birth
        r0 = source.flow.radius_at_age(birth_age)
        if hasattr(source.flow, "azimuth_at_age"):
            a0 = source.flow.azimuth_at_age(birth_age)
        else:
            a0 = source.flow.streamline(r0)["azimuth_rad"]
        phi0 = np.remainder(self.phi-self.azimuth+a0, 2*math.pi)
        log_r0 = np.log(r0)
        bandwidth = np.interp(birth_age, source._bandwidth_age, source._bandwidth)
        spiral_phase = 2*phi0-3*(log_r0-math.log(source.rin))+source.spiral_phases[birth_index]
        spiral_gain = self._smoothstep((bandwidth-2)/.25)*source.spiral_contrast
        numerator = 1+source.parameters.contrast*source.spiral_weight*(1+spiral_gain*np.cos(spiral_phase))
        denominator = np.full_like(numerator, 1+source.parameters.contrast*source.spiral_weight)
        # Each hit intersects only a few compact radial supports. No dense
        # N_hits x N_packets array is formed or retained in the prepared state.
        for j, cap in enumerate(source.packet_caps):
            center = source.packet_centers[birth_index, j]
            x = (log_r0-center)/source.packet_widths[birth_index, j]
            active = np.abs(x) < 1
            if not np.any(active):
                continue
            idx = birth_index[active]
            radial = np.exp(1-1/(1-x[active]**2))
            weight = source.parameters.contrast*source.packet_amplitudes[idx, j]*radial
            theta = (phi0[active]+source.packet_slopes[idx, j]*(log_r0[active]-center[active])
                     +source.packet_phases[idx, j])
            kernel = self._kernel(theta, bandwidth[active], int(cap))
            numerator[active] += weight*kernel
            denominator[active] += weight
        return numerator/denominator

    def evaluate(self, observer_time_tg, structured=True):
        time = float(observer_time_tg)
        if not math.isfinite(time):
            raise ValueError("observer time must be finite")
        result = np.zeros(math.prod(self.shape), dtype=np.float64)
        if not self.valid_indices.size:
            return result.reshape(self.shape)
        p = self.source.parameters
        if not structured or p.contrast == 0:
            result[self.valid_indices] = self.background
            return result.reshape(self.shape)
        # Reduce time before selecting births so negative propagation delays
        # and negative emitted times use the identical periodic forcing law.
        t_emit = np.remainder(time-self.delay, p.loop_period_tg)
        # Floating-point remainder of a negative infinitesimal can round up
        # to the period itself. It denotes the same birth seam as zero.
        t_emit = np.where(t_emit >= p.loop_period_tg, 0., t_emit)
        interval = self.source.birth_interval_tg
        index = np.floor(t_emit/interval).astype(np.int64)
        age = t_emit-index*interval
        weight = np.sin(.5*math.pi*age/interval)**2
        modulation = weight*self._birth_modulation(age, index)
        modulation += (1-weight)*self._birth_modulation(age+interval, (index-1) % p.births_per_loop)
        result[self.valid_indices] = self.background*modulation
        return result.reshape(self.shape)


# Drop-in name for callers transitioning from the v1.1 module.
ShearingDiskEmission = DynamicDiskEmission
