"""Neutral equatorial accretion and geodesic plunge in Kerr--Newman.

Units G=M=c=1, metric signature -+++.  The outer disk is a *prescribed*
slowly accreting flow, not a viscous or magnetohydrodynamic solution:
v_r=-alpha*(H/r)^2*v_phi,ZAMO.  This familiar thin-disk scaling is used as
a small radial boost in the circular emitter's orthonormal radial frame.
The existing circular angular frequency is retained and u.u=-1 exactly.
alpha and H/r set a dimensionless kinematic drift; no disk thickness is
traced and no thermodynamic consistency with Page--Thorne is claimed.

The ISCO launch has the same finite radial four-velocity as the outer
flow. Its E=-u_t and L=u_phi are then conserved along the *exact neutral
equatorial timelike geodesic* inside the ISCO.  Starting with circular
ISCO E,L and zero radial speed would instead approach the ISCO only
asymptotically, so it cannot give a finite continuous launch.

Sources for the calculated metric and geodesic equations:
Wang, Lee & Lin (2022), Eqs. 9--15, https://arxiv.org/abs/2208.11906
Ko, Lee & Lin (2024), Eqs. 1--12, https://arxiv.org/abs/2309.13832
The phenomenological drift scaling is discussed by Penna et al. (2013),
Sec. 1, https://doi.org/10.1093/mnras/sts185
The source brightness of plunging matter is a separate prescription.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import numpy as np

from kn_orbits import horizon_radius, isco_radius, orbital_quantities


def equatorial_metric(radius_rg, spin_a_star=.9, charge_q=.4):
    """Nonzero equatorial Boyer--Lindquist metric components, outside r+."""
    r = np.asarray(radius_rg, dtype=np.float64)
    a, q = float(spin_a_star), float(charge_q)
    rp = horizon_radius(a, q)
    if np.any(~np.isfinite(r)) or np.any(r <= rp):
        raise ValueError("finite radii outside the outer horizon required")
    # Factoring Delta is accurate close to the horizon.
    delta = (r-rp)*(r-(2-rp))
    f = (2*r-q*q)/(r*r)
    return {"tt": -1+f, "tphi": -a*f,
            "phiphi": r*r+a*a+a*a*f, "rr": r*r/delta,
            "delta": delta}


@dataclass(frozen=True)
class FlowParameters:
    spin_a_star: float = .9
    charge_q: float = .4
    alpha_viscosity: float = .1
    aspect_ratio: float = .15
    outer_radius_rg: float = 40.
    horizon_offset_rg: float = 1e-5
    radial_samples: int = 65537

    def __post_init__(self):
        ri = isco_radius(self.spin_a_star, self.charge_q)
        rp = horizon_radius(self.spin_a_star, self.charge_q)
        vals = [self.alpha_viscosity, self.aspect_ratio, self.outer_radius_rg,
                self.horizon_offset_rg]
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("flow parameters must be finite")
        if not (0 < self.alpha_viscosity <= 1 and 0 < self.aspect_ratio <= .5
                and self.outer_radius_rg > ri and 0 < self.horizon_offset_rg < ri-rp
                and isinstance(self.radial_samples, int) and self.radial_samples >= 1025):
            raise ValueError("invalid flow scales or streamline resolution")


class NeutralAccretionFlow:
    """Normalized weak outer drift, smoothly matched to a neutral ISCO plunge.

    ``four_velocity(r)`` gives contravariant ut, ur, uphi, coordinate dr_dt,
    omega, energy and angular_momentum. ``streamline(r)`` gives travel time,
    proper time and accumulated azimuth relative to the ISCO launch. These
    increase along the inward flow; phi-azimuth_rad and t-age_tg are
    material invariants. They are tabulated by convergent quadrature, not
    by rounding the orbital frequencies or by moving already rendered pixels.
    """

    def __init__(self, parameters=None, **parameter_options):
        if parameters is not None and parameter_options:
            raise ValueError("supply a parameter object or keyword options, not both")
        self.parameters = parameters or FlowParameters(**parameter_options)
        p = self.parameters
        self.spin_a_star, self.charge_q = p.spin_a_star, p.charge_q
        self.isco_rg = isco_radius(p.spin_a_star, p.charge_q)
        self.horizon_rg = horizon_radius(p.spin_a_star, p.charge_q)
        self.inner_radius_rg = self.horizon_rg+p.horizon_offset_rg
        self.outer_radius_rg = p.outer_radius_rg
        launch = self._outer_velocity(np.asarray(self.isco_rg))
        self.launch_radial_speed = float(launch["radial_speed_comoving"])
        self.launch_gamma = float(launch["radial_gamma"])
        self.plunge_energy = float(launch["energy"])
        self.plunge_angular_momentum = float(launch["angular_momentum"])
        circ = orbital_quantities(self.isco_rg, p.spin_a_star, p.charge_q)
        self._isco_circular_energy = float(circ["energy"])
        self._isco_circular_momentum = float(circ["angular_momentum"])
        self._gamma2_minus_one = self.launch_radial_speed**2/(1-self.launch_radial_speed**2)
        self._build_streamline()

    def _outer_velocity(self, r):
        p = self.parameters
        metric = equatorial_metric(r, p.spin_a_star, p.charge_q)
        circ = orbital_quantities(r, p.spin_a_star, p.charge_q)
        drag = -metric["tphi"]/metric["phiphi"]
        v_phi_zamo = (circ["omega"]-drag)*metric["phiphi"]/np.sqrt(metric["delta"])
        v_radial = -p.alpha_viscosity*p.aspect_ratio**2*v_phi_zamo
        gamma = 1/np.sqrt(1-v_radial*v_radial)
        ut = gamma*circ["ut"]
        ur = gamma*v_radial*np.sqrt(metric["delta"])/r
        return {"ut": ut, "ur": ur, "uphi": circ["omega"]*ut,
                "omega": circ["omega"], "dr_dt": ur/ut,
                "energy": gamma*circ["energy"],
                "angular_momentum": gamma*circ["angular_momentum"],
                "radial_speed_comoving": v_radial, "radial_gamma": gamma}

    def plunge_potential(self, radius_rg):
        """Stable R for the exact boosted-ISCO geodesic, including its launch.

        The zero-speed circular ISCO potential has a triple root. Expanding
        about that root avoids cancellation of O(1) terms when the finite
        launch speed is only O(1e-3). The two discarded lower Taylor
        coefficients vanish analytically at the independently solved ISCO.
        """
        r = np.asarray(radius_rg, dtype=np.float64)
        delta = equatorial_metric(r, self.spin_a_star, self.charge_q)["delta"]
        x = r-self.isco_rg
        c4 = self._isco_circular_energy**2-1
        base = x*x*x*(c4*x+4*c4*self.isco_rg+2)
        return (1+self._gamma2_minus_one)*base+self._gamma2_minus_one*delta*r*r

    def four_velocity(self, radius_rg):
        r = np.asarray(radius_rg, dtype=np.float64)
        metric = equatorial_metric(r, self.spin_a_star, self.charge_q)
        outer = r >= self.isco_rg
        # The circular formula is never evaluated inside the ISCO/plunge.
        safe = np.maximum(r, self.isco_rg)
        velocity = self._outer_velocity(safe)
        a, E, L = self.spin_a_star, self.plunge_energy, self.plunge_angular_momentum
        P = E*(r*r+a*a)-a*L
        R = self.plunge_potential(r)
        if np.any((~outer) & (R <= 0)):
            raise ValueError("plunge radial potential must remain positive")
        ut = ((r*r+a*a)*P/metric["delta"]+a*(L-a*E))/(r*r)
        uphi = (a*P/metric["delta"]+(L-a*E))/(r*r)
        ur = -np.sqrt(np.maximum(R, 0))/(r*r)
        ut = np.where(outer, velocity["ut"], ut)
        uphi = np.where(outer, velocity["uphi"], uphi)
        ur = np.where(outer, velocity["ur"], ur)
        return {"ut": ut, "ur": ur, "uphi": uphi, "omega": uphi/ut,
                "dr_dt": ur/ut,
                "energy": np.where(outer, velocity["energy"], E),
                "angular_momentum": np.where(outer, velocity["angular_momentum"], L)}

    def redshift(self, radius_rg, photon_lambda, photon_kr_over_E):
        """Frequency ratio E_observer/E_emitter for the actual flowing emitter.

        ``photon_kr_over_E`` is the FORWARD, COVARIANT radial photon momentum
        divided by conserved photon energy at infinity, with either radial
        sign retained: +/-sqrt(R_photon)/Delta. A backward-traced ray with
        positive du/ds has positive forward radial momentum.
        Bolometric intensity transforms by g^4; this method returns g.
        """
        r, lam, kr = np.broadcast_arrays(np.asarray(radius_rg, dtype=np.float64),
                                         np.asarray(photon_lambda, dtype=np.float64),
                                         np.asarray(photon_kr_over_E, dtype=np.float64))
        v = self.four_velocity(r)
        denominator = v["ut"]-lam*v["uphi"]-kr*v["ur"]
        if np.any(~np.isfinite(denominator)) or np.any(denominator <= 0):
            raise ValueError("future-directed null photon must have positive emitter energy")
        return 1/denominator

    def _build_streamline(self):
        p = self.parameters
        # Independent grids avoid subtracting large outer travel times from
        # the short inner plunge. Log(r-r+) resolves the horizon asymptote.
        inner = self.horizon_rg+np.geomspace(p.horizon_offset_rg,
            self.isco_rg-self.horizon_rg, p.radial_samples)
        inner[-1] = self.isco_rg
        outer = np.geomspace(self.isco_rg, p.outer_radius_rg, p.radial_samples)
        outer[0] = self.isco_rg
        def primitive(rr):
            v = self.four_velocity(rr)
            integrands = np.stack((v["ut"], v["uphi"], np.ones_like(rr)))/(-v["ur"])
            return np.concatenate((np.zeros((3, 1)), np.cumsum(
                .5*(integrands[:, 1:]+integrands[:, :-1])*np.diff(rr), axis=1)), axis=1)
        pin, pout = primitive(inner), primitive(outer)
        pin = pin[:, -1:] - pin
        self._r_grid = np.concatenate((inner[:-1], outer))
        self._integrals = np.concatenate((pin[:, :-1], -pout), axis=1)
        self._age_grid = self._integrals[0]
        if np.any(np.diff(self._age_grid) >= 0):
            raise ValueError("flow travel-time coordinate must be strictly monotone")

    def streamline(self, radius_rg):
        r = np.asarray(radius_rg, dtype=np.float64)
        if np.any(~np.isfinite(r)) or np.any(r < self.inner_radius_rg) or np.any(r > self.outer_radius_rg):
            raise ValueError("streamline radius outside tabulated flow domain")
        vals = [np.interp(r, self._r_grid, row) for row in self._integrals]
        return dict(zip(("age_tg", "azimuth_rad", "proper_age_tg"), vals))

    def radius_at_age(self, age_tg, *, clip=False):
        age = np.asarray(age_tg, dtype=np.float64)
        if np.any(~np.isfinite(age)):
            raise ValueError("finite flow age required")
        if not clip and (np.any(age < self._age_grid[-1]) or np.any(age > self._age_grid[0])):
            raise ValueError("age outside tabulated flow domain")
        return np.interp(age, self._age_grid[::-1], self._r_grid[::-1])

    def azimuth_at_age(self, age_tg, *, clip=False):
        """Accumulated flow angle at a material age, without a second r lookup."""
        age = np.asarray(age_tg, dtype=np.float64)
        if np.any(~np.isfinite(age)):
            raise ValueError("finite flow age required")
        if not clip and (np.any(age < self._age_grid[-1]) or np.any(age > self._age_grid[0])):
            raise ValueError("age outside tabulated flow domain")
        return np.interp(age, self._age_grid[::-1], self._integrals[1, ::-1])

    def metadata(self):
        p = self.parameters
        return {**asdict(p), "isco_rg": self.isco_rg,
            "horizon_rg": self.horizon_rg, "inner_radius_rg": self.inner_radius_rg,
            "neutral_emitter_specific_charge": 0.,
            "outer_motion": "normalized circular four-velocity plus comoving radial tetrad boost; prescribed alpha(H/r)^2 ZAMO circular-speed scaling",
            "outer_drift_is_prescribed": True,
            "launch_radial_speed_comoving_c": self.launch_radial_speed,
            "launch_radial_lorentz_factor": self.launch_gamma,
            "plunge_energy_per_rest_mass": self.plunge_energy,
            "plunge_angular_momentum_per_rest_mass": self.plunge_angular_momentum,
            "plunge_motion": "neutral equatorial timelike geodesic with E,L fixed by the finite-drift ISCO launch; Carter constant zero",
            "plunge_radial_potential": "R=[E(r^2+a^2)-aL]^2-Delta[r^2+(L-aE)^2], u^r=-sqrt(R)/r^2",
            "frequency_shift": "g=1/[u^t-lambda*u^phi-(k_r/E)*u^r]; forward covariant photon k_r retained with its radial sign",
            "streamline_coordinate": "age and accumulated azimuth relative to ISCO; t-age(r) and phi-azimuth(r) advect along flow",
            "streamline_quadrature": "cumulative trapezoid on separate exterior geometric and inner log(r-r+) grids; float64 linear lookup",
            "event_horizon": "no emission or timelike emitter is defined on or inside r+",
            "not_solved": ["viscous hydrodynamics", "magnetohydrodynamics", "charged-particle forces", "disk thickness", "thermodynamics", "plunge emissivity"],
            "sources": [
                {"arxiv": "2208.11906", "role": "Kerr-Newman metric and neutral timelike first integrals"},
                {"arxiv": "2309.13832", "role": "neutral Kerr-Newman inspiral/plunge equations and interpretation"},
                {"doi": "10.1093/mnras/sts185", "role": "motivation only for the prescribed alpha(H/r)^2 radial drift scaling"}]}
