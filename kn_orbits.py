"""Neutral equatorial circular geodesics in Kerr--Newman, G=M=c=1.

The positive-angular-velocity branch is used. Positive a* is consequently
aligned/prograde; negative a* is antialigned/retrograde. q is the BLACK HOLE'S
dimensionless charge Q/(M sqrt(4 pi epsilon_0 G)), not the particle charge,
which is zero. Only subextremal black holes a*^2+q^2<1 are supported.

Metric and orbit conventions: Pugliese, Quevedo & Ruffini (2013),
https://arxiv.org/abs/1303.6250, Eqs. 1--5 and Sec. II.
Closed neutral E,L expressions: Ivanov & Prodanov (2005),
https://arxiv.org/abs/gr-qc/0504025, Eqs. 19--20, positive branch.
Neutral radial potential: Liu, Lee & Lin (2017),
https://arxiv.org/abs/1706.05466, Eq. 11 with Carter constant zero. Their
Eq. 22 provides an independent ISCO check in verify_kn_orbit_emission.py.
ISCO is found from R''=0 at fixed E,L using the neutral timelike radial
potential, independently of any tabulated Kerr or charged-particle ISCO.
"""
from __future__ import annotations

from functools import lru_cache
import math
import numpy as np


def validate_parameters(spin_a_star=0.9, charge_q=0.4):
    a, q = float(spin_a_star), float(charge_q)
    if not math.isfinite(a) or not math.isfinite(q) or a*a+q*q >= 1:
        raise ValueError("finite subextremal parameters a*^2+q^2<1 required")
    return a, q


def horizon_radius(spin_a_star=0.9, charge_q=0.4):
    a, q = validate_parameters(spin_a_star, charge_q)
    return 1.0+math.sqrt(1.0-a*a-q*q)


def orbital_quantities(radius_rg, spin_a_star=0.9, charge_q=0.4):
    """E,L,u^t,Omega and analytic radial derivatives, positive Omega branch.

    Timelike circular radii must be outside the horizon with D>0. Circular
    orbits inside the ISCO are permitted here for stability verification;
    the emitting disk separately excludes them.
    """
    a, q = validate_parameters(spin_a_star, charge_q)
    r = np.asarray(radius_rg, dtype=np.float64)
    if np.any(~np.isfinite(r)) or np.any(r <= horizon_radius(a, q)):
        raise ValueError("finite circular radius outside the horizon required")
    q2 = q*q
    s = np.sqrt(r-q2)
    sp = 0.5/s
    D = r*r-3*r+2*q2+2*a*s
    if np.any(D <= 0):
        raise ValueError("radius does not admit a timelike circular orbit")
    denominator = r*np.sqrt(D)
    dlogden = 1/r+0.5*(2*r-3+2*a*sp)/D
    NE = r*r-2*r+q2+a*s
    NL = (r*r+a*a)*s-a*(2*r-q2)
    E, L = NE/denominator, NL/denominator
    Ep = ((2*r-2+a*sp)-NE*dlogden)/denominator
    Lp = ((2*r*s+(r*r+a*a)*sp-2*a)-NL*dlogden)/denominator
    W = r*r+a*s
    omega = s/W
    omegap = r*(4*q2-3*r)/(2*s*W*W)
    inverse_ut = denominator/W
    return {"energy": E, "angular_momentum": L,
            "denergy_dr": Ep, "dangular_momentum_dr": Lp,
            "omega": omega, "domega_dr": omegap,
            "ut": 1/inverse_ut, "energy_minus_omega_l": inverse_ut}


def radial_potential_derivatives(radius_rg, energy, angular_momentum,
                                 spin_a_star=0.9, charge_q=0.4):
    """Return R,R',R'' with E,L HELD FIXED (neutral equatorial unit mass).

    R=[E(r^2+a^2)-aL]^2-Delta[r^2+(L-aE)^2]; Sigma^2 (dr/dtau)^2=R.
    Thus R''<0 means radial stability, and ISCO has R=R'=R''=0.
    """
    a, q = validate_parameters(spin_a_star, charge_q)
    r, E, L = np.broadcast_arrays(np.asarray(radius_rg, dtype=float),
                                 np.asarray(energy, dtype=float),
                                 np.asarray(angular_momentum, dtype=float))
    delta = r*r-2*r+a*a+q*q
    dp = 2*r-2
    P = E*(r*r+a*a)-a*L
    B = r*r+(L-a*E)**2
    R = P*P-delta*B
    Rp = 4*r*E*P-dp*B-2*r*delta
    Rpp = 4*E*P+8*r*r*E*E-2*B-4*r*dp-2*delta
    return R, Rp, Rpp


@lru_cache(maxsize=128)
def isco_radius(spin_a_star=0.9, charge_q=0.4):
    """Outermost marginal stable neutral circular orbit on positive branch.

    Bracket the transition from stable to unstable moving inward from 32 M.
    This avoids selecting a horizon/inner root of a squared ISCO polynomial.
    """
    a, q = validate_parameters(spin_a_star, charge_q)
    horizon = horizon_radius(a, q)
    radii = horizon+np.geomspace(1e-7, 32-horizon, 8192)
    D = radii*radii-3*radii+2*q*q+2*a*np.sqrt(radii-q*q)
    radii = radii[D > 0]
    orbit = orbital_quantities(radii, a, q)
    curvature = radial_potential_derivatives(
        radii, orbit["energy"], orbit["angular_momentum"], a, q)[2]
    transitions = np.flatnonzero((curvature[:-1] > 0) & (curvature[1:] <= 0))
    if not transitions.size:
        raise ValueError("failed to bracket a subextremal neutral ISCO")
    j = int(transitions[-1])
    lo, hi = float(radii[j]), float(radii[j+1])
    for _ in range(80):
        mid = (lo+hi)/2
        orb = orbital_quantities(mid, a, q)
        rpp = float(radial_potential_derivatives(mid, orb["energy"], orb["angular_momentum"], a, q)[2])
        if rpp > 0:
            lo = mid
        else:
            hi = mid
        if hi-lo <= 4*np.finfo(float).eps*mid:
            break
    return (lo+hi)/2
