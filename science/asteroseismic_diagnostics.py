"""Independent, reproducible asteroseismic diagnostics for the illustrative browser lab.

Units: all frequencies in microHz, effective temperature in kelvin. The
scaling relations are empirical approximations, not stellar-evolution fits.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

NU_MAX_SUN = 3090.0
DELTA_NU_SUN = 135.1
TEFF_SUN = 5772.0


@dataclass(frozen=True)
class ScalingResult:
    mass_solar: float
    radius_solar: float
    covariance: np.ndarray  # 2x2 in linear solar units, first-order delta method


def asymptotic_modes(orders, delta_nu: float, epsilon: float, d0: float,
                     degrees=(0, 1, 2, 3)) -> np.ndarray:
    """Return structured (n,l,nu_microhz) p-mode comb using Tassoul's approximation.

    Small separation nu(n,0)-nu(n-1,2) = 6*d0, for the same epsilon.
    This is a pedagogical high-radial-order asymptotic relation, not a
    peak-bagged frequency set for KIC 5184732.
    """
    if not (np.isfinite(delta_nu) and delta_nu > 0 and np.isfinite(epsilon)
            and np.isfinite(d0) and d0 >= 0):
        raise ValueError("delta_nu must be positive; epsilon and d0 finite, d0 nonnegative")
    ns = np.asarray(orders)
    ls = np.asarray(degrees)
    if ns.ndim != 1 or ls.ndim != 1 or ns.size == 0 or ls.size == 0:
        raise ValueError("orders and degrees must be nonempty one-dimensional arrays")
    if np.any(ns != np.floor(ns)) or np.any(ns < 1) or np.any(ls != np.floor(ls)) or np.any((ls < 0) | (ls > 3)):
        raise ValueError("radial orders must be positive integers and degrees integers 0..3")
    n, l = np.meshgrid(ns.astype(int), ls.astype(int), indexing="ij")
    nu = delta_nu * (n + l / 2 + epsilon) - l * (l + 1) * d0
    if not np.all(np.isfinite(nu)) or np.any(nu <= 0):
        raise ValueError("nonphysical frequencies; choose appropriate high radial orders")
    out = np.empty(n.size, dtype=[("n", "i4"), ("l", "i4"), ("nu_microhz", "f8")])
    out["n"], out["l"], out["nu_microhz"] = n.ravel(), l.ravel(), nu.ravel()
    return out


def scaling_with_covariance(nu_max: float, delta_nu: float, teff: float,
                            input_covariance: np.ndarray) -> ScalingResult:
    """Mass/radius scaling with full correlated 3x3 input covariance.

    Inputs are ordered (nu_max [microHz], delta_nu [microHz], Teff [K]).
    Propagation is first order in log observables; only suitable for
    small relative uncertainties and positive measurements.
    """
    x = np.asarray([nu_max, delta_nu, teff], dtype=float)
    cov = np.asarray(input_covariance, dtype=float)
    if not np.all(np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("observables must be positive and finite")
    if cov.shape != (3, 3) or not np.all(np.isfinite(cov)) or not np.allclose(cov, cov.T, rtol=1e-12, atol=1e-12):
        raise ValueError("input covariance must be a finite symmetric 3x3 matrix")
    if np.linalg.eigvalsh(cov).min() < -1e-10 * max(1., np.max(np.abs(cov))):
        raise ValueError("input covariance must be positive semidefinite")
    m = (nu_max / NU_MAX_SUN)**3 * (delta_nu / DELTA_NU_SUN)**-4 * (teff / TEFF_SUN)**1.5
    r = (nu_max / NU_MAX_SUN) * (delta_nu / DELTA_NU_SUN)**-2 * (teff / TEFF_SUN)**0.5
    exponent = np.array([[3., -4., 1.5], [1., -2., 0.5]])
    jacobian = np.diag([m, r]) @ exponent @ np.diag(1 / x)
    out_cov = jacobian @ cov @ jacobian.T
    return ScalingResult(float(m), float(r), out_cov)


def estimate_large_separation(frequency_microhz, power, expected_microhz: float,
                              relative_window: float = 0.25) -> dict:
    """Find the strongest normalized comb autocorrelation near a prior Delta-nu.

    Requires uniformly sampled frequency bins; the expected separation
    is a *prior* used to avoid half-spacing aliases. Reports no uncertainty.
    The spectrum is synthetic or externally supplied; no observation is fetched.
    """
    f = np.asarray(frequency_microhz, dtype=float)
    p = np.asarray(power, dtype=float)
    if f.ndim != 1 or p.shape != f.shape or f.size < 32:
        raise ValueError("frequency and power must be equal-length 1D arrays with >=32 bins")
    step = np.diff(f)
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(p)) or np.any(p < 0) or np.any(step <= 0):
        raise ValueError("finite ascending frequencies and nonnegative finite power required")
    if not np.allclose(step, step[0], rtol=1e-7, atol=1e-10):
        raise ValueError("frequency grid must be uniform")
    if not (np.isfinite(expected_microhz) and expected_microhz > 0 and 0 < relative_window < 0.5):
        raise ValueError("invalid expected separation or window")
    lo = max(1, int(np.ceil(expected_microhz * (1 - relative_window) / step[0])))
    hi = min(len(f) - 3, int(np.floor(expected_microhz * (1 + relative_window) / step[0])))
    if lo > hi:
        raise ValueError("requested separation is outside the sampled frequency range")
    q = p - np.median(p)
    if np.max(q) <= 0 or np.linalg.norm(q) == 0:
        raise ValueError("spectrum has no positive contrast above its median")
    lags = np.arange(lo, hi + 1)
    scores = np.array([np.dot(q[:-k], q[k:]) / (np.linalg.norm(q[:-k]) * np.linalg.norm(q[k:]))
                       for k in lags])
    best = int(np.argmax(scores))
    return {"delta_nu_microhz": float(lags[best] * step[0]),
            "correlation": float(scores[best]),
            "frequency_bin_microhz": float(step[0])}
