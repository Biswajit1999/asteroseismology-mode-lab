"""Independent, dependency-free diagnostics for the illustrative p-mode lab.

The global asteroseismic scaling relations are approximate, especially for
non-solar stars; none of the routines is a stellar-evolution inference engine.
"""
from __future__ import annotations

import math
import random
from typing import Sequence

NUMAX_SUN_UHZ = 3090.0
DNU_SUN_UHZ = 135.1
TEFF_SUN_K = 5772.0


def _check_positive(*values: float) -> None:
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError("all physical inputs must be positive and finite")


def scaling(nu_max_uhz: float, dnu_uhz: float, teff_k: float) -> tuple[float, float]:
    """Return dimensionless M/Msun and R/Rsun from uncorrected scaling laws."""
    _check_positive(nu_max_uhz, dnu_uhz, teff_k)
    n, d, t = nu_max_uhz / NUMAX_SUN_UHZ, dnu_uhz / DNU_SUN_UHZ, teff_k / TEFF_SUN_K
    return n**3 * d**-4 * t**1.5, n * d**-2 * math.sqrt(t)


def _cholesky(a: Sequence[Sequence[float]]) -> list[list[float]]:
    n = len(a)
    if n != 3 or any(len(row) != 3 for row in a):
        raise ValueError("covariance must be 3x3 in [nu_max, dnu, teff] order")
    if any(not math.isfinite(v) for row in a for v in row):
        raise ValueError("nonfinite covariance")
    if any(abs(a[i][j] - a[j][i]) > 1e-9 * max(1.0, abs(a[i][j]), abs(a[j][i])) for i in range(n) for j in range(n)):
        raise ValueError("covariance must be symmetric")
    l = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = a[i][j] - sum(l[i][k] * l[j][k] for k in range(j))
            if i == j:
                if s <= 0:
                    raise ValueError("covariance must be positive definite")
                l[i][j] = math.sqrt(s)
            else:
                l[i][j] = s / l[j][j]
    return l


def propagate_scaling(nu_max_uhz: float, dnu_uhz: float, teff_k: float,
                      covariance: Sequence[Sequence[float]]) -> dict:
    """Delta-method covariance, retaining input correlations.

    Valid for small uncertainties; systematic corrections to seismic scaling
    relations are *not* included.
    """
    m, r = scaling(nu_max_uhz, dnu_uhz, teff_k)
    _cholesky(covariance)
    jac = ([3*m/nu_max_uhz, -4*m/dnu_uhz, 1.5*m/teff_k],
           [r/nu_max_uhz, -2*r/dnu_uhz, 0.5*r/teff_k])
    cov_out = [[sum(jac[i][k] * covariance[k][l] * jac[j][l]
                    for k in range(3) for l in range(3)) for j in range(2)] for i in range(2)]
    return {"mass_msun": m, "radius_rsun": r, "sigma_mass_msun": math.sqrt(cov_out[0][0]),
            "sigma_radius_rsun": math.sqrt(cov_out[1][1]), "covariance_m_r": cov_out,
            "corr_m_r": cov_out[0][1]/math.sqrt(cov_out[0][0]*cov_out[1][1])}


def monte_carlo_scaling(nu_max_uhz: float, dnu_uhz: float, teff_k: float,
                        covariance: Sequence[Sequence[float]], trials: int = 5000,
                        seed: int = 20261015) -> dict:
    """Reproducible Gaussian measurement-error validation of the delta method."""
    if trials < 2:
        raise ValueError("at least two trials are required")
    l = _cholesky(covariance)
    rng = random.Random(seed)
    values = []
    for _ in range(trials):
        z = [rng.gauss(0, 1) for _ in range(3)]
        obs = [v + sum(l[i][j] * z[j] for j in range(i+1))
               for i, v in enumerate((nu_max_uhz, dnu_uhz, teff_k))]
        if min(obs) <= 0:
            raise ValueError("uncertainties too large: sampled nonphysical measurement")
        values.append(scaling(*obs))
    means = [sum(v[k] for v in values)/trials for k in range(2)]
    variances = [sum((v[k]-means[k])**2 for v in values)/(trials-1) for k in range(2)]
    return {"trials": trials, "seed": seed, "mean_mass_msun": means[0],
            "mean_radius_rsun": means[1], "sigma_mass_msun": math.sqrt(variances[0]),
            "sigma_radius_rsun": math.sqrt(variances[1])}


def _solve_3x3(a: list[list[float]], b: list[float]) -> list[float]:
    """Pivoted Gaussian elimination for a small weighted normal system."""
    m = [list(row) + [v] for row, v in zip(a, b)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda k: abs(m[k][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise ValueError("mode design is singular; vary radial orders and degrees")
        m[col], m[pivot] = m[pivot], m[col]
        factor = m[col][col]
        for j in range(col, 4):
            m[col][j] /= factor
        for i in range(3):
            if i != col:
                factor = m[i][col]
                for j in range(col, 4):
                    m[i][j] -= factor*m[col][j]
    return [m[i][3] for i in range(3)]


def fit_asymptotic_modes(modes: Sequence[tuple[int, int, float, float]]) -> dict:
    """Fit nu_nl = Dnu(n+l/2+epsilon)-l(l+1)D0 with known mode IDs.

    Inputs: (radial order n, degree l, frequency in microHz, sigma in microHz).
    The model omits glitches, curvature and rotational splitting.
    """
    if len(modes) < 4:
        raise ValueError("need at least four identified modes")
    a = [[0.0]*3 for _ in range(3)]
    b = [0.0]*3
    for n, degree, frequency, sigma in modes:
        if not isinstance(n, int) or not isinstance(degree, int) or degree < 0 or degree > 3:
            raise ValueError("invalid radial order or degree")
        _check_positive(frequency, sigma)
        x = [n + degree/2, 1.0, -degree*(degree+1)]
        w = 1.0/sigma**2
        for i in range(3):
            b[i] += w*x[i]*frequency
            for j in range(3):
                a[i][j] += w*x[i]*x[j]
    dnu, intercept, d0 = _solve_3x3(a, b)
    if dnu <= 0:
        raise ValueError("fitted large separation is nonphysical")
    chi2 = sum(((frequency - (dnu*(n+degree/2)+intercept-degree*(degree+1)*d0))/sigma)**2
               for n, degree, frequency, sigma in modes)
    return {"delta_nu_uhz": dnu, "epsilon": intercept/dnu, "D0_uhz": d0,
            "chi2": chi2, "dof": len(modes)-3, "reduced_chi2": chi2/(len(modes)-3)}


def observing_window(baseline_days: float, cadence_seconds: float, nu_max_uhz: float,
                     dnu_uhz: float, linewidth_uhz: float) -> dict:
    """Rectangular-window Rayleigh scale and uniform-sampling Nyquist check.

    Not a TESS/Kepler spacecraft window-function simulation; gaps, apodization,
    and actual instrumental transfer functions are excluded.
    """
    _check_positive(baseline_days, cadence_seconds, nu_max_uhz, dnu_uhz, linewidth_uhz)
    if baseline_days*86400 < 2*cadence_seconds:
        raise ValueError("baseline must contain at least two cadences")
    rayleigh = 1e6/(baseline_days*86400)
    nyquist = 1e6/(2*cadence_seconds)
    folded = abs((nu_max_uhz + nyquist) % (2*nyquist) - nyquist)
    return {"rayleigh_uhz": rayleigh, "nyquist_uhz": nyquist,
            "nu_max_alias_uhz": folded, "nu_max_below_nyquist": nu_max_uhz < nyquist,
            "linewidth_resolved_heuristic": linewidth_uhz >= rayleigh,
            "large_separation_resolved": dnu_uhz >= rayleigh,
            "n_independent_frequency_bins_per_dnu": dnu_uhz/rayleigh}
