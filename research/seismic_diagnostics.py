"""Independent, temperature-aware solar-like oscillation scaling diagnostics.

Reference νmax and Δν values are measured; uncertainties in the pilot are
illustrative and are *not* attributed to the Kepler LEGACY catalogue.
"""

from __future__ import annotations

import numpy as np

NUMAX_SUN = 3090.0  # microHz; matches the browser's adopted solar anchor
DNU_SUN = 135.1  # microHz; matches the browser's adopted solar anchor
TEFF_SUN = 5772.0  # K; IAU nominal solar effective temperature

# Jacobian of [ln(M/Msun), ln(R/Rsun)] w.r.t. [ln(nu_max), ln(Delta_nu), ln(Teff)].
LOG_JACOBIAN = np.array([[3.0, -4.0, 1.5], [1.0, -2.0, 0.5]])


def mass_radius(numax: float, delta_nu: float, teff: float,
                *, solar_numax: float = NUMAX_SUN,
                solar_delta_nu: float = DNU_SUN,
                solar_teff: float = TEFF_SUN) -> tuple[float, float]:
    """Kjeldsen-Bedding global scaling, with no stellar-model correction."""
    values = np.array([numax, delta_nu, teff, solar_numax, solar_delta_nu,
                       solar_teff], dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("all frequencies and temperatures must be positive and finite")
    x = numax / solar_numax
    y = delta_nu / solar_delta_nu
    t = teff / solar_teff
    return float(x**3 * y**-4 * t**1.5), float(x * y**-2 * t**0.5)


def temperature_omission_bias(numax: float, delta_nu: float,
                              teff: float) -> dict[str, float]:
    """Compare the full scaling to the browser's solar-Teff proxy."""
    mass, radius = mass_radius(numax, delta_nu, teff)
    old_mass, old_radius = mass_radius(numax, delta_nu, TEFF_SUN)
    return {"mass_with_teff": mass, "radius_with_teff": radius,
            "mass_solar_teff_proxy": old_mass, "radius_solar_teff_proxy": old_radius,
            "fractional_mass_proxy_bias": old_mass / mass - 1,
            "fractional_radius_proxy_bias": old_radius / radius - 1}


def log_covariance_propagation(observables: tuple[float, float, float],
                               covariance: np.ndarray) -> dict[str, object]:
    """Exact log-linear covariance propagation for fixed solar reference values.

    `covariance` is the 3x3 covariance of natural logarithms of the three
    observables. It must be symmetric positive semidefinite. Returns the
    covariance of [ln M, ln R], not a Gaussian approximation in M and R.
    """
    mass, radius = mass_radius(*observables)
    cov = np.asarray(covariance, dtype=float)
    if cov.shape != (3, 3) or not np.all(np.isfinite(cov)):
        raise ValueError("expected a finite 3x3 covariance")
    if not np.allclose(cov, cov.T, rtol=0, atol=1e-12):
        raise ValueError("covariance must be symmetric")
    cov = (cov + cov.T) / 2
    if np.linalg.eigvalsh(cov).min() < -1e-12:
        raise ValueError("covariance must be positive semidefinite")
    propagated = LOG_JACOBIAN @ cov @ LOG_JACOBIAN.T
    return {"mass": mass, "radius": radius,
            "log_covariance": propagated.tolist(),
            "fractional_log_sigma_mass": float(np.sqrt(max(0, propagated[0, 0]))),
            "fractional_log_sigma_radius": float(np.sqrt(max(0, propagated[1, 1])))}


def sample_scaling(observables: tuple[float, float, float], covariance: np.ndarray,
                   *, samples: int = 20000, seed: int = 20261026) -> np.ndarray:
    """Reproducible correlated lognormal Monte Carlo; columns are mass, radius."""
    if samples < 2:
        raise ValueError("samples must be >= 2")
    log_covariance_propagation(observables, covariance)  # validate inputs
    rng = np.random.default_rng(seed)
    log_inputs = rng.multivariate_normal(np.log(observables), covariance,
                                         size=samples, check_valid="raise")
    log_ratios = log_inputs - np.log([NUMAX_SUN, DNU_SUN, TEFF_SUN])
    return np.exp(log_ratios @ LOG_JACOBIAN.T)


def asymptotic_modes(delta_nu: float, epsilon: float, d0: float,
                     *, orders=range(15, 27), degrees=(0, 1, 2)) -> list[dict]:
    """Tassoul high-order asymptotic p-mode ridge positions (microHz).

    Does not model mode excitation, linewidths, rotation, or mixed modes.
    """
    if not np.isfinite([delta_nu, epsilon, d0]).all() or delta_nu <= 0 or d0 < 0:
        raise ValueError("require finite delta_nu > 0, epsilon and d0 >= 0")
    result = []
    for n in orders:
        if not isinstance(n, (int, np.integer)) or n < 0:
            raise ValueError("orders must be nonnegative integers")
        for ell in degrees:
            if ell not in (0, 1, 2, 3):
                raise ValueError("degrees must be 0, 1, 2 or 3")
            frequency = delta_nu * (n + ell / 2 + epsilon) - ell * (ell + 1) * d0
            if frequency <= 0:
                raise ValueError("mode frequency must be positive")
            result.append({"n": int(n), "ell": int(ell), "frequency_microhz": float(frequency),
                           "echelle_phase_microhz": float(frequency % delta_nu)})
    return result


def small_separation_02(delta_nu: float, epsilon: float, d0: float,
                        *, n: int = 20) -> float:
    """nu(n,0)-nu(n-1,2) = 6 D0 under the adopted asymptotic formula."""
    modes = asymptotic_modes(delta_nu, epsilon, d0,
                            orders=[n - 1, n], degrees=(0, 2))
    return (next(m["frequency_microhz"] for m in modes if m["n"] == n and m["ell"] == 0)
            - next(m["frequency_microhz"] for m in modes
                   if m["n"] == n - 1 and m["ell"] == 2))
