import numpy as np
import pytest

from research.seismic_diagnostics import (
    DNU_SUN, NUMAX_SUN, TEFF_SUN, asymptotic_modes,
    log_covariance_propagation, mass_radius, sample_scaling,
    small_separation_02, temperature_omission_bias,
)


def test_solar_reference_is_one():
    assert mass_radius(NUMAX_SUN, DNU_SUN, TEFF_SUN) == pytest.approx((1, 1))


def test_temperature_exponents():
    a = mass_radius(2100, 96, TEFF_SUN)
    b = mass_radius(2100, 96, TEFF_SUN * 1.04)
    assert b[0] / a[0] == pytest.approx(1.04**1.5)
    assert b[1] / a[1] == pytest.approx(1.04**0.5)


def test_browser_proxy_equals_solar_temperature_case():
    d = temperature_omission_bias(2089, 95.5, TEFF_SUN)
    assert d["fractional_mass_proxy_bias"] == pytest.approx(0)
    assert d["fractional_radius_proxy_bias"] == pytest.approx(0)


def test_covariance_exact_for_independent_terms():
    s = np.array([0.01, 0.02, 0.03])
    r = log_covariance_propagation((2089, 95.5, 5846), np.diag(s**2))
    expected = 9 * s[0]**2 + 16 * s[1]**2 + 2.25 * s[2]**2
    assert r["log_covariance"][0][0] == pytest.approx(expected)


def test_correlated_mc_matches_log_analytic():
    s = np.array([0.01, 0.003, 0.015])
    corr = np.array([[1, 0.35, 0], [0.35, 1, -0.15], [0, -0.15, 1]])
    cov = np.outer(s, s) * corr
    analytic = np.asarray(log_covariance_propagation((2089, 95.5, 5846), cov)["log_covariance"])
    draws = sample_scaling((2089, 95.5, 5846), cov, samples=30000, seed=11)
    observed = np.cov(np.log(draws).T)
    assert np.allclose(observed, analytic, rtol=0.05, atol=2e-5)


def test_mc_seed_reproducible():
    cov = np.diag([0.01, 0.01, 0.01]) ** 2
    a = sample_scaling((2089, 95.5, 5846), cov, samples=10, seed=7)
    b = sample_scaling((2089, 95.5, 5846), cov, samples=10, seed=7)
    assert np.array_equal(a, b)


def test_asymptotic_small_separation_exact():
    assert small_separation_02(95.5, 1.2, 1.3) == pytest.approx(7.8)


def test_mode_order_translation_by_delta_nu():
    a = asymptotic_modes(95.5, 1.2, 1.3, orders=[19])
    b = asymptotic_modes(95.5, 1.2, 1.3, orders=[20])
    for x, y in zip(a, b):
        assert y["frequency_microhz"] - x["frequency_microhz"] == pytest.approx(95.5)
        assert y["echelle_phase_microhz"] == pytest.approx(x["echelle_phase_microhz"])


@pytest.mark.parametrize("values", [(0, 100, 5772), (2000, -1, 5772),
                                     (2000, 100, float("nan"))])
def test_reject_invalid_inputs(values):
    with pytest.raises(ValueError):
        mass_radius(*values)


def test_reject_indefinite_covariance():
    with pytest.raises(ValueError):
        log_covariance_propagation((2089, 95.5, 5846),
                                   np.array([[1, 2, 0], [2, 1, 0], [0, 0, 1]]))


def test_reject_invalid_mode_degree():
    with pytest.raises(ValueError):
        asymptotic_modes(95.5, 1.2, 1.3, degrees=(4,))
