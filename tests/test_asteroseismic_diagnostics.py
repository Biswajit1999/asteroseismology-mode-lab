"""Run with python -m unittest discover -s tests -v."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "science"))
from asteroseismic_diagnostics import (asymptotic_modes, scaling_with_covariance,
                                       estimate_large_separation, NU_MAX_SUN, DELTA_NU_SUN, TEFF_SUN)


class AsteroseismicDiagnosticsTest(unittest.TestCase):
    def test_small_separation_six_d0(self):
        modes = asymptotic_modes(range(12, 23), 100, 1.2, 1.4)
        for n in range(13, 23):
            radial = modes[(modes["n"] == n) & (modes["l"] == 0)]["nu_microhz"][0]
            quadrupole = modes[(modes["n"] == n - 1) & (modes["l"] == 2)]["nu_microhz"][0]
            self.assertAlmostEqual(radial - quadrupole, 8.4, places=9)

    def test_radial_large_separation(self):
        a = asymptotic_modes([18, 19], 95.5, 1.0, 1.0, degrees=[0])
        self.assertAlmostEqual(a["nu_microhz"][1] - a["nu_microhz"][0], 95.5)

    def test_bad_mode_parameters(self):
        with self.assertRaises(ValueError):
            asymptotic_modes([0], 100, 1, 1)

    def test_solar_scaling(self):
        r = scaling_with_covariance(NU_MAX_SUN, DELTA_NU_SUN, TEFF_SUN, np.zeros((3, 3)))
        self.assertAlmostEqual(r.mass_solar, 1)
        self.assertAlmostEqual(r.radius_solar, 1)
        np.testing.assert_allclose(r.covariance, np.zeros((2, 2)), atol=1e-16)

    def test_scaling_covariance_monte_carlo(self):
        x = np.array([2089., 95.5, 5800.])
        cov = np.diag([8.**2, .15**2, 25.**2])
        result = scaling_with_covariance(*x, cov)
        rng = np.random.default_rng(20261016)
        s = rng.multivariate_normal(x, cov, size=20000)
        mass = (s[:, 0] / NU_MAX_SUN)**3 * (s[:, 1] / DELTA_NU_SUN)**-4 * (s[:, 2] / TEFF_SUN)**1.5
        radius = (s[:, 0] / NU_MAX_SUN) * (s[:, 1] / DELTA_NU_SUN)**-2 * (s[:, 2] / TEFF_SUN)**.5
        np.testing.assert_allclose(np.cov([mass, radius]), result.covariance, rtol=.06)

    def test_correlated_covariance(self):
        x = [3000., 140., 5700.]
        cov = np.diag([30.**2, 1.**2, 80.**2])
        cov[0, 1] = cov[1, 0] = 15.
        a = scaling_with_covariance(*x, cov)
        b = scaling_with_covariance(*x, np.diag(np.diag(cov)))
        self.assertNotAlmostEqual(a.covariance[0, 0], b.covariance[0, 0])

    def test_invalid_covariance(self):
        with self.assertRaises(ValueError):
            scaling_with_covariance(3090, 135.1, 5772, np.diag([-1., 1., 1.]))

    def test_recover_synthetic_comb(self):
        f = np.arange(1000., 4000., .25)
        p = .01 + sum(np.exp(-.5 * ((f - (1300 + n * 95.5)) / 1.6)**2)
                        for n in range(25))
        r = estimate_large_separation(f, p, expected_microhz=95.5)
        self.assertAlmostEqual(r["delta_nu_microhz"], 95.5, delta=.25)
        self.assertGreater(r["correlation"], .9)

    def test_reject_nonuniform_grid(self):
        with self.assertRaises(ValueError):
            estimate_large_separation(np.arange(64, dtype=float)**1.1,
                                      np.ones(64), 5)


if __name__ == "__main__":
    unittest.main()
