import math
import unittest
from analysis.seismic_diagnostics import (scaling, propagate_scaling, monte_carlo_scaling,
                                          fit_asymptotic_modes, observing_window)

class SeismicTests(unittest.TestCase):
    def test_solar_normalization(self):
        self.assertEqual(scaling(3090,135.1,5772),(1.0,1.0))
    def test_scaling_radius_mass(self):
        m,r=scaling(3090,135.1/2**1.5,5772)
        self.assertAlmostEqual(m,64.0,9)
        self.assertAlmostEqual(r,8.0,9)
    def test_covariance_propagation_diagonal(self):
        nu,d,t=2089,95.5,5846
        s=[0.01*nu,0.003*d,0.015*t]
        result=propagate_scaling(nu,d,t,[[s[i]**2 if i==j else 0 for j in range(3)] for i in range(3)])
        self.assertAlmostEqual(result['sigma_mass_msun']/result['mass_msun'],
                               math.sqrt((3*.01)**2+(4*.003)**2+(1.5*.015)**2),12)
    def test_covariance_correlation_matters(self):
        nu,d,t=2089,95.5,5846
        s=[nu*.01,d*.003,t*.015]
        c=[[s[i]*s[j]*(1 if i==j else (0.4 if {i,j}=={0,1} else 0)) for j in range(3)] for i in range(3)]
        corr=propagate_scaling(nu,d,t,c)
        diag=propagate_scaling(nu,d,t,[[s[i]**2 if i==j else 0 for j in range(3)] for i in range(3)])
        self.assertLess(corr['sigma_mass_msun'],diag['sigma_mass_msun'])
    def test_mc_agrees_with_linear_error(self):
        nu,d,t=2089,95.5,5846
        s=[nu*.004,d*.002,t*.006]
        c=[[s[i]**2 if i==j else 0 for j in range(3)] for i in range(3)]
        analytic=propagate_scaling(nu,d,t,c)
        mc=monte_carlo_scaling(nu,d,t,c,trials=4000)
        self.assertLess(abs(mc['sigma_mass_msun']/analytic['sigma_mass_msun']-1),.08)
        self.assertLess(abs(mc['sigma_radius_rsun']/analytic['sigma_radius_rsun']-1),.08)
    def test_reject_invalid_covariance(self):
        with self.assertRaises(ValueError):
            propagate_scaling(2089,95.5,5846,[[1,2,0],[2,1,0],[0,0,1]])
    def test_fit_recovers_exact_comb(self):
        modes=[(n,l,95.5*(n+l/2+1.2)-l*(l+1)*1.5,.05)
               for n in range(18,25) for l in (0,1,2)]
        fit=fit_asymptotic_modes(modes)
        self.assertAlmostEqual(fit['delta_nu_uhz'],95.5,8)
        self.assertAlmostEqual(fit['epsilon'],1.2,8)
        self.assertAlmostEqual(fit['D0_uhz'],1.5,8)
        self.assertAlmostEqual(fit['chi2'],0,7)
    def test_fit_rejects_degenerate_design(self):
        with self.assertRaises(ValueError):
            fit_asymptotic_modes([(20,0,2100,.1)]*5)
    def test_window_alias(self):
        r=observing_window(27.4,1800,2089,95.5,1.0)
        self.assertFalse(r['nu_max_below_nyquist'])
        self.assertAlmostEqual(r['nyquist_uhz'],277.7777777778,6)
        self.assertLess(r['nu_max_alias_uhz'],r['nyquist_uhz'])
    def test_window_high_cadence(self):
        r=observing_window(90,120,2089,95.5,1.0)
        self.assertTrue(r['nu_max_below_nyquist'])
        self.assertTrue(r['linewidth_resolved_heuristic'])
    def test_bad_inputs(self):
        with self.assertRaises(ValueError): observing_window(0,120,2089,95.5,1)
        with self.assertRaises(ValueError): scaling(3090,-1,5772)

if __name__ == '__main__': unittest.main()
