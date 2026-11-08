import math
import unittest
from robustness import scaling, log_uncertainty, window_response

class SeismicTests(unittest.TestCase):
    def test_solar(self):
        self.assertEqual(scaling(3090,135.1,5772),(1,1))
    def test_covariance(self):
        cov=[[30**2,0,0],[0,1.0,0],[0,0,70**2]]
        result=log_uncertainty(3090,135.1,5772,cov)
        expected=(3*30/3090)**2+(4/135.1)**2+(1.5*70/5772)**2
        self.assertAlmostEqual(result[0][0],expected)
    def test_correlated_covariance(self):
        cov=[[900,15,0],[15,1,0],[0,0,4900]]
        self.assertLess(log_uncertainty(3090,135.1,5772,cov)[0][0],log_uncertainty(3090,135.1,5772,[[900,0,0],[0,1,0],[0,0,4900]])[0][0])
    def test_alias(self):
        times=[i/24 for i in range(72)]
        result=window_response(times,[0,1e6/3600,1e6/7200])
        self.assertAlmostEqual(result[0],1)
        self.assertAlmostEqual(result[1],1)
        self.assertLess(result[2],1e-20)
    def test_invalid(self):
        with self.assertRaises(ValueError): scaling(0,135,5772)

if __name__=='__main__': unittest.main()
