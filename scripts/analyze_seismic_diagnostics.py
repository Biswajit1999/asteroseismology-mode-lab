"""Run with: python -m scripts.analyze_seismic_diagnostics [path/to/reference.json]"""
import json
import math
import random
import sys
from pathlib import Path
from analysis.seismic_diagnostics import (propagate_scaling, monte_carlo_scaling,
                                          fit_asymptotic_modes, observing_window)


def run(reference: Path) -> dict:
    measured=json.loads(reference.read_text())['measured_parameters']
    nu,d,t=measured['numax_microHz'],measured['Dnu_microHz'],measured['Teff_K']
    # Hypothetical measurement errors for demonstration; the reference does not
    # provide a full published covariance matrix.
    sigma=[nu*.01,d*.003,t*.015]
    corr=[[1,.25,0],[.25,1,0],[0,0,1]]
    cov=[[sigma[i]*sigma[j]*corr[i][j] for j in range(3)] for i in range(3)]
    rng=random.Random(20261015)
    # Synthetic identified mode frequencies: NOT measured KIC 5184732 frequencies.
    modes=[(n,l,d*(n+l/2+1.2)-l*(l+1)*1.5+rng.gauss(0,.08),.08)
           for n in range(18,25) for l in (0,1,2)]
    return {'star':measured['name'], 'data_provenance':'real global parameters; synthetic modes',
            'assumed_errors':{'fractional_numax':.01,'fractional_dnu':.003,
                              'fractional_teff':.015,'corr_numax_dnu':.25},
            'scaling_delta_method':propagate_scaling(nu,d,t,cov),
            'scaling_mc':monte_carlo_scaling(nu,d,t,cov,trials=5000),
            'synthetic_mode_fit':fit_asymptotic_modes(modes),
            'windows':{'90d_120s':observing_window(90,120,nu,d,1.0),
                       '27d_120s':observing_window(27.4,120,nu,d,1.0),
                       '27d_1800s':observing_window(27.4,1800,nu,d,1.0)}}

if __name__=='__main__':
    source=Path(sys.argv[1]) if len(sys.argv)>1 else Path('data/reference.json')
    print(json.dumps(run(source),indent=2))
