"""Reference-anchored asteroseismic pilot with clearly illustrative uncertainties."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from research.seismic_diagnostics import (  # noqa: E402
    asymptotic_modes, log_covariance_propagation, sample_scaling,
    small_separation_02, temperature_omission_bias,
)


def main() -> None:
    reference = json.loads((ROOT / "data/reference.json").read_text(encoding="utf-8"))
    params = reference["measured_parameters"]
    observables = (float(params["numax_microHz"]), float(params["Dnu_microHz"]),
                   float(params["Teff_K"]))
    # Hypothetical, not catalogued: 1% nu_max, 0.3% Delta_nu, 1.5% Teff;
    # corr(nu_max, Delta_nu)=0.35, corr(Delta_nu, Teff)=-0.15.
    s = np.array([0.01, 0.003, 0.015])
    corr = np.array([[1, 0.35, 0], [0.35, 1, -0.15], [0, -0.15, 1]])
    cov = np.outer(s, s) * corr
    analytic = log_covariance_propagation(observables, cov)
    draws = sample_scaling(observables, cov, samples=30000, seed=20261026)
    modes = asymptotic_modes(observables[1], epsilon=1.2, d0=1.3)
    report = {
        "reference_object": params["name"],
        "source": reference["source_url"],
        "measured_inputs": {"numax_microHz": observables[0],
                            "delta_nu_microHz": observables[1], "teff_K": observables[2]},
        "temperature_bias": temperature_omission_bias(*observables),
        "uncertainty_assumptions": "illustrative lognormal 1%, 0.3%, 1.5%; correlations 0.35 and -0.15; not LEGACY errors",
        "analytic_log_covariance": analytic["log_covariance"],
        "mc_log_covariance": np.cov(np.log(draws).T).tolist(),
        "mc_draws": len(draws),
        "asymptotic_assumptions": "epsilon=1.2, D0=1.3 microHz, not fitted to KIC 5184732",
        "small_separation_02_microHz": small_separation_02(observables[1], 1.2, 1.3),
        "mode_count": len(modes),
        "limitations": "global scaling only; no stellar-model corrections, observed mode fitting, or age determination",
    }
    out = ROOT / "results/seismic_pilot.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
