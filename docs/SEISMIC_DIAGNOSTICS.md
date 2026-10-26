# Independent asteroseismic scaling and mode diagnostics (2026-10-26)

This optional Python research layer uses the **published global inputs** for
16 Cygni A (KIC 5184732) already pinned in `data/reference.json`: νmax = 2089
μHz, Δν = 95.5 μHz and Teff = 5846 K. The five envelope points in that JSON
are **illustrative**, not measured individual oscillation-mode frequencies.

1. **Temperature-aware global mass/radius scaling** implements the standard
   Kjeldsen–Bedding exponents and compares against the browser's implicit
   solar-Teff approximation. Solar anchors match the browser (3090 and 135.1
   μHz); Teff,Sun = 5772 K. This is not a calibrated stellar-model mass.
2. **Covariance propagation** uses the exact log-linear Jacobian for M and R,
   including correlations between observables, and verifies it with a seeded
   lognormal Monte Carlo. The 1%, 0.3%, 1.5% input uncertainties and the two
   correlations in the pilot are **hypothetical**, not errors reported by Lund
   et al. (2017). Solar-reference uncertainties are held fixed.
3. **Asymptotic p-mode ridges** implement ν(n,ℓ) ≈ Δν(n+ℓ/2+ε)−ℓ(ℓ+1)D0.
   The identity δν02=6D0 is tested, and frequencies can be folded modulo Δν
   for an échelle diagram. ε=1.2 and D0=1.3 μHz are demonstration choices,
   not a fit to observed individual modes. This is distinct from the browser's
   hand-tuned Gaussian-comb visualisation.

Run `python scripts/run_seismic_diagnostics.py` with NumPy installed. It reads
`data/reference.json` and writes `results/seismic_pilot.json`. Run
`python -m pytest -q tests/test_seismic_diagnostics.py` for regression checks.

No new Kepler light curves were downloaded, no stellar age was inferred, and
no claim of independent asteroseismic parameter measurement is made.

Install the optional research environment using
`python -m pip install -r requirements-science.txt`.
`python scripts/plot_seismic_diagnostics.py` creates two labelled figures
(illustrative asymptotic ridges and hypothetical correlated mass-radius errors).
The committed `results/seismic_pilot.json` was generated using the pinned
reference inputs and seed 20261026.
