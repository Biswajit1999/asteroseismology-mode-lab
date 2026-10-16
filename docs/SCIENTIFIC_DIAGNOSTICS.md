# Asteroseismic diagnostics (Python, independent of browser UI)

These functions supplement the browser's illustrative p-mode spectrum. They do **not** fit real Kepler/TESS light curves and do not derive stellar ages or evolutionary states. No observational frequency lists are bundled or inferred from the visual reference anchors.

## Reproduce

```bash
python -m pip install -r requirements-science.txt
python -m unittest discover -s tests -v
```

### Three scientific analyses

1. **Asymptotic p-mode comb** — `asymptotic_modes(orders, delta_nu, epsilon, d0)` evaluates `nu(n,l) = DeltaNu(n+l/2+epsilon)-l(l+1)D0` for degrees 0–3. The unit test verifies the `l=0/2` small separation `nu(n,0)-nu(n-1,2)=6D0`. It is a high-radial-order approximation, not a stellar-structure eigenfrequency solver.
2. **Correlated scaling-relation uncertainty** — `scaling_with_covariance` evaluates solar-normalized mass and radius using `nu_max`, `DeltaNu`, and `Teff`. It propagates the **full** 3x3 input covariance through an analytic Jacobian, returning a 2x2 mass-radius covariance. It assumes small errors, and excludes known scaling-relation systematics and corrections (e.g. `f_DeltaNu`). The seeded Monte Carlo test validates the linear approximation in the small-error regime.
3. **Comb spacing recovery** — `estimate_large_separation` computes a normalized autocorrelation of a uniformly sampled nonnegative power spectrum in a prior-limited window. It is validated on a deterministic synthetic 95.5 microHz comb. Half-spacing aliases, mode linewidth, envelope curvature, finite resolution and stochastic excitation can bias real-data estimates; this is **not** a measurement of KIC 5184732.

### Units and provenance

Frequencies are microHz, temperatures kelvin, output mass/radius in solar units. Constants used: `nu_max,sun=3090 microHz`, `DeltaNu_sun=135.1 microHz`, `Teff_sun=5772 K`. Reference: Kjeldsen & Bedding (1995), *A&A* 293, 87; Chaplin & Miglio (2013), *ARA&A* 51, 353. The seeded tests use **synthetic** data only. Nothing in these tests confirms the UI's simplified inclination model or its radius-only evolutionary-state label.
