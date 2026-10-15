# Reproducible seismic diagnostics (Python standard library)

Run from the repository root:

```sh
python -m unittest discover -s tests -p 'test_seismic_diagnostics.py' -v
python -m scripts.analyze_seismic_diagnostics > seismic_diagnostics.json
```

Three independent analyses supplement (and do not alter) the existing illustrative browser model:

1. **Correlated-error propagation for global seismic scaling relations.** Compute mass/radius estimates and their full 2×2 covariance from a 3×3 covariance of (`nu_max`, `Delta_nu`, `Teff`). Analytic Jacobian (delta method) is cross-checked with seeded Gaussian Monte Carlo. This quantifies *measurement* error only: systematic corrections to Δν scaling, metallicity and stellar evolution are outside scope.
2. **Weighted fit to an identified asymptotic mode comb.** Estimate Δν, phase ε and D₀ from input tuples `(radial_order, angular_degree, frequency_µHz, sigma_µHz)`. The demo mode list is **synthetic**, not measured Kepler individual-mode data. It omits rotational splitting, glitches and mixed modes; mode identification is assumed known.
3. **Observing-baseline and cadence limits.** Compute Rayleigh resolution `1/T`, uniform-sampling Nyquist `1/(2Δt)`, folded alias frequency and basic resolvability flags. These are necessary sampling diagnostics, not a complete Kepler/TESS window function; gaps and aperture/shot-noise effects are excluded.

The bundled `data/reference.json` provides published-style *global* KIC 5184732 values but no full parameter covariance or individual-mode frequency table. The script therefore labels its illustrative uncertainty assumptions and simulated mode frequencies explicitly. The script prints JSON for inspection; it does not silently replace published values or add synthetic measurements to the reference dataset.
