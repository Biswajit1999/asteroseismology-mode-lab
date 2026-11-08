# Seismic uncertainty and window diagnostics (8 November 2026)

This standalone standard-library Python module is an **independent numerical audit**, not a replacement for the browser model and not a fit to Kepler/TESS observations.

1. The uncorrected Kjeldsen–Bedding solar scalings give mass and radius from nu_max, Delta_nu and effective temperature.
2. The 3x3 input covariance propagates to the joint 2x2 log-mass/log-radius covariance via the Jacobian. Correlated measurement errors can either increase or reduce derived uncertainties.
3. The normalized complex Fourier sampling window exposes aliases; hourly sampling aliases at 1/hour. It is not a power-spectrum detection significance test.

Run: `cd research/seismic_robustness && python -m unittest -v test_robustness.py`.

Assumptions: first-order linear propagation, uncorrected scaling laws, finite observation times, no granulation, no stochastic mode excitation or surface corrections. Validation is on analytic/synthetic cases only. AI-assisted implementation; Biswajit Jana is the requested sole Git commit author.
