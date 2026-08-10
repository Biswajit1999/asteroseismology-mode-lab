# Asteroseismology Mode Lab

Solar-like p-mode combs, large separations and mode-envelope diagnostics.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread. It is a teaching/exploration tool for the shape of a solar-like oscillation power spectrum, not a fitting pipeline for real time-series data.

## Background: stellar oscillations as structural probes

Stars like the Sun are not acoustically quiet. Convective turbulence in the outer envelope stochastically excites a spectrum of resonant standing waves that ring throughout the stellar interior. Observing these oscillations — asteroseismology — gives an otherwise inaccessible view of stellar interior structure, in the same way helioseismology mapped the Sun's rotation profile and convection zone.

**p-modes vs g-modes.** Two broad oscillation families exist:

- **p-modes (pressure/acoustic modes):** the restoring force is the pressure gradient. They are essentially trapped sound waves, most sensitive to the outer envelope and well described (for main-sequence and subgiant stars) by an *asymptotic* frequency pattern.
- **g-modes (gravity modes):** the restoring force is buoyancy. They propagate in stably stratified, non-convective regions (typically the deep radiative interior) and are the primary probe of core structure, but they are evanescent in convective envelopes and hard to detect from the surface in Sun-like stars.

This lab models the **p-mode** regime, which is what dominates the observed power spectra of the Sun and other solar-like oscillators.

**The asymptotic relation.** For high-order, low-degree p-modes, Tassoul's asymptotic theory predicts frequencies organized into a near-regular comb:

```
nu(n, l) ≈ Δν (n + l/2 + ε) − l(l+1) D0
```

where `n` is the radial order, `l` the angular degree, `Δν` the large frequency separation, `ε` a surface-layer phase offset, and `D0` a small curvature term related to the sound-speed gradient near the core. Modes of a given `l` form a near-equally-spaced comb in frequency with spacing `Δν`; modes of adjacent even/odd degree (l=0/2 and l=1/3) interleave, and the *small separation* `δν` between them is sensitive to the core's evolutionary state.

**Large separation Δν and nu_max.** Two headline observables summarize the whole comb:

- **Δν**, the large frequency separation, is set by the sound-crossing time of the star and scales with mean stellar density: `Δν ∝ sqrt(M / R^3)`.
- **nu_max**, the frequency of maximum oscillation power, tracks the acoustic cutoff frequency of the atmosphere, which depends on surface gravity and effective temperature.

**Scaling relations (Kjeldsen & Bedding 1995).** Combining the two observables with an effective temperature gives near-direct estimates of stellar mass and radius relative to solar values:

```
nu_max / nu_max,sun ≈ (g / g_sun) * (Teff / Teff,sun)^(-1/2)
Δν / Δν,sun         ≈ sqrt(M/Msun) * (R/Rsun)^(-3/2)
```

Solving the pair for mass and radius (with `Teff` from spectroscopy or photometry) yields:

```
M/Msun ≈ (nu_max/nu_max,sun)^3 (Δν/Δν,sun)^-4 (Teff/Teff,sun)^1.5
R/Rsun ≈ (nu_max/nu_max,sun)   (Δν/Δν,sun)^-2 (Teff/Teff,sun)^0.5
```

Because Δν tracks mean density and nu_max tracks surface gravity, and both correlate with stellar age through the star's evolving radius and luminosity, these two numbers — extracted from nothing more than a power spectrum — give asteroseismic mass, radius, and (via stellar-evolution-model grids) age, all to a few-percent precision for bright solar-like targets.

**Kepler and TESS context.** The *Kepler* mission's four-year, ultra-precise photometry made asteroseismology routine for tens of thousands of solar-type and red-giant stars, anchoring exoplanet host-star parameters (and hence planet radii) and testing stellar-evolution physics across the HR diagram. *TESS*, with shorter but nearly all-sky coverage, extends this to the brightest, closest solar-like oscillators, prioritizing targets bright enough for granulation and p-mode signal to rise above photon noise in a ~27-day (or longer, in the continuous viewing zones) baseline.

## How it works

This is a zero-build, dependency-free static site: `index.html` + vanilla JS, no bundler.

1. **`index.html` / `styles.css`** lay out a dark "mission control" dashboard: a controls panel, a spectrum plot canvas, a parameter-response heatmap canvas, and a telemetry/notes panel.
2. **`app.js`** owns UI state. On load it fetches `data/reference.json` — the published low-degree solar oscillation anchor points — and renders them as fixed reference markers on the plot. It builds slider controls from the `LAB.controls` definition (`deltaNu`, `nuMax`, `linewidth`, `inclination`), and on every slider change posts the current parameter set to a `Worker`.
3. **`physicsWorker.js`** runs off the UI thread. Its `asteroseismology(p)` function builds a synthetic p-mode power spectrum over 1000–5000 microHz: a broad Gaussian mode envelope centered on `nuMax` (width fixed at 650 microHz, standing in for the observed power envelope), a comb of Lorentzian-like Gaussian peaks spaced by `deltaNu/2` (approximating the interleaved even/odd-degree ridges of the true `Δν`-spaced comb), each peak's amplitude modulated by its position in the envelope and by a cosine term parameterized by `inclination` (a simplified stand-in for how rotational splitting/visibility varies with the star's inclination to the line of sight), and peak sharpness controlled by `linewidth`. It also produces a companion "echelle response" heatmap — a 2D field whose oscillatory ridge mimics the way an échelle diagram (power folded modulo Δν) shows near-vertical mode ridges.
4. The worker posts back `{series, metrics, heatmap}`. `app.js` draws the model spectrum plus the fixed reference anchors on a Canvas 2D plot with auto-scaled axes, renders the heatmap to a second canvas via `ImageData`, and lists live metrics in the telemetry panel: `nu_max`, `large_separation`, `peak_power`, and two scaling-relation proxies, `mass_proxy_Msun` and `radius_proxy_Rsun`, computed from the current `nu_max`/`Δν` sliders relative to solar reference values (`nu_max,sun` = 3090 microHz, `Δν,sun` = 135.1 microHz) via the Kjeldsen & Bedding (1995) relations, assuming solar `Teff` since the lab has no effective-temperature control.
5. **`data/reference.json`** is a small, auditable set of published-style anchor points (envelope edges, p-mode ridge, nu_max) with a citation, used purely as a visual reference overlay — not as fit targets, since the model is illustrative rather than a least-squares fit to data.
6. **`research-overlay.js`** injects a small non-invasive status panel (bottom-right) confirming the research-quality validation layer is active; it is cosmetic and does not affect the simulation.
7. **`scripts/validate.js`** (`npm run check`) is a no-dependency sanity check: confirms required files exist, the reference JSON has finite points and required citations appear in the README, worker/app syntax is valid (`node --check`), and no unfinished scaffold markers remain in the tracked source. `scripts/validate_repository.mjs` (`npm run validate:research`) checks the smaller `data/research-reference.json` metadata bundle in the same spirit.

Everything the model computes is a hand-tuned illustrative approximation of the true asymptotic p-mode pattern — good for building intuition about how Δν, nu_max, linewidth and inclination reshape a solar-like spectrum, not a substitute for a real peak-bagging pipeline (e.g. fitting Lorentzians to *Kepler*/TESS power spectra).

## Usage

Zero build step — just serve the static files (fetch of `data/reference.json` requires an HTTP origin, so opening `index.html` directly via `file://` will not work in most browsers).

```bash
python -m http.server 8080
```

Open `http://localhost:8080`. Drag the `Large separation [microHz]`, `nu_max [microHz]`, `Line width [microHz]`, and `Inclination [deg]` sliders and watch the spectrum, échelle-style heatmap, and telemetry panel update live via the worker.

## Validate

```bash
npm run check
npm run validate:research
```

`npm run check` verifies required files, reference-data integrity, worker/app syntax, README citation coverage, and the absence of unfinished scaffold markers. `npm run validate:research` checks the smaller research-quality metadata bundle in `data/research-reference.json`.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation (`asteroseismology` lab function, plus other unused lab functions retained from the shared worker template).
- `data/reference.json`: small auditable reference-data bundle (published-style p-mode anchor points).
- `data/research-reference.json`: small metadata bundle used by the repository-level validation script.
- `research-overlay.js`: cosmetic research-quality status panel.
- `scripts/validate.js`: no-dependency repository validation.
- `scripts/validate_repository.mjs`: no-dependency research-metadata validation.

## Reference Data

Low-degree solar-like oscillation anchors compressed for validating the mode-comb renderer.

## Math appendix

Symbols: `n` radial order, `l` angular degree, `Δν` large separation, `δν` small separation, `ε` surface phase term, `D0` core-sensitive curvature term, `nu_max` frequency of maximum oscillation power, `g` surface gravity, `Teff` effective temperature.

**Asymptotic p-mode relation:**

```
nu(n, l) ≈ Δν (n + l/2 + ε) − l(l+1) D0
```

**Large separation from mean density:**

```
Δν ∝ sqrt(M / R^3) = sqrt(rho_mean)
```

**Frequency of maximum power (Kjeldsen & Bedding 1995):**

```
nu_max ∝ g * Teff^(-1/2)
```

**Mass and radius from the scaling relations:**

```
M/Msun ≈ (nu_max/nu_max,sun)^3 (Δν/Δν,sun)^-4 (Teff/Teff,sun)^1.5
R/Rsun ≈ (nu_max/nu_max,sun)   (Δν/Δν,sun)^-2 (Teff/Teff,sun)^0.5
```

**This repository's illustrative spectrum model** (`physicsWorker.js`, `asteroseismology(p)`), for parameters `deltaNu`, `nuMax`, `linewidth`, `inclination`:

```
envelope(f)   = exp( -0.5 * ((f − nuMax) / 650)^2 )
mode centers  = nuMax + n * (deltaNu / 2),   n = −12 … 14
spectrum(f)   = 0.02 + Σ_n envelope(f0_n) * exp( -0.5 * ((f − f0_n) / linewidth)^2 )
                          * (1 + 0.35 * cos(inclination_rad + n))
```

This is a hand-built approximation of the asymptotic comb (fixed 650 microHz envelope width, `deltaNu/2` spacing standing in for interleaved degrees, and a cosine visibility term standing in for inclination/rotational-splitting effects) — chosen for interactive responsiveness and visual clarity, not derived from a stellar-structure model.

## References

- Aerts, C., Christensen-Dalsgaard, J. and Kurtz, D.W., 2010. Asteroseismology. Springer.
- Chaplin, W.J. and Miglio, A., 2013. Asteroseismology of solar-type and red-giant stars. Annual Review of Astronomy and Astrophysics, 51, pp.353-392.
- Kjeldsen, H. and Bedding, T.R., 1995. Amplitudes of stellar oscillations: the implications for asteroseismology. Astronomy and Astrophysics, 293, pp.87-106.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors, equations and research boundaries added to this repository.
