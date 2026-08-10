# Asteroseismology Mode Lab

Solar-like p-mode combs, large separations and mode-envelope diagnostics.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `scripts/validate.js`: no-dependency repository validation.

## Run

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and absence of unfinished scaffold tokens.

## Reference Data

Low-degree solar-like oscillation anchors compressed for validating the mode-comb renderer.

## References

- Aerts, C., Christensen-Dalsgaard, J. and Kurtz, D.W., 2010. Asteroseismology. Springer.
- Chaplin, W.J. and Miglio, A., 2013. Asteroseismology of solar-type and red-giant stars. Annual Review of Astronomy and Astrophysics, 51, pp.353-392.
