"""Plot illustrative asymptotic echelle ridges and lognormal scaling uncertainty."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from research.seismic_diagnostics import asymptotic_modes, sample_scaling  # noqa: E402


def main() -> None:
    ref = json.loads((ROOT / "data/reference.json").read_text(encoding="utf-8"))
    p = ref["measured_parameters"]
    dnu = float(p["Dnu_microHz"])
    modes = asymptotic_modes(dnu, epsilon=1.2, d0=1.3)
    fig, ax = plt.subplots(figsize=(6.8, 4.5), layout="constrained")
    for ell in (0, 1, 2):
        selected = [m for m in modes if m["ell"] == ell]
        ax.scatter([m["echelle_phase_microhz"] for m in selected],
                   [m["frequency_microhz"] for m in selected],
                   label=f"degree l={ell}", s=22)
    ax.set(xlabel="Frequency modulo Delta nu [microHz]",
           ylabel="Asymptotic mode frequency [microHz]",
           title="Illustrative Tassoul echelle ridges (not observed modes)")
    ax.legend()
    dest = ROOT / "results/asymptotic_echelle.png"
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=180)
    plt.close(fig)

    s = np.array([0.01, 0.003, 0.015])
    corr = np.array([[1, 0.35, 0], [0.35, 1, -0.15], [0, -0.15, 1]])
    samples = sample_scaling((float(p["numax_microHz"]), dnu, float(p["Teff_K"])),
                             np.outer(s, s) * corr, samples=10000)
    fig, ax = plt.subplots(figsize=(6.8, 4.5), layout="constrained")
    ax.hexbin(samples[:, 0], samples[:, 1], gridsize=45, mincnt=1)
    ax.set(xlabel="Scaling mass [solar units]", ylabel="Scaling radius [solar units]",
           title="Illustrative correlated uncertainties (not LEGACY errors)")
    dest = ROOT / "results/mass_radius_uncertainty.png"
    fig.savefig(dest, dpi=180)
    plt.close(fig)
    print("Created two labelled illustrative diagnostics")


if __name__ == "__main__":
    main()
