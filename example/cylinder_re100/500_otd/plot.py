#!/usr/bin/env python3
"""plot.py: OTD about the steady cylinder wake at Re = 100, seeded with the leading eigenmode.

(a) finite-time Lyapunov exponents of the two OTD modes against time, with the
    growth rate sigma of 310_stability_direct as a guide line;
(b) drift of the exponents per time step (otd_residuals.dat), the quantity the
    convergence test of the run reads.

OUTPUT: plot_otd.png
USAGE:  python plot.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
SPECTRUM = CASE_DIR.parent / "310_stability_direct" / "direct" / "Spectre_NSd.dat"


def sigma_from_310() -> float | None:
    """Growth rate of the leading eigenvalue of 310_stability_direct, or None without the file."""
    if not SPECTRUM.exists():
        return None
    return float(SPECTRUM.read_text().split()[0].replace("D", "E"))


def main() -> None:
    ftle = np.loadtxt(CASE_DIR / "otd_ftle.dat")
    resid = np.loadtxt(CASE_DIR / "otd_residuals.dat")
    sigma = sigma_from_310()

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)

    ax = axes[0]
    styles = [(3.0, "o", "FTLE 1"), (1.6, "s", "FTLE 2")]
    for col, (lw, marker, label) in enumerate(styles, start=1):
        ax.plot(
            ftle[:, 0],
            ftle[:, col],
            lw=lw,
            marker=marker,
            ms=3.5,
            markevery=max(1, len(ftle) // 12),
            label=label,
        )
    if sigma is not None:
        ax.axhline(
            sigma, color="k", lw=0.8, ls="--", label=rf"$\sigma$ of 310 ({sigma:.4f})"
        )
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"finite-time exponent ($1/t_c$)")
    ax.set_title("(a) Exponents of the OTD modes", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    ax = axes[1]
    for col, (lw, marker, label) in enumerate(styles, start=1):
        ax.semilogy(
            resid[:, 0],
            resid[:, col],
            lw=lw,
            marker=marker,
            ms=3.5,
            markevery=max(1, len(resid) // 12),
            label=label,
        )
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"$|\Delta\,\mathrm{FTLE}|$ per print step ($1/t_c$)")
    ax.set_title("(b) Drift of the exponents", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    out = CASE_DIR / "plot_otd.png"
    fig.savefig(out, dpi=400)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
