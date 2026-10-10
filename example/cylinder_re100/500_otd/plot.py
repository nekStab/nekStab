#!/usr/bin/env python3
"""plot.py: OTD about the steady cylinder wake at Re = 100, seeded with the leading eigenmode.

(a) eigenvalues of the reduced operator of the two OTD modes against time
    (otd_eigenvalues.dat), with the growth rate sigma of 310_stability_direct as a guide line;
(b) finite-time Lyapunov exponents of the two modes (otd_ftle.dat), with the same guide line;
(c) drift of the exponents per print step (otd_residuals.dat), the quantity the
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
    eig = np.loadtxt(CASE_DIR / "otd_eigenvalues.dat")
    ftle = np.loadtxt(CASE_DIR / "otd_ftle.dat")
    resid = np.loadtxt(CASE_DIR / "otd_residuals.dat")
    sigma = sigma_from_310()

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2), constrained_layout=True)
    styles = [(3.0, "o", "mode 1"), (1.6, "s", "mode 2")]

    # The values differ from sigma by about 3e-5, so the tick labels carry all digits.
    for ax, data, title in (
        (axes[0], eig, "(a) Eigenvalues of the reduced operator"),
        (axes[1], ftle, "(b) Exponents of the OTD modes"),
    ):
        for col, (lw, marker, label) in enumerate(styles, start=1):
            ax.plot(
                data[:, 0],
                data[:, col],
                lw=lw,
                marker=marker,
                ms=3.5,
                markevery=max(1, len(data) // 12),
                label=label,
            )
        if sigma is not None:
            ax.axhline(
                sigma,
                color="k",
                lw=0.8,
                ls="--",
                label=rf"$\sigma$ of 310 ({sigma:.7f})",
            )
        ax.ticklabel_format(axis="y", useOffset=False)
        ax.set_xlabel(r"$t$ (convective time units)")
        ax.set_ylabel(r"exponent ($1/t_c$)")
        ax.set_title(title, loc="left", fontsize=9)
        ax.legend(frameon=True, fontsize=8)

    ax = axes[2]
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
    ax.set_title("(c) Drift of the exponents", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    out = CASE_DIR / "plot_otd.png"
    fig.savefig(out, dpi=400)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
