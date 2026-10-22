#!/usr/bin/env python3
"""Plot the direct and adjoint spectra of the cubic cavity at Re=1914.

The script reads the converged eigenvalues of this stage (Spectre_NSd_conv.dat)
and of ../320_stability_adjoint (Spectre_NSa_conv.dat). Each file lists
sigma and omega in units of U/L. It writes plot_spectrum.png.

Usage: python plot.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def main() -> None:
    """Draw both spectra on one set of axes, direct thick at the back."""
    direct = np.loadtxt(CASE_DIR / "Spectre_NSd_conv.dat")
    adjoint_file = CASE_DIR.parent / "320_stability_adjoint" / "Spectre_NSa_conv.dat"
    adjoint = np.loadtxt(adjoint_file)
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.plot(direct[:, 0], direct[:, 1], "o", c="0.6", ms=9, label="direct")
    ax.plot(adjoint[:, 0], adjoint[:, 1], "s", c="k", ms=4, mfc="none", label="adjoint")
    ax.axvline(0.0, c="0.4", lw=0.6)
    ax.set_xlabel(r"growth rate $\sigma$ ($U/L$)")
    ax.set_ylabel(r"frequency $\omega$ ($U/L$)")
    ax.set_title("(a) Spectrum at Re = 1914", loc="left", fontsize=8)
    ax.set_xticks([-0.02, -0.01, 0.0])
    ax.legend(loc="center")
    output = CASE_DIR / "plot_spectrum.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
