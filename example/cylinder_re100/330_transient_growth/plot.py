#!/usr/bin/env python3
"""Plot the optimal energy gain G(tau) of the cylinder wake at Re=40.

The script reads growth_sweep.csv, which holds one run of this stage per
horizon tau, and writes plot_envelope.png.

Usage: python plot.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def main() -> None:
    """Draw G(tau) with one marker for each computed horizon."""
    tau, gain = np.loadtxt(
        CASE_DIR / "growth_sweep.csv", delimiter=",", skiprows=3, unpack=True
    )
    fig, ax = plt.subplots(figsize=(3.5, 2.45))
    ax.plot(tau, gain, "-s", c="k", lw=0.8, ms=4, mfc="none", label="nekStab, Re = 40")
    ax.set_xlabel(r"horizon $\tau$ ($D/U_\infty$)")
    ax.set_ylabel(r"$G(\tau)$")
    ax.legend(loc="lower right")
    output = CASE_DIR / "plot_envelope.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
