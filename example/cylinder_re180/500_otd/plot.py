#!/usr/bin/env python3
"""plot.py: OTD about the Newton orbit of the cylinder at Re = 180 (userParam01 = 5.1, lpert = 1).

(a) exponent of the OTD mode for each period of the orbit (otd_ftle_period.dat), with the
    exponent of the phase mode of the orbit, 0, as a guide line;
(b) the absolute value of the same exponents on a logarithmic axis, with the tolerance of
    ref/reference.json as a guide line.

OUTPUT: plot_otd.png
USAGE:  python plot.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
TOLERANCE = 1.0e-4  # tol_abs of otd_period_exponent_final in ref/reference.json


def main() -> None:
    data = np.loadtxt(CASE_DIR / "otd_ftle_period.dat")
    t, sigma = data[:, 0], data[:, 1]
    label = f"OTD mode, {len(t)} periods"

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)

    ax = axes[0]
    ax.plot(t, sigma, lw=1.6, marker="o", ms=3.5, label=label)
    ax.axhline(0.0, color="k", lw=0.8, ls="--", label="phase mode of the orbit (0)")
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"exponent per period ($1/t_c$)")
    ax.set_title("(a) Exponent of the OTD mode", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    ax = axes[1]
    ax.semilogy(t, np.abs(sigma), lw=1.6, marker="o", ms=3.5, label=label)
    ax.axhline(
        TOLERANCE,
        color="k",
        lw=0.8,
        ls="--",
        label=f"tolerance of the reference ({TOLERANCE:.0e})",
    )
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"$|$exponent per period$|$ ($1/t_c$)")
    ax.set_title("(b) Absolute value", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    out = CASE_DIR / "plot_otd.png"
    fig.savefig(out, dpi=400)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
