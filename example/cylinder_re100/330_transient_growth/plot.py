#!/usr/bin/env python3
"""Plot the optimal energy gain G(tau) of the cylinder wake at Re=40.

The script reads growth_sweep.csv, which holds one run of this stage per
horizon tau, and writes plot_envelope.png. When the field files of the
tau = 40 run are present it also writes plot_optimal_perturbation.png
(pRe1cyl0.f00001, the initial perturbation at t = 0) and
plot_optimal_response.png (ore1cyl0.f00001, its response at t = tau).

Usage: python plot.py
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(next(d for d in CASE_DIR.parents if (d / "nekplot.py").exists())))
import nekplot as nk

# Window of each field: the perturbation sits upstream of and around the cylinder,
# the response (t = tau = 40) has been carried into the wake.
WINDOWS = {"perturbation": ((-15, 10), (-6, 6)), "response": ((-2, 42), (-8, 8))}


def plot_optimal(kind, fname, output):
    """Draw vx (a) and vy (b) of an optimal perturbation or response, one colour scale each."""
    x, y, fields, _time = nk.read_field(fname)
    triang = nk.make_triangulation(x, y)
    xlim, ylim = WINDOWS[kind]
    width = nk.COL_WIDTH * (1.4 if kind == "perturbation" else 2.0)
    fig, axes = plt.subplots(2, 1, figsize=(width, width * (ylim[1] - ylim[0]) / (xlim[1] - xlim[0]) * 2.3))
    for ax, comp, tag in zip(axes, ("vx", "vy"), ("a", "b")):
        q = fields[comp]
        bound = round(float(np.nanpercentile(np.abs(q), 99)), 2)
        cf = nk.tricontourf(ax, triang, q, levels=257, cmap="RdBu", vmin=-bound, vmax=bound, extend="both")
        cbar = nk.inset_colorbar(ax, cf, orientation="horizontal", width="20%", height="6%", loc=1,
                                 ticks=[-bound, bound], tick_labels=[f"{-bound:g}", f"{bound:g}"])
        cbar.set_label(rf"$u_{comp[1]}$", fontsize=5, labelpad=1)
        nk.add_cylinder_patches(ax)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$", labelpad=-1)
        ax.set_ylabel(r"$y$", labelpad=1)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        what = "perturbation at $t = 0$" if kind == "perturbation" else r"response at $t = \tau = 40$"
        ax.set_title(rf"({tag}) optimal {what}, $u_{comp[1]}$", loc="left", fontsize=7)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def main() -> None:
    """Draw G(tau) with one marker for each computed horizon."""
    tau, gain = np.loadtxt(
        CASE_DIR / "growth_sweep.csv", delimiter=",", skiprows=3, unpack=True
    )
    fig, ax = plt.subplots(figsize=(3.5, 2.45))
    ax.plot(tau, gain, "s", c="k", ls="none", ms=4, mfc="none", label="nekStab, Re = 40")
    ax.set_xlabel(r"horizon $\tau$ ($D/U_\infty$)")
    ax.set_ylabel(r"$G(\tau)$")
    ax.legend(loc="lower right")
    output = CASE_DIR / "plot_envelope.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)

    nk.configure_style()
    for kind, name in (("perturbation", "pRe1cyl0.f00001"), ("response", "ore1cyl0.f00001")):
        if (CASE_DIR / name).exists():
            plot_optimal(kind, CASE_DIR / name, CASE_DIR / f"plot_optimal_{kind}.png")


if __name__ == "__main__":
    main()
