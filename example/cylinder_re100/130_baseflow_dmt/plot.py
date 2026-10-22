#!/usr/bin/env python3
"""Plot the DMT base-flow stage of the Re = 100 cylinder wake.

Writes plot_ic.png (the seed), plot_residual.png (the residual of the run
against the convergence gate) and plot_bf.png (the converged base flow).

Usage: python plot.py
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def par_value(name):
    """Return the text after `name =` in the .par file, without the comment."""
    text = next(CASE_DIR.glob("*.par")).read_text()
    match = re.search(rf"^\s*{name}\s*=\s*([^\s#]+)", text, re.M | re.I)
    return match.group(1) if match else None


def render_field(nk, ax, path, vmax=1.5):
    """Draw the velocity magnitude of a field file around the cylinder."""
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    speed = np.hypot(fields["vx"], fields["vy"])
    cf = nk.tricontourf(
        ax, triang, speed, levels=257, cmap="Blues", vmin=0, vmax=vmax, extend="max"
    )
    nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=[0, vmax],
        tick_labels=["0", f"{vmax}"],
    )
    ax.set_aspect("equal")
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_xlabel(r"$x$", labelpad=-1)
    ax.set_ylabel(r"$y$", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def save_field(nk, path, output):
    """Render one field file to a wide figure."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
    render_field(nk, ax, path)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_residual(nk, path, gate, output):
    """Plot the residual history of the run with the gate as a horizontal line."""
    data = np.genfromtxt(path)
    time, residual = data[:, 0], data[:, 1]
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
    ax.semilogy(
        time,
        residual,
        color="C0",
        lw=0.8,
        marker="o",
        ms=2,
        markevery=max(1, len(time) // 40),
        label=f"DMT ({time[-1]:.1f} t.u.)",
    )
    ax.axhline(gate, color="k", ls=":", lw=0.5, label=f"gate {gate:.0e}")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\|u - \bar{u}\| / \Delta t$")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, which="both", ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    nk = load_nekplot()
    nk.configure_style()

    seed = CASE_DIR / par_value("startFrom")
    if seed.exists():
        save_field(nk, seed, CASE_DIR / "plot_ic.png")

    residual_file = CASE_DIR / "residu_dmt.dat"
    if residual_file.exists():
        save_residual(
            nk,
            residual_file,
            float(par_value("userParam15")),
            CASE_DIR / "plot_residual.png",
        )

    base_flow = CASE_DIR / "BF_1cyl0.f00001"
    if base_flow.exists():
        save_field(nk, base_flow, CASE_DIR / "plot_bf.png")


if __name__ == "__main__":
    main()
