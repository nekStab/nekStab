#!/usr/bin/env python3
"""plot.py: cylinder Re = 100 DMD panels.

  plot_spectrum.png - growth rate against Strouhal number of the DMD
                      eigenvalues, with the multiples of the shedding frequency
  plot_svd.png      - singular values of the snapshot matrix and the rank used
  plot_mode.png     - vy of the real part of the shedding mode

The saved modes are the first modal_nsave rows of dmd_spectrum.dat
(dm11cyl0.f00001, dm11cyl0.f00002, ...). The shedding mode is the saved mode
with the lowest positive Strouhal number.

OUTPUTS: plot_spectrum.png, plot_svd.png, plot_mode.png
USAGE:   python plot.py
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
MODE_GLOB = "dm11cyl0.f*"


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    root = next(d for d in CASE_DIR.parents if (d / "nekplot.py").exists())
    sys.path.insert(0, str(root))
    import nekplot

    return nekplot


def read_spectrum(path):
    """Return the table of dmd_spectrum.dat.

    The columns are mode, |mu|, sigma, omega, St, Re mu and Im mu.
    """
    rows = [
        [float(value) for value in line.split()[:7]]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    return np.array(rows)


def read_svd(path):
    """Return the singular values and the rank used from dmd_svd.dat."""
    text = path.read_text()
    rank = int(re.search(r"rank_used\s*=\s*(\d+)", text).group(1))
    rows = [
        [float(value) for value in line.split()]
        for line in text.splitlines()
        if line.strip() and not line.startswith("#")
    ]
    return np.array(rows)[:, 1], rank


def shedding_index(table, n_saved):
    """Return the row (0-based) of the saved mode with the lowest positive St."""
    saved = table[:n_saved]
    positive = np.flatnonzero(saved[:, 4] > 0)
    return int(positive[np.argmin(saved[positive, 4])])


def save_spectrum(nk, table, fundamental, output):
    """Draw growth rate against Strouhal number for the eigenvalues with St >= 0."""
    keep = table[:, 4] >= 0
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.7))
    for multiple in range(1, 6):
        ax.axvline(
            multiple * fundamental,
            color="0.6",
            lw=0.6,
            ls=":",
            label=rf"multiples of $St$ = {fundamental:.4f}" if multiple == 1 else None,
        )
    ax.axhline(0, color="k", lw=0.5, ls="--", label=r"$\sigma = 0$")
    ax.plot(
        table[keep, 4],
        table[keep, 2],
        ls="none",
        marker="o",
        ms=4,
        color="C0",
        label="DMD eigenvalue",
    )
    ax.set_xlabel(r"$St$")
    ax.set_ylabel(r"$\sigma$")
    ax.legend(fontsize=7, loc="lower right")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_svd(nk, path, output):
    """Draw the singular values of the snapshot matrix and mark the rank used."""
    sigma, rank = read_svd(path)
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.7))
    ax.axvline(rank, color="0.6", lw=0.6, ls=":", label=f"rank used = {rank}")
    ax.plot(
        np.arange(1, len(sigma) + 1),
        sigma,
        color="C0",
        lw=0.7,
        marker="o",
        ms=2.5,
        label="singular value",
    )
    ax.set_yscale("log")
    ax.set_xlabel(r"$i$")
    ax.set_ylabel(r"$\sigma_i$")
    ax.legend(fontsize=7, loc="upper right")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_mode(nk, path, output):
    """Draw vy of a mode file with a colour range that is symmetric about zero."""
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    values = fields["vy"]
    bound = round(float(np.nanpercentile(np.abs(values), 99)), 3) or 1.0
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
    cf = nk.tricontourf(
        ax,
        triang,
        values,
        levels=257,
        cmap="RdBu_r",
        vmin=-bound,
        vmax=bound,
        extend="both",
    )
    nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=[-bound, bound],
        tick_labels=[f"{-bound:g}", f"{bound:g}"],
    )
    nk.add_cylinder_patches(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_xlabel(r"$x$", labelpad=-1)
    ax.set_ylabel(r"$y$", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    nk = load_nekplot()
    nk.configure_style()

    spectrum_file = CASE_DIR / "dmd_spectrum.dat"
    modes = sorted(CASE_DIR.glob(MODE_GLOB))
    if spectrum_file.exists():
        table = read_spectrum(spectrum_file)
        n_saved = len(modes) or 1
        fundamental = table[shedding_index(table, n_saved), 4]
        save_spectrum(nk, table, fundamental, CASE_DIR / "plot_spectrum.png")
        if modes:
            index = shedding_index(table, len(modes))
            save_mode(nk, modes[index], CASE_DIR / "plot_mode.png")

    svd_file = CASE_DIR / "dmd_svd.dat"
    if svd_file.exists():
        save_svd(nk, svd_file, CASE_DIR / "plot_svd.png")


if __name__ == "__main__":
    main()
