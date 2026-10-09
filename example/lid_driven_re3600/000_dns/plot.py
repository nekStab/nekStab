#!/usr/bin/env python3
"""plot.py: lid-driven cavity Re = 3600 DNS.

  plot_snapshot.png - velocity magnitude of the last field file (cav0.f*)

OUTPUTS: plot_snapshot.png
USAGE:   python plot.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    root = next(d for d in CASE_DIR.parents if (d / "nekplot.py").exists())
    sys.path.insert(0, str(root))
    import nekplot

    return nekplot


def last_match(pattern):
    """Return the match of the glob with the highest field number, or None."""
    matches = sorted(CASE_DIR.glob(pattern))
    return matches[-1] if matches else None


def save_snapshot(nk, path, output):
    """Draw velocity magnitude of a field file."""
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    values = np.hypot(fields["vx"], fields["vy"])
    vmin, vmax = 0.0, round(float(np.nanpercentile(values, 99)), 3)
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.8, nk.COL_WIDTH * 0.8))
    cf = nk.tricontourf(
        ax,
        triang,
        values,
        levels=257,
        cmap="Blues",
        vmin=vmin,
        vmax=vmax,
        extend="max",
    )
    nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=[vmin, vmax],
        tick_labels=[f"{vmin:g}", f"{vmax:g}"],
    )
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_aspect("equal")
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

    snapshot = last_match("cav0.f*")
    if snapshot:
        save_snapshot(nk, snapshot, CASE_DIR / "plot_snapshot.png")


if __name__ == "__main__":
    main()
