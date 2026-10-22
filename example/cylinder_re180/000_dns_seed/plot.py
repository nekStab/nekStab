#!/usr/bin/env python3
"""plot.py: cylinder Re = 180 DNS seed.

  plot_snapshot.png - velocity magnitude of the last field file (1cyl0.f*)
  plot_signal.png   - vy(t) at the first probe of 1cyl.his

OUTPUTS: plot_snapshot.png, plot_signal.png
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
    vmin, vmax = 0.0, 1.5
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
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


def save_signal(nk, path, output):
    """Plot vy of the first probe of a history file against time."""
    with path.open() as handle:
        n_probes = int(handle.readline())
        probes = [
            tuple(map(float, handle.readline().split()[:3])) for _ in range(n_probes)
        ]
        data = np.loadtxt(handle, ndmin=2)
    if data.size == 0:
        print(f"No history rows in {path.name}: signal figure skipped")
        return
    time, signal = data[:, 0], data[:, 2]
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
    ax.plot(
        time,
        signal,
        color="C0",
        lw=0.7,
        marker="o",
        ms=1.5,
        markevery=max(1, len(time) // 60),
        label=rf"$v_y$ at ({probes[0][0]:g}, {probes[0][1]:g})",
    )
    low, high = ax.get_ylim()
    ax.set_ylim(low, high + 0.3 * (high - low))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$v_y$")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    nk = load_nekplot()
    nk.configure_style()

    snapshot = last_match("1cyl0.f*")
    if snapshot:
        save_snapshot(nk, snapshot, CASE_DIR / "plot_snapshot.png")

    history = last_match("1cyl.his")
    if history:
        save_signal(nk, history, CASE_DIR / "plot_signal.png")


if __name__ == "__main__":
    main()
