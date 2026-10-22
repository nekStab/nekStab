#!/usr/bin/env python3
"""plot.py: harmonically forced jet at Re = 2000, time-delayed feedback (mode 1.4).

Reads residu.dat (column 1 time, column 2 root mean square over the domain of q(t) - q(t - T), column 3 change of
the residual per step) and the last snapshot tpjet0.f0000N, the flow at the end of the run.

The residual rises while the feedback takes hold of the seed and then falls as a power law of the time. The script
fits residual = c t^p over the window t = 100 to the end of the run and writes tdf_decay.dat with one row:
p, the residual at t = 100, the residual at t = 300, the residual at the end, the largest residual of the run
and the time of that largest residual. It writes plot_residual.png with the residual against time on a logarithmic
and on a double logarithmic scale, and the axial velocity of the last snapshot below them.

Lengths are in units of the nozzle diameter D, speeds in units of the centre-line speed U, time in D / U.

USAGE: python plot.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
FIT_FROM = 100.0  # start of the power-law window
TOLERANCE = 1.0e-8  # residualTol of tpjet.par: the run writes BF_tpjet0.f00001 below it
XLIM = (0.0, 40.0)
YLIM = (0.0, 2.0)


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(
        0,
        str(next(d for d in CASE_DIR.parents if (d / "nekplot.py").exists())),
    )
    import nekplot

    return nekplot


def last_snapshot() -> Path | None:
    """Return the snapshot with the highest number, or None when the run wrote none."""
    snapshots = sorted(p for p in CASE_DIR.glob("tpjet0.f0*") if p.suffix[2:].isdigit())
    return snapshots[-1] if snapshots else None


def value_at(time: np.ndarray, values: np.ndarray, target: float) -> float:
    """Return the value of the first sample at or after the target time."""
    return float(values[np.searchsorted(time, target)])


def main() -> None:
    """Fit the decay of the residual, write tdf_decay.dat and draw the figure."""
    nk = load_nekplot()
    nk.configure_style()
    history = np.loadtxt(CASE_DIR / "residu.dat")
    time, residual = history[:, 0], history[:, 1]

    window = time >= FIT_FROM
    slope, intercept = np.polyfit(np.log(time[window]), np.log(residual[window]), 1)
    peak = int(np.argmax(residual))
    row = [
        slope,
        value_at(time, residual, 100.0),
        value_at(time, residual, 300.0),
        float(residual[-1]),
        float(residual[peak]),
        float(time[peak]),
    ]
    np.savetxt(
        CASE_DIR / "tdf_decay.dat",
        np.array([row]),
        header="p  residual(t=100)  residual(t=300)  residual(end)  residual(max)  t(max)    "
        f"(residual = c t^p fitted over t >= {FIT_FROM:g})",
        fmt="%.6e",
    )
    print(
        f"fit over t >= {FIT_FROM:g}: residual ~ t^{slope:.3f}; residual {row[1]:.3e} at t = 100, "
        f"{row[2]:.3e} at t = 300, {row[3]:.3e} at t = {time[-1]:.1f}; largest {row[4]:.3e} at t = {row[5]:.1f}"
    )

    figure = plt.figure(
        figsize=(2 * nk.COL_WIDTH, 1.55 * nk.COL_WIDTH), layout="constrained"
    )
    grid = figure.add_gridspec(2, 2, height_ratios=[1.0, 0.45])
    every = max(1, len(time) // 30)
    fitted = np.exp(intercept) * time[window] ** slope

    for index, (name, logx) in enumerate((("a", False), ("b", True))):
        axes = figure.add_subplot(grid[0, index])
        axes.plot(
            time[window],
            fitted,
            c="0.7",
            lw=3.0,
            zorder=1,
            label=rf"fit $t^{{{slope:.2f}}}$",
        )
        axes.plot(
            time,
            residual,
            c="C0",
            lw=1.2,
            marker="o",
            ms=3,
            markevery=every,
            zorder=2,
            label="residual",
        )
        axes.axhline(TOLERANCE, c="k", lw=0.8, ls="--", label=r"residualTol $10^{-8}$")
        axes.axvline(
            FIT_FROM, c="C3", lw=0.8, ls=":", label=f"fit from $t={FIT_FROM:g}$"
        )
        axes.set_yscale("log")
        if logx:
            axes.set_xscale("log")
        axes.set_xlabel(r"$t$ ($D/U$)")
        axes.set_ylabel(r"$\|q(t)-q(t-T)\|_{\mathrm{rms}}$")
        axes.set_title(
            f"({name}) Residual, " + ("logarithmic time" if logx else "linear time"),
            loc="left",
            fontsize=7.5,
        )
        axes.legend(loc="upper right", fontsize=6.0)

    snapshot = last_snapshot()
    if snapshot is not None:
        axes = figure.add_subplot(grid[1, :])
        x, y, fields, field_time = nk.read_field(str(snapshot))
        triangulation = nk.make_triangulation(x, y)
        nk.tricontourf(
            axes,
            triangulation,
            fields["vx"],
            levels=np.linspace(0.0, 1.0, 41),
            cmap="Blues",
            vmin=0.0,
            vmax=1.0,
            extend="both",
        )
        axes.set_xlim(XLIM)
        axes.set_ylim(YLIM)
        axes.set_aspect("auto")
        axes.set_xlabel(r"$z$ ($D$)")
        axes.set_ylabel(r"$r$ ($D$)")
        axes.set_title(
            rf"(c) Axial velocity at $t={float(field_time):.0f}$",
            loc="left",
            fontsize=7.5,
        )

    output = CASE_DIR / "plot_residual.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
