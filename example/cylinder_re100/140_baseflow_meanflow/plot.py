#!/usr/bin/env python3
"""Plot the mean-flow stage of the Re = 100 cylinder wake.

The script reads the probe history 1cyl.his, the averaging window meanflow_window.dat and the
mean flow BFA1cyl0.f00001. When the Newton solution of ../210_baseflow_newton/fp is present, it
draws the steady base flow of the same wake next to the mean flow.

It writes meanflow_axis.dat (recirculation length and minimum velocity of the mean flow on the
axis y = 0) and plot_meanflow.png. Lengths are in units of the cylinder diameter, velocities in
units of the free-stream speed, times in units of D/U.

Usage: python plot.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
STEADY_FIELD = CASE_DIR.parent / "210_baseflow_newton" / "fp" / "BF_1cyl0.f00001"
CYLINDER_RADIUS = 0.5
AXIS_TOLERANCE = 1e-5


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def read_probe(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return time and v of the probe in a Nek5000 history file.

    The first line holds the number of probes, followed by one line per probe point.
    Each later line holds the time and then u, v, p at the probe.
    """
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    n_probes = int(lines[0][0])
    rows = np.array(
        [[float(v) for v in line] for line in lines[1 + n_probes :] if len(line) >= 3]
    )
    return rows[:, 0], rows[:, 2]


def axis_profile(
    x: np.ndarray, y: np.ndarray, u: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return x and u on the axis y = 0 behind the cylinder, sorted by x.

    The axis is a line of element faces, so every GLL node of the axis is a mesh node. The field
    file stores the coordinates in single precision, so y of these nodes is zero only to about 1e-7.
    A node on a face belongs to two elements; it appears once in the result.
    """
    on_axis = (np.abs(y) < AXIS_TOLERANCE) & (x >= CYLINDER_RADIUS - AXIS_TOLERANCE)
    x_nodes, first = np.unique(x[on_axis], return_index=True)
    return x_nodes, u[on_axis][first]


def recirculation_length(x_axis: np.ndarray, u_axis: np.ndarray) -> float:
    """Return the distance from the rear of the cylinder to the end of the reversed flow.

    The end is the last point where u on the axis goes from negative to positive,
    found by linear interpolation between two mesh nodes.
    """
    reversed_flow = np.flatnonzero(u_axis < 0.0)
    if reversed_flow.size == 0 or reversed_flow[-1] + 1 >= u_axis.size:
        raise ValueError("no end of the recirculation bubble on the axis")
    i = reversed_flow[-1]
    x_end = x_axis[i] - u_axis[i] * (x_axis[i + 1] - x_axis[i]) / (
        u_axis[i + 1] - u_axis[i]
    )
    return float(x_end - CYLINDER_RADIUS)


def main() -> None:
    """Write meanflow_axis.dat and draw the averaging window, the axis profile and the mean flow."""
    nk = load_nekplot()
    nk.configure_style()

    mean_field = CASE_DIR / "BFA1cyl0.f00001"
    x, y, fields, _time = nk.read_field(str(mean_field))
    x_axis, u_axis = axis_profile(x, y, fields["vx"])
    length_mean = recirculation_length(x_axis, u_axis)
    np.savetxt(
        CASE_DIR / "meanflow_axis.dat",
        [[length_mean, u_axis.min(), x_axis[np.argmin(u_axis)]]],
        header="recirculation_length  u_min_on_axis  x_of_u_min",
        fmt="%.8g",
    )
    print(
        f"mean flow: recirculation length {length_mean:.4f} behind the rear, u_min on axis {u_axis.min():.4f}"
    )

    window = np.loadtxt(CASE_DIR / "meanflow_window.dat", ndmin=2)[0]
    t_open, t_close, n_periods, period, strouhal = window
    print(
        f"window t = {t_open:.3f} to {t_close:.3f}: {int(n_periods)} periods of {period:.5f} (St = {strouhal:.6f})"
    )

    steady = None
    if STEADY_FIELD.is_file():
        xs, ys, steady_fields, _ = nk.read_field(str(STEADY_FIELD))
        steady = (xs, ys, steady_fields["vx"])
        x_steady, u_steady = axis_profile(xs, ys, steady_fields["vx"])
        length_steady = recirculation_length(x_steady, u_steady)
        print(f"steady base flow: recirculation length {length_steady:.4f}")

    fig = plt.figure(figsize=(7.0, 5.0), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0], width_ratios=[1.0, 1.5])
    ax_probe = fig.add_subplot(grid[0, :])
    ax_axis = fig.add_subplot(grid[1, 0])
    ax_field = fig.add_subplot(grid[1, 1])

    t_probe, v_probe = read_probe(CASE_DIR / "1cyl.his")
    ax_probe.plot(
        t_probe,
        v_probe,
        c="0.45",
        lw=0.8,
        marker="o",
        ms=1.5,
        markevery=max(1, t_probe.size // 400),
        label="$v$ at $(x, y) = (5, 0)$",
    )
    ax_probe.axvline(
        t_open, c="C0", lw=1.0, ls="--", label=f"average opens, $t = {t_open:.2f}$"
    )
    ax_probe.axvline(
        t_close, c="C3", lw=1.0, ls="--", label=f"average closes, $t = {t_close:.2f}$"
    )
    ax_probe.set_xlabel(r"$t$ ($D/U$)")
    ax_probe.set_ylabel(r"$v$ ($U$)")
    ax_probe.set_title(
        f"(a) Probe signal and the average over {int(n_periods)} periods",
        loc="left",
        fontsize=8,
    )
    ax_probe.set_ylim(1.5 * v_probe.min(), 1.5 * v_probe.max())
    ax_probe.legend(loc="upper center", fontsize=6.5, ncol=3)

    if steady is not None:
        ax_axis.plot(
            x_steady,
            u_steady,
            c="0.7",
            lw=2.5,
            marker="s",
            ms=2.5,
            markevery=8,
            label=f"steady base flow, $L$ = {length_steady:.2f}",
        )
    ax_axis.plot(
        x_axis,
        u_axis,
        c="C0",
        lw=1.2,
        marker="o",
        ms=2.5,
        markevery=8,
        label=f"mean flow, $L$ = {length_mean:.2f}",
    )
    ax_axis.axhline(0.0, c="k", lw=0.6, ls=":", label="$u = 0$")
    ax_axis.set_xlim(CYLINDER_RADIUS, 10.0)
    ax_axis.set_xlabel(r"$x$ ($D$)")
    ax_axis.set_ylabel(r"$u(x, 0)$ ($U$)")
    ax_axis.set_title("(b) Velocity on the axis", loc="left", fontsize=8)
    ax_axis.set_ylim(-0.25, 1.15)
    ax_axis.legend(loc="upper left", fontsize=6.5)

    triang = nk.make_triangulation(x, y)
    speed_levels = np.linspace(-0.2, 1.2, 29)
    filled = ax_field.tricontourf(
        triang, fields["vx"], levels=speed_levels, cmap="RdBu_r", extend="both"
    )
    mean_line = ax_field.tricontour(
        triang, fields["vx"], levels=[0.0], colors="k", linewidths=1.2
    )
    handles = [mean_line.legend_elements()[0][0]]
    labels = ["$u = 0$, mean flow"]
    if steady is not None:
        steady_line = ax_field.tricontour(
            nk.make_triangulation(steady[0], steady[1]),
            steady[2],
            levels=[0.0],
            colors="k",
            linewidths=1.2,
            linestyles="--",
        )
        handles.append(steady_line.legend_elements()[0][0])
        labels.append("$u = 0$, steady base flow")
    nk.add_cylinder_patches(ax_field)
    ax_field.set_aspect("equal")
    ax_field.set_xlim(-1.5, 10.0)
    ax_field.set_ylim(-3.0, 3.0)
    ax_field.set_xlabel(r"$x$ ($D$)")
    ax_field.set_ylabel(r"$y$ ($D$)")
    ax_field.set_title(
        "(c) Streamwise velocity of the mean flow", loc="left", fontsize=8
    )
    ax_field.legend(handles, labels, loc="upper right", fontsize=6.5)
    fig.colorbar(
        filled, ax=ax_field, label=r"$u$ ($U$)", ticks=[-0.2, 0.0, 0.5, 1.0], shrink=0.8
    )

    output = CASE_DIR / "plot_meanflow.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
