#!/usr/bin/env python3
"""Plot the Newton iteration for the forced periodic orbit of the oscillating cylinder at Re = 100.

The seed is the last field of 010_dns_forced (BF_seed_1cyl0.f00001, phase 0 of the forcing). Newton-GMRES
(userParam01 = 2.2) writes the orbit as BF_1cyl0.f00001, also at phase 0. The script reads residu_newton.dat and
both fields.

It writes newton_change.dat and plot_newton.png. newton_change.dat holds one row with:
  change_near   ||u_orbit - u_seed|| / ||u_orbit|| over the nodes with -2 < x < 10 and |y| < 4
                (u = the two velocity components)
  change_full   the same ratio over all nodes
The DNS seed has not reached the orbit: the distance is what the Newton iteration removes.

Lengths are in units of the cylinder diameter, speeds in units of the free-stream speed.

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
NEAR_X = (-2.0, 10.0)
NEAR_Y = 4.0


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def read_velocity(
    nk, path: Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the node coordinates x, y and the velocity components vx, vy of a field file."""
    x, y, fields, _time = nk.read_field(str(path))
    return x, y, fields["vx"], fields["vy"]


def main() -> None:
    """Measure the distance between the seed and the orbit and draw the residual and the change."""
    nk = load_nekplot()
    nk.configure_style()

    history = np.atleast_2d(np.loadtxt(CASE_DIR / "residu_newton.dat"))
    iteration = history[:, 0]
    vectors = history[:, 3]  # cumulative number of Arnoldi vectors
    residual = history[:, 6]
    tolerance = history[0, 7]

    x, y, vx_seed, vy_seed = read_velocity(nk, CASE_DIR / "BF_seed_1cyl0.f00001")
    _, _, vx_orbit, vy_orbit = read_velocity(nk, CASE_DIR / "BF_1cyl0.f00001")
    near = (x > NEAR_X[0]) & (x < NEAR_X[1]) & (np.abs(y) < NEAR_Y)
    delta = np.hypot(vx_orbit - vx_seed, vy_orbit - vy_seed)
    speed = np.hypot(vx_orbit, vy_orbit)
    change_near = float(np.linalg.norm(delta[near]) / np.linalg.norm(speed[near]))
    change_full = float(np.linalg.norm(delta) / np.linalg.norm(speed))
    np.savetxt(
        CASE_DIR / "newton_change.dat",
        [[change_near, change_full]],
        header="change_near  change_full   (||u_orbit - u_seed|| / ||u_orbit||)",
        fmt="%.8g",
    )
    print(
        f"{int(iteration[-1])} Newton iterations, {int(vectors[-1])} Arnoldi vectors; residual "
        f"{residual[0]:.3e} -> {residual[-1]:.3e}; distance of the seed to the orbit {change_near:.3e} "
        f"in the near field, {change_full:.3e} over all nodes"
    )

    fig = plt.figure(figsize=(7.0, 3.2), layout="constrained")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.5])
    ax_res = fig.add_subplot(grid[0, 0])
    ax_diff = fig.add_subplot(grid[0, 1])

    ax_res.semilogy(
        iteration,
        residual,
        c="C0",
        lw=1.4,
        marker="o",
        ms=3.5,
        label=f"residual ({int(vectors[-1])} Arnoldi vectors)",
    )
    ax_res.axhline(
        tolerance, c="k", lw=0.7, ls="--", label=f"tolerance {tolerance:.0e}"
    )
    ax_res.set_xlabel("Newton iteration")
    ax_res.set_ylabel(r"$\|F(q) - q\|^2$")
    ax_res.set_title("(a) Residual", loc="left", fontsize=8)
    ax_res.legend(loc="upper right", fontsize=6.5)

    bound = float(np.nanpercentile(delta[near], 99.0))
    filled = ax_diff.tricontourf(
        nk.make_triangulation(x, y),
        delta,
        levels=np.linspace(0.0, bound, 25),
        cmap="Blues",
        extend="max",
    )
    nk.add_cylinder_patches(ax_diff)
    ax_diff.set_aspect("equal")
    ax_diff.set_xlim(-2.0, 12.0)
    ax_diff.set_ylim(-3.0, 3.0)
    ax_diff.set_xlabel(r"$x$ ($D$)")
    ax_diff.set_ylabel(r"$y$ ($D$)")
    ax_diff.set_title(
        r"(b) Change of the velocity by Newton, $|\Delta u|$ at phase 0",
        loc="left",
        fontsize=8,
    )
    fig.colorbar(
        filled,
        ax=ax_diff,
        label=r"$|\Delta u|$ ($U$)",
        shrink=0.8,
        ticks=np.linspace(0.0, bound, 5),
        format="%.1e",
    )

    output = CASE_DIR / "plot_newton.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
