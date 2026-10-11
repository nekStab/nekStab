#!/usr/bin/env python3
"""Plot the DMT filter bank on the forced periodic orbit of the oscillating cylinder at Re = 100.

The run starts from the orbit of 210_baseflow_newton (BF_1cyl0.f00001, phase 0 of the forcing, time 0) and
integrates 30 periods T = 1 / f of the forcing with four filters: k = 1, 2, 3 at the frequency of the forcing
and its second and third harmonic, and k = 4 at frequency 0 (the mean flow). The files bRe1cyl0.f0000k and
bIm1cyl0.f0000k hold Re and Im of the amplitude uhat_k of each filter, averaged over the last two periods, with
u(x, t) = Re[uhat_k(x) exp(i k omega t)] and t the time of the run (phase 0 at t = 0).

The orbit at phase 0 is the sum of the mean and of Re uhat_k over all harmonics. The script rebuilds it
from the filters, one harmonic at a time, and compares with the orbit of stage 210:

  e_m = || u_orbit - (mean + sum_{k <= m} Re uhat_k) || / || u_orbit - mean ||,   m = 0, 1, 2, 3

over the nodes with -2 < x < 10 and |y| < 4 (u = the two velocity components). It writes bank_check.dat with
one row per m: m, e_m (near field), e_m (all nodes), and || uhat_m || / || mean || (near field; for m = 0 the
norm of the mean over itself, 1). It also writes plot_dmt_bank.png.

On the wall of the cylinder the flow has the wall velocity, A 2 pi f cos(omega t) along y, so uhat_1 = (0, A omega)
and the other amplitudes are 0 there. bank_wall.dat holds, for each filter, the mean of Re uhat_y and Im uhat_y over
the nodes on the wall.

Lengths are in units of the cylinder diameter D, speeds in units of the free-stream speed U, time in D / U.

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
FORCING_FREQUENCY = 0.1643  # userParam06 of 1cyl.par
NEAR_X = (-2.0, 10.0)
NEAR_Y = 4.0
HARMONICS = (1, 2, 3)  # filters 1 to 3; filter 4 is the mean flow
MEAN = 4
AMPLITUDE = 1.08  # userParam05 of 1cyl.par, in D
WALL_TOLERANCE = 1.0e-4  # in D: distance of a wall node from the circle of radius 0.5 D at phase 0
# The run ends one step (2.4e-4 of a period) after 30 periods, and the cylinder moves 1.6e-3 D in it.
PHASE_TOLERANCE = 1.0e-3  # in periods
MESH_TOLERANCE = 5.0e-3  # in D


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def read_field(path: Path) -> dict[str, np.ndarray]:
    """Return x, y, vx, vy, pr and the time stamp of a field file; a field that is zero stays zero."""
    from pymech.neksuite import readnek

    data = readnek(str(path))
    out = {
        "x": np.concatenate([e.pos[0].ravel() for e in data.elem]),
        "y": np.concatenate([e.pos[1].ravel() for e in data.elem]),
        "vx": np.concatenate([e.vel[0].ravel() for e in data.elem]),
        "vy": np.concatenate([e.vel[1].ravel() for e in data.elem]),
        "time": np.array(data.time),
    }
    return out


def l2(mask: np.ndarray, *components: np.ndarray) -> float:
    """Norm of the stacked components over the nodes selected by the mask."""
    return float(np.sqrt(sum(np.sum(c[mask] ** 2) for c in components)))


def main() -> None:
    """Rebuild the orbit from the filters, write bank_check.dat and draw the figure."""
    nk = load_nekplot()
    nk.configure_style()

    period = 1.0 / FORCING_FREQUENCY
    history = np.loadtxt(CASE_DIR / "dmt_bank.dat")
    time = history[:, 0]

    orbit_file = CASE_DIR / "BF_1cyl0.f00001"
    if not orbit_file.exists():
        orbit_file = CASE_DIR.parent / "210_baseflow_newton" / "BF_1cyl0.f00001"
    orbit = read_field(orbit_file)
    amp = {}
    for k in (*HARMONICS, MEAN):
        re = read_field(CASE_DIR / f"bRe1cyl0.f{k:05d}")
        im = read_field(CASE_DIR / f"bIm1cyl0.f{k:05d}")
        # The time loop takes one step past the end time, and the cylinder moves during it, so the
        # mesh of the amplitude files differs from the mesh of the orbit by less than MESH_TOLERANCE.
        phase = float(re["time"]) / period % 1.0
        if min(phase, 1.0 - phase) > PHASE_TOLERANCE:
            raise SystemExit(
                f"filter {k}: the run ended at phase {phase:.4f} of the forcing, not at phase 0"
            )
        same_mesh = np.allclose(
            re["x"], orbit["x"], rtol=0.0, atol=MESH_TOLERANCE
        ) and np.allclose(re["y"], orbit["y"], rtol=0.0, atol=MESH_TOLERANCE)
        if not same_mesh:
            raise SystemExit(
                f"filter {k}: the mesh of bRe is not the mesh of the orbit"
            )
        amp[k] = (re["vx"], re["vy"], im["vx"], im["vy"])

    x, y = orbit["x"], orbit["y"]
    near = (x > NEAR_X[0]) & (x < NEAR_X[1]) & (np.abs(y) < NEAR_Y)
    everything = np.ones_like(near)
    mean_x, mean_y = amp[MEAN][0], amp[MEAN][1]
    fluct_x, fluct_y = orbit["vx"] - mean_x, orbit["vy"] - mean_y
    norm_mean = l2(near, mean_x, mean_y)

    rows = []
    sum_x = np.zeros_like(x)
    sum_y = np.zeros_like(y)
    for m in (0, *HARMONICS):
        if m > 0:
            sum_x = sum_x + amp[m][0]
            sum_y = sum_y + amp[m][1]
        err = []
        for mask in (near, everything):
            err.append(
                l2(mask, fluct_x - sum_x, fluct_y - sum_y) / l2(mask, fluct_x, fluct_y)
            )
        a_norm = (
            1.0
            if m == 0
            else l2(near, amp[m][0], amp[m][1], amp[m][2], amp[m][3]) / norm_mean
        )
        rows.append([m, err[0], err[1], a_norm])
    np.savetxt(
        CASE_DIR / "bank_check.dat",
        np.array(rows),
        header="m  e_m_near  e_m_all  norm_uhat_m/norm_mean    (e_m = |u_orbit - mean - sum_k<=m Re uhat_k| / |u_orbit - mean|)",
        fmt=["%d", "%.8g", "%.8g", "%.8g"],
    )
    print("m  e_m(near)  e_m(all)  |uhat_m|/|mean|")
    for r in rows:
        print(f"{int(r[0])}  {r[1]:.4e}  {r[2]:.4e}  {r[3]:.4e}")

    wall = np.abs(np.hypot(x, y) - 0.5) < WALL_TOLERANCE
    wall_rows = [
        [k, float(amp[k][1][wall].mean()), float(amp[k][3][wall].mean())]
        for k in (*HARMONICS, MEAN)
    ]
    np.savetxt(
        CASE_DIR / "bank_wall.dat",
        np.array(wall_rows),
        header="k  re_uhat_y  im_uhat_y    (mean over the wall nodes; exact: re_uhat_y = A 2 pi f for k = 1, 0 for the others)",
        fmt=["%d", "%.8g", "%.8g"],
    )
    exact = AMPLITUDE * 2.0 * np.pi * FORCING_FREQUENCY
    print(f"wall: Re uhat_1y = {wall_rows[0][1]:.5f} (exact A 2 pi f = {exact:.5f})")

    fig = plt.figure(figsize=(7.2, 6.2), layout="constrained")
    grid = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.9])
    triang = nk.make_triangulation(x, y)

    for col, k in enumerate(HARMONICS):
        ax = fig.add_subplot(grid[0, 2 * col : 2 * col + 2])
        magnitude = np.sqrt(sum(c**2 for c in amp[k]))
        top = float(np.nanpercentile(magnitude[near], 99.5))
        filled = ax.tricontourf(
            triang,
            magnitude,
            levels=np.linspace(0.0, top, 25),
            cmap="viridis",
            extend="max",
        )
        nk.add_cylinder_patches(ax)
        ax.set_aspect("equal")
        ax.set_xlim(-2.0, 12.0)
        ax.set_ylim(-3.0, 3.0)
        ax.set_xlabel(r"$x$ ($D$)")
        if col == 0:
            ax.set_ylabel(r"$y$ ($D$)")
        label = ("a", "b", "c")[col]
        ax.set_title(
            rf"({label}) $|\hat u_{k}|$, $k\omega$ with $k={k}$",
            loc="left",
            fontsize=7.5,
        )
        fig.colorbar(
            filled,
            ax=ax,
            orientation="horizontal",
            shrink=0.9,
            pad=0.02,
            ticks=np.linspace(0.0, top, 3),
            format="%.2f",
        )

    ax_hist = fig.add_subplot(grid[1, 0:4])
    for k, lw in zip(HARMONICS, (3.0, 2.0, 1.2)):
        column = history[
            :, 2 + k
        ]  # |uhat_k| of the history: columns time, |u|, |u - sum b|, then the filters
        ax_hist.plot(
            time / period,
            column / column[-1],
            lw=lw,
            marker="o",
            ms=3,
            markevery=max(1, len(time) // 25),
            label=rf"$k={k}$, $\omega$ = {k * 2 * np.pi * FORCING_FREQUENCY:.3f}",
        )
    ax_hist.axvline(28.0, c="k", lw=0.7, ls="--", label="average from $28\\,T$")
    ax_hist.set_xlabel(r"$t / T$")
    ax_hist.set_ylabel(r"$\|\hat u_k\|$ / value at $t=30\,T$")
    ax_hist.set_title(
        "(d) The amplitude of each filter settles", loc="left", fontsize=7.5
    )
    ax_hist.legend(loc="lower right", fontsize=6.5)

    ax_err = fig.add_subplot(grid[1, 4:6])
    ax_err.semilogy(
        [r[0] for r in rows],
        [r[1] for r in rows],
        c="C3",
        lw=1.6,
        marker="s",
        ms=4,
        label="near field",
    )
    ax_err.semilogy(
        [r[0] for r in rows],
        [r[2] for r in rows],
        c="C0",
        lw=1.0,
        marker="o",
        ms=3,
        label="all nodes",
    )
    ax_err.set_xticks(range(0, 4))
    ax_err.set_xlabel("harmonics in the sum")
    ax_err.set_ylabel(r"$e_m$")
    ax_err.set_title("(e) Orbit from the filters", loc="left", fontsize=7.5)
    ax_err.legend(loc="upper right", fontsize=6.5)

    output = CASE_DIR / "plot_dmt_bank.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
