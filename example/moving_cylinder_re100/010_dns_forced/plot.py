#!/usr/bin/env python3
"""Plot the forced-oscillation DNS of the Re = 100 cylinder and measure its lock-in to the forcing.

The cylinder moves transversely with the amplitude A = 1.08 and the frequency f = 0.1643 (userParam05 and
userParam06). The DNS writes one field per forcing period, all at phase 0 of the forcing, and the probe
history 1cyl.his. The script reads both.

It writes lockin.dat and plot_lockin.png. lockin.dat holds one row with:
  probe_period        mean spacing of the upward zero crossings of v at the probe over the last 10 periods
  forcing_period      1 / userParam06
  last_change         ||u_n - u_(n-1)|| / ||u_n|| of the last two fields in the near field (u = the two velocity
                      components at the nodes with -2 < x < 10 and |y| < 4)
  last_change_full    the same ratio over all nodes
  probe_amplitude     half of the range of v at the probe over the last 5 periods
  n_fields            number of fields
and lockin_history.dat, the change from period n - 1 to n in the near field and over all nodes, for n = 2 to the last.

The near field measures the lock-in of the flow around the cylinder. The far wake, which the start-up transient
leaves behind, travels out of the domain slowly and keeps the change over all nodes larger.

Lengths are in units of the cylinder diameter, speeds in units of the free-stream speed, times in units of D/U.

Usage: python plot.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
LAST_PERIODS = 10
AMPLITUDE_PERIODS = 5
NEAR_X = (-2.0, 10.0)
NEAR_Y = 4.0
TOLERANCE = 1e-3  # tol_abs of last_change in ref/reference.json


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def par_value(name: str) -> float:
    """Return the number after `name =` in the .par file."""
    text = next(CASE_DIR.glob("*.par")).read_text()
    match = re.search(rf"^\s*{name}\s*=\s*([0-9.eE+-]+)", text, re.M | re.I)
    if match is None:
        raise ValueError(f"no {name} in the par file")
    return float(match.group(1))


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


def upward_crossings(t: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Return the times at which v crosses zero upward, by linear interpolation."""
    i = np.flatnonzero((v[:-1] < 0.0) & (v[1:] >= 0.0))
    return t[i] - v[i] * (t[i + 1] - t[i]) / (v[i + 1] - v[i])


def read_velocity(nk, path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the node coordinates x, y and the velocity vector (u values, then v values) of a field file."""
    x, y, fields, _time = nk.read_field(str(path))
    return x, y, np.concatenate([fields["vx"], fields["vy"]])


def relative_change(fields: list[np.ndarray], keep: np.ndarray) -> np.ndarray:
    """Return ||u_n - u_(n-1)|| / ||u_n|| over the nodes in keep, for n = 2 to the number of fields."""
    both = np.concatenate([keep, keep])
    return np.array(
        [
            np.linalg.norm((b - a)[both]) / np.linalg.norm(b[both])
            for a, b in zip(fields[:-1], fields[1:])
        ]
    )


def main() -> None:
    """Measure the lock-in, write lockin.dat and draw the probe signal, the period-to-period change and the last field."""
    nk = load_nekplot()
    nk.configure_style()

    forcing_period = 1.0 / par_value("userParam06")
    t_probe, v_probe = read_probe(CASE_DIR / "1cyl.his")
    t_end = t_probe[-1]

    window = t_probe >= t_end - LAST_PERIODS * forcing_period
    crossings = upward_crossings(t_probe[window], v_probe[window])
    probe_period = float((crossings[-1] - crossings[0]) / (crossings.size - 1))
    late = t_probe >= t_end - AMPLITUDE_PERIODS * forcing_period
    amplitude = 0.5 * float(v_probe[late].max() - v_probe[late].min())

    files = sorted(CASE_DIR.glob("1cyl0.f*"))
    x, y, _ = read_velocity(nk, files[0])
    near = (x > NEAR_X[0]) & (x < NEAR_X[1]) & (np.abs(y) < NEAR_Y)
    fields = [read_velocity(nk, f)[2] for f in files]
    change_near = relative_change(fields, near)
    change_full = relative_change(fields, np.ones_like(near))
    np.savetxt(
        CASE_DIR / "lockin.dat",
        [
            [
                probe_period,
                forcing_period,
                change_near[-1],
                change_full[-1],
                amplitude,
                len(files),
            ]
        ],
        header="probe_period  forcing_period  last_change  last_change_full  probe_amplitude  n_fields",
        fmt="%.8g",
    )
    np.savetxt(
        CASE_DIR / "lockin_history.dat",
        np.column_stack([np.arange(2, len(files) + 1), change_near, change_full]),
        header="period n  change_near  change_full   (||u_n - u_(n-1)|| / ||u_n||)",
        fmt=["%d", "%.8g", "%.8g"],
    )
    print(
        f"probe period {probe_period:.5f} (forcing {forcing_period:.5f}); last change {change_near[-1]:.3e} "
        f"in the near field, {change_full[-1]:.3e} over all nodes; probe amplitude {amplitude:.4f}; "
        f"{len(files)} fields"
    )

    fig = plt.figure(figsize=(7.0, 5.0), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0], width_ratios=[1.0, 1.5])
    ax_probe = fig.add_subplot(grid[0, :])
    ax_change = fig.add_subplot(grid[1, 0])
    ax_field = fig.add_subplot(grid[1, 1])

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
        t_end - LAST_PERIODS * forcing_period,
        c="C0",
        lw=1.0,
        ls="--",
        label=f"last {LAST_PERIODS} periods: probe period {probe_period:.4f}",
    )
    ax_probe.set_ylim(1.5 * v_probe.min(), 1.5 * v_probe.max())
    ax_probe.set_xlabel(r"$t$ ($D/U$)")
    ax_probe.set_ylabel(r"$v$ ($U$)")
    ax_probe.set_title("(a) Probe signal", loc="left", fontsize=8)
    ax_probe.legend(loc="upper center", fontsize=6.5, ncol=2)

    periods = np.arange(2, len(files) + 1)
    ax_change.semilogy(
        periods,
        change_full,
        c="0.55",
        lw=2.2,
        marker="s",
        ms=3,
        label="all nodes",
    )
    ax_change.semilogy(
        periods,
        change_near,
        c="C0",
        lw=1.2,
        marker="o",
        ms=3,
        label=f"near field ({change_near[-1]:.1e})",
    )
    ax_change.axhline(TOLERANCE, c="k", lw=0.6, ls=":", label=f"{TOLERANCE:.0e}")
    ax_change.set_xlabel("forcing period $n$")
    ax_change.set_ylabel(r"$\|u_n - u_{n-1}\| / \|u_n\|$")
    ax_change.set_title("(b) Change from period to period", loc="left", fontsize=8)
    ax_change.legend(loc="upper right", fontsize=6.5)

    x, y, last, time = nk.read_field(str(files[-1]))
    vorticity = last["t"]
    bound = float(np.nanpercentile(np.abs(vorticity), 99.0))
    filled = ax_field.tricontourf(
        nk.make_triangulation(x, y),
        vorticity,
        levels=np.linspace(-bound, bound, 25),
        cmap="RdBu_r",
        extend="both",
    )
    nk.add_cylinder_patches(ax_field)
    ax_field.set_aspect("equal")
    ax_field.set_xlim(-2.0, 12.0)
    ax_field.set_ylim(-3.0, 3.0)
    ax_field.set_xlabel(r"$x$ ($D$)")
    ax_field.set_ylabel(r"$y$ ($D$)")
    ax_field.set_title(
        f"(c) Vorticity at $t = {time:.2f}$, phase 0", loc="left", fontsize=8
    )
    colorbar = fig.colorbar(
        filled, ax=ax_field, label=r"$\omega_z$ ($U/D$)", shrink=0.8
    )
    colorbar.set_ticks([-4.0, 0.0, 4.0])

    output = CASE_DIR / "plot_lockin.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
