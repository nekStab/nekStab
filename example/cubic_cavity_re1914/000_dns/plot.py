#!/usr/bin/env python3
"""Plot the probe signal and the mid-plane speed of the cubic cavity DNS.

The script reads the probe history cav.his and the last snapshot cav0.f0*
of this stage. Lengths are in units of the cavity edge L, speeds in units of
the lid speed U, time in units of L/U. It writes plot_dns.png.

Usage: python plot.py
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from pymech.neksuite import readnek
from scipy.spatial import Delaunay

CASE_DIR = Path(__file__).resolve().parent
PROBE_HEADER_LINES = 2  # point count and probe position
START_OF_PLOT = 10.0  # the start-up transient before this time hides the later signal


def reynolds_number(par_file: Path) -> float:
    """Return the Reynolds number: the viscosity in the par file is -Re."""
    found = re.search(r"^\s*viscosity\s*=\s*-?([0-9.eE+]+)", par_file.read_text(), re.M)
    if found is None:
        raise ValueError(f"no viscosity in {par_file}")
    return float(found.group(1))


def midplane_speed(field_file: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and the in-plane speed on the plane z = 0 of a field file.

    The plane z = 0 is an element face of the cavity mesh. Only the elements
    that start at z = 0 contribute, so each node appears once per element.
    """
    field = readnek(str(field_file))
    x_nodes, y_nodes, speed = [], [], []
    for element in field.elem:
        z_first = element.pos[2][0]
        if np.abs(z_first).max() > 1e-9:
            continue
        x_nodes.append(element.pos[0][0].ravel())
        y_nodes.append(element.pos[1][0].ravel())
        speed.append(np.hypot(element.vel[0][0], element.vel[1][0]).ravel())
    return np.concatenate(x_nodes), np.concatenate(y_nodes), np.concatenate(speed)


def main() -> None:
    """Draw the probe history (a) and the speed on the mid-plane z = 0 (b)."""
    reynolds = reynolds_number(CASE_DIR / "cav.par")
    probe = np.loadtxt(CASE_DIR / "cav.his", skiprows=PROBE_HEADER_LINES)
    snapshot = sorted(CASE_DIR.glob("cav0.f0*"))[-1]
    x, y, speed = midplane_speed(snapshot)
    triangulation = Triangulation(
        x, y, triangles=Delaunay(np.column_stack((x, y))).simplices
    )

    fig, (ax_probe, ax_field) = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.0),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.2, 1.0]},
    )
    ax_probe.plot(
        probe[:, 0],
        probe[:, 2],
        "-o",
        c="k",
        lw=1.0,
        ms=3,
        markevery=max(len(probe) // 20, 1),
        label="probe",
    )
    ax_probe.axhline(probe[-1, 2], c="0.5", lw=0.8, ls="--", label="final value")
    ax_probe.set_xlim(START_OF_PLOT, probe[-1, 0])
    late = probe[:, 0] >= START_OF_PLOT
    ax_probe.set_ylim(probe[late, 2].min() - 0.002, probe[late, 2].max() + 0.002)
    ax_probe.set_xlabel(r"time $t$ ($L/U$)")
    ax_probe.set_ylabel(r"$v$ at $(0.25, 0.25, 0)$ ($U$)")
    ax_probe.set_title(f"(a) Probe signal, Re = {reynolds:.0f}", loc="left", fontsize=8)
    ax_probe.legend(loc="best")

    levels = np.linspace(0.0, 1.0, 21)
    filled = ax_field.tricontourf(
        triangulation, speed, levels=levels, cmap="Blues", extend="max"
    )
    colorbar = fig.colorbar(filled, ax=ax_field, label=r"in-plane speed ($U$)")
    colorbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_field.set_aspect("equal")
    ax_field.set_xticks([-0.5, 0.0, 0.5])
    ax_field.set_xlabel(r"$x$ ($L$)")
    ax_field.set_ylabel(r"$y$ ($L$)")
    ax_field.set_title("(b) Plane z = 0, last snapshot", loc="left", fontsize=8)

    output = CASE_DIR / "plot_dns.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
