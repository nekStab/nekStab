#!/usr/bin/env python3
"""Plot the Newton residual history and the mid-plane speed of the cavity base flow.

The script reads residu_newton.dat and the base flow BF_cav0.f00001 of this
stage. Lengths are in units of the cavity edge L, speeds in units of the lid
speed U. It writes plot_newton.png.

Usage: python plot.py
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib.tri import Triangulation
from pymech.neksuite import readnek
from scipy.spatial import Delaunay

CASE_DIR = Path(__file__).resolve().parent


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
    """Draw the Newton residual (a) and the speed on the mid-plane z = 0 (b)."""
    reynolds = reynolds_number(CASE_DIR / "cav.par")
    history = np.atleast_2d(np.loadtxt(CASE_DIR / "residu_newton.dat"))
    x, y, speed = midplane_speed(CASE_DIR / "BF_cav0.f00001")
    triangulation = Triangulation(
        x, y, triangles=Delaunay(np.column_stack((x, y))).simplices
    )

    fig, (ax_res, ax_field) = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.0),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.2, 1.0]},
    )
    ax_res.plot(
        history[:, 0], history[:, 6], "-o", c="k", lw=1.0, ms=4, label="residual"
    )
    ax_res.axhline(history[-1, 7], c="0.5", lw=0.8, ls="--", label="tolerance")
    ax_res.set_yscale("log")
    ax_res.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax_res.set_xlabel("Newton iteration")
    ax_res.set_ylabel(r"squared residual $\|F(q) - q\|^2$")
    ax_res.set_title(
        f"(a) Newton residual, Re = {reynolds:.0f}", loc="left", fontsize=8
    )
    ax_res.legend(loc="best")

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
    ax_field.set_title("(b) Base flow, plane z = 0", loc="left", fontsize=8)

    output = CASE_DIR / "plot_newton.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
