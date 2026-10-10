#!/usr/bin/env python3
"""Plot the wavemaker and the leading direct and adjoint modes of the cubic cavity.

The script reads the wavemaker wm_cav0.f00001 of this stage and the leading
modes dRecav0.f00001, dImcav0.f00001, aRecav0.f00001 and aImcav0.f00001 that
the run read. Each field is shown on the plane z = 0 and scaled to its maximum
on that plane. Lengths are in units of the cavity edge L. It writes
plot_wavemaker.png.

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


def reynolds_number(par_file: Path) -> float:
    """Return the Reynolds number: the viscosity in the par file is -Re."""
    found = re.search(r"^\s*viscosity\s*=\s*-?([0-9.eE+]+)", par_file.read_text(), re.M)
    if found is None:
        raise ValueError(f"no viscosity in {par_file}")
    return float(found.group(1))


def midplane(
    field_file: Path, imag_file: Path | None = None, scalar: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and a field amplitude on the plane z = 0.

    A file with a velocity field gives the root of the sum of the squared
    velocity components. With a second file, which holds the imaginary part
    of a mode, the squares of both parts are added. With scalar = True the
    function returns the scalar of the first temperature slot, where the run
    stores the wavemaker. The plane z = 0 is an element face of the cavity
    mesh. Only the elements that start at z = 0 contribute, so each node
    appears once per element.
    """
    first = readnek(str(field_file))
    second = readnek(str(imag_file)) if imag_file is not None else first
    x_nodes, y_nodes, values = [], [], []
    for element, other in zip(first.elem, second.elem, strict=True):
        if np.abs(element.pos[2][0]).max() > 1e-9:
            continue
        if scalar:
            amplitude = element.temp[0][0]
        else:
            squares = sum(element.vel[k][0] ** 2 for k in range(3))
            if imag_file is not None:
                squares = squares + sum(other.vel[k][0] ** 2 for k in range(3))
            amplitude = np.sqrt(squares)
        x_nodes.append(element.pos[0][0].ravel())
        y_nodes.append(element.pos[1][0].ravel())
        values.append(amplitude.ravel())
    return np.concatenate(x_nodes), np.concatenate(y_nodes), np.concatenate(values)


def draw(ax: plt.Axes, data: tuple[np.ndarray, np.ndarray, np.ndarray], title: str):
    """Draw one scaled field on the axes and return the filled contour set."""
    x, y, values = data
    values = values / values.max()
    triangulation = Triangulation(
        x, y, triangles=Delaunay(np.column_stack((x, y))).simplices
    )
    filled = ax.tricontourf(
        triangulation, values, levels=np.linspace(0.0, 1.0, 21), cmap="Blues"
    )
    ax.set_aspect("equal")
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_xlabel(r"$x$ ($L$)")
    ax.set_title(title, loc="left", fontsize=8)
    return filled


def main() -> None:
    """Draw the wavemaker (a), the direct mode (b) and the adjoint mode (c) on z = 0."""
    reynolds = reynolds_number(CASE_DIR / "cav.par")
    fig, axes = plt.subplots(
        1, 3, figsize=(7.0, 2.6), layout="constrained", sharey=True
    )
    draw(
        axes[0],
        midplane(CASE_DIR / "wm_cav0.f00001", scalar=True),
        f"(a) Wavemaker, Re = {reynolds:.0f}",
    )
    draw(
        axes[1],
        midplane(CASE_DIR / "dRecav0.f00001", CASE_DIR / "dImcav0.f00001"),
        "(b) Direct mode",
    )
    filled = draw(
        axes[2],
        midplane(CASE_DIR / "aRecav0.f00001", CASE_DIR / "aImcav0.f00001"),
        "(c) Adjoint mode",
    )
    axes[0].set_ylabel(r"$y$ ($L$)")
    colorbar = fig.colorbar(filled, ax=axes, label="amplitude (maximum = 1)")
    colorbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])

    output = CASE_DIR / "plot_wavemaker.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
