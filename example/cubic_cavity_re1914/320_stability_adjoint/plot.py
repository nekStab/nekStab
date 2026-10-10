#!/usr/bin/env python3
"""Plot the adjoint spectrum and the leading adjoint mode of the cubic cavity.

The script reads the converged eigenvalues Spectre_NSa_conv.dat and the leading
mode aRecav0.f00001 (real part) and aImcav0.f00001 (imaginary part) of this
stage. Growth rates and frequencies are in units of U/L, lengths in units of
the cavity edge L. The mode amplitude is scaled to its maximum on the plane
z = 0. It writes plot_adjoint.png.

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


def midplane_amplitude(
    real_file: Path, imag_file: Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and the mode amplitude on the plane z = 0.

    The amplitude is the root of the sum of the squared real and imaginary
    parts of the three velocity components. It does not depend on the phase of
    the mode. The plane z = 0 is an element face of the cavity mesh. Only the
    elements that start at z = 0 contribute, so each node appears once per element.
    """
    real_part = readnek(str(real_file))
    imag_part = readnek(str(imag_file))
    x_nodes, y_nodes, amplitude = [], [], []
    for element_re, element_im in zip(real_part.elem, imag_part.elem, strict=True):
        if np.abs(element_re.pos[2][0]).max() > 1e-9:
            continue
        squares = sum(
            element_re.vel[k][0] ** 2 + element_im.vel[k][0] ** 2 for k in range(3)
        )
        x_nodes.append(element_re.pos[0][0].ravel())
        y_nodes.append(element_re.pos[1][0].ravel())
        amplitude.append(np.sqrt(squares).ravel())
    return np.concatenate(x_nodes), np.concatenate(y_nodes), np.concatenate(amplitude)


def main() -> None:
    """Draw the spectrum (a) and the leading mode on the plane z = 0 (b)."""
    reynolds = reynolds_number(CASE_DIR / "cav.par")
    spectrum = np.atleast_2d(np.loadtxt(CASE_DIR / "Spectre_NSa_conv.dat"))
    x, y, amplitude = midplane_amplitude(
        CASE_DIR / "aRecav0.f00001", CASE_DIR / "aImcav0.f00001"
    )
    amplitude = amplitude / amplitude.max()
    triangulation = Triangulation(
        x, y, triangles=Delaunay(np.column_stack((x, y))).simplices
    )

    fig, (ax_spec, ax_mode) = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.0),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.0, 1.0]},
    )
    ax_spec.axvline(0.0, c="0.5", lw=0.8, ls="--", label=r"$\sigma = 0$")
    ax_spec.plot(
        spectrum[:, 0], spectrum[:, 1], "o", c="k", ms=4, ls="none", label="adjoint"
    )
    ax_spec.set_xlabel(r"growth rate $\sigma$ ($U/L$)")
    ax_spec.set_ylabel(r"frequency $\omega$ ($U/L$)")
    ax_spec.set_title(
        f"(a) Adjoint spectrum, Re = {reynolds:.0f}", loc="left", fontsize=8
    )
    ax_spec.legend(loc="best")

    levels = np.linspace(0.0, 1.0, 21)
    filled = ax_mode.tricontourf(triangulation, amplitude, levels=levels, cmap="Blues")
    colorbar = fig.colorbar(filled, ax=ax_mode, label="mode amplitude (maximum = 1)")
    colorbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_mode.set_aspect("equal")
    ax_mode.set_xticks([-0.5, 0.0, 0.5])
    ax_mode.set_xlabel(r"$x$ ($L$)")
    ax_mode.set_ylabel(r"$y$ ($L$)")
    ax_mode.set_title("(b) Leading adjoint mode, plane z = 0", loc="left", fontsize=8)

    output = CASE_DIR / "plot_adjoint.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
