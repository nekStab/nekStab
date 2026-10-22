#!/usr/bin/env python3
"""Plot the optimal transient growth of the cubic cavity at Re = 1914 and check its optimal response.

The script reads the gains of this stage (Spectre_Hp_conv.dat), the horizon tau of cav.par, the gains
of the other horizons (growth_sweep.csv), the optimal initial perturbation pRecav0.f00001 and its
response orecav0.f00001. It computes the kinetic energy of both fields with GLL quadrature on the
brick elements of the mesh. The energy of the response divided by the energy of the perturbation is
the gain of that perturbation, and equals G when the adjoint is consistent with the direct map.
The script writes growth_check.dat and plot_growth.png. Lengths are in units of the cavity edge L,
speeds in units of the lid speed U, times in units of L/U.

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
from numpy.polynomial.legendre import Legendre
from pymech.neksuite import readnek
from scipy.spatial import Delaunay

CASE_DIR = Path(__file__).resolve().parent


def horizon(par_file: Path) -> float:
    """Return the horizon tau of the optimisation: the end time of the par file."""
    found = re.search(r"^\s*endTime\s*=\s*([0-9.eE+-]+)", par_file.read_text(), re.M)
    if found is None:
        raise ValueError(f"no endTime in {par_file}")
    return float(found.group(1))


def gll_weights(n_points: int) -> np.ndarray:
    """Return the Gauss-Lobatto-Legendre quadrature weights on [-1, 1] for n_points nodes."""
    degree = n_points - 1
    legendre = Legendre.basis(degree)
    nodes = np.concatenate(([-1.0], np.sort(legendre.deriv().roots().real), [1.0]))
    return 2.0 / (degree * (degree + 1) * legendre(nodes) ** 2)


def kinetic_energy(field_file: Path) -> float:
    """Return the integral of the squared velocity over the cavity.

    Each element is a brick, so its Jacobian is constant: one eighth of its volume.
    """
    total = 0.0
    for element in readnek(str(field_file)).elem:
        n_z, n_y, n_x = element.pos.shape[1:]
        volume = np.prod([np.ptp(element.pos[axis]) for axis in range(3)])
        weight = (
            gll_weights(n_z)[:, None, None]
            * gll_weights(n_y)[None, :, None]
            * gll_weights(n_x)[None, None, :]
            * volume
            / 8.0
        )
        total += float(np.sum(weight * np.sum(element.vel**2, axis=0)))
    return total


def midplane_speed(field_file: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and the speed on the plane z = 0 of a field file.

    The plane z = 0 is an element face of the cavity mesh. Only the elements
    that start at z = 0 contribute, so each node appears once per element.
    """
    x_nodes, y_nodes, speed = [], [], []
    for element in readnek(str(field_file)).elem:
        if np.abs(element.pos[2][0]).max() > 1e-9:
            continue
        x_nodes.append(element.pos[0][0].ravel())
        y_nodes.append(element.pos[1][0].ravel())
        speed.append(np.sqrt(np.sum(element.vel[:, 0] ** 2, axis=0)).ravel())
    return np.concatenate(x_nodes), np.concatenate(y_nodes), np.concatenate(speed)


def triangulation(x: np.ndarray, y: np.ndarray) -> Triangulation:
    """Return the Delaunay triangulation of the points (x, y)."""
    return Triangulation(x, y, triangles=Delaunay(np.column_stack((x, y))).simplices)


def main() -> None:
    """Check the optimal response, write growth_check.dat and draw G(tau) and the two fields."""
    tau = horizon(CASE_DIR / "cav.par")
    gain = float(np.loadtxt(CASE_DIR / "Spectre_Hp_conv.dat", ndmin=2)[0, 0])
    perturbation = sorted(CASE_DIR.glob("pRecav0.f*"))[0]
    response = sorted(CASE_DIR.glob("orecav0.f*"))[0]
    energy_in = kinetic_energy(perturbation)
    energy_out = kinetic_energy(response)
    ratio = energy_out / energy_in
    print(f"tau = {tau:g}: G = {gain:.4f}")
    print(
        f"energy of the perturbation {energy_in:.6f}, of the response {energy_out:.4f}"
    )
    print(
        f"gain of the stored perturbation {ratio:.4f}; divided by G: {ratio / gain:.6f}"
    )
    np.savetxt(
        CASE_DIR / "growth_check.dat",
        [[tau, gain, energy_in, energy_out, ratio, ratio / gain]],
        header="tau  G  energy_perturbation  energy_response  energy_ratio  ratio_over_G",
        fmt="%.8g",
    )

    sweep_tau, sweep_gain = np.loadtxt(
        CASE_DIR / "growth_sweep.csv", delimiter=",", skiprows=3, ndmin=2, unpack=True
    )
    fig, (ax_gain, ax_in, ax_out) = plt.subplots(
        1,
        3,
        figsize=(7.0, 2.6),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]},
    )
    ax_gain.plot(
        sweep_tau,
        sweep_gain,
        "s-",
        c="0.6",
        lw=2.0,
        ms=5,
        label="one run per horizon",
    )
    ax_gain.plot([tau], [gain], "o", c="k", ms=3, label=f"this run, G = {gain:.1f}")
    ax_gain.axvline(tau, c="0.5", lw=0.8, ls="--", label=r"horizon $\tau$ of this run")
    ax_gain.set_yscale("log")
    ax_gain.set_xlabel(r"horizon $\tau$ ($L/U$)")
    ax_gain.set_ylabel(r"optimal gain $G(\tau)$")
    ax_gain.set_title("(a) Optimal gain, Re = 1914", loc="left", fontsize=8)
    ax_gain.legend(loc="upper left", fontsize=6.5)

    levels = np.linspace(0.0, 1.0, 21)
    for ax, field_file, title in (
        (ax_in, perturbation, "(b) Perturbation, t = 0"),
        (ax_out, response, rf"(c) Response, t = {tau:g}"),
    ):
        x, y, speed = midplane_speed(field_file)
        filled = ax.tricontourf(
            triangulation(x, y), speed / speed.max(), levels=levels, cmap="Blues"
        )
        ax.set_aspect("equal")
        ax.set_xticks([-0.5, 0.0, 0.5])
        ax.set_xlabel(r"$x$ ($L$)")
        ax.set_title(title, loc="left", fontsize=8)
    ax_in.set_ylabel(r"$y$ ($L$)")
    colorbar = fig.colorbar(
        filled, ax=[ax_in, ax_out], label="speed / maximum, plane z = 0"
    )
    colorbar.set_ticks([0.0, 0.5, 1.0])

    output = CASE_DIR / "plot_growth.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
