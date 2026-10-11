#!/usr/bin/env python3
"""Plot the energy budget of the leading Floquet mode of the cubic cavity at Re = 1950.

Mode 4.11 marches the real and the imaginary part of the leading Floquet mode of
../311_stability_direct_floquet over one period T = 10.78909 of the limit cycle, on a stored copy of the
orbit, and averages the terms of the energy equation with the weight exp(-2 sigma t) dt / T. The run writes
PKE_floquet_cav0.f00001, one number per line: lines 1 to 9 are the integrals of the nine production fields,
line 10 is minus the dissipation, line 11 the production P, line 12 the net P - D, line 13 the energy E, line
14 the growth rate sigma from the budget, line 15 the sponge damping (zero here) and line 16 the growth rate
from the energy of the marched mode, ln(E(T) / E(0)) / (2 T). It also writes the orbit-averaged production
fields F01cav0.f00001 to F01cav0.f00003: one file for each direction of the gradient of the base flow, with the
three momentum components in vx, vy and vz.

The script compares the three growth rates and writes budget_closure.dat with one row: sigma of the Floquet
stage (row 1 of Spectre_NSd_conv.dat), sigma from the budget, sigma from the energy of the marched mode, the
relative difference of each to the first, P, D, E, and the gap of the balance P - D = sigma E in units of D:
(sigma from the budget - sigma of the Floquet stage) E / D. It writes plot_budget.png: the production density, the
sum of the nine fields, on the three mid-planes, and the three growth rates.

Lengths are in units of the cavity edge L, speeds in units of the lid speed U, time in L/U.

Usage: python plot.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
from pymech.neksuite import readnek

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
PLANES = (
    ("z", 0, 1),
    ("y", 0, 2),
    ("x", 1, 2),
)  # normal axis name, and the two axes drawn


def read_summary(path: Path) -> np.ndarray:
    """Return the numbers of PKE_floquet_cav0.f00001, one per line."""
    return np.loadtxt(path)


def production_density() -> tuple[np.ndarray, np.ndarray]:
    """Return the node positions (n, 3) and the sum of the nine production fields at them."""
    total = None
    positions = None
    for number in (1, 2, 3):
        data = readnek(str(CASE_DIR / f"F01cav0.f{number:05d}"))
        part = np.concatenate(
            [e.vel.sum(axis=0).ravel() for e in data.elem]
        )  # vx + vy + vz of the three components
        if positions is None:
            positions = np.stack(
                [
                    np.concatenate([e.pos[i].ravel() for e in data.elem])
                    for i in range(3)
                ],
                axis=1,
            )
        total = part if total is None else total + part
    return positions, total


def plane_values(
    positions: np.ndarray, values: np.ndarray, normal: int, axes: tuple[int, int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the two coordinates and the values at the nodes of the mid-plane normal to one axis.

    A node on an element edge belongs to several elements; its entries are averaged.
    """
    on_plane = np.isclose(positions[:, normal], 0.0, atol=1e-9)
    keys = np.round(positions[on_plane][:, list(axes)], 9)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    average = np.bincount(inverse, weights=values[on_plane]) / np.bincount(inverse)
    return unique[:, 0], unique[:, 1], average


def main() -> None:
    """Write budget_closure.dat and draw the figure."""
    summary = read_summary(CASE_DIR / "PKE_floquet_cav0.f00001")
    production, net, energy = summary[10], summary[11], summary[12]
    sigma_budget, sigma_energy = summary[13], summary[15]
    dissipation = production - net
    sigma_floquet = float(np.loadtxt(CASE_DIR / "Spectre_NSd_conv.dat", ndmin=2)[0, 0])
    closure = (sigma_budget - sigma_floquet) / sigma_floquet
    energy_gap = (sigma_energy - sigma_floquet) / sigma_floquet
    # The balance P - D = sigma E has the error (sigma_budget - sigma_floquet) E; in units of D it does not depend on
    # how small sigma is.
    imbalance = (sigma_budget - sigma_floquet) * energy / dissipation
    np.savetxt(
        CASE_DIR / "budget_closure.dat",
        np.array(
            [
                [
                    sigma_floquet,
                    sigma_budget,
                    sigma_energy,
                    closure,
                    energy_gap,
                    production,
                    dissipation,
                    energy,
                    imbalance,
                ]
            ]
        ),
        header=(
            "sigma_floquet  sigma_budget  sigma_energy  (budget-floquet)/floquet  (energy-floquet)/floquet  P  D  E  "
            "(budget-floquet)*E/D"
        ),
        fmt="%.8g",
    )
    print(
        f"sigma: Floquet {sigma_floquet:.6e}, budget {sigma_budget:.6e} ({closure:+.2%}), "
        f"marched energy {sigma_energy:.6e} ({energy_gap:+.2%})"
    )
    print(
        f"P = {production:.6e}, D = {dissipation:.6e}, E = {energy:.6e}; "
        f"the balance P - D - sigma E is {imbalance:.2%} of D"
    )

    positions, density = production_density()
    density = density / np.percentile(np.abs(density), 99.9)

    figure = plt.figure(figsize=(7.2, 5.4), layout="constrained")
    grid = figure.add_gridspec(2, 3, height_ratios=[1.25, 0.8])
    labels = ("x", "y", "z")
    for index, (name, first, second) in enumerate(PLANES):
        normal = labels.index(name)
        axes = figure.add_subplot(grid[0, index])
        u, v, values = plane_values(positions, density, normal, (first, second))
        filled = axes.tricontourf(
            u,
            v,
            values,
            levels=np.linspace(-1.0, 1.0, 25),
            cmap="RdBu_r",
            extend="both",
        )
        axes.set_aspect("equal")
        axes.set_xlabel(rf"${labels[first]}$ ($L$)")
        axes.set_ylabel(rf"${labels[second]}$ ($L$)")
        axes.set_title(rf"({'abc'[index]}) ${name}=0$", loc="left", fontsize=7.5)
        figure.colorbar(
            filled,
            ax=axes,
            orientation="horizontal",
            shrink=0.9,
            pad=0.02,
            ticks=[-1.0, 0.0, 1.0],
            format="%.0f",
        )

    # The growth rate is the small difference of two larger terms, P/E and D/E; the legend gives both.
    axes = figure.add_subplot(grid[1, :])
    axes.axhline(
        sigma_floquet,
        c="k",
        lw=1.0,
        ls="--",
        label=rf"$\sigma$ of the Floquet stage: {sigma_floquet:.4e}",
    )
    axes.plot(
        [1],
        [sigma_budget],
        c="C3",
        lw=0.0,
        marker="o",
        ms=7,
        label=(
            rf"budget $(\bar P-\bar D)/\bar E$: {sigma_budget:.4e} ({closure:+.1%})"
            "\n"
            rf"$\bar P/\bar E$ = {production / energy:.5f}, $\bar D/\bar E$ = "
            rf"{dissipation / energy:.5f}; the gap is {imbalance:.2%} of $\bar D/\bar E$"
        ),
    )
    axes.plot(
        [2],
        [sigma_energy],
        c="C0",
        lw=0.0,
        marker="s",
        ms=7,
        label=rf"energy of the marched mode: {sigma_energy:.4e} ({energy_gap:+.1%})",
    )
    axes.set_xlim(0.5, 2.5)
    axes.set_xticks([1, 2], ["budget", "energy"])
    axes.set_ylabel(r"$\sigma$ ($U/L$)")
    axes.set_title("(d) Growth rate by two routes", loc="left", fontsize=7.5)
    low, high = axes.get_ylim()
    axes.set_ylim(low - 0.05 * (high - low), high + 0.5 * (high - low))
    axes.legend(loc="upper center", fontsize=6.5)

    output = CASE_DIR / "plot_budget.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
