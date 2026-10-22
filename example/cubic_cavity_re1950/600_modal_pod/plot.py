#!/usr/bin/env python3
"""Plot the POD of the limit cycle of the cubic cavity at Re = 1950.

The stage reads 100 snapshots of the DNS of ../000_dns, one every 0.5 time units, from t = 1740 to 1789.5.
The snapshots lie on the saturated cycle of period T = 10.78909 (../000_dns_period). The files are:
pod_spectrum.dat (mode, eigenvalue, energy in percent, cumulative energy in percent), pod_coefficients.dat
(snapshot number and the temporal coefficients a_1 to a_10) and podcav0.f00001, podcav0.f00002 (spatial modes 1
and 2, normalised to unit norm).

For a cycle the two leading modes form a pair: their coefficients (a_1, a_2) turn once per period, and the next
pairs, (a_3, a_4) to (a_9, a_10), turn at the second to the fifth harmonic. The script measures the angular
frequency of each pair as the rate of the unwrapped phase atan2(a_2, a_1) over the snapshots, by a linear fit,
and divides it by 2 pi / T. It writes pod_summary.dat with one row: E_1, E_2, E_1 + E_2 (in percent) and the
ratios of the five pairs. It writes plot_pod.png.

Lengths are in units of the cavity edge L, speeds in units of the lid speed U, time in L/U. The modes have no unit.

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
CYCLE_PERIOD = 10.78909  # 000_dns_period
SNAPSHOT_DT = 0.5  # writeInterval of 000_dns, modal_dt of cav.usr
MODES_SHOWN = 2  # spatial modes drawn
MODES_LISTED = 12  # modes of the spectrum drawn; the rest is rounding error


def mid_plane(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and the speed |u| at the nodes of the plane z = 0 of a field file, one entry per distinct point.

    The plane is a face of the elements whose lower edge is z = 0. A node on an element edge belongs to several
    elements; its entries are averaged.
    """
    data = readnek(str(path))
    x_all, y_all, z_all, speed = [], [], [], []
    for element in data.elem:
        x_all.append(element.pos[0].ravel())
        y_all.append(element.pos[1].ravel())
        z_all.append(element.pos[2].ravel())
        speed.append(np.sqrt(np.sum(element.vel**2, axis=0)).ravel())
    x, y, z, speed = (np.concatenate(a) for a in (x_all, y_all, z_all, speed))
    on_plane = np.isclose(z, 0.0, atol=1e-9)
    keys = np.round(np.stack([x[on_plane], y[on_plane]], axis=1), 9)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    average = np.bincount(inverse, weights=speed[on_plane]) / np.bincount(inverse)
    return unique[:, 0], unique[:, 1], average


def main() -> None:
    """Measure the pair frequency, write pod_summary.dat and draw the figure."""
    spectrum = np.loadtxt(
        CASE_DIR / "pod_spectrum.dat", comments="#", usecols=(0, 1, 2)
    )
    coefficients = np.loadtxt(CASE_DIR / "pod_coefficients.dat", comments="#")
    energy = spectrum[:, 2]
    a_1, a_2 = coefficients[:, 1], coefficients[:, 2]

    time = np.arange(len(a_1)) * SNAPSHOT_DT
    ratios = []
    for first in (
        1,
        3,
        5,
        7,
        9,
    ):  # column of a_1, a_3, ... a_9: the first mode of each pair
        phase = np.unwrap(
            np.arctan2(coefficients[:, first + 1], coefficients[:, first])
        )
        rate = abs(np.polyfit(time, phase, 1)[0])
        ratios.append(rate / (2.0 * np.pi / CYCLE_PERIOD))
    ratio = ratios[0]
    np.savetxt(
        CASE_DIR / "pod_summary.dat",
        np.array([[energy[0], energy[1], energy[0] + energy[1], *ratios]]),
        header="E_1  E_2  E_1+E_2 (percent of the energy)  omega_pair / (2 pi / T) for the pairs (1,2) (3,4) (5,6) (7,8) (9,10)",
        fmt="%.8g",
    )
    print(
        f"energy {energy[0]:.3f} % and {energy[1]:.3f} %, sum {energy[0] + energy[1]:.3f} %; "
        f"the pair turns at {ratio:.5f} of the frequency of the cycle"
    )

    figure = plt.figure(figsize=(7.2, 6.4), layout="constrained")
    grid = figure.add_gridspec(2, 2)

    axes = figure.add_subplot(grid[0, 0])
    shown = slice(0, MODES_LISTED)
    axes.semilogy(
        spectrum[shown, 0],
        energy[shown],
        c="C0",
        lw=1.4,
        marker="o",
        ms=4,
        label="mode",
    )
    axes.set_xticks(range(2, MODES_LISTED + 1, 2))
    axes.set_xlabel("mode")
    axes.set_ylabel("energy (%)")
    axes.set_title(
        f"(a) Energy of the POD modes, {energy[0] + energy[1]:.1f} % in two",
        loc="left",
        fontsize=7.5,
    )
    axes.legend(loc="upper right", fontsize=6.5)

    axes = figure.add_subplot(grid[0, 1])
    axes.plot(a_1, a_2, c="C3", lw=1.0, marker="o", ms=2.5)
    axes.set_aspect("equal", adjustable="datalim")
    axes.locator_params(axis="x", nbins=4)
    axes.set_xlabel(r"$a_1$")
    axes.set_ylabel(r"$a_2$")
    axes.set_title(rf"(b) Pair turns at {ratio:.4f} $2\pi/T$", loc="left", fontsize=7.5)

    for index in range(MODES_SHOWN):
        axes = figure.add_subplot(grid[1, index])
        x, y, speed = mid_plane(CASE_DIR / f"podcav0.f{index + 1:05d}")
        top = float(np.percentile(speed, 99.5))
        filled = axes.tricontourf(
            x, y, speed, levels=np.linspace(0.0, top, 25), cmap="viridis", extend="max"
        )
        axes.set_aspect("equal")
        axes.set_xlabel(r"$x$ ($L$)")
        axes.set_ylabel(r"$y$ ($L$)")
        axes.set_title(
            rf"({'cd'[index]}) $|\Phi_{index + 1}|$ at $z=0$, {energy[index]:.1f} %",
            loc="left",
            fontsize=7.5,
        )
        figure.colorbar(
            filled,
            ax=axes,
            orientation="horizontal",
            shrink=0.9,
            pad=0.02,
            ticks=np.linspace(0.0, top, 3),
            format="%.2f",
        )

    output = CASE_DIR / "plot_pod.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
