#!/usr/bin/env python3
"""Plot the DMD of the limit cycle of the cubic cavity at Re = 1950.

The stage reads 100 snapshots of the DNS of ../000_dns, one every 0.5 time units, from t = 1740 to 1789.5, on the
saturated cycle of period T = 10.78909 (../000_dns_period). The files are dmd_spectrum.dat (mode, |mu|, sigma,
omega, St, Re mu, Im mu and the norm of the real part of the mode, per time unit, for the 20 modes of the rank)
and the mode fields dm1cav0.f000mm (real part) and dm2cav0.f000mm (imaginary part) of mode mm of that list.

The list is sorted by |mu|, so the order of the modes follows their damping and not their frequency. A cycle of
period T has the frequencies n 2 pi / T, n = 1, 2, 3, ... For each n up to HARMONICS the script takes the mode with
the angular frequency omega > 0 nearest to n 2 pi / T. It writes dmd_harmonics.dat with one row per n: n, the number
of the mode in the list, omega / (2 pi / T) (n if the mode is on the harmonic), the growth rate sigma of the mode in
1 per time unit, and its norm. It writes plot_dmd.png.

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
HARMONICS = 6  # harmonics of the cycle listed in dmd_harmonics.dat
HARMONICS_SHOWN = 3  # harmonics drawn as fields


def read_spectrum(path: Path) -> np.ndarray:
    """Return the columns mode, |mu|, sigma, omega, St, Re mu, Im mu and the norm (NaN where the file has none)."""
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        values = [float(word) for word in line.split()]
        rows.append(values + [np.nan] * (8 - len(values)))
    return np.array(rows)


def mid_plane_amplitude(
    real_path: Path, imaginary_path: Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and the modulus of the complex velocity mode at the nodes of the plane z = 0.

    The modulus is sqrt(|Re u|^2 + |Im u|^2) at each node. A node on an element edge belongs to several elements;
    its entries are averaged.
    """
    real_data, imaginary_data = readnek(str(real_path)), readnek(str(imaginary_path))
    x_all, y_all, z_all, modulus = [], [], [], []
    for real_element, imaginary_element in zip(
        real_data.elem, imaginary_data.elem, strict=True
    ):
        x_all.append(real_element.pos[0].ravel())
        y_all.append(real_element.pos[1].ravel())
        z_all.append(real_element.pos[2].ravel())
        modulus.append(
            np.sqrt(
                np.sum(real_element.vel**2, axis=0)
                + np.sum(imaginary_element.vel**2, axis=0)
            ).ravel()
        )
    x, y, z, modulus = (np.concatenate(a) for a in (x_all, y_all, z_all, modulus))
    on_plane = np.isclose(z, 0.0, atol=1e-9)
    keys = np.round(np.stack([x[on_plane], y[on_plane]], axis=1), 9)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    average = np.bincount(inverse, weights=modulus[on_plane]) / np.bincount(inverse)
    return unique[:, 0], unique[:, 1], average


def main() -> None:
    """Match the modes to the harmonics, write dmd_harmonics.dat and draw the figure."""
    spectrum = read_spectrum(CASE_DIR / "dmd_spectrum.dat")
    number, magnitude, growth, omega = (
        spectrum[:, 0],
        spectrum[:, 1],
        spectrum[:, 2],
        spectrum[:, 3],
    )
    norm = spectrum[:, 7]
    omega_cycle = 2.0 * np.pi / CYCLE_PERIOD
    ratio = omega / omega_cycle

    rows = []
    for n in range(1, HARMONICS + 1):
        candidates = np.where(omega > 0.0)[0]
        best = candidates[np.argmin(np.abs(ratio[candidates] - n))]
        rows.append([n, number[best], ratio[best], growth[best], norm[best]])
    np.savetxt(
        CASE_DIR / "dmd_harmonics.dat",
        np.array(rows),
        header="n  mode  omega/(2 pi/T)  sigma  norm   (the mode with omega > 0 nearest to n 2 pi/T)",
        fmt=["%d", "%d", "%.8g", "%.8g", "%.8g"],
    )
    print("n  mode  omega/(2 pi/T)  sigma  norm")
    for row in rows:
        print(f"{int(row[0])}  {int(row[1])}  {row[2]:.5f}  {row[3]:.3e}  {row[4]:.4f}")
    print(f"largest |mu| of the {len(number)} modes: {magnitude.max():.6f}")

    figure = plt.figure(figsize=(7.2, 4.9), layout="constrained")
    grid = figure.add_gridspec(2, 6, height_ratios=[0.75, 1.0])
    positive = omega > 0.0

    axes = figure.add_subplot(grid[0, 0:3])
    for n in range(1, HARMONICS + 1):
        axes.axvline(n, c="0.75", lw=0.8, ls=":")
    axes.plot(ratio[positive], norm[positive], c="C0", lw=0.0, marker="o", ms=5)
    axes.set_xlabel(r"$\omega\,/\,(2\pi/T)$")
    axes.set_ylabel("norm of the mode")
    axes.set_title("(a) Modes sit on the harmonics", loc="left", fontsize=7.5)

    axes = figure.add_subplot(grid[0, 3:6])
    axes.axhline(0.0, c="k", lw=0.8, ls="--", label="neutral")
    axes.plot(
        ratio[positive],
        growth[positive],
        c="C3",
        lw=0.0,
        marker="o",
        ms=5,
        label=r"$\sigma$ of a mode",
    )
    axes.set_xlabel(r"$\omega\,/\,(2\pi/T)$")
    axes.set_ylabel(r"$\sigma$ ($U/L$)")
    axes.set_title("(b) No mode grows", loc="left", fontsize=7.5)
    axes.legend(loc="lower left", fontsize=6.5)

    for index in range(HARMONICS_SHOWN):
        mode = int(rows[index][1])
        axes = figure.add_subplot(grid[1, 2 * index : 2 * index + 2])
        x, y, modulus = mid_plane_amplitude(
            CASE_DIR / f"dm1cav0.f{mode:05d}", CASE_DIR / f"dm2cav0.f{mode:05d}"
        )
        top = float(np.percentile(modulus, 99.5))
        filled = axes.tricontourf(
            x,
            y,
            modulus,
            levels=np.linspace(0.0, top, 25),
            cmap="viridis",
            extend="max",
        )
        axes.set_aspect("equal")
        axes.set_xlabel(r"$x$ ($L$)")
        if index == 0:
            axes.set_ylabel(r"$y$ ($L$)")
        axes.set_title(
            rf"({'cde'[index]}) $|\Phi|$, $n={index + 1}$, $z=0$",
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
            format="%.1f",
        )

    output = CASE_DIR / "plot_dmd.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
