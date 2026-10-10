#!/usr/bin/env python3
"""Plot the leading eigenvalue and eigenmode of the Re = 100 cylinder mean flow.

The script reads Spectre_NSd.dat and the real part of the leading mode (dRe1cyl0.f00001) of this
stage. When the neighbouring stages are present, it adds the leading eigenvalue of the steady
base flow (../310_stability_direct/direct) and the shedding frequency of the DNS
(../140_baseflow_meanflow/meanflow_window.dat) as references.

It writes plot_meanflow_stability.png. Growth rates are in units of U/D, angular frequencies in
units of U/D, lengths in units of the cylinder diameter.

Usage: python plot.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.contour import QuadContourSet

CASE_DIR = Path(__file__).resolve().parent
STEADY_SPECTRUM = (
    CASE_DIR.parent / "310_stability_direct" / "direct" / "Spectre_NSd.dat"
)
WINDOW_FILE = CASE_DIR.parent / "140_baseflow_meanflow" / "meanflow_window.dat"


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def leading_eigenvalue(spectrum_file: Path) -> tuple[float, float]:
    """Return growth rate and angular frequency of the first row of a spectrum file."""
    row = np.loadtxt(spectrum_file, ndmin=2)[0]
    return float(row[0]), float(row[1])


def draw_mode(nk, ax, field_file: Path, component: str, title: str) -> QuadContourSet:
    """Draw one velocity component of a mode on a symmetric colour scale and return the contour set."""
    x, y, fields, _time = nk.read_field(str(field_file))
    values = fields[component]
    bound = float(np.nanpercentile(np.abs(values), 99))
    filled = ax.tricontourf(
        nk.make_triangulation(x, y),
        values / bound,
        levels=np.linspace(-1.0, 1.0, 21),
        cmap="RdBu_r",
        extend="both",
    )
    nk.add_cylinder_patches(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, 12.0)
    ax.set_ylim(-3.0, 3.0)
    ax.set_ylabel(r"$y$ ($D$)")
    ax.set_title(title, loc="left", fontsize=8)
    return filled


def main() -> None:
    """Draw the mean-flow eigenvalue with the references and the leading mode."""
    nk = load_nekplot()
    nk.configure_style()

    growth, frequency = leading_eigenvalue(CASE_DIR / "Spectre_NSd.dat")
    print(f"mean flow: growth rate {growth:.6f}, angular frequency {frequency:.6f}")

    fig = plt.figure(figsize=(7.0, 3.9), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.45])
    ax_eig = fig.add_subplot(grid[:, 0])
    ax_vx = fig.add_subplot(grid[0, 1])
    ax_vy = fig.add_subplot(grid[1, 1])

    if STEADY_SPECTRUM.is_file():
        growth_steady, frequency_steady = leading_eigenvalue(STEADY_SPECTRUM)
        print(
            f"steady base flow: growth rate {growth_steady:.6f}, angular frequency {frequency_steady:.6f}"
        )
        ax_eig.plot(
            [growth_steady],
            [frequency_steady],
            "s",
            c="0.55",
            ms=7,
            label=rf"steady flow: $\sigma$ = {growth_steady:.4f}, $\omega$ = {frequency_steady:.4f}",
        )
    ax_eig.plot(
        [growth],
        [frequency],
        "o",
        c="C0",
        ms=4,
        label=rf"mean flow: $\sigma$ = {growth:.4f}, $\omega$ = {frequency:.4f}",
    )
    ax_eig.axvline(0.0, c="k", lw=0.8, ls=":", label=r"$\sigma = 0$")
    if WINDOW_FILE.is_file():
        period = float(np.loadtxt(WINDOW_FILE, ndmin=2)[0, 3])
        omega_dns = 2.0 * np.pi / period
        print(f"DNS: period {period:.5f}, angular frequency {omega_dns:.6f}")
        ax_eig.axhline(
            omega_dns,
            c="C3",
            lw=0.9,
            ls="--",
            label=rf"DNS, $2\pi/T$ = {omega_dns:.4f}",
        )
    ax_eig.set_xlabel(r"growth rate $\sigma$ ($U/D$)")
    ax_eig.set_ylabel(r"angular frequency $\omega$ ($U/D$)")
    ax_eig.set_title("(a) Leading eigenvalue", loc="left", fontsize=8)
    ax_eig.set_xlim(-0.01, 0.14)
    ax_eig.set_ylim(0.70, 1.10)
    ax_eig.legend(loc="center right", fontsize=5.5)

    mode_file = CASE_DIR / "dRe1cyl0.f00001"
    filled = draw_mode(nk, ax_vx, mode_file, "vx", r"(b) Eigenmode, $u$ (real part)")
    draw_mode(nk, ax_vy, mode_file, "vy", r"(c) Eigenmode, $v$ (real part)")
    ax_vx.set_xticklabels([])
    ax_vy.set_xlabel(r"$x$ ($D$)")
    fig.colorbar(
        filled,
        ax=[ax_vx, ax_vy],
        label="mode / its 99th percentile",
        ticks=[-1, 0, 1],
        shrink=0.8,
    )

    output = CASE_DIR / "plot_meanflow_stability.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
