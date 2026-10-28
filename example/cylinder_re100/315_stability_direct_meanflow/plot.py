#!/usr/bin/env python3
"""Plot the leading direct eigenmode of the Re = 100 cylinder about the mean flow and about the steady flow.

The script reads Spectre_NSd.dat and the real part of the leading mode (dRe1cyl0.f00001) of this
stage (mean flow) and of ../310_stability_direct/direct (steady flow). Both modes are drawn as the
transverse velocity v on the same window, with the same normalisation (maximum of |v| equal to 1,
sign chosen so that v is positive where |v| is maximal) and the same colour scale.

It writes plot_meanflow_stability.png, plot_bf.png (velocity magnitude of the mean flow) and
plot_spectrum.png (full Ritz spectrum of the mean flow, with the steady-flow spectrum of
../310_stability_direct/direct in gray when present). Growth rates are in units of U/D, angular frequencies in
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
STEADY_DIR = CASE_DIR.parent / "310_stability_direct" / "direct"
STEADY_SPECTRUM = STEADY_DIR / "Spectre_NSd.dat"
WINDOW_FILE = CASE_DIR.parent / "140_baseflow_meanflow" / "meanflow_window.dat"
BF_FILE = CASE_DIR / "BF_1cyl0.f00001"
MEAN_FLOW_FILE = CASE_DIR.parent / "140_baseflow_meanflow" / "BFA1cyl0.f00001"
MEAN_SPECTRUM = CASE_DIR / "Spectre_NSd.dat"


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
    """Draw one velocity component of a mode, scaled to a maximum |q| of 1 with a positive sign there, on [-1, 1]."""
    x, y, fields, _time = nk.read_field(str(field_file))
    values = fields[component]
    values = values / values[np.argmax(np.abs(values))]
    filled = ax.tricontourf(
        nk.make_triangulation(x, y),
        values,
        levels=np.linspace(-1.0, 1.0, 21),
        cmap="RdBu_r",
        extend="both",
    )
    nk.add_cylinder_patches(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_ylabel(r"$y$ ($D$)")
    ax.set_title(title, loc="left", fontsize=8)
    return filled


def draw_base_flow(nk, ax, field_file: Path) -> None:
    """Draw the velocity magnitude of a base flow the way stage 310 draws its steady base flow."""
    x, y, fields, _time = nk.read_field(str(field_file))
    triang = nk.make_triangulation(x, y)
    speed = np.sqrt(fields["vx"] ** 2 + fields["vy"] ** 2)
    vmax = 1.5
    cf = nk.tricontourf(
        ax, triang, speed, levels=257, cmap="Blues", vmin=0, vmax=vmax, extend="max"
    )
    nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=[0, vmax],
        tick_labels=["0", f"{vmax}"],
    )
    ax.set_aspect("equal")
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_xlabel(r"$x$ ($D$)", labelpad=-1)
    ax.set_ylabel(r"$y$ ($D$)", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def draw_spectrum(ax, spectrum_file: Path, name: str, color: str, marker: str) -> None:
    """Plot the Ritz values of one spectrum file as sigma against f; converged values are stars with a black edge and the legend gives sigma and f of the leading one."""
    data = np.loadtxt(spectrum_file, ndmin=2)
    sigma = data[:, 0]
    freq = data[:, 1] / (2.0 * np.pi)
    conv = data[:, 2] < 1e-5
    print(f"{name}: {len(sigma)} values, {int(conv.sum())} converged")
    ax.plot(freq[~conv], sigma[~conv], ls="none", marker=marker, ms=3, color=color, label=name)
    label = f"{name}, converged"
    lead = np.flatnonzero(conv & (freq > 0))
    if len(lead):
        k = lead[np.argmax(sigma[lead])]
        label = rf"{name}, converged: $\sigma$ = {sigma[k]:.4f}, $f$ = {freq[k]:.4f}"
    ax.plot(
        freq[conv],
        sigma[conv],
        ls="none",
        marker="*",
        ms=10,
        color=color,
        mec="k",
        mew=0.5,
        label=label,
    )


def plot_base_flow(nk) -> None:
    """Draw the velocity magnitude of the mean flow as plot_bf.png."""
    field_file = BF_FILE if BF_FILE.is_file() else MEAN_FLOW_FILE
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
    draw_base_flow(nk, ax, field_file)
    output = CASE_DIR / "plot_bf.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def plot_spectrum(nk) -> None:
    """Draw the full Ritz spectrum of the mean flow, with the steady-flow spectrum when present, as plot_spectrum.png."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.6, nk.COL_WIDTH * 1.0))
    draw_spectrum(ax, MEAN_SPECTRUM, "mean flow", "C0", "o")
    if STEADY_SPECTRUM.is_file():
        draw_spectrum(ax, STEADY_SPECTRUM, "steady flow", "0.55", "s")
    ax.axhline(0, color="k", lw=0.6, ls="--")
    ax.set_xlabel(r"$f = \omega/2\pi$ ($U/D$)")
    ax.set_ylabel(r"$\sigma$ ($U/D$)")
    ax.set_title("(a) Ritz spectrum, mean flow and steady flow", loc="left", fontsize=8)
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)  # below the axes: inside it hides Ritz values
    output = CASE_DIR / "plot_spectrum.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def main() -> None:
    """Draw the mean-flow eigenvalue with the references and the leading mode."""
    nk = load_nekplot()
    nk.configure_style()

    growth, frequency = leading_eigenvalue(CASE_DIR / "Spectre_NSd.dat")
    print(f"mean flow: growth rate {growth:.6f}, angular frequency {frequency:.6f}")

    fig, (ax_mean, ax_steady) = plt.subplots(
        2, 1, figsize=(7.0, 3.9), sharex=True, layout="constrained"
    )
    mean_file = CASE_DIR / "dRe1cyl0.f00001"
    steady_file = STEADY_DIR / "dRe1cyl0.f00001"
    filled = draw_mode(
        nk,
        ax_mean,
        mean_file,
        "vy",
        rf"(a) Mean flow: $\sigma$ = {growth:.4f}, $f$ = {frequency / (2.0 * np.pi):.4f}",
    )
    if steady_file.is_file() and STEADY_SPECTRUM.is_file():
        growth_steady, frequency_steady = leading_eigenvalue(STEADY_SPECTRUM)
        print(
            f"steady base flow: growth rate {growth_steady:.6f}, angular frequency {frequency_steady:.6f}"
        )
        draw_mode(
            nk,
            ax_steady,
            steady_file,
            "vy",
            rf"(b) Steady flow: $\sigma$ = {growth_steady:.4f}, $f$ = {frequency_steady / (2.0 * np.pi):.4f}",
        )
    ax_mean.set_xlabel("")
    ax_steady.set_xlabel(r"$x$ ($D$)")
    fig.colorbar(
        filled,
        ax=[ax_mean, ax_steady],
        label=r"$v$ of the mode / its maximum",
        ticks=[-1, 0, 1],
        shrink=0.8,
    )

    output = CASE_DIR / "plot_meanflow_stability.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)

    plot_base_flow(nk)
    plot_spectrum(nk)


if __name__ == "__main__":
    main()
