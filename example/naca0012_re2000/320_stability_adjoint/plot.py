#!/usr/bin/env python3
"""plot.py: NACA 0012 Re = 2000 adjoint stability panels.

Panel sequence (stability lane: BF -> SPECTRUM -> MODE):
  plot_bf.png       - base flow velocity magnitude (BF_naca00120.f00001)
  plot_spectrum.png - converged eigenvalues of Spectre_NSa.dat
                      (growth rate against frequency)
  plot_mode.png     - vy of the leading adjoint mode (aRenaca00120.f00001)

Geometry: NACA 0012 airfoil, xlim = (-1, 5), ylim = (-2, 2).

OUTPUTS: plot_bf.png, plot_spectrum.png, plot_mode.png
USAGE:   python plot.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
XLIM = (-1, 5)
YLIM = (-2, 2)


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    root = next(d for d in CASE_DIR.parents if (d / "nekplot.py").exists())
    sys.path.insert(0, str(root))
    import nekplot

    return nekplot


def first_match(*patterns):
    """Return the first file of the case folder that matches one of the patterns."""
    for pattern in patterns:
        matches = sorted(CASE_DIR.glob(pattern))
        if matches:
            return matches[0]
    return None


def render_field(nk, ax, path, what, cmap, symmetric):
    """Draw the velocity magnitude (what = 'umag') or one field of a file."""
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    if what == "umag":
        values = np.hypot(fields["vx"], fields["vy"])
    else:
        values = fields[what]
    if symmetric:
        bound = float(np.nanpercentile(np.abs(values), 99))
        vmin, vmax, extend = -bound, bound, "both"
    else:
        vmin, vmax, extend = 0.0, 1.5, "max"
    cf = nk.tricontourf(
        ax, triang, values, levels=257, cmap=cmap, vmin=vmin, vmax=vmax, extend=extend
    )
    ticks = [round(vmin, 3), round(vmax, 3)]
    nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=ticks,
        tick_labels=[f"{tick:g}" for tick in ticks],
    )
    nk.add_naca0012_patch(ax)
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$", labelpad=-1)
    ax.set_ylabel(r"$y$", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def save_field(nk, path, what, cmap, symmetric, output):
    """Render one field file to a wide figure."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.8, nk.COL_WIDTH * 0.7))
    render_field(nk, ax, path, what, cmap, symmetric)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_spectrum(nk, path, output):
    """Draw the converged eigenvalues; the stable half-plane is shaded."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
    nk.plot_spectrum_NS_paper(ax, path, label=r"$Re = 2000$")
    ax.axhspan(min(ax.get_ylim()[0], -0.2), 0, facecolor="gray", alpha=0.3, zorder=-1)
    ax.legend(fontsize=7, loc="lower center")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    nk = load_nekplot()
    nk.configure_style()

    base_flow = first_match("BF_naca00120.f00001", "BF_*naca*0.f00001")
    if base_flow:
        save_field(nk, base_flow, "umag", "Blues", False, CASE_DIR / "plot_bf.png")

    spectrum = first_match("Spectre_NSa.dat", "Spectre_NSa_conv.dat")
    if spectrum:
        save_spectrum(nk, spectrum, CASE_DIR / "plot_spectrum.png")

    mode = first_match("aRenaca00120.f00001")
    if mode:
        save_field(nk, mode, "vy", "RdBu_r", True, CASE_DIR / "plot_mode.png")


if __name__ == "__main__":
    main()
