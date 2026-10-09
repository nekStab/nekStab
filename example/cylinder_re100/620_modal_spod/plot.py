#!/usr/bin/env python3
"""plot.py: cylinder Re = 100 streaming SPOD panels.

  plot_spectrum.png - SPOD eigenvalues against Strouhal number
  plot_mode.png     - vy of the real part of the saved SPOD mode with the
                      largest leading eigenvalue

The code saves the leading modes of the frequency bins 1, 1 + nfreq // 8,
1 + 2 * (nfreq // 8), ... and of the last bin (spod_save_modes in
src/modal_spod.f90). The saved bins are St = 0, 0.125, ..., 1 for nfft = 64
and dt = 0.5. The bin of the spectrum peak is not always one of them.

OUTPUTS: plot_spectrum.png, plot_mode.png
USAGE:   python plot.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

CASE_DIR = Path(__file__).resolve().parent
MODE_GLOB = "sRe1cyl0.f*"


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    root = next(d for d in CASE_DIR.parents if (d / "nekplot.py").exists())
    sys.path.insert(0, str(root))
    import nekplot

    return nekplot


def read_spectrum(path):
    """Return the table of spod_stream_spectrum.dat.

    The columns are St, lambda_1, lambda_2, and so on.
    """
    rows = [
        [float(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    return np.array(rows)


def saved_bins(nfreq):
    """Return the 0-based bins whose modes the code saves, in file order."""
    step = max(1, nfreq // 8)
    return [
        ifreq - 1
        for ifreq in range(1, nfreq + 1)
        if ifreq in (1, nfreq) or (ifreq - 1) % step == 0
    ]


def save_spectrum(nk, table, output):
    """Draw the eigenvalues of every bin; the thickest line is at the back."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.7))
    widths = [1.8, 1.2, 0.7]
    markers = ["o", "s", "^"]
    for index in range(table.shape[1] - 1):
        ax.plot(
            table[:, 0],
            np.maximum(table[:, index + 1], 1e-20),
            color=f"C{index}",
            lw=widths[index % 3],
            marker=markers[index % 3],
            ms=3,
            zorder=3 - index,
            label=rf"$\lambda_{index + 1}$",
        )
    ax.set_yscale("log")
    ax.set_xlabel(r"$St$")
    ax.set_ylabel(r"$\lambda$")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_mode(nk, path, strouhal, eigenvalue, output):
    """Draw vy of a mode file; the legend gives its bin and eigenvalue."""
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    values = fields["vy"]
    bound = round(float(np.nanpercentile(np.abs(values), 99)), 3) or 1.0
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
    cf = nk.tricontourf(
        ax,
        triang,
        values,
        levels=257,
        cmap="RdBu_r",
        vmin=-bound,
        vmax=bound,
        extend="both",
    )
    nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=[-bound, bound],
        tick_labels=[f"{-bound:g}", f"{bound:g}"],
    )
    nk.add_cylinder_patches(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_xlabel(r"$x$", labelpad=-1)
    ax.set_ylabel(r"$y$", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    label = rf"$St$ = {strouhal:g}, $\lambda_1$ = {eigenvalue:.3g}"
    ax.legend(
        handles=[Line2D([], [], ls="none", label=label)],
        fontsize=7,
        loc="lower right",
        handlelength=0,
    )
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    nk = load_nekplot()
    nk.configure_style()

    spectrum_file = CASE_DIR / "spod_stream_spectrum.dat"
    if not spectrum_file.exists():
        return
    table = read_spectrum(spectrum_file)
    save_spectrum(nk, table, CASE_DIR / "plot_spectrum.png")

    modes = sorted(CASE_DIR.glob(MODE_GLOB))
    bins = saved_bins(len(table))
    if len(modes) != len(bins):
        print(
            f"{len(modes)} mode files for {len(bins)} saved bins: mode figure skipped"
        )
        return
    best = max(range(len(bins)), key=lambda n: table[bins[n], 1])
    row = table[bins[best]]
    save_mode(nk, modes[best], row[0], row[1], CASE_DIR / "plot_mode.png")


if __name__ == "__main__":
    main()
