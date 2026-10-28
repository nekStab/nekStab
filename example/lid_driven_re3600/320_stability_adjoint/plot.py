#!/usr/bin/env python3
"""plot.py: lid-driven cavity Re = 3600 adjoint stability panels.

Panel sequence (stability lane: BF -> SPECTRUM -> MODE):
  plot_bf.png       - base flow velocity magnitude (BF_cav0.f00001)
  plot_spectrum.png - converged eigenvalues of Spectre_NSa.dat
                      (growth rate against frequency)
  plot_mode.png     - vx (a) and vy (b), sign fixed, of the leading adjoint mode (aRecav0.f00001)

Geometry: the cavity spans x in [-0.5, 0.5] and y in [0, userParam10]; the
limits follow the coordinates of the field file.

OUTPUTS: plot_bf.png, plot_spectrum.png, plot_mode.png
USAGE:   python plot.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
KIND = "adjoint"


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


def render_field(nk, ax, path, what, cmap, symmetric, fix_sign=False, label=None):
    """Draw the velocity magnitude (what = 'umag') or one field of a file.

    fix_sign: an eigenvector has no sign. Multiply the field by +1 or -1 so
    that vy is positive where |vy| is largest. The direct and the adjoint
    mode then share the sign convention and their colours can be compared.
    """
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    if what == "umag":
        values = np.hypot(fields["vx"], fields["vy"])
    else:
        values = fields[what]
    if fix_sign:
        vy = fields["vy"]
        values = values * np.sign(vy[np.argmax(np.abs(vy))])
    bound = round(float(np.nanpercentile(np.abs(values), 99)), 3) or 1.0
    if symmetric:
        vmin, vmax, extend = -bound, bound, "both"
    else:
        vmin, vmax, extend = 0.0, bound, "max"
    cf = nk.tricontourf(
        ax, triang, values, levels=257, cmap=cmap, vmin=vmin, vmax=vmax, extend=extend
    )
    cbar = nk.inset_colorbar(
        ax,
        cf,
        orientation="horizontal",
        width="50%",
        height="5%",
        loc=9,
        ticks=[vmin, vmax],
        tick_labels=[f"{vmin:g}", f"{vmax:g}"],
    )
    if label:
        cbar.set_label(label, fontsize=6, labelpad=1)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$", labelpad=-1)
    ax.set_ylabel(r"$y$", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def save_field(nk, path, what, cmap, symmetric, output):
    """Render one field file to a square figure."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.8, nk.COL_WIDTH * 0.8))
    render_field(nk, ax, path, what, cmap, symmetric)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_mode(nk, path, output):
    """Draw vx (a) and vy (b) of the leading mode with the sign fixed.

    The direct and the adjoint mode are 90 to 95 % correlated but differ by
    0.3 in vx and 0.44 in vy (relative L2 norm, norm-matched). Showing both
    components with one sign convention makes that difference visible.
    """
    fig, axes = plt.subplots(1, 2, figsize=(nk.COL_WIDTH * 1.6, nk.COL_WIDTH * 0.9))
    for ax, what, tag in zip(axes, ("vx", "vy"), ("a", "b")):
        render_field(nk, ax, path, what, "RdBu_r", True, fix_sign=True,
                     label=rf"${what[0]}_{what[1]}$ ({KIND} mode, $\max|v_y|>0$)")
        ax.set_title(rf"({tag}) ${what[0]}_{what[1]}$ of the {KIND} mode", loc="left", fontsize=7)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_spectrum(nk, path, output):
    """Draw the converged eigenvalues; the stable half-plane is shaded."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
    nk.plot_spectrum_NS_paper(ax, path, label=r"$Re = 3600$")
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

    base_flow = first_match("BF_cav0.f00001", "BF_*cav*0.f00001")
    if base_flow:
        save_field(nk, base_flow, "umag", "Blues", False, CASE_DIR / "plot_bf.png")

    spectrum = first_match("Spectre_NSa.dat", "Spectre_NSa_conv.dat")
    if spectrum:
        save_spectrum(nk, spectrum, CASE_DIR / "plot_spectrum.png")

    mode = first_match("aRecav0.f00001")
    if mode:
        save_mode(nk, mode, CASE_DIR / "plot_mode.png")


if __name__ == "__main__":
    main()
