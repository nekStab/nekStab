#!/usr/bin/env python3
"""plot.py: NACA 0012 Re = 2000 adjoint stability panels.

Panel sequence (stability lane: BF -> SPECTRUM -> MODE):
  plot_bf.png       - base flow velocity magnitude (BF_naca00120.f00001)
  plot_spectrum.png - converged eigenvalues of Spectre_NSa.dat (growth rate
                      against frequency); the pair sigma = -0.557 +- 9.43 i is marked
  plot_mode.png     - vy of the leading adjoint mode (aRenaca00120.f00001) in (a) and of the
                      pair sigma = -0.557 +- 9.43 i in (b). The mode file number is the row
                      number of Spectre_NSa_conv.dat.
  ref/mode.gif      - vy of the leading mode over one period,
                      Re[(qr + i qi) exp(i omega t)] = qr cos(omega t) - qi sin(omega t)

Geometry: NACA 0012 airfoil, xlim = (-1, 5), ylim = (-2, 2).

OUTPUTS: plot_bf.png, plot_spectrum.png, plot_mode.png, ref/mode.gif
USAGE:   python plot.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

CASE_DIR = Path(__file__).resolve().parent
XLIM = (-1, 5)
YLIM = (-2, 2)
N_FRAMES = 24


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


def draw_field(nk, ax, triang, values, bound, cmap, levels=257):
    """Draw values on a fixed colour scale [-bound, bound] (or [0, 1.5] when bound is None)."""
    if bound is None:
        vmin, vmax, extend = 0.0, 1.5, "max"
    else:
        vmin, vmax, extend = -bound, bound, "both"
    cf = nk.tricontourf(
        ax, triang, values, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, extend=extend
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


def render_field(nk, ax, path, what, cmap, symmetric):
    """Draw the velocity magnitude (what = 'umag') or one field of a file."""
    x, y, fields, _time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    if what == "umag":
        values = np.hypot(fields["vx"], fields["vy"])
    else:
        values = fields[what]
    bound = float(np.nanpercentile(np.abs(values), 99)) if symmetric else None
    draw_field(nk, ax, triang, values, bound, cmap)


def save_field(nk, path, what, cmap, symmetric, output):
    """Render one field file to a wide figure."""
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.8, nk.COL_WIDTH * 0.7))
    render_field(nk, ax, path, what, cmap, symmetric)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def pair_row(conv):
    """Index of the conv row with omega > 0 closest to the pair sigma = -0.557 +- 9.43 i."""
    data = np.genfromtxt(conv)
    return int(np.argmin(np.abs(data[:, 1] - 9.433) + (data[:, 1] < 0) * 1e3)), data


def save_spectrum(nk, path, conv, output):
    """Draw every converged eigenvalue (flag 1); the stable half-plane is shaded."""
    sigma, omega, _res, flag = nk._read_spectrum(path)
    f = omega / (2.0 * np.pi)
    ok = flag > 0
    k, cdata = pair_row(conv)
    lead = ok & (np.abs(np.abs(omega) - cdata[0, 1]) < 1e-3)
    pair = ok & (np.abs(np.abs(omega) - cdata[k, 1]) < 1e-3)
    rest = ok & ~lead & ~pair
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.8))
    ax.axhspan(-1.6, 0, facecolor="gray", alpha=0.3, zorder=-1)
    ax.axhline(0, lw=0.5, c="gray", ls="dotted")
    ax.axvline(0, lw=0.5, c="gray", ls="dotted")
    for sel, label, mk, col, ms in (
        (rest, "other converged", "d", "k", 3.5),
        (lead, rf"$\sigma={cdata[0, 0]:.4f}\pm{cdata[0, 1]:.3f}i$", "o", "tab:red", 5),
        (pair, rf"$\sigma={cdata[k, 0]:.3f}\pm{cdata[k, 1]:.2f}i$", "s", "tab:blue", 5),
    ):
        ax.plot(f[sel], sigma[sel], ls="none", marker=mk, mfc=col, mec="k", mew=0.3,
                ms=ms, label=label, zorder=3)
    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$\sigma$")
    ax.set_ylim(-1.6, 0.1)
    ax.legend(fontsize=6, loc="lower center", framealpha=0.9)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_modes(nk, conv, output):
    """(a) leading mode, (b) mode of the pair sigma = -0.557 +- 9.43 i; vy of the real part."""
    k, cdata = pair_row(conv)
    rows = ((0, "leading mode"), (k, "stable pair"))
    fig, axes = plt.subplots(2, 1, figsize=(nk.COL_WIDTH * 1.8, nk.COL_WIDTH * 1.5))
    for ax, (row, name), letter in zip(axes, rows, "ab"):
        path = CASE_DIR / f"aRenaca00120.f{row + 1:05d}"
        render_field(nk, ax, path, "vy", "RdBu_r", True)
        ax.set_title(
            rf"({letter}) {name}, $\sigma={cdata[row, 0]:.4f}$, $\omega={cdata[row, 1]:.3f}$",
            loc="left", fontsize=8,
        )
    fig.subplots_adjust(hspace=0.45)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def save_gif(nk, conv, output):
    """vy of the leading mode over one period: qr cos(omega t) - qi sin(omega t)."""
    omega = abs(np.genfromtxt(conv)[0, 1])
    x, y, fr, _ = nk.read_field(str(CASE_DIR / "aRenaca00120.f00001"))
    _, _, fi, _ = nk.read_field(str(CASE_DIR / "aImnaca00120.f00001"))
    qr, qi = np.nan_to_num(fr["vy"]), np.nan_to_num(fi["vy"])
    triang = nk.make_triangulation(x, y)
    # fixed symmetric scale: the peak of the amplitude |qr + i qi| over the period
    bound = float(np.nanpercentile(np.hypot(qr, qi), 99))
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.8, nk.COL_WIDTH * 0.7))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.17, top=0.9)

    def update(i):
        ax.clear()
        phase = 2.0 * np.pi * i / N_FRAMES
        draw_field(nk, ax, triang, qr * np.cos(phase) - qi * np.sin(phase), bound,
                   "RdBu_r", levels=129)
        ax.set_title(rf"$vy$, $t={i / N_FRAMES:.2f}\,T$, $T=2\pi/\omega={2 * np.pi / omega:.3f}$",
                     loc="left", fontsize=8)
        return []

    anim = FuncAnimation(fig, update, frames=N_FRAMES, blit=False)
    anim.save(output, writer=PillowWriter(fps=12), dpi=150)
    plt.close(fig)
    print(f"Saved {output}")


def main():
    nk = load_nekplot()
    nk.configure_style()

    base_flow = first_match("BF_naca00120.f00001", "BF_*naca*0.f00001")
    if base_flow:
        save_field(nk, base_flow, "umag", "Blues", False, CASE_DIR / "plot_bf.png")

    spectrum = first_match("Spectre_NSa.dat", "Spectre_NSa_conv.dat")
    conv = first_match("Spectre_NSa_conv.dat")
    if spectrum and conv:
        save_spectrum(nk, spectrum, conv, CASE_DIR / "plot_spectrum.png")

    if conv and first_match("aRenaca00120.f0000*"):
        save_modes(nk, conv, CASE_DIR / "plot_mode.png")
        save_gif(nk, conv, CASE_DIR / "ref" / "mode.gif")


if __name__ == "__main__":
    main()
