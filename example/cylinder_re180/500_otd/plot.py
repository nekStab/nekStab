#!/usr/bin/env python3
"""plot.py: OTD about the Newton orbit of the cylinder at Re = 180 (userParam01 = 5.1, lpert = 1).

Phase 1 (files p1_otd_*.dat; plain otd_*.dat when the run has no phase 1): white-noise start.
Phase 2 (otd_*.dat, r01 fields): one period from the basis at 40 T, r011cyl0.f00001..8 every
T/8, and the orbit 1cyl0.f00001..8 at the same times.

plot_otd.png          (a) exponent of the OTD mode for each period, guide 0 = phase mode;
                      (b) its absolute value, guide = tolerance of ref/reference.json
plot_otd_history.png  (a) instantaneous growth rate (eigenvalue of the symmetric part of the
                      reduced operator) over the periods; (b) absolute per-period exponent
plot_otd_mode.png     vy of the OTD mode at phases 0, T/4, T/2, 3T/4
plot_otd_vs_floquet.png  OTD mode at phase 0 against the phase mode dU/dt of the orbit
plot_otd_mode.gif     vy of the OTD mode over the 8 snapshots, fixed colour scale
USAGE:  python plot.py     (a figure whose data files are missing is skipped)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nekplot as nk
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

CASE_DIR = Path(__file__).resolve().parent
TOLERANCE = 1.0e-4  # tol_abs of otd_period_exponent_final in ref/reference.json
T = 5.1896284145  # period of the orbit
XLIM, YLIM = (-2.0, 20.0), (-4.0, 4.0)  # wake window of 311 mode.png
NF = 8  # snapshots per period


def p1(name):
    """Phase-1 file when it exists, else the plain name."""
    f = CASE_DIR / f"p1_{name}"
    return f if f.is_file() else CASE_DIR / name


def exponent_plot() -> None:
    data = np.loadtxt(p1("otd_ftle_period.dat"))
    t, sigma = data[:, 0], data[:, 1]
    label = f"OTD mode, {len(t)} periods"

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)

    ax = axes[0]
    ax.plot(t, sigma, lw=1.6, marker="o", ms=3.5, label=label)
    ax.axhline(0.0, color="k", lw=0.8, ls="--", label="phase mode of the orbit (0)")
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"exponent per period ($1/t_c$)")
    ax.set_title("(a) Exponent of the OTD mode", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    ax = axes[1]
    ax.semilogy(t, np.abs(sigma), lw=1.6, marker="o", ms=3.5, label=label)
    ax.axhline(
        TOLERANCE,
        color="k",
        lw=0.8,
        ls="--",
        label=f"tolerance of the reference ({TOLERANCE:.0e})",
    )
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"$|$exponent per period$|$ ($1/t_c$)")
    ax.set_title("(b) Absolute value", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8)

    out = CASE_DIR / "plot_otd.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


def history_plot() -> None:
    g = np.loadtxt(p1("otd_growth_rates.dat"))
    f = np.loadtxt(p1("otd_ftle_period.dat"))
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)

    ax = axes[0]
    ax.plot(g[:, 0], g[:, 1], lw=1.0, marker="o", ms=3.0, markevery=1500,
            label=r"$\sigma_1$ of $L_s$")
    ax.axhline(0.0, color="k", lw=0.8, ls="--", label="phase mode of the orbit (0)")
    for k in range(1, int(round(g[-1, 0] / T)) + 1):
        ax.axvline(k * T, color="0.8", lw=0.4, zorder=0)
    ax.set_yscale("symlog", linthresh=1.0e-2)
    ax.set_ylim(top=3.0)  # room for the legend above the data
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"growth rate ($1/t_c$)")
    ax.set_title("(a) Instantaneous growth rate", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8, loc="upper left")

    ax = axes[1]
    ax.semilogy(f[:, 0], np.abs(f[:, 1]), lw=1.6, marker="o", ms=3.5,
                label=f"OTD mode, {len(f)} periods")
    ax.axhline(TOLERANCE, color="k", lw=0.8, ls="--",
               label=f"tolerance of the reference ({TOLERANCE:.0e})")
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"$|$exponent per period$|$ ($1/t_c$)")
    ax.set_title("(b) Exponent per period", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=8, loc="lower left")

    out = CASE_DIR / "plot_otd_history.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


def frames(prefix):
    """x, y and {vx, vy} of the 8 snapshots; snapshot k is at t = k T/8, snapshot 8 at phase 0."""
    files = [CASE_DIR / f"{prefix}1cyl0.f{k:05d}" for k in range(1, NF + 1)]
    if not all(p.is_file() for p in files):
        return None
    out = [nk.read_field(p) for p in files]
    x, y = out[0][0], out[0][1]
    for o, p in zip(out, files):
        assert np.array_equal(o[0], x) and np.array_equal(o[1], y), f"mesh of {p.name} differs"
    return x, y, np.array([o[2]["vx"] for o in out]), np.array([o[2]["vy"] for o in out]), np.array([o[3] for o in out])


def window_p99(x, y, q):
    m = (x >= XLIM[0]) & (x <= XLIM[1]) & (y >= YLIM[0]) & (y <= YLIM[1])
    return float(np.nanpercentile(np.abs(q[..., m]), 99))


def draw(ax, triang, q, bd, levels=257):
    cf = nk.tricontourf(ax, triang, q, levels=levels, cmap="RdBu_r", vmin=-bd, vmax=bd, extend="both")
    nk.add_cylinder_patches(ax)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_aspect("equal")
    ax.set_ylabel(r"$y$", labelpad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    return cf


def stack(items, bd, name):
    """Panels (a), (b), ... of fields on one mesh, one shared colour scale."""
    x, y, triang = items[0]
    n = len(items) - 1
    fig, axes = plt.subplots(n, 1, figsize=(6.2, 1.55 * n + 0.5), constrained_layout=True, sharex=True)
    for ax, (q, sub), letter in zip(axes, items[1:], "abcdefgh"):
        cf = draw(ax, triang, q, bd)
        ax.set_title(f"({letter}) {sub}", loc="left", fontsize=9)
    axes[-1].set_xlabel(r"$x$", labelpad=-1)
    cb = fig.colorbar(cf, ax=axes, location="right", shrink=0.8, aspect=30, ticks=[-bd, 0, bd])
    cb.ax.tick_params(labelsize=6)
    cb.ax.set_yticklabels([f"{-bd:.2g}", "0", f"{bd:.2g}"])
    out = CASE_DIR / name
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


def mode_plots() -> None:
    r = frames("r01")
    if r is None:
        print("INFO: no r01 snapshots, mode figures skipped")
        return
    x, y, vx, vy, t = r
    triang = nk.make_triangulation(x, y)
    bd = window_p99(x, y, vy)
    # snapshot k is at t = k T/8: phases 0, T/4, T/2, 3T/4 are snapshots 8, 2, 4, 6
    sel = [(NF, "0"), (2, "T/4"), (4, "T/2"), (6, "3T/4")]
    stack([(x, y, triang)] + [(vy[k - 1], rf"OTD mode $v_y$, phase $t={p}$") for k, p in sel],
          bd, "plot_otd_mode.png")

    fig, ax = plt.subplots(figsize=(6.0, 2.4), constrained_layout=True)
    cf = draw(ax, triang, vy[0], bd, levels=65)
    ax.set_xlabel(r"$x$", labelpad=-1)
    cb = fig.colorbar(cf, ax=ax, location="right", shrink=0.8, aspect=20, ticks=[-bd, 0, bd])
    cb.ax.tick_params(labelsize=6)

    def upd(i):
        for c in ax.collections:
            c.remove()
        draw(ax, triang, vy[i], bd, levels=65)
        ax.set_title(rf"(a) OTD mode $v_y$, $t={(i + 1) * T / NF:.3f}$ ({(i + 1)}T/8)", loc="left", fontsize=9)
        return []

    anim = FuncAnimation(fig, upd, frames=NF)
    out = CASE_DIR / "plot_otd_mode.gif"
    anim.save(out, writer=PillowWriter(fps=3), dpi=90)
    plt.close(fig)
    print(f"Saved {out}")


def phase_mode_plot() -> None:
    r, o = frames("r01"), frames("")
    if r is None or o is None:
        print("INFO: no OTD or orbit snapshots, comparison skipped")
        return
    x, y, vx, vy, t = r
    assert np.array_equal(x, o[0]) and np.array_equal(y, o[1]), "mesh of orbit and OTD fields differs"
    assert np.allclose(t, o[4]), "snapshot times of orbit and OTD fields differ"
    # dU/dt at phase 0 (snapshot 8) by the Fourier derivative of the 8 samples of one period,
    # kept up to harmonic kc. Harmonics 3 and 4 carry about the same energy, so 8 samples
    # resolve dU/dt only coarsely: the difference is printed for kc = 1, 2, 3 to show that.
    unit = lambda v: v / np.linalg.norm(v)
    otd = unit(np.concatenate([vx[NF - 1], vy[NF - 1]]))
    F = [np.fft.fft(c, axis=0) for c in (o[2], o[3])]
    h = np.fft.fftfreq(NF, d=1.0 / NF)
    res = {}
    for kc in (1, 2, 3):
        k = np.where(np.abs(h) <= kc, h * 2 * np.pi / T, 0.0)
        ph = unit(np.concatenate([np.fft.ifft(1j * k[:, None] * f, axis=0).real[NF - 1] for f in F]))
        s = np.sign(otd @ ph)
        res[kc] = (s * ph, np.linalg.norm(otd - s * ph))
    print("OTD mode vs phase mode dU/dt, relative L2 difference for harmonics up to 1, 2, 3: "
          + ", ".join(f"{res[kc][1]:.3f}" for kc in (1, 2, 3)))
    ph, rel = res[3]
    diff = otd - ph
    n = len(x)
    triang = nk.make_triangulation(x, y)
    a, b, c = otd[n:], ph[n:], diff[n:]
    bd = window_p99(x, y, np.stack([a, b]))
    stack([(x, y, triang), (a, r"OTD mode $v_y$, phase 0"), (b, r"phase mode $\partial_t v_y$ of the orbit"),
           (c, "difference of the two")], bd, "plot_otd_vs_floquet.png")


def main() -> None:
    nk.configure_style()
    plt.rcParams["savefig.dpi"] = 400
    exponent_plot()
    history_plot()
    mode_plots()
    phase_mode_plot()


if __name__ == "__main__":
    main()
