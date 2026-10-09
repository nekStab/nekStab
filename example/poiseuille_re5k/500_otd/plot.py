#!/usr/bin/env python3
"""plot.py: OTD of plane Poiseuille flow at Re = 5000, three modes about the parabolic profile.

plot_otd_eigenvalues.png
    (a) real parts of the eigenvalues of Lr against time, with the two exact values as guide lines:
        -nu (pi/2)^2 for the x-independent mode and the Orr-Sommerfeld real part for alpha = 1;
    (b) distance of the first eigenvalue from -nu (pi/2)^2, with the tolerance of ref/reference.json.
plot_otd_modes.png
    the three OTD basis fields at the last output: u of mode 1, v of modes 2 and 3.

USAGE: python plot.py
"""

import sys
from pathlib import Path

sys.path.insert(
    0,
    str(
        next(d for d in Path(__file__).resolve().parents if (d / "nekplot.py").exists())
    ),
)
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nekplot as nk
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
from orr_sommerfeld import spectrum  # noqa: E402

RE = 5000.0
EXACT_1 = -(np.pi**2) / (
    4.0 * RE
)  # u = cos(pi y / 2) decays by viscosity alone: -nu (pi/2)^2
TOL_1 = 1.0e-3 * abs(EXACT_1)  # tol_rel of the first quantity in ref/reference.json
T_FROM = 200.0  # start of the eigenvalue curves


def spread(n_rows, n_marks, shift=0):
    """Indices of n_marks samples spread evenly over n_rows rows, moved by shift rows."""
    return sorted(set(int(i) for i in (np.linspace(0, n_rows - 1, n_marks) + shift).clip(0, n_rows - 1)))


def eigenvalue_figure():
    data = np.loadtxt(CASE_DIR / "otd_eigenvalues.dat")
    # The start from white noise first selects the subspace; its eigenvalues swing by orders of
    # magnitude over the first 200 time units and are left out.
    data = data[data[:, 0] >= T_FROM]
    t = data[:, 0]
    os_real = spectrum(RE, 1.0)[0].real
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)

    ax = axes[0]
    # the thickest line is drawn first, at the back
    styles = [(3.4, "o", 1), (2.2, "s", 2), (1.2, "^", 3)]
    for (lw, marker, col), zorder in zip(styles, (1, 2, 3)):
        if col >= data.shape[1]:
            continue
        ax.plot(
            t,
            data[:, col],
            lw=lw,
            marker=marker,
            ms=3.5,
            markevery=spread(len(t), 14, shift=(col - 1) * len(t) // 42),
            zorder=zorder,
            label=rf"$\Re\,\lambda_{col}$ (end value {data[-1, col]:.4e})",
        )
    ax.axhline(
        EXACT_1, color="k", lw=0.8, ls="--", label=rf"$-\nu(\pi/2)^2$ ({EXACT_1:.4e})"
    )
    ax.axhline(
        os_real,
        color="0.4",
        lw=0.8,
        ls=":",
        label=rf"Orr-Sommerfeld, $\alpha=1$ ({os_real:.4e})",
    )
    ax.set_xlim(T_FROM, t[-1])
    ax.set_ylim(-3.2e-3, 1.7e-3)
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"real part of the eigenvalue of $L_r$ ($1/t_c$)")
    ax.set_title(r"(a) Eigenvalues of $L_r$", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=6.5, loc="upper right")

    ax = axes[1]
    err = np.abs(data[:, 1] - EXACT_1)
    keep = err > 0
    ax.semilogy(
        t[keep],
        err[keep],
        lw=2.0,
        marker="o",
        ms=3.5,
        markevery=spread(int(keep.sum()), 14),
        label=r"$|\Re\,\lambda_1 + \nu(\pi/2)^2|$",
    )
    ax.axhline(
        TOL_1, color="k", lw=0.8, ls="--", label=rf"reference tolerance ({TOL_1:.2e})"
    )
    ax.set_xlim(T_FROM, t[-1])
    ax.set_xlabel(r"$t$ (convective time units)")
    ax.set_ylabel(r"error of the first eigenvalue ($1/t_c$)")
    ax.set_title(r"(b) Convergence of $\lambda_1$", loc="left", fontsize=9)
    ax.legend(frameon=True, fontsize=6.5, loc="upper right")
    out = CASE_DIR / "plot_otd_eigenvalues.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


def last_field(prefix):
    files = nk.find_fields(f"{prefix}poiseuille_OTD0.f*", CASE_DIR)
    return files[-1] if files else None


def mode_figure():
    jobs = [(1, "vx", "u"), (2, "vy", "v"), (3, "vy", "v")]
    panels = []
    for mode, key, name in jobs:
        f = last_field(f"r{mode:02d}")
        if f is not None:
            panels.append((mode, key, name, f))
    if not panels:
        return
    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(6.4, 1.5 * len(panels) + 0.4),
        constrained_layout=True,
        squeeze=False,
    )
    for ax, (mode, key, name, f), letter in zip(axes[:, 0], panels, "abc"):
        x, y, fields, _ = nk.read_field(f)
        q = fields[key]
        bd = float(np.nanmax(np.abs(q)))
        triang = nk.make_triangulation(x, y)
        cf = nk.tricontourf(
            ax, triang, q, levels=65, cmap="RdBu", vmin=-bd, vmax=bd, extend="both"
        )
        nk.inset_colorbar(
            ax,
            cf,
            orientation="horizontal",
            width="30%",
            height="6%",
            loc=9,
            ticks=[-bd, bd],
            tick_labels=[f"{-bd:.2f}", f"{bd:.2f}"],
        )
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$", labelpad=-1)
        ax.set_ylabel(r"$y$", labelpad=1)
        ax.set_title(rf"({letter}) OTD mode {mode}, ${name}$", loc="left", fontsize=9)
    out = CASE_DIR / "plot_otd_modes.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    nk.configure_style()
    if (CASE_DIR / "otd_eigenvalues.dat").exists():
        eigenvalue_figure()
    mode_figure()


if __name__ == "__main__":
    main()
