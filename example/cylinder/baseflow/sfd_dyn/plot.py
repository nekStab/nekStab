#!/usr/bin/env python3
"""plot.py: Visualize the SFD dynamic-tolerance base-flow testcase.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import matplotlib.pyplot as plt
import numpy as np

import nekplot as nk

CASE_DIR = Path(__file__).resolve().parent
OUTPUT = CASE_DIR / "plot.png"


def find_result_field():
    checkpoint_files = nk.find_fields("1cyl0.f*", CASE_DIR)
    if checkpoint_files:
        return checkpoint_files[-1]

    result_files = nk.find_fields("BF_1cyl0.f*", CASE_DIR)
    if result_files:
        return result_files[-1]

    fallback_files = nk.find_fields("BF*1cyl0.f*", CASE_DIR)
    return fallback_files[-1] if fallback_files else None


def main():
    nk.configure_style()

    residu = CASE_DIR / "residu.dat"
    baseline_residu = CASE_DIR.parent / "sfd" / "residu.dat"
    dyn_tol = CASE_DIR / "dyn_tol.dat"
    has_conv = residu.exists()
    has_dyn_tol = dyn_tol.exists()
    bf_file = find_result_field()

    if has_conv and has_dyn_tol and bf_file is not None:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 3, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 0.8, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_tol = fig.add_subplot(gs[1])
        ax_bf = fig.add_subplot(gs[2])
    elif has_conv and bf_file is not None:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_tol = None
        ax_bf = fig.add_subplot(gs[1])
    elif has_conv and has_dyn_tol:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 0.8], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_tol = fig.add_subplot(gs[1])
        ax_bf = None
    elif bf_file is not None:
        fig, ax_bf = plt.subplots(
            1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45)
        )
        ax_conv = None
        ax_tol = None
    else:
        fig, ax_conv = plt.subplots(
            1, 1, figsize=(nk.COL_WIDTH * 0.5, nk.COL_WIDTH * 0.67)
        )
        ax_tol = None
        ax_bf = None

    labels = iter("abcdefgh")

    if ax_conv is not None and has_conv:
        if baseline_residu.exists():
            nk.plot_residuals(ax_conv, baseline_residu, label="SFD fixed",
                              color="0.45")
        nk.plot_residuals(ax_conv, residu, label="SFD dyn", color="m")
        ax_conv.axhline(1e-9, color="k", lw=0.4, ls=":")
        ax_conv.set_title("SFD With Dynamic Tolerances", fontsize=8)
        nk.panel_label(ax_conv, rf"$\bf{{({next(labels)})}}$")

    if ax_tol is not None and has_dyn_tol:
        residu_data = np.genfromtxt(str(residu))
        dyn_data = np.genfromtxt(str(dyn_tol), comments="#")
        if residu_data.ndim == 1:
            residu_data = residu_data.reshape(1, -1)
        if dyn_data.ndim == 1:
            dyn_data = dyn_data.reshape(1, -1)

        ax_tol.semilogy(residu_data[:, 0], residu_data[:, 3], color="0.35",
                        lw=0.7, label="current tol")
        ax_tol.semilogy(dyn_data[:, 0], dyn_data[:, 3], "o-", ms=2.0,
                        lw=0.6, color="C1", label="requested")
        ax_tol.semilogy(dyn_data[:, 0], dyn_data[:, 4], "o-", ms=2.0,
                        lw=0.6, color="C2", label="used")
        if np.any(dyn_data[:, 5] > 0.0):
            ax_tol.semilogy(dyn_data[:, 0], dyn_data[:, 5], "--",
                            lw=0.6, color="C3", label="cap")
        ax_tol.set_xlabel(r"$t$")
        ax_tol.set_ylabel("solver tol")
        ax_tol.set_title("Tolerance Scheduler", fontsize=8)
        ax_tol.legend(fontsize=5.5, loc="best")
        nk.panel_label(ax_tol, rf"$\bf{{({next(labels)})}}$")

    if ax_bf is not None and bf_file is not None:
        x, y, fields, time = nk.read_field(bf_file)
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields["vx"] ** 2 + fields["vy"] ** 2)
        cf = nk.tricontourf(
            ax_bf,
            triang,
            umag,
            levels=257,
            cmap="Blues",
            vmin=0,
            vmax=1.5,
            extend="max",
        )
        nk.inset_colorbar(
            ax_bf,
            cf,
            orientation="horizontal",
            width="50%",
            height="5%",
            loc=9,
            ticks=[0, 1.5],
            tick_labels=["0", "1.5"],
        )
        ax_bf.set_aspect("equal")
        nk.add_cylinder_patches(ax_bf)
        ax_bf.set_xlim(-2, 20)
        ax_bf.set_ylim(-4, 4)
        ax_bf.set_xlabel(r"$x$", labelpad=-1)
        ax_bf.set_ylabel(r"$y$", labelpad=1)
        ax_bf.spines["right"].set_visible(False)
        ax_bf.spines["top"].set_visible(False)
        nk.panel_label(ax_bf, rf"$\bf{{({next(labels)})}}$")

    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight")
    print(f"Saved {OUTPUT}")
    plt.close()


if __name__ == "__main__":
    main()
