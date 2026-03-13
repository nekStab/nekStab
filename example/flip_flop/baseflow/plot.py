#!/usr/bin/env python3
"""plot.py: Visualize flip-flop Newton base flow (UPO) results.

Matches AMR_Krylov_V5 paper style: PiYG colormap for vx (±1.5),
dual cylinder patches (gap=0.7, radius=0.5), Delaunay tricontourf.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUTPUT = CASE_DIR / 'plot.png'

# Flip-flop dual cylinders
GAP = 0.7
RADIUS = 0.5


def main():
    nk.configure_style()

    # --- Detect available data ---
    has_conv = ((CASE_DIR / 'residu_newton.dat').exists() or
                (CASE_DIR / 'residu_arnoldi.dat').exists())
    bf_files = nk.find_fields('BF_Re60_2cyl0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_*2cyl0.f00001', CASE_DIR)

    if has_conv and bf_files:
        fig = plt.figure(figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1,
                                   figsize=(nk.COL_WIDTH * 0.6,
                                            nk.COL_WIDTH * 0.5))
        ax_conv = None
    else:
        fig, ax_conv = plt.subplots(1, 1,
                                     figsize=(nk.COL_WIDTH * 0.5,
                                              nk.COL_WIDTH * 0.67))
        ax_bf = None

    labels = iter('abcdefgh')

    # --- Panel: Newton convergence ---
    if ax_conv is not None and has_conv:
        nk.plot_newton_convergence(ax_conv, CASE_DIR)
        ax_conv.set_title('Convergence', fontsize=8)
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    # --- Panel: UPO snapshot vx (PiYG, ±1.5) ---
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='PiYG', vmin=-1.5, vmax=1.5,
                                extend='both')
            nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                              width="50%", height="6%", loc=9,
                              ticks=[-1.5, 1.5],
                              tick_labels=['-1.5', '1.5'])
        ax_bf.set_aspect('equal')
        nk.add_dual_cylinder_patches(ax_bf, gap=GAP, radius=RADIUS)
        ax_bf.set_xlim(-5, 30)
        ax_bf.set_ylim(-8, 8)
        ax_bf.set_xlabel(r'$x$', labelpad=-1)
        ax_bf.set_ylabel(r'$y$', labelpad=1)
        ax_bf.spines['right'].set_visible(False)
        ax_bf.spines['top'].set_visible(False)
        nk.panel_label(ax_bf, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
