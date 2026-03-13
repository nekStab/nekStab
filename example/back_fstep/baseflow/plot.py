#!/usr/bin/env python3
"""plot.py: Visualize backward-facing step Newton base flow results.

Matches AMR_Krylov_V5 paper style: cividis colormap for base flow vx,
step geometry patch, Delaunay tricontourf. Double-column width for
the elongated [-5,40]×[-1,2] domain.

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

XLIM = (-5, 40)
YLIM = (-1, 2)
FW = 2 * nk.COL_WIDTH  # double-column width for elongated domain


def main():
    nk.configure_style()

    # --- Detect available data ---
    has_conv = ((CASE_DIR / 'residu_newton.dat').exists() or
                (CASE_DIR / 'residu_arnoldi.dat').exists())
    bf_files = nk.find_fields('BF_bfs0.f00001', CASE_DIR)

    if has_conv and bf_files:
        fig = plt.figure(figsize=(FW, FW * 0.22))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.6, 1.4], wspace=0.3)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(FW, FW * 0.15))
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
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    # --- Panel: Base flow vx (cividis, 0→3.5) ---
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='cividis', vmin=0., vmax=3.5,
                                extend='both')
            nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                              width="25%", height="12%", loc=1,
                              ticks=[0, 3.5],
                              tick_labels=['0', '3.5'])
        nk.add_step_patch(ax_bf, x0=-5, y0=-1, w=5, h=1)
        ax_bf.set_xlim(XLIM)
        ax_bf.set_ylim(YLIM)
        ax_bf.set_aspect('equal')
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
