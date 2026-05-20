#!/usr/bin/env python3
"""plot.py: Visualize thermosyphon Newton base flow results.

Matches AMR_Krylov_V5 paper style: annular domain, temperature field,
inferno colormap, equal-aspect circular geometry.

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

# Thermosyphon geometry: R1=1, R2=2
R_INNER = 1.0
R_OUTER = 2.0


def main():
    nk.configure_style()

    # --- Detect available data ---
    residu_nwt = CASE_DIR / 'residu_newton.dat'
    residu_arn = CASE_DIR / 'residu_arnoldi.dat'
    has_conv = residu_nwt.exists() or residu_arn.exists()

    bf_files = nk.find_fields('BF_tsyphon0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_Ra*_tsyphon0.f00001', CASE_DIR)

    if has_conv and bf_files:
        fig = plt.figure(figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                   nk.COL_WIDTH * 0.5))
        ax_conv = None
    else:
        fig, ax_conv = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.67))
        ax_bf = None

    # --- Panel: Newton convergence ---
    if ax_conv is not None and has_conv:
        nk.plot_newton_convergence(ax_conv, CASE_DIR)
        ax_conv.set_title('Convergence', fontsize=8)
        nk.panel_label(ax_conv, r'$\bf{(a)}$')

    # --- Panel: Base flow temperature ---
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)

        q = fields.get('t')
        if q is not None:
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='inferno', vmin=0., vmax=1.,
                                extend='neither')
            nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[0, 1],
                              tick_labels=['0', '1'])
        ax_bf.axis('equal')
        nk.add_annulus_patches(ax_bf, R_INNER, R_OUTER)
        _setup_annulus_axes(ax_bf)
        lbl = r'$\bf{(b)}$' if ax_conv else r'$\bf{(a)}$'
        nk.panel_label(ax_bf, lbl)

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


def _setup_annulus_axes(ax):
    """Configure axes for circular thermosyphon domain."""
    ticks = [-2, 0, 2]
    labels = ['-2', '0', '2']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
