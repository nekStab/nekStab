#!/usr/bin/env python3
"""plot.py: Visualize Newton UPO results.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUTPUT = CASE_DIR / 'plot.png'


def main():
    nk.configure_style()

    has_conv = ((CASE_DIR / 'residu_newton.dat').exists() or
                (CASE_DIR / 'residu_arnoldi.dat').exists())
    bf_files = nk.find_fields('BF_1cyl0.f*', CASE_DIR)

    if has_conv and bf_files:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45))
        ax_conv = None
    else:
        fig, ax_conv = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5, nk.COL_WIDTH * 0.67))
        ax_bf = None

    labels = iter('abcdefgh')

    if ax_conv is not None and has_conv:
        nk.plot_newton_convergence(ax_conv, CASE_DIR)
        ax_conv.set_title('Convergence', fontsize=8)
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax_bf, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax_bf.set_aspect('equal')
        nk.add_cylinder_patches(ax_bf)
        ax_bf.set_xlim(-2, 20)
        ax_bf.set_ylim(-4, 4)
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
