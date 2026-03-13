#!/usr/bin/env python3
"""plot.py: Visualize Newton (dynamic + temperature) base flow results.

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
    bf_files = nk.find_fields('BF_*0.f*', CASE_DIR)

    ncols = int(has_conv) + 2  # velocity + temperature panels
    if not bf_files:
        ncols = max(int(has_conv), 1)

    ratios = []
    if has_conv:
        ratios.append(0.8)
    if bf_files:
        ratios.extend([1, 1])

    fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
    gs = fig.add_gridspec(1, len(ratios), width_ratios=ratios, wspace=0.35)
    labels = iter('abcdefgh')
    col = 0

    if has_conv:
        ax_conv = fig.add_subplot(gs[col]); col += 1
        nk.plot_newton_convergence(ax_conv, CASE_DIR)
        ax_conv.set_title('Convergence', fontsize=8)
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    if bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)

        # Velocity magnitude
        ax_vel = fig.add_subplot(gs[col]); col += 1
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax_vel, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax_vel, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax_vel.set_aspect('equal')
        nk.add_cylinder_patches(ax_vel)
        ax_vel.set_xlim(-2, 20)
        ax_vel.set_ylim(-4, 4)
        ax_vel.set_xlabel(r'$x$', labelpad=-1)
        ax_vel.set_ylabel(r'$y$', labelpad=1)
        ax_vel.spines['right'].set_visible(False)
        ax_vel.spines['top'].set_visible(False)
        nk.panel_label(ax_vel, rf'$\bf{{({next(labels)})}}$')

        # Temperature
        ax_t = fig.add_subplot(gs[col])
        if 't' in fields:
            cf2 = nk.tricontourf(ax_t, triang, fields['t'], levels=257,
                                 cmap='inferno', extend='neither')
            nk.inset_colorbar(ax_t, cf2, orientation='horizontal',
                              width="50%", height="5%", loc=9)
        ax_t.set_aspect('equal')
        nk.add_cylinder_patches(ax_t)
        ax_t.set_xlim(-2, 20)
        ax_t.set_ylim(-4, 4)
        ax_t.set_xlabel(r'$x$', labelpad=-1)
        ax_t.spines['right'].set_visible(False)
        ax_t.spines['top'].set_visible(False)
        nk.panel_label(ax_t, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
