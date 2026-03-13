#!/usr/bin/env python3
"""plot.py: Visualize CI test results.

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


def main():
    nk.configure_style()

    ff = nk.find_fields('*0.f00001', CASE_DIR)

    if ff:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        x, y, fields, time = nk.read_field(ff[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        nk.panel_label(ax, r'$\bf{(a)}$')
    else:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
        ax.text(0.5, 0.5, 'No field files found',
                transform=ax.transAxes, ha='center')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
