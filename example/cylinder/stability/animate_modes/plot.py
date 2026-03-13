#!/usr/bin/env python3
"""plot.py: Visualize animate modes results.

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

    bf_files = sorted(CASE_DIR.glob('BF_*0.f0*'))
    if not bf_files:
        bf_files = nk.find_fields('*1cyl0.f00001', CASE_DIR)

    if bf_files:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vy', fields.get('vx'))
        if q is not None:
            bd = np.nanpercentile(np.abs(q), 99)
            cf = nk.tricontourf(ax, triang, q, levels=257,
                                cmap='RdBu', vmin=-bd, vmax=bd, extend='both')
            bdr = round(bd, 2)
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
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
