#!/usr/bin/env python3
"""plot.py: Visualize slot jet with FST results.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
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
        if 'vx' in fields and 'vy' in fields:
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            vmax = np.nanpercentile(umag, 99)
            cf = nk.tricontourf(ax, triang, umag, levels=257,
                                cmap='Blues', vmin=0, vmax=vmax, extend='max')
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[0, round(vmax, 1)],
                              tick_labels=['0', f'{round(vmax, 1)}'])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        nk.panel_label(ax, r'$\bf{(a)}$')
    else:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
        ax.text(0.5, 0.5, 'No field files found\n(run simulation first)',
                transform=ax.transAxes, ha='center', va='center')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
