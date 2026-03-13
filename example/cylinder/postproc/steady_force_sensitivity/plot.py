#!/usr/bin/env python3
"""plot.py: Visualize steady force sensitivity fields.

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

    sr_files = sorted(CASE_DIR.glob('sr_*0.f0*'))
    si_files = sorted(CASE_DIR.glob('si_*0.f0*'))

    panels = []
    for files, label in [(sr_files, 'Real'), (si_files, 'Imag')]:
        if files:
            x, y, fields, time = nk.read_field(files[0])
            if 'vx' in fields and 'vy' in fields:
                umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
                panels.append((label, x, y, umag))

    n = max(len(panels), 1)
    fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4 * n))
    gs = fig.add_gridspec(n, 1, hspace=0.4) if n > 1 else fig.add_gridspec(1, 1)

    labels = iter('abcdefgh')

    if panels:
        for idx, (title, x, y, q) in enumerate(panels):
            ax = fig.add_subplot(gs[idx])
            triang = nk.make_triangulation(x, y)
            vmax = np.nanpercentile(q, 99)
            cf = nk.tricontourf(ax, triang, q, levels=257,
                                cmap='hot_r', vmin=0, vmax=vmax, extend='max')
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[0, round(vmax, 2)],
                              tick_labels=['0', f'{round(vmax, 2)}'])
            ax.set_aspect('equal')
            nk.add_cylinder_patches(ax)
            ax.set_xlim(-2, 20)
            ax.set_ylim(-4, 4)
            ax.set_xlabel(r'$x$', labelpad=-1)
            ax.set_ylabel(r'$y$', labelpad=1)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            nk.panel_label(ax, rf'$\bf{{({next(labels)})}}$')
    else:
        ax = fig.add_subplot(gs[0])
        ax.text(0.5, 0.5, 'No sensitivity files found',
                transform=ax.transAxes, ha='center')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
