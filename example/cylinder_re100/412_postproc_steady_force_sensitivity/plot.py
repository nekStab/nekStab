#!/usr/bin/env python3
"""plot.py: Visualize the steady force sensitivity field (velocity magnitude of fsr*0.f00001).

OUTPUTS: plot_real.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def _plot_sensitivity(x, y, q, output):
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
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
    fig.savefig(output, dpi=600, bbox_inches='tight')
    print(f'Saved {output}')
    plt.close(fig)


def main():
    nk.configure_style()

    fsr_files = sorted(CASE_DIR.glob('fsr*0.f0*'))

    if fsr_files:
        x, y, fields, time = nk.read_field(fsr_files[0])
        if 'vx' in fields and 'vy' in fields:
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            _plot_sensitivity(x, y, umag, CASE_DIR / 'plot_real.png')


if __name__ == '__main__':
    main()
