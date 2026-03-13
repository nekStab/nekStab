#!/usr/bin/env python3
"""plot.py: Visualize OTD analysis results.

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

    has_dat = any((CASE_DIR / f).exists()
                  for f in ('otd_ftle.dat', 'otd_eigenvalues.dat', 'otd_growth_rates.dat'))
    bf_files = sorted(CASE_DIR.glob('bf01*0.f0*'))

    if has_dat and bf_files:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_otd = fig.add_subplot(gs[0])
        ax_mode = fig.add_subplot(gs[1])
    elif has_dat:
        fig, ax_otd = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        ax_mode = None
    elif bf_files:
        fig, ax_mode = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        ax_otd = None
    else:
        fig, ax_otd = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        ax_otd.text(0.5, 0.5, 'No OTD data found',
                    transform=ax_otd.transAxes, ha='center')
        ax_mode = None

    labels = iter('abcdefgh')

    if has_dat and ax_otd is not None:
        nk.plot_otd_exponents(ax_otd, CASE_DIR)
        ax_otd.set_title('Lyapunov exponents', fontsize=8)
        nk.panel_label(ax_otd, rf'$\bf{{({next(labels)})}}$')

    if bf_files and ax_mode is not None:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vy', fields.get('vx'))
        if q is not None:
            bd = np.nanpercentile(np.abs(q), 99)
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='RdBu', vmin=-bd, vmax=bd, extend='both')
            bdr = round(bd, 2)
            nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
        ax_mode.set_aspect('equal')
        nk.add_cylinder_patches(ax_mode)
        ax_mode.set_xlim(-2, 20)
        ax_mode.set_ylim(-4, 4)
        ax_mode.set_xlabel(r'$x$', labelpad=-1)
        ax_mode.set_ylabel(r'$y$', labelpad=1)
        ax_mode.spines['right'].set_visible(False)
        ax_mode.spines['top'].set_visible(False)
        nk.panel_label(ax_mode, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
