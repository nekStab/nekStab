#!/usr/bin/env python3
"""plot.py: Visualize structural sensitivity / wavemaker field.

OUTPUTS: plot_wavemaker.png, plot_direct_mode.png, plot_adjoint_mode.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def _plot_panel(kind, x, y, q, output):
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
    triang = nk.make_triangulation(x, y)
    if kind == 'wavemaker':
        vmax = np.nanpercentile(q, 99)
        cf = nk.tricontourf(ax, triang, q, levels=257,
                            cmap='hot_r', vmin=0, vmax=vmax, extend='max')
        nk.inset_colorbar(ax, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, round(vmax, 3)],
                          tick_labels=['0', f'{round(vmax, 3)}'])
    else:
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
    fig.savefig(output, dpi=600, bbox_inches='tight')
    print(f'Saved {output}')
    plt.close(fig)


def main():
    nk.configure_style()

    dRe_files = sorted(CASE_DIR.glob('dRe*0.f0*'))
    aRe_files = sorted(CASE_DIR.glob('aRe*0.f0*'))

    if dRe_files and aRe_files:
        xd, yd, fd, td = nk.read_field(dRe_files[0])
        xa, ya, fa, ta = nk.read_field(aRe_files[0])

        if 'vx' in fd and 'vy' in fd and 'vx' in fa and 'vy' in fa:
            mag_d = np.sqrt(fd['vx']**2 + fd['vy']**2)
            mag_a = np.sqrt(fa['vx']**2 + fa['vy']**2)
            wavemaker = mag_d * mag_a
            _plot_panel('wavemaker', xd, yd, wavemaker, CASE_DIR / 'plot_wavemaker.png')
            _plot_panel('direct', xd, yd, fd.get('vy', fd['vx']), CASE_DIR / 'plot_direct_mode.png')
            _plot_panel('adjoint', xa, ya, fa.get('vy', fa['vx']), CASE_DIR / 'plot_adjoint_mode.png')


if __name__ == '__main__':
    main()
