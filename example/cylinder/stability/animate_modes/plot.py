#!/usr/bin/env python3
"""plot.py: Visualize animated mode snapshots over one oscillation period.

Reads dQ_*0.f0* files (mode superimposed on base flow) and plots
16 evenly spaced frames showing the transverse velocity (v_y).

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

NFRAMES = 16


def main():
    nk.configure_style()

    dq_files = sorted(CASE_DIR.glob('dQ_*0.f0*'))
    if not dq_files:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
        ax.text(0.5, 0.5, 'No dQ_ animation files found',
                transform=ax.transAxes, ha='center')
        fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
        plt.close()
        return

    nfiles = len(dq_files)
    indices = np.linspace(0, nfiles - 1, NFRAMES, dtype=int)

    ncols = 4
    nrows = NFRAMES // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4 * nrows))
    axes_flat = axes.ravel()

    labels = iter('abcdefghijklmnopqrst')

    # Global color range
    bd_global = 0
    for idx in indices:
        _, _, fields, _ = nk.read_field(dq_files[idx])
        q = np.nan_to_num(fields['vy'], nan=0.0)
        bd_global = max(bd_global, np.nanpercentile(np.abs(q), 99))

    cf = None
    for i, idx in enumerate(indices):
        ax = axes_flat[i]
        x, y, fields, _ = nk.read_field(dq_files[idx])
        triang = nk.make_triangulation(x, y, fields)
        q = np.nan_to_num(fields['vy'], nan=0.0)
        cf = nk.tricontourf(ax, triang, q, levels=257,
                            cmap='RdBu_r', vmin=-bd_global, vmax=bd_global,
                            extend='both')
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-1, 12)
        ax.set_ylim(-3, 3)

        row, col = divmod(i, ncols)
        if row == nrows - 1:
            ax.set_xlabel(r'$x$', labelpad=-1)
        else:
            ax.set_xticklabels([])
        if col == 0:
            ax.set_ylabel(r'$y$', labelpad=1)
        else:
            ax.set_yticklabels([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        phase = idx / max(nfiles - 1, 1)
        lbl = next(labels)
        nk.panel_label(ax, rf'$\bf{{({lbl})}}$  $\phi={phase:.2f}T$')

    if cf is not None:
        bdr = round(bd_global, 2)
        cbar = fig.colorbar(cf, ax=axes_flat.tolist(), orientation='horizontal',
                            fraction=0.03, pad=0.08, aspect=50)
        cbar.set_ticks([-bdr, 0, bdr])
        cbar.set_ticklabels([f'{-bdr}', '0', f'{bdr}'])
        cbar.set_label(r'$v_y$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
