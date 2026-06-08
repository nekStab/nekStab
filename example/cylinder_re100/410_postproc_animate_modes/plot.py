#!/usr/bin/env python3
"""plot.py: Visualize animated mode snapshots over one oscillation period.

Reads dQ_*0.f0* files (mode superimposed on base flow) and animates
all frames showing the transverse velocity (v_y).

OUTPUTS: plot.gif
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

CASE_DIR = Path(__file__).resolve().parent


def main():
    nk.configure_style()

    dq_files = sorted(CASE_DIR.glob('dQ_*0.f0*'))
    if not dq_files:
        print('No dQ_ animation files found')
        return

    nfiles = len(dq_files)

    # Global color range
    bd_global = 0
    sample_files = dq_files[::4] or dq_files
    for fpath in sample_files:
        _, _, fields, _ = nk.read_field(fpath)
        q = np.nan_to_num(fields['vy'], nan=0.0)
        bd_global = max(bd_global, np.nanpercentile(np.abs(q), 99))

    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))

    def update(idx):
        ax.clear()
        x, y, fields, _ = nk.read_field(dq_files[idx])
        triang = nk.make_triangulation(x, y)
        q = np.nan_to_num(fields['vy'], nan=0.0)
        nk.tricontourf(ax, triang, q, levels=257,
                       cmap='RdBu_r', vmin=-bd_global, vmax=bd_global,
                       extend='both')
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-1, 12)
        ax.set_ylim(-3, 3)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        phase = idx / max(nfiles - 1, 1)
        ax.set_title(rf'$\phi={phase:.2f}T$', fontsize=8)
        return []

    anim = FuncAnimation(fig, update, frames=nfiles, blit=False)
    output = CASE_DIR / 'plot.gif'
    anim.save(output, writer=PillowWriter(fps=10), dpi=150)
    print(f'Saved {output}')
    plt.close(fig)


if __name__ == '__main__':
    main()
