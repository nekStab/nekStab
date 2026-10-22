#!/usr/bin/env python3
"""plot.py: Animate mode snapshots for the UPO case (cylinder Re=180).

Reads dQ_*0.f0* field files (or dQb/dQm patterns) and creates an animated
GIF showing the transverse velocity (v_y) over one oscillation period.

OUTPUTS: plot.gif  (only if animation data files are present)
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

CASE_DIR = Path(__file__).resolve().parent


def find_animation_files():
    """Try multiple glob patterns for mode animation files."""
    for pat in ('dQ_*0.f0*', 'dQb*0.f0*', 'dQm*0.f0*', '*0.f0*'):
        files = sorted(CASE_DIR.glob(pat))
        if files:
            return files
    return []


def main():
    nk.configure_style()

    dq_files = find_animation_files()
    if not dq_files:
        print('No animation data found — skipping plot.gif')
        return

    nfiles = len(dq_files)

    # Global color range: sample every 4th file for speed
    bd_global = 0.0
    sample_files = dq_files[::4] or dq_files
    for fpath in sample_files:
        _, _, fields, _ = nk.read_field(fpath)
        q = np.nan_to_num(fields.get('vy', fields.get('vx', np.array([0.0]))), nan=0.0)
        bd_global = max(bd_global, float(np.nanpercentile(np.abs(q), 99)))

    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))

    def update(idx):
        ax.clear()
        x, y, fields, _ = nk.read_field(dq_files[idx])
        triang = nk.make_triangulation(x, y)
        q = np.nan_to_num(fields.get('vy', fields.get('vx', np.zeros(len(x)))), nan=0.0)
        nk.tricontourf(ax, triang, q, levels=257,
                       cmap='RdBu_r', vmin=-bd_global, vmax=bd_global,
                       extend='both')
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        phase = idx / max(nfiles - 1, 1)
        ax.set_title(rf'$\phi={phase:.2f}T$', fontsize=8)
        return []

    anim = FuncAnimation(fig, update, frames=nfiles, blit=False)
    output = CASE_DIR / 'plot.gif'
    anim.save(output, writer=PillowWriter(fps=10), dpi=100)
    print(f'Saved {output}')
    plt.close(fig)


if __name__ == '__main__':
    main()
