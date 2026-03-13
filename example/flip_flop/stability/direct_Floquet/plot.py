#!/usr/bin/env python3
"""plot.py: Visualize flip-flop direct Floquet results.

Matches AMR_Krylov_V5 paper style: Floquet spectrum on unit circle,
UPO snapshot vx (PiYG ±1.5) and leading mode vx (PiYG ±0.5),
dual cylinder patches (gap=0.7, radius=0.5).

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
BF_DIR = CASE_DIR.parent.parent / 'baseflow'
OUTPUT = CASE_DIR / 'plot.png'

# Flip-flop dual cylinders
GAP = 0.7
RADIUS = 0.5


def main():
    nk.configure_style()

    # --- Detect available data ---
    spec_h = CASE_DIR / 'Spectre_Hd.dat'
    bf_files = nk.find_fields('BF_*2cyl0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_*2cyl0.f00001', BF_DIR)
    mode_files = nk.find_fields('dRe2cyl0.f*', CASE_DIR)

    has_fields = bool(bf_files) or bool(mode_files)

    if has_fields:
        ncols = int(spec_h.exists()) + int(bool(bf_files)) + int(bool(mode_files))
        ratios = []
        if spec_h.exists():
            ratios.append(0.8)
        if bf_files:
            ratios.append(1)
        if mode_files:
            ratios.append(1)
        fig = plt.figure(figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.4))
        gs = fig.add_gridspec(1, ncols, width_ratios=ratios, wspace=0.3)
        col = 0
        ax_spec = None
        if spec_h.exists():
            ax_spec = fig.add_subplot(gs[col]); col += 1
        ax_upo = None
        if bf_files:
            ax_upo = fig.add_subplot(gs[col]); col += 1
        ax_mode = None
        if mode_files:
            ax_mode = fig.add_subplot(gs[col])
    else:
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.5))
        ax_upo = ax_mode = None

    labels = iter('abcdefgh')

    # --- Panel: Floquet spectrum ---
    if ax_spec is not None and spec_h.exists():
        nk.setup_unit_circle_axes(ax_spec, lim=1.5)
        nk.plot_spectrum_H_paper(ax_spec, spec_h, label=r'$Re = 60$')
        nk.panel_label(ax_spec, rf'$\bf{{({next(labels)})}}$')

    # --- Panel: UPO snapshot (vx, PiYG ±1.5) ---
    if ax_upo is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_upo, triang, q, levels=257,
                                cmap='PiYG', vmin=-1.5, vmax=1.5,
                                extend='both')
            nk.inset_colorbar(ax_upo, cf, orientation='horizontal',
                              width="50%", height="6%", loc=9,
                              ticks=[-1.5, 1.5],
                              tick_labels=['-1.5', '1.5'])
        ax_upo.set_aspect('equal')
        nk.add_dual_cylinder_patches(ax_upo, gap=GAP, radius=RADIUS)
        _setup_flipflop_axes(ax_upo, show_ylabel=True)
        nk.panel_label(ax_upo, rf'$\bf{{({next(labels)})}}$')

    # --- Panel: Leading Floquet mode (vx, PiYG ±0.5) ---
    if ax_mode is not None and mode_files:
        x, y, fields, time = nk.read_field(mode_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='PiYG', vmin=-0.5, vmax=0.5,
                                extend='both')
            nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                              width="50%", height="6%", loc=9,
                              ticks=[-0.5, 0.5],
                              tick_labels=['-0.5', '0.5'])
        ax_mode.set_aspect('equal')
        nk.add_dual_cylinder_patches(ax_mode, gap=GAP, radius=RADIUS)
        _setup_flipflop_axes(ax_mode, show_ylabel=False)
        nk.panel_label(ax_mode, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


def _setup_flipflop_axes(ax, show_ylabel=True):
    """Configure axes for flip-flop domain."""
    ax.set_xlim(-5, 30)
    ax.set_ylim(-8, 8)
    ax.set_xlabel(r'$x$', labelpad=-1)
    if show_ylabel:
        ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
