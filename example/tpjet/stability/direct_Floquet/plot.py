#!/usr/bin/env python3
"""plot.py: Visualize turbulent pulsed jet direct Floquet results.

Matches AMR_Krylov_V5 paper style: Floquet spectrum on unit circle,
base flow vx (Blues, 0→5) and leading mode vx (RdBu, ±2),
axisymmetric domain [0,40]×[0,2] with z/r axis labels.
Double-column width for elongated domain.

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

XLIM = (0, 40)
YLIM = (0, 2)
FW = 2 * nk.COL_WIDTH  # double-column width


def main():
    nk.configure_style()

    # --- Detect available data ---
    spec_h = CASE_DIR / 'Spectre_Hd.dat'
    bf_files = nk.find_fields('BF_*tpjet0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_*tpjet0.f00001', BF_DIR / 'newton')
    if not bf_files:
        bf_files = nk.find_fields('BF_*tpjet0.f00001', BF_DIR / 'tdf')
    mode_files = nk.find_fields('dRetpjet0.f*', CASE_DIR)

    has_fields = bool(bf_files) or bool(mode_files)

    if has_fields and spec_h.exists():
        fig = plt.figure(figsize=(FW, FW * 0.45))
        gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.4, 0.4], hspace=0.5)
        ax_spec = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
        ax_mode = fig.add_subplot(gs[2])
    elif has_fields:
        fig = plt.figure(figsize=(FW, FW * 0.2))
        gs = fig.add_gridspec(2, 1, hspace=0.4)
        ax_spec = None
        ax_bf = fig.add_subplot(gs[0])
        ax_mode = fig.add_subplot(gs[1])
    else:
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.5))
        ax_bf = ax_mode = None

    labels = iter('abcdefgh')

    # --- Panel: Floquet spectrum ---
    if ax_spec is not None and spec_h.exists():
        nk.setup_unit_circle_axes(ax_spec, lim=1.5)
        nk.plot_spectrum_H_paper(ax_spec, spec_h)
        nk.panel_label(ax_spec, rf'$\bf{{({next(labels)})}}$')

    # --- Panel: Base flow vx (Blues, 0→5) ---
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='Blues', vmin=0., vmax=5.,
                                extend='both')
            nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                              width="20%", height="14%", loc=1,
                              ticks=[0, 5],
                              tick_labels=['0', '5'])
        ax_bf.set_xlim(XLIM)
        ax_bf.set_ylim(YLIM)
        ax_bf.set_aspect('auto')
        _setup_tpjet_axes(ax_bf)
        nk.panel_label(ax_bf, rf'$\bf{{({next(labels)})}}$')

    # --- Panel: Leading mode vx (RdBu, ±2) ---
    if ax_mode is not None and mode_files:
        x, y, fields, time = nk.read_field(mode_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='RdBu', vmin=-2., vmax=2.,
                                extend='both')
            nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                              width="20%", height="14%", loc=1,
                              ticks=[-2, 2],
                              tick_labels=['-2', '2'])
        ax_mode.set_xlim(XLIM)
        ax_mode.set_ylim(YLIM)
        ax_mode.set_aspect('auto')
        _setup_tpjet_axes(ax_mode)
        nk.panel_label(ax_mode, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


def _setup_tpjet_axes(ax):
    """Configure axes for axisymmetric tpjet domain (z/r labels)."""
    ax.set_xlabel(r'$z$', labelpad=-1)
    ax.set_ylabel(r'$r$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
