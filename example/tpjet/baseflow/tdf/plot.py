#!/usr/bin/env python3
"""plot.py: Visualize turbulent pulsed jet TDF base flow results.

Matches AMR_Krylov_V5 paper style: Blues colormap for base flow vx (0→5),
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
OUTPUT = CASE_DIR / 'plot.png'

XLIM = (0, 40)
YLIM = (0, 2)
FW = 2 * nk.COL_WIDTH  # double-column width


def main():
    nk.configure_style()

    # --- Detect available data ---
    res_file = CASE_DIR / 'residu_tdf.dat'
    if not res_file.exists():
        res_file = CASE_DIR / 'residu.dat'
    has_conv = res_file.exists()
    bf_files = nk.find_fields('BF_Re1900_tpjet0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_*tpjet0.f00001', CASE_DIR)

    if has_conv and bf_files:
        fig = plt.figure(figsize=(FW, FW * 0.35))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 0.5], hspace=0.45)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(FW, FW * 0.12))
        ax_conv = None
    else:
        fig, ax_conv = plt.subplots(1, 1,
                                     figsize=(nk.COL_WIDTH * 0.5,
                                              nk.COL_WIDTH * 0.67))
        ax_bf = None

    labels = iter('abcdefgh')

    # --- Panel: TDF convergence ---
    if ax_conv is not None and has_conv:
        nk.plot_residuals(ax_conv, res_file, label='TDF')
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

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
        ax_bf.set_xlabel(r'$z$', labelpad=-1)
        ax_bf.set_ylabel(r'$r$', labelpad=1)
        ax_bf.spines['right'].set_visible(False)
        ax_bf.spines['top'].set_visible(False)
        nk.panel_label(ax_bf, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
