#!/usr/bin/env python3
"""plot.py: Visualize turbulent pulsed jet direct Floquet results.

Matches AMR_Krylov_V5 paper style: Floquet spectrum on unit circle,
base flow vx (Blues, 0→5) and leading mode vx (RdBu, ±2),
axisymmetric domain [0,40]×[0,2] with z/r axis labels.
Double-column width for elongated domain.

OUTPUTS: plot_spectrum.png, plot_baseflow.png, plot_mode.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
BF_DIR = CASE_DIR.parent.parent / 'baseflow'

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

    # --- Panel: Floquet spectrum ---
    if spec_h.exists():
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.5))
        nk.setup_unit_circle_axes(ax_spec, lim=1.5)
        nk.plot_spectrum_H_paper(ax_spec, spec_h)
        output = CASE_DIR / 'plot_spectrum.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # --- Panel: Base flow vx (Blues, 0→5) ---
    if bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(FW, FW * 0.15))
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vx')
        if q is not None:
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='Blues', vmin=0., vmax=1.,
                                extend='both')
            nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                              width="20%", height="14%", loc=1,
                              ticks=[0, 1],
                              tick_labels=['0', '1'])
        ax_bf.set_xlim(XLIM)
        ax_bf.set_ylim(YLIM)
        ax_bf.set_aspect('auto')
        _setup_tpjet_axes(ax_bf)
        output = CASE_DIR / 'plot_baseflow.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # --- Panel: Leading mode vx (RdBu, ±2) ---
    # Use last mode file: spurious modes (|mu|>>1) converge first,
    # physical mode (|mu| ~ 1) is written last
    if mode_files:
        fig, ax_mode = plt.subplots(1, 1, figsize=(FW, FW * 0.15))
        x, y, fields, time = nk.read_field(mode_files[-1])
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
        output = CASE_DIR / 'plot_mode.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)


def _setup_tpjet_axes(ax):
    """Configure axes for axisymmetric tpjet domain (z/r labels)."""
    ax.set_xlabel(r'$z$', labelpad=-1)
    ax.set_ylabel(r'$r$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
