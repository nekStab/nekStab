#!/usr/bin/env python3
"""plot.py: Visualize forced jet adjoint Floquet results.

Matches AMR_Krylov_V5 paper style: Floquet spectrum on unit circle,
base flow vx (Blues, 0→1) and leading adjoint mode vx (RdBu, ±6),
axisymmetric domain [0,40]×[0,2] with z/r axis labels. The adjoint mode
lives at the nozzle, so its panel shows [0,4]×[0,1].
Double-column width for elongated domain.

OUTPUTS: plot_spectrum.png, plot_baseflow.png, plot_mode.png (a: unfiltered, b: post-filtered)
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(next(d for d in Path(__file__).resolve().parents if (d / 'nekplot.py').exists())))
import nekplot as nk
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
BF_DIR = CASE_DIR.parent.parent / 'baseflow'

XLIM = (0, 40)
YLIM = (0, 2)
FW = 2 * nk.COL_WIDTH  # double-column width
MODE_XLIM = (0, 4)  # the adjoint mode lives at the nozzle
MODE_YLIM = (0, 1)


def main():
    nk.configure_style()

    # --- Detect available data ---
    spec_h = CASE_DIR / 'Spectre_Ha.dat'
    bf_files = nk.find_fields('BF_*tpjet0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_*tpjet0.f00001', BF_DIR / 'newton')
    if not bf_files:
        bf_files = nk.find_fields('BF_*tpjet0.f00001', BF_DIR / 'tdf')
    mode_files = nk.find_fields('aRetpjet0.f*', CASE_DIR)

    # --- Panel: Floquet spectrum ---
    if spec_h.exists():
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.5))
        nk.setup_unit_circle_axes(ax_spec, lim=1.5)
        nk.plot_spectrum_H_paper(ax_spec, spec_h)
        output = CASE_DIR / 'plot_spectrum.png'
        fig.savefig(output, dpi=400, bbox_inches='tight')
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
        fig.savefig(output, dpi=400, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # --- Panel: Leading adjoint mode vx (RdBu, ±6), zoomed on the nozzle ---
    # The mode files hold the converged multipliers by decreasing modulus: the first is the leading mode.
    if mode_files:
        # (a) the file as written by the solver, (b) the same field after the post-filter of
        # nekplot.modal_lowpass (top Legendre modes damped in each element, shared points averaged).
        # The filter is applied here only; the field file is not changed.
        x, y, fields, time = nk.read_field(mode_files[0])
        # Strong filter for the adjoint: the lip streaks span several elements (reviewer asked for a stronger one).
        filtered = nk.modal_lowpass(mode_files[0], ('vx',), top=(0.5, 0.2, 0.0), passes=3)
        triang = nk.make_triangulation(x, y)
        fig, axes = plt.subplots(2, 1, figsize=(FW * 0.5, FW * 0.5 * 0.3 * 2), sharex=True)
        fig.subplots_adjust(hspace=0.5)  # room for the (b) subtitle
        for ax_mode, q, tag, name in ((axes[0], fields.get('vx'), 'a', 'unfiltered'),
                                      (axes[1], filtered['vx'], 'b', 'filtered')):
            if q is None:
                continue
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='RdBu', vmin=-6., vmax=6.,
                                extend='both')
            cbar = nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                                     width="20%", height="14%", loc=1,
                                     ticks=[-6, 6],
                                     tick_labels=['-6', '6'])
            cbar.set_label(r'$v_z$', fontsize=5, labelpad=1)
            ax_mode.set_xlim(MODE_XLIM)
            ax_mode.set_ylim(MODE_YLIM)
            ax_mode.set_aspect('auto')
            ax_mode.set_title(f'({tag}) leading adjoint mode, {name}', loc='left', fontsize=7)
            _setup_tpjet_axes(ax_mode)
        axes[0].set_xlabel('')  # shared z axis: label only under (b)
        output = CASE_DIR / 'plot_mode.png'
        fig.savefig(output, dpi=400, bbox_inches='tight')
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
