#!/usr/bin/env python3
"""plot.py: Visualize thermosyphon direct stability results.

Matches AMR_Krylov_V5 paper style: annular domain, temperature fields,
inferno (base flow) + RdBu (eigenmode), equal-aspect circular geometry.

OUTPUTS: plot_spectrum.png, plot_baseflow.png, plot_mode.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
BF_DIR = CASE_DIR.parent / 'baseflow'

# Thermosyphon geometry: R1=1, R2=2
R_INNER = 1.0
R_OUTER = 2.0


def main():
    nk.configure_style()

    spec_ns = CASE_DIR / 'Spectre_NSd.dat'
    bf_files = nk.find_fields('BF_tsyphon0.f00001', BF_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_tsyphon0.f00001', CASE_DIR)
    mode_files = nk.find_fields('dRetsyphon0.f*', CASE_DIR)

    # --- Panel: NS spectrum ---
    if spec_ns.exists():
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_spectrum_NS_paper(ax_spec, spec_ns, label=r'$Ra=500$')
        # Shade stable half-plane
        ylim = ax_spec.get_ylim()
        ax_spec.axhspan(min(ylim[0], -0.2), 0, facecolor='gray',
                        alpha=0.3, zorder=-1)
        ax_spec.set_title(r'$\sigma$ vs $f$', fontsize=8)
        output = CASE_DIR / 'plot_spectrum.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # --- Panel: Base flow temperature ---
    if bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)

        # Temperature is field index 5 in XUPT layout → pymech 't'
        q = fields.get('t')
        if q is not None:
            fig, ax_bf = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.5))
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='inferno', vmin=0., vmax=1.,
                                extend='neither')
            nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[0, 1],
                              tick_labels=['0', '1'])
            ax_bf.axis('equal')
            nk.add_annulus_patches(ax_bf, R_INNER, R_OUTER)
            _setup_annulus_axes(ax_bf, show_ylabel=True)
            output = CASE_DIR / 'plot_baseflow.png'
            fig.savefig(output, dpi=600, bbox_inches='tight')
            print(f'Saved {output}')
            plt.close(fig)

    # --- Panel: Eigenmode temperature (real part) ---
    if mode_files:
        x, y, fields, time = nk.read_field(mode_files[0])
        triang = nk.make_triangulation(x, y)

        q = fields.get('t')
        if q is not None:
            fig, ax_mode = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                       nk.COL_WIDTH * 0.5))
            bd = np.nanpercentile(np.abs(q), 99)
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='RdBu', vmin=-bd, vmax=bd,
                                extend='both')
            bdr = round(bd, 3)
            nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
            ax_mode.axis('equal')
            nk.add_annulus_patches(ax_mode, R_INNER, R_OUTER)
            _setup_annulus_axes(ax_mode, show_ylabel=True)
            output = CASE_DIR / 'plot_mode.png'
            fig.savefig(output, dpi=600, bbox_inches='tight')
            print(f'Saved {output}')
            plt.close(fig)


def _setup_annulus_axes(ax, show_ylabel=True):
    """Configure axes for circular thermosyphon domain."""
    ticks = [-2, 0, 2]
    labels = ['-2', '0', '2']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
