#!/usr/bin/env python3
"""plot.py: Visualize thermosyphon direct stability results.

Matches AMR_Krylov_V5 paper style: annular domain, temperature fields,
inferno (base flow) + RdBu (eigenmode), equal-aspect circular geometry.

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
BF_DIR = CASE_DIR.parent / 'baseflow'
OUTPUT = CASE_DIR / 'plot.png'

# Thermosyphon geometry: R1=1, R2=2
R_INNER = 1.0
R_OUTER = 2.0


def main():
    nk.configure_style()

    # --- Determine layout based on available files ---
    spec_ns = CASE_DIR / 'Spectre_NSd.dat'
    bf_files = nk.find_fields('BF_tsyphon0.f00001', BF_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_tsyphon0.f00001', CASE_DIR)
    mode_files = nk.find_fields('dRetsyphon0.f*', CASE_DIR)

    has_fields = bool(bf_files) or bool(mode_files)

    if has_fields:
        fig = plt.figure(figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.572))
        if spec_ns.exists():
            # 3 panels: spectrum + base flow temp + mode temp
            gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.35)
            ax_spec = fig.add_subplot(gs[0])
            ax_bf = fig.add_subplot(gs[1])
            ax_mode = fig.add_subplot(gs[2])
            axes_field = [ax_bf, ax_mode]
        else:
            gs = fig.add_gridspec(1, 2, wspace=0.15)
            ax_spec = None
            ax_bf = fig.add_subplot(gs[0])
            ax_mode = fig.add_subplot(gs[1])
            axes_field = [ax_bf, ax_mode]
    else:
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH,
                                                     nk.COL_WIDTH * 0.67))
        ax_bf = ax_mode = None
        axes_field = []

    # --- Panel: NS spectrum ---
    if ax_spec is not None and spec_ns.exists():
        nk.plot_spectrum_NS_paper(ax_spec, spec_ns, label=r'$Ra=500$')
        # Shade stable half-plane
        ylim = ax_spec.get_ylim()
        ax_spec.axhspan(min(ylim[0], -0.2), 0, facecolor='gray',
                        alpha=0.3, zorder=-1)
        ax_spec.set_title(r'$\sigma$ vs $f$', fontsize=8)
        nk.panel_label(ax_spec, r'$\bf{(a)}$')

    # --- Panel: Base flow temperature ---
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)

        # Temperature is field index 5 in XUPT layout → pymech 't'
        q = fields.get('t')
        if q is not None:
            cf = nk.tricontourf(ax_bf, triang, q, levels=257,
                                cmap='inferno', vmin=0., vmax=1.,
                                extend='neither')
            cbar = nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                                     width="50%", height="5%", loc=9,
                                     ticks=[0, 1],
                                     tick_labels=['0', '1'])
        ax_bf.axis('equal')
        nk.add_annulus_patches(ax_bf, R_INNER, R_OUTER)
        _setup_annulus_axes(ax_bf, show_ylabel=True)
        nk.panel_label(ax_bf, r'$\bf{(b)}$' if ax_spec else r'$\bf{(a)}$')

    # --- Panel: Eigenmode temperature (real part) ---
    if ax_mode is not None and mode_files:
        x, y, fields, time = nk.read_field(mode_files[0])
        triang = nk.make_triangulation(x, y)

        q = fields.get('t')
        if q is not None:
            bd = np.nanpercentile(np.abs(q), 99)
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='RdBu', vmin=-bd, vmax=bd,
                                extend='both')
            bdr = round(bd, 3)
            cbar = nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                                     width="50%", height="5%", loc=9,
                                     ticks=[-bdr, bdr],
                                     tick_labels=[f'{-bdr}', f'{bdr}'])
        ax_mode.axis('equal')
        nk.add_annulus_patches(ax_mode, R_INNER, R_OUTER)
        _setup_annulus_axes(ax_mode, show_ylabel=False)
        lbl = r'$\bf{(c)}$' if ax_spec else r'$\bf{(b)}$'
        nk.panel_label(ax_mode, lbl)

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


def _setup_annulus_axes(ax, show_ylabel=True):
    """Configure axes for circular thermosyphon domain."""
    ticks = [-2, 0, 2]
    labels = ['-2', '0', '2']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels if show_ylabel else [])
    ax.set_xlabel(r'$x$', labelpad=-1)
    if show_ylabel:
        ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
