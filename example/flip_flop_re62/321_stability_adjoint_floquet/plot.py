#!/usr/bin/env python3
"""plot.py: Visualize flip-flop adjoint Floquet results.

Matches AMR_Krylov_V5 paper style: Floquet spectrum on unit circle,
UPO snapshot vx (PiYG ±1.5) and leading mode vx (PiYG ±0.5),
dual cylinder patches (gap=0.7, radius=0.5).

OUTPUTS: ref/spectrum.png, ref/upo_snapshot.png, ref/mode.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
BF_DIR = CASE_DIR.parent.parent / 'baseflow'
REF_DIR = CASE_DIR / 'ref'

# Flip-flop dual cylinders
GAP = 0.7
RADIUS = 0.5


def main():
    nk.configure_style()
    REF_DIR.mkdir(exist_ok=True)

    # --- Detect available data ---
    spec_h = CASE_DIR / 'Spectre_Ha.dat'
    bf_files = nk.find_fields('BF_*2cyl0.f00001', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_*2cyl0.f00001', BF_DIR)
    mode_files = nk.find_fields('aRe2cyl0.f*', CASE_DIR)

    # --- Panel: Floquet spectrum ---
    if spec_h.exists():
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                                     nk.COL_WIDTH * 0.5))
        nk.setup_unit_circle_axes(ax_spec, lim=1.5)
        nk.plot_spectrum_H_paper(ax_spec, spec_h, label=r'$Re = 60$')
        output = REF_DIR / 'spectrum.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # --- Panel: UPO snapshot (vx, PiYG ±1.5) ---
    if bf_files:
        fig, ax_upo = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2,
                                                   nk.COL_WIDTH * 0.4))
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
        _setup_flipflop_axes(ax_upo)
        output = REF_DIR / 'upo_snapshot.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # --- Panel: Leading Floquet mode (vx, PiYG ±0.5) ---
    if mode_files:
        fig, ax_mode = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2,
                                                    nk.COL_WIDTH * 0.4))
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
        _setup_flipflop_axes(ax_mode)
        output = REF_DIR / 'mode.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)


def _setup_flipflop_axes(ax):
    """Configure axes for flip-flop domain."""
    ax.set_xlim(-5, 30)
    ax.set_ylim(-8, 8)
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


if __name__ == '__main__':
    main()
