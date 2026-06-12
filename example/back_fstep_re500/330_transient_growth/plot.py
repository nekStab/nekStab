#!/usr/bin/env python3
"""plot.py: Visualize backward-facing step transient growth results.

Matches AMR_Krylov_V5 paper style: stacked panels showing BF vx (cividis),
optimal perturbation vx (RdBu), and optimal response vx (seismic),
with step geometry patch. Double-column width for elongated domain.

OUTPUTS: plot_baseflow.png, plot_optimal_perturbation.png, plot_optimal_response.png, plot_envelope.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
BF_DIR = CASE_DIR.parent / 'baseflow'

XLIM = (-5, 40)
YLIM = (-1, 2)
FW = 2 * nk.COL_WIDTH  # double-column width


def _plot_field(kind, fpath, output):
    x, y, fields, time = nk.read_field(fpath)
    q = fields.get('vx')
    if q is None:
        return

    ph = 0.15  # panel height ratio
    fig, ax = plt.subplots(1, 1, figsize=(FW, FW * ph))
    triang = nk.make_triangulation(x, y)

    cmaps = {'bf': 'cividis', 'pert': 'RdBu', 'resp': 'seismic'}
    if kind == 'bf':
        cf = nk.tricontourf(ax, triang, q, levels=257,
                            cmap=cmaps[kind], vmin=0., vmax=3.5,
                            extend='both')
        ticks = [0, 3.5]
        tick_labels = ['0', '3.5']
    else:
        bd = np.nanpercentile(np.abs(q), 99)
        cf = nk.tricontourf(ax, triang, q, levels=257,
                            cmap=cmaps[kind], vmin=-bd, vmax=bd,
                            extend='both')
        bdr = round(bd, 2)
        ticks = [-bdr, bdr]
        tick_labels = [f'{-bdr}', f'{bdr}']
    nk.inset_colorbar(ax, cf, orientation='horizontal',
                      width="20%", height="14%", loc=1,
                      ticks=ticks, tick_labels=tick_labels)

    nk.add_step_patch(ax, x0=-5, y0=-1, w=5, h=1)
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches='tight')
    print(f'Saved {output}')
    plt.close(fig)


def main():
    nk.configure_style()

    # --- Detect available field files ---
    bf_files = nk.find_fields('BF_bfs0.f00001', BF_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_bfs0.f00001', CASE_DIR)
    pert_files = nk.find_fields('pRebfs0.f*', CASE_DIR)
    resp_files = nk.find_fields('rRebfs0.f*', CASE_DIR)

    panels = []
    if bf_files:
        panels.append(('bf', bf_files[0], CASE_DIR / 'plot_baseflow.png'))
    if pert_files:
        panels.append(('pert', pert_files[0], CASE_DIR / 'plot_optimal_perturbation.png'))
    if resp_files:
        panels.append(('resp', resp_files[0], CASE_DIR / 'plot_optimal_response.png'))

    if not panels:
        # Fall back to G(t) envelope only
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        ref = CASE_DIR.parent / 'barkley2008_fig5.ref'
        ref_path = str(ref) if ref.exists() else None
        nk.plot_transient_growth(ax, CASE_DIR.parent, ref_file=ref_path)
        output = CASE_DIR / 'plot_envelope.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)
        return

    for kind, fpath, output in panels:
        _plot_field(kind, fpath, output)


if __name__ == '__main__':
    main()
