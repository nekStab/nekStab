#!/usr/bin/env python3
"""plot.py: Visualize Newton base flow results.

Convergence panel overlays Newton vs Newton_dyn to show that dynamic
solver tolerances reach the same residual in fewer linear solver calls.

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
NEWTON_DYN_DIR = CASE_DIR.parent / 'newton_dyn'
OUTPUT = CASE_DIR / 'plot.png'


def _read_newton(d):
    """Read Newton + Arnoldi residuals from a directory."""
    arn, nwt = None, None
    f = d / 'residu_arnoldi.dat'
    if f.exists():
        data = np.genfromtxt(str(f))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        arn = data
    f = d / 'residu_newton.dat'
    if f.exists():
        data = np.genfromtxt(str(f))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        nwt = data
    return arn, nwt


def main():
    nk.configure_style()

    arn_nwt, nwt_nwt = _read_newton(CASE_DIR)
    arn_dyn, nwt_dyn = _read_newton(NEWTON_DYN_DIR)
    has_conv = any(x is not None for x in [arn_nwt, nwt_nwt, arn_dyn, nwt_dyn])
    bf_files = nk.find_fields('BF*1cyl0.f*', CASE_DIR)

    if has_conv and bf_files:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45))
        ax_conv = None
    else:
        fig, ax_conv = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5, nk.COL_WIDTH * 0.67))
        ax_bf = None

    labels = iter('abcdefgh')

    # ── Convergence: Newton vs Newton_dyn ──────────────────────────────
    if ax_conv is not None and has_conv:
        for arn, nwt, color, tag in [
            (arn_nwt, nwt_nwt, 'C0', 'Newton'),
            (arn_dyn, nwt_dyn, 'C1', 'Newton dyn.'),
        ]:
            # Arnoldi residual (faint background curve)
            if arn is not None:
                k = np.arange(len(arn))
                ax_conv.semilogy(k, arn[:, 3], color=color, lw=0.4, alpha=0.3)

            # Newton residual markers at cumulative iteration positions
            if nwt is not None:
                k_sum = nwt[:, 3]
                res = nwt[:, 6]
                n_calls = int(k_sum[-1]) if len(k_sum) > 0 else 0
                ax_conv.semilogy(k_sum, res, 's', ms=4, mfc='none',
                                 color=color, lw=0.8,
                                 label=f'{tag} ({n_calls} calls)')
                # Convergence threshold
                if len(nwt) > 0:
                    dtol = nwt[0, 7]
                    ax_conv.axhline(dtol, color='r', ls='--', lw=0.3,
                                    zorder=1)

        ax_conv.set_xlabel('linear solver calls')
        ax_conv.set_ylabel(r'$\|r\|^2$')
        ax_conv.legend(fontsize=5.5, loc='upper right')
        ax_conv.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    # ── Base flow field ────────────────────────────────────────────────
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax_bf, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax_bf.set_aspect('equal')
        nk.add_cylinder_patches(ax_bf)
        ax_bf.set_xlim(-2, 20)
        ax_bf.set_ylim(-4, 4)
        ax_bf.set_xlabel(r'$x$', labelpad=-1)
        ax_bf.set_ylabel(r'$y$', labelpad=1)
        ax_bf.spines['right'].set_visible(False)
        ax_bf.spines['top'].set_visible(False)
        nk.panel_label(ax_bf, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
