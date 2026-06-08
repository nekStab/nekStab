#!/usr/bin/env python3
"""plot.py: Flip-flop dual-cylinder Newton base flow — canonical baseflow panels.

Panel sequence (canonical baseflow lane: IC -> RESIDUAL -> BF):
  plot_ic.png       - initial condition (seed file from startFrom)
  plot_residual.png - Newton convergence history
  plot_bf.png       - converged base flow (vx, PiYG diverging)

Geometry: dual cylinders, gap=0.7, radius=0.5.

OUTPUTS: plot_ic.png, plot_residual.png, plot_bf.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys, re
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent

# Flip-flop dual cylinders
GAP = 0.7
RADIUS = 0.5

XLIM = (-5, 30)
YLIM = (-8, 8)


def parse_startfrom(par_path):
    if not par_path.exists():
        return None
    txt = par_path.read_text()
    m = re.search(r'^\s*startFrom\s*=\s*([^\s#]+)', txt, re.M | re.I)
    if not m:
        return None
    val = m.group(1).strip().strip('"').strip("'")
    if val in ('0', ''):
        return None
    return val


def find_ic(case_dir):
    par_files = list(case_dir.glob('*.par'))
    if not par_files:
        return None
    ic_name = parse_startfrom(par_files[0])
    if not ic_name:
        return None
    for d in (case_dir, case_dir.parent, case_dir.parent.parent):
        p = d / ic_name
        if p.exists():
            return p
    return None


def find_bf(case_dir):
    # Try specific names first, then fallback patterns
    for pat in ('BF_2cyl0.f00001', 'BF_Re60_2cyl0.f00001',
                'BF_*2cyl0.f*', 'BF_*0.f00001'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def render_flipflop_field(ax, field_path):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    q = fields.get('vx')
    if q is None:
        return
    cf = nk.tricontourf(ax, triang, q, levels=257,
                        cmap='PiYG', vmin=-1.5, vmax=1.5, extend='both')
    nk.inset_colorbar(ax, cf, orientation='horizontal',
                      width='50%', height='6%', loc=9,
                      ticks=[-1.5, 1.5],
                      tick_labels=['-1.5', '1.5'])
    nk.add_dual_cylinder_patches(ax, gap=GAP, radius=RADIUS)
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def main():
    nk.configure_style()

    # PANEL 1 — IC (seed from startFrom)
    ic = find_ic(CASE_DIR)
    if ic:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.5))
        render_flipflop_field(ax, ic)
        fig.savefig(CASE_DIR / 'plot_ic.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_ic.png"}')
        plt.close(fig)

    # PANEL 2 — RESIDUAL (Newton convergence)
    residu_nwt = CASE_DIR / 'residu_newton.dat'
    residu_arn = CASE_DIR / 'residu_arnoldi.dat'
    if residu_nwt.exists() or residu_arn.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_newton_convergence(ax, CASE_DIR)
        handles, labs = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labs, fontsize=6, loc='best')
        ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_residual.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_residual.png"}')
        plt.close(fig)

    # PANEL 3 — BF (converged base flow)
    bf = find_bf(CASE_DIR)
    if bf:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.5))
        render_flipflop_field(ax, bf)
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
