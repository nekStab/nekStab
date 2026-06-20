#!/usr/bin/env python3
"""plot.py: Turbulent pulsed jet Newton base flow — canonical baseflow panels.

Panel sequence (canonical baseflow lane: IC -> RESIDUAL -> BF):
  plot_ic.png       - initial condition (seed file from startFrom)
  plot_residual.png - Newton convergence history
  plot_bf.png       - converged base flow (vx, Blues)

Geometry: axisymmetric, xlim=(0,40), ylim=(0,2), axes labeled z/r.

OUTPUTS: plot_ic.png, plot_residual.png, plot_bf.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys, re
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent

XLIM = (0, 40)
YLIM = (0, 2)
FW = 2 * nk.COL_WIDTH  # double-column for elongated domain


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
    # Check all par files (skip commented-out startFrom)
    for par in par_files:
        ic_name = parse_startfrom(par)
        if ic_name:
            for d in (case_dir, case_dir.parent, case_dir.parent.parent):
                p = d / ic_name
                if p.exists():
                    return p
    return None


def find_bf(case_dir):
    for pat in ('BF_tpjet0.f00001', 'BF_Re*tpjet0.f00001',
                'BF_*tpjet0.f*', 'BF_*0.f00001'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def render_tpjet_field(ax, field_path, vmin=0, vmax=1):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    q = fields.get('vx')
    if q is None:
        return
    cf = nk.tricontourf(ax, triang, q, levels=257,
                        cmap='Blues', vmin=vmin, vmax=vmax, extend='both')
    nk.inset_colorbar(ax, cf, orientation='horizontal',
                      width='20%', height='14%', loc=1,
                      ticks=[vmin, vmax],
                      tick_labels=[str(vmin), str(vmax)])
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect('auto')
    ax.set_xlabel(r'$z$', labelpad=-1)
    ax.set_ylabel(r'$r$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def main():
    nk.configure_style()

    # PANEL 1 — IC (seed from startFrom)
    ic = find_ic(CASE_DIR)
    if ic:
        fig, ax = plt.subplots(1, 1, figsize=(FW, FW * 0.12))
        render_tpjet_field(ax, ic)
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
        fig, ax = plt.subplots(1, 1, figsize=(FW, FW * 0.12))
        render_tpjet_field(ax, bf)
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
