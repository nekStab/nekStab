#!/usr/bin/env python3
"""plot.py: Lid-driven cavity Re=3600 Newton base flow — canonical baseflow panels.

Panel sequence (canonical baseflow lane: IC -> RESIDUAL -> BF):
  plot_ic.png       - initial condition (seed file from startFrom)
  plot_residual.png - Newton convergence history (renamed from plot_convergence)
  plot_bf.png       - converged base flow (renamed from plot_baseflow)

Geometry: cavity on x in [-0.5, 0.5], y in [0, uparam10]; axis limits are
derived from the field coordinates (not hardcoded) so the full domain shows.

OUTPUTS: plot_ic.png, plot_residual.png, plot_bf.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys, re
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


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
    for pat in ('BF_cav0.f*', 'BF_*cav0.f*', 'BF_*0.f00001'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def render_cavity_field(ax, field_path):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
    vmax = float(np.nanpercentile(umag, 99))
    vmax = round(vmax, 2) if vmax > 0 else 1.0
    cf = nk.tricontourf(ax, triang, umag, levels=257,
                        cmap='Blues', vmin=0, vmax=vmax, extend='max')
    nk.inset_colorbar(ax, cf, orientation='horizontal',
                      width='50%', height='5%', loc=9,
                      ticks=[0, vmax],
                      tick_labels=['0', f'{vmax}'])
    ax.set_xlim(x.min(), x.max())   # domain is x in [-0.5, 0.5], not [0, 1]
    ax.set_ylim(y.min(), y.max())   # y rescaled to [0, uparam10] in usrdat2
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
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.6, nk.COL_WIDTH * 0.6))
        render_cavity_field(ax, ic)
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
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.6, nk.COL_WIDTH * 0.6))
        render_cavity_field(ax, bf)
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
