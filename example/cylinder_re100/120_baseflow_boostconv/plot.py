#!/usr/bin/env python3
"""plot.py: Cylinder SFD baseline (Akervik) — canonical baseflow panels.

Panel sequence (canonical baseflow lane: IC -> RESIDUAL -> BF):
  plot_ic.png       - initial condition (seed file from startFrom)
  plot_residual.png - SFD residual history with target tol
  plot_bf.png       - converged base flow (velocity magnitude)

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


def parse_startfrom(par_path):
    if not par_path.exists():
        return None
    txt = par_path.read_text()
    m = re.search(r'^\s*startFrom\s*=\s*([^\s#]+)', txt, re.M)
    if not m:
        return None
    val = m.group(1).strip().strip('"').strip("'")
    if val in ('0', ''):
        return None
    return val


def render_cyl_field(ax, field_path, what='umag', vmax=1.5, cmap='Blues'):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    if what == 'umag':
        q = np.sqrt(fields['vx']**2 + fields['vy']**2)
    else:
        q = fields.get(what, fields.get('vx'))
    cf = nk.tricontourf(ax, triang, q, levels=257,
                        cmap=cmap, vmin=0, vmax=vmax, extend='max')
    nk.inset_colorbar(ax, cf, orientation='horizontal',
                      width='50%', height='5%', loc=9,
                      ticks=[0, vmax], tick_labels=['0', f'{vmax}'])
    ax.set_aspect('equal')
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


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
    for pat in ('BF_*0.f00001', 'BF*0.f00001', '1cyl0.f*'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            # final converged BF is the highest-index file
            return matches[-1]
    return None


def main():
    nk.configure_style()

    # PANEL 1 — IC (seed from startFrom)
    ic = find_ic(CASE_DIR)
    if ic:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        render_cyl_field(ax, ic)
        fig.savefig(CASE_DIR / 'plot_ic.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_ic.png"}')
        plt.close(fig)

    # PANEL 2 — RESIDUAL (SFD residu.dat OR Newton residu_newton.dat)
    residu_sfd = CASE_DIR / 'residu.dat'
    residu_nwt = CASE_DIR / 'residu_newton.dat'
    fig = None
    if residu_sfd.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_residuals(ax, residu_sfd, label='SFD', color='C0')
        ax.axhline(1e-9, color='k', ls=':', lw=0.5, label=r'tol $10^{-9}$')
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r'$\|r\|$')
        ax.legend(fontsize=7, loc='upper right')
    elif residu_nwt.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_newton_convergence(ax, CASE_DIR)
        handles, labs = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labs, fontsize=6, loc='best')
    if fig is not None:
        ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_residual.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_residual.png"}')
        plt.close(fig)

    # PANEL 3 — BF (converged base flow)
    bf = find_bf(CASE_DIR)
    if bf:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        render_cyl_field(ax, bf)
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
