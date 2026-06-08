#!/usr/bin/env python3
"""plot.py: Cylinder Re=100 adjoint stability — canonical stability panels.

Panel sequence (stability lane: BF -> SPECTRUM -> MODE):
  plot_bf.png       - base flow velocity magnitude (BF_1cyl0.f00001)
  plot_spectrum.png - NS eigenvalue spectrum (Spectre_NSa.dat or Spectre_Ha.dat)
  plot_mode.png     - adjoint mode (aRe1cyl0.f00001, vy field, RdBu)

Geometry: cylinder, xlim=(-2,20), ylim=(-4,4).

OUTPUTS: plot_bf.png, plot_spectrum.png, plot_mode.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def find_bf(case_dir):
    for pat in ('BF_1cyl0.f00001', 'BF_*1cyl0.f00001', 'BF_*0.f00001'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def find_spectrum(case_dir):
    """Prefer NS spectrum; fall back to H spectrum."""
    for name in ('Spectre_NSa.dat', 'Spectre_Ha.dat',
                 'Spectre_NSa_conv.dat', 'Spectre_Ha_conv.dat'):
        p = case_dir / name
        if p.exists():
            return p, name
    return None, None


def find_mode(case_dir):
    """Find adjoint direct-mode real part."""
    matches = sorted(case_dir.glob('aRe*0.f00001'))
    if not matches:
        matches = sorted(case_dir.glob('aRe*0.f0*'))
    return matches[0] if matches else None


def render_cyl_field(ax, field_path, what='umag', cmap='Blues',
                     vmin=0, vmax=1.5, extend='max', symmetric=False):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    if what == 'umag':
        q = np.sqrt(fields['vx']**2 + fields['vy']**2)
    else:
        q = fields.get(what, fields.get('vy', fields.get('vx')))
    if q is None:
        return None
    if symmetric:
        bd = float(np.nanpercentile(np.abs(q), 99))
        vmin, vmax = -bd, bd
        extend = 'both'
    cf = nk.tricontourf(ax, triang, q, levels=257,
                        cmap=cmap, vmin=vmin, vmax=vmax, extend=extend)
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return cf


def main():
    nk.configure_style()

    # PANEL 1 — BF (base flow velocity magnitude)
    bf = find_bf(CASE_DIR)
    if bf:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        cf = render_cyl_field(ax, bf, what='umag', cmap='Blues',
                              vmin=0, vmax=1.5, extend='max')
        if cf is not None:
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width='50%', height='5%', loc=9,
                              ticks=[0, 1.5], tick_labels=['0', '1.5'])
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)

    # PANEL 2 — SPECTRUM
    spec_path, spec_name = find_spectrum(CASE_DIR)
    if spec_path is not None:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.8, nk.COL_WIDTH * 0.8))
        if 'H' in spec_name:
            nk.setup_unit_circle_axes(ax, lim=1.5)
            nk.plot_spectrum_H_paper(ax, spec_path, label=r'$Re=100$')
        else:
            nk.plot_spectrum_NS_paper(ax, spec_path, label=r'$Re=100$')
            ylim = ax.get_ylim()
            ax.axhspan(min(ylim[0], -0.2), 0, facecolor='gray', alpha=0.3, zorder=-1)
            ax.set_title(r'$\sigma$ vs $f$ (adjoint)', fontsize=8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_spectrum.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_spectrum.png"}')
        plt.close(fig)

    # PANEL 3 — MODE (adjoint vy field, RdBu)
    mode = find_mode(CASE_DIR)
    if mode:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        cf = render_cyl_field(ax, mode, what='vy', cmap='RdBu',
                              symmetric=True)
        if cf is not None:
            bd = float(max(abs(cf.get_clim()[0]), abs(cf.get_clim()[1])))
            bdr = round(bd, 2)
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width='50%', height='5%', loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
        fig.savefig(CASE_DIR / 'plot_mode.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_mode.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
