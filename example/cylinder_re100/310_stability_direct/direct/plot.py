#!/usr/bin/env python3
"""plot.py: Cylinder direct stability — canonical 5-panel set.

Panel sequence (canonical stability lane with twin leading-mode components):
  plot_bf.png               - converged base flow (velocity magnitude)
  plot_spectrum_circle.png  - Hessenberg eigenvalues on the unit circle
  plot_spectrum_ns_log.png  - Navier-Stokes spectrum, growth rate sigma against frequency f
  plot_mode_vx.png          - leading eigenmode, streamwise component (v_x)
  plot_mode_vy.png          - leading eigenmode, transverse component (v_y)

OUTPUTS: plot_bf.png, plot_spectrum_circle.png, plot_spectrum_ns_log.png,
         plot_mode_vx.png, plot_mode_vy.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def render_field(ax, field_path, q_name='umag', vmax=1.5, cmap='Blues', signed=False):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    if q_name == 'umag':
        q = np.sqrt(fields['vx']**2 + fields['vy']**2)
    else:
        q = fields.get(q_name, fields.get('vx'))
    if signed:
        bd = float(np.nanpercentile(np.abs(q), 99))
        cf = nk.tricontourf(ax, triang, q, levels=257, cmap=cmap,
                            vmin=-bd, vmax=bd, extend='both')
        bdr = round(bd, 2)
        nk.inset_colorbar(ax, cf, orientation='horizontal', width='50%',
                          height='5%', loc=9, ticks=[-bdr, bdr],
                          tick_labels=[f'{-bdr}', f'{bdr}'])
    else:
        cf = nk.tricontourf(ax, triang, q, levels=257, cmap=cmap,
                            vmin=0, vmax=vmax, extend='max')
        nk.inset_colorbar(ax, cf, orientation='horizontal', width='50%',
                          height='5%', loc=9, ticks=[0, vmax],
                          tick_labels=['0', f'{vmax}'])
    ax.set_aspect('equal')
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20); ax.set_ylim(-4, 4)
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def main():
    nk.configure_style()

    # PANEL 1 — BF
    bf_files = sorted(CASE_DIR.glob('BF_*0.f00001'))
    if bf_files:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        render_field(ax, bf_files[0])
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)

    # PANEL 2 — Hessenberg eigenvalues on unit circle
    spec_h = CASE_DIR / 'Spectre_Hd.dat'
    if spec_h.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.7, nk.COL_WIDTH * 0.7))
        nk.setup_unit_circle_axes(ax, lim=1.5)
        nk.plot_spectrum_H_paper(ax, spec_h, label=r'$Re=100$')
        ax.set_title(r'Hessenberg spectrum (unit circle)', fontsize=8)
        fig.savefig(CASE_DIR / 'plot_spectrum_circle.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_spectrum_circle.png"}')
        plt.close(fig)

    # PANEL 3 — Navier-Stokes spectrum, growth rate against frequency
    spec_ns = CASE_DIR / 'Spectre_NSd.dat'
    if spec_ns.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        # Load + plot manually so we can convert omega to f on the frequency axis.
        data = np.genfromtxt(str(spec_ns))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        # Columns: sigma (growth), omega (angular frequency, in U/D), residual; check headers
        if data.shape[1] >= 3:
            sig = data[:, 0]
            freq = data[:, 1] / (2 * np.pi)
            res = data[:, 2]
            converged = res < 1e-5
            ax.plot(freq[~converged], sig[~converged], ls='none', marker='o', ms=3,
                    color='lightgray', label='unconverged')
            # Leading converged value: the one with the largest growth rate and f > 0.
            lead = np.flatnonzero(converged & (freq > 0))
            lead = lead[np.argmax(sig[lead])] if len(lead) else None
            lbl = 'converged'
            if lead is not None:
                lbl = (rf'converged: $\sigma$ = {sig[lead]:.4f}, $\omega$ = {data[lead, 1]:.4f}, '
                       rf'$f$ = {freq[lead]:.4f}')
                ax.axvline(freq[lead], color='C3', lw=0.6, ls=':')
            ax.plot(freq[converged], sig[converged], ls='none', marker='*', ms=9,
                    color='C3', mec='k', mew=0.4, label=lbl)
            ax.axhline(0, color='k', lw=0.6, ls='--')
            ax.set_xlabel(r'$f = \omega/2\pi$')
            ax.set_ylabel(r'$\sigma$')
            ax.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=1)  # below the axes: inside it hides Ritz values
            ax.grid(True, ls=':', lw=0.3, alpha=0.5)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            fig.savefig(CASE_DIR / 'plot_spectrum_ns_log.png', dpi=600, bbox_inches='tight')
            print(f'Saved {CASE_DIR / "plot_spectrum_ns_log.png"}')
        plt.close(fig)

    # PANELS 4 + 5 — leading mode components (vx + vy)
    dRe_files = sorted(CASE_DIR.glob('dRe*0.f0*'))
    if dRe_files:
        for component, fname in (('vx', 'plot_mode_vx.png'),
                                 ('vy', 'plot_mode_vy.png')):
            fig, ax = plt.subplots(1, 1,
                                   figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
            render_field(ax, dRe_files[0], q_name=component,
                         cmap='RdBu', signed=True)
            fig.savefig(CASE_DIR / fname, dpi=600, bbox_inches='tight')
            print(f'Saved {CASE_DIR / fname}')
            plt.close(fig)


if __name__ == '__main__':
    main()
