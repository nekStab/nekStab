#!/usr/bin/env python3
"""plot.py: Visualize POD modal analysis results for cylinder Re=100.

Panels (emitted only when data is present):
  plot_snapshot.png  - middle DNS snapshot velocity magnitude
  plot_signal.png    - POD mode-1 temporal coefficient vs snapshot index
  plot_eigvalues.png - POD eigenvalue spectrum (energy per mode)
  plot_mode1.png     - POD spatial mode 1 (v_y)
  plot_mode2.png     - POD spatial mode 2 (v_y), if available

OUTPUTS: plot_snapshot.png, plot_signal.png, plot_eigvalues.png, plot_mode1.png [, plot_mode2.png]
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def main():
    nk.configure_style()

    # PANEL: SNAPSHOT — middle DNS snapshot
    snap_files = sorted(CASE_DIR.glob('1cyl0.f*'))
    if snap_files:
        mid = len(snap_files) // 2
        snap_file = snap_files[mid]
        x, y, fields, _ = nk.read_field(snap_file)
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        cf = nk.tricontourf(ax, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(f'DNS snapshot (frame {mid + 1}/{len(snap_files)})', fontsize=8)
        output = CASE_DIR / 'plot_snapshot.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # PANEL: SIGNAL — POD mode-1 temporal coefficient vs snapshot index
    coeff_file = CASE_DIR / 'pod_coefficients.dat'
    if coeff_file.exists():
        data = np.genfromtxt(str(coeff_file), comments='#')
        if data.ndim == 2 and data.shape[1] >= 2:
            snap_idx = data[:, 0]
            a1 = data[:, 1]  # mode-1 temporal coefficient
            fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
            ax.plot(snap_idx, a1, 'b-', lw=0.8)
            if data.shape[1] >= 3:
                a2 = data[:, 2]
                ax.plot(snap_idx, a2, 'r-', lw=0.8, alpha=0.7, label='mode 2')
                ax.legend(fontsize=6)
            ax.set_xlabel('snapshot index')
            ax.set_ylabel('POD coefficient')
            ax.set_title('POD temporal coefficients (modes 1–2)', fontsize=8)
            ax.grid(True, ls=':', lw=0.3, alpha=0.5)
            output = CASE_DIR / 'plot_signal.png'
            fig.savefig(output, dpi=600, bbox_inches='tight')
            print(f'Saved {output}')
            plt.close(fig)

    # PANEL: EIGVALUES — POD eigenvalue spectrum
    spec_file = CASE_DIR / 'pod_energy.dat'
    if not spec_file.exists():
        spec_file = CASE_DIR / 'pod_spectrum.dat'
    if spec_file.exists():
        data = np.genfromtxt(str(spec_file), comments='#')
        if data.ndim == 2 and data.shape[1] >= 3:
            modes = data[:, 0].astype(int)
            energy_pct = data[:, 2]
            cumul_pct = data[:, 3] if data.shape[1] >= 4 else np.cumsum(energy_pct)
            n_show = min(20, len(modes))
            fig, ax1 = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.6))
            ax1.bar(modes[:n_show], energy_pct[:n_show],
                    color='steelblue', edgecolor='k', linewidth=0.4, alpha=0.8)
            ax1.set_xlabel('POD mode')
            ax1.set_ylabel('Energy (%)', color='steelblue')
            ax1.tick_params(axis='y', labelcolor='steelblue')
            ax2 = ax1.twinx()
            ax2.plot(modes[:n_show], cumul_pct[:n_show], 'o-',
                     color='firebrick', ms=3, lw=0.8)
            ax2.axhline(99, color='gray', ls='--', lw=0.5, alpha=0.7)
            ax2.set_ylabel('Cumulative energy (%)', color='firebrick')
            ax2.tick_params(axis='y', labelcolor='firebrick')
            ax2.set_ylim(0, 105)
            ax1.set_title('POD eigenvalue spectrum', fontsize=8)
            output = CASE_DIR / 'plot_eigvalues.png'
            fig.savefig(output, dpi=600, bbox_inches='tight')
            print(f'Saved {output}')
            plt.close(fig)

    # PANEL: MODE 1 — POD spatial mode 1 (v_y)
    pod_mode_files = sorted(CASE_DIR.glob('pod*1cyl0.f*'))
    for idx, (mode_file, title, fname) in enumerate([
        (pod_mode_files[0] if len(pod_mode_files) >= 1 else None,
         'POD mode 1 ($v_y$)', 'plot_mode1.png'),
        (pod_mode_files[1] if len(pod_mode_files) >= 2 else None,
         'POD mode 2 ($v_y$)', 'plot_mode2.png'),
    ]):
        if mode_file is None or not mode_file.exists():
            continue
        x, y, fields, _ = nk.read_field(mode_file)
        triang = nk.make_triangulation(x, y)
        q = fields.get('vy', fields.get('vx'))
        if q is None:
            continue
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        bd = np.nanpercentile(np.abs(q), 99)
        if bd < 1e-15:
            bd = 1.0
        cf = nk.tricontourf(ax, triang, q, levels=257,
                            cmap='RdBu', vmin=-bd, vmax=bd, extend='both')
        bdr = round(bd, 2)
        nk.inset_colorbar(ax, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[-bdr, bdr],
                          tick_labels=[f'{-bdr}', f'{bdr}'])
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(title, fontsize=8)
        output = CASE_DIR / fname
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)


if __name__ == '__main__':
    main()
