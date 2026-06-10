#!/usr/bin/env python3
"""plot.py: Slot jet with FST Re=495 — canonical DNS panels.

Canonical panel set (DNS lane):
  plot_snapshot.png - snapshot of the DNS field (|u|, Blues)
  plot_signal.png   - vy(t) probe time series (if .his data with >= 2 rows)
  plot_fft.png      - power spectrum of wake signal (if .his data exists)
  plot_phase.png    - phase portrait vx vs vy (if .his data exists)

Geometry: slot jet, elongated domain.

OUTPUTS: plot_snapshot.png [, plot_signal.png, plot_fft.png, plot_phase.png]
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def read_his(his_path):
    """Parse Nek5000 .his file. Returns (n_probes, probes[(x,y,z)], data)."""
    with open(his_path) as f:
        n_probes = int(f.readline().strip())
        probes = []
        for _ in range(n_probes):
            coords = [float(x) for x in f.readline().split()]
            probes.append(tuple(coords[:3]))
        data = np.loadtxt(f)
    return n_probes, probes, data


def find_snapshot(case_dir):
    """Find best snapshot field file."""
    for pat in ('BF_*0.f00001', '*0.f00001'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def main():
    nk.configure_style()

    # PANEL 1 — SNAPSHOT
    snap = find_snapshot(CASE_DIR)
    if snap:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        x, y, fields, time = nk.read_field(str(snap))
        triang = nk.make_triangulation(x, y)
        if 'vx' in fields and 'vy' in fields:
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            vmax = float(np.nanpercentile(umag, 99))
            vmax = round(vmax, 1) if vmax > 0 else 1.0
            cf = nk.tricontourf(ax, triang, umag, levels=257,
                                cmap='Blues', vmin=0, vmax=vmax, extend='max')
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width='50%', height='5%', loc=9,
                              ticks=[0, vmax],
                              tick_labels=['0', f'{vmax}'])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_snapshot.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_snapshot.png"}')
        plt.close(fig)
    else:
        print('WARNING: no snapshot field found — plot_snapshot.png skipped')

    # SIGNAL / FFT / PHASE — only if .his has sufficient data
    his_candidates = sorted(CASE_DIR.glob('*.his'))
    his_path = his_candidates[0] if his_candidates else None

    data_ok = False
    if his_path is not None and his_path.exists():
        try:
            n_probes, probes, data = read_his(his_path)
            if data.ndim == 2 and data.shape[0] >= 10 and data.shape[1] >= 3:
                data_ok = True
        except Exception as e:
            print(f'WARNING: could not parse {his_path}: {e}')

    if data_ok:
        t = data[:, 0]
        vx = data[:, 1]
        vy = data[:, 2]
        probe_xyz = probes[0] if probes else (0.0, 0.0, 0.0)

        # Drop first 30% to remove initial transients
        n_skip = max(int(0.3 * len(t)), 0)
        t_s, vx_s, vy_s = t[n_skip:], vx[n_skip:], vy[n_skip:]

        # PANEL 2 — SIGNAL
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
        ax.plot(t_s, vy_s, lw=0.7, color='C0')
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r'$v_y$' + f' at ({probe_xyz[0]:.1f}, {probe_xyz[1]:.1f})')
        ax.grid(True, ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_signal.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_signal.png"}')
        plt.close(fig)

        # PANEL 3 — FFT
        # Use median of strictly positive diffs to handle duplicate timestamps
        diffs = np.diff(t_s)
        pos_diffs = diffs[diffs > 0]
        dt = float(np.median(pos_diffs)) if len(pos_diffs) > 0 else 0.0
        if dt > 0:
            fs = 1.0 / dt
            try:
                from scipy.signal import welch
                nperseg = min(len(vy_s), 2048)
                freqs, psd = welch(vy_s - vy_s.mean(), fs=fs, nperseg=nperseg)
            except ImportError:
                n = len(vy_s)
                freqs = np.fft.rfftfreq(n, d=dt)
                psd = np.abs(np.fft.rfft(vy_s - vy_s.mean()))**2 / n

            fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
            ax.semilogy(freqs, psd, lw=0.7, color='C0')
            valid = (freqs > 0.01) & (freqs < 0.5)
            if valid.any():
                i_peak = int(np.argmax(psd[valid]))
                f_peak = freqs[valid][i_peak]
                ax.axvline(f_peak, color='C3', ls='--', lw=0.6,
                           label=f'St $\\approx$ {f_peak:.3f}')
                ax.legend(fontsize=7, loc='upper right')
            ax.set_xlim(0, 0.5)
            ax.set_xlabel(r'frequency')
            ax.set_ylabel(r'PSD $|\hat{v}_y|^2$')
            ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            fig.savefig(CASE_DIR / 'plot_fft.png', dpi=600, bbox_inches='tight')
            print(f'Saved {CASE_DIR / "plot_fft.png"}')
            plt.close(fig)
        else:
            print('WARNING: uniform time spacing (all zero diffs) — plot_fft.png skipped')

        # PANEL 4 — PHASE portrait
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.7, nk.COL_WIDTH * 0.7))
        ax.plot(vx_s, vy_s, lw=0.4, color='C0', alpha=0.7)
        ax.set_xlabel(r'$v_x$' + f' at ({probe_xyz[0]:.1f}, {probe_xyz[1]:.1f})')
        ax.set_ylabel(r'$v_y$' + f' at ({probe_xyz[0]:.1f}, {probe_xyz[1]:.1f})')
        ax.set_aspect('equal')
        ax.grid(True, ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_phase.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_phase.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
