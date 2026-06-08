#!/usr/bin/env python3
"""plot.py: Visualize cylinder Re=100 DNS results.

Canonical DNS panel set:
  plot_snapshot.png  - velocity magnitude snapshot of the final field
  plot_signal.png    - vy(t) probe time series in the wake
  plot_fft.png       - power spectrum of the wake signal, Strouhal peak marked
  plot_phase.png     - phase portrait vx vs vy at the same probe
  plot.gif           - animation of vortex shedding (only if multiple snapshots)

Probe location is read from the .his header (line 2: "x y z" of the first
monitor point); defaults to (3.0, 0.0) if no .his file is present.

OUTPUTS: plot_snapshot.png, plot_signal.png, plot_fft.png, plot_phase.png, plot.gif
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def read_his(his_path):
    """Parse Nek5000 .his file. Returns (n_probes, probes[(x,y,z)], data array (t, vx_p0, vy_p0, vz_p0, vx_p1, ...))"""
    with open(his_path) as f:
        n_probes = int(f.readline().strip())
        probes = []
        for _ in range(n_probes):
            coords = [float(x) for x in f.readline().split()]
            probes.append(tuple(coords[:3]))
        data = np.loadtxt(f)
    return n_probes, probes, data


def main():
    nk.configure_style()

    # ============================================================
    # PANEL 1 — SNAPSHOT (final field, velocity magnitude)
    # ============================================================
    ff = nk.find_fields('*1cyl0.f00001', CASE_DIR)
    if not ff:
        ff = nk.find_fields('*0.f00001', CASE_DIR)

    if ff:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        x, y, fields, time = nk.read_field(ff[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
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
        re = nk.re_from_par(CASE_DIR)
        nk.re_label(ax, re)
        fig.savefig(CASE_DIR / 'plot_snapshot.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_snapshot.png"}')
        plt.close(fig)

    # ============================================================
    # SIGNAL / FFT / PHASE — read probe history
    # ============================================================
    his_path = CASE_DIR / '1cyl.his'
    if not his_path.exists():
        # Try generic glob fallback
        candidates = sorted(CASE_DIR.glob('*.his'))
        his_path = candidates[0] if candidates else None

    if his_path is not None and his_path.exists():
        n_probes, probes, data = read_his(his_path)
        # Per Nek5000 .his layout, columns after time are vx,vy,vz per probe.
        # Use the first probe by default — matches the wake-monitor convention.
        t = data[:, 0]
        vx = data[:, 1]
        vy = data[:, 2]
        probe_xyz = probes[0] if probes else (3.0, 0.0, 0.0)

        # Drop the first 30% of the signal so initial transients don't
        # dominate the spectrum or the phase portrait.
        n_keep = max(int(0.3 * len(t)), 0)
        t_s, vx_s, vy_s = t[n_keep:], vx[n_keep:], vy[n_keep:]

        # --------- PANEL 2: SIGNAL (vy time series) ---------
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

        # --------- PANEL 3: FFT (power spectrum of vy) ---------
        # Estimate sampling rate from time spacing
        dt = float(np.median(np.diff(t_s))) if len(t_s) > 1 else 1.0
        fs = 1.0 / dt
        # Use Welch's PSD for cleaner peak
        try:
            from scipy.signal import welch
            nperseg = min(len(vy_s), 2048)
            freqs, psd = welch(vy_s - vy_s.mean(), fs=fs, nperseg=nperseg)
        except ImportError:
            # numpy FFT fallback
            n = len(vy_s)
            freqs = np.fft.rfftfreq(n, d=dt)
            psd = np.abs(np.fft.rfft(vy_s - vy_s.mean()))**2 / n

        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.5))
        ax.semilogy(freqs, psd, lw=0.7, color='C0')
        # Mark the dominant peak (Strouhal frequency); restrict to physical band
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

        # --------- PANEL 4: PHASE portrait vx vs vy ---------
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.7,
                                              nk.COL_WIDTH * 0.7))
        ax.plot(vx_s, vy_s, lw=0.4, color='C0', alpha=0.7)
        ax.set_xlabel(r'$v_x$' + f' at ({probe_xyz[0]:.1f}, {probe_xyz[1]:.1f})')
        ax.set_ylabel(r'$v_y$' + f' at ({probe_xyz[0]:.1f}, {probe_xyz[1]:.1f})')
        ax.set_box_aspect(1)
        ax.grid(True, ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_phase.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_phase.png"}')
        plt.close(fig)

    # ============================================================
    # PANEL 5 — ANIM (vortex shedding GIF, only if multiple snapshots)
    # ============================================================
    snaps = sorted(CASE_DIR.glob('1cyl0.f000[1-9][0-9]'))
    snaps += sorted(CASE_DIR.glob('1cyl0.f00[1-9][0-9][0-9]'))
    snaps = [s for s in snaps if s.name >= '1cyl0.f00002']  # skip the f00001 we already plotted

    if len(snaps) >= 8:
        from matplotlib.animation import FuncAnimation, PillowWriter

        # Pre-compute once
        x0, y0, fields0, _ = nk.read_field(snaps[0])
        triang = nk.make_triangulation(x0, y0)

        fig, ax = plt.subplots(1, 1,
                               figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        def update(i):
            ax.collections.clear()
            x, y, fields, t = nk.read_field(snaps[i])
            tri = nk.make_triangulation(x, y) if (x.shape != x0.shape) else triang
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            nk.tricontourf(ax, tri, umag, levels=64,
                           cmap='Blues', vmin=0, vmax=1.5, extend='max')
            ax.set_title(f't = {t:.2f}', fontsize=7)
            return ax,

        n_frames = min(32, len(snaps))
        anim = FuncAnimation(fig, update, frames=n_frames,
                             interval=100, blit=False)
        anim.save(str(CASE_DIR / 'plot.gif'),
                  writer=PillowWriter(fps=10), dpi=120)
        print(f'Saved {CASE_DIR / "plot.gif"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
