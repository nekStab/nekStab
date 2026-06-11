#!/usr/bin/env python3
"""Analyze the moving-cylinder Re=100 DNS stage."""
from __future__ import annotations

import json
import math
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk


CASE_DIR = Path(__file__).resolve().parent
LATE_FRACTION = 0.5


def read_his(path: Path) -> tuple[list[tuple[float, float, float]], np.ndarray]:
    with path.open() as handle:
        n_probes = int(handle.readline().strip())
        probes = []
        for _ in range(n_probes):
            fields = handle.readline().split()
            probes.append(tuple(float(v) for v in fields[:3]))
        data = np.loadtxt(handle)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return probes, data


def final_snapshot() -> Path:
    matches = sorted(CASE_DIR.glob("1cyl0.f*"))
    if not matches:
        raise FileNotFoundError("no 1cyl0.f* field output found")
    return matches[-1]


def dominant_frequency(t: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    dt = float(np.median(np.diff(t)))
    nfft = 1
    while nfft < len(y) * 8:
        nfft *= 2
    window = np.hanning(len(y))
    freqs = np.fft.rfftfreq(nfft, d=dt)
    power = np.abs(np.fft.rfft((y - y.mean()) * window, n=nfft)) ** 2
    band = (freqs > 0.03) & (freqs < 0.5)
    if not np.any(band):
        raise ValueError("no FFT bins in the expected shedding band")
    idx = np.where(band)[0][np.argmax(power[band])]
    return float(freqs[idx]), freqs, power


def render_snapshot(path: Path) -> None:
    nk.configure_style()
    x, y, fields, time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    if "t" in fields:
        values = fields["t"]
        cmap = "RdBu_r"
        bound = float(np.nanpercentile(np.abs(values), 99.0))
        vmin, vmax, extend = -bound, bound, "both"
        ticks = [round(-bound, 3), round(bound, 3)]
    else:
        values = np.sqrt(fields["vx"]**2 + fields["vy"]**2)
        cmap = "viridis"
        vmax = float(np.nanpercentile(values, 99.0))
        vmin, extend = 0.0, "max"
        ticks = [0.0, round(vmax, 3)]

    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.6, nk.COL_WIDTH * 0.7))
    cf = nk.tricontourf(ax, triang, values, levels=257, cmap=cmap,
                        vmin=vmin, vmax=vmax, extend=extend)
    nk.inset_colorbar(ax, cf, orientation="horizontal", width="46%",
                      height="5%", loc=9, ticks=ticks,
                      tick_labels=[f"{tick:g}" for tick in ticks])
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 8)
    ax.set_ylim(-3, 3)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(f"t = {time:.2f}")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(CASE_DIR / "snapshot.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_reference(frequency: float, amplitude: float) -> None:
    spec = {
        "case": "moving_cylinder_re100/000_dns",
        "description": "Moving-cylinder Re=100 DNS late-window dominant frequency and amplitude from 1cyl.his.",
        "produced": {
            "date": "2026-06-10",
            "ranks": 8,
            "mode": "userParam01=0",
        },
        "quantities": {
            "dominant_frequency": {
                "value": frequency,
                "source": "analysis_metrics.dat:col1:first",
                "tol_rel": 1.0e-3,
            },
            "late_amplitude": {
                "value": amplitude,
                "source": "analysis_metrics.dat:col2:first",
                "tol_rel": 0.02,
            },
        },
    }
    ref_dir = CASE_DIR / "ref"
    ref_dir.mkdir(exist_ok=True)
    (ref_dir / "reference.json").write_text(json.dumps(spec, indent=2) + "\n")


def main() -> None:
    nk.configure_style()
    probes, data = read_his(CASE_DIR / "1cyl.his")
    t = data[:, 0]
    vy = data[:, 2]
    start = t[0] + LATE_FRACTION * (t[-1] - t[0])
    late = t >= start
    t_late = t[late]
    vy_late = vy[late]
    frequency, freqs, power = dominant_frequency(t_late, vy_late)
    amplitude = float(0.5 * (np.max(vy_late) - np.min(vy_late)))

    np.savetxt(
        CASE_DIR / "analysis_metrics.dat",
        np.array([[frequency, amplitude]]),
        header="dominant_frequency late_amplitude",
    )

    probe = probes[0]
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.35, nk.COL_WIDTH * 0.65))
    ax.plot(t, vy, lw=0.6)
    ax.axvspan(t_late[0], t_late[-1], color="0.9", zorder=-1)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$v_y$ at " + f"({probe[0]:.1f}, {probe[1]:.1f})")
    ax.grid(True, ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(CASE_DIR / "signal.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.35, nk.COL_WIDTH * 0.65))
    ax.semilogy(freqs, power, lw=0.7)
    ax.axvline(frequency, color="C3", ls="--", lw=0.8,
               label=f"St = {frequency:.5f}")
    ax.set_xlim(0.03, 0.5)
    ax.set_xlabel("frequency")
    ax.set_ylabel(r"$|\hat{v}_y|^2$")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(CASE_DIR / "spectrum.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    render_snapshot(final_snapshot())
    write_reference(frequency, amplitude)

    ref_dir = CASE_DIR / "ref"
    for name in ("signal.png", "spectrum.png", "snapshot.png"):
        shutil.copy2(CASE_DIR / name, ref_dir / name)

    print(f"dominant_frequency={frequency:.8g}")
    print(f"late_amplitude={amplitude:.8g}")
    print(f"late_window={t_late[0]:.8g}:{t_late[-1]:.8g}")


if __name__ == "__main__":
    main()
