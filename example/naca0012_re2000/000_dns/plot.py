#!/usr/bin/env python3
"""Analyze the NACA 0012 Re=2000 DNS stage."""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.signal import hilbert

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk

CASE_DIR = Path(__file__).resolve().parent
LINEAR_SIGMA = -0.1127203
LINEAR_OMEGA = 8.036574
RINGDOWN_WINDOW = (30.0, 70.0)
AMPLITUDE_RATIO_WINDOWS = ((30.0, 40.0), (60.0, 70.0))


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


def find_snapshot() -> Path:
    patterns = ("naca00120.f*", "resnaca00120.f*", "*naca00120.f*")
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(CASE_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError("no NACA field file found")
    return sorted(set(matches))[-1]


def peak_frequency(t: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    dt = float(np.median(np.diff(t)))
    nfft = 1
    while nfft < len(y) * 16:
        nfft *= 2
    window = np.hanning(len(y))
    freqs = np.fft.rfftfreq(nfft, d=dt)
    power = np.abs(np.fft.rfft((y - y.mean()) * window, n=nfft)) ** 2
    band = (freqs > 0.5) & (freqs < 2.0)
    if not np.any(band):
        raise ValueError("no FFT bins in the expected wake band")
    idx = np.where(band)[0][np.argmax(power[band])]
    if 0 < idx < len(power) - 1 and power[idx - 1] > 0 and power[idx + 1] > 0:
        alpha = math.log(power[idx - 1])
        beta = math.log(power[idx])
        gamma = math.log(power[idx + 1])
        denom = alpha - 2.0 * beta + gamma
        delta = 0.5 * (alpha - gamma) / denom if denom != 0.0 else 0.0
        f_peak = float(freqs[idx] + delta * (freqs[1] - freqs[0]))
    else:
        f_peak = float(freqs[idx])
    return f_peak, freqs, power


def fit_decay(t: np.ndarray, y: np.ndarray) -> float:
    envelope = np.abs(hilbert(y))
    keep = envelope > max(float(np.max(envelope)) * 1.0e-3, 1.0e-12)
    slope, _ = np.polyfit(t[keep], np.log(envelope[keep]), 1)
    return float(slope)


def estimate_frequency_from_phase(t: np.ndarray, y: np.ndarray) -> float:
    phase = np.unwrap(np.angle(hilbert(y)))
    slope, _ = np.polyfit(t, phase, 1)
    return float(slope / (2.0 * math.pi))


def peak_to_peak_amplitude(y: np.ndarray) -> float:
    return float(0.5 * (np.max(y) - np.min(y)))


def window_amplitude_ratio_sigma(t: np.ndarray, y: np.ndarray) -> float:
    (a0, a1), (b0, b1) = AMPLITUDE_RATIO_WINDOWS
    first = (t >= a0) & (t <= a1)
    last = (t >= b0) & (t <= b1)
    amp_first = peak_to_peak_amplitude(y[first])
    amp_last = peak_to_peak_amplitude(y[last])
    t_first = 0.5 * (a0 + a1)
    t_last = 0.5 * (b0 + b1)
    return float(math.log(amp_last / amp_first) / (t_last - t_first))


def render_snapshot(path: Path) -> None:
    nk.configure_style()
    x, y, fields, time = nk.read_field(str(path))
    triang = nk.make_triangulation(x, y)
    umag = np.sqrt(fields["vx"] ** 2 + fields["vy"] ** 2)

    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.55))
    cf = nk.tricontourf(ax, triang, umag, levels=257, cmap="Blues",
                        vmin=0, vmax=1.5, extend="max")
    nk.inset_colorbar(ax, cf, orientation="horizontal", width="50%",
                      height="5%", loc=9, ticks=[0, 1.5],
                      tick_labels=["0", "1.5"])
    nk.add_naca0012_patch(ax)
    ax.set_xlim(-1, 5)
    ax.set_ylim(-2, 2)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(f"t = {time:.2f}")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(CASE_DIR / "snapshot.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_reference(outcome: str, frequency: float, omega: float,
                    sigma_hilbert: float, sigma_ratio: float,
                    amplitude: float) -> None:
    quantities = {
        "dominant_frequency": {
            "source": "analysis_metrics.dat:col1:first",
            "value": frequency,
            "tol_rel": 0.02,
        },
        "late_amplitude": {
            "source": "analysis_metrics.dat:col2:first",
            "value": amplitude,
            "tol_rel": 0.02,
        },
        "fitted_sigma": {
            "source": "analysis_metrics.dat:col3:first",
            "value": sigma_hilbert,
            "tol_rel": 0.02,
        },
        "fitted_omega": {
            "source": "analysis_metrics.dat:col4:first",
            "value": omega,
            "tol_rel": 0.02,
        },
        "sigma_window_ratio": {
            "source": "analysis_metrics.dat:col5:first",
            "value": sigma_ratio,
            "tol_rel": 0.02,
        },
    }
    description = (
        "NACA 0012 Re=2000 DNS from finite-amplitude Re=2000 wake restart; "
        f"outcome={outcome}. The late-window DNS frequency agrees with the "
        "linear eigenvalue omega=8.036574 to about 0.4%, but the decay rate "
        "differs from the linear prediction sigma=-0.1127203 and remains under "
        "investigation. The DNS steady state differs from the Newton base flow "
        "by O(1e-2) in the wake and O(0.2) in the sponge zone."
    )
    spec = {
        "case": "naca0012_re2000/000_dns",
        "description": description,
        "produced": {
            "date": "2026-06-10",
            "ranks": "8",
            "nekstab_mode": "0 (DNS)",
            "initial_condition": "resnaca00120.f00010",
        },
        "quantities": quantities,
    }
    ref_dir = CASE_DIR / "ref"
    ref_dir.mkdir(exist_ok=True)
    (ref_dir / "reference.json").write_text(json.dumps(spec, indent=2) + "\n")


def main() -> None:
    nk.configure_style()
    probes, data = read_his(CASE_DIR / "naca0012.his")
    t = data[:, 0]
    vy = data[:, 2]
    final_mean = float(np.mean(vy[t >= 65.0]))
    residual = vy - final_mean
    fit_mask = (t >= RINGDOWN_WINDOW[0]) & (t <= RINGDOWN_WINDOW[1])
    t_late = t[fit_mask]
    vy_late = residual[fit_mask]

    fft_frequency, freqs, power = peak_frequency(t_late, vy_late)
    frequency = estimate_frequency_from_phase(t_late, vy_late)
    omega = 2.0 * math.pi * frequency
    sigma_hilbert = fit_decay(t_late, vy_late)
    sigma_ratio = window_amplitude_ratio_sigma(t, residual)
    amplitude = peak_to_peak_amplitude(residual[t >= 65.0])
    omega_err = abs((omega - LINEAR_OMEGA) / LINEAR_OMEGA)
    sigma_estimates_agree = abs((sigma_hilbert - sigma_ratio) / sigma_hilbert) <= 0.05
    outcome = "ringdown" if sigma_hilbert < 0.0 and sigma_ratio < 0.0 and sigma_estimates_agree else "limit-cycle"

    np.savetxt(
        CASE_DIR / "analysis_metrics.dat",
        np.array([[frequency, amplitude, sigma_hilbert, omega, sigma_ratio, fft_frequency]]),
        header="frequency amplitude sigma_hilbert omega sigma_window_ratio fft_frequency",
    )

    probe = probes[0]
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.4, nk.COL_WIDTH * 0.65))
    ax.plot(t, vy, lw=0.7)
    ax.axvspan(t_late[0], t_late[-1], color="0.9", zorder=-1)
    ax.text(0.02, 0.95, f"fit window: {RINGDOWN_WINDOW[0]:.0f}-{RINGDOWN_WINDOW[1]:.0f}",
            transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$v_y$ at " + f"({probe[0]:.1f}, {probe[1]:.1f})")
    ax.grid(True, ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(CASE_DIR / "signal.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.4, nk.COL_WIDTH * 0.65))
    ax.semilogy(freqs, power, lw=0.7)
    ax.axvline(frequency, color="C3", ls="--", lw=0.8,
               label=f"DNS f = {frequency:.4f}")
    ax.axvline(LINEAR_OMEGA / (2.0 * math.pi), color="C2", ls=":", lw=0.8,
               label=f"linear f = {LINEAR_OMEGA / (2.0 * math.pi):.4f}")
    ax.set_xlim(0, 2.0)
    ax.set_xlabel("frequency")
    ax.set_ylabel(r"$|\hat{v}_y|^2$")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", ls=":", lw=0.3, alpha=0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(CASE_DIR / "spectrum.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    render_snapshot(find_snapshot())
    write_reference(outcome, frequency, omega, sigma_hilbert, sigma_ratio, amplitude)

    for name in ("signal.png", "spectrum.png", "snapshot.png"):
        (CASE_DIR / "ref" / name).write_bytes((CASE_DIR / name).read_bytes())

    print(f"outcome={outcome}")
    print(f"frequency={frequency:.8g}")
    print(f"fft_frequency={fft_frequency:.8g}")
    print(f"omega={omega:.8g}")
    print(f"sigma_hilbert={sigma_hilbert:.8g}")
    print(f"sigma_window_ratio={sigma_ratio:.8g}")
    print(f"late_amplitude={amplitude:.8g}")
    print(f"sigma_relerr_vs_linear={abs((sigma_hilbert - LINEAR_SIGMA) / LINEAR_SIGMA):.6g}")
    print(f"omega_relerr_vs_linear={omega_err:.6g}")


if __name__ == "__main__":
    main()
