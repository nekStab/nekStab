#!/usr/bin/env python3
"""Measure the period of the Re = 1950 limit cycle of the cubic cavity and plot the probe signal.

The script reads the probe file cav.his of this stage (vy at the cavity centre, written every
time step). The oscillation amplitude overshoots by 2 % near t = 1100 and settles by t = 1300,
so every estimate uses the window FIT_START to the end of the run. The period comes from a
least-squares fit of a mean and three harmonics with a free frequency. The script also gives the
FFT peak and the zero-crossing count of the same window. It writes period_fit.dat and
plot_period.png. Times are in units of L/U.

Usage: python plot.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

CASE_DIR = Path(__file__).resolve().parent
FIT_START = 1400.0  # the amplitude has settled by t = 1300; the window starts later
START_UP = 50.0  # the seed flow of Re = 1914 adjusts to Re = 1950 first; panel (a) starts after it
HARMONICS = 3


def read_probe(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return the time and the vy signal of the first probe of a Nek5000 .his file."""
    raw = np.loadtxt(path, skiprows=2)
    return raw[:, 0], raw[:, 2]


def harmonic_design(t: np.ndarray, frequency: float) -> np.ndarray:
    """Return the matrix of a mean and cosines and sines of the first HARMONICS multiples."""
    columns = [np.ones_like(t)]
    for k in range(1, HARMONICS + 1):
        columns.append(np.cos(2.0 * np.pi * k * frequency * t))
        columns.append(np.sin(2.0 * np.pi * k * frequency * t))
    return np.column_stack(columns)


def fit_residual(frequency: float, t: np.ndarray, v: np.ndarray) -> float:
    """Return the rms residual of the linear least-squares fit at a given frequency, in the units of v."""
    design = harmonic_design(t, frequency)
    coefficients, *_ = np.linalg.lstsq(design, v, rcond=None)
    return float(np.sqrt(np.mean((design @ coefficients - v) ** 2)))


def fit_period(t: np.ndarray, v: np.ndarray, guess: float) -> tuple[float, float]:
    """Return the period and the rms residual of the harmonic fit within 1 % of the guessed period."""
    found = minimize_scalar(
        fit_residual,
        args=(t, v),
        bounds=(0.99 / guess, 1.01 / guess),
        method="bounded",
        options={"xatol": 1e-10},
    )
    return 1.0 / found.x, float(found.fun)


def fft_peak_period(t: np.ndarray, v: np.ndarray) -> float:
    """Return the period of the FFT peak of the uniformly resampled signal (parabolic refinement)."""
    dt = float(np.median(np.diff(t)))
    uniform = np.arange(t[0], t[-1], dt)
    signal = np.interp(uniform, t, v)
    signal -= signal.mean()
    n_fft = 1 << int(np.ceil(np.log2(len(signal) * 4)))
    power = np.abs(np.fft.rfft(signal * np.hanning(len(signal)), n=n_fft)) ** 2
    frequencies = np.fft.rfftfreq(n_fft, d=dt)
    k = int(np.argmax(power[1:])) + 1
    a, b, c = power[k - 1], power[k], power[k + 1]
    shift = 0.5 * (a - c) / (a - 2.0 * b + c)
    return 1.0 / ((k + shift) * (frequencies[1] - frequencies[0]))


def crossing_period(t: np.ndarray, v: np.ndarray) -> float:
    """Return the mean time between upward zero crossings of the signal without its mean."""
    dt = float(np.median(np.diff(t)))
    uniform = np.arange(t[0], t[-1], dt)
    signal = np.interp(uniform, t, v)
    signal -= signal.mean()
    negative = np.signbit(signal)
    up = np.where((~negative[1:]) & negative[:-1])[0]
    times = uniform[up] - signal[up] * (uniform[up + 1] - uniform[up]) / (
        signal[up + 1] - signal[up]
    )
    return float((times[-1] - times[0]) / (len(times) - 1))


def period_envelope(
    t: np.ndarray, v: np.ndarray, width: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the centre, the maximum and the minimum of the signal in consecutive complete windows of one width."""
    edges = np.arange(t[0], t[-1], width)
    bounds = np.searchsorted(t, edges)
    centres, upper, lower = [], [], []
    for first, last in zip(bounds[:-1], bounds[1:], strict=True):
        if last - first > 10:
            centres.append(t[first:last].mean())
            upper.append(v[first:last].max())
            lower.append(v[first:last].min())
    return np.array(centres), np.array(upper), np.array(lower)


def main() -> None:
    """Fit the period, write period_fit.dat and draw the probe signal and the fit."""
    t_all, v_all = read_probe(CASE_DIR / "cav.his")
    end_time = float(t_all[-1])
    window = t_all >= FIT_START
    t_fit, v_fit = t_all[window], v_all[window]
    period_fft = fft_peak_period(t_fit, v_fit)
    period_crossing = crossing_period(t_fit, v_fit)
    period, rms = fit_period(t_fit, v_fit, period_fft)
    design_fit = harmonic_design(t_fit, 1.0 / period)
    coefficients, *_ = np.linalg.lstsq(design_fit, v_fit, rcond=None)
    mean = float(coefficients[0])
    amplitude = float(0.5 * (v_fit.max() - v_fit.min()))
    print(f"end time {end_time:.1f}; window {FIT_START:.0f} to {end_time:.0f}")
    print(
        f"period: harmonic fit {period:.5f}, FFT peak {period_fft:.5f}, zero crossings {period_crossing:.5f}"
    )
    print(
        f"St = {1.0 / period:.5f}; fit rms residual {rms:.2e}; probe mean {mean:.6f}; half-range {amplitude:.3e}"
    )

    header = "period_fit  St_fit  period_fft  period_crossings  fit_rms  probe_mean  probe_half_range  end_time"
    np.savetxt(
        CASE_DIR / "period_fit.dat",
        [
            [
                period,
                1.0 / period,
                period_fft,
                period_crossing,
                rms,
                mean,
                amplitude,
                end_time,
            ]
        ],
        header=header,
        fmt="%.8g",
    )

    centres, upper, lower = period_envelope(t_all, v_all, period)
    after_start_up = centres > START_UP
    centres, upper, lower = (a[after_start_up] for a in (centres, upper, lower))
    t_zoom = (t_all >= end_time - 5.0 * period) & (t_all <= end_time)
    design_zoom = harmonic_design(t_all[t_zoom], 1.0 / period)

    fig, (ax_all, ax_zoom) = plt.subplots(
        1, 2, figsize=(7.0, 3.0), layout="constrained"
    )
    ax_all.plot(
        centres,
        upper,
        "o-",
        c="k",
        lw=0.8,
        ms=2.5,
        markevery=3,
        label="maximum per period",
    )
    ax_all.plot(
        centres,
        lower,
        "s-",
        c="0.5",
        lw=0.8,
        ms=2.5,
        markevery=3,
        label="minimum per period",
    )
    ax_all.axvline(FIT_START, c="0.5", lw=0.8, ls="--", label="fit window starts")
    ax_all.set_xlabel(r"time $t$ ($L/U$)")
    ax_all.set_ylabel(r"$v_y$ at the cavity centre ($U$)")
    ax_all.set_title(
        "(a) Envelope of the probe signal, Re = 1950", loc="left", fontsize=8
    )
    ax_all.legend(loc="lower left", fontsize=6.5)

    ax_zoom.plot(
        t_all[t_zoom],
        v_all[t_zoom],
        "o-",
        c="0.6",
        lw=2.0,
        ms=3.0,
        markevery=max(1, int(t_zoom.sum() // 60)),
        label="DNS",
    )
    ax_zoom.plot(
        t_all[t_zoom],
        design_zoom @ coefficients,
        "s-",
        c="k",
        lw=0.8,
        ms=2.5,
        markevery=max(1, int(t_zoom.sum() // 45)),
        label=f"fit, T = {period:.5f}",
    )
    ax_zoom.set_xlabel(r"time $t$ ($L/U$)")
    ax_zoom.set_ylabel(r"$v_y$ ($U$)")
    ax_zoom.set_title("(b) Last five periods", loc="left", fontsize=8)
    low, high = ax_zoom.get_ylim()
    ax_zoom.set_ylim(low, high + 0.3 * (high - low))  # headroom for the legend
    ax_zoom.legend(loc="upper center", ncols=2, fontsize=7)

    output = CASE_DIR / "plot_period.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
