#!/usr/bin/env python3
"""Plot the probe signal of the DNS on the limit cycle of the cubic cavity at Re = 1950.

The stage starts from the on-orbit field of ../000_dns_period at t = 1700 and runs to t = 1800. It writes 200
snapshots for the POD and DMD stages and the signal of one probe at (0.25, 0.25, 0), cav.his: time, v_x, v_y, v_z
and the pressure, one line per time step. The script checks that the run stays on the cycle:

- the period of a least-squares fit of a mean and three harmonics to v_y, with a free frequency, against the
  period of the DNS of ../000_dns_period, 10.78909;
- the amplitude of v_y in each complete period, which does not grow or fall when the start field is on the orbit.

It writes cycle_check.dat with one row: the fitted period, the rms residual of the fit, the mean of v_y, the half
range of v_y, the half range in the last complete period over the half range in the first, and the first and last
time of the signal. It writes plot_probe.png.

Lengths are in units of the cavity edge L, speeds in units of the lid speed U, time in L/U.

Usage: python plot.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
from scipy.optimize import minimize_scalar

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
CYCLE_PERIOD = 10.78909  # 000_dns_period
HARMONICS = 3


def read_probe(path: Path) -> np.ndarray:
    """Return the columns time, v_x, v_y of the first probe of a Nek5000 .his file."""
    raw = np.loadtxt(path, skiprows=2)
    return raw[:, :3]


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


def period_half_range(
    t: np.ndarray, v: np.ndarray, width: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return the centre time and the half range of v in each complete window of one width from the first sample."""
    centres, half_ranges = [], []
    start = t[0]
    while start + width <= t[-1]:
        inside = (t >= start) & (t < start + width)
        centres.append(start + 0.5 * width)
        half_ranges.append(0.5 * (v[inside].max() - v[inside].min()))
        start += width
    return np.array(centres), np.array(half_ranges)


def main() -> None:
    """Fit the period, write cycle_check.dat and draw the probe signal."""
    time, v_x, v_y = read_probe(CASE_DIR / "cav.his").T
    period, rms = fit_period(time, v_y, CYCLE_PERIOD)
    design = harmonic_design(time, 1.0 / period)
    coefficients, *_ = np.linalg.lstsq(design, v_y, rcond=None)
    mean = float(coefficients[0])
    half_range = float(0.5 * (v_y.max() - v_y.min()))
    centres, per_period = period_half_range(time, v_y, period)
    ratio = float(per_period[-1] / per_period[0])
    np.savetxt(
        CASE_DIR / "cycle_check.dat",
        np.array([[period, rms, mean, half_range, ratio, time[0], time[-1]]]),
        header="period_fit  fit_rms  probe_mean  probe_half_range  half_range(last period)/half_range(first)  t_first  t_last",
        fmt="%.8g",
    )
    print(
        f"signal from t = {time[0]:.2f} to {time[-1]:.2f}; period {period:.5f} "
        f"({period / CYCLE_PERIOD - 1.0:+.2e} against {CYCLE_PERIOD}); fit rms {rms:.2e}"
    )
    print(
        f"probe mean {mean:.6f}, half range {half_range:.3e}; half range of the last period over the first {ratio:.4f}"
    )

    figure = plt.figure(figsize=(7.2, 2.8), layout="constrained")
    grid = figure.add_gridspec(1, 3, width_ratios=[1.5, 1.0, 1.0])

    axes = figure.add_subplot(grid[0, 0])
    every = max(1, len(time) // 60)
    axes.plot(
        time, v_y, c="0.6", lw=2.6, marker="o", ms=3, markevery=every, label="DNS"
    )
    axes.plot(
        time,
        design @ coefficients,
        c="k",
        lw=0.9,
        marker="s",
        ms=2.5,
        markevery=every,
        label=f"fit, $T$ = {period:.4f}",
    )
    axes.set_xlabel(r"$t$ ($L/U$)")
    axes.set_ylabel(r"$v_y$ at (0.25, 0.25, 0) ($U$)")
    axes.set_title("(a) Probe signal", loc="left", fontsize=7.5)
    low, high = axes.get_ylim()
    axes.set_ylim(low, high + 0.3 * (high - low))
    axes.legend(loc="upper center", ncols=2, fontsize=6.5)

    axes = figure.add_subplot(grid[0, 1])
    axes.plot(v_x, v_y, c="C3", lw=1.0, marker="o", ms=2, markevery=every)
    axes.locator_params(axis="x", nbins=4)
    axes.set_xlabel(r"$v_x$ ($U$)")
    axes.set_ylabel(r"$v_y$ ($U$)")
    axes.set_title("(b) Closed orbit", loc="left", fontsize=7.5)

    axes = figure.add_subplot(grid[0, 2])
    relative = per_period / per_period.mean()
    axes.plot(centres, relative, c="C0", lw=1.4, marker="o", ms=4)
    axes.set_ylim(min(0.99, relative.min() - 0.002), max(1.01, relative.max() + 0.002))
    axes.set_xticks(np.round(np.linspace(centres[0], centres[-1], 2), -1))
    axes.set_xlabel(r"$t$ ($L/U$)")
    axes.set_ylabel(r"half range / its mean")
    axes.set_title("(c) Amplitude per period", loc="left", fontsize=7.5)

    output = CASE_DIR / "plot_probe.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
