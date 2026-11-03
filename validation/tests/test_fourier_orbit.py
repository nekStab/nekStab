"""Fourier series of one scalar record (src/fourier_series.f90).

A closed sum of harmonics returns its coefficients and rebuilds the
signal. The same sum, one sample at a time, matches an FFT of the
record. The mode count follows the doubling rule in the source.
Evaluation at a time between samples uses the period. A jump decays
as 1/k and a kink as 1/k^2.
Needs gfortran. The tests skip when it is missing.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"

REL_COEFF = 1e-13
ABS_EVAL = 1e-12
DECAY_TOL = 0.1


@pytest.fixture(scope="module")
def driver(build_driver) -> Path:
    return build_driver(
        "fourier_series",
        [ROOT / "src" / "fourier_series.f90", FORTRAN / "fourier_series_driver.f90"],
    )


def run_driver(exe: Path, text: str) -> list[str]:
    run = subprocess.run([str(exe)], input=text, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    return run.stdout.splitlines()


def fmt(values) -> str:
    return " ".join(repr(float(v)) for v in values)


def accum(exe: Path, samples: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    text = f"accum\n{len(samples)} {k}\n{fmt(samples)}\n"
    rows = run_driver(exe, text)
    assert len(rows) == k + 1
    cre = np.empty(k + 1)
    cim = np.empty(k + 1)
    for i, row in enumerate(rows):
        cre[i], cim[i] = (float(x) for x in row.split())
    return cre, cim


def evaluate(exe: Path, cre: np.ndarray, cim: np.ndarray, period: float, times: np.ndarray) -> np.ndarray:
    k = len(cre) - 1
    text = f"eval\n{k} {period!r} {len(times)}\n{fmt(cre)}\n{fmt(cim)}\n{fmt(times)}\n"
    rows = run_driver(exe, text)
    assert len(rows) == len(times)
    return np.array([float(row) for row in rows])


def select(exe: Path, ek: np.ndarray, kmax_fixed: int):
    kacc = len(ek) - 1
    text = f"select\n{kacc} {kmax_fixed}\n{fmt(ek)}\n"
    rows = run_driver(exe, text)
    k_sel, missing, ntrace = rows[0].split()
    traces = []
    for row in rows[1:]:
        a, b, c = row.split()
        traces.append((int(a), int(b), float(c)))
    assert len(traces) == int(ntrace)
    return int(k_sel), float(missing), traces


def kaccum(exe: Path, nsamples: int, kmax_fixed: int) -> int:
    rows = run_driver(exe, f"kaccum\n{nsamples} {kmax_fixed}\n")
    assert len(rows) == 1
    return int(rows[0])


def harmonic_samples(n: int, period: float, a0: float, amp_cos: list[float], amp_sin: list[float]) -> np.ndarray:
    """Samples at t = dt, 2 dt, ..., N dt. N dt = period."""
    step = np.arange(1, n + 1)
    t = step * period / n
    signal = np.full(n, a0, dtype=float)
    for k, (ak, bk) in enumerate(zip(amp_cos, amp_sin), start=1):
        signal += ak * np.cos(2.0 * np.pi * k * t / period) + bk * np.sin(2.0 * np.pi * k * t / period)
    return signal


def assert_coeff(got: float, expected: float) -> None:
    if expected == 0.0:
        assert abs(got) <= REL_COEFF
    else:
        assert abs(got - expected) <= REL_COEFF * abs(expected)


def fft_reference(samples: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """FFT of the record, with the last sample standing for time 0."""
    ordered = np.empty(len(samples))
    ordered[0] = samples[-1]
    ordered[1:] = samples[:-1]
    coeff = np.fft.fft(ordered) / len(samples)
    return coeff.real[: k + 1].copy(), coeff.imag[: k + 1].copy()


def test_three_harmonics_match_and_rebuild(driver):
    period = 5.0
    n = 64
    a0 = 0.3
    amp_cos = [1.25, -0.4, 0.15]
    amp_sin = [-0.7, 0.2, 0.85]
    samples = harmonic_samples(n, period, a0, amp_cos, amp_sin)
    cre, cim = accum(driver, samples, 8)

    assert_coeff(cre[0], a0)
    assert_coeff(cim[0], 0.0)
    for k in range(1, 4):
        assert_coeff(cre[k], amp_cos[k - 1] / 2.0)
        assert_coeff(cim[k], -amp_sin[k - 1] / 2.0)
    for k in range(4, 9):
        assert_coeff(cre[k], 0.0)
        assert_coeff(cim[k], 0.0)

    times = np.arange(1, n + 1) * period / n
    rebuilt = evaluate(driver, cre[:4], cim[:4], period, times)
    assert np.max(np.abs(rebuilt - samples)) <= ABS_EVAL

    t = 0.37 * period
    exact_off = a0 + sum(
        ak * np.cos(2.0 * np.pi * k * t / period) + bk * np.sin(2.0 * np.pi * k * t / period)
        for k, (ak, bk) in enumerate(zip(amp_cos, amp_sin), start=1)
    )
    got = evaluate(driver, cre[:4], cim[:4], period, np.array([t]))
    assert abs(got[0] - exact_off) <= ABS_EVAL


def test_one_sample_sum_matches_fft(driver):
    rng = np.random.default_rng(1)
    samples = rng.normal(size=64)
    k = (len(samples) - 1) // 2
    cre, cim = accum(driver, samples, k)
    ref_re, ref_im = fft_reference(samples, k)
    assert np.max(np.abs(cre - ref_re)) <= ABS_EVAL
    assert np.max(np.abs(cim - ref_im)) <= ABS_EVAL


def test_mode_count_follows_the_doubling_rule(driver):
    # The count never includes the Nyquist mode. One step stores no mode.
    assert kaccum(driver, 1, 0) == 0
    k_sel, missing, traces = select(driver, np.array([1.0]), 0)
    assert k_sel == 0
    assert missing == 0.0
    assert traces == []
    assert kaccum(driver, 64, 0) == 31
    assert kaccum(driver, 64, 2) == 2
    assert kaccum(driver, 64, 1000) == 31
    assert kaccum(driver, 16, 0) == 7
    assert kaccum(driver, 5, 0) == 2

    # Two harmonics on 64 samples. Modes above 2 are empty, but the
    # rule starts at 8 and does not drop a mode. K stays 8.
    samples = harmonic_samples(64, 1.0, 0.2, [1.0, -0.3], [0.5, 0.4])
    cre, cim = accum(driver, samples, 31)
    ek = np.zeros(32)
    ek[0] = cre[0] ** 2
    ek[1:] = cre[1:] ** 2 + cim[1:] ** 2
    k_sel, missing, traces = select(driver, ek, 0)
    assert k_sel == 8
    assert missing < 1e-10
    assert traces == []

    # Five samples hold only modes 0..2. The same signal then returns K = 2.
    samples = harmonic_samples(5, 1.0, 0.2, [1.0, -0.3], [0.5, 0.4])
    cre, cim = accum(driver, samples, 2)
    ek = np.zeros(3)
    ek[0] = cre[0] ** 2
    ek[1:] = cre[1:] ** 2 + cim[1:] ** 2
    k_sel, missing, traces = select(driver, ek, 0)
    assert k_sel == 2
    assert missing < 1e-10
    assert traces == []

    # Flat fluctuation energy. 16 keeps 0.75, and 0.75 > 0.1 * 0.875,
    # so the count stays at 8. The mean in ek[0] does not enter.
    ek = np.ones(65)
    ek[0] = 1.0e6
    k_sel, missing, traces = select(driver, ek, 0)
    assert k_sel == 8
    assert abs(missing - 0.875) <= 1e-14
    assert len(traces) == 1
    assert traces[0][0] == 8
    assert traces[0][1] == 16
    assert abs(traces[0][2] - 0.75) <= 1e-14

    # Modes 1..10 hold equal energy. The step 8 -> 16 captures the
    # rest, so K becomes 16.
    ek = np.zeros(33)
    ek[1:11] = 1.0
    k_sel, missing, traces = select(driver, ek, 0)
    assert k_sel == 16
    assert abs(missing) <= 1e-14
    assert len(traces) == 1
    assert (traces[0][0], traces[0][1]) == (8, 16)
    assert abs(traces[0][2]) <= 1e-14

    # A fixed cap keeps every accumulated mode. No doubling is printed.
    ek = np.ones(33)
    ek[0] = 0.0
    k_sel, missing, traces = select(driver, ek, 32)
    assert k_sel == 32
    assert abs(missing) <= 1e-14
    assert traces == []


def test_off_grid_time_uses_the_period(driver):
    period = 5.0
    n = 64
    a0 = 0.3
    amp_cos = [1.25, -0.4, 0.15]
    amp_sin = [-0.7, 0.2, 0.85]
    samples = harmonic_samples(n, period, a0, amp_cos, amp_sin)
    cre, cim = accum(driver, samples, 3)
    cre, cim = cre[:4], cim[:4]

    def exact(t: float) -> float:
        return a0 + sum(
            ak * np.cos(2.0 * np.pi * k * t / period) + bk * np.sin(2.0 * np.pi * k * t / period)
            for k, (ak, bk) in enumerate(zip(amp_cos, amp_sin), start=1)
        )

    # Halfway between sample 10 and sample 11. A line between those
    # samples misses the signal; the series does not.
    t_mid = (10.5) * period / n
    got = evaluate(driver, cre, cim, period, np.array([t_mid, t_mid + period, t_mid - period]))
    assert abs(got[0] - exact(t_mid)) <= ABS_EVAL
    assert abs(got[1] - exact(t_mid)) <= ABS_EVAL
    assert abs(got[2] - exact(t_mid)) <= ABS_EVAL
    linear = 0.5 * (samples[9] + samples[10])
    assert abs(linear - exact(t_mid)) > 1e-4
    assert abs(got[0] - linear) > 1e-4


def decay_slope(mag: np.ndarray) -> float:
    modes = np.arange(4, 33)
    return float(np.polyfit(np.log(modes), np.log(mag[modes]), 1)[0])


def test_jump_and_kink_decay(driver):
    n = 512
    step = np.arange(1, n + 1)
    saw = step / n - 0.5
    kink = (step / n - 0.5) ** 2
    cre, cim = accum(driver, saw, 40)
    saw_mag = np.sqrt(cre**2 + cim**2)
    saw_slope = decay_slope(saw_mag)
    # The jump is a sine series. |c_4| / |c_32| is 8 for a 1/k decay.
    # The cosine part of this record is the jump leakage, about 1/(2N).
    saw_ratio = abs(cim[4]) / abs(cim[32])
    cre, cim = accum(driver, kink, 40)
    kink_mag = np.sqrt(cre**2 + cim**2)
    kink_slope = decay_slope(kink_mag)
    # The kink is a cosine series. |c_4| / |c_32| is 64 for a 1/k^2 decay.
    kink_ratio = abs(cre[4]) / abs(cre[32])
    assert abs(saw_slope - (-1.0)) <= DECAY_TOL
    assert abs(kink_slope - (-2.0)) <= DECAY_TOL
    assert abs(saw_ratio - 8.0) <= 0.5
    assert abs(kink_ratio - 64.0) <= 4.0
