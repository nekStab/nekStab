"""Period estimate from upward crossings (src/zero_cross.f90).

A sine of period T, sampled at a fixed dt, returns T. The tolerance is
the interpolation bound for that sine, not a fitted number. Dropping the
linear interpolation leaves an error of order dt, above the bound.
Needs gfortran. The tests skip when it is missing.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"

# dt does not divide T, so a crossing falls between two samples.
SINE_CASES = [
    (1.0, 0.017, 400),
    (1.0, 0.023, 280),
    (1.0, 0.031, 220),
]


def period_bound(period: float, dt: float) -> float:
    """Bound on |T_hat - T| for sin(2*pi*t/T).

    On a bracket of width dt about a zero, |f''| <= omega**3*dt.
    Linear interpolation then misses each crossing by at most
    omega**2*dt**3 / (8*cos(omega*dt)). A period has two crossings.
    """
    omega = 2.0 * math.pi / period
    return omega**2 * dt**3 / (4.0 * math.cos(omega * dt))


def mean_bound(period: float, dt: float) -> float:
    """Bound on the trapezoid error of the mean of one period of a sine."""
    omega = 2.0 * math.pi / period
    return omega**2 * dt**2 / 12.0


@pytest.fixture(scope="module")
def driver(build_driver) -> Path:
    return build_driver(
        "zero_cross",
        [ROOT / "src" / "zero_cross.f90", FORTRAN / "zero_cross_driver.f90"],
    )


def run_sine(exe: Path, n: int, dt: float, period: float, offset: float, v_init: float):
    text = f"sine\n{n} {dt!r} {period!r} {offset!r} {v_init!r}\n"
    run = subprocess.run([str(exe)], input=text, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    rows = []
    done = None
    for line in run.stdout.splitlines():
        parts = line.split()
        if parts[0] == "done":
            done = (int(parts[1]), float(parts[2]))
        else:
            rows.append(tuple(float(x) for x in parts))
    assert done is not None
    return rows, done


@pytest.mark.parametrize("period,dt,n", SINE_CASES)
def test_sine_period_matches_to_interpolation_order(driver, period, dt, n):
    bound = period_bound(period, dt)
    assert bound < 0.25 * dt
    rows, (n_cross, _) = run_sine(driver, n, dt, period, 0.0, 0.0)
    t_last = (n - 1) * dt
    n_period = math.floor(t_last / period) - 1
    assert n_period >= 3
    assert len(rows) == n_period
    assert n_cross == n_period + 1
    for estimated, v_used, v_new in rows:
        assert abs(estimated - period) <= bound
        assert abs(v_used) <= mean_bound(period, dt)
        assert abs(v_new) <= mean_bound(period, dt)


def test_mean_over_the_first_period_removes_the_offset(driver):
    period, dt, n = 1.0, 0.017, 160
    offset = -0.25
    rows, _ = run_sine(driver, n, dt, period, offset, 0.0)
    assert rows
    estimated, v_used, v_new = rows[0]
    # The first period still uses the initial level, so its crossings
    # are not the zeros of the sine. The generic interpolation bound applies.
    omega = 2.0 * math.pi / period
    slope = math.sqrt(1.0 - offset**2)
    bound = omega * dt**2 / (4.0 * slope)
    assert abs(estimated - period) <= bound
    assert v_used == 0.0
    assert abs(v_new - offset) <= mean_bound(period, dt)
