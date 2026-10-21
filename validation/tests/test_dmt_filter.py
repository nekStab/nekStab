"""Oscillator-pair filter of the DMT bank (src/dmt_filter.f90).

The module is compiled alone (it needs no Nek5000) and driven by
fortran/dmt_filter_driver.f90. The references are computed here with numpy only:
the step coefficients against a fourth-order Runge-Kutta integration of the
oscillator with a linearly varying input, the amplitude coefficients against the
formula in the file header, and the filter on a tone with a large mean against
the amplitude of the tone.
Needs gfortran; the tests skip when it is missing.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"


@pytest.fixture(scope="module")
def driver(build_driver) -> Path:
    return build_driver(
        "dmt_filter",
        [ROOT / "src" / "dmt_filter.f90", FORTRAN / "dmt_filter_driver.f90"],
    )


def ask(exe: Path, *lines: str) -> list[list[float]]:
    run = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True
    )
    assert run.returncode == 0, run.stdout + run.stderr
    return [[float(x) for x in row.split()] for row in run.stdout.strip().split("\n")]


def rk4_step(omega: float, beta: float, h: float, z0, u0: float, u1: float, n: int = 4000):
    """Integrate p' = q, q' = u - beta*q - omega^2*p over h with u linear from u0 to u1."""
    z = np.array(z0, dtype=float)
    dt = h / n

    def f(t, z):
        u = u0 + (u1 - u0) * t / h
        return np.array([z[1], u - beta * z[1] - omega**2 * z[0]])

    t = 0.0
    for _ in range(n):
        k1 = f(t, z)
        k2 = f(t + dt / 2, z + dt / 2 * k1)
        k3 = f(t + dt / 2, z + dt / 2 * k2)
        k4 = f(t + dt, z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    return z


CASES = [
    (1.032, 0.0516, 1.47e-3),
    (1.032, 0.0516, 0.5),
    (3.0, 0.3, 2.0),
    (0.2, 0.01, 1.0e-4),
    (0.0, 0.05, 0.3),
    (0.0, 0.05, 5.0),
]


@pytest.mark.parametrize("omega,beta,h", CASES)
def test_step_coefficients_match_a_runge_kutta_integration(driver, omega, beta, h):
    (row,) = ask(driver, f"coef {omega!r} {beta!r} {h!r}")
    phi = np.array([[row[0], row[1]], [row[2], row[3]]])
    a = np.array(row[4:6])
    b = np.array(row[6:8])
    rng = np.random.default_rng(7)
    for _ in range(3):
        z0 = rng.normal(size=2)
        u0, u1 = rng.normal(size=2)
        got = phi @ z0 + a * u0 + b * u1
        want = rk4_step(omega, beta, h, z0, u0, u1)
        assert np.allclose(got, want, rtol=1e-10, atol=1e-12), (got, want)


def test_a_constant_input_leaves_the_equilibrium_unchanged(driver):
    omega, beta, h = 1.032, 0.0516, 0.01
    (row,) = ask(driver, f"coef {omega!r} {beta!r} {h!r}")
    phi = np.array([[row[0], row[1]], [row[2], row[3]]])
    a, b = np.array(row[4:6]), np.array(row[6:8])
    z = np.array([2.5 / omega**2, 0.0])
    assert np.allclose(phi @ z + (a + b) * 2.5, z, rtol=1e-13, atol=1e-14)


@pytest.mark.parametrize("omega,beta", [(1.032, 0.0516), (0.5, 0.2), (0.0, 0.05)])
@pytest.mark.parametrize("t", [0.0, 0.7, 123.456])
def test_amplitude_coefficients_follow_the_header_formula(driver, omega, beta, t):
    (row,) = ask(driver, f"amp {omega!r} {beta!r} {t!r}")
    cre, cim = np.array(row[:3]), np.array(row[3:])
    rng = np.random.default_rng(3)
    q, p, u = rng.normal(size=3)
    qd = u - beta * q - omega**2 * p
    if omega == 0.0:
        re_want, im_want = beta * q, 0.0
    else:
        s, c = math.sin(omega * t), math.cos(omega * t)
        re_want = beta * (c * q - s * qd / omega)
        im_want = -beta * (s * q + c * qd / omega)
    assert math.isclose(cre @ np.array([q, p, u]), re_want, rel_tol=1e-12, abs_tol=1e-14)
    assert math.isclose(cim @ np.array([q, p, u]), im_want, rel_tol=1e-12, abs_tol=1e-14)


@pytest.mark.parametrize("beta_over_omega", [0.05, 0.1])
def test_the_amplitude_of_a_tone_on_a_large_mean(driver, beta_over_omega):
    # A mean of 1 under a tone of 0.01. The estimate beta*(q + i*omega*p) passes the mean with
    # the gain i*beta/omega and has an error of 5 in units of the tone amplitude; the
    # derivative form has none, only the start-up transient exp(-beta*t/2).
    omega = 1.032
    beta = beta_over_omega * omega
    period = 2 * math.pi / omega
    h = period / 4133
    nper = 40
    n = round(nper * period / h)
    phase, amp = 0.7, 0.01
    (row,) = ask(driver, f"tone {omega!r} {beta!r} {h!r} {n} 1.0 {amp!r} {phase!r}")
    got = complex(row[0], row[1])
    want = amp * complex(math.cos(phase), math.sin(phase))
    transient = math.exp(-beta * row[2] / 2)
    assert abs(got - want) / amp < 3 * transient + 1e-6


def test_the_mean_filter_follows_the_mean(driver):
    # omega = 0: u = 3 + 0.5*cos(0*t) is constant 3.5, the exponential mean converges to it.
    (row,) = ask(driver, "tone 0.0 0.05 0.05 4000 3.0 0.5 0.0")
    assert math.isclose(row[0], 3.5, rel_tol=1e-4)
    assert row[1] == 0.0
