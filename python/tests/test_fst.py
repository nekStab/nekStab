"""Physics + I/O regression tests for the nekStab FST inflow generator.

Run:  python -m pytest python/tests -q     (from the repo root)

These guard the pure-Python port against regressions: spectral accuracy of the
Chebyshev operators, the Blasius wall shear (hardcoded constant from the MATLAB),
the wavenumber-shell geometry, the OSS mode physics (Newton far-field BC, unit
energy, divergence-free reconstruction), the consolidated-file round-trip, and
parallel/serial reproducibility.
"""
import numpy as np
import pytest

# np.trapz removed in NumPy 2.0 (renamed np.trapezoid); support both.
try:
    from numpy import trapezoid as _trapz
except ImportError:  # NumPy < 2.0
    from numpy import trapz as _trapz

from nekstab.fst.cheb import cheb_collocation
from nekstab.fst.blasius import blasius_profile, C_WALL
from nekstab.fst.wavenumber import generate_wavenumbers
from nekstab.fst.oss import solve_oss_mode
from nekstab.fst.generator import (generate_fst_inflow, write_inflow_file,
                                   read_inflow_binary)


# ----------------------------------------------------------------- cheb
def test_cheb_spectral_differentiation():
    res = cheb_collocation(48, 0.3, 0.1)
    D1, D2, y = res["D1"], res["D2"], res["y"]
    k = np.pi
    f = np.sin(k * y)
    interior = slice(1, -1)
    assert np.max(np.abs((D1 @ f - k * np.cos(k * y))[interior])) < 1e-6
    assert np.max(np.abs((D2 @ f + k**2 * np.sin(k * y))[interior])) < 1e-5


# -------------------------------------------------------------- blasius
def test_blasius_wall_shear_matches_matlab_constant():
    y = np.linspace(0.0, 20.0, 200)
    _U, Up, _Upp = blasius_profile(20.0, y)
    # f''(0) must equal the hardcoded constant the MATLAB actually runs with.
    assert abs(Up[0] - C_WALL) < 1e-9


def test_blasius_reaches_freestream():
    y = np.linspace(0.0, 20.0, 200)
    U, _Up, _Upp = blasius_profile(20.0, y)
    assert abs(U[-1] - 1.0) < 1e-3
    assert np.all(np.diff(U) >= -1e-12)  # monotone increasing


# ----------------------------------------------------------- wavenumber
def test_wavenumber_shape_signs_radii():
    w = generate_wavenumbers(numk=5, kkini=0.5, kkfin=5.0, seed=42)
    assert w.shape == (50, 3)
    assert np.all(w[:, 0] > 0) and np.all(w[:, 1] > 0)  # omega>0, gamma>0
    r = np.sqrt((w**2).sum(axis=1))
    kk = np.linspace(0.5, 5.0, 5)
    for i in range(5):
        assert np.max(np.abs(r[i * 10:(i + 1) * 10] - kk[i])) < 1e-9


def test_wavenumber_reproducible():
    a = generate_wavenumbers(5, 0.5, 5.0, seed=7)
    b = generate_wavenumbers(5, 0.5, 5.0, seed=7)
    assert np.allclose(a, b)


# ---------------------------------------------------------------- oss
@pytest.mark.parametrize("omega,gamma,beta", [
    (0.8, 0.5, 0.6),
    (1.2, 0.3, 0.9),
    (0.5, 0.7, 0.4),
])
def test_oss_mode_physics(omega, gamma, beta):
    Re, Ny, Ly = 495.0, 35 * 6, 20.0
    rng = np.random.default_rng(0)
    y, U, V, W, diag = solve_oss_mode(omega, gamma, beta, Re, Ny, Ly, rng)
    assert np.all(np.isfinite(U)) and np.all(np.isfinite(V)) and np.all(np.isfinite(W))
    assert diag["newton_res"] < 1e-8            # far-field BC satisfied
    assert diag["spectral_div"] < 1e-9          # raw reconstruction div-free
    # unit free-stream energy
    from scipy.interpolate import CubicSpline
    ydm = Ly - 0.2 * Ly
    y1 = np.linspace(5.0, ydm, Ny)
    e = (np.abs(CubicSpline(y, U)(y1))**2 + np.abs(CubicSpline(y, V)(y1))**2
         + np.abs(CubicSpline(y, W)(y1))**2)
    e_norm = 0.5 * abs(_trapz(e, y1)) / abs(y1[0] - y1[-1])
    assert abs(e_norm - 1.0) < 1e-6


# --------------------------------------------------------- generator/IO
def test_parallel_equals_serial():
    r1 = generate_fst_inflow(495, 60, 20, 2, 0.8, 8, seed=0, jobs=1, verbose=False)
    r2 = generate_fst_inflow(495, 60, 20, 2, 0.8, 8, seed=0, jobs=2, verbose=False)
    for a, b in zip(r1["modes"], r2["modes"]):
        assert np.allclose(a["U"], b["U"])
        assert np.allclose(a["V"], b["V"])
        assert np.allclose(a["W"], b["W"])


def test_binary_roundtrip(tmp_path):
    r = generate_fst_inflow(495, 60, 20, 2, 0.8, 8, seed=0, jobs=1, verbose=False)
    p = tmp_path / "fst.bin"
    write_inflow_file(str(p), r, fmt="bin")
    numk, nmodes, Ny, Re, okini, okfin, length, tu, y, modes = read_inflow_binary(str(p))
    assert numk == r["numk"]
    assert nmodes == r["nmodes"]
    assert Ny == len(r["y"])
    assert Re == r["Re"]
    assert okini == r["okini"]
    assert okfin == r["okfin"]
    assert abs(length - r["length"]) < 1e-12
    assert abs(tu - r["tu"]) < 1e-12
    assert np.allclose(y, r["y"])
    om, ga, be, block = modes[0]
    assert np.allclose(block[:, 0], r["modes"][0]["U"].real)
    assert np.allclose(block[:, 3], r["modes"][0]["V"].imag)


def test_formatted_writes_y_grid_once(tmp_path):
    r = generate_fst_inflow(495, 60, 20, 2, 0.8, 8, seed=0, jobs=1, verbose=False)
    p = tmp_path / "fst.dat"
    write_inflow_file(str(p), r, fmt="dat")
    Ny = len(r["y"])
    # exactly one non-comment line should hold the full Ny-value y-grid
    lines = [ln for ln in p.read_text().splitlines()
             if ln and not ln.startswith("#")]
    y_lines = [ln for ln in lines if len(ln.split()) == Ny]
    assert len(y_lines) == 1
