"""
test_parsers.py — unit tests for validation/parsers/ (residual + spectrum).

All test data is synthesised in-memory with tempfile; no real run outputs
required.  Only numpy and the standard library are used.
"""

from __future__ import annotations

import math
import tempfile
import os

import numpy as np
import pytest

from validation.parsers import (
    parse_residual_log,
    parse_force_history,
    detect_convergence,
    parse_spectrum_file,
    detect_unit_circle_modes,
    leading_modes,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tmp(content: str, suffix: str = ".dat") -> str:
    """Write *content* to a named temp file and return its path."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as fh:
        fh.write(content)
    return path


# ---------------------------------------------------------------------------
# parse_residual_log — boostconv / SFD style (3-col)
# ---------------------------------------------------------------------------

class TestParseResidualBoostconvStyle:
    """3-col: time, residual, rate — residual should be monotone decreasing."""

    def setup_method(self):
        # Build a synthetic boostconv-style residual file.
        # Row 0 has rate = Infinity (Nek5000 outputs this literal string).
        lines = ["  0.0000000E+00  1.0000000E+00       Infinity\n"]
        t = 0.05
        r = 1.0
        for _ in range(9):
            r *= 0.5  # strict monotone decrease
            rate = -r  # negative rate (convergence)
            lines.append(f"  {t:.7E}  {r:.7E}  {rate:.7E}\n")
            t += 0.05
        self.path = _write_tmp("".join(lines))

    def teardown_method(self):
        os.unlink(self.path)

    def test_shape(self):
        time, resid, rate = parse_residual_log(self.path)
        assert time.shape == (10,), f"expected (10,), got {time.shape}"
        assert resid.shape == (10,), f"expected (10,), got {resid.shape}"
        assert rate is not None, "rate should not be None for 3-col file"
        assert rate.shape == (10,)

    def test_monotone_decrease(self):
        _, resid, _ = parse_residual_log(self.path)
        diffs = np.diff(resid)
        assert np.all(diffs < 0), "residual should be strictly decreasing"

    def test_first_rate_is_inf(self):
        _, _, rate = parse_residual_log(self.path)
        assert np.isinf(rate[0]), "first rate row (Infinity) should map to np.inf"

    def test_time_monotone(self):
        time, _, _ = parse_residual_log(self.path)
        assert np.all(np.diff(time) > 0), "time should be strictly increasing"


# ---------------------------------------------------------------------------
# parse_residual_log — Newton style (8-col)
# ---------------------------------------------------------------------------

class TestParseResidualNewtonStyle:
    """8-col: step cumul iter inner time res_prev residual tol."""

    def setup_method(self):
        rows = []
        for i in range(1, 6):
            step = i
            cumul = i * 100
            it = 100
            inner = i * 10
            t = i * 0.5
            res_prev = 1e-3 / (i)
            res = 1e-3 / (i + 1)
            tol = 1e-9
            rows.append(
                f"  {step:6d}  {cumul:6d}  {it:6d}  {inner:6d}"
                f"  {t:.6E}  {res_prev:.6E}  {res:.6E}  {tol:.6E}\n"
            )
        self.path = _write_tmp("".join(rows))
        self.expected_times = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
        self.expected_resids = np.array(
            [1e-3 / 2, 1e-3 / 3, 1e-3 / 4, 1e-3 / 5, 1e-3 / 6]
        )

    def teardown_method(self):
        os.unlink(self.path)

    def test_time_from_col4(self):
        time, _, _ = parse_residual_log(self.path)
        np.testing.assert_allclose(
            time, self.expected_times, rtol=1e-6,
            err_msg="Newton: time should come from col4 (simulation time)"
        )

    def test_residual_from_col6(self):
        _, resid, _ = parse_residual_log(self.path)
        np.testing.assert_allclose(
            resid, self.expected_resids, rtol=1e-6,
            err_msg="Newton: residual should come from col6"
        )

    def test_rate_is_none(self):
        _, _, rate = parse_residual_log(self.path)
        assert rate is None, "Newton 8-col file: rate should be None"


# ---------------------------------------------------------------------------
# parse_residual_log — GMRES style (7-col)
# ---------------------------------------------------------------------------

class TestParseResidualGmresStyle:
    """7-col: step ? iter cumul res_prev residual tol."""

    def setup_method(self):
        rows = []
        for i in range(1, 5):
            step = i
            something = 1
            it = i * 5
            cumul = i * 5
            res_prev = 1e-3 / i
            res = 1e-3 / (i + 0.5)
            tol = 1e-9
            rows.append(
                f"  {step:6d}  {something:6d}  {it:6d}  {cumul:6d}"
                f"  {res_prev:.6E}  {res:.6E}  {tol:.6E}\n"
            )
        self.path = _write_tmp("".join(rows))

    def teardown_method(self):
        os.unlink(self.path)

    def test_shape(self):
        time, resid, rate = parse_residual_log(self.path)
        assert time.shape == (4,)
        assert resid.shape == (4,)
        assert rate is None

    def test_time_is_step(self):
        time, _, _ = parse_residual_log(self.path)
        np.testing.assert_array_equal(
            time, [1.0, 2.0, 3.0, 4.0],
            err_msg="GMRES: time should be step index from col0"
        )

    def test_residual_positive(self):
        _, resid, _ = parse_residual_log(self.path)
        assert np.all(resid > 0), "all residuals should be positive"


# ---------------------------------------------------------------------------
# detect_convergence
# ---------------------------------------------------------------------------

class TestDetectConvergenceConverged:
    def test_converged_true(self):
        n = 20
        resid = np.logspace(0, -8, n)  # decays from 1 to 1e-8
        time = np.linspace(0, 1, n)
        result = detect_convergence(resid, time, tol=1e-6)
        assert result["converged"] is True, "should detect convergence"
        assert isinstance(result["converge_time"], float), (
            "converge_time should be a float"
        )
        assert result["converge_time"] <= time[-1]

    def test_final_residual_value(self):
        resid = np.array([1.0, 0.1, 0.01, 1e-7])
        time = np.array([0.0, 1.0, 2.0, 3.0])
        result = detect_convergence(resid, time, tol=1e-6)
        assert math.isclose(result["final_residual"], 1e-7, rel_tol=1e-9)


class TestDetectConvergenceNotConverged:
    def test_not_converged(self):
        resid = np.array([1.0, 0.5, 0.1, 1e-5])  # above tol=1e-6
        time = np.array([0.0, 1.0, 2.0, 3.0])
        result = detect_convergence(resid, time, tol=1e-6)
        assert result["converged"] is False, "should not be converged"
        assert result["converge_time"] is None, (
            "converge_time should be None when residual never drops to tol"
        )

    def test_final_residual_above_tol(self):
        resid = np.array([0.5, 0.3, 0.01])
        time = np.array([0.0, 1.0, 2.0])
        result = detect_convergence(resid, time, tol=1e-6)
        assert result["final_residual"] == pytest.approx(0.01)


# ---------------------------------------------------------------------------
# parse_force_history
# ---------------------------------------------------------------------------

class TestParseForceHistory:
    """Nek5000 lift_drag.dat: step, time, col1..col8."""

    def setup_method(self):
        lines = []
        for i in range(5):
            cols = [i, i * 0.1] + [float(j + i) for j in range(8)]
            lines.append("  " + "  ".join(f"{c:.7E}" for c in cols) + "\n")
        self.path = _write_tmp("".join(lines))

    def teardown_method(self):
        os.unlink(self.path)

    def test_keys_present(self):
        result = parse_force_history(self.path)
        for key in ("step", "time", "col_2"):
            assert key in result, f"expected key {key!r} in result"

    def test_step_shape(self):
        result = parse_force_history(self.path)
        assert result["step"].shape == (5,)

    def test_time_values(self):
        result = parse_force_history(self.path)
        expected = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
        np.testing.assert_allclose(result["time"], expected, atol=1e-10)

    def test_force_column_count(self):
        result = parse_force_history(self.path)
        # 10 cols total: step, time, col_2 .. col_9
        force_keys = [k for k in result if k.startswith("col_")]
        assert len(force_keys) == 8, (
            f"expected 8 force columns (col_2..col_9), got {force_keys}"
        )


# ---------------------------------------------------------------------------
# parse_spectrum_file — generic 2-col
# ---------------------------------------------------------------------------

class TestParseSpectrumGeneric2Col:
    def setup_method(self):
        lines = []
        for i in range(6):
            r = 0.9 + 0.05 * i
            im = 0.1 * i - 0.25
            lines.append(f"  {r:.7E}  {im:.7E}\n")
        self.path = _write_tmp("".join(lines))
        self.n = 6

    def teardown_method(self):
        os.unlink(self.path)

    def test_shape(self):
        ev = parse_spectrum_file(self.path)
        assert ev.shape == (self.n, 2), f"expected ({self.n}, 2), got {ev.shape}"

    def test_values(self):
        ev = parse_spectrum_file(self.path)
        np.testing.assert_allclose(ev[0, 0], 0.9, atol=1e-7)


# ---------------------------------------------------------------------------
# parse_spectrum_file — Spectre_*.dat (4-col)
# ---------------------------------------------------------------------------

class TestParseSpectre4Col:
    """Spectre_Hp.dat style: real, imag, magnitude, flag."""

    def setup_method(self):
        pairs = [(3.7e5, 9.0e4), (3.1e5, -1.2e5), (2.7e5, 1.4e5)]
        lines = []
        for r, im in pairs:
            mag = math.hypot(r, im)
            lines.append(f"  {r:.7E}  {im:.7E}  {mag:.7E} 0\n")
        self.path = _write_tmp("".join(lines))
        self.pairs = pairs

    def teardown_method(self):
        os.unlink(self.path)

    def test_shape(self):
        ev = parse_spectrum_file(self.path)
        assert ev.shape == (3, 2)

    def test_real_imag_correct(self):
        ev = parse_spectrum_file(self.path)
        for idx, (r, im) in enumerate(self.pairs):
            np.testing.assert_allclose(ev[idx, 0], r, rtol=1e-6,
                err_msg=f"row {idx}: real part mismatch")
            np.testing.assert_allclose(ev[idx, 1], im, rtol=1e-6,
                err_msg=f"row {idx}: imag part mismatch")


# ---------------------------------------------------------------------------
# parse_spectrum_file — DMD with '# DMD' header (8-col)
# ---------------------------------------------------------------------------

class TestParseSpectrumDMDHeader:
    """dmd_spectrum.dat: mu_real in col5, mu_imag in col6."""

    def setup_method(self):
        # Columns: mode |mu| sigma omega St mu_real mu_imag ||Phi||
        mu_reals = [-0.4896, 0.8673, 0.5047]
        mu_imags = [0.8717, 0.4973, 0.8629]
        lines = [
            "# DMD Eigenvalue Spectrum\n",
            "# mode   |mu|   sigma   omega   St   mu_real   mu_imag   ||Phi||\n",
        ]
        for i, (mr, mi) in enumerate(zip(mu_reals, mu_imags), start=1):
            mag = math.hypot(mr, mi)
            lines.append(
                f"  {i:4d}  {mag:.6E}  0.0  0.0  0.0  {mr:.6E}  {mi:.6E}  1.0\n"
            )
        self.path = _write_tmp("".join(lines))
        self.mu_reals = mu_reals
        self.mu_imags = mu_imags

    def teardown_method(self):
        os.unlink(self.path)

    def test_shape(self):
        ev = parse_spectrum_file(self.path)
        assert ev.shape == (3, 2)

    def test_mu_real_used(self):
        ev = parse_spectrum_file(self.path)
        np.testing.assert_allclose(ev[:, 0], self.mu_reals, rtol=1e-5,
            err_msg="DMD: real part should come from col5 (mu_real)")

    def test_mu_imag_used(self):
        ev = parse_spectrum_file(self.path)
        np.testing.assert_allclose(ev[:, 1], self.mu_imags, rtol=1e-5,
            err_msg="DMD: imag part should come from col6 (mu_imag)")


# ---------------------------------------------------------------------------
# detect_unit_circle_modes
# ---------------------------------------------------------------------------

class TestDetectUnitCircleModes:
    def _make_eigvals(self):
        """Mix of on-circle, near-circle, and far-from-circle."""
        data = np.array([
            [1.0, 0.0],        # idx 0 — exactly on circle, magnitude = 1
            [0.0, 1.0],        # idx 1 — exactly on circle
            [0.7071, 0.7071],  # idx 2 — magnitude ≈ 1.0  (on circle)
            [0.5, 0.0],        # idx 3 — magnitude = 0.5, off circle
            [2.0, 0.0],        # idx 4 — magnitude = 2.0, off circle
            [0.9995, 0.0],     # idx 5 — within tol=1e-3 of unit circle
        ])
        return data

    def test_returns_on_circle_indices(self):
        ev = self._make_eigvals()
        indices = detect_unit_circle_modes(ev, tol=1e-3)
        # Indices 0, 1, 2, 5 should be on or near unit circle
        for expected_idx in [0, 1, 5]:
            assert expected_idx in indices, (
                f"index {expected_idx} (|mu|≈1) should be in unit-circle set"
            )

    def test_off_circle_excluded(self):
        ev = self._make_eigvals()
        indices = detect_unit_circle_modes(ev, tol=1e-3)
        for off_idx in [3, 4]:
            assert off_idx not in indices, (
                f"index {off_idx} (|mu| far from 1) should not be in unit-circle set"
            )

    def test_returns_ndarray(self):
        ev = self._make_eigvals()
        indices = detect_unit_circle_modes(ev)
        assert isinstance(indices, np.ndarray)


# ---------------------------------------------------------------------------
# leading_modes
# ---------------------------------------------------------------------------

class TestLeadingModes:
    def _make_eigvals(self):
        """Magnitudes: 0.5, 3.0, 1.0, 2.5, 0.1, 4.0."""
        data = np.array([
            [0.5,  0.0],   # |0| = 0.5
            [3.0,  0.0],   # |1| = 3.0
            [1.0,  0.0],   # |2| = 1.0
            [2.5,  0.0],   # |3| = 2.5
            [0.1,  0.0],   # |4| = 0.1
            [4.0,  0.0],   # |5| = 4.0
        ])
        return data

    def test_top3_correct(self):
        ev = self._make_eigvals()
        idx = leading_modes(ev, k=3)
        # Expected order: 5 (4.0), 1 (3.0), 3 (2.5)
        np.testing.assert_array_equal(
            idx, [5, 1, 3],
            err_msg="top-3 leading modes should be indices 5, 1, 3 (by magnitude)"
        )

    def test_top1(self):
        ev = self._make_eigvals()
        idx = leading_modes(ev, k=1)
        assert idx[0] == 5, "leading mode should be index 5 (magnitude 4.0)"

    def test_k_larger_than_n(self):
        ev = self._make_eigvals()
        idx = leading_modes(ev, k=100)
        assert len(idx) == len(ev), "k > N should be clamped to N"

    def test_descending_order(self):
        ev = self._make_eigvals()
        idx = leading_modes(ev, k=4)
        mags = np.hypot(ev[idx, 0], ev[idx, 1])
        assert np.all(np.diff(mags) <= 0), "leading modes should be in descending magnitude order"

    def test_returns_ndarray(self):
        ev = self._make_eigvals()
        idx = leading_modes(ev, k=2)
        assert isinstance(idx, np.ndarray)
