"""Negative controls for scripts/check_against_ref.py.

A reference check is worth something only if it fails when the run is wrong. Each test
builds a small stage folder (a reference, a column file, a log, a DONE file) and checks
that the verdict follows: pass for a good run; fail for a shifted value, a missing output,
a bad exit code, an ERROR line in the log, a NaN in the log, and an empty list of quantities;
the advice that names a non-finite norm in quotes is not a failure.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "check_against_ref.py"


def make_stage(
    tmp_path: Path,
    value: float = 0.125,
    log: str = "total elapsed time : 1.0E+00 sec\n",
    done: str = "0",
    quantities: dict | None = None,
) -> Path:
    """A stage folder with residu.dat (last column-2 value = `value`) and a reference of 0.125."""
    stage = tmp_path / "stage"
    (stage / "ref").mkdir(parents=True)
    (stage / "residu.dat").write_text(f"1 0.5\n2 {value}\n")
    (stage / "logfile").write_text(log)
    if done is not None:
        (stage / "DONE").write_text(done + "\n")
    if quantities is None:
        quantities = {
            "last_residual": {
                "source": "residu.dat:col2:last",
                "value": 0.125,
                "tol_rel": 1e-3,
            }
        }
    (stage / "ref" / "reference.json").write_text(
        json.dumps({"case": "test", "quantities": quantities})
    )
    return stage


def run_check(stage: Path) -> tuple[int, str]:
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(stage)], capture_output=True, text=True
    )
    return proc.returncode, proc.stdout + proc.stderr


def test_good_run_passes(tmp_path: Path) -> None:
    code, out = run_check(make_stage(tmp_path))
    assert code == 0, out
    assert "RESULT: PASS" in out


def test_shifted_value_fails(tmp_path: Path) -> None:
    code, out = run_check(make_stage(tmp_path, value=0.13))
    assert code == 1, out
    assert "FAIL" in out


def test_missing_output_fails(tmp_path: Path) -> None:
    stage = make_stage(tmp_path)
    (stage / "residu.dat").unlink()
    code, out = run_check(stage)
    assert code == 1, out
    assert "ERROR" in out


def test_bad_exit_code_fails(tmp_path: Path) -> None:
    code, out = run_check(make_stage(tmp_path, done="3"))
    assert code == 1, out
    assert "DONE holds exit code 3" in out


def test_error_line_in_log_fails(tmp_path: Path) -> None:
    log = "ERROR: ABORT: Could not open Spectre_NSd_conv.dat\n"
    code, out = run_check(make_stage(tmp_path, log=log))
    assert code == 1, out
    assert "run gate" in out


def test_nan_in_log_fails(tmp_path: Path) -> None:
    code, out = run_check(make_stage(tmp_path, log="residual = NaN\n"))
    assert code == 1, out
    assert "run gate" in out


def test_non_finite_norm_message_fails(tmp_path: Path) -> None:
    log = "ERROR [k_normalize]: non-finite norm; refusing to scale Krylov vector.\n"
    code, out = run_check(make_stage(tmp_path, log=log))
    assert code == 1, out
    assert "run gate" in out


def test_advice_about_non_finite_norm_passes(tmp_path: Path) -> None:
    log = (
        " nekStab: NOTE: [VELOCITY] residualProj = yes in a Newton or\n"
        '   Krylov mode. If the run stops with "non-finite norm", set\n'
        "   residualProj = no in [VELOCITY]; naca0012 needed it at 24+ ranks.\n"
    )
    code, out = run_check(make_stage(tmp_path, log=log))
    assert code == 0, out


def test_unset_time_step_line_passes(tmp_path: Path) -> None:
    # Nek5000 prints DTCFL and DTFS unset when the time step is fixed; a build can show NaN there.
    log = "     DT/DTCFL/DTFS/DTINIT   0.250E-01         NaN         NaN   0.000E+00\n"
    code, out = run_check(make_stage(tmp_path / "clean", log=log))
    assert code == 0, out
    code, out = run_check(make_stage(tmp_path / "nan", log=log + "residual = NaN\n"))
    assert code == 1, out


def test_empty_quantity_list_fails(tmp_path: Path) -> None:
    code, out = run_check(make_stage(tmp_path, quantities={}))
    assert code == 1, out
    assert "checks nothing" in out
