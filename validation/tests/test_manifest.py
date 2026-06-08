"""Tests the standalone per-case run provenance manifest writer."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "write_manifest.py"


def test_manifest_happy_path(tmp_path: Path) -> None:
    case_dir = tmp_path / "test_case"
    case_dir.mkdir()
    (case_dir / "test_case.par").write_text(
        "[GENERAL]\n"
        "reynoldsNumber = 47.5\n"
        "userParam07 = 200\n"
    )
    (case_dir / "logfile").write_text("total elapsed time in nekton =    123.45\n")
    (case_dir / "residu.dat").write_text("1 9.99e-01\n2 1.23e-08\n")
    (case_dir / "BF_test0.f00001").write_text("")

    out_path = tmp_path / "manifest.json"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(case_dir),
            "--ranks",
            "4",
            "--jobid",
            "42",
            "--out",
            str(out_path),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    manifest = json.loads(out_path.read_text())

    assert manifest["schema_version"] == "0.1.0"
    assert isinstance(manifest["case_id"], str)
    assert manifest["case_id"]
    assert isinstance(manifest["run_date_utc"], str)
    assert manifest["run_date_utc"]
    assert manifest["sbatch_jobid"] == 42
    assert manifest["ranks"] == 4
    assert manifest["final_residual"] == pytest.approx(1.23e-08)
    assert manifest["parameters"]["reynolds"] == pytest.approx(47.5)
    assert str(manifest["parameters"]["userparam07"]).lower() == "200"
    assert any(artifact["role"] == "CHECKPOINT" for artifact in manifest["artifacts"])
    assert manifest["slurm_log"] == "logfile"
    assert manifest["wallclock_seconds"] == pytest.approx(123.45)
