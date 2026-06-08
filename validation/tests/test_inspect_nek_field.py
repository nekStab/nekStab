from __future__ import annotations

import json
import math
import os
import re
import struct
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "inspect_nek_field.py"
REAL_FIELD = ROOT / "example" / "cylinder" / "dns" / "1cyl0.f00001"


def write_fake_fld(path: Path) -> Path:
    nx = 2
    ny = 2
    nz = 1
    nelt = 2
    time = 12.5
    istep = 25
    rdcode = "XUP"
    header = (
        f"#std 4 {nx} {ny} {nz} {nelt} {nelt} {time:.6e} "
        f"{istep} 0 1 {rdcode} 1 6.543210e-01"
    )
    raw_header = header.encode("ascii")
    assert len(raw_header) <= 132

    coords = [0.0] * (nelt * 3 * nx * ny * nz)
    vx_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    vy_values = [-1.0] * 8
    vz_values = [0.5] * 8
    pressure = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    values = coords + vx_values + vy_values + vz_values + pressure

    with path.open("wb") as fh:
        fh.write(raw_header.ljust(132, b" "))
        fh.write(struct.pack(f"{len(values)}f", *values))
    return path


def run_script(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["INSPECT_NEK_FIELD_NO_PYMECH"] = "1"
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=str(ROOT),
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def fake_file(tmp_path: Path) -> Path:
    return write_fake_fld(tmp_path / "fake0.f00001")


def parse_field_line(output: str, field: str) -> tuple[float, float, float]:
    match = re.search(
        rf"^\s*{field}\s+min=\s*([-+0-9.eE]+)\s+max=\s*([-+0-9.eE]+)\s+norm=\s*([-+0-9.eE]+)",
        output,
        re.MULTILINE,
    )
    assert match, output
    return tuple(float(group) for group in match.groups())


def test_header_parsed(tmp_path: Path) -> None:
    result = run_script(str(fake_file(tmp_path)))

    assert result.returncode == 0, result.stderr
    assert "nelt=2" in result.stdout
    assert "time=1.250000e+01" in result.stdout


def test_min_max_vx(tmp_path: Path) -> None:
    result = run_script(str(fake_file(tmp_path)))

    assert result.returncode == 0, result.stderr
    vmin, vmax, norm = parse_field_line(result.stdout, "vx")
    assert vmin == pytest.approx(1.0)
    assert vmax == pytest.approx(8.0)
    assert norm == pytest.approx(math.sqrt(sum(value * value for value in range(1, 9))))


def test_json_output(tmp_path: Path) -> None:
    result = run_script("--json", str(fake_file(tmp_path)))

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["nelt"] == 2
    assert set(payload["fields"]["vx"]) == {"min", "max", "norm"}
    assert payload["fields"]["vx"]["min"] == pytest.approx(1.0)
    assert payload["fields"]["vx"]["max"] == pytest.approx(8.0)


def test_missing_file(tmp_path: Path) -> None:
    result = run_script(str(tmp_path / "missing.f00001"))

    assert result.returncode == 1
    assert "UNREADABLE" in result.stdout


@pytest.mark.skipif(not REAL_FIELD.exists(), reason="real cylinder DNS field is unavailable")
def test_real_file_smoke() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(REAL_FIELD)],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "vx" in result.stdout
