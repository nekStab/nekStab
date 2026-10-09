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


def write_fake_fld(path: Path, nz: int = 1, wdsize: int = 4) -> Path:
    """Write a small field file in the layout of Nek5000: header, test float, element map, data.

    The coordinates and the velocity are stored element by element (all components of element 1,
    then those of element 2), as pressure comes as one block. nz = 1 gives a 2D file, nz = 2 a 3D file.
    """
    nx = 2
    ny = 2
    nelt = 2
    time = 12.5
    istep = 25
    rdcode = "XUP"
    header = (
        f"#std {wdsize} {nx} {ny} {nz} {nelt} {nelt} {time:.6e} "
        f"{istep} 0 1 {rdcode}"
    )
    raw_header = header.encode("ascii")
    assert len(raw_header) <= 132

    ndim = 2 if nz == 1 else 3
    nxyz = nx * ny * nz
    vx_values = [[float(nxyz * e + i + 1) for i in range(nxyz)] for e in range(nelt)]
    vy_values = [-1.0] * nxyz
    vz_values = [0.5] * nxyz
    pressure = [10.0 * (nxyz * e + i + 1) for e in range(nelt) for i in range(nxyz)]

    values: list[float] = []
    for _ in range(nelt):
        values += [0.0] * (ndim * nxyz)
    for e in range(nelt):
        values += vx_values[e] + vy_values + (vz_values if ndim == 3 else [])
    values += pressure
    code = "f" if wdsize == 4 else "d"

    with path.open("wb") as fh:
        fh.write(raw_header.ljust(132, b" "))
        fh.write(struct.pack("<f", 6.54321))
        fh.write(struct.pack(f"<{nelt}i", *range(1, nelt + 1)))
        fh.write(struct.pack(f"<{len(values)}{code}", *values))
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


def test_3d_double_precision(tmp_path: Path) -> None:
    result = run_script("--json", str(write_fake_fld(tmp_path / "fake3d0.f00001", nz=2, wdsize=8)))

    assert result.returncode == 0, result.stderr
    fields = json.loads(result.stdout)["fields"]
    assert fields["vx"]["min"] == pytest.approx(1.0)
    assert fields["vx"]["max"] == pytest.approx(16.0)
    assert fields["vy"]["min"] == pytest.approx(-1.0)
    assert fields["vz"]["max"] == pytest.approx(0.5)
    assert fields["pr"]["max"] == pytest.approx(160.0)


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
