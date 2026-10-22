"""Tests of bin/c0check, the continuity check of Nek5000 field files."""

from __future__ import annotations

import importlib.machinery
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "c0check"
NX = 3


def load_tool() -> ModuleType:
    """Import bin/c0check, which has no .py suffix."""
    loader = importlib.machinery.SourceFileLoader("c0check", str(SCRIPT))
    spec = importlib.util.spec_from_loader("c0check", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["c0check"] = module  # dataclasses look the module up here
    loader.exec_module(module)
    return module


def two_elements() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return the coordinates and a smooth u, v of two 2D elements.

    Element 1 covers 0 <= x <= 1, element 2 covers 1 <= x <= 2, both
    0 <= y <= 1. They share the edge x = 1 (3 points).
    """
    g = np.linspace(0.0, 1.0, NX)
    gx, gy = np.meshgrid(g, g)
    coords = np.stack(
        [np.stack([gx.ravel() + e, gy.ravel()]) for e in range(2)]
    )  # (nel, 2, nxyz)
    x, y = coords[:, 0, :], coords[:, 1, :]
    vel = np.stack([np.sin(x) * np.cos(y), x * y], axis=1)  # (nel, 2, nxyz)
    return coords, vel


def write_field(
    path: Path,
    coords: npt.NDArray[np.float64] | None,
    vel: npt.NDArray[np.float64],
    elmap: list[int],
    wdsize: int = 8,
) -> Path:
    """Write a 2D Nek5000 field file: header, test float, element map, data."""
    nel = vel.shape[0]
    rdcode = ("X" if coords is not None else "") + "UP"
    hdr = f"#std {wdsize} {NX} {NX} 1 {nel} {nel} 0.0 1 0 1 {rdcode}"
    real = np.float64 if wdsize == 8 else np.float32
    with path.open("wb") as fh:
        fh.write(hdr.ljust(132).encode("ascii"))
        fh.write(np.array([6.54321], "<f4").tobytes())
        fh.write(np.array(elmap, "<i4").tobytes())
        if coords is not None:
            fh.write(coords.astype(real).tobytes())
        fh.write(vel.astype(real).tobytes())
        fh.write(np.zeros((nel, NX * NX), real).tobytes())
    return path


def test_continuous_field_passes(tmp_path: Path) -> None:
    tool = load_tool()
    coords, vel = two_elements()
    path = write_field(tmp_path / "ok0.f00001", coords, vel, [2, 1])
    field = tool.read_field(path)
    report = tool.check(field, field.coords, 1e-10)
    assert report.nranks == 2
    assert report.shared_by_nranks == {2: (3, 0)}
    assert report.max_rel_jump == 0.0
    assert tool.main([str(path)]) == 0


def test_jump_on_one_copy_is_found(tmp_path: Path) -> None:
    tool = load_tool()
    coords, vel = two_elements()
    vel[1, 0, 3] += 0.5  # element 2, point (x = 1, y = 0.5): one copy only
    path = write_field(tmp_path / "bad0.f00001", coords, vel, [2, 1])
    field = tool.read_field(path)
    report = tool.check(field, field.coords, 1e-10)
    assert report.shared_by_nranks == {2: (3, 1)}
    x, y, _, jump = report.worst[0]
    assert (x, y) == pytest.approx((1.0, 0.5))
    assert jump == pytest.approx(0.5 / np.abs(vel[:, 0, :]).max())
    assert tool.main([str(path)]) == 1


def test_same_rank_copies_count_as_one_rank(tmp_path: Path) -> None:
    tool = load_tool()
    coords, vel = two_elements()
    path = write_field(tmp_path / "one0.f00001", coords, vel, [1, 2])
    report = tool.check(tool.read_field(path), coords, 1e-10)
    assert report.nranks == 1
    assert report.shared_by_nranks == {1: (3, 0)}


def test_32bit_file_uses_the_single_precision_limit(tmp_path: Path) -> None:
    tool = load_tool()
    coords, vel = two_elements()
    path = write_field(tmp_path / "sp0.f00001", coords, vel, [2, 1], wdsize=4)
    assert tool.read_field(path).wdsize == 4
    assert tool.main([str(path)]) == 0


def test_field_without_mesh_reads_coords_from_another_file(tmp_path: Path) -> None:
    tool = load_tool()
    coords, vel = two_elements()
    mesh = write_field(tmp_path / "bf0.f00001", coords, vel, [2, 1])
    vel[0, 1, 5] -= 1.0  # element 1, point (x = 1, y = 0.5)
    mode = write_field(tmp_path / "dRe0.f00002", None, vel, [2, 1])
    assert tool.main([str(mode)]) == 2
    assert tool.main(["--coords", str(mesh), str(mode)]) == 1
