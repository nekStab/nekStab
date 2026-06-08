#!/usr/bin/env python3
"""Inspect Nek5000 field files and report basic per-field statistics.

Usage:

    python3 scripts/inspect_nek_field.py [--json] FIELD [FIELD ...]

The reader tries pymech first. If pymech is unavailable or cannot read the
file, it falls back to a small raw reader for standard packed Nek5000 field
files with a 132-byte ASCII header.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import struct
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

HDR_BYTES = 132


Stats = Dict[str, float]


class InspectError(RuntimeError):
    """Raised when a field file cannot be inspected."""


def _empty_stats() -> Stats:
    return {"min": math.inf, "max": -math.inf, "norm": 0.0}


def _update_stats(stats: Stats, values: Iterable[float]) -> None:
    total = stats["norm"] * stats["norm"]
    seen = False
    for value in values:
        value = float(value)
        seen = True
        if value < stats["min"]:
            stats["min"] = value
        if value > stats["max"]:
            stats["max"] = value
        total += value * value
    if seen:
        stats["norm"] = math.sqrt(total)


def _finalize_stats(stats: Stats) -> Stats:
    if math.isinf(stats["min"]):
        return {"min": 0.0, "max": 0.0, "norm": 0.0}
    return stats


def _array_values(array: object) -> Iterable[float]:
    try:
        return array.ravel()  # type: ignore[attr-defined]
    except AttributeError:
        return _flatten(array)


def _flatten(value: object) -> Iterable[float]:
    if isinstance(value, (list, tuple)):
        for item in value:
            yield from _flatten(item)
    else:
        yield float(value)  # type: ignore[arg-type]


def _stats_from_array(array: object) -> Stats:
    stats = _empty_stats()
    _update_stats(stats, _array_values(array))
    return _finalize_stats(stats)


def _inspect_with_pymech(path: Path) -> Dict[str, object]:
    if os.environ.get("INSPECT_NEK_FIELD_NO_PYMECH"):
        raise InspectError("pymech disabled by INSPECT_NEK_FIELD_NO_PYMECH")

    try:
        import pymech  # type: ignore
    except Exception as exc:
        raise InspectError(f"pymech unavailable: {exc}") from exc

    try:
        ds = pymech.readnek(str(path))
    except Exception as exc:
        raise InspectError(f"pymech failed: {exc}") from exc

    try:
        lx1, ly1, lz1 = int(ds.nx1), int(ds.ny1), int(ds.nz1)
    except AttributeError:
        try:
            lx1, ly1, lz1 = (int(value) for value in ds.lr1)
        except Exception as exc:
            raise InspectError(f"pymech metadata missing polynomial order: {exc}") from exc

    fields: Dict[str, Stats] = {}

    for name, component in (("vx", 0), ("vy", 1), ("vz", 2)):
        stats = _empty_stats()
        seen = False
        for elem in ds.elem:
            if getattr(elem, "vel", None) is None or elem.vel.shape[0] <= component:
                continue
            _update_stats(stats, _array_values(elem.vel[component]))
            seen = True
        if seen:
            fields[name] = _finalize_stats(stats)

    stats = _empty_stats()
    seen = False
    for elem in ds.elem:
        if getattr(elem, "pres", None) is None or elem.pres.shape[0] < 1:
            continue
        _update_stats(stats, _array_values(elem.pres[0]))
        seen = True
    if seen:
        fields["pr"] = _finalize_stats(stats)

    nscalars = 0
    for elem in ds.elem:
        temp = getattr(elem, "temp", None)
        if temp is not None:
            nscalars = max(nscalars, int(temp.shape[0]))
    for scalar in range(nscalars):
        stats = _empty_stats()
        seen = False
        for elem in ds.elem:
            temp = getattr(elem, "temp", None)
            if temp is None or temp.shape[0] <= scalar:
                continue
            _update_stats(stats, _array_values(temp[scalar]))
            seen = True
        if seen:
            fields["t" if scalar == 0 else f"s{scalar + 1}"] = _finalize_stats(stats)

    return {
        "path": str(path),
        "nelt": int(ds.nel),
        "lx1": lx1,
        "ly1": ly1,
        "lz1": lz1,
        "time": float(ds.time),
        "istep": int(ds.istep),
        "fields": fields,
    }


def _parse_header(path: Path) -> Tuple[Dict[str, object], str]:
    with path.open("rb") as fh:
        raw = fh.read(HDR_BYTES)
    if len(raw) != HDR_BYTES:
        raise InspectError("file is shorter than Nek5000 header")
    text = raw.decode("ascii", errors="replace").rstrip("\x00 ")
    parts = text.split()
    if len(parts) < 12:
        raise InspectError(f"unreadable header: {text!r}")
    try:
        meta = {
            "path": str(path),
            "wdsize": int(parts[1]),
            "lx1": int(parts[2]),
            "ly1": int(parts[3]),
            "lz1": int(parts[4]),
            "nelt": int(parts[6]),
            "time": float(parts[7]),
            "istep": int(parts[8]),
            "rdcode": parts[11],
        }
    except (ValueError, IndexError) as exc:
        raise InspectError(f"invalid header: {text!r}") from exc
    return meta, text


def _read_floats(fh: object, count: int, fmt: str) -> List[float]:
    if count <= 0:
        return []
    size = struct.calcsize(fmt)
    data = fh.read(count * size)  # type: ignore[attr-defined]
    if len(data) != count * size:
        raise InspectError("field data ended unexpectedly")
    return list(struct.unpack(f"{count}{fmt}", data))


def _inspect_raw(path: Path) -> Dict[str, object]:
    meta, _ = _parse_header(path)
    wdsize = int(meta["wdsize"])
    if wdsize == 4:
        fmt = "f"
    elif wdsize == 8:
        fmt = "d"
    else:
        raise InspectError(f"unsupported wdsize={wdsize}")

    nelt = int(meta["nelt"])
    nxyz = int(meta["lx1"]) * int(meta["ly1"]) * int(meta["lz1"])
    total_points = nelt * nxyz
    rdcode = str(meta["rdcode"]).upper()
    fields: Dict[str, Stats] = {}

    with path.open("rb") as fh:
        fh.seek(HDR_BYTES)
        if "X" in rdcode:
            _read_floats(fh, total_points * 3, fmt)
        if "U" in rdcode:
            for name in ("vx", "vy", "vz"):
                fields[name] = _stats_from_array(_read_floats(fh, total_points, fmt))
        if "P" in rdcode:
            fields["pr"] = _stats_from_array(_read_floats(fh, total_points, fmt))
        if "T" in rdcode:
            scalar_index = 0
            item_size = struct.calcsize(fmt)
            while True:
                chunk = fh.read(total_points * item_size)
                if not chunk:
                    break
                if len(chunk) != total_points * item_size:
                    raise InspectError("scalar data ended unexpectedly")
                values = struct.unpack(f"{total_points}{fmt}", chunk)
                fields["t" if scalar_index == 0 else f"s{scalar_index + 1}"] = _stats_from_array(values)
                scalar_index += 1

    return {
        "path": str(path),
        "nelt": nelt,
        "lx1": int(meta["lx1"]),
        "ly1": int(meta["ly1"]),
        "lz1": int(meta["lz1"]),
        "time": float(meta["time"]),
        "istep": int(meta["istep"]),
        "fields": fields,
    }


def inspect(path: Path) -> Dict[str, object]:
    if not path.exists():
        raise InspectError("missing file")
    if not path.is_file():
        raise InspectError("not a regular file")

    try:
        return _inspect_with_pymech(path)
    except InspectError as pymech_error:
        try:
            return _inspect_raw(path)
        except InspectError as raw_error:
            raise InspectError(f"{pymech_error}; raw reader failed: {raw_error}") from raw_error


def _format_stats(name: str, stats: Stats) -> str:
    return (
        f"  {name:<3} min={stats['min']: .6e}  "
        f"max={stats['max']: .6e}  norm={stats['norm']: .6e}"
    )


def print_text(results: Sequence[Dict[str, object]]) -> None:
    for result in results:
        print(f"FILE: {result['path']}")
        print(
            f"  nelt={result['nelt']}  lx1={result['lx1']}  ly1={result['ly1']}  "
            f"lz1={result['lz1']}  time={float(result['time']):.6e}  "
            f"istep={result['istep']}"
        )
        fields = result["fields"]
        assert isinstance(fields, dict)
        for name in ("vx", "vy", "vz", "pr", "t"):
            if name in fields:
                print(_format_stats(name, fields[name]))  # type: ignore[arg-type]
        for name in sorted(k for k in fields if k.startswith("s")):
            print(_format_stats(name, fields[name]))  # type: ignore[arg-type]


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="write machine-readable JSON")
    parser.add_argument("fields", nargs="+", metavar="FIELD")
    args = parser.parse_args(argv)

    results: List[Dict[str, object]] = []
    for field in args.fields:
        path = Path(field)
        try:
            results.append(inspect(path))
        except InspectError as exc:
            print(f"UNREADABLE {path}: {exc}", file=sys.stderr if args.json else sys.stdout)

    if args.json:
        payload: object = results[0] if len(results) == 1 else results
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_text(results)

    return 0 if results else 1


if __name__ == "__main__":
    sys.exit(main())
