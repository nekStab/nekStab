#!/usr/bin/env python3
"""Inspect Nek5000 field-file headers to verify IC provenance before runs.

Reads the 132-byte ASCII header of one or more `.f0000*` files and prints
nelt, time, polynomial order, and the scalar-field code. Use this before
launching a Newton/SFD/BoostConv case to confirm:

- the IC mesh matches the case mesh (nelt agrees with `SIZE`/`.re2`);
- the IC time is what the algorithm expects (0 for fp Newton; on-orbit
  for UPO Newton; matches `startTime` for restarts);
- thermal cases load XUPT (not XUP), which proves the scalar field is
  populated and coupled.

Header format (Nek5000 fld v1, ASCII):

    #std <wdsize> <lx1> <ly1> <lz1> <nelt_per_pid> <nelt_total> <time>
         <istep> <fid0> <nfileoo> <rdcode> <ranks> <fp_test_value>

Usage:

    python3 scripts/inspect_nek_ic.py PATH [PATH ...]

Exits 0 on success even if some files are missing (logged inline).
"""
from __future__ import annotations

import sys
from pathlib import Path

HDR_BYTES = 132


def inspect(path: Path) -> int:
    if not path.exists():
        print(f"  MISSING  {path}")
        return 1
    with path.open("rb") as fh:
        raw = fh.read(HDR_BYTES)
    text = raw.decode("ascii", errors="replace").rstrip("\x00 ")
    parts = text.split()
    if len(parts) < 12:
        print(f"  UNREADABLE  {path}  ({text!r})")
        return 1
    _, wdsize, lx, ly, lz, _, nelt, time_s, istep = parts[:9]
    rdcode = parts[11]
    size = path.stat().st_size
    print(f"  {path}")
    print(
        f"    nelt={nelt}  lx={lx}  ly={ly}  lz={lz}  "
        f"wdsize={wdsize}  time={time_s}  istep={istep}  "
        f"rdcode={rdcode}  size_bytes={size}"
    )
    return 0


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(f"usage: {argv[0]} PATH [PATH ...]")
        return 2
    rc = 0
    for arg in argv[1:]:
        rc |= inspect(Path(arg))
    return 0 if rc == 0 else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
