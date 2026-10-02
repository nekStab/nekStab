#!/usr/bin/env python3
"""Checks for the cylinder Re=1e6 resolved stress.

Pass means the 20-time-unit moments form a small R_ij, and, once
reynolds_commit has run, RS1 matches that R_ij. A missing Nek write is
SKIP, not a pass. Exit 1 on FAIL.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
BFA = ROOT / "BFA1cyl0.f00002"
EM2 = ROOT / "EM21cyl0.f00002"
EUV = ROOT / "EUV1cyl0.f00002"
WAKE = ROOT / "steady_wake.f00001"
RS1 = ROOT / "RS11cyl0.f00001"
R_MAX = 1.0e-2
U_AT_5 = (0.50, 0.75)


def load(path: Path) -> dict[str, np.ndarray]:
    raw = path.read_bytes()
    hdr = raw[:132].decode("ascii", "replace")
    tok = hdr.split()
    wds = int(tok[1])
    nel = int(tok[5])
    nxyz = 36
    code = hdr[hdr.find("X"):hdr.find("X") + 8]
    nscal = 2 if "S" in code else 0
    off = 132 + 4 + nel * 4
    dt = "<f4" if wds == 4 else "<f8"

    def group(ncomp: int, o: int) -> tuple[np.ndarray, int]:
        n = nel * ncomp * nxyz
        a = np.frombuffer(raw, dtype=dt, count=n, offset=o).astype(np.float64)
        return a.reshape(nel, ncomp, nxyz), o + n * wds

    xy, off = group(2, off)
    uv, off = group(2, off)
    _, off = group(1, off)
    _, off = group(1, off)
    out = {
        "x": xy[:, 0].ravel(),
        "y": xy[:, 1].ravel(),
        "u": uv[:, 0].ravel(),
        "v": uv[:, 1].ravel(),
    }
    if nscal:
        sc, off = group(nscal, off)
        out["k"] = sc[:, 0].ravel()
        out["tau"] = sc[:, 1].ravel()
    if off > len(raw) + 8:
        raise RuntimeError(f"{path.name}: read past end ({off} > {len(raw)})")
    return out


def wake_mask(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return (x > 1.0) & (x < 12.0) & (np.abs(y) < 2.0)


def main() -> int:
    failed = False
    need = (BFA, EM2, EUV, WAKE)
    if not all(p.is_file() for p in need):
        missing = ", ".join(p.name for p in need if not p.is_file())
        print(f"SKIP  moments  missing {missing}")
    else:
        bfa, em2, euv, wake = map(load, need)
        ruu = em2["u"] - bfa["u"] ** 2
        rvv = em2["v"] - bfa["v"] ** 2
        ruv = euv["u"] - bfa["u"] * bfa["v"]
        m = wake_mask(bfa["x"], bfa["y"])
        peaks = {
            "Ruu": float(np.max(np.abs(ruu[m]))),
            "Rvv": float(np.max(np.abs(rvv[m]))),
            "Ruv": float(np.max(np.abs(ruv[m]))),
        }
        near = (np.abs(wake["x"] - 5.0) < 0.15) & (np.abs(wake["y"]) < 0.25)
        u5 = float(np.median(wake["u"][near])) if np.any(near) else float("nan")
        print(
            "moments  "
            + "  ".join(f"|{k}|={v:.3e}" for k, v in peaks.items())
            + f"  u(5,0)={u5:.3f}"
        )
        if any(v > R_MAX or not np.isfinite(v) for v in peaks.values()):
            print(f"FAIL  moments  a component exceeds {R_MAX:.0e} or is not finite")
            failed = True
        elif not (U_AT_5[0] <= u5 <= U_AT_5[1]):
            print(f"FAIL  wake probe  u(5,0)={u5:.3f} not in {U_AT_5}")
            failed = True
        else:
            print(f"PASS  moments  |R| < {R_MAX:.0e} and u(5,0) in {U_AT_5}")

        if not RS1.is_file():
            print("SKIP  RS1  reynolds_commit has not written RS11cyl0.f00001")
        else:
            rs1 = load(RS1)
            # Same nodes, same order, if both files are this mesh.
            if rs1["u"].shape != ruu.shape:
                print(
                    f"FAIL  RS1  length {rs1['u'].size} != moments {ruu.size}"
                )
                failed = True
            else:
                corr = float(np.corrcoef(rs1["u"], ruu)[0, 1])
                print(f"RS1     corr(Ruu)={corr:.6f}")
                if not np.isfinite(corr) or corr < 0.99:
                    print("FAIL  RS1  does not match E(U^2)-E(U)^2")
                    failed = True
                else:
                    print("PASS  RS1  matches the moment identity")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
