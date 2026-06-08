#!/usr/bin/env python3
"""check_against_ref.py — verify a freshly run case against its ref/reference.json.

Usage:
    scripts/check_against_ref.py example/<case>/<stage>

Reads <stage>/ref/reference.json, re-extracts each listed quantity from the
freshly produced run output, and reports PASS/FAIL per quantity. Exits non-zero
if any quantity is outside tolerance. See example/REFERENCE_RESULTS.md.

A quantity's "source" tells the extractor where to read the value from a fresh
run, in the form "<file>:col<N>:<reducer>" with reducer in {last,first,min,max},
e.g. "residu.dat:col2:last". Tolerance is "tol_rel" (relative) or "tol_abs".
"""
from __future__ import annotations
import json
import sys
from pathlib import Path


def _column(path: Path, col: int) -> list[float]:
    """All numeric values in a 0-indexed column, skipping header/non-numeric rows."""
    values = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) <= col:
            continue
        try:
            values.append(float(fields[col].replace("D", "E")))
        except ValueError:
            continue
    return values


def _strouhal(path: Path, signal_col: int, time_col: int = 0,
              band=(0.05, 0.5), settle=0.3) -> float:
    """Dominant shedding frequency from FFT of a probe-history signal column."""
    import numpy as np
    t = np.array(_column(path, time_col))
    y = np.array(_column(path, signal_col))
    n = min(len(t), len(y))
    t, y = t[:n], y[:n]
    m = int(n * settle)  # drop initial transient
    t, y = t[m:], y[m:]
    dt = float(np.median(np.diff(t)))
    y = y - y.mean()
    freqs = np.fft.rfftfreq(len(y), d=dt)
    psd = np.abs(np.fft.rfft(y)) ** 2
    sel = (freqs > band[0]) & (freqs < band[1])
    return float(freqs[sel][np.argmax(psd[sel])])


def extract(stage: Path, source: str) -> float:
    """Read a scalar from a fresh run output.

    Spec forms:
      "<file>:col<N>:<reducer>"   reducer in {last,first,min,max}
      "<file>:strouhal:col<N>"    FFT dominant frequency of column N (probe vy)
    """
    parts = source.split(":")
    if len(parts) != 3:
        raise ValueError(f"unsupported source spec: {source!r}")
    fname, op, arg = parts
    path = stage / fname
    if not path.exists():
        raise FileNotFoundError(f"run output not found: {path} (run the case first?)")

    if op == "strouhal":
        return _strouhal(path, int(arg[3:]) - 1)

    if not op.startswith("col"):
        raise ValueError(f"unsupported source spec: {source!r}")
    col = int(op[3:]) - 1
    values = _column(path, col)
    if not values:
        raise ValueError(f"no numeric data in column {col + 1} of {path}")
    reducers = {"last": values[-1], "first": values[0],
                "min": min(values), "max": max(values)}
    if arg not in reducers:
        raise ValueError(f"unknown reducer {arg!r} in {source!r}")
    return reducers[arg]


def check_quantity(name: str, spec: dict, actual: float) -> tuple[bool, str]:
    expected = spec["value"]
    if "tol_rel" in spec:
        denom = abs(expected) if expected != 0 else 1.0
        err = abs(actual - expected) / denom
        tol = spec["tol_rel"]
        kind = "rel"
    elif "tol_abs" in spec:
        err = abs(actual - expected)
        tol = spec["tol_abs"]
        kind = "abs"
    else:
        raise ValueError(f"{name}: need tol_rel or tol_abs")
    ok = err <= tol
    return ok, (f"{name:24s} {'PASS' if ok else 'FAIL'}  "
                f"expected={expected:.6g} actual={actual:.6g} "
                f"{kind}err={err:.3g} tol={tol:.3g}")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    stage = Path(argv[1])
    ref = stage / "ref" / "reference.json"
    if not ref.exists():
        print(f"no reference file: {ref}")
        return 2

    spec = json.loads(ref.read_text())
    print(f"case: {spec.get('case', stage)}")
    prod = spec.get("produced", {})
    if prod:
        print("  reference produced with: "
              + ", ".join(f"{k}={v}" for k, v in prod.items()))

    all_ok = True
    for name, q in spec.get("quantities", {}).items():
        try:
            actual = extract(stage, q["source"])
        except (FileNotFoundError, ValueError) as exc:
            print(f"{name:24s} ERROR  {exc}")
            all_ok = False
            continue
        ok, line = check_quantity(name, q, actual)
        print("  " + line)
        all_ok = all_ok and ok

    print("RESULT:", "PASS — matches reference" if all_ok
          else "FAIL — out of tolerance")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
