#!/usr/bin/env python3
"""Cross-method baseflow convergence comparison for cylinder_re100.

Overlays Newton, all 6 SFD variants, and BoostConv on a single semilog plot.
Run from repo root: python3 scripts/plot_cross_method_convergence.py
"""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FAMILY = ROOT / "example" / "cylinder_re100"
OUTDIR = ROOT / "validation" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

METHODS = [
    # label, file, color, linestyle, time-col, residu-col
    ("Newton",                "210_baseflow_newton/fp/residu_newton.dat",            "k",         "-",  4, 6),
    ("SFD Akervik static",    "110_baseflow_sfd/akervik/residu.dat",                 "#d62728",   "-",  0, 1),
    ("SFD Akervik + dyn-tol", "110_baseflow_sfd/akervik_dyn/residu.dat",             "#d62728",   "--", 0, 1),
    ("SFD Akervik + OIFS",    "110_baseflow_sfd/akervik_dyn_oifs/residu.dat",        "#d62728",   ":",  0, 1),
    ("SFD Casacuberta static","110_baseflow_sfd/casacuberta/residu.dat",             "#1f77b4",   "-",  0, 1),
    ("SFD Casacub. + dyn-tol","110_baseflow_sfd/casacuberta_dyn/residu.dat",         "#1f77b4",   "--", 0, 1),
    ("SFD Casacub. + OIFS",   "110_baseflow_sfd/casacuberta_dyn_oifs/residu.dat",    "#1f77b4",   ":",  0, 1),
    ("BoostConv",             "120_baseflow_boostconv/residu.dat",                   "#2ca02c",   "-",  0, 1),
]

fig, ax = plt.subplots(figsize=(8, 5.2))

for label, fname, color, ls, tcol, rcol in METHODS:
    p = FAMILY / fname
    if not p.exists():
        print(f"  skip (missing): {fname}")
        continue
    data = np.loadtxt(p)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    t = data[:, tcol]
    r = np.abs(data[:, rcol])
    r[r < 1e-15] = 1e-15
    marker = "o" if label == "Newton" else None
    ms = 7 if label == "Newton" else 0
    lw = 1.6 if label == "Newton" else 1.0
    ax.plot(t, r, color=color, linestyle=ls, marker=marker, markersize=ms,
            linewidth=lw, label=label)
    print(f"  loaded: {fname}  ({len(t)} pts, final residu={r[-1]:.2e})")

ax.set_yscale("log")
ax.set_xlabel("simulated time")
ax.set_ylabel("baseflow residual")
ax.set_title("Cylinder Re=100 -- baseflow convergence across methods (Story 2)")
ax.axhline(1e-9, color="gray", linestyle=":", linewidth=0.5, label="tol=1e-9")
ax.set_ylim(1e-11, 10)
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="upper right", fontsize=8, ncol=2)

fig.tight_layout()
out = OUTDIR / "cylinder_re100_baseflow_convergence.png"
fig.savefig(out, dpi=150)
print(f"saved: {out}")

# Companion table: final-state agreement (vy norm)
print()
print("Final-state vy-norm agreement vs Newton:")
import subprocess, re
def inspect(p):
    o = subprocess.run(["python3", "scripts/inspect_nek_field.py", str(p)],
                       capture_output=True, text=True, cwd=ROOT).stdout
    m = re.search(r"vy\s+min=\S+\s+max=\s*\S+\s+norm=\s*([\d.eE+-]+)", o)
    return float(m.group(1)) if m else None
ref = inspect(FAMILY / "210_baseflow_newton/fp/BF_1cyl0.f00001")
print(f"  Newton (reference): vy_norm = {ref:.6f}")
for label, fname, *_ in METHODS[1:]:
    case_dir = FAMILY / Path(fname).parent
    bf = case_dir / "BF_1cyl0.f00001"
    if bf.exists():
        v = inspect(bf)
        diff = abs(v - ref) / ref * 100 if ref else float("nan")
        print(f"  {label:30s}  vy_norm = {v:.6f}  (Delta = {diff:.3f}%)")
