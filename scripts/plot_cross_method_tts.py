#!/usr/bin/env python3
"""Cross-method baseflow TTS (time-to-solution) comparison for cylinder_re100.

Plots residual trajectories, with vertical annotations at t@1e-6 (the unified
convergence threshold) showing which method reaches the BF first.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FAMILY = ROOT / "example" / "cylinder_re100"
OUTDIR = ROOT / "validation" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

CASES = [
    # label, file, color, ls, tcol, rcol
    ("Newton",                  "210_baseflow_newton/fp/residu_newton.dat",         "k",       "-",  4, 6),
    ("SFD Akervik static",      "110_baseflow_sfd/akervik/residu.dat",              "#d62728", "-",  0, 1),
    ("SFD Akervik + dyn",       "110_baseflow_sfd/akervik_dyn/residu.dat",          "#d62728", "--", 0, 1),
    ("SFD Akervik + OIFS",      "110_baseflow_sfd/akervik_dyn_oifs/residu.dat",     "#d62728", ":",  0, 1),
    ("SFD Casacuberta static",  "110_baseflow_sfd/casacuberta/residu.dat",          "#1f77b4", "-",  0, 1),
    ("SFD Casacub. + dyn",      "110_baseflow_sfd/casacuberta_dyn/residu.dat",      "#1f77b4", "--", 0, 1),
    ("SFD Casacub. + OIFS",     "110_baseflow_sfd/casacuberta_dyn_oifs/residu.dat", "#1f77b4", ":",  0, 1),
    ("BoostConv",               "120_baseflow_boostconv/residu.dat",                "#2ca02c", "-",  0, 1),
]

THRESH = 1e-6  # unified convergence threshold

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [3, 1.4]})

tts_data = []

for label, fname, color, ls, tcol, rcol in CASES:
    p = FAMILY / fname
    if not p.exists():
        continue
    data = np.loadtxt(p)
    if data.ndim == 1: data = data.reshape(1, -1)
    t = data[:, tcol]
    r = np.abs(data[:, rcol])
    r[r < 1e-15] = 1e-15

    # TTS = first time residu < THRESH
    below = np.where(r < THRESH)[0]
    tts = t[below[0]] if len(below) else None
    tts_data.append((label, tts, r[-1], color))

    marker = "o" if label == "Newton" else None
    ms = 7 if label == "Newton" else 0
    lw = 1.6 if label == "Newton" else 1.0
    ax1.plot(t, r, color=color, linestyle=ls, marker=marker, markersize=ms,
             linewidth=lw, label=label, alpha=0.85)

ax1.set_yscale("log")
ax1.set_xlabel("simulated time")
ax1.set_ylabel("residual")
ax1.set_title("Baseflow convergence trajectories (residualTol unified to 1e-6)")
ax1.axhline(THRESH, color="gray", linestyle="-.", linewidth=0.8, label=f"tol = {THRESH:.0e}")
ax1.set_ylim(1e-11, 10)
ax1.grid(True, which="both", alpha=0.3)
ax1.legend(loc="upper right", fontsize=8, ncol=2)

# Right panel: TTS bar chart
labels = [t[0] for t in tts_data]
tts_vals = [(t[1] if t[1] is not None else float("nan")) for t in tts_data]
colors = [t[3] for t in tts_data]
y_pos = np.arange(len(labels))[::-1]
bars = ax2.barh(y_pos, tts_vals, color=colors, alpha=0.8, edgecolor="black", linewidth=0.4)
for i, (label, tts, fin, color) in enumerate(tts_data):
    yp = y_pos[i]
    if tts is None:
        ax2.text(2, yp, f"  (no cross, plateau {fin:.1e})", va="center", fontsize=7, style="italic", color="gray")
    else:
        ax2.text(tts + 15, yp, f"  t = {tts:.0f}", va="center", fontsize=8)
ax2.set_yticks(y_pos)
ax2.set_yticklabels(labels, fontsize=8)
ax2.set_xlabel("simulated time to reach 1e-6")
ax2.set_title("TTS @ 1e-6")
ax2.set_xlim(0, 1100)
ax2.grid(True, axis="x", alpha=0.3)
ax2.invert_yaxis()

fig.suptitle("Cylinder Re=100 -- Story 2: which baseflow method gets there first?", fontsize=11, y=1.02)
fig.tight_layout()
out = OUTDIR / "cylinder_re100_baseflow_tts.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved: {out}")

# Print summary table
print()
print(f"TTS summary (unified threshold = {THRESH:.0e}):")
print(f"{'Method':<27} {'t@1e-6':>9}  {'final_res':>10}")
print("-" * 50)
for label, tts, fin, _ in tts_data:
    tts_s = f"{tts:.0f}" if tts is not None else "did not cross"
    print(f"{label:<27} {tts_s:>9}  {fin:>10.2e}")
