#!/usr/bin/env python3
"""Coupled operator: the near-wall k-tau continuum buries the wake mode.

Sorted growth-rate spectra at k_dim=100 and 200 both plateau in sigma~2.3-3.9 and
never approach the physical wake mode (the quasilaminar leading sigma=0.40).
Doubling the subspace barely lowers the floor -> the wake mode is unreachable
by direct Arnoldi on the coupled operator; the quasilaminar operator
(prescribed base mu_t) is required, not merely preferred.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A100 = np.loadtxt("coupled/Spectre_NSd_kdim100.dat")[:, 0]
A200 = np.loadtxt("coupled/Spectre_NSd.dat")[:, 0]
B = np.loadtxt("quasilaminar/Spectre_NSd.dat")[:, 0]
sB = B.max()

fig, ax = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)
ax.plot(range(1, len(A100)+1), sorted(A100, reverse=True), "o-", ms=3, lw=0.8,
        color="salmon", label=f"coupled, k_dim=100  (floor σ={A100.min():.2f})")
ax.plot(range(1, len(A200)+1), sorted(A200, reverse=True), "s-", ms=3, lw=0.8,
        color="C3", label=f"coupled, k_dim=200  (floor σ={A200.min():.2f})")
ax.axhline(sB, color="C0", ls="--", lw=2,
           label=f"quasilaminar wake mode  σ={sB:.2f} (St=0.206)")
ax.axhline(0, color="k", lw=0.8)
ax.axhspan(sB, A200.min(), alpha=0.12, color="grey")
ax.text(100, 0.5*(sB + A200.min()), "GAP — coupled has NO modes here\n(wake mode unreachable)",
        ha="center", va="center", fontsize=9, color="dimgray")
ax.set_xlabel("mode index (sorted by growth rate, descending)")
ax.set_ylabel("σ  (growth rate)")
ax.set_title("Coupled k-τ operator: near-wall continuum buries the wake mode\n"
             "doubling k_dim only lowers the floor 2.50→2.27 — the physical mode never surfaces")
ax.legend(loc="upper right"); ax.grid(alpha=0.3); ax.set_ylim(-1.6, 4.1)
fig.savefig("coupled_continuum.png", dpi=200)
print(f"A100 floor={A100.min():.3f}  A200 floor={A200.min():.3f}  quasilaminar wake σ={sB:.3f}")
print("saved coupled_continuum.png")
