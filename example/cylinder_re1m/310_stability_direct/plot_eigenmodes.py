#!/usr/bin/env python3
"""Spatial structure of the leading direct-stability eigenmodes, coupled vs quasilaminar.

Reads the nekStab eigenvector fields dRe<session>0.f0000N (real part of mode N,
in spectrum order) and plots the transverse velocity perturbation v' in the
wake. A physical wake mode = coherent alternating rollers shed downstream; a
spurious/contaminated mode = small-scale or near-body localized structure.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymech import readnek

def load(path):
    f = readnek(path)
    x = np.concatenate([e.pos[0].ravel() for e in f.elem])
    y = np.concatenate([e.pos[1].ravel() for e in f.elem])
    u = np.concatenate([e.vel[0].ravel() for e in f.elem])
    v = np.concatenate([e.vel[1].ravel() for e in f.elem])
    return x, y, u, v

# (dir, title, [(sigma, St) for modes 1..3])
operators = [
    ("coupled", "coupled",  [(3.920, 0.000), (3.897, 0.216), (3.682, 0.054)]),
    ("quasilaminar", "quasilaminar", [(0.395, 0.206), (0.351, 0.397), (0.302, 0.279)]),
]

fig, axes = plt.subplots(2, 3, figsize=(16, 7.5), constrained_layout=True)
for r, (d, title, eigs) in enumerate(operators):
    for c in range(3):
        ax = axes[r, c]
        x, y, u, v = load(f"{d}/dRe1cyl0.f{c+1:05d}")
        m = (x > -3) & (x < 22) & (np.abs(y) < 6)
        vmax = np.max(np.abs(v[m])) or 1.0
        vn = v[m] / vmax
        ax.tricontourf(x[m], y[m], vn, levels=np.linspace(-1, 1, 23),
                       cmap="RdBu_r", extend="both")
        sig, St = eigs[c]
        ax.set_title(f"{title} — mode {c+1}: $\\sigma$={sig}, St={St}", fontsize=9)
        ax.set_aspect("equal"); ax.set_xlim(-3, 22); ax.set_ylim(-6, 6)
        ax.set_xlabel("x"); ax.set_ylabel("y")
fig.suptitle("Direct-stability eigenmode structure (v', normalized) — cylinder_re1m 2D k-tau RANS, Re=1e6",
             fontsize=12)
fig.savefig("eigenmodes.png", dpi=180)
print("saved eigenmodes.png")
# quick spatial-extent diagnostic: where is |v'| concentrated? (wake vs near-body)
for d, title, _ in operators:
    x, y, u, v = load(f"{d}/dRe1cyl0.f00001")
    a = np.abs(v); thr = a > 0.5 * a.max()
    print(f"{title} mode1: |v'|>50%% peak at  x in [{x[thr].min():.1f},{x[thr].max():.1f}]  "
          f"y in [{y[thr].min():.1f},{y[thr].max():.1f}]  (centroid x={np.average(x,weights=a):.2f})")
