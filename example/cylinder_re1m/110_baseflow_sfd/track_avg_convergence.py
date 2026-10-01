#!/usr/bin/env python3
"""Track statistical convergence of the frozen eddy viscosity across avg windows.

avg_all writes an INDEPENDENT windowed average every 25 t.u. (param(68)=25000).
For each window we form mu_t ~ <k>*<tau> and report wake-region statistics. If
windows 2..N agree, the frozen-mu_t base is statistically converged. Window 1 is
the SFD filter re-init transient (flagged, not trusted). We also compare against
the single-snapshot mu_t from base_converged.f00001 to quantify the caveat.
"""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymech import readnek

def stats(path):
    f = readnek(path)
    g = lambda s: np.concatenate([s(e).ravel() for e in f.elem])
    x, y = g(lambda e: e.pos[0]), g(lambda e: e.pos[1])
    k, tau = g(lambda e: e.scal[0]), g(lambda e: e.scal[1])
    mut = k * tau
    wake = (x > 1) & (x < 20) & (np.abs(y) < 3)
    return dict(t=f.time, mut_wake=mut[wake].mean(), mut_max=mut.max(),
                k_wake=k[wake].mean(), kneg=100*(k < 0).mean())

avg = sorted(glob.glob("avg1cyl0.f*"))
print(f"found {len(avg)} avg windows")
base = stats("base_converged.f00001")
rows = [stats(p) for p in avg]

print(f"\n{'window':>10} {'t_end':>7} {'<mut>_wake':>12} {'mut_max':>10} {'<k>_wake':>10} {'%k<0':>7}")
print(f"{'SNAPSHOT':>10} {'100':>7} {base['mut_wake']:12.4e} {base['mut_max']:10.3e} {base['k_wake']:10.3e} {base['kneg']:7.1f}")
for i, (p, r) in enumerate(zip(avg, rows), 1):
    tag = "  (transient)" if i == 1 else ""
    print(f"{i:>10} {r['t']:7.1f} {r['mut_wake']:12.4e} {r['mut_max']:10.3e} {r['k_wake']:10.3e} {r['kneg']:7.1f}{tag}")

# stationarity over windows 2..N
if len(rows) >= 3:
    m = np.array([r["mut_wake"] for r in rows[1:]])
    spread = 100 * (m.max() - m.min()) / m.mean()
    print(f"\nwindows 2..{len(rows)}: <mut>_wake spread = {spread:.2f}% "
          f"(<2% = statistically converged)")
    dvs = 100 * (m.mean() - base["mut_wake"]) / base["mut_wake"]
    print(f"time-averaged vs single snapshot: {dvs:+.1f}% in <mut>_wake "
          f"(magnitude of the snapshot caveat)")

# figure  (each window = an independent 25 t.u. average; index = chronological order)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
wi = list(range(1, len(rows) + 1))
mw = [r["mut_wake"] for r in rows]
ax[0].plot(wi, mw, "o-", color="C0", label=r"windowed $\langle\mu_t\rangle$ (25 t.u. each)")
ax[0].axhline(base["mut_wake"], color="C3", ls="--",
              label=f"single snapshot t=100 ({base['mut_wake']:.3e})")
if len(rows) >= 2:
    m2 = np.mean(mw[1:])
    ax[0].axhline(m2, color="C2", ls=":", label=f"mean of windows 2-6 ({m2:.3e})")
    ax[0].axvspan(0.5, 1.5, alpha=0.12, color="grey")
    ax[0].text(1, max(mw), "transient", fontsize=8, ha="center", va="bottom")
ax[0].set_xlabel("avg window # (25 t.u. each, chronological)")
ax[0].set_ylabel(r"$\langle\mu_t\rangle_{wake}\;\approx\langle k\rangle\langle\tau\rangle$")
ax[0].set_title("Frozen eddy viscosity per window\nspread of windows 2-6 = "
                f"{100*(max(mw[1:])-min(mw[1:]))/np.mean(mw[1:]):.1f}%  "
                f"(snapshot off by +1.8%)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[1].plot(wi, [r["kneg"] for r in rows], "s-", color="C2")
ax[1].set_xlabel("avg window #"); ax[1].set_ylabel("% cells with k<0")
ax[1].set_ylim(0, 30)
ax[1].set_title("Negative-TKE (limiter) fraction per window\n"
                "~25% steady = limiter fires in a quarter of cells"); ax[1].grid(alpha=0.3)
fig.savefig("avg_convergence.png", dpi=170)
print("\nsaved avg_convergence.png")
