#!/usr/bin/env python3
"""Build a PROPER long-time average base from the avg_all windows + prove convergence.

avg_all writes an independent 25 t.u. windowed average each dump (it resets atime).
Because the windows are equal-duration at constant dt, the mean of windows 2..n is
EXACTLY the cumulative time-average over that span -- so the running cumulative
<mu_t> vs n is the publishable convergence curve, and the final cumulative field is
the converged frozen-mu_t base. Window 1 is discarded (SFD filter re-init transient).

Outputs:
  cumulative_convergence.png  -- running cumulative <mu_t>_wake vs averaging time
  avg_base.f00001             -- the converged averaged field (stage as BF_ for quasilaminar)
"""
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymech import readnek, writenek

WIN_TU = 25.0  # t.u. per window (param(68)=25000, dt=1e-3)
files = sorted(glob.glob("avg1cyl0.f*"))
print(f"{len(files)} windows; discarding window 1 (transient)")
clean = files[1:]                       # drop transient
fields = []                             # skip any file still mid-write (concurrent job)
for p in clean:
    try:
        fields.append(readnek(p))
    except Exception as e:
        print(f"skip {p} (incomplete: {e})")
print(f"using {len(fields)} complete windows")

def cat(f, sel):
    return np.concatenate([sel(e).ravel() for e in f.elem])

x = cat(fields[0], lambda e: e.pos[0]); y = cat(fields[0], lambda e: e.pos[1])
wake = (x > 1) & (x < 20) & (np.abs(y) < 3)

# per-window mu_t = <k>*<tau>
kk = np.array([cat(f, lambda e: e.scal[0]) for f in fields])
tt = np.array([cat(f, lambda e: e.scal[1]) for f in fields])
mut_win = kk * tt

# drop any window read partially (concurrent write) -> wake-mean far from median
wm = (mut_win[:, wake]).mean(axis=1)
med = np.median(wm)
keep = np.abs(wm - med) < 0.05 * med
if (~keep).any():
    print(f"dropped {int((~keep).sum())} partial/outlier window(s)")
fields = [f for f, k in zip(fields, keep) if k]
mut_win = mut_win[keep]
n_keep = len(fields)

# running cumulative: C_n = mean of the first n kept windows
n_win = n_keep
t_end = np.array([(i + 1) * WIN_TU for i in range(n_win)])   # cumulative averaging duration
cum_mut_wake = np.array([mut_win[: i + 1].mean(axis=0)[wake].mean() for i in range(n_win)])

# convergence metric: spread over the last third of the running curve
tail = cum_mut_wake[max(1, 2 * n_win // 3):]
spread = 100 * (tail.max() - tail.min()) / tail.mean()
final = cum_mut_wake[-1]
print(f"final cumulative <mu_t>_wake = {final:.5e} over {t_end[-1]-WIN_TU:.0f} t.u.")
print(f"running-cumulative spread over last third = {spread:.2f}%  (<1% = publishable)")

# figure
fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
ax.plot(t_end, cum_mut_wake, "o-", color="C0", ms=4,
        label="running cumulative average")
ax.axhline(final, color="C2", ls=":", label=f"final = {final:.4e}")
ax.fill_between(t_end[max(1, 2*n_win//3):], tail.min(), tail.max(),
                alpha=0.15, color="C2")
ax.set_xlabel("averaging window end-time  [t.u.]")
ax.set_ylabel(r"cumulative $\langle\mu_t\rangle_{wake}$")
ax.set_title("Proper long-time average of the frozen eddy viscosity\n"
             f"running cumulative flattens to {spread:.2f}% over the last third "
             f"({t_end[-1]-WIN_TU:.0f} t.u. total)")
ax.legend(); ax.grid(alpha=0.3)
fig.savefig("cumulative_convergence.png", dpi=170)
print("saved cumulative_convergence.png")

# build the averaged field: overwrite window-2 handle with the mean across clean windows
out = fields[0]
for ie, e in enumerate(out.elem):
    for c in range(e.vel.shape[0]):
        e.vel[c] = np.mean([fields[w].elem[ie].vel[c] for w in range(n_win)], axis=0)
    e.pres[0] = np.mean([fields[w].elem[ie].pres[0] for w in range(n_win)], axis=0)
    for s in range(e.scal.shape[0]):
        e.scal[s] = np.mean([fields[w].elem[ie].scal[s] for w in range(n_win)], axis=0)
out.time = t_end[-1] - WIN_TU
writenek("avg_base.f00001", out)
print("wrote avg_base.f00001 (converged averaged base; stage as BF_1cyl0.f00001)")
