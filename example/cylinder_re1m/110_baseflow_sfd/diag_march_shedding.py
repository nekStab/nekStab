#!/usr/bin/env python3
"""Forward RANS march (dt=5e-3, no SFD): show the wake sheds strongly => URANS.

Columns in 1cyl.his (2D, no w): time, u, v, p, temp, tke, tau. Probe order is
x=2 (near), x=5 (mid), x=10 (far). We plot the mid-wake probe (index 1).
The SFD run held this same probe's v at +/-0.006 (near-steady); the bare march
lets it grow to its saturated shedding amplitude.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

with open("1cyl.his") as f:
    nprobe = int(f.readline().split()[0])
    [f.readline() for _ in range(nprobe)]
data = np.loadtxt("1cyl.his", skiprows=1 + nprobe)
ntime = data.shape[0] // nprobe
data = data[: ntime * nprobe].reshape(ntime, nprobe, -1)
pj = 1  # mid-wake x=5
t, u, v = data[:, pj, 0], data[:, pj, 1], data[:, pj, 2]

fig, ax = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True, sharex=True)
ax[0].plot(t, u, lw=0.6, color="C0")
ax[0].set_ylabel("u  (mid-wake x=5)")
ax[0].set_title("Bare forward k-tau RANS march, Re=1e6, dt=5e-3 (no SFD) — wake sheds (URANS)")
ax[0].grid(alpha=0.3)

ax[1].plot(t, v, lw=0.6, color="C3", label="bare march (URANS)")
ax[1].axhline(0.006, ls="--", lw=1, color="grey")
ax[1].axhline(-0.006, ls="--", lw=1, color="grey",
              label="SFD-suppressed level (+/-0.006 = steady RANS mean)")
ax[1].set_ylabel("v  (mid-wake x=5)")
ax[1].set_xlabel("time [convective units]")
ax[1].grid(alpha=0.3)
ax[1].legend(loc="upper left", fontsize=8)

fig.savefig("diag_march_shedding.png", dpi=300)
m = (t >= 40)
print("over t>=40:  v in [%.3f, %.3f]   u in [%.3f, %.3f]" %
      (v[m].min(), v[m].max(), u[m].min(), u[m].max()))
print("v_rms (t>=40) = %.3f  vs SFD-suppressed ~0.006  => factor ~%.0f" %
      (np.std(v[m]), np.std(v[m]) / 0.006))
print("saved diag_march_shedding.png")
