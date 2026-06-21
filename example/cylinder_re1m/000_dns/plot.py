#!/usr/bin/env python3
"""cylinder_re1m RANS (k-tau) forward-solve diagnostic.

Reads the wake-probe history (1cyl.his) and plots the near-wake velocity and
turbulent kinetic energy vs time, showing the relaxation from the inflow IC
onto the (unsteady, shedding) RANS state. Saves plot.png (gallery figure).
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt

HIS = "1cyl.his"
NP = 3  # probes: (2,0,10), (5,0,10), (10,0,10)
PROBE_LABELS = ["x=2 (near wake)", "x=5", "x=10 (far wake)"]

# Parse: line 1 = npoints, next NP lines = coords, then groups of NP data rows.
with open(HIS) as f:
    lines = [ln.split() for ln in f if ln.strip()]
npoints = int(lines[0][0])
data_rows = lines[1 + npoints:]  # skip count + coord lines

# columns: t, vx, vy, vz, pr, temp, tke, tau
series = {p: {"t": [], "vx": [], "vy": [], "tke": []} for p in range(npoints)}
for i in range(0, len(data_rows) - npoints + 1, npoints):
    for p in range(npoints):
        r = data_rows[i + p]
        if len(r) < 8:
            continue
        series[p]["t"].append(float(r[0]))
        series[p]["vx"].append(float(r[1]))
        series[p]["vy"].append(float(r[2]))
        series[p]["tke"].append(float(r[6]))

fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
for p in range(npoints):
    ax[0].plot(series[p]["t"], series[p]["vx"], lw=0.8, label=PROBE_LABELS[p])
    ax[1].plot(series[p]["t"], series[p]["vy"], lw=0.8, label=PROBE_LABELS[p])
    ax[2].plot(series[p]["t"], series[p]["tke"], lw=0.8, label=PROBE_LABELS[p])
ax[0].set_ylabel(r"$u_x$  [-]")
ax[1].set_ylabel(r"$u_y$  [-]")
ax[2].set_ylabel(r"tke $k$  [-]")
ax[2].set_xlabel(r"time  [convective units]")
ax[0].set_title(r"cylinder_re1m k-$\tau$ RANS forward solve (Re=40000): wake-probe history")
for a in ax:
    a.grid(True, alpha=0.3)
    a.legend(fontsize=8, loc="best")
fig.tight_layout()
fig.savefig("plot.png", dpi=300, bbox_inches="tight")
print("wrote plot.png")
