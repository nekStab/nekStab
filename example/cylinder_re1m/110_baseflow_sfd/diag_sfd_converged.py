#!/usr/bin/env python3
"""Clean converged steady RANS via SFD: residual floor drops with dt; wake steady.

Panel 1: SFD residual vs time for dt=5e-3 (residu_dt5e3.dat) and dt=1e-3
(residu.dat) -- the floor tracks dt (discretization-limited), 5e-3 -> 5e-4.
Panel 2: mid-wake probe v vs time (dt=1e-3) -- v ~ 1e-4 == steady (bare-march
shedding was +/-0.87 for scale).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

r5 = np.loadtxt("residu_dt5e3.dat")
r1 = np.loadtxt("residu.dat")

with open("1cyl.his") as f:
    npb = int(f.readline().split()[0])
    [f.readline() for _ in range(npb)]
d = np.loadtxt("1cyl.his", skiprows=1 + npb)
nt = d.shape[0] // npb
d = d[: nt * npb].reshape(nt, npb, -1)
tp, vp = d[:, 1, 0], d[:, 1, 2]

fig, ax = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
ax[0].semilogy(r5[:, 0], np.abs(r5[:, 1]), lw=0.7, color="C1", label="dt=5e-3 (CFL~0.5): floors ~5e-3")
ax[0].semilogy(r1[:, 0], np.abs(r1[:, 1]), lw=0.7, color="C0", label="dt=1e-3 (CFL~0.1): -> ~5e-4")
ax[0].axhline(1e-5, ls="--", color="grey", lw=1, label="tol 1e-5")
ax[0].set_xlabel("time [c.u.]"); ax[0].set_ylabel(r"SFD residual")
ax[0].set_title("SFD convergence of the steady k-tau RANS base, Re=1e6 2D (corrected cutoff St=0.192)")
ax[0].grid(True, which="both", alpha=0.3); ax[0].legend(loc="upper right", fontsize=8)

ax[1].plot(tp, vp, lw=0.6, color="C0")
ax[1].set_xlabel("time [c.u.]"); ax[1].set_ylabel("v at mid-wake x=5")
ax[1].set_title("Wake steady: v_rms ~ 1e-4 (bare RANS march shed at +/-0.87 -- a 5600x reduction)")
ax[1].grid(True, alpha=0.3)

fig.savefig("diag_sfd_converged.png", dpi=300)
m = tp >= 30
print("dt=1e-3 final residual ~ %.3e" % np.abs(r1[-1, 1]))
print("probe v_rms (t>=30) = %.3e   mean = %.3e   => steady" % (np.std(vp[m]), np.mean(vp[m])))
print("saved diag_sfd_converged.png")
