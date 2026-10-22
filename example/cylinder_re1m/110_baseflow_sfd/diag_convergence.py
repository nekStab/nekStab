#!/usr/bin/env python3
"""Honest convergence diagnostic for the cylinder_re1m 2D k-tau RANS SFD base.

Two figures:
  (1) sfd_residual.png  -- SFD convergence residual res=||x_n - x_{n-1}||_knorm
      and its rate vs time, for dt=1e-3 (residu.dat) and dt=5e-3 (residu_dt5e3.dat).
      The base NEVER hit the SFD tolerance gate; this shows where/why it floored.
  (2) base_fields.png   -- all base-flow quantities of base_converged.f00001:
      u, v, p, k(tke), tau, and mu_t ~ k*tau (standard k-tau eddy-viscosity scaling).
      v should be ~0 everywhere for a clean steady base; residual wake structure in v
      is the visual signature of an incompletely-damped (limit-cycling) fixed point.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymech import readnek

# ---------- (1) residual history ----------
r1 = np.loadtxt("residu.dat")        # dt=1e-3
r5 = np.loadtxt("residu_dt5e3.dat")  # dt=5e-3
vel_tol = 1.0e-6                      # [VELOCITY] residualTol

fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
ax[0].semilogy(r5[:, 0], r5[:, 1], lw=0.7, color="goldenrod",
               label=f"dt=5e-3  (floor {r5[:,1].min():.2e})")
ax[0].semilogy(r1[:, 0], r1[:, 1], lw=0.7, color="C0",
               label=f"dt=1e-3  (floor {r1[:,1].min():.2e})")
ax[0].axhline(vel_tol, color="k", ls=":", lw=1, label="velocity solver tol 1e-6")
imin = r1[:, 1].argmin()
ax[0].plot(r1[imin, 0], r1[imin, 1], "v", color="C3", ms=9,
           label=f"min {r1[imin,1]:.2e} @ t={r1[imin,0]:.0f} (then rises)")
ax[0].set_ylabel(r"residual  $\|x_n-x_{n-1}\|_{k\mathrm{-norm}}$")
ax[0].set_title("SFD convergence residual — RANS base NEVER reached the tolerance gate\n"
                "floors ~5e-4, ~500x above the velocity solver tol = NOT a clean fixed point")
ax[0].legend(loc="upper right", fontsize=9); ax[0].grid(alpha=0.3, which="both")

ax[1].plot(r1[:, 0], r1[:, 2], lw=0.5, color="C0", label="rate dt=1e-3")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_ylabel("rate  d(res)/dt"); ax[1].set_xlabel("time")
ax[1].set_ylim(-0.15, 0.15)
ax[1].set_title("rate oscillates around zero (not steadily negative) "
                "= low-amplitude limit cycle, not still converging", fontsize=10)
ax[1].legend(loc="upper right", fontsize=9); ax[1].grid(alpha=0.3)
fig.savefig("sfd_residual.png", dpi=200)
print(f"min res dt1e-3 = {r1[:,1].min():.4e}  final = {r1[-1,1]:.4e}")
print(f"min res dt5e-3 = {r5[:,1].min():.4e}  final = {r5[-1,1]:.4e}")
print("saved sfd_residual.png")

# ---------- (2) base-flow fields ----------
f = readnek("base_converged.f00001")
x = np.concatenate([e.pos[0].ravel() for e in f.elem])
y = np.concatenate([e.pos[1].ravel() for e in f.elem])
u = np.concatenate([e.vel[0].ravel() for e in f.elem])
v = np.concatenate([e.vel[1].ravel() for e in f.elem])
p = np.concatenate([e.pres[0].ravel() for e in f.elem])
k = np.concatenate([e.scal[0].ravel() for e in f.elem])
tau = np.concatenate([e.scal[1].ravel() for e in f.elem])
mut = k * tau  # standard k-tau: nu_t ~ k*tau (tau = 1/omega)

m = (x > -4) & (x < 25) & (np.abs(y) < 8)
panels = [
    ("u  (streamwise)", u, "RdBu_r", None),
    ("v  (transverse) — should be ~0 for steady base", v, "RdBu_r", "sym"),
    ("p  (pressure)", p, "viridis", None),
    ("k  (tke)", k, "inferno", None),
    (r"$\tau$  (=1/$\omega$)", tau, "inferno", None),
    (r"$\mu_t \approx k\,\tau$  (eddy viscosity)", mut, "magma", None),
]
fig, axes = plt.subplots(3, 2, figsize=(15, 11), constrained_layout=True)
for ax, (title, q, cmap, mode) in zip(axes.ravel(), panels):
    qm = q[m]
    if mode == "sym":
        a = np.percentile(np.abs(qm), 99.5)
        lv = np.linspace(-a, a, 25)
        cf = ax.tricontourf(x[m], y[m], qm, levels=lv, cmap=cmap, extend="both")
        title += f"\n(|v|max={np.abs(qm).max():.2e}, rms={np.sqrt(np.mean(qm**2)):.2e})"
    else:
        lo, hi = np.percentile(qm, [1, 99.5])
        cf = ax.tricontourf(x[m], y[m], qm, levels=np.linspace(lo, hi, 25),
                            cmap=cmap, extend="both")
    fig.colorbar(cf, ax=ax, shrink=0.8)
    ax.set_title(title, fontsize=10); ax.set_aspect("equal")
    ax.set_xlim(-4, 25); ax.set_ylim(-8, 8)
fig.suptitle("cylinder_re1m base_converged.f00001 — all base-flow quantities "
             "(2D k-tau RANS, Re=1e6, SFD)", fontsize=13)
fig.savefig("base_fields.png", dpi=170)
print(f"v: |max|={np.abs(v[m]).max():.3e}  rms={np.sqrt(np.mean(v[m]**2)):.3e}")
print(f"u range [{u[m].min():.3f},{u[m].max():.3f}]  k range [{k[m].min():.2e},{k[m].max():.2e}]")
print("saved base_fields.png")
