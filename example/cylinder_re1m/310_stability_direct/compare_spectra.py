#!/usr/bin/env python3
"""Coupled (full k-tau Frechet) vs quasilaminar (prescribed mu_t) stability spectra.

Spectre_NSd.dat columns: sigma (growth rate), omega (frequency), residual, idx.
St = omega/(2 pi). Re=1e6, 2D, converged steady RANS base (SFD).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = np.loadtxt("coupled/Spectre_NSd.dat")
B = np.loadtxt("quasilaminar/Spectre_NSd.dat")
StA, sigA = A[:, 1] / (2*np.pi), A[:, 0]
StB, sigB = B[:, 1] / (2*np.pi), B[:, 0]
St_obs = 0.192

fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

# full range — shows the coupled operator's contaminated high-sigma cluster
ax[0].scatter(StA, sigA, s=28, c="C3", marker="o", label="coupled (full k-tau)")
ax[0].scatter(StB, sigB, s=28, c="C0", marker="s", label="quasilaminar (prescribed mu_t)")
ax[0].axhline(0, color="k", lw=0.8)
ax[0].axvline(St_obs, ls="--", color="grey", lw=1, label=f"observed St={St_obs}")
ax[0].set_xlabel("St = $\\omega/2\\pi$"); ax[0].set_ylabel("$\\sigma$ (growth rate)")
ax[0].set_title("Full spectrum: coupled leading cluster at $\\sigma\\approx3.9$ (k-tau-contaminated)")
ax[0].grid(alpha=0.3); ax[0].legend(loc="upper right", fontsize=8)

# zoom to the physical range (repo convention sigma in [-0.4,0.3])
ax[1].scatter(StA, sigA, s=36, c="C3", marker="o", label="coupled")
ax[1].scatter(StB, sigB, s=36, c="C0", marker="s", label="quasilaminar")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].axvline(St_obs, ls="--", color="grey", lw=1, label=f"observed St={St_obs}")
ax[1].set_xlim(-0.5, 0.5); ax[1].set_ylim(-0.5, 0.6)
ax[1].set_xlabel("St = $\\omega/2\\pi$"); ax[1].set_ylabel("$\\sigma$")
ax[1].set_title("Physical range: quasilaminar leading mode $\\sigma$=0.40 @ St=0.206")
ax[1].grid(alpha=0.3); ax[1].legend(loc="upper right", fontsize=8)

fig.suptitle("cylinder_re1m 2D k-tau RANS — direct stability, coupled vs quasilaminar", fontsize=11)
fig.savefig("compare_spectra.png", dpi=300)

def lead(St, sig):
    i = np.argmax(sig)
    return sig[i], St[i]
sA, fA = lead(StA, sigA); sB, fB = lead(StB, sigB)
print("coupled leading: sigma=%.4f  St=%.4f  (#modes=%d, sigma range [%.3f,%.3f])"
      % (sA, fA, len(sigA), sigA.min(), sigA.max()))
print("quasilaminar leading: sigma=%.4f  St=%.4f  (#modes=%d, sigma range [%.3f,%.3f])"
      % (sB, fB, len(sigB), sigB.min(), sigB.max()))
print("coupled has a mode near (St=0.2, sigma<0.6)? ",
      np.any((np.abs(StA)-0.2 < 0.05) & (sigA < 0.6) & (sigA > 0.1)))
print("saved compare_spectra.png")
