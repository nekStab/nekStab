#!/usr/bin/env python3
"""Figures for the cylinder Re=1e6 moments check.

Two states, same mesh: the SFD base and the unforced steady wake.
Then R_ij = E(u_i u_j) - E(u_i)E(u_j) from the 20-time-unit window.
This is not the k-tau closure.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from check_force import EM2, EUV, BFA, WAKE, load

ROOT = Path(__file__).resolve().parent
SFD = ROOT.parent / "110_baseflow_sfd" / "base_converged.f00001"
XI, YI = np.meshgrid(np.linspace(-2.0, 22.0, 480), np.linspace(-6.0, 6.0, 280))
HOLE = XI**2 + YI**2 < 0.55**2


def on_grid(x, y, q):
    m = (x > -3.0) & (x < 23.0) & (np.abs(y) < 7.0)
    zi = griddata(np.column_stack((x[m], y[m])), q[m], (XI, YI), method="linear")
    zi = np.ma.masked_where(HOLE | ~np.isfinite(zi), zi)
    return zi


def panel(ax, x, y, q, levels, cmap, title):
    zi = on_grid(x, y, q)
    cf = ax.contourf(XI, YI, zi, levels=levels, cmap=cmap, extend="both")
    ax.add_patch(plt.Circle((0, 0), 0.5, color="k", zorder=3))
    ax.set_aspect("equal")
    ax.set_xlim(-2, 22)
    ax.set_ylim(-6, 6)
    ax.set_title(title, fontsize=10)
    return cf


def main():
    sfd = load(SFD)
    wake = load(WAKE)
    bfa = load(BFA)
    em2 = load(EM2)
    euv = load(EUV)
    ruu = em2["u"] - bfa["u"] ** 2
    rvv = em2["v"] - bfa["v"] ** 2
    ruv = euv["u"] - bfa["u"] * bfa["v"]

    u_levels = np.linspace(-0.05, 1.15, 25)
    v_levels = np.linspace(-0.4, 0.4, 21)
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 6.2), constrained_layout=True)
    panel(axes[0, 0], sfd["x"], sfd["y"], sfd["u"], u_levels, "RdBu_r", "SFD base: u")
    panel(axes[0, 1], sfd["x"], sfd["y"], sfd["v"], v_levels, "RdBu_r", "SFD base: v")
    cf_u = panel(axes[1, 0], wake["x"], wake["y"], wake["u"], u_levels, "RdBu_r", "unforced wake: u")
    cf_v = panel(axes[1, 1], wake["x"], wake["y"], wake["v"], v_levels, "RdBu_r", "unforced wake: v")
    fig.colorbar(cf_u, ax=axes[:, 0], shrink=0.85, label="u")
    fig.colorbar(cf_v, ax=axes[:, 1], shrink=0.85, label="v")
    fig.suptitle("Same mesh. SFD base and the unforced steady wake. Not a shedding street.", fontsize=12)
    fig.savefig(ROOT / "states.png", dpi=160)
    print("saved states.png")

    m = (bfa["x"] > 1.0) & (bfa["x"] < 12.0) & (np.abs(bfa["y"]) < 2.0)
    peaks = [float(np.max(np.abs(q[m]))) for q in (ruu, rvv, ruv)]
    lim = max(peaks)
    r_levels = np.linspace(-lim, lim, 21)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8), constrained_layout=True)
    names = (r"$R_{uu}$", r"$R_{vv}$", r"$R_{uv}$")
    fields = (ruu, rvv, ruv)
    for ax, name, q, peak in zip(axes, names, fields, peaks):
        cf = panel(ax, bfa["x"], bfa["y"], q, r_levels, "RdBu_r", f"{name}   max |R|={peak:.2e}")
    fig.colorbar(cf, ax=axes, shrink=0.85, label=r"$R_{ij}$")
    fig.suptitle(
        r"20 t.u. window. $R_{ij}=E(u_i u_j)-E(u_i)E(u_j)$. All three are $\sim 10^{-4}$, not $O(1)$.",
        fontsize=11,
    )
    fig.savefig(ROOT / "stress.png", dpi=160)
    print("saved stress.png")


if __name__ == "__main__":
    main()
