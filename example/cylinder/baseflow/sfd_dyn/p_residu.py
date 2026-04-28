#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

params = {
    "text.usetex": False,
    "font.size": 8,
    "legend.fontsize": 8,
    "legend.handlelength": 2.5,
}
plt.rcParams.update(params)
try:
    plt.style.use("seaborn-v0_8-white")
except OSError:
    plt.style.use("default")

FORMT = "png"
ADJUST = "tight"
QUAL = 500
FIG_WIDTH = 6.0
FIG_HEIGHT = 2.45


def load_table(filename):
    path = Path(filename)
    print(f"Reading {path.name}")
    return np.genfromtxt(path, comments="#")


def as_table(data):
    if data.ndim == 1:
        return data.reshape(1, -1)
    return data


if __name__ == "__main__":
    residu = as_table(load_table("residu.dat"))
    baseline_path = Path("../sfd/residu.dat")
    baseline = as_table(load_table(baseline_path)) if baseline_path.exists() else None
    dyn_tol_path = Path("dyn_tol.dat")
    dyn_tol = as_table(load_table("dyn_tol.dat")) if dyn_tol_path.exists() else None

    fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH, FIG_HEIGHT))

    ax_res = axes[0]
    ax_res.set_yscale("log")
    ax_res.set_xlabel(r"$t$")
    ax_res.set_ylabel(r"$\|r\|$")
    ax_res.axhline(y=1e-9, lw=0.3, c="k", ls="dotted", label=r"$10^{-9}$")
    if baseline is not None:
        ax_res.plot(baseline[:, 0], baseline[:, 1], c="0.45", lw=0.5,
                    label="SFD fixed")
    ax_res.plot(residu[:, 0], residu[:, 1], c="m", lw=0.45, label="SFD dyn")
    ax_res.axhline(y=np.nanmin(residu[:, 1]), c="m", lw=0.3, ls="--")
    ax_res.axvline(x=np.nanmax(residu[:, 0]), c="m", lw=0.3, ls="--")
    ax_res.legend(loc="best", fontsize=6)

    ax_tol = axes[1]
    ax_tol.set_yscale("log")
    ax_tol.set_xlabel(r"$t$")
    ax_tol.set_ylabel("solver tolerance")
    ax_tol.plot(residu[:, 0], residu[:, 3], c="0.35", lw=0.45, label="current tol")

    if dyn_tol is not None and np.size(dyn_tol) > 0:
        ax_tol.plot(dyn_tol[:, 0], dyn_tol[:, 3], "o-", ms=2.0, lw=0.45,
                    c="C1", label="requested")
        ax_tol.plot(dyn_tol[:, 0], dyn_tol[:, 4], "o-", ms=2.0, lw=0.45,
                    c="C2", label="used")
        if np.any(dyn_tol[:, 5] > 0.0):
            ax_tol.plot(dyn_tol[:, 0], dyn_tol[:, 5], "--", lw=0.45,
                        c="C3", label="cap")

    ax_tol.legend(loc="best", fontsize=6)

    fname = "residu." + FORMT
    plt.savefig(fname, format=FORMT, dpi=QUAL, bbox_inches=ADJUST)
    print(f"Saving {fname}")
    plt.close()
    print("------------------------------------------")
