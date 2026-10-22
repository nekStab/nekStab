#!/usr/bin/env python3
"""Plot flip-flop Floquet post-processing fields."""
from pathlib import Path
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk


CASE_DIR = Path(__file__).resolve().parent
GAP = 0.7
RADIUS = 0.5


def _read_required(pattern):
    matches = sorted(CASE_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"no files match {pattern!r}")
    return matches[0]


def _velocity(path):
    x, y, fields, _time = nk.read_field(path)
    if "vx" not in fields or "vy" not in fields:
        raise KeyError(f"{path.name} does not contain vx/vy")
    return x, y, fields["vx"], fields["vy"]


def _production_density():
    """Orbit-averaged production density of mode 1, the sum of the stored components.

    F012cyl0.f00001 holds the components for the x-derivative of the base flow, F012cyl0.f00002 those for the
    y-derivative (the z components of a 2D run are not written). Their volume integral is P_bar of the log.
    """
    total = None
    for number in (1, 2):
        x, y, fields, _time = nk.read_field(_read_required(f"F012cyl0.f0000{number}"))
        part = sum(fields[name] for name in ("vx", "vy", "vz") if name in fields)
        total = part if total is None else total + part
    return x, y, total


def _plot(kind, x, y, values, output):
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.35, nk.COL_WIDTH * 0.55))
    triang = nk.make_triangulation(x, y)
    if kind == "wavemaker":
        vmax = float(np.nanpercentile(values, 99.5))
        cf = nk.tricontourf(ax, triang, values, levels=257, cmap="hot_r",
                            vmin=0.0, vmax=vmax, extend="max")
        ticks = [0.0, round(vmax, 3)]
    else:
        bound = float(np.nanpercentile(np.abs(values), 99.0))
        cf = nk.tricontourf(ax, triang, values, levels=257, cmap="RdBu_r",
                            vmin=-bound, vmax=bound, extend="both")
        ticks = [round(-bound, 3), round(bound, 3)]
    nk.inset_colorbar(ax, cf, orientation="horizontal", width="48%",
                      height="5%", loc=9, ticks=ticks,
                      tick_labels=[f"{tick:g}" for tick in ticks])
    nk.add_dual_cylinder_patches(ax, gap=GAP, radius=RADIUS)
    ax.set_aspect("equal")
    ax.set_xlim(-3, 14)
    ax.set_ylim(-3, 3)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output.name}")


def main():
    nk.configure_style()

    xd, yd, dvx, dvy = _velocity(_read_required("dRe2cyl0.f00001"))
    xa, ya, avx, avy = _velocity(_read_required("aRe2cyl0.f00001"))
    xb, yb, production = _production_density()

    direct_mag = np.sqrt(dvx**2 + dvy**2)
    adjoint_mag = np.sqrt(avx**2 + avy**2)
    overlap = direct_mag * adjoint_mag
    _plot("wavemaker", xd, yd, overlap, CASE_DIR / "plot_wavemaker.png")
    _plot("direct", xd, yd, dvy, CASE_DIR / "plot_direct_mode.png")
    _plot("adjoint", xa, ya, avy, CASE_DIR / "plot_adjoint_mode.png")
    _plot("production", xb, yb, production, CASE_DIR / "plot_budget_field.png")

    ref_dir = CASE_DIR / "ref"
    ref_dir.mkdir(exist_ok=True)
    for png in CASE_DIR.glob("plot_*.png"):
        shutil.copy2(png, ref_dir / png.name)


if __name__ == "__main__":
    main()
