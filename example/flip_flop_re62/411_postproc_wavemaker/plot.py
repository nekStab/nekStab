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


def _budget_field():
    for pattern in ("F022cyl0.f00001", "F012cyl0.f00001"):
        path = _read_required(pattern)
        x, y, fields, _time = nk.read_field(path)
        if "vx" in fields and "vy" in fields:
            return x, y, np.sqrt(fields["vx"]**2 + fields["vy"]**2)
        for key in ("t", "p"):
            if key in fields:
                return x, y, np.abs(fields[key])
    raise KeyError("no budget field with vx/vy, t, or p found")


def _plot(kind, x, y, values, output):
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 1.35, nk.COL_WIDTH * 0.55))
    triang = nk.make_triangulation(x, y)
    if kind in {"wavemaker", "budget"}:
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


def _write_metrics(x, y, values):
    idx = int(np.nanargmax(values))
    metrics = {
        "wm_max": float(values[idx]),
        "wm_x_at_max": float(x[idx]),
        "wm_y_at_max": float(y[idx]),
        "wm_mean": float(np.nanmean(values)),
    }
    with (CASE_DIR / "wm_metrics.dat").open("w") as fh:
        fh.write("# wm_max wm_x_at_max wm_y_at_max wm_mean\n")
        fh.write("{wm_max:.16e} {wm_x_at_max:.16e} "
                 "{wm_y_at_max:.16e} {wm_mean:.16e}\n".format(**metrics))
    return metrics


def main():
    nk.configure_style()

    xd, yd, dvx, dvy = _velocity(_read_required("dRe2cyl0.f00001"))
    xa, ya, avx, avy = _velocity(_read_required("aRe2cyl0.f00001"))
    xb, yb, budget_mag = _budget_field()

    direct_mag = np.sqrt(dvx**2 + dvy**2)
    adjoint_mag = np.sqrt(avx**2 + avy**2)
    overlap = direct_mag * adjoint_mag
    _plot("wavemaker", xd, yd, overlap, CASE_DIR / "plot_wavemaker.png")
    _plot("direct", xd, yd, dvy, CASE_DIR / "plot_direct_mode.png")
    _plot("adjoint", xa, ya, avy, CASE_DIR / "plot_adjoint_mode.png")
    _plot("budget", xb, yb, budget_mag, CASE_DIR / "plot_budget_field.png")
    metrics = _write_metrics(xd, yd, overlap)

    ref_dir = CASE_DIR / "ref"
    ref_dir.mkdir(exist_ok=True)
    for png in CASE_DIR.glob("plot_*.png"):
        shutil.copy2(png, ref_dir / png.name)
    print("Wrote wm_metrics.dat")
    for name, value in metrics.items():
        print(f"{name} {value:.16e}")


if __name__ == "__main__":
    main()
