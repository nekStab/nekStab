#!/usr/bin/env python3
"""Plot thermosyphon wavemaker and mode fields."""
from pathlib import Path
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk


CASE_DIR = Path(__file__).resolve().parent
R_INNER = 1.0
R_OUTER = 2.0


def _read_required(pattern):
    matches = sorted(CASE_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"no files match {pattern!r}")
    return matches[0]


def _field_component(path, preferred):
    x, y, fields, _time = nk.read_field(path)
    for key in preferred:
        if key in fields:
            return x, y, fields[key]
    raise KeyError(f"{path.name} has none of {preferred}")


def _plot(kind, x, y, values, output):
    fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.75, nk.COL_WIDTH * 0.75))
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
    nk.inset_colorbar(ax, cf, orientation="horizontal", width="50%",
                      height="5%", loc=9, ticks=ticks,
                      tick_labels=[f"{tick:g}" for tick in ticks])
    ax.set_aspect("equal")
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    nk.add_annulus_patches(ax, R_INNER, R_OUTER)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output.name}")


def _write_metrics(x, y, wavemaker):
    idx = int(np.nanargmax(wavemaker))
    metrics = {
        "wm_max": float(wavemaker[idx]),
        "wm_x_at_max": float(x[idx]),
        "wm_y_at_max": float(y[idx]),
        "wm_mean": float(np.nanmean(wavemaker)),
    }
    with (CASE_DIR / "wm_metrics.dat").open("w") as fh:
        fh.write("# wm_max wm_x_at_max wm_y_at_max wm_mean\n")
        fh.write("{wm_max:.16e} {wm_x_at_max:.16e} "
                 "{wm_y_at_max:.16e} {wm_mean:.16e}\n".format(**metrics))
    return metrics


def main():
    nk.configure_style()
    x, y, wavemaker = _field_component(_read_required("wm_tsyphon0.f00001"), ["t"])
    xd, yd, direct = _field_component(_read_required("dRetsyphon0.f00001"), ["t", "vy", "vx"])
    xa, ya, adjoint = _field_component(_read_required("aRetsyphon0.f00002"), ["t", "vy", "vx"])

    _plot("wavemaker", x, y, wavemaker, CASE_DIR / "plot_wavemaker.png")
    _plot("direct", xd, yd, direct, CASE_DIR / "plot_direct_mode.png")
    _plot("adjoint", xa, ya, adjoint, CASE_DIR / "plot_adjoint_mode.png")
    metrics = _write_metrics(x, y, wavemaker)

    ref_dir = CASE_DIR / "ref"
    ref_dir.mkdir(exist_ok=True)
    for png in CASE_DIR.glob("plot_*.png"):
        shutil.copy2(png, ref_dir / png.name)
    print("Wrote wm_metrics.dat")
    for name, value in metrics.items():
        print(f"{name} {value:.16e}")


if __name__ == "__main__":
    main()
