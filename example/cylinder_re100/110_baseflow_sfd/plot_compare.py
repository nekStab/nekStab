#!/usr/bin/env python3
"""plot_compare.py: residual of the six SFD variants of the Re = 100 cylinder wake against time and against wall time.

plot_compare.png
    (a) residual against the time of the simulation, with the two stop gates (1e-10 and 1e-6) as guide lines;
    (b) the same residual against the wall time of the run, read from the `Step` lines of each logfile.
The legend gives the wall time of each run.

USAGE: python plot_compare.py [PATTERN]
    PATTERN is the folder of a variant with {v} for its name; the default is ./{v}. The folder holds residu.dat and logfile.
"""

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
VARIANTS = [
    # name, line width, marker, colour; the thickest line is drawn first, at the back
    ("casacuberta_dyn_oifs", 4.4, "v", "C5"),
    ("akervik_dyn_oifs", 3.6, "^", "C4"),
    ("casacuberta_dyn", 3.0, "D", "C3"),
    ("casacuberta", 2.4, "s", "C2"),
    ("akervik_dyn", 1.8, "o", "C1"),
    ("akervik", 1.2, "x", "C0"),
]
STEP = re.compile(r"^Step\s+\d+, t=\s*\S+, DT=\s*\S+, C=\s*\S+\s+(\S+)\s+\S+")
ELAPSED = re.compile(r"total elapsed time\s*:\s*(\S+)")
GATES = [(1e-10, "gate 1e-10 (BDF3 variants)"), (1e-6, "gate 1e-6 (OIFS variants)")]


def spread(n_rows, n_marks):
    """Indices of n_marks samples spread evenly over n_rows rows."""
    return sorted({int(i) for i in np.linspace(0, n_rows - 1, n_marks)})


def load(folder):
    """Time, residual, wall time at each time step and the total elapsed time of one run."""
    data = np.loadtxt(folder / "residu.dat")
    t, res = data[:, 0], data[:, 1]
    wall, elapsed = np.zeros(len(t)), None
    steps = []
    for line in (folder / "logfile").read_text(errors="replace").splitlines():
        m = STEP.match(line)
        if m:
            steps.append(float(m.group(1)))
        elif elapsed is None:
            m = ELAPSED.search(line)
            if m:
                elapsed = float(m.group(1))
    if len(steps) >= len(t):
        wall = np.array(steps[: len(t)])
    else:
        wall = np.linspace(elapsed / len(t), elapsed, len(t))
    return t, res, wall, elapsed


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else str(HERE / "{v}")
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4), constrained_layout=True)
    for zorder, (name, lw, marker, colour) in enumerate(VARIANTS):
        folder = Path(pattern.format(v=name))
        if not (folder / "residu.dat").exists():
            print(f"skipped {name}: no residu.dat in {folder}")
            continue
        t, res, wall, elapsed = load(folder)
        for ax, x in zip(axes, (t, wall)):
            ax.plot(
                x,
                res,
                lw=lw,
                marker=marker,
                ms=4.0,
                markevery=spread(len(x), 12),
                color=colour,
                zorder=zorder + 1,
                label=f"{name} ({elapsed:.0f} s)",
            )
    for ax in axes:
        for value, label in GATES:
            ax.axhline(
                value, color="k", lw=0.8, ls=":" if value < 1e-8 else "--", label=label
            )
        ax.set_yscale("log")
        ax.set_ylabel(r"rms change of $u$ per step ($U$)")
    axes[0].set_xlabel(r"$t$ ($D/U$)")
    axes[1].set_xlabel("wall time (s)")
    axes[0].set_title("(a) Residual against time", loc="left", fontsize=9)
    axes[1].set_title("(b) Residual against wall time", loc="left", fontsize=9)
    axes[0].legend(frameon=True, fontsize=6.5, loc="upper right")
    out = HERE / "plot_compare.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
