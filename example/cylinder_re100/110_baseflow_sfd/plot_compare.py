#!/usr/bin/env python3
"""plot_compare.py: residual of the base-flow methods of the cylinder wake against time and against wall time.

The curves are the six SFD variants of this stage, BoostConv (120_baseflow_boostconv, at Re = 50, not 100),
DMT (130_baseflow_dmt) and Newton (210_baseflow_newton/fp, one point per Newton iteration).
SFD, BoostConv and DMT write the root mean square of the change of the velocity in one time step;
Newton writes the norm of F(q) = U(T) - q of each iterate. The residuals are not the same quantity.

plot_compare.png
    (a) residual against the time of the simulation, with the stop gates as guide lines;
        for Newton the time is the total time integrated up to that iteration;
    (b) the same residual against the wall time of the run, read from the `Step` lines of each logfile.
        Newton prints no `Step` line per iteration: its wall time is the total elapsed time
        scaled by the fraction of the total time integrated.
The legend gives the wall time of each run.
plot_compare_akervik.png and plot_compare_casacuberta.png: the fixed, dyn and dyn OIFS runs of one formulation,
    against the time step (a) and the wall time (b); see formulation_figure().

USAGE: python plot_compare.py [PATTERN]
    PATTERN is the folder of a run, with {v} for its path from the family folder (for example 110_baseflow_sfd/akervik)
    or {d} for the same path with "-" in place of "/". The default is ../{v}. The folder holds the residual file and logfile.
"""

import re
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = [
    # path from the family folder, label, residual file, time and residual columns, line width, marker, colour;
    # the thickest line is drawn first, at the back
    ("210_baseflow_newton/fp", "Newton", "residu_newton.dat", 4, 6, 3.4, "*", "k"),
    ("130_baseflow_dmt", "DMT", "residu_dmt.dat", 0, 1, 3.0, "P", "C7"),
    ("120_baseflow_boostconv", "BoostConv, Re = 50", "residu.dat", 0, 1, 2.7, "h", "C6"),
    ("110_baseflow_sfd/casacuberta_dyn_oifs", "casacuberta_dyn_oifs", "residu.dat", 0, 1, 2.4, "v", "C5"),
    ("110_baseflow_sfd/akervik_dyn_oifs", "akervik_dyn_oifs", "residu.dat", 0, 1, 2.1, "^", "C4"),
    ("110_baseflow_sfd/casacuberta_dyn", "casacuberta_dyn", "residu.dat", 0, 1, 1.8, "D", "C3"),
    ("110_baseflow_sfd/casacuberta", "casacuberta", "residu.dat", 0, 1, 1.5, "s", "C2"),
    ("110_baseflow_sfd/akervik_dyn", "akervik_dyn", "residu.dat", 0, 1, 1.2, "o", "C1"),
    ("110_baseflow_sfd/akervik", "akervik", "residu.dat", 0, 1, 0.9, "x", "C0"),
]
STEP = re.compile(r"^Step\s+\d+, t=\s*\S+, DT=\s*\S+, C=\s*\S+\s+(\S+)\s+\S+")
ELAPSED = re.compile(r"total elapsed time\s*:\s*(\S+)")
GATES = [
    (1e-10, ":", "gate 1e-10 (SFD, BDF3)"),
    (1e-9, "-.", "gate 1e-9 (Newton)"),
    (1e-7, (0, (5, 1, 1, 1, 1, 1)), "gate 1e-7 (DMT)"),
    (1e-6, "--", "gate 1e-6 (SFD OIFS, BoostConv)"),
]


def spread(n_rows, n_marks):
    """Indices of n_marks samples spread evenly over n_rows rows."""
    return sorted({int(i) for i in np.linspace(0, n_rows - 1, n_marks)})


def load(folder, residual_file, tcol, rcol):
    """Time, residual, wall time of each row and the total elapsed time of one run."""
    data = np.loadtxt(folder / residual_file, ndmin=2)
    t, res = data[:, tcol], data[:, rcol]
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
    if len(steps) >= len(t) and residual_file != "residu_newton.dat":
        wall = np.array(steps[: len(t)])
    else:
        # No wall time per row: the cost of a run grows with the time integrated.
        wall = elapsed * t / t[-1]
    return t, res, wall, elapsed


def formulation_figure(name, pattern):
    """Residual of the fixed, dyn and dyn OIFS runs of one formulation, against time step and wall time.

    Writes HERE / f"plot_compare_{name}.png". A run without residu.dat or logfile is skipped with a message.
    The gate lines are the largest residualTol of each run's own 1cyl.par.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), layout="constrained")
    gates = set()
    # suffix, line width, marker, colour; the thickest line is drawn first, at the back
    runs = (("", 2.4, "x", "C0"), ("_dyn", 1.6, "o", "C1"), ("_dyn_oifs", 0.9, "^", "C4"))
    for zorder, (suffix, lw, marker, colour) in enumerate(runs):
        run = f"{name}{suffix}"
        folder = Path(pattern.format(v=f"110_baseflow_sfd/{run}", d=f"110_baseflow_sfd-{run}"))
        if not (folder / "residu.dat").exists() or not (folder / "logfile").exists():
            print(f"skipped {run}: no residu.dat or logfile in {folder}")
            continue
        t, res, wall, elapsed = load(folder, "residu.dat", 0, 1)
        steps = np.arange(1, len(res) + 1)
        for ax, x in zip(axes, (steps, wall)):
            ax.plot(
                x,
                res,
                lw=lw,
                marker=marker,
                ms=5.0,
                markevery=spread(len(x), 12),
                color=colour,
                zorder=zorder + 1,
                label=f"{run} ({len(res)} steps, {elapsed:.0f} s)",
            )
        par = (HERE / run / "1cyl.par").read_text(errors="replace")
        gates.add(max(float(v) for v in re.findall(r"^\s*residualTol\s*=\s*(\S+)", par, re.M)))
    for value in sorted(gates, reverse=True):
        for ax in axes:
            ax.axhline(value, color="k", lw=0.8, ls=":", label=f"gate {value:.0e}")
    for ax in axes:
        ax.set_yscale("log")
        ax.set_ylabel("residual (nondimensional)")
    axes[0].set_xlabel("time step")
    # Four ticks at most: five-digit step counts overlap otherwise.
    axes[0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    axes[1].set_xlabel("wall time (s)")
    axes[0].set_title("(a) Residual against time step", loc="left", fontsize=9)
    axes[1].set_title("(b) Residual against wall time", loc="left", fontsize=9)
    # One legend for both panels, outside the axes: inside it hides the curves.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=True, fontsize=6.5, loc="outside right upper")
    out = HERE / f"plot_compare_{name}.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else str(HERE.parent / "{v}")
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6), layout="constrained")
    for zorder, (path, name, residual_file, tcol, rcol, lw, marker, colour) in enumerate(RUNS):
        folder = Path(pattern.format(v=path, d=path.replace("/", "-")))
        if not (folder / residual_file).exists():
            print(f"skipped {name}: no {residual_file} in {folder}")
            continue
        t, res, wall, elapsed = load(folder, residual_file, tcol, rcol)
        for ax, x in zip(axes, (t, wall)):
            ax.plot(
                x,
                res,
                lw=lw,
                marker=marker,
                ms=5.0,
                markevery=spread(len(x), 12),
                color=colour,
                zorder=zorder + 1,
                label=f"{name} ({elapsed:.0f} s)",
            )
    for ax in axes:
        for value, ls, label in GATES:
            ax.axhline(value, color="k", lw=0.8, ls=ls, label=label)
        ax.set_yscale("log")
        ax.set_ylabel("residual (nondimensional)")
    axes[0].set_xlabel(r"$t$ ($D/U$)")
    axes[1].set_xlabel("wall time (s)")
    axes[0].set_title("(a) Residual against time", loc="left", fontsize=9)
    axes[1].set_title("(b) Residual against wall time", loc="left", fontsize=9)
    # One legend for both panels, outside the axes: inside it hides the curves.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=True, fontsize=6.5, loc="outside right upper")
    out = HERE / "plot_compare.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"Saved {out}")
    for name in ("akervik", "casacuberta"):
        formulation_figure(name, pattern)


if __name__ == "__main__":
    main()
