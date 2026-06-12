#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASE = Path(__file__).resolve().parent
TAUS = [1.0, 1.67, 2.78, 4.64, 7.74, 12.9]
TAU_DIR = {
    1.0: "t_1.0",
    1.67: "t_1.67",
    2.78: "t_2.78",
    4.64: "t_4.64",
    7.74: "t_7.74",
    12.9: "t_12.9",
}

# Only the old tau=2.78 scalar was preserved in docs/goals.md.  The old
# Spectre_Hp.dat files for the other completed runs were replaced by reruns.
OLD_GAIN = {
    1.0: "",
    1.67: "",
    2.78: "4796.6",
    4.64: "",
    7.74: "",
    12.9: "NaN",
}
OLD_STATUS = {
    1.0: "clean; old scalar overwritten",
    1.67: "clean; old scalar overwritten",
    2.78: "clean; docs/goals.md preserved scalar",
    4.64: "clean; old scalar overwritten",
    7.74: "clean; old scalar overwritten",
    12.9: "NaN; old log exited in k_normalize",
}


def first_gain(path: Path) -> float:
    with path.open() as handle:
        first = handle.readline().split()[0]
    return float(first.replace("D", "E"))


def barkley_interp(tau: float, data: np.ndarray) -> float:
    return float(np.interp(tau, data[:, 0], data[:, 1]))


def main() -> None:
    barkley = np.loadtxt(CASE / "barkley2008_fig5.ref")
    rows = []
    for tau in TAUS:
        run_dir = CASE / TAU_DIR[tau]
        new_gain = first_gain(run_dir / "Spectre_Hp.dat")
        new_log = run_dir / f"logfile.t{tau}.v2"
        old_log = run_dir / f"logfile.t{tau}"
        rows.append(
            {
                "tau": f"{tau:g}",
                "G_oldBF": OLD_GAIN[tau],
                "G_newBF": f"{new_gain:.7g}",
                "Barkley_interp": f"{barkley_interp(tau, barkley):.7g}",
                "status": f"new clean; {OLD_STATUS[tau]}",
                "old_log": old_log.relative_to(CASE).as_posix(),
                "new_log": new_log.relative_to(CASE).as_posix(),
            }
        )

    csv_path = CASE / "tg_six_point_table.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    md_path = CASE / "tg_six_point_table.md"
    with md_path.open("w") as handle:
        handle.write("| tau | G_oldBF | G_newBF | Barkley | clean/NaN | logs |\n")
        handle.write("| ---: | ---: | ---: | ---: | --- | --- |\n")
        for row in rows:
            handle.write(
                f"| {row['tau']} | {row['G_oldBF'] or 'not recovered'} | "
                f"{row['G_newBF']} | {row['Barkley_interp']} | {row['status']} | "
                f"old: {row['old_log']}; new: {row['new_log']} |\n"
            )

    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    ax.plot(barkley[:, 0], barkley[:, 1], color="0.15", lw=1.0, label="Barkley 2008 Fig. 5")
    ax.scatter(TAUS, [float(r["G_newBF"]) for r in rows], s=32, marker="o",
               facecolors="none", edgecolors="#0066cc", label="new converged BF")
    ax.scatter([2.78], [13.12], s=42, marker="x", color="#cc5500",
               label="main-BF tau=2.78, G=13.12")
    ax.scatter([2.78], [4796.6], s=30, marker="s", facecolors="none",
               edgecolors="#aa0000", label="old stale BF tau=2.78")
    ax.set_yscale("log")
    ax.set_xlim(0.8, 14.5)
    ax.set_ylim(1.0, 1.0e4)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$G(\tau)$")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.65)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig_path = CASE / "tg_six_point_logscale.png"
    fig.savefig(fig_path, dpi=300)
    print(csv_path)
    print(md_path)
    print(fig_path)


if __name__ == "__main__":
    main()
