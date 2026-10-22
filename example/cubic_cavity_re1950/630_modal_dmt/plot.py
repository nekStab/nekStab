#!/usr/bin/env python3
"""Plot the DMT filter bank on the limit cycle of the cubic cavity at Re = 1950.

The run starts from the on-orbit checkpoint of 000_dns_period (t = 1700, time set to 0) and
integrates 15 periods T = 10.78909 of the cycle with four filters: k = 1, 2, 3 at the frequency of
the cycle and its second and third harmonic, and k = 4 at frequency 0 (the mean flow). The files
bRecav0.f0000k and bImcav0.f0000k hold Re and Im of the amplitude uhat_k of each filter, averaged over
the last 4 periods, with u(x, t) = Re[uhat_k(x) exp(i k omega t)] and t the time of the run.

The last field of the run, cav0.f00001, is the flow at t = 15 T. The script rebuilds it from the
filters, one harmonic at a time, and measures what is left:

  e_m = || u_end - (mean + sum_{k <= m} Re[uhat_k exp(i k omega t_end)]) || / || u_end - mean ||

for m = 0, 1, 2, 3, with the velocity in the kinetic-energy norm over the whole cavity (Gauss-Lobatto-
Legendre quadrature on the brick elements). The cycle is unstable to a slow mode of period 43.5, so
the field also holds a small part that no harmonic of the cycle carries; e_3 does not reach zero.

It writes bank_check.dat with one row per m = 0..3: m, e_m, and || uhat_m || / || mean || (for m = 0
the norm of the mean over itself, 1), and the file bank_settle.dat with one row per filter k = 1..3:
k, the amplitude of the filter at the start of the averaging window over its value at the end. It also
writes plot_dmt_bank.png.

Lengths are in units of the cavity edge L, speeds in units of the lid speed U, time in L/U.

Usage: python plot.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
from numpy.polynomial.legendre import Legendre
from pymech.neksuite import readnek

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
PERIOD = 10.78909  # 000_dns_period; the same value is in cav.usr
AVERAGE_FROM = 11.0 * PERIOD  # avg_from of cav.usr
HARMONICS = (1, 2, 3)  # filters 1 to 3; filter 4 is the mean flow
MEAN = 4


def gll_weights(n_points: int) -> np.ndarray:
    """Return the Gauss-Lobatto-Legendre quadrature weights on [-1, 1] for n_points nodes."""
    degree = n_points - 1
    legendre = Legendre.basis(degree)
    nodes = np.concatenate(([-1.0], np.sort(legendre.deriv().roots().real), [1.0]))
    return 2.0 / (degree * (degree + 1) * legendre(nodes) ** 2)


def read_field(path: Path) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Return the velocity, the node positions, the time and the quadrature weights of a field file.

    The velocity has the shape (elements, 3, nz, ny, nx), the positions the same shape. Each element
    is a brick, so its Jacobian is constant: one eighth of its volume. The weights have the shape
    (elements, nz, ny, nx).
    """
    data = readnek(str(path))
    velocity = np.stack([element.vel for element in data.elem])
    position = np.stack([element.pos for element in data.elem])
    weight = np.empty(velocity.shape[:1] + velocity.shape[2:])
    for index, element in enumerate(data.elem):
        n_z, n_y, n_x = element.pos.shape[1:]
        volume = np.prod([np.ptp(element.pos[axis]) for axis in range(3)])
        weight[index] = (
            gll_weights(n_z)[:, None, None]
            * gll_weights(n_y)[None, :, None]
            * gll_weights(n_x)[None, None, :]
            * volume
            / 8.0
        )
    return velocity, position, float(data.time), weight


def norm(field: np.ndarray, weight: np.ndarray) -> float:
    """Return the square root of the integral of |field|^2 over the cavity."""
    return float(np.sqrt(np.sum(weight * np.sum(field**2, axis=1))))


def mid_plane(
    position: np.ndarray, value: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y and value at the nodes of the plane z = 0, one entry per distinct point.

    The plane is the face of the elements whose lower edge is z = 0. A node on an element edge
    belongs to several elements; its entries are averaged.
    """
    on_plane = np.isclose(position[:, 2], 0.0, atol=1e-9)
    selected_x = position[:, 0][on_plane]
    selected_y = position[:, 1][on_plane]
    selected_v = value[on_plane]
    keys = np.round(np.stack([selected_x, selected_y], axis=1), 9)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    counts = np.bincount(inverse)
    mean_value = np.bincount(inverse, weights=selected_v) / counts
    return unique[:, 0], unique[:, 1], mean_value


def main() -> None:
    """Rebuild the final field from the filters, write the check files and draw the figure."""
    history = np.loadtxt(CASE_DIR / "dmt_bank.dat")
    time = history[:, 0]
    omega1 = 2.0 * np.pi / PERIOD

    final, position, t_end, weight = read_field(CASE_DIR / "cav0.f00001")
    amplitude: dict[int, np.ndarray] = {}
    for k in (*HARMONICS, MEAN):
        real, _, _, _ = read_field(CASE_DIR / f"bRecav0.f{k:05d}")
        imaginary, _, _, _ = read_field(CASE_DIR / f"bImcav0.f{k:05d}")
        amplitude[k] = np.stack([real, imaginary])  # (2, elements, 3, nz, ny, nx)

    mean_flow = amplitude[MEAN][0]
    fluctuation = final - mean_flow
    norm_mean = norm(mean_flow, weight)
    norm_fluctuation = norm(fluctuation, weight)

    rows = []
    partial = np.zeros_like(final)
    for m in (0, *HARMONICS):
        if m > 0:
            phase = m * omega1 * t_end
            partial = partial + (
                amplitude[m][0] * np.cos(phase) - amplitude[m][1] * np.sin(phase)
            )
        error = norm(fluctuation - partial, weight) / norm_fluctuation
        size = (
            1.0
            if m == 0
            else float(
                np.sqrt(
                    norm(amplitude[m][0], weight) ** 2
                    + norm(amplitude[m][1], weight) ** 2
                )
            )
            / norm_mean
        )
        rows.append([m, error, size])
    np.savetxt(
        CASE_DIR / "bank_check.dat",
        np.array(rows),
        header="m  e_m  norm_uhat_m/norm_mean    (e_m = |u_end - mean - sum_k<=m Re uhat_k exp(i k omega t_end)| / |u_end - mean|)",
        fmt=["%d", "%.8g", "%.8g"],
    )
    print(f"final field at t = {t_end:.5f} = {t_end / PERIOD:.4f} T")
    print("m  e_m  |uhat_m|/|mean|")
    for row in rows:
        print(f"{int(row[0])}  {row[1]:.4e}  {row[2]:.4e}")

    settle = []
    for k in HARMONICS:
        column = history[:, 2 + k]  # columns: time, |u|, |u - sum b|, then |uhat_k|
        start = column[np.searchsorted(time, AVERAGE_FROM)]
        settle.append([k, start / column[-1]])
    np.savetxt(
        CASE_DIR / "bank_settle.dat",
        np.array(settle),
        header="k  |uhat_k|(start of the averaging window) / |uhat_k|(end)",
        fmt=["%d", "%.8g"],
    )
    print("settling ratio at the start of the averaging window:", settle)

    figure = plt.figure(figsize=(7.2, 6.2), layout="constrained")
    grid = figure.add_gridspec(2, 6, height_ratios=[1.0, 0.9])
    for column_index, k in enumerate(HARMONICS):
        axes = figure.add_subplot(grid[0, 2 * column_index : 2 * column_index + 2])
        magnitude = np.sqrt(
            np.sum(amplitude[k] ** 2, axis=(0, 2))
        )  # (elements, nz, ny, nx)
        plane_x, plane_y, plane_value = mid_plane(position, magnitude)
        top = float(np.percentile(plane_value, 99.5))
        filled = axes.tricontourf(
            plane_x,
            plane_y,
            plane_value,
            levels=np.linspace(0.0, top, 25),
            cmap="viridis",
            extend="max",
        )
        axes.set_aspect("equal")
        axes.set_xlabel(r"$x$ ($L$)")
        if column_index == 0:
            axes.set_ylabel(r"$y$ ($L$)")
        label = ("a", "b", "c")[column_index]
        axes.set_title(
            rf"({label}) $|\hat u_{k}|$, $k\omega$ with $k={k}$, $z=0$",
            loc="left",
            fontsize=7.5,
        )
        figure.colorbar(
            filled,
            ax=axes,
            orientation="horizontal",
            shrink=0.9,
            pad=0.02,
            ticks=np.linspace(0.0, top, 3),
            format="%.1e",
        )

    axes_history = figure.add_subplot(grid[1, 0:4])
    for k, width in zip(HARMONICS, (3.0, 2.0, 1.2), strict=True):
        column = history[:, 2 + k]
        axes_history.plot(
            time / PERIOD,
            column / column[-1],
            lw=width,
            marker="o",
            ms=3,
            markevery=max(1, len(time) // 25),
            label=rf"$k={k}$, $\omega$ = {k * omega1:.3f}",
        )
    axes_history.axvline(
        AVERAGE_FROM / PERIOD, c="k", lw=0.7, ls="--", label="average from $11\\,T$"
    )
    axes_history.set_xlabel(r"$t / T$")
    axes_history.set_ylabel(r"$\|\hat u_k\|$ / value at $t=15\,T$")
    axes_history.set_title(
        "(d) The amplitude of each filter settles", loc="left", fontsize=7.5
    )
    axes_history.legend(loc="lower right", fontsize=6.5)

    axes_error = figure.add_subplot(grid[1, 4:6])
    axes_error.semilogy(
        [row[0] for row in rows],
        [row[1] for row in rows],
        c="C3",
        lw=1.6,
        marker="s",
        ms=4,
        label="final field",
    )
    axes_error.set_xticks(range(0, 4))
    axes_error.set_xlabel("harmonics in the sum")
    axes_error.set_ylabel(r"$e_m$")
    axes_error.set_title("(e) Flow from the filters", loc="left", fontsize=7.5)
    axes_error.legend(loc="upper right", fontsize=6.5)

    output = CASE_DIR / "plot_dmt_bank.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
