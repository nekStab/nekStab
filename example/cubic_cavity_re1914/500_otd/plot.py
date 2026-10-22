#!/usr/bin/env python3
"""Plot the OTD basis of the cubic cavity at Re = 1914 and check that it stays on the leading eigenplane.

The base flow is the steady flow of 210_baseflow_newton, frozen. The two OTD modes start on the plane
spanned by the real and imaginary parts of the leading eigenmode of 310_stability_direct. That plane is
invariant under the linearized operator, so the basis must stay on it, and the mean of the two diagonal
entries of the reduced operator Lr (half the trace of Lr on the plane) must equal the growth rate sigma
of the eigenpair. The entries alone do not: the basis rotates in the plane at the frequency omega of the
pair, and each entry follows the rotation.

The script reads otd_eigenvalues.dat (time and the two entries), otd_ftle.dat (running exponents) and the
basis written during the run (r01cav0.f*, r02cav0.f*). The plane of the eigenpair is read from
../310_stability_direct (dRecav0.f00001, dImcav0.f00001), its growth rate from Spectre_NSd.dat there.
The inner product is the kinetic energy, computed by GLL quadrature on the brick elements of the mesh.
For each written basis it computes the sine of the largest principal angle between the basis and the plane.

It writes otd_summary.dat and plot_otd.png. Lengths are in units of the cavity edge L, speeds in units of
the lid speed U, times in units of L/U and rates in units of U/L.

Usage: python plot.py [directory of the 310 stage]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import Legendre
from pymech.neksuite import readnek

CASE_DIR = Path(__file__).resolve().parent
EIGEN_DIR = (
    Path(sys.argv[1]) if len(sys.argv) > 1 else CASE_DIR.parent / "310_stability_direct"
)
PLANE_TOLERANCE = 1e-5  # tol_abs of plane_sine_final in ref/reference.json


def gll_weights(n_points: int) -> np.ndarray:
    """Return the Gauss-Lobatto-Legendre quadrature weights on [-1, 1] for n_points nodes."""
    degree = n_points - 1
    legendre = Legendre.basis(degree)
    nodes = np.concatenate(([-1.0], np.sort(legendre.deriv().roots().real), [1.0]))
    return 2.0 / (degree * (degree + 1) * legendre(nodes) ** 2)


def read_velocity(field_file: Path) -> tuple[np.ndarray, float, np.ndarray]:
    """Return the velocity of a field file, its time and the quadrature weight of every node.

    The velocity has the shape (elements, 3, nz, ny, nx). Each element is a brick, so its Jacobian is
    constant: one eighth of its volume. The weights have the shape (elements, nz, ny, nx).
    """
    data = readnek(str(field_file))
    velocity = np.stack([element.vel for element in data.elem])
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
    return velocity, float(data.time), weight


def inner(a: np.ndarray, b: np.ndarray, weight: np.ndarray) -> float:
    """Return the integral of a . b over the cavity."""
    return float(np.sum(weight * np.sum(a * b, axis=1)))


def orthonormalize(fields: list[np.ndarray], weight: np.ndarray) -> list[np.ndarray]:
    """Return the modified Gram-Schmidt orthonormalisation of the fields."""
    basis: list[np.ndarray] = []
    for field in fields:
        field = field.copy()
        for other in basis:
            field -= inner(field, other, weight) * other
        basis.append(field / np.sqrt(inner(field, field, weight)))
    return basis


def plane_sine(
    basis: list[np.ndarray], plane: list[np.ndarray], weight: np.ndarray
) -> float:
    """Return the sine of the largest principal angle between span(basis) and span(plane).

    Both lists hold orthonormal fields. The sine is the largest singular value of the part of the
    basis that lies outside the plane.
    """
    outside = []
    for field in basis:
        part = field.copy()
        for other in plane:
            part -= inner(field, other, weight) * other
        outside.append(part)
    gram = np.array([[inner(a, b, weight) for b in outside] for a in outside])
    return float(np.sqrt(np.linalg.eigvalsh(gram).max()))


def main() -> None:
    """Write otd_summary.dat and draw the entries of Lr, their mean and the distance to the eigenplane."""
    sigma = float(
        (EIGEN_DIR / "Spectre_NSd.dat").read_text().split()[0].replace("D", "E")
    )
    entries = np.loadtxt(CASE_DIR / "otd_eigenvalues.dat")
    ftle = np.loadtxt(CASE_DIR / "otd_ftle.dat")
    entries_mean = 0.5 * (entries[:, 1] + entries[:, 2])
    ftle_mean = 0.5 * (ftle[:, 1] + ftle[:, 2])

    real_part, _, weight = read_velocity(EIGEN_DIR / "dRecav0.f00001")
    imaginary_part, _, _ = read_velocity(EIGEN_DIR / "dImcav0.f00001")
    plane = orthonormalize([real_part, imaginary_part], weight)
    times, sines = [], []
    for first in sorted(CASE_DIR.glob("r01cav0.f*")):
        second = first.with_name(first.name.replace("r01", "r02"))
        mode_1, time, _ = read_velocity(first)
        mode_2, _, _ = read_velocity(second)
        times.append(time)
        sines.append(
            plane_sine(orthonormalize([mode_1, mode_2], weight), plane, weight)
        )
    times, sines = np.array(times), np.array(sines)

    after_start = entries[:, 0] > 0.0
    spread = float(np.ptp(entries_mean[after_start]))
    np.savetxt(
        CASE_DIR / "otd_summary.dat",
        [
            [
                entries[-1, 0],
                sigma,
                entries_mean[0],
                entries_mean[-1],
                spread,
                ftle_mean[-1],
                sines[-1],
                sines.max(),
            ]
        ],
        header=(
            "t_end  sigma_310  entries_mean_start  entries_mean_final  entries_mean_spread  "
            "ftle_mean_final  plane_sine_final  plane_sine_max"
        ),
        fmt="%.8g",
    )
    print(
        f"sigma of 310 {sigma:.6e}; mean of the entries of Lr at t = 0 {entries_mean[0]:.6e}, "
        f"at t = {entries[-1, 0]:.1f} {entries_mean[-1]:.6e} (spread after t = 0: {spread:.2e}); "
        f"mean exponent {ftle_mean[-1]:.6e}; sine of the largest angle to the eigenplane "
        f"{sines[-1]:.2e} (largest over {sines.size} bases {sines.max():.2e})"
    )

    fig, (ax_entries, ax_mean, ax_plane) = plt.subplots(
        1, 3, figsize=(8.0, 2.8), layout="constrained"
    )
    marker_gap = max(1, len(entries) // 12)
    for column, (lw, marker, label) in enumerate(
        [(2.4, "o", "entry 1"), (1.6, "s", "entry 2")], start=1
    ):
        ax_entries.plot(
            entries[:, 0],
            entries[:, column],
            lw=lw,
            marker=marker,
            ms=3,
            markevery=marker_gap,
            label=label,
        )
    low, high = entries[:, 1:].min(), entries[:, 1:].max()
    ax_entries.set_ylim(low - 0.45 * (high - low), high + 0.05 * (high - low))
    ax_entries.set_xlabel(r"$t$ ($L/U$)")
    ax_entries.set_ylabel(r"diagonal of $L_r$ ($U/L$)")
    ax_entries.set_title("(a) Entries of $L_r$", loc="left", fontsize=8)
    ax_entries.legend(loc="lower center", ncol=2, fontsize=6.5)

    ax_mean.plot(
        entries[:, 0],
        entries_mean,
        c="C0",
        lw=2.4,
        marker="o",
        ms=3,
        markevery=marker_gap,
        label=f"mean of the entries ({entries_mean[-1]:.2e})",
    )
    ax_mean.plot(
        ftle[:, 0],
        ftle_mean,
        c="C1",
        lw=1.4,
        marker="s",
        ms=3,
        markevery=max(1, len(ftle) // 12),
        label=f"mean running exponent ({ftle_mean[-1]:.2e})",
    )
    ax_mean.axhline(
        sigma, c="k", lw=0.8, ls="--", label=rf"$\sigma$ of 310 ({sigma:.2e})"
    )
    low = min(entries_mean.min(), ftle_mean.min(), sigma)
    high = max(entries_mean.max(), ftle_mean.max(), sigma)
    ax_mean.set_ylim(low - 2.4 * (high - low), high + 0.3 * (high - low))
    ax_mean.set_xlabel(r"$t$ ($L/U$)")
    ax_mean.set_ylabel(r"mean over the two modes ($U/L$)")
    ax_mean.set_title("(b) Mean of the two modes", loc="left", fontsize=8)
    ax_mean.legend(loc="lower center", fontsize=6.5)

    ax_plane.semilogy(
        times,
        np.maximum(sines, 1e-9),
        c="C2",
        lw=1.6,
        marker="o",
        ms=4,
        label="OTD basis",
    )
    ax_plane.axhline(
        PLANE_TOLERANCE, c="k", lw=0.8, ls="--", label="tolerance of the reference"
    )
    ax_plane.set_ylim(1e-8, 1.0)
    ax_plane.set_xlabel(r"$t$ ($L/U$)")
    ax_plane.set_ylabel(r"$\sin\theta_{\max}$ to the plane")
    ax_plane.set_title("(c) Distance to the plane", loc="left", fontsize=8)
    ax_plane.legend(loc="lower right", fontsize=6.5)

    output = CASE_DIR / "plot_otd.png"
    fig.savefig(output, dpi=400)
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
