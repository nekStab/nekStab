#!/usr/bin/env python3
"""Plot the Floquet analysis of the forced periodic orbit of the oscillating cylinder at Re = 100.

The orbit is the output of 210_baseflow_newton, the monodromy map is one period T = 1 / userParam06 of the forcing.
The script reads Spectre_Hd.dat (Floquet multipliers: real part, imaginary part, residual, converged flag),
Spectre_NSd.dat (the same values as sigma = ln|mu| / T and omega = arg(mu) / T), the leading mode
dRe1cyl0.f00001 and dIm1cyl0.f00001, and, when it exists, ../010_dns_forced/lockin_history.dat.

The forced orbit has no phase mode: every multiplier lies inside the unit circle when the wake locks to the forcing.
The forced DNS of 010_dns_forced approaches the orbit with the leading multipliers. Its distance to the orbit falls
by |mu_1| per period once the other multipliers have died out, so the change from period to period of the near
field falls by the same factor. The script compares |mu_1| with that factor.

It writes floquet.dat and plot_floquet.png. floquet.dat holds one row with:
  modulus      |mu_1|, the largest converged multiplier
  argument     arg(mu_1) in rad
  sigma        ln|mu_1| / T
  n_converged  number of converged multipliers
  dns_ratio    (c_N / c_(N-5))^(1/5), the mean factor per period of the near-field change c of the DNS over its
               last 5 periods (NaN without lockin_history.dat)
  orbit_period time stamp of BF_1cyl0.f00001, the length of the map (NaN without that file)

Usage: python plot.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
DNS_HISTORY = CASE_DIR.parent / "010_dns_forced" / "lockin_history.dat"
RATIO_PERIODS = 5
FIRST_FIT_PERIOD = 12


def load_nekplot():
    """Import the plot helpers shared by the examples."""
    sys.path.insert(0, str(CASE_DIR.parents[1]))
    import nekplot

    return nekplot


def par_value(name: str) -> float:
    """Return the number after `name =` in the .par file."""
    text = next(CASE_DIR.glob("*.par")).read_text()
    match = re.search(rf"^\s*{name}\s*=\s*([0-9.eE+-]+)", text, re.M | re.I)
    if match is None:
        raise ValueError(f"no {name} in the par file")
    return float(match.group(1))


def dns_ratio(change: np.ndarray) -> float:
    """Return the mean factor per period of the last RATIO_PERIODS periods of the change."""
    return float((change[-1] / change[-1 - RATIO_PERIODS]) ** (1.0 / RATIO_PERIODS))


def mode_amplitude(nk, real_file: Path, imag_file: Path):
    """Return x, y and sqrt(vx_re^2 + vy_re^2 + vx_im^2 + vy_im^2), scaled to a maximum of 1.

    A real multiplier has a real mode. The file of its imaginary part holds zeros, and the reader
    returns no field for it.
    """
    x, y, real_part, _ = nk.read_field(str(real_file))
    _, _, imag_part, _ = nk.read_field(str(imag_file))
    zero = np.zeros_like(real_part["vx"])
    amplitude = np.sqrt(
        real_part["vx"] ** 2
        + real_part["vy"] ** 2
        + imag_part.get("vx", zero) ** 2
        + imag_part.get("vy", zero) ** 2
    )
    return x, y, amplitude / amplitude.max()


def main() -> None:
    """Write floquet.dat and draw the multipliers, the DNS relaxation and the leading mode."""
    nk = load_nekplot()
    nk.configure_style()

    period = 1.0 / par_value("userParam06")
    ritz = np.atleast_2d(np.loadtxt(CASE_DIR / "Spectre_Hd.dat"))
    converged = ritz[:, 3] > 0.5
    multiplier = ritz[:, 0] + 1j * ritz[:, 1]
    leading = int(np.argmax(np.where(converged, np.abs(multiplier), -1.0)))
    modulus = float(np.abs(multiplier[leading]))
    argument = float(np.angle(multiplier[leading]))
    sigma = float(np.log(modulus) / period)

    have_dns = DNS_HISTORY.is_file()
    if have_dns:
        history = np.loadtxt(DNS_HISTORY)
        periods, change = history[:, 0], history[:, 1]
        ratio = dns_ratio(change)
    else:
        ratio = float("nan")
    orbit_file = CASE_DIR / "BF_1cyl0.f00001"
    orbit_period = (
        float(nk.read_field(orbit_file)[3]) if orbit_file.is_file() else float("nan")
    )
    np.savetxt(
        CASE_DIR / "floquet.dat",
        [[modulus, argument, sigma, int(converged.sum()), ratio, orbit_period]],
        header="modulus  argument  sigma  n_converged  dns_ratio  orbit_period",
        fmt=["%.8g", "%.8g", "%.8g", "%d", "%.8g", "%.12g"],
    )
    print(
        f"leading multiplier |mu| = {modulus:.6f}, arg = {argument:.4f} rad, sigma = {sigma:.5f}; "
        f"{int(converged.sum())} converged; DNS factor per period {ratio:.4f}"
    )

    fig = plt.figure(figsize=(7.2, 5.0), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 0.8])
    ax_circle = fig.add_subplot(grid[0, 0])
    ax_dns = fig.add_subplot(grid[0, 1])
    ax_mode = fig.add_subplot(grid[1, :])

    angle = np.linspace(0.0, 2.0 * np.pi, 361)
    ax_circle.plot(np.cos(angle), np.sin(angle), c="0.45", lw=0.8, label="unit circle")
    ax_circle.plot(
        multiplier.real[~converged],
        multiplier.imag[~converged],
        ls="none",
        marker="o",
        ms=3.0,
        mfc="none",
        c="0.5",
        label="not converged",
    )
    ax_circle.plot(
        multiplier.real[converged],
        multiplier.imag[converged],
        ls="none",
        marker="o",
        ms=4.5,
        c="C0",
        label=f"converged ({int(converged.sum())})",
    )
    ax_circle.plot(
        [multiplier[leading].real],
        [multiplier[leading].imag],
        ls="none",
        marker="s",
        ms=6.0,
        mfc="none",
        c="C3",
        label=rf"$|\mu_1|$ = {modulus:.4f}",
    )
    ax_circle.set_aspect("equal")
    ax_circle.set_xlim(-1.15, 1.15)
    ax_circle.set_ylim(-1.15, 1.15)
    ax_circle.set_xlabel(r"$\mathrm{Re}\,\mu$")
    ax_circle.set_ylabel(r"$\mathrm{Im}\,\mu$")
    ax_circle.set_title("(a) Floquet multipliers", loc="left", fontsize=8)
    ax_circle.legend(loc="lower left", fontsize=6.0)

    if have_dns:
        ax_dns.semilogy(
            periods,
            change,
            c="C0",
            lw=1.2,
            marker="o",
            ms=3.0,
            label="forced DNS, near field",
        )
        fit = periods >= FIRST_FIT_PERIOD
        ax_dns.semilogy(
            periods[fit],
            change[-1] * modulus ** (periods[fit] - periods[-1]),
            c="C3",
            lw=1.0,
            ls="--",
            marker="s",
            ms=2.5,
            label=rf"$|\mu_1|^n$ = {modulus:.4f}$^n$, from the last period",
        )
        ax_dns.axvline(
            FIRST_FIT_PERIOD, c="k", lw=0.5, ls=":", label=f"period {FIRST_FIT_PERIOD}"
        )
        ax_dns.legend(loc="upper right", fontsize=6.0)
    ax_dns.set_xlabel("forcing period $n$")
    ax_dns.set_ylabel(r"$\|u_n - u_{n-1}\| / \|u_n\|$")
    ax_dns.set_title("(b) Relaxation of the DNS", loc="left", fontsize=8)

    mode_files = (CASE_DIR / "dRe1cyl0.f00001", CASE_DIR / "dIm1cyl0.f00001")
    if all(path.is_file() for path in mode_files):
        x, y, amplitude = mode_amplitude(nk, *mode_files)
        filled = ax_mode.tricontourf(
            nk.make_triangulation(x, y),
            amplitude,
            levels=np.linspace(0.0, 1.0, 21),
            cmap="Blues",
        )
        nk.add_cylinder_patches(ax_mode)
        ax_mode.set_aspect("equal")
        ax_mode.set_xlim(-2.0, 60.0)
        ax_mode.set_ylim(-3.0, 3.0)
        fig.colorbar(filled, ax=ax_mode, label="amplitude (maximum = 1)", shrink=0.8)
    ax_mode.set_xlabel(r"$x$ ($D$)")
    ax_mode.set_ylabel(r"$y$ ($D$)")
    ax_mode.set_title(
        "(c) Leading mode at phase 0 of the forcing", loc="left", fontsize=8
    )

    output = CASE_DIR / "plot_floquet.png"
    fig.savefig(output, dpi=400, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
