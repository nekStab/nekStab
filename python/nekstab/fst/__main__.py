"""CLI entry point for the nekStab FST inflow generator.

Examples
--------
Generate the slot_FST Re=495 inflow (matches example/slot_FST defaults):
    PYTHONPATH=python python -m nekstab.fst \
        --re 495 --ny 210 --ly 20 --numk 25 \
        --kmin 1.122 --kmax 22.44 --seed 0 \
        --out example/slot_fst_re495/FST_inflow.dat

Notes
-----
For a Blasius inflow the wavenumber band is set by the domain geometry:
    kmin = 2*pi / (Ly * delta_star)        (largest resolvable wavelength)
    kmax = 2*pi / (dy_min * delta_star)     (smallest resolvable wavelength)
Pass --kmin/--kmax directly (case-agnostic), or --delta-star/--dy-min to derive
them with the slot_FST convention.
"""

import argparse
import math

from .generator import generate_fst_inflow, write_inflow_file


def _build_parser():
    p = argparse.ArgumentParser(
        prog="python -m nekstab.fst",
        description="Generate Free-Stream-Turbulence inflow (continuous-spectrum "
                    "Orr-Sommerfeld-Squire modes) as one consolidated file.")
    p.add_argument("--re", type=float, required=True,
                   help="Reynolds number based on displacement thickness delta*.")
    p.add_argument("--ny", type=int, required=True,
                   help="Wall-normal collocation points (e.g. Ney*lx1).")
    p.add_argument("--ly", type=float, required=True,
                   help="Domain height in delta* units.")
    p.add_argument("--numk", type=int, required=True,
                   help="Number of spherical wavenumber shells (n_modes = numk*10).")
    # Wavenumber band: either give kmin/kmax directly, or derive from geometry.
    p.add_argument("--kmin", type=float, default=None, help="Smallest wavenumber.")
    p.add_argument("--kmax", type=float, default=None, help="Largest wavenumber.")
    p.add_argument("--delta-star", type=float, default=None,
                   help="delta* (used with --dy-min to derive kmin/kmax).")
    p.add_argument("--dy-min", type=float, default=None,
                   help="Minimum near-wall spacing (used to derive kmax).")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed (shell rotations + phases) for reproducibility.")
    p.add_argument("--tu", type=float, default=0.010,
                   help="Turbulence intensity Tu.")
    p.add_argument("--length", type=float, default=None,
                   help="von Karman length scale; default 1.80/kmax.")
    p.add_argument("--out", required=True, help="Output consolidated inflow file.")
    p.add_argument("--format", choices=("dat", "bin"), default="dat",
                   help="Output format: 'dat' (formatted ASCII, default) or "
                        "'bin' (little-endian stream binary, faster Fortran read).")
    p.add_argument("--jobs", type=int, default=1,
                   help="Parallel worker processes for the independent OSS solves "
                        "(default 1 = serial).")
    p.add_argument("--quiet", action="store_true", help="Suppress per-mode progress.")
    return p


def _resolve_band(args):
    """Return (kmin, kmax), either explicit or derived from geometry."""
    if args.kmin is not None and args.kmax is not None:
        return args.kmin, args.kmax
    if args.delta_star is not None and args.dy_min is not None:
        kmin = 2.0 * math.pi / (args.ly * args.delta_star)
        kmax = 2.0 * math.pi / (args.dy_min * args.delta_star)
        return kmin, kmax
    raise SystemExit("error: provide either --kmin and --kmax, "
                     "or --delta-star and --dy-min to derive them.")


def main(argv=None):
    args = _build_parser().parse_args(argv)
    kmin, kmax = _resolve_band(args)
    print(f"FST inflow: Re={args.re} Ny={args.ny} Ly={args.ly} "
          f"numk={args.numk} (n_modes={args.numk * 10}) "
          f"kmin={kmin:.4f} kmax={kmax:.4f}")
    result = generate_fst_inflow(args.re, args.ny, args.ly, args.numk,
                                 kmin, kmax, seed=args.seed, jobs=args.jobs,
                                 verbose=not args.quiet, tu=args.tu, length=args.length)
    write_inflow_file(args.out, result, fmt=args.format)
    print(f"Wrote {args.out}  ({len(result['modes'])} modes, Ny={args.ny}, "
          f"format={args.format})")


if __name__ == "__main__":
    main()
