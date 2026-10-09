# Flip-Flop UPO at Re=62 — Newton-GMRES

## Physics
2D flow around two side-by-side circular cylinders at Re = 62 (slightly above
Re_c ≈ 61.17). Computes the unstable periodic orbit (UPO) of the flip-flop
instability using Newton-GMRES. The converged orbit has period T = 8.815177.

## nekStab Mode
- `userParam01 = 2.1` — Newton-GMRES for UPOs
- `userParam07 = 100` — Krylov subspace dimension (GMRES for Newton linear solve)
- `endTime = 8.73356` — starting value of the period. The Newton iteration changes the period; the converged value is the time stamp
  of `BF_2cyl0.f00001` (T = 8.815177). The stages 311, 321 and 411 read the period from that time stamp.
- `viscosity = -62.0` — Re = 62

## Initial condition
`startFrom = BF_Re60_2cyl0.f00001` — a nearby Re=60 orbit snapshot used as the
Newton initial guess (shipped here as the seed). Newton converges the orbit at
Re = 62 from this guess.

## Run
```bash
mks 2cyl
sbatch run.local.slurm     # 8 ranks
```

## Result (2026-10-09, Jean Zay, `ifort`, 16 ranks)
- Ten Newton iterations. The residual (column 7 of `residu_newton.dat`) falls from 5.77e-2 to 8.05e-9; the tolerance is 1e-8.
- The period changes from 8.73356 to 8.815177.
- Wall time = 922.6 s.

## Reference (ref/)
- `reference.json` — the final Newton residual (`residu_newton.dat`, column 7, last row, below 1e-8), the period of the
  orbit (time stamp of `BF_2cyl0.f00001`, within 1e-5) and five statistics of its fields (min, max and norm of `vx`, min and max of `vy`, within
  1e-4 to 5e-4). Verify a fresh run with `scripts/check_against_ref.py example/flip_flop_re62/210_baseflow_newton`.
- `ref/bf.png` — converged UPO base flow
- `ref/ic.png` — initial guess
- `ref/residual.png` — Newton residual decay

The converged orbit `BF_2cyl0.f00001` is the input to
`../311_stability_direct_floquet`.

## Reference
Carini et al. (2015), JFM 778.
