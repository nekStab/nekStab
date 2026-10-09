# Triple-port jet — forced periodic orbit (Newton-GMRES), Re=1900

> **Operating point: Re = 1900** (`viscosity = -1900`). See the family-level
> note in `../README.md` about the `_re2005` directory naming.

## Physics
Axisymmetric triple-port jet at Re = 1900 (above the critical Re_c ≈ 1371).
Computes the forced periodic orbit (harmonic forcing at St = 0.6, period
T = 1/St = 1.66667) with Newton-GMRES.

## nekStab Mode
- `userParam01 = 2.2` — Newton-GMRES for forced periodic orbits
- `userParam07 = 100` — Krylov dimension `k_dim`
- `endTime = 1.666667` — forcing period (St = 0.6)
- `viscosity = -1900.0` — Re = 1900

## Initial condition
`startFrom = BF_Re1900_tpjet0.f00001` is a field of the Re = 1900 jet at
polynomial order 5 (`lx1 = 6`). `SIZE` has `lx1 = 5`, so Nek5000 interpolates
the field to the order of the run. The field is not on the periodic orbit of
this mesh and order. The first Newton residual is 1.13e-4, which is 11000 times
the tolerance of 1e-8.

Newton writes the orbit to `BF_tpjet0.f00001`. The stage does not ship this
file, because a shipped orbit would let the run stop after one iteration.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 16 ranks
```

## Convergence record (run of 2026-10-10, Jean Zay, ifort 2021.9, 16 ranks)
The run takes 2370 s. The residual is the squared norm of F(q) − q
(`residu_newton.dat`, column 7). Each linear solve stops at one fifth of the
residual at the start of its iteration (`ifdyntol` in `tpjet.usr`), so the
residual falls by a factor of 4 to 6 per iteration:

| Newton iter | residual  |
|-------------|-----------|
| 1           | 1.130e-4  |
| 2           | 2.615e-5  |
| 3           | 6.417e-6  |
| 4           | 1.228e-6  |
| 5           | 2.277e-7  |
| 6           | 3.921e-8  |
| 7           | 9.253e-9  → **converged** (target 1e-8) |

The six linear solves use 7, 17, 12, 7, 9 and 17 Arnoldi vectors (69 in total).

## Reference (ref/)
- `reference.json` checks six quantities. Three come from `residu_newton.dat`:
  the first residual (1.13e-4, tolerance 2 %), the iteration count (7,
  tolerance 1) and the final residual (below the 1e-8 target). Three come from
  `BF_tpjet0.f00001`: the norms of `vx` and `vy` and the largest `vx`
  (tolerance 1e-3). The first residual and the iteration count fail if the run
  starts from a field that is already on the orbit.
  Verify with
  `scripts/check_against_ref.py example/tpjet_re2005/210_baseflow_newton`.
- `ref/bf.png`, `ref/ic.png`, `ref/residual.png`.

## Related
`../311_stability_direct_floquet` ships its own copy of the orbit. That copy
differs from the orbit of this stage by 4.7e-4 (relative l2 norm of the
velocity). Both reach the tolerance of 1e-8 on the squared residual.
