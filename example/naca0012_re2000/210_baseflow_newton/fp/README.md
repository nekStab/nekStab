# NACA 0012 base flow at Re=2000 — Newton-GMRES

## Physics
2D flow around a NACA 0012 airfoil at Re = 2000. Computes the steady (unstable)
base flow using Newton-GMRES with sponge regions.

## nekStab Mode
- `userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam07 = 220` — Krylov subspace dimension
- `userParam10 = 1.7` — sponge strength (left = 1, right = 5)
- `viscosity = -2000.0` — Re = 2000

## Initial condition
`startFrom = BF_seed_naca00120.f00001` is the last snapshot (t = 70) of
`../../000_dns`. The flow at that time is not steady: the first Newton residual
is 6.08e-7, 600 times the tolerance of 1e-9. `residu_newton.dat` holds the squared
norm of F(q) - q. Newton writes the base flow to `BF_naca00120.f00001` and leaves
the seed unchanged. `SPGnaca00120.f00001` is the sponge field.

## The base flow depends on the seed
The sponge forces the velocity toward the field that the run holds when it starts.
In this stage that field is the seed. The Newton solution therefore solves the
Navier-Stokes equations plus a forcing toward the seed, and the field inside the
sponge zone depends on the seed. The right sponge starts at x = 3.4, reaches full
strength at x = 5.1 and ends at the outflow, x = 6.75. Three Newton runs in a row, each seeded
with the output of the one before (8 ranks, Jean Zay), and the direct stability
stage run on each base flow:

| pass | seed | first residual | Newton iterations | time (s) | field change | sigma | omega |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | DNS snapshot, t = 70 | 6.08e-7 | 5 | 956 | | -0.06427 | 8.0244 |
| 2 | output of pass 1 | 3.76e-8 | 4 | 494 | 3.5e-3 | -0.06610 | 8.0169 |
| 3 | output of pass 2 | 9.69e-9 | 3 | 227 | 2.7e-3 | | |

The field change is the relative l2 norm of the velocity difference to the previous
pass. For x < 3 it falls by a factor of 4 to 5 per pass (6e-4 to 8e-4, then 1e-4 to
4e-4); for 5 < x < 6 it falls by a factor of 1.5 (1.2e-2, then 8.1e-3). The base
flow that `../../310_stability_direct` ships differs from the pass-1 field by 1.1e-2
(7e-4 to 8e-4 for x < 1, 3e-2 for 5 < x < 6) and gives sigma = -0.06654 and
omega = 8.0058. The stage runs one pass. The growth rate of pass 1 lies 0.0023 above the
shipped value (3.4 % of its magnitude), and that of pass 2 lies 0.0004 above (0.7 %).

## Mesh
4368 elements (`lelg=4368`), polynomial order N=5 (`lx1=6`).

## Run
```bash
mks naca0012
sbatch run.local.slurm     # 8 ranks
```
The run works at 8 to 40 ranks. `[VELOCITY] residualProj = no` is needed: with the
velocity residual projection on, the Newton run stopped with
`k_normalize non-finite norm` at 24, 32 and 40 ranks on Jean Zay (ifort 2021.9) and
converged at 8 and 16. With it off, the Newton residuals at 8, 16 and 40 ranks are
equal to the seven printed digits at every iteration, and the base flow at 16 and 40
ranks differs from the one at 8 ranks by 1e-11 (relative l2 norm of the velocity).
The run takes 956 s on 8 ranks, 477 s on 16 and 220 s on 40.

## Reference (ref/)
- `reference.json` checks six quantities. Three come from `residu_newton.dat`: the first
  residual (6.08e-7, tolerance 2 %), the iteration count (5, tolerance 1) and the final
  residual (below the 1e-9 tolerance). Three come from `BF_naca00120.f00001`: the norms
  of `vx` and `vy` and the largest `vx` (tolerance 1e-3). The first residual and the
  iteration count fail if the run starts from a field that is already converged.
  Verify with
  `scripts/check_against_ref.py example/naca0012_re2000/210_baseflow_newton/fp`.
- `ref/ic.png`, `ref/residual.png`, `ref/bf.png`: seed, residual history and base flow.

## Related
The stability stages `../../310_stability_direct` and `../../320_stability_adjoint`
ship their own base flow. Their references hold for that field; a base flow made by
this stage gives a growth rate that differs by the amounts in the table.
