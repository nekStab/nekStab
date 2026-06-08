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
`startFrom = BF_naca00120.f00001` — the base-flow guess (shipped). Newton
converges the steady state in place. `SPGnaca00120.f00001` is the sponge field.

## Mesh
4368 elements (`lelg=4368`), polynomial order N=5 (`lx1=6`).

## Run
```bash
mks naca0012
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — converged Newton residual (`residu_newton.dat` col 6, last
  iteration → 1e-9). Verify with
  `scripts/check_against_ref.py example/naca0012_re2000/210_baseflow_newton/fp`.

The converged `BF_naca00120.f00001` is the input to the stability stages
`../../310_stability_direct/direct` and `../../320_stability_adjoint`.
