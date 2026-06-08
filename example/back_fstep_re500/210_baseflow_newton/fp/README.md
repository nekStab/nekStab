# Backward-facing step base flow at Re=500 — Newton-GMRES

## Physics
2D backward-facing step at Re = 500. Computes the steady base flow with
Newton-GMRES (the classic Barkley/Blackburn/Sherwin separated-flow benchmark).

## nekStab Mode
- `userParam01 = 2.0` — Newton-GMRES for fixed points
- `viscosity = -500.0` — Re = 500

## Initial condition
`startFrom = BF_bfs0.f00001` — the base-flow guess (shipped). `rstbfs0.f00001`
is the upstream restart/seed.

## Run
```bash
mks bfs
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
Figure/IC reference only — the shipped residual files were empty/corrupt, so
`reference.json` asserts no numeric scalar. A fresh run regenerates
`residu_newton.dat` (convergence target 1e-8); compare the converged base flow
visually against `ref/bf.png`.
- `ref/bf.png` base flow, `ref/ic.png` initial guess, `ref/residual.png`.

## Related
The `330_transient_growth` stage for this case is under investigation (TG
diverged to NaN) and is not yet migrated.
