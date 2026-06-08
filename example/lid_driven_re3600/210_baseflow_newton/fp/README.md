# Lid-driven cavity base flow at Re=3600 — Newton-GMRES

## Physics
2D lid-driven cavity (aspect ratio 1.5) at Re = 3600. Computes the steady
(unstable) base flow with Newton-GMRES.

## nekStab Mode
- `userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam10 = 1.5` — cavity aspect ratio
- `viscosity = -3600.0` — Re = 3600

## Initial condition
`startFrom = BF_cav0.f00001` — the base-flow guess (shipped). Newton converges
the steady state in place.

## Run
```bash
mks cav
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — converged Newton residual (`residu_newton.dat` col 6, last
  iteration → 1e-9). Verify with
  `scripts/check_against_ref.py example/lid_driven_re3600/210_baseflow_newton/fp`.
- `ref/bf.png` base flow, `ref/ic.png` initial guess, `ref/residual.png` decay.
