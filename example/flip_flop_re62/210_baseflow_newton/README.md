# Flip-Flop UPO at Re=62 — Newton-GMRES

## Physics
2D flow around two side-by-side circular cylinders at Re = 62 (slightly above
Re_c ≈ 61.17). Computes the unstable periodic orbit (UPO) of the flip-flop
instability using Newton-GMRES. Converged orbit period T ≈ 8.73356.

## nekStab Mode
- `userParam01 = 2.1` — Newton-GMRES for UPOs
- `userParam07 = 30` — Krylov subspace dimension (GMRES for Newton linear solve)
- `endTime = 8.73356` — approximate orbit period
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

## Reference (ref/)
- `reference.json` — converged Newton residual (`residu_newton.dat` col 6, last
  iteration → 1e-8). Verify a fresh run with
  `scripts/check_against_ref.py example/flip_flop_re62/210_baseflow_newton`.
- `ref/bf.png` — converged UPO base flow
- `ref/ic.png` — initial guess
- `ref/residual.png` — Newton residual decay

The converged orbit `BF_2cyl0.f00001` is the input to
`../311_stability_direct_floquet`.

## Reference
Carini et al. (2015), JFM 778.
