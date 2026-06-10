# Lid-Driven Cavity — Newton-GMRES base flow

## Physics
2D lid-driven cavity at Re = 3600 with aspect ratio 1.5 (`userParam10 = 1.5`). Computes the steady base flow using Newton-GMRES.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 90` — Krylov subspace dimension
- `userParam10 = 1.5` — cavity aspect ratio

## Prerequisites
- Initial guess: `BF_cav0.f00001`

## Mesh
100 elements (`lelg=100` in `SIZE`). Compact cavity mesh — runs fast.

## Run
```bash
mks cav
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Verified 2026-05-17 on 8 ranks (Slurm):

- 1 Newton outer iteration (already near converged BF), |F| = 2.5e-10
- Wall time 0.58 s
- Output: `BF_cav0.f00001` (144 KB)

## Reference
Standard lid-driven cavity benchmark.
