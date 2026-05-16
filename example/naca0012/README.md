# NACA 0012 — Newton-GMRES base flow

## Physics
2D flow around a NACA 0012 airfoil at Re = 2000. Computes the steady base flow using Newton-GMRES with sponge regions.

## nekStab Mode
`userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam07 = 220` — Krylov subspace dimension
- Sponge: left = 1, right = 5, strength = 1.7

## Prerequisites
- Initial guess: `BF_naca00120.f00001`

## Run
```bash
mks naca0012        # Compile
nekbmpi naca0012 N  # Run on N MPI ranks
```

## Mesh
4368 elements (`lelg=4368` in `SIZE`, `lpmin=8`).  Polynomial order N=5
(`lx1=6`).

## Run
```bash
mks naca0012
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Verified 2026-05-17 on 8 ranks (Slurm):

- 10 Newton outer iterations to |F| = 9.4e-10
- 27957 inner solver steps, 194 cumulative GMRES iters
- Wall time 722 s (≈ 12 min)
- Outputs: converged `BF_naca00120.f00001` (6.3 MB)

## Reference
Internal example.
