# Cylinder Newton UPO — Periodic orbit via Newton-GMRES

## Physics
2D flow around a circular cylinder at Re = 50. Computes the periodic orbit (stable limit cycle from the supercritical Hopf bifurcation at Re ≈ 46.6) using Newton-GMRES for time-periodic solutions.

## nekStab Mode
`userParam01 = 2.1` — Newton-GMRES for UPOs
- `userParam07 = 100` — Krylov subspace dimension
- `endTime = 7.9` — approximate orbit period

## Prerequisites
- Restart file: `rstcyl0.f00001` (initial guess near the periodic orbit)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged periodic orbit. Newton residual and period correction decrease to machine precision.

## Reference
Loiseau et al. (2019), in Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics.
