# Cylinder Newton Dynamic — Newton-GMRES with dynamic solver tolerances

## Physics
2D flow around a circular cylinder at Re = 50. Computes the unstable steady base flow using Newton-GMRES with dynamic solver tolerances (`ifdyntol = .true.` in `.usr`) and tighter PDE solver tolerances (1e-11).

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension
- `ifdyntol = .true.` — set in `.usr` for adaptive Nek5000 solver tolerances

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow at Re = 50 with high-precision solver tolerances.

## Reference
Loiseau et al. (2019), in Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics.
