# Cylinder Newton Dynamic — Newton-GMRES with dynamic solver tolerances

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady base flow using Newton-GMRES with dynamic solver tolerances (`ifdyntol = .true.`). Early Newton iterations use relaxed GMRES tolerances (since the outer residual is large), tightening as Newton converges. This reduces total linear solver calls compared to the fixed-tolerance Newton case in `newton/`.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension
- `ifdyntol = .true.` — dynamic solver tolerances (set in `.usr`)

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow at Re = 100. Fewer total GMRES calls than the fixed-tolerance case (`newton/`). Comparison plot (`plot.py`) overlays both convergence histories.

## Reference
Loiseau et al. (2019), in Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics.
