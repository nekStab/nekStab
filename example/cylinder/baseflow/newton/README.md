# Cylinder Newton — Base flow via Newton-GMRES

## Physics
2D flow around a circular cylinder at Re = 50. Computes the unstable steady base flow using Newton-GMRES iteration.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 30` — Krylov subspace dimension

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow. Newton residual decreases quadratically near convergence.

## Reference
Loiseau et al. (2019), in Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics.
