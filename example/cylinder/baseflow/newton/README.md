# Cylinder Newton — Base flow via Newton-GMRES

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady base flow using Newton-GMRES iteration with fixed solver tolerances. Compare with `newton_dyn/` (dynamic tolerances) to see the reduction in total linear solver calls.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension (same as newton_dyn)
- `ifdyntol = .false.` — fixed solver tolerances (baseline for comparison)

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow. Newton residual decreases quadratically near convergence. Comparison plot (`plot.py`) overlays Newton vs Newton_dyn convergence.

## Reference
Loiseau et al. (2019), in Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics.
