# Cylinder Steady Force Sensitivity — Sensitivity to steady forcing

## Physics
2D flow around a circular cylinder at Re = 50. Computes the sensitivity of the eigenvalue to a steady body force, identifying where a small steady force most effectively modifies the instability.

## nekStab Mode
`userParam01 = 4.41` — Sensitivity to steady forcing
- `userParam07 = 128` — Krylov subspace dimension

## Prerequisites
- Direct and adjoint eigenmode files (from `../../stability/`)
- No restart file needed

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Sensitivity field showing optimal placement for steady forcing to control the instability.

## Reference
Marquet et al. (2008), JFM 615.
