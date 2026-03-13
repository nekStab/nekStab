# Cylinder Adjoint Stability — Adjoint eigenvalue problem

## Physics
2D flow around a circular cylinder at Re = 50. Computes the leading adjoint eigenvalues, identifying regions of maximum receptivity.

## nekStab Mode
`userParam01 = 3.2` — Adjoint stability
- `userParam07 = 48` — Krylov subspace dimension
- Sponge: left = 5, right = 5, strength = 1.7

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (from `../baseflow/`)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Adjoint eigenvalues and modes. Adjoint modes peak upstream of the cylinder (receptivity region).

## Reference
Giannetti & Luchini (2007), JFM 581.
