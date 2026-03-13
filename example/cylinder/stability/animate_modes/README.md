# Cylinder Animate Modes — Eigenmode reconstruction

## Physics
2D flow around a circular cylinder at Re = 50. Reconstructs and animates the leading eigenmode over one oscillation period (T = 8.305).

## nekStab Mode
`userParam01 = 4.5` — Animate modes
- `userParam07 = 92` — Krylov subspace dimension
- `endTime = 8.305` — period of leading mode

## Prerequisites
- Base flow: `BF_1cyl0.f00001`
- Converged eigenmode files (from `../direct/`)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Sequence of field files showing the eigenmode oscillation over one period.

## Reference
Internal example.
