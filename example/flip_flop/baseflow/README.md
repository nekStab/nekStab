# Flip-Flop Baseflow — UPO via Newton-GMRES

## Physics
2D flow around two side-by-side circular cylinders at Re = 62. Computes the unstable periodic orbit (flip-flop instability) using Newton-GMRES (T ~ 8.734).

## nekStab Mode
`userParam01 = 2.1` — Newton-GMRES for UPOs
- `userParam07 = 30` — Krylov subspace dimension (GMRES for Newton linear solve)
- `endTime = 8.73356` — approximate orbit period

## Prerequisites
- Initial guess: `BF_Re60_2cyl0.f00001`

## Run
```bash
mks 2cyl        # Compile
nekbmpi 2cyl N  # Run on N MPI ranks
```

## Expected Output
Converged periodic orbit at Re = 62 (slightly above Re_c ~ 61.17).

## Reference
Carini et al. (2015), JFM 778.
