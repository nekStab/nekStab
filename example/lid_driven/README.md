# Lid-Driven Cavity — Newton-GMRES base flow

## Physics
2D lid-driven cavity at Re = 3600 with aspect ratio 1.5 (`userParam10 = 1.5`). Computes the steady base flow using Newton-GMRES.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 90` — Krylov subspace dimension
- `userParam10 = 1.5` — cavity aspect ratio

## Prerequisites
- Initial guess: `BF_cav0.f00001`

## Run
```bash
mks cav        # Compile
nekbmpi cav N  # Run on N MPI ranks
```

## Expected Output
Converged 2D steady base flow at Re = 3600, aspect ratio 1.5.

## Reference
Standard lid-driven cavity benchmark.
