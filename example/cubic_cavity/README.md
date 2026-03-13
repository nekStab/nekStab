# Cubic Cavity — Newton-GMRES base flow

## Physics
3D lid-driven cubic cavity at Re = 2500. Computes the steady base flow using Newton-GMRES.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension

## Prerequisites
- Initial guess: `BF_cav0.f00001`

## Run
```bash
mks cav        # Compile
nekbmpi cav N  # Run on N MPI ranks
```

## Expected Output
Converged 3D steady base flow at Re = 2500.

## Reference
Standard 3D lid-driven cavity benchmark.
