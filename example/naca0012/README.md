# NACA 0012 — Newton-GMRES base flow

## Physics
2D flow around a NACA 0012 airfoil at Re = 2000. Computes the steady base flow using Newton-GMRES with sponge regions.

## nekStab Mode
`userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam07 = 220` — Krylov subspace dimension
- Sponge: left = 1, right = 5, strength = 1.7

## Prerequisites
- Initial guess: `BF_naca00120.f00001`

## Run
```bash
mks naca0012        # Compile
nekbmpi naca0012 N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow at Re = 2000.

## Reference
Internal example.
