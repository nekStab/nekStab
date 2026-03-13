# Torus — Direct stability of toroidal pipe flow

## Physics
Flow in a toroidal pipe at Re = 1650 with cyclic boundary conditions. Direct stability analysis to find the leading eigenvalues.

## nekStab Mode
`userParam01 = 3.1` — Direct stability
- `userParam07 = 200` — Krylov subspace dimension
- `cyclicBoundaries = yes`

## Prerequisites
- Base flow: `BF_torus0.f00001`

## Run
```bash
mks torus        # Compile
nekbmpi torus N  # Run on N MPI ranks
```

## Expected Output
Leading eigenvalues for the toroidal pipe flow at Re = 1650.

## Reference
Internal example.
