# Torus Full — Direct stability of full toroidal pipe

## Physics
Flow in a full toroidal pipe at Re = 1650 with cyclic boundary conditions. Direct stability analysis on the complete torus geometry.

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
Leading eigenvalues for the full toroidal pipe flow at Re = 1650.

## Reference
Internal example.
