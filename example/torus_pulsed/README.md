# Torus Pulsed — Newton-GMRES for pulsatile toroidal pipe flow

## Physics
Pulsatile flow in a toroidal pipe at Re = 1700 with cyclic boundary conditions. Newton-GMRES for forced periodic orbits (T ~ 8.72).

## nekStab Mode
`userParam01 = 2.2` — Newton-GMRES for forced periodic orbits
- `userParam07 = 200` — Krylov subspace dimension
- `endTime = 8.7195` — forcing period
- `cyclicBoundaries = yes`

## Prerequisites
- Base flow: `BF_torus0.f00001`

## Run
```bash
mks torus        # Compile
nekbmpi torus N  # Run on N MPI ranks
```

## Expected Output
Converged pulsatile periodic orbit at Re = 1700.

## Reference
Internal example.
