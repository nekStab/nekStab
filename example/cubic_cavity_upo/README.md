# Cubic Cavity UPO — Unstable periodic orbit

## Physics
3D lid-driven cubic cavity at Re = 1950. Computes an unstable periodic orbit using Newton-GMRES for time-periodic solutions (T ~ 10.77).

## nekStab Mode
`userParam01 = 2.1` — Newton-GMRES for UPOs
- `userParam07 = 200` — Krylov subspace dimension
- `endTime = 10.7738` — approximate orbit period

## Prerequisites
- Initial guess: `1960_PO_cav0.f00001` (near-periodic solution at Re = 1960)

## Run
```bash
mks cav        # Compile
nekbmpi cav N  # Run on N MPI ranks
```

## Expected Output
Converged 3D periodic orbit at Re = 1950 (Re_c ~ 1916).

## Reference
Standard 3D lid-driven cavity benchmark.
