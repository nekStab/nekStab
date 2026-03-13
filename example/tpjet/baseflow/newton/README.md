# Triple-Port Jet Newton — Forced periodic orbit via Newton-GMRES

## Physics
Axisymmetric triple-port jet at Re = 1900. Computes the forced periodic orbit using Newton-GMRES at the forcing frequency St = 0.6 (T = 1.667).

## nekStab Mode
`userParam01 = 2.2` — Newton-GMRES for forced periodic orbits
- `userParam05 = 0.6` — forcing frequency (St_D)
- `userParam07 = 30` — Krylov subspace dimension (GMRES for Newton linear solve)
- `axiSymmetry = yes`

## Prerequisites
- Initial guess: `BF_tpjet0.f00001`

## Run
```bash
mks tpjet        # Compile
nekbmpi tpjet N  # Run on N MPI ranks
```

## Expected Output
Converged forced periodic orbit at Re = 1900 (above Re_c ~ 1371).

## Reference
Internal example.
