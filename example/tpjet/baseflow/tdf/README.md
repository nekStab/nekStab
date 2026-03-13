# Triple-Port Jet TDF — Time-delayed feedback base flow

## Physics
Axisymmetric triple-port jet at Re = 2005. Computes the periodic base flow using time-delayed feedback (TDF) at the forced frequency St = 0.60.

## nekStab Mode
`userParam01 = 1.4` — TDF (time-delayed feedback)
- `userParam05 = 0.60` — forcing frequency (St_D)
- `axiSymmetry = yes`

## Prerequisites
- Initial guess: `BF_Re1900_tpjet0.f00001`

## Run
```bash
mks tpjet        # Compile
nekbmpi tpjet N  # Run on N MPI ranks
```

## Expected Output
Converged periodic base flow stabilized by time-delayed feedback.

## Reference
Internal example.
