# Cylinder SFD OIFS — SFD with operator-integration-factor splitting

## Physics
2D flow around a circular cylinder at Re = 50. SFD base flow computation using OIFS for larger time steps (CFL = 5).

## nekStab Mode
`userParam01 = 1.1` — SFD
- `userParam04 = 0.12` — filter frequency
- `userParam05 = 0.05` — gain
- `extrapolation = OIFS`, `targetCFL = 5`

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow. OIFS allows CFL ~ 5, significantly faster than standard BDF3.

## Reference
Akervik et al. (2006), Phys. Fluids 18.
