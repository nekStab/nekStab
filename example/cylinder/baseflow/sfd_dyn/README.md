# Cylinder SFD Dynamic — SFD with dynamic solver tolerances and OIFS

## Physics
2D flow around a circular cylinder at Re = 50. Steady base flow via SFD with dynamic solver tolerances (`ifdyntol = .true.` in `.usr`) and OIFS time stepping (CFL = 5).

## nekStab Mode
`userParam01 = 1.1` — SFD
- `userParam04 = 0.12` — filter frequency
- `userParam05 = 0.05` — gain
- `ifdyntol = .true.` — set in `.usr` for adaptive Nek5000 solver tolerances
- `extrapolation = OIFS`, `targetCFL = 5`

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow. Solver tolerances adapt during convergence. OIFS allows larger time steps than standard BDF3.

## Reference
Akervik et al. (2006), Phys. Fluids 18.
