# Thermosyphon Baseflow — Newton-GMRES with buoyancy

## Physics
Axisymmetric thermosyphon at Ra = 500, Pr = 5. Computes the steady base flow (conductive state) using Newton-GMRES with coupled velocity-temperature.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam06 = 500.0` — Rayleigh number
- `userParam07 = 100` — Krylov subspace dimension

## Prerequisites
- Initial guess: `BF_Ra400_tsyphon0.f00001` (base flow at Ra = 400)

## Run
```bash
mks tsyphon        # Compile
nekbmpi tsyphon N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow at Ra = 500 (slightly above Ra_c ~ 494).

## Plot Notes
- Current `plot.py` shows velocity magnitude; paper shows base flow **temperature** field
- Add temperature panel alongside velocity magnitude to match paper figure

## Reference
Soucasse et al. (2019), PRF 4.
