# Thermosyphon Baseflow — Newton-GMRES with buoyancy

## Physics
Axisymmetric thermosyphon at Ra = 500, Pr = 5. Computes the steady base flow (conductive state) using Newton-GMRES with coupled velocity-temperature.

The thermosyphon is not scaled like the thermal cylinder case. Here the buoyancy
forcing entering momentum is proportional to `Pr * Ra * T`, so the automatic
thermal norm weighting must use that combined coefficient rather than `Ra` alone.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam06 = 500.0` — Rayleigh number
- `userParam07 = 100` — Krylov subspace dimension
- thermal norm weighting uses:
  `thermal_norm_mode = auto_clipped`
  `thermal_buoyancy_coeff = Pr * Ra`

## Prerequisites
- Initial guess: `BF_Ra400_tsyphon0.f00001` (base flow at Ra = 400)

## Run
```bash
mks tsyphon        # Compile
nekbmpi tsyphon N  # Run on N MPI ranks
```

The directory must contain a standard `SESSION.NAME` file. Without it, Nek5000
fails before reading the case with:
`ERROR:  Error while reading SESSION.NAME!`

## Expected Output
Converged steady base flow at Ra = 500 (slightly above Ra_c ~ 494).

## Latest Result

Verified 2026-05-17 on 4 MPI ranks (Slurm):

- automatic thermal norm setup active (`thermal_norm_mode = TN_CLIP`)
- `thermal_buoyancy_coeff = Pr * Ra = 2.5e3`
- 10 Newton outer iterations, |F| = 7.5e-12
- Wall time 43 s
- Outputs: `BF_tsyphon0.f00001`, `restsyphon0.f00001`

## Plot Notes
- Current `plot.py` shows velocity magnitude; paper shows base flow **temperature** field
- Add temperature panel alongside velocity magnitude to match paper figure

## Reference
Soucasse et al. (2019), PRF 4.
