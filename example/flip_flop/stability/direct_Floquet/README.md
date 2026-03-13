# Flip-Flop Direct Floquet — Floquet stability of flip-flop orbit

## Physics
2D flow around two side-by-side cylinders at Re = 62. Floquet analysis of the periodic flip-flop orbit (T = 8.734) to detect secondary instabilities.

## nekStab Mode
`userParam01 = 3.11` — Direct Floquet
- `userParam07 = 48` — Krylov subspace dimension
- `endTime = 8.73356` — orbit period

## Prerequisites
- UPO file: `BF_2cyl0.f00001` (from `../../baseflow/`)

## Run
```bash
mks 2cyl        # Compile
nekbmpi 2cyl N  # Run on N MPI ranks
```

## Expected Output
Floquet multipliers of the flip-flop periodic orbit.

## Reference
Carini et al. (2015), JFM 778.
