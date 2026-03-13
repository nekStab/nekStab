# Cylinder Direct Floquet — Floquet stability of periodic orbit

## Physics
2D flow around a circular cylinder at Re = 50. Floquet analysis of the vortex shedding limit cycle to detect secondary instabilities.

## nekStab Mode
`userParam01 = 3.11` — Direct Floquet
- `userParam07 = 100` — Krylov subspace dimension
- Sponge: left = 5, right = 5, strength = 1.7

## Prerequisites
- UPO file: `BF_1cyl0.f00001` (periodic orbit, endTime adjusted from file)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Floquet multipliers. Mode A (|mu| > 1) appears near Re ~ 189, mode B near Re ~ 259.

## Reference
Barkley & Henderson (1996), JFM 322.
