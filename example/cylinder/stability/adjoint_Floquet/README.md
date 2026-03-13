# Cylinder Adjoint Floquet — Adjoint Floquet stability

## Physics
2D flow around a circular cylinder at Re = 50. Adjoint Floquet analysis of the periodic orbit for receptivity of secondary instabilities.

## nekStab Mode
`userParam01 = 3.21` — Adjoint Floquet
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
Adjoint Floquet multipliers and modes. Identifies receptivity regions for secondary instabilities.

## Reference
Barkley & Henderson (1996), JFM 322.
