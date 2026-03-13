# Cylinder SFD — Selective frequency damping

## Physics
2D flow around a circular cylinder at Re = 50. Computes the unstable steady base flow using SFD.

## nekStab Mode
`userParam01 = 1.1` — SFD
- `userParam04 = 0.12` — filter frequency (Strouhal of leading mode)
- `userParam05 = 0.05` — gain (twice growth rate of leading mode)

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001` (base flow at nearby Re)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow `BF_1cyl0.f00001` at Re = 50. Residual decreases monotonically.

## Reference
Akervik et al. (2006), Phys. Fluids 18.
