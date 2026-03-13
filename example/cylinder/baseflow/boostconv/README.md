# Cylinder BoostConv — Base flow via BoostConv acceleration

## Physics
2D flow around a circular cylinder at Re = 50. Computes the unstable steady base flow using BoostConv.

## nekStab Mode
`userParam01 = 1.2` — BoostConv

## Prerequisites
- Initial guess: `BFRe40_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow at Re = 50.

## Reference
Citro et al. (2017), JFM 813.
