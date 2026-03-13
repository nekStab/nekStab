# Cylinder OTD — Optimally time-dependent modes

## Physics
2D flow around a circular cylinder at Re = 180. OTD decomposition tracks instantaneous flow instabilities.

## nekStab Mode
`userParam01 = 5` — OTD

## Prerequisites
- Base flow: `BF_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
OTD modes and instantaneous Lyapunov exponents. Krylov dimension k = 92.

## Reference
Babaee & Sapsis (2016), Proc. R. Soc. A.
