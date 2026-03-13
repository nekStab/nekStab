# Cylinder Direct Stability — Direct eigenvalue problem

## Physics
2D flow around a circular cylinder at Re = 50. Computes the leading eigenvalues of the linearized Navier-Stokes operator around the steady base flow.

## nekStab Mode
`userParam01 = 3.1` — Direct stability
- `userParam03 = 2` — number of target Schur eigenvalues
- `userParam07 = 40` — Krylov subspace dimension

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (from `../baseflow/`)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Leading eigenvalues. The dominant mode crosses the imaginary axis near Re_c ~ 46.6 with St ~ 0.12.

## Reference
Barkley & Henderson (1996), JFM 322.
