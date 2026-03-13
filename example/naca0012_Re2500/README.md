# NACA 0012 Re = 2500 — Direct stability analysis

## Physics
2D flow around a NACA 0012 airfoil at Re = 2500. Computes the leading eigenvalues of the linearized operator around the steady base flow.

## nekStab Mode
`userParam01 = 3.1` — Direct stability
- `userParam07 = 220` — Krylov subspace dimension
- Sponge: left = 1, right = 5, strength = 1.7

## Prerequisites
- Base flow: `BF_naca00120.f00001`

## Run
```bash
mks naca0012        # Compile
nekbmpi naca0012 N  # Run on N MPI ranks
```

## Expected Output
Leading eigenvalues at Re = 2500.

## Reference
Internal example.
