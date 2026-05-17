# Time-Periodic Jet Newton — Forced periodic orbit via Newton-GMRES

## Physics
Axisymmetric time-periodic jet at Re = 1900. Computes the forced periodic orbit using Newton-GMRES at the forcing frequency St = 0.6 (T = 1.667).

## nekStab Mode
`userParam01 = 2.2` — Newton-GMRES for forced periodic orbits
- `userParam05 = 0.6` — forcing frequency (St_D)
- `userParam07 = 30` — Krylov subspace dimension (GMRES for Newton linear solve)
- `axiSymmetry = yes`

## Prerequisites
- Initial guess: `BF_tpjet0.f00001`

## Run
```bash
mks tpjet        # Compile
nekbmpi tpjet N  # Run on N MPI ranks
```

## Expected Output
Converged forced periodic orbit at Re = 1900 (above Re_c ~ 1371).

## Latest Result
Verified 2026-05-17 on 8 ranks (Slurm):

- 1 Newton outer iteration (already near converged UPO), |F| = 9.3e-9
- Wall time 8.5 s
- BF_tpjet0.f00001 refreshed (6.9 MB) and committed as proof

## Plot Notes
- Current `plot.py` shows velocity magnitude; paper plots **vx** (axial velocity, Blues cmap 0 → 5)
- Axisymmetric domain: axis labels should be z, r (not x, y); plot bounds: z ∈ [0, 40], r ∈ [0, 2]

## Reference
Internal example.
