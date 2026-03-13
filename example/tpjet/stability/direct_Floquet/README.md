# Triple-Port Jet Direct Floquet — Floquet stability of forced jet

## Physics
Axisymmetric triple-port jet at Re = 1900. Direct Floquet analysis of the forced periodic orbit (St = 0.6, T = 1.667) to detect secondary instabilities.

## nekStab Mode
`userParam01 = 3.11` — Direct Floquet
- `userParam05 = 0.6` — forcing frequency (St_D)
- `userParam07 = 48` — Krylov subspace dimension
- `axiSymmetry = yes`

## Prerequisites
- Periodic base flow: `BF_tpjet0.f00001` (from `../../baseflow/`)

## Run
```bash
mks tpjet        # Compile
nekbmpi tpjet N  # Run on N MPI ranks
```

## Expected Output
Floquet multipliers at Re = 1900 (above Re_c ~ 1371).

## Reference
Internal example.
