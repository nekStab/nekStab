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

## Plot Notes
- Current `plot.py` shows velocity magnitude; paper plots **vx** (axial velocity, streamwise)
- Axisymmetric domain: axis labels should be z, r (not x, y); plot bounds: z ∈ [0, 40], r ∈ [0, 2]
- Base flow: Blues cmap 0 → 5; Floquet mode: RdBu ±2
- Leading Floquet multiplier μ ≈ −1 → period-doubling bifurcation at Re_c ≈ 1371
- Full paper figure requires multi-Re sweep (Re = 1300, 1370, 1375, 2000)
- Reference data: Shabani et al. (2017) — `shab_2000.dat`

## Reference
Internal example.
