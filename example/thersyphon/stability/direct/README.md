# Thermosyphon Direct Stability — Eigenvalue problem with buoyancy

## Physics
Axisymmetric thermosyphon at Ra = 500, Pr = 5. Direct stability analysis of the conductive base flow with coupled velocity-temperature perturbations.

## nekStab Mode
`userParam01 = 3.1` — Direct stability
- `userParam06 = 500.0` — Rayleigh number
- `userParam07 = 120` — Krylov subspace dimension

## Prerequisites
- Base flow: `BF_tsyphon0.f00001` (from `../../baseflow/`)

## Run
```bash
mks tsyphon        # Compile
nekbmpi tsyphon N  # Run on N MPI ranks
```

## Expected Output
Leading eigenvalues at Ra = 500 (slightly above Ra_c ~ 494). Onset of convective instability.

## Plot Notes
- Paper plots **temperature** field (index 5) for both base flow and eigenmode — current `plot.py` correctly uses temperature
- Pitchfork: real eigenvalue crosses σ = 0 at Ra_c,1 ≈ 494 (symmetry-breaking in θ)
- Full paper figure requires multi-Ra sweep (Ra = 480, 490, 500, 510) to show crossing
- Hopf bifurcation at Ra_c,2 ≈ 16081 (oscillating convection cell) — not covered here

## Reference
Soucasse et al. (2019), PRF 4.
