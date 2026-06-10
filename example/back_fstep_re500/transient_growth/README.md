# Backward-Facing Step Transient Growth — Optimal perturbation

## Physics
2D flow over a backward-facing step at Re = 500. Computes the optimal initial perturbation for maximum transient energy growth.

## nekStab Mode
`userParam01 = 3.3` — Transient growth
- `userParam07 = 64` — Krylov subspace dimension
- Sponge: left = 5, right = 10, strength = 2

## Prerequisites
- Base flow: `BF_bfs0.f00001` (from `../baseflow/`)

## Run
```bash
mks bfs        # Compile
nekbmpi bfs N  # Run on N MPI ranks
```

## Expected Output
Optimal gain G(t) and the associated initial/response perturbation fields.

## Plot Notes
- Current `plot.py` shows **vy**; paper plots **vx** (shows Kelvin-Helmholtz structure along shear layer)
- Paper uses 3 panels: (a) base flow vx (cividis), (b) optimal perturbation vx (RdBu), (c) optimal response vx (seismic)
- Plot bounds: x ∈ [−5, 25], y ∈ [−1, 1]; step body: rectangle at (−5, −1) size 5 × 1
- G(τ) envelope from `Spectre_Hp.dat` files; reference data: Barkley et al. (2008) Fig. 5
- Full paper figure requires multi-τ sweep (τ = 8, 18, 38, 58, 88); peak at τ = 58

## Reference
Blackburn et al. (2008), JFM 603.
