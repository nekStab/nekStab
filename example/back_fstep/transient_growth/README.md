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

## Reference
Blackburn et al. (2008), JFM 603.
