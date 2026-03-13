# Backward-Facing Step Baseflow — Newton-GMRES

## Physics
2D flow over a backward-facing step at Re = 500. Computes the steady base flow using Newton-GMRES.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 30` — Krylov subspace dimension (GMRES for Newton linear solve)

## Prerequisites
- Initial guess: `BF_bfs0.f00001`

## Run
```bash
mks bfs        # Compile
nekbmpi bfs N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow at Re = 500. Recirculation region behind the step.

## Reference
Barkley et al. (2002), JFM 473; Blackburn et al. (2008), JFM 603.
