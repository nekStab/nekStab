# NACA 0012 adjoint global stability — Re=2000

## Physics
Adjoint global linear stability of the NACA 0012 steady base flow at Re = 2000
via Krylov-Schur. The adjoint eigenvalue matches the direct one (σ, ω agree to
the convergence tolerance), as expected; the adjoint mode localizes the receptivity.

## nekStab Mode
- `userParam01 = 3.2` — adjoint stability
- `viscosity = -2000.0` — Re = 2000

## Initial condition
`startFrom = BF_naca00120.f00001` — the steady base flow from
`../210_baseflow_newton/fp` (shipped here).

## Run
```bash
mks naca0012
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — leading adjoint eigenvalue (σ, ω) from `Spectre_NSa.dat`
  row 1 (σ = -0.11367, ω = 8.03673). Verify with
  `scripts/check_against_ref.py example/naca0012_re2000/320_stability_adjoint`.

Outputs: leading adjoint eigenmodes `aRenaca00120.f0000{1,2}` (shipped).
