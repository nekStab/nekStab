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

## Latest Result
- 2026-06-11, 8 ranks, `k_dim = 220`, sponge strength 1.7: σ = -0.0665164, ω = 8.00577.
- Sponge-fix provenance: recomputed after the sponge-strength correction.

## Reference (ref/)
- `reference.json` — leading adjoint eigenvalue (σ, ω) from `Spectre_NSa.dat`
  row 1 (σ = -0.0665164, ω = 8.00577). Verify with
  `scripts/check_against_ref.py example/naca0012_re2000/320_stability_adjoint`.

Outputs: leading adjoint eigenmodes `aRenaca00120.f0000{1,2}` (shipped).

Gate run of 2026-10-09 (Jean Zay, `ifort`, 8 ranks): σ = -0.0665164, ω = 8.00577, wall time 1334 s.
`plot.py` writes `plot_bf.png` (velocity magnitude of the base flow), `plot_spectrum.png` (converged
eigenvalues, growth rate against frequency) and `plot_mode.png` (vy of the leading adjoint mode, which
sits at the trailing edge); `ref/` holds copies as `bf.png`, `spectrum.png` and `mode.png`.
