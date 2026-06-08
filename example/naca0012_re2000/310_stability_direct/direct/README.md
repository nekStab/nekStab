# NACA 0012 direct global stability — Re=2000

## Physics
Direct global linear stability of the NACA 0012 steady base flow at Re = 2000
via Krylov-Schur on the time-stepper. The leading eigenvalue is stable at this
Reynolds number.

## nekStab Mode
- `userParam01 = 3.1` — direct stability
- `viscosity = -2000.0` — Re = 2000

## Initial condition
`startFrom = BF_naca00120.f00001` — the steady base flow from
`../../210_baseflow_newton/fp` (shipped here).

## Run
```bash
mks naca0012
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — leading eigenvalue (σ, ω) from `Spectre_NSd.dat` row 1
  (σ = -0.11272, ω = 8.03657). Verify with
  `scripts/check_against_ref.py example/naca0012_re2000/310_stability_direct/direct`.

Outputs: leading direct eigenmodes `dRenaca00120.f0000{1,2}` (shipped).
