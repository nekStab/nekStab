# Thermosyphon direct stability — Ra=500

## Physics
Direct global linear stability of the thermosyphon steady base flow at
Ra = 500 via Krylov-Schur. Just above Ra_c ≈ 494 the leading eigenvalue is a
real, steady (ω = 0) monotonically growing mode.

## nekStab Mode
- `userParam01 = 3.1` — direct stability
- `userParam06 = 500.0` — Rayleigh number Ra

## Initial condition
`startFrom = BF_tsyphon0.f00001` — the steady base flow from
`../../210_baseflow_newton` (shipped here).

## Run
```bash
mks tsyphon
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — leading eigenvalue (σ = 0.10222, ω = 0) from
  `Spectre_NSd.dat` row 1. Verify with
  `scripts/check_against_ref.py example/thermosyphon_ra500/310_stability_direct/direct`.
- `ref/baseflow.png`, `ref/mode.png`, `ref/spectrum.png`.

Outputs: leading direct eigenmodes `dRetsyphon0.f0000{1,2}` (shipped).
