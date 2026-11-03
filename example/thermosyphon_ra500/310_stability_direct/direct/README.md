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
sbatch run.local.slurm     # 4 ranks
```

## Reference (ref/)
- `reference.json` — leading eigenvalue (σ = 0.10222, ω = 0) from
  `Spectre_NSd.dat` row 1. Verify with
  `scripts/check_against_ref.py example/thermosyphon_ra500/310_stability_direct/direct`.
- `ref/baseflow.png`, `ref/mode.png`, `ref/spectrum.png`.

## How many modes converge
The run stops on the leading value (`schur_tgt = 1`) after 14 matrix-vector products, about
15 s on 4 ranks. The next value is real too:
sigma = -0.2751999, measured with `schur_tgt = 2`, at a residual of 9.7e-7 against the
tolerance of 1e-6. That sits so close to the tolerance that the same run took 19 restarts
twice (33 min) and was not done after 49 (85 min), so this stage does not ask for it.

Versions before 2.0.0 also returned a value at sigma = 0. It was not a flow: a constant
pressure in a closed domain, with multiplier exactly 1. `DOC.md` (Troubleshooting) has the
cause and the fix.

Outputs: leading direct eigenmode `dRetsyphon0.f00001` (shipped).
