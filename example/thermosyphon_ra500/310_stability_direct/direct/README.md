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
Two values converge, both real: sigma = 0.10222 and sigma = 0 to round-off
(1e-13; -2e-8 in another run). A run with `k_dim = 300` and `schur_tgt = 6` did 130
restarts of the 300-vector basis (4.3 h on 4 ranks) and still had these two only. The
next value, near sigma = -0.276, stays at a residual near 2e-3. This run was stopped and
is not a stage of this folder.

The value at 0 is not explained. Its mode is a continuous field (`bin/c0check` finds no
jump between elements), so it is not the element-jump mode of `DOC.md`. It is mostly
velocity (rms 0.20 for each component against a temperature rms of 1.5e-4). It is not a
rigid rotation of the annulus. The reference checks the leading value only.

Outputs: leading direct eigenmodes `dRetsyphon0.f0000{1,2}` (shipped).
