# Flip-Flop Adjoint Floquet - Re=62

## Physics
2D flow around two side-by-side cylinders at Re = 62. This stage computes the
adjoint Floquet spectrum of the periodic flip-flop orbit with period T = 8.73356.
The leading adjoint conjugate pair matches the direct Floquet multiplier from
stage 311.

## nekStab Mode
- `userParam01 = 3.21` - adjoint Floquet analysis
- `userParam07 = 90` - Krylov subspace dimension, bumped from 50 for adjoint convergence
- `endTime = 8.73356` - orbit period from the Newton UPO
- `viscosity = -62.0` - Re = 62
- `startFrom = BF_2cyl0.f00001` - UPO base-flow orbit

## Base-Flow Replay
The base-flow orbit is built once and then replayed time-reversed for every
adjoint matvec. This keeps the adjoint Floquet operator aligned with the stored
periodic orbit used by the direct Floquet stage.

## Run
```bash
mks 2cyl
sbatch run.local.slurm
```

`run.local.slurm` refreshes `SESSION.NAME` for the current directory and writes
stdout/stderr to `logfile`.

## Latest Result (verified 2026-06-10)
- Reference produced with 16 ranks.
- Leading adjoint pair: sigma = 0.0078864, omega = +/-0.1407585.
- Direct 311 leading pair: sigma = 0.0078930, omega = +/-0.1407626.

## Reference (ref/)
- `reference.json` - leading adjoint Floquet pair from `Spectre_NSa.dat` row 1.
- `spectrum.png` - adjoint Floquet spectrum.
- `mode.png` - leading adjoint Floquet mode.
- `upo_snapshot.png` - UPO base-flow snapshot.

Verify a fresh run from the repository root with:

```bash
python3 scripts/check_against_ref.py example/flip_flop_re62/321_stability_adjoint_floquet
```
