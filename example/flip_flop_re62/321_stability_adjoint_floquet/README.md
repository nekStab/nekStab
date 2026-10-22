# Flip-Flop Adjoint Floquet - Re=62

## Physics
2D flow around two side-by-side cylinders at Re = 62. This stage computes the
adjoint Floquet spectrum of the periodic flip-flop orbit with period T = 8.815172.
The leading adjoint conjugate pair matches the direct Floquet multiplier from
stage 311.

## nekStab Mode
- `userParam01 = 3.21` - adjoint Floquet analysis
- `userParam07 = 90` - Krylov subspace dimension, bumped from 50 for adjoint convergence
- `endTime = 8.815172` - orbit period. The run takes the period from the time stamp of `BF_2cyl0.f00001`
  and replaces `endTime` with it. The first version of this stage carried the starting value of the Newton iteration, 8.73356.
- `viscosity = -62.0` - Re = 62
- `startFrom = BF_2cyl0.f00001` - UPO base-flow orbit

## Base-Flow Replay
The base-flow orbit is built once and then replayed time-reversed for every
adjoint matvec. This keeps the adjoint Floquet operator aligned with the stored
periodic orbit used by the direct Floquet stage.

## Run
```bash
mks 2cyl
sbatch --ntasks=16 --mem-per-cpu=2000 run.local.slurm
```

`run.local.slurm` refreshes `SESSION.NAME` for the current directory and writes
stdout/stderr to `logfile`.

## Result (2026-10-09, Jean Zay, `ifort`, 16 ranks)
- Leading adjoint pair: sigma = 0.007889578, omega = +/-0.1407610. The run of 2026-06-11 gave the same values.
- Direct 311 leading pair: sigma = 0.007896230, omega = +/-0.1407652. The adjoint sigma is 6.7e-6 (0.08 %) lower.
- Wall time = 1453 s.

## Reference (ref/)
- `reference.json` - leading adjoint Floquet pair from `Spectre_NSa.dat` row 1.
- `spectrum.png` - adjoint Floquet spectrum.
- `mode.png` - leading adjoint Floquet mode.
- `upo_snapshot.png` - UPO base-flow snapshot.

Verify a fresh run from the repository root with:

```bash
python3 scripts/check_against_ref.py example/flip_flop_re62/321_stability_adjoint_floquet
```
