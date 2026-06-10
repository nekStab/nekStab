# Thermosyphon Adjoint Stability - Ra=500

## Physics
Thermal convection in the annular thermosyphon at Ra = 500 and Pr = 5. This
stage computes the leading steady adjoint eigenvalue of the Newton base flow.
The adjoint eigenvalue matches the direct 310 stability value, as expected for
the direct and adjoint operators.

## nekStab Mode
- `userParam01 = 3.2` - steady adjoint stability
- `userParam06 = 500.0` - Rayleigh number
- `userParam07 = 120` - Krylov subspace dimension
- `startFrom = BF_tsyphon0.f00001` - steady Newton base flow
- `endTime = 0.5`

## Run
```bash
mks tsyphon
sbatch run.local.slurm
```

`run.local.slurm` refreshes `SESSION.NAME` for the current directory and writes
stdout/stderr to `logfile`.

## Latest Result (verified 2026-06-10)
- Reference produced with 8 ranks.
- Leading adjoint eigenvalue: sigma = 0.1022272, omega = 0.0.
- Direct 310 leading eigenvalue: sigma = 0.1022191.

## Reference (ref/)
- `reference.json` - leading adjoint eigenvalue from `Spectre_NSa.dat` row 1.
- `spectrum.png` - adjoint spectrum.
- `mode.png` - leading adjoint mode.
- `baseflow.png` - steady base-flow temperature.

Verify a fresh run from the repository root with:

```bash
python3 scripts/check_against_ref.py example/thermosyphon_ra500/320_stability_adjoint
```
