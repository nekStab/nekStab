# Thermosyphon Ra=500 wavemaker, base-flow sensitivity, and PKE budget

This stage post-processes the committed direct and adjoint thermal stability
modes for the annular thermosyphon at Ra = 500. It is assembled from the
validated `320_stability_adjoint` case files and uses `userParam01 = 4`, the
same post-processing dispatch as the cylinder 411 model.

## Inputs

- Base flow: `BF_tsyphon0.f00001`
- Direct modes: `dRetsyphon0.f0000{1,2}` and `dImtsyphon0.f0000{1,2}`
- Adjoint modes: `aRetsyphon0.f0000{1,2}` and `aImtsyphon0.f0000{1,2}`

The direct and adjoint mode fields include velocity, pressure, and temperature.
The stage-local `maxmodes` is set to 2 because those are the shipped mode
snapshots available to the PKE budget loop.

## Run

```bash
mks tsyphon
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
python3 plot.py
python3 scripts/check_against_ref.py example/thermosyphon_ra500/411_postproc_wavemaker
```

## Outputs

- `wm_tsyphon0.f00001`: structural-sensitivity wavemaker field
- `sr_*`, `si_*`, `pr_*`, `pi_*`, `tr_*`, `ti_*`: base-flow sensitivity fields
- `K01*`, `K02*`, `PKE_dRe*`: PKE budget fields and integrals
- `plot_wavemaker.png`, `plot_direct_mode.png`, `plot_adjoint_mode.png`
- `wm_metrics.dat`: scalar reference inputs used by `ref/reference.json`
