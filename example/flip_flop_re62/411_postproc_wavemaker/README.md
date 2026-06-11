# Flip-flop Re=62 Floquet post-processing

This stage runs the time-periodic post-processing path for the flip-flop UPO at
Re = 62. It is assembled from the committed `321_stability_adjoint_floquet`
case files and uses `userParam01 = 4.11`, which dispatches
`energy_budget_floquet` over the UPO period `T = 8.73356`.

## Inputs

- Base-flow orbit seed: `BF_2cyl0.f00001`
- Sponge field: `SPG2cyl0.f00001`
- Direct Floquet modes: `dRe2cyl0.f0000{1,2}` and `dIm2cyl0.f0000{1,2}`
- Adjoint Floquet modes: `aRe2cyl0.f0000{1,2}` and `aIm2cyl0.f0000{1,2}`

The stage-local `maxmodes` is set to 2 because those are the direct Floquet
mode snapshots available to the budget loop.

## Run

```bash
mks 2cyl
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
python3 plot.py
python3 scripts/check_against_ref.py example/flip_flop_re62/411_postproc_wavemaker
```

## Outputs

- `F01*`, `F02*`: orbit-averaged Floquet budget fields
- `PKE_floquet_2cyl0.f0000{1,2}`: orbit-averaged budget integrals
- `plot_wavemaker.png`: direct/adjoint overlap proxy for the leading Floquet mode
- `plot_direct_mode.png`, `plot_adjoint_mode.png`, `plot_budget_field.png`
- `wm_metrics.dat`: scalar reference inputs used by `ref/reference.json`

## Latest Result

Verified 2026-06-11 on 8 ranks (Slurm), using the stage-local committed
direct/adjoint Floquet mode inputs:

- Nek run completed one UPO period, T = 8.73356; wall time: 40.37 s.
- Wavemaker metrics: max = 0.03216049123256942 at x = 2.4090390205383301,
  y = 0.6052520275115967; mean = 0.0005153653561364195.
- `python3 plot.py` regenerated `plot_wavemaker.png`,
  `plot_direct_mode.png`, `plot_adjoint_mode.png`, `plot_budget_field.png`,
  `wm_metrics.dat`, and the matching `ref/plot_*.png` copies.
- Upstream fresh `311`/`321` Floquet mode reruns are pending `SIZE` approval;
  this result refreshes the post-processing stage from its existing mode inputs.
