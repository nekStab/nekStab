# NACA 0012 wavemaker, base-flow sensitivity, and PKE budget

This stage post-processes the Re = 2000 NACA 0012 direct and adjoint global
stability modes. It follows the cylinder `411_postproc_wavemaker` pattern:
`userParam01 = 4` runs the PKE budget, wavemaker, and base-flow sensitivity
paths from existing mode files.

## Inputs

- Base flow: `BF_naca00120.f00001`
- Direct modes: `dRenaca00120.f0000{1,2}` and `dImnaca00120.f0000{1,2}`
- Adjoint modes: `aRenaca00120.f0000{1,2}` and `aImnaca00120.f0000{1,2}`

The stage-local `maxmodes` is set to 2 because only the two shipped direct
mode snapshots are present.

## Run

```bash
mks naca0012
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
python3 plot.py
python3 scripts/check_against_ref.py example/naca0012_re2000/411_postproc_wavemaker
```

## Outputs

- `wm_naca00120.f00001`: structural-sensitivity wavemaker field
- `sr_*`, `si_*`, `pr_*`, `pi_*`, `tr_*`, `ti_*`: base-flow sensitivity fields
- `K01*`, `K02*`, `PKE_dRe*`: PKE budget fields and integrals
- `plot_wavemaker.png`, `plot_direct_mode.png`, `plot_adjoint_mode.png`
- `wm_metrics.dat`: scalar reference inputs used by `ref/reference.json`
