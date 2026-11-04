# NACA 0012 wavemaker, base-flow sensitivity, and PKE budget

This stage post-processes the Re = 2000 NACA 0012 direct and adjoint global
stability modes. It follows the cylinder `411_postproc_wavemaker` pattern:
`userParam01 = 4` runs the PKE budget, wavemaker, and base-flow sensitivity
paths from existing mode files.

## Inputs

- Base flow: `BF_naca00120.f00001`
- Direct modes: `dRenaca00120.f00001` and `dImnaca00120.f00001` (the wavemaker and the base-flow sensitivity read these)
- Adjoint modes: `aRenaca00120.f00002` and `aImnaca00120.f00002`

The code reads the direct modes from `f00001` and the adjoint modes from `f00002`. The stage ships exactly these five
files, so the wavemaker and the base-flow sensitivity run without another stage. The energy budget loops over
`maxmodes = 2` (in `naca0012.usr`) and reads `dRenaca00120.f00001` and `f00002` and the two `dIm` files; the stage
does not ship the `f00002` direct files (`../310_stability_direct/` ships only `f00001`), so a run stops with
`Error reading restart header` for `dRenaca00120.f00002` unless you copy `dRenaca00120.f00002` and
`dImnaca00120.f00002` from a run of `../310_stability_direct/`. Stage 320 ships the adjoint modes `f00001`;
a run of `../320_stability_adjoint/` writes the `f00002` files too.

## Run

```bash
mks naca0012
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
python3 plot.py
python3 scripts/check_against_ref.py example/naca0012_re2000/411_postproc_wavemaker
```

## Latest Result
- 2026-06-11, 8 ranks, refreshed 310/320 modes with sponge strength 1.7:
  `wm_max = 57.89675140` at `(x, y) = (0.90037745, -0.07858546)`,
  `wm_mean = 0.98948703`.
- Sponge-fix provenance: recomputed after the sponge-strength correction.

## Outputs

- `wm_naca00120.f00001`: structural-sensitivity wavemaker field
- `sr_*`, `si_*`, `pr_*`, `pi_*`, `tr_*`, `ti_*`: base-flow sensitivity fields
- `K01*`, `K02*`, `PKE_dRe*`: PKE budget fields and integrals
- `plot_wavemaker.png`, `plot_direct_mode.png`, `plot_adjoint_mode.png`
- `wm_metrics.dat`: scalar reference inputs used by `ref/reference.json`
