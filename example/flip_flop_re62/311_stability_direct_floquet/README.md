# Flip-Flop Direct Floquet — Re=62

## Physics
2D flow around two side-by-side cylinders at Re = 62. Floquet analysis of the
periodic flip-flop orbit (T = 8.815172) to detect secondary instabilities. The
Neimark-Sacker (flip-flop) bifurcation sits at Re_c,2 ≈ 61.17, so at Re = 62 the
leading multiplier is just above the imaginary axis.

## nekStab Mode
- `userParam01 = 3.11` — direct Floquet analysis
- `userParam07 = 100` — Krylov subspace dimension
- `endTime = 8.815172` — orbit period. The run takes the period from the time stamp of `BF_2cyl0.f00001`, which holds the same number,
  and replaces `endTime` with it. The first version of this stage carried the starting value of the Newton iteration, 8.73356.
- `maxmodes = 10` in `2cyl.usr` — the run writes ten converged modes
- `viscosity = -62.0` — Re = 62

## Initial condition
`startFrom = BF_2cyl0.f00001` — an orbit of `../210_baseflow_newton`, shipped here
(T = 8.815172; the committed run of 210 gives T = 8.815177).

## Run
```bash
mks 2cyl
sbatch --ntasks=16 --mem-per-cpu=2000 run.local.slurm
```

## Result (2026-10-09, Jean Zay, `ifort`, 16 ranks)
- `Spectre_NSd_conv.dat` has ten rows, the ten leading converged multipliers (`maxmodes = 10`).
- Leading multiplier: sigma = 0.007896230, omega = +/-0.1407652 (just above zero - consistent with
  Re = 62 barely above the secondary bifurcation). The run of 2026-06-11 gave sigma = 0.007896231.
- Rows 3 and 4 are real multipliers at sigma = 3.3e-6 and 6.5e-11, that is |mu| = 1 to within 3e-5.
- Wall time = 1131 s

## Reference (ref/)
- `reference.json` — leading Floquet multiplier (σ, ω) from `Spectre_NSd.dat`
  row 1. Verify a fresh run with
  `scripts/check_against_ref.py example/flip_flop_re62/311_stability_direct_floquet`.
- `ref/mode.png` — leading Floquet mode
- `ref/spectrum.png` — Floquet spectrum
- `ref/upo_snapshot.png` — UPO base-flow snapshot

## Plot notes
- `plot.py` plots the **vx** component with the PiYG diverging colormap: ±1.5 for
  the UPO base flow, ±0.5 for the Floquet mode. Bounds x ∈ [−5, 30], y ∈ [−8, 8];
  cylinders at (0, ±0.85), radius 0.5.
- The spectrum shows the multipliers of `Spectre_Hd.dat` with a residual below 1e-5:
  44 of the 100 Ritz values. `Spectre_NSd_conv.dat` lists the ten that the run writes as modes.
- A full paper figure needs a multi-Re sweep (Re = 60, 62, 63, 67); reference
  data Carini et al. (2014).

## Reference
Carini et al. (2015), JFM 778.
