# Flip-Flop Direct Floquet — Re=62

## Physics
2D flow around two side-by-side cylinders at Re = 62. Floquet analysis of the
periodic flip-flop orbit (T = 8.73356) to detect secondary instabilities. The
Neimark-Sacker (flip-flop) bifurcation sits at Re_c,2 ≈ 61.17, so at Re = 62 the
leading multiplier is just above the imaginary axis.

## nekStab Mode
- `userParam01 = 3.11` — direct Floquet analysis
- `userParam07 = 50` — Krylov subspace dimension
- `endTime = 8.73356` — orbit period (from the Newton UPO)
- `viscosity = -62.0` — Re = 62

## Initial condition
`startFrom = BF_2cyl0.f00001` — the converged UPO base flow from
`../210_baseflow_newton/upo` (shipped here).

## Run
```bash
mks 2cyl
sbatch --ntasks=16 --mem-per-cpu=2000 run.local.slurm
```

## Latest result (verified 2026-06-11, 16 ranks)
- 3 multipliers converged with residual ≤ 1e-7
- Leading multiplier: sigma = 0.007896231, omega = +/-0.1407652 (just above zero - consistent with
  Re = 62 barely above the secondary bifurcation)
- Wall time = 652.93 s

## Reference (ref/)
- `reference.json` — leading Floquet multiplier (σ, ω) from `Spectre_NSd.dat`
  row 1. Verify a fresh run with
  `scripts/check_against_ref.py example/flip_flop_re62/311_stability_direct_floquet`.
- `ref/mode.png` — leading Floquet mode
- `ref/spectrum.png` — Floquet spectrum
- `ref/upo_snapshot.png` — UPO base-flow snapshot

## Plot notes
- `plot.py` shows velocity magnitude; paper plots the **vx** component (asymmetry
  visible). PiYG diverging colormap: ±1.5 clim for the UPO base flow, ±0.5 for
  the Floquet mode. Bounds x ∈ [−3, 20], y ∈ [−3, 3]; cylinders at (0, ±0.85),
  radius 0.5.
- A full paper figure needs a multi-Re sweep (Re = 60, 62, 63, 67); reference
  data Carini et al. (2014).

## Reference
Carini et al. (2015), JFM 778.
