# Cylinder Sensitivity, Budget & Wavemaker — Post-processing

## Physics
2D flow around a circular cylinder at Re = 100. Computes structural sensitivity
(the wavemaker field), perturbation kinetic energy budget terms (production,
dissipation, transport), and sensitivity-to-base-flow modifications, all from
pre-computed direct and adjoint eigenmodes.

## nekStab Mode
`userParam01 = 4` — Full post-processing (includes 4.1 + 4.2 + 4.3)

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (from `../../baseflow/newton/`)
- Direct modes: `dRe1cyl0.f0000{1,2}`, `dIm1cyl0.f0000{1,2}` (from `../../stability/direct/`)
- Adjoint modes: `aRe1cyl0.f0000{1,2}`, `aIm1cyl0.f0000{1,2}` (from `../../stability/adjoint/`)
- No restart-from-time-history needed (the case reads the eigenmodes directly).

## Mesh
2128 elements (`lelg=2128` in `SIZE`). Same mesh as the rest of the cylinder
isothermal suite.

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Outputs (regenerated per run)
- `wm_1cyl0.f00001` — wavemaker field (structural sensitivity)
- `sr_1cyl0.f00001`, `si_1cyl0.f00001` — sensitivity-to-base-flow real/imag parts
- `pr_1cyl0.f00001`, `pi_1cyl0.f00001` — production real/imag parts
- `tr_1cyl0.f00001`, `ti_1cyl0.f00001` — transport real/imag parts

## Latest Result
Verified 2026-05-15 on 8 ranks (Slurm):

- Run completes in 0.7 s (pure field arithmetic on pre-computed modes)
- All 7 output fields written
- Plot shows the canonical Giannetti-Luchini structure: wavemaker concentrated
  in the near-wake recirculation region (≈ 1 < x < 8), where direct mode
  amplitude and adjoint mode amplitude overlap most strongly.

## Plot
```bash
python3 plot.py
```

Three-panel figure:
- (a) Wavemaker field |dRe · aRe + dIm · aIm| (structural sensitivity magnitude)
- (b) Leading direct eigenmode (real part)
- (c) Leading adjoint eigenmode (real part)

## Reference
Giannetti & Luchini (2007), J. Fluid Mech. 581.
