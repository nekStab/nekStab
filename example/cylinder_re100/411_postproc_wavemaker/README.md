# Cylinder Sensitivity, Budget & Wavemaker — Post-processing

## Physics
2D flow around a circular cylinder at Re = 100. Computes structural sensitivity
(the wavemaker field), perturbation kinetic energy budget terms (production,
dissipation, transport), and sensitivity-to-base-flow modifications, all from
pre-computed direct and adjoint eigenmodes.

## nekStab Mode
`userParam01 = 4` — Full post-processing (includes 4.1 + 4.2 + 4.3)

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (shipped; the base flow of `../310_stability_direct/direct`)
- Direct modes: `dRe1cyl0.f0000{1,2}`, `dIm1cyl0.f0000{1,2}` (copy from `../310_stability_direct/direct/`)
- Adjoint modes: `aRe1cyl0.f0000{1,2}`, `aIm1cyl0.f0000{1,2}` (copy from `../320_stability_adjoint/`)
- The code reads `dRe`, `dIm` from file `f00001` and `aRe`, `aIm` from file `f00002`. Only the `f00002` files are
  tracked; run the two stability stages first and copy all eight mode files.
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
Gate run 2026-10-09 on 8 ranks (Jean Zay, ifort 2021.9), with the modes of the gate runs of the two
stability stages:

- Run completes in 1.4 s (pure field arithmetic on pre-computed modes)
- All 7 output fields written
- The wavemaker maximum is 0.1709 at x = 3.84, y = -0.64. The field exceeds half of this maximum for
  1.7 < x < 5.7 (|y| < 1.0), so it lies in the near-wake recirculation region, where direct mode
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
Giannetti & Luchini (2007), J. Fluid Mech. 581 place the wavemaker inside the recirculation bubble. Here the
maximum of the wavemaker field lies at x = 3.84 diameters from the cylinder centre, inside the bubble of the
steady Re=100 base flow.

`ref/reference.json` holds regression values of this repository: maximum and norm of the wavemaker field
`wm_1cyl0.f00001` (stored in the temperature slot) and the streamwise position of its maximum. Check a run with
`uv run --with numpy --with pymech python scripts/check_against_ref.py example/cylinder_re100/411_postproc_wavemaker`.
