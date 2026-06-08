# Poiseuille OTD — optimally time-dependent modes at Re=5000

## Physics
Optimally time-dependent (OTD) decomposition of 2D plane Poiseuille flow at
Re = 5000. The OTD method evolves an orthonormal basis that continuously tracks
the most unstable instantaneous directions of the linearized Navier-Stokes
operator, yielding finite-time Lyapunov exponents (FTLE).

## nekStab Mode
- OTD mode (`nekstab_mode = otd`, code 500)
- `viscosity = -5000.0` — Re = 5000

## Initial condition
`poiseuille_OTD0.f00001` — base trajectory IC (shipped). The OTD basis fields
`bf0poiseuille_OTD0.f0000{1,2}` and modes `ip{0,1,2}poiseuille_OTD0.f00001` are
the shipped results.

## Run
```bash
mks poiseuille_OTD
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — converged FTLE of the two leading OTD modes
  (`otd_ftle.dat` last row: mode 1 ≈ -1.511e-3, mode 2 ≈ -1.772e-3, both stable).
  Verify with `scripts/check_against_ref.py example/poiseuille_re5k/500_otd`.
- `ref/lyapunov_exponents.png`, `ref/otd_mode1.png`, `ref/otd_mode2.png`,
  `ref/otd_residuals.png`.
