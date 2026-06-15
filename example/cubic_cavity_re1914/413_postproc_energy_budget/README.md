# cubic_cavity_re1914 / 413_postproc_energy_budget

**Stage**: 413 (postprocessing — perturbation kinetic energy budget)
**uparam01**: 4.1 (`nekstab_mode = energy_budget`, steady base flow)
**Re**: 1914 (steady, below the Hopf at Re≈1916.6)
**Inputs**: `BF_cav0.f00001` (converged base flow) + `dRe`/`dIm` eigenmodes
(leading 2), copied from `310_stability_direct/direct`.
**Ranks**: 16

## What it computes

The perturbation kinetic energy (PKE) budget for the leading steady eigenmodes.
For a complex eigenmode `q = q_r + i q_i`, the Reynolds–Orr balance closes as

```
σ_r = ( P − D ) / ( 2 E )
P = − ∫ (q_r,i q_r,j + q_i,i q_i,j) ∂U_i/∂x_j dV   (production)
D = (1/Re) ∫ (∇q_r:∇q_r + ∇q_i:∇q_i) dV             (dissipation)
E = ½ ∫ (|q_r|² + |q_i|²) dV                          (energy)
```

The run prints `σ_r(budget)` vs the eigenvalue and writes spatial budget fields
(`K01…`) plus integrals (`PKE_dRe…`). Closure (`σ_r(budget) ≈ σ_r(eig)`) verifies
the budget. This is the **steady** counterpart to the periodic-orbit
(complex-Floquet) budget — full derivation in `docs/pke-periodic-orbit.md`.

Run with `sbatch run.local.slurm`.
