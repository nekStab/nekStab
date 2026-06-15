# cubic_cavity_re1950 / 413_postproc_energy_budget_floquet

**Stage**: 413 (postprocessing — periodic-orbit PKE budget, complex Floquet)
**uparam01**: 4.11 (`nekstab_mode = energy_budget_floquet`)
**Re**: 1950 (above the Hopf — admits an unstable periodic orbit)
**Ranks**: 16

## What it computes

The orbit-averaged perturbation kinetic energy budget for a **complex** Floquet
mode `q = e^{μt} ũ`, `μ = σ_r + iω`, of the periodic-orbit (UPO) base flow. The
Hermitian product removes the `e^{iωt}` phase, so the budget closes with the same
form as the steady case but averaged over one orbit period `T`:

```
σ_r = ( P̄ − D̄ ) / ( 2 Ē )            (orbit averages ⟨·⟩ = T⁻¹∫₀ᵀ)
```

Writing the complex evolution as `q = a + i b` (`a` = evolution of the mode real
part `dRe`, `b` = of the imaginary part `dIm`), the Hermitian terms are additive:
`a_i a_j + b_i b_j`, `∇a:∇a + ∇b:∇b`, `|a|² + |b|²`, with growth weight
`exp(−2σ_r t) dt/T` (σ_r only — ω cancels). The driver integrates two passes over
the stored base orbit (one per mode part). Full derivation:
`docs/pke-periodic-orbit.md`.

## Inputs (gated — wire before running)

This stage needs, copied/linked into this directory:
- `BF_cav0.f00001` — the **periodic-orbit** base flow IC from
  `../210_baseflow_newton/<variant>/` (Newton UPO).
- `dRecav0.f0000N`, `dImcav0.f0000N` — the leading **complex** Floquet eigenmode
  (real + imaginary parts) from `../311_stability_direct_floquet/`.

Confirm the multiplier is complex (`ω ≠ 0` in the eigenvalue file); if it is real,
the dIm pass is skipped and this reduces to the real-multiplier budget. Then
`mks cav --fresh` and `sbatch run.local.slurm`.
