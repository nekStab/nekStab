# cylinder_re1m/310_stability_direct

**Stage**: stability direct
**Variants**: coupled, quasilaminar

This stage has multiple variants that share the same dispatch subroutine but
use different parameter choices. Each variant lives in its own subdirectory.

The two variants are the two RANS stability operators, run on the same
imperfect SFD field (`base_converged.f00001`). Both spectra are base-flow
spectra of that field, not mean-flow spectra. Coupled versus quasilaminar is
a choice about the closure linearization, not a Reynolds stress:

- `coupled/` — full RANS Frechet. Explicit-filter solve (`filterWeight = 1`,
  `filterCutoffRatio = 0.67`): leading mode σ = 3.41 at St = 0.491, on the
  cylinder. A converged mode sits at St = 0.193 with σ = 3.35.
- `quasilaminar/` — hydrodynamic operator. Same explicit-filter solve:
  leading mode σ = 1.15 at St = 0, a stationary ring on the cylinder.
  A wake-frequency mode remains at St = 0.189, σ = 0.261, but it is not
  leading and was not among the ten written fields.

The hpfrt spectra (quasilaminar σ = 0.237 at St = 0.194) are in
`spectre_hpfrt.tar` in each directory.

Cross-variant residual decay plot: `scripts/nstab_compare.py example/cylinder_re1m/310_stability_direct/`
