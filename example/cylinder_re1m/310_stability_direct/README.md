# cylinder_re1m/310_stability_direct

**Stage**: stability direct
**Variants**: coupled, quasilaminar

This stage has multiple variants that share the same dispatch subroutine but
use different parameter choices. Each variant lives in its own subdirectory.

The two variants are the two RANS stability operators, run on the same SFD
steady RANS base (`base_converged.f00001`):

- `coupled/` — full RANS Frechet (`iffindiff` alone). Leading cluster
  sigma = 3.920 at St = 0, floor 2.27 at k_dim=200. The eigenfunctions
  are on the cylinder: all of v'^2 is inside r < 1.5.
- `quasilaminar/` — hydrodynamic operator (`iffindiff` + `ifquasilaminar`).
  Leading mode sigma = 0.395 at St = 0.206, vs observed St = 0.192.
  84% of v'^2 is downstream of x = 2 (centroid x = 5.78). The tau in the
  eigenmode file is the base flow, not the perturbation.

Cross-variant residual decay plot: `scripts/nstab_compare.py example/cylinder_re1m/310_stability_direct/`
