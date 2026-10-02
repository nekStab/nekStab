# cylinder_re1m / 310_stability_direct/coupled

**Stage**: direct stability via finite-difference Frechet — the **coupled**
RANS operator (`iffindiff=.true.` in .usr): the perturbation carries velocity
*and* the turbulence scalars (k', tau'), so the eddy viscosity responds to the
perturbation — the full RANS Jacobian (Crouch, Garbaruk & Magidov 2007; Sarras
et al. 2024).
**Re**: 1e6 (k-tau RANS, user-directed).

## Two RANS stability operators

- `coupled` (this directory, `iffindiff` alone): full RANS Frechet — velocity,
  k and tau are all perturbed and the eddy viscosity responds.
- `../quasilaminar` (`iffindiff` + `ifquasilaminar`): the hydrodynamic
  operator — velocity-only perturbation with the eddy viscosity prescribed
  from the base flow (the frozen-eddy-viscosity / quasi-laminar linearisation;
  Mettot & Sipp 2014; Pickering et al. 2021).

This run does not show that agreement. The coupled leading eigenfunctions,
including the St = 0.216 eigenvalue, sit on the cylinder. The quasilaminar
leading velocity field is the wake.

## Status

Run on the SFD field `base_converged.f00001`. Explicit-filter solve
(`filterWeight = 1`, `filterCutoffRatio = 0.67`): leading mode σ = 3.41
at St = 0.491, on the cylinder. A converged mode is at St = 0.193 with
σ = 3.35. The hpfrt spectrum is in `spectre_hpfrt.tar`.
