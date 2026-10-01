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

Run on the SFD field `base_converged.f00001`. Leading cluster sigma = 3.920
at St = 0, floor 2.27 at k_dim=200 (`Spectre_NSd.dat`). The eigenfunctions
are not a buried wake roller: in `dRe1cyl0.f00001` and `f00002`, all of
v'^2 is inside r < 1.5 and the centroid is on the body. Velocity is 83% of
the mean square, k 15%, tau 2%. The wake mode is the quasilaminar result
(sigma = 0.395 at St = 0.206, vs observed St = 0.192).
