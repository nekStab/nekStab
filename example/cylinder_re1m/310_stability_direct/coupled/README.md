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

Run both from the *same* base flow to compare them: the leading eigenvalue
typically agrees to a few percent, while the structural sensitivity /
wavemaker and the predicted onset differ.

## Status

Run on the SFD steady RANS field `base_converged.f00001`. Documented
**negative result**: the coupled operator's leading cluster sits at
sigma ~ 3.9 and its near-wall k-tau continuum floors at sigma ~ 2.27 even at
k_dim=200, burying the physical wake mode — see `../coupled_continuum.py`.
The physical leading mode is recovered by the quasilaminar operator in
`../quasilaminar` (sigma = 0.395 at St = 0.206, vs observed St = 0.192).
