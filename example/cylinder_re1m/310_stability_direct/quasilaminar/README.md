# cylinder_re1m / 310_stability_direct/quasilaminar

**Stage**: direct stability via finite-difference Frechet (RANS),
**quasilaminar** — the hydrodynamic operator: velocity-only perturbation with
the eddy viscosity prescribed from the base flow (the frozen-eddy-viscosity /
quasi-laminar linearisation). Sibling of `../coupled` (the full RANS Frechet).
**Re**: 1e6 (k-tau RANS, user-directed).

## Two RANS stability operators

- `../coupled` (`iffindiff` alone): the perturbation carries velocity *and*
  the turbulence scalars (k', tau'), so the eddy viscosity responds to the
  perturbation — the full RANS Jacobian (Crouch, Garbaruk & Magidov 2007;
  Sarras et al. 2024).
- `quasilaminar` (this directory): the perturbation is velocity-only and the
  eddy viscosity is held at its base value — the hydrodynamic-only operator
  (Mettot & Sipp 2014; Meliga et al. 2012; Pickering et al. 2021).

Run both from the *same* base flow to compare them: the leading eigenvalue
typically agrees to a few percent, while the structural sensitivity /
wavemaker and the predicted onset differ.

## How this case selects the quasilaminar operator

- `1cyl.usr`: `ifquasilaminar=.true.` (velocity-only Frechet seed — k', tau'
  are zeroed in the perturbation) and `ifKEnorm=.true.` (velocity
  kinetic-energy norm — the turbulence scalars carry no weight in the inner
  product).
- The base eddy viscosity is snapshotted once, on the first finite-difference
  substep (`nekStab_usrchk` calls `base_mut_capture`, the general routine in
  `src/base_mut.f90`), and replayed in `uservp` via `base_mut_get` for every
  later evaluation — so mu_t never responds to the perturbation. The explicit
  snapshot makes the prescription model-agnostic (correct for a
  strain-dependent SST model too, where `solver=none` alone would still let
  mu_t track the perturbed velocity) and immune to in-place `limit_ktau`
  clipping of the base k, tau.
- `1cyl.par`: `[SCALAR01]`/`[SCALAR02]` `solver = none` keeps k, tau pinned at
  the base flow (belt-and-suspenders with the zeroed seed). For the standard
  k-tau model (`m_id=4`, `mu_t = rho*alp_str*k*tau`, velocity-independent) the
  captured snapshot and `solver=none` coincide — a useful consistency check.

## Result

Run on the SFD steady RANS field `base_converged.f00001`. Leading mode
sigma = 0.395 at St = 0.206 — the von Karman wake mode, against the observed
St = 0.192.
