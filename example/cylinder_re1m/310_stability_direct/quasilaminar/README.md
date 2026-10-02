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

This run does not show that agreement. The coupled leading eigenfunctions
sit on the cylinder. The quasilaminar leading velocity field is the wake.

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

Explicit-filter solve (`filterWeight = 1`, `filterCutoffRatio = 0.67`)
on `base_converged.f00001`. Leading mode σ = 1.15 at St = 0, a stationary
ring on the cylinder. A wake-frequency mode remains at St = 0.189,
σ = 0.261, but it is not leading and was not written out. The hpfrt
result, σ = 0.237 at St = 0.194, is in `spectre_hpfrt.tar`.
