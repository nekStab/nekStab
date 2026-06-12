# Changelog

Differences between tagged versions. `RELEASE.md` is the complete change
description for the 2.0 series; entries here are the per-tag deltas.

## [Unreleased]

### Fixed
- Newton backtracking with periodic-orbit modes evaluates each trial at
  its damped period: the time horizon (`nsteps`) is re-prepared per
  trial and the orbit storage grows when a trial orbit is longer than
  the allocation. This closes the rc4 known limitation; the flag is now
  safe for UPO runs.

## [2.0.0-rc4] - 2026-06-12

Newton robustness. 15 commits since rc3.

### Added
- Newton residual backtracking globalization (opt-in via
  `ifnewton_backtrack` in the `.usr`, default off): a non-monotone
  Armijo line search wraps the Newton update, rejecting overshooting
  or NaN steps and damping the periodic-orbit period correction with
  the same factor as the state. The non-monotone window (worst of the
  last three accepted residuals) deliberately tolerates the transient
  residual growth seen with rough initial guesses. Validated on the
  backward-facing-step baseflow: a configuration whose full Newton
  step previously grew the residual now converges, with the
  overshooting step rejected at full length and accepted at half.
- Backward-facing-step Re=500 transient-growth stage validated against
  Barkley, Blackburn & Sherwin (2008) at six time horizons.
- Moving-cylinder DNS stage with reference data; NACA0012 DNS ringdown
  and wavemaker stages; thermosyphon and flip-flop wavemaker stages,
  all with reference data.

### Fixed
- Sponge strength is applied to every forcing branch (perturbation and
  temperature, not just the velocity DNS branch); sponge-affected
  stability references regenerated.
- Newton aborts loudly when the first residual is exactly zero (a
  degenerate forward map previously reported instant fake convergence).
- The cached periodic-orbit boundary vector carries the orbit period,
  fixing the period column of the UPO Jacobian.
- Krylov vector normalization fails loudly on a NaN norm instead of
  silently zeroing the vector.
- Cleaner Krylov seed initialization and baseflow reload in the
  eigensolvers.

### Known limitations
- With backtracking enabled for periodic-orbit Newton, trial residuals
  are still evaluated at the undamped period; keep the flag off for
  UPO runs until the follow-up fix lands.

## [2.0.0-rc3] - 2026-06-10

Adjoint correctness and reproducibility. ~95 commits since rc2.

### Added
- General Boussinesq buoyancy API (`nekStab_set_buoyancy` + `nekStab_qvol`)
  carrying the adjoint transpose terms — thermal adjoint stability from the
  `.usr` file alone.
- Dynamic Mode Tracking (DMT) base-flow method (mode 1.3 / `dmt`).
- FST (free-stream turbulence) inflow as a reusable `src/` module with a
  self-describing binary mode file.
- Reproducibility framework: per-stage `ref/` reference results
  (`reference.json` + figures) with `scripts/check_against_ref.py`, plus
  parsers, a Slurm submit/check driver, and field-inspection utilities.
- Validation gallery deployed to GitHub Pages, catalog-driven and
  Playwright-tested.
- Numbered stage layout (`000_dns` through `6x0_modal_*`) tracked for all
  16 example families; validated flagship stages committed with reference
  data (cylinder Re=100 family, flip-flop Re=62 UPO + direct/adjoint
  Floquet, thermosyphon Ra=500 direct/adjoint, and others).

### Fixed
- Adjoint Floquet analysis replays the stored base-flow orbit
  time-reversed; adjoint Floquet spectra now match their direct
  counterparts (flip-flop Re=62: 0.0078864 +/- 0.1407585i vs direct
  0.0078930 +/- 0.1407626i).
- Thermal (Boussinesq) adjoint: buoyancy transpose enters the adjoint
  scalar source and base temperature gradients feed the adjoint momentum
  production; adjoint eigenvalues match the direct ones (thermosyphon
  Ra=500: 0.1022272 vs 0.1022191).
- BoostConv QR re-orthogonalization restored to modified Gram-Schmidt.
- POD/DMD/SPOD and mean-field outposts write the scalar field again.
- Drag computed directly on wall (`W`) faces with the correct scaling.
- Floquet mode animation writes the correct deformed output.
- Builds survive configuration changes in non-interactive shells
  (Slurm, CI); include dependency edges added so stale objects rebuild
  after `NEKSTAB.inc`/`SIZE` changes.
- Orbit allocation failures report the requested size and exit cleanly;
  time-step and energy-budget divisions guarded; per-step solver logging
  throttled (orbit-replay logfiles shrink by two orders of magnitude).

### Changed
- Solver dispatch gated by integer mode-code constants behind the string
  selector (`nekstab_mode`) and the numeric `userParam01` codes.
- Legacy flat-layout example directories removed in favor of the numbered
  stages.

## [2.0.0-rc2] - 2026-05-19

Validation campaign. 675 commits since alpha.1.

### Added
- Per-case verification on a local Slurm queue: 24 cases re-verified at
  their reference parameters (SFD ladder, BoostConv, Newton fixed points
  and UPOs, direct/adjoint stability, Floquet, thermal Newton-GMRES),
  each with base flow and leading-mode evidence committed.
- String-based mode selection (`nekstab_mode = 'direct'` etc.) with
  3-priority resolution: string > flags > `uparam(1)`.

### Changed
- All isothermal cylinder cases unified on a 2128-element mesh with a
  shared seed initial condition for cross-case comparison.

### Fixed
- `makeneks` no longer rewrites `NEKSTAB.inc` spuriously on every build.
- Mode flags cleared before the string-mode override is applied.

## [2.0.0-alpha.1] - 2026-04-28

The 2.0 rewrite. ~660 commits since 1.0.1; see `RELEASE.md` for the full
description.

- All sources restructured into free-form Fortran modules;
  `nekstab_nek_bridge` isolates the Nek5000 common-block imports.
- Multiple passive scalars, RANS (finite-difference Fréchet operator),
  and conjugate heat transfer (`lt != lv`) supported through the whole
  stability/modal pipeline.
- MPI reduction batching: Gram matrices, projections, and cross-spectral
  densities computed via BLAS with a single allreduce; CGS2
  orthogonalization (2 allreduces per Arnoldi step, independent of k).
- Krylov-Schur improvements: heap allocation, eigenvalue locking,
  adaptive restart, conjugate-pair handling, combined abs/rel
  convergence.
- Inexact Newton-GMRES with residual-proportional adaptive tolerance and
  stagnation guard.
- Modal analysis framework: POD, projected DMD, and SPOD (batch and
  streaming) sharing the batch inner-product infrastructure.
- Floquet energy budget (mode 4.11) for time-periodic base flows.

## [1.0.1] - 2021-03-20

- Fixed an extra outpost when using Krylov-Schur.
- CFL limiter threshold raised from 1 to 10.
- Sensitivity subroutine with normalization; adjoint reference modes.
- New example cases (porous, MFM); GitHub Actions compiler check.

## [1.0] - 2021-01-22

First public release: steady-state computation (SFD, BoostConv, Newton),
direct and adjoint eigenvalue problems, transient growth, and
post-processing tools for Nek5000.
