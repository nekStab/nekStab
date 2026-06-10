# Changelog

All notable changes to nekStab are documented here.

## [2.0.0-rc3] - 2026-06-10

### Added
- General Boussinesq buoyancy API (`nekStab_set_buoyancy` + `nekStab_qvol`
  hooks) carrying the adjoint transpose terms, so thermal adjoint stability
  works in any buoyant case from the `.usr` file alone.
- Reproducibility framework: per-stage `ref/` reference results
  (`reference.json` + figures) with `scripts/check_against_ref.py`, plus
  parsers, a Slurm submit/check driver, and field-inspection utilities.
- Validation gallery deployed to GitHub Pages, catalog-driven and
  Playwright-tested.
- FST (free-stream turbulence) inflow consolidated into a reusable `src/`
  module with a self-describing binary format.
- Dynamic Mode Tracking (DMT) baseflow method.
- Numbered stage layout (`000_dns` through `6xx_modal`) tracked for all
  example families, including the 3D cubic cavity, NACA0012, lid-driven
  cavity, backward-facing step, thermosyphon, flip-flop, tpjet, slot FST,
  and Poiseuille cases.
- Validated flagship stages with committed reference data, among them the
  cylinder Re=100 family, flip-flop Re=62 Newton UPO with direct and adjoint
  Floquet analysis, and thermosyphon Ra=500 baseflow with direct and adjoint
  stability.

### Fixed
- Adjoint Floquet analysis replays the stored base-flow orbit time-reversed;
  adjoint Floquet spectra now match their direct counterparts.
- Thermal adjoint for Boussinesq flows: the buoyancy transpose enters through
  the scalar source and the base temperature gradients feed the adjoint
  momentum production; adjoint eigenvalues now match the direct ones.
- BoostConv QR re-orthogonalization restored to modified Gram-Schmidt.
- POD/DMD/SPOD and mean-field outposts write the scalar field again.
- Drag computed directly on wall (`W`) faces with the correct scaling.
- Floquet mode animation writes the correct deformed output.
- Builds survive configuration changes in non-interactive shells (Slurm, CI);
  missing include dependency edges added so stale objects are rebuilt after
  `NEKSTAB.inc`/`SIZE` changes.
- Orbit allocation failures report the requested size and exit cleanly;
  time-step and energy-budget divisions guarded; per-step solver logging
  throttled (orbit-replay logfiles shrink by two orders of magnitude).

### Changed
- Solver dispatch gated by integer mode-code constants with a string selector
  (`nekstab_mode`) in addition to the numeric `userParam01` codes.
- Legacy flat-layout example directories removed in favor of the numbered
  stages.

## [2.0.0-rc2] - 2026-05-19

Per-case verification campaign on a local Slurm queue (Newton, SFD,
BoostConv, Floquet, and thermal cases re-verified at their reference
parameters), cylinder examples unified on a 2128-element mesh, and build
fixes for spurious full rebuilds.

## [2.0.0-alpha.1] - 2026-04-28

First 2.0 preview: sources reorganized into free-form `.f90` modules,
Krylov routines batched through BLAS, Krylov-Schur improvements (locking,
adaptive restart), and the example tree reworked around per-case directories.

## [1.0.1] - 2021-03-20

Maintenance release of the original toolbox.

## [1.0] - 2021-01-22

First public release: steady-state computation (SFD, BoostConv, Newton),
direct/adjoint eigenvalue problems, transient growth, and post-processing
for Nek5000.
