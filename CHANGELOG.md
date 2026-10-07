# Changelog

Differences between tagged versions. `RELEASE.md` is the complete change
description for the 2.0 series; entries here are the per-tag deltas.

## [Unreleased]

### Added
- Setup checks at start-up. A Floquet run (direct, adjoint or transient
  growth) now stops with a message when the restart file (the orbit) was
  written at another polynomial order than `lx1` in `SIZE`; Nek5000
  interpolates such a file, and the multiplier 1 of the phase mode is
  lost. The steady stability modes print a warning for the same mismatch.
  Newton runs are not checked: they may start from a seed of another
  order. A stability mode with an active sponge prints a note that the
  sponge also forces the base flow. A Floquet run warns when no converged
  multiplier lies within 1e-3 of 1.

### Changed
- Newton residual backtracking is now enabled by default; a case can
  still disable it by setting `ifnewton_backtrack = .false.` in the
  `.usr`. A stalling Newton run now leaves the loop gracefully and
  outposts its best state instead of aborting: when no damped step can
  reduce the residual the current state is held, and three iterations
  without progress end the loop.
- The two RANS stability operators in `example/cylinder_re1m/310_stability_direct/`
  are renamed: `findiff/` becomes `coupled/` (the full RANS Frechet operator —
  velocity, k and tau perturbed, eddy viscosity responding, still selected by
  `iffindiff` alone) and `findiff_frozen/` becomes `quasilaminar/` (the
  hydrodynamic operator — velocity-only perturbation, eddy viscosity
  prescribed from the base flow). The `.usr` flag `iffrozenEV` is renamed to
  `ifquasilaminar` in the same common-block slot, and `src/frozen_mut.f90`
  becomes `src/base_mut.f90` (module `nekstab_base_mut`, procedures
  `base_mut_capture`/`base_mut_get`/`base_mut_ready`/`base_mut_reset`, object
  `base_mut.o`). The "School A"/"School B" names are retired; no aliases or
  compatibility shims are kept.

- `example/cylinder_re180/000_dns_seed/` is one DNS at Re = 180 from rest.
  It replaces the Re = 150, 170 and 175 seed folders. The DNS lands on the
  limit cycle, so its last field and its wake-probe period start the Newton
  stage directly. The Newton stage runs without the explicit filter.

### Fixed
- Transient growth writes the true optimal response. The `ore` output was
  the optimal perturbation multiplied by G, because the response call went
  through `matvec`, which dispatches on the mode flags and so still applied
  the direct-adjoint map. The call now switches to the direct map for that
  one step and restores the flags and the `evop` tag. Check on the back
  step: the mass-weighted energy ratio of response to perturbation equals G
  (7.69523 against 7.695277). The gains were never affected.

- `cylinder_re180` stages run at one mesh order. The Floquet stages were
  compiled at polynomial order 5 and read the order-7 Newton orbit by
  interpolation. The orbit then did not close, and the multiplier of the
  phase mode missed 1 by 6.8e-4 whatever the time step. With the `SIZE` of
  the Newton stage the miss is 1.8e-6. The Floquet stages run with the
  sponge off, because the sponge also forces the base flow.

- `cubic_cavity_re1914` runs with the uniform lid `u = 1`, as in the
  thesis. Since 2026-06-16 the cases used a regularized lid, under which the
  flow at Re = 1914 is more strongly stable and the Hopf pair is absent.
  With the uniform lid the pair is back: sigma = -1.8e-4, omega = 0.5868
  (St = 0.0934); direct and adjoint agree to 2e-4 in sigma.

- Example plot scripts run again with matplotlib 3.9 and later.
  `nekplot.discrete_cmap` looked up colormaps through `plt.cm.get_cmap`,
  which matplotlib removed. Nine scripts also found `nekplot.py` through a
  fixed folder depth that broke when the stage folders were flattened;
  they now search the parent folders.

- `moving_cylinder_re100/000_dns` compiles again. Its `.usr` used the old
  module name `krylov_subspace` instead of `nekstab_krylov_subspace`.

- Example run scripts no longer name the Slurm partition `local`. They
  use the default partition of the queue they are submitted to.

- `nekStab_avg` now averages every scalar slot, not only the
  temperature. The mean and the second moment loop over `i = 1..ldimt`,
  so k and tau of a RANS model are averaged too. The scalar length is
  `nelt`, as in Nek5000 `avg_all`. Before, slots 2 to `ldimt` stayed
  zero.

- Eigenmode outpost no longer writes leftover base scalars into slots
  the Krylov product did not solve. `solver = none` leaves `ifpsco`
  true, so a slot is copied only when `idpss` is not negative.
  Unsolved slots stay zero. A quasilaminar run had been writing the
  base tau into every mode file.

- GMRES no longer reads one past `yvec` when Arnoldi uses every
  column without reaching its tolerance. The column count is stored
  inside the loop and is what `k_matmul` uses. A finished DO index is
  not read: the standard leaves it undefined, and gfortran and ifort
  set it one past the last column. `Q` is one longer than `yvec`, so
  the illegal access was the coefficient, and the restart matvec then
  died in `k_normalize`. An early exit still uses that column.
  `arnoldi_factorization` had the same pattern on `mstep`, which is
  why the residual matvec logged `from 101/100`; `mstep` now stays on
  the column just built.
- Arnoldi keeps the Hessenberg columns already built. `H` was
  `intent(out)`, so a Krylov-Schur restart formally discarded the
  Schur block it had just condensed. The factorization only writes
  the new columns.
- The Newton residual log starts `k_sum` at 0. `k_out` is a saved
  count of the previous GMRES columns and is added before GMRES runs;
  without an initial value the first row was undefined.
- Biorthogonalization passes scalar fields with shape `(lt, ldimt)`.
  The old length-`lv` arrays were short by `ldimt` on every RANS
  `SIZE`, so ifx rejected the build (error #7983) before it could
  link.
- The Newton closing message reports success only when the residual
  is below the target. A finished DO index is not read. A run that
  hits the iteration cap, stagnates, or diverges is not reported as
  a successful finish.
- Krylov normalization rejects a non-finite norm (NaN or Inf) before
  scaling. An infinite norm used to pass the NaN test, become a zero
  scale factor, and write NaN into the next Krylov vector.
- Newton aborts when the forward residual is not finite, at
  `Computed residual:`, instead of falling through into GMRES. The
  comparisons `residual < dtol` and `residual > 1e8 * initial` are both
  false for NaN, so the failure was previously reported one step later
  from `k_normalize`.
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
