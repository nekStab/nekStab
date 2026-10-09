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
  sponge also forces the base flow. A Floquet run prints the closure of its
  orbit, `||U(T)-U(0)|| / ||U(0)||` of the base flow integrated over the period
  of the base-flow time stamp, and warns above 1e-3. The orbit of
  `flip_flop_re62` closes to 8e-6; with a period 9 % too long it closes to
  8e-3. A Floquet run also warns when no converged multiplier lies within 1e-3
  of 1, but that test is weaker: with the wrong period the converged list still
  held two multipliers within 1e-11 of 1, and their eigenvectors sat on a few
  mesh vertices instead of on the wake. A Newton or Krylov mode with
  `residualProj = yes` in `[VELOCITY]` prints a note that the velocity
  residual projection can stop the run with a non-finite norm (see Fixed).
- `validation/tests/test_mode_selector.py` compiles `mode_config.f90` against a
  stand-in for the Nek5000 bridge and checks every legacy `userParam01`
  decimal, every mode string, the flag path, the Floquet modifier, the
  conflicts and the invalid inputs.
- OTD about a base flow that the DNS advances: `userParam01 = 5.1`, string
  `'otd_dns'`, flag `ifotdDNS` (mode code 510). The code for the frozen base
  flow stays `5` (`'otd'`). Before, the branch for an evolving base flow in
  `otd.f90` could not be reached, because the mode table had one OTD code.
  With `otd_FTLEPeriod` set, `otd_ftle_period.dat` holds one line per
  completed period, the exponents averaged over exactly that period.
- `ts_gmres` returns the squared residual of its solution in the optional
  argument `final_res`. The steady force sensitivity (modes 4.41 and 4.42) uses
  it to print a warning when GMRES ends above its target of 1e-6; the field is
  still written. Before, the run ended without a message.
- `example/cylinder_re100/130_baseflow_dmt`: the Re = 100 base flow by Dynamic
  Mode Tracking (mode 1.3). With omega_c = 0.7337 (the leading mode of stage
  310), beta = 1 and chi = 0.5 the residual falls below 1e-7 at t = 766, and
  the field agrees with the Newton solution to 4e-5. The README states that DMT
  filters one frequency band, so a steady flow is the only state without force.
- `scripts/compiler_battery.py` exports the tracked files of an example stage,
  builds it under nine compiler profiles (gfortran at -O0, -O2 and -O3, a
  portable build, floating-point traps and signalling NaN; ifx at -O0, the
  default level and a portable build), runs it, and lists every reference
  quantity under each profile with the spread. Seven stages (cylinder DNS,
  thermosyphon Newton and direct, BoostConv, flip-flop wavemaker and Floquet,
  dynamic Newton) pass the run gate and the reference check under every
  profile. The relative spread of the field statistics is 1.7e-4 or less; the
  cylinder DNS agrees to 5e-11 and gives the same Strouhal number to all
  printed digits.
- `validation/tests/test_references_not_vacuous.py` exports `example/` of HEAD
  without running anything and requires that no quantity of any
  `ref/reference.json` passes. A second test requires every quantity to have a
  numeric value, a source and exactly one tolerance.

### Changed
- A `userParam01` that is not a mode code stops the run with an error that
  names the value. Before, `3.15` or `7.3` matched no mode and the run
  fell back to DNS with a warning, and `3.104` was read as `3.1`. A value
  must lie within 1e-4 of a code in the mode table. The float-tolerance
  decoder that followed the integer decoder could not match anything the
  integer decoder had missed, and is removed. Existing `.par` files are
  unaffected: every `userParam01` in the examples decodes.
- The seed loader of the eigensolver and the steady force sensitivity choose
  their files from the mode flags (`isDirect`, `isAdjoint`, `ifForceSensReal`
  and so on) instead of comparing `uparam(01)` to floats. A seed load in a mode
  with no seed file, and a force sensitivity without either flag, now stop the
  run with a message (the file name used to stay undefined).
- Newton residual backtracking is now enabled by default; a case can
  still disable it by setting `ifnewton_backtrack = .false.` in the
  `.usr`. A stalling Newton run now leaves the loop gracefully and
  outposts its best state instead of aborting: when no damped step can
  reduce the residual the current state is held, and six iterations
  without progress end the loop. Three ended a converging run too early:
  `lid_driven_re3600/210_baseflow_newton/fp` from a DNS seed stalled for
  three iterations at iteration 6, wrote no base flow, and converged at
  iteration 15 once the limit was lifted.
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

- `example/cylinder_re100/412_postproc_steady_force_sensitivity/` is removed.
  Its GMRES solve ended at a squared relative residual of 1.35e-2 against a
  target of 1e-6 (T = 1, `k_dim` = 50, 10 restarts), and the stage reported a
  result without saying so. The residual stayed above 7e-3 for T = 1 to 3 and
  for `k_dim` up to 400, and the adjoint map grows to a squared norm of 3.9e25
  at T = 10. The cause is not found. The modes 4.41 and 4.42 stay in the code
  and are documented as not validated in 2.0.0 (`DOC.md`, "Steady force
  sensitivity").

- `example/cylinder_re180/000_dns_seed/` is one DNS at Re = 180 from rest.
  It replaces the Re = 150, 170 and 175 seed folders. The DNS lands on the
  limit cycle, so its last field and its wake-probe period start the Newton
  stage directly. The Newton stage runs without the explicit filter.

### Fixed
- The convergence residual of SFD (`residu.dat`, column 2) and of the TDF mode is
  the root mean square over the domain again. The multi-field change of
  2026-06-21 replaced `normvc` by `k_norm`, the integral norm, which grows with
  the square root of the domain volume: it is 79 times larger on the cylinder
  mesh (volume 6250). The gate, the larger `residualTol` of the `.par`
  (1e-10 in `cylinder_re100/110_baseflow_sfd`), then lay below the noise of the
  solves: the residual stopped falling at 6e-10, SFD never reported
  convergence, and it never wrote `BF_*`. The dynamic solver tolerance, set to
  0.9 times the residual, was 79 times looser than before. BoostConv still
  used `normvc` and was not affected. A cylinder SFD run (`akervik`, 8 ranks)
  converges again at t = 668 with residual 9.998e-11, as before the change.
- The first residual of a steady Newton run was computed with the explicit
  velocity filter and every later residual without it. The steady linearized
  solver runs with `ifbase = .false.`, `q_filter` skips the velocity filter in
  that state, and the first linearized sweep left the flag false for the later
  nonlinear evaluations. For a seed from a filtered DNS the first residual was
  therefore smaller than the later ones (7e-7, then 3e-5 for the same state in
  `lid_driven_re3600`), the first Newton step was wasted, and the stagnation
  count started with a false stall. The nonlinear map of a steady Newton run now
  runs unfiltered from its first evaluation, like the Jacobian. The converged
  base flows are the same: they were already fixed points of the unfiltered map.
- A run that stops on an error ends on all ranks. `nekStab_error` called the
  collective `nek_end`; several callers sit inside `if (nid == 0)` (a missing
  `Spectre_NSd_conv.dat`, a residual file that does not open), so rank 0 waited
  in `nek_end` for ranks that waited in another collective, and the job hung
  until its time limit. `nekStab_error` now prints the message and calls
  `MPI_Abort`: the run ends within seconds with exit code 1.
- `scripts/check_against_ref.py` fails a run whose log holds a line that
  starts with `ERROR:`.
- `naca0012_re2000/210_baseflow_newton` runs at any rank count up to 40. With
  the velocity residual projection on, the Newton run stopped with
  `k_normalize non-finite norm` at 24, 32 and 40 ranks on Jean Zay (ifort
  2021.9) and converged at 8 and 16. The projection is now off in that `.par`
  (`[VELOCITY] residualProj = no`, as in the other Newton and Krylov stages); the
  residuals are identical to seven digits at 8, 16, 24 and 32 ranks, and
  agree with the 16-rank run with the projection on.
- The OTD stages `cylinder_re100/500_otd` and `cubic_cavity_re1914/500_otd`
  set the initial perturbation fields to zero in `useric`, so the OTD basis was
  zero, the FTLE stayed zero and the run reported convergence after 300 steps.
  The `cylinder_re100`, `cylinder_re180` and `cubic_cavity_re1914` stages also
  returned a uniform flow for the base flow in `useric`; `setics` overwrote the
  base flow that the driver had loaded, and the perturbations evolved about
  u = (1, 0) instead of the steady wake (growth rates of -0.014 for all modes).
  All three now read `vxpic`, `vypic`, `vzpic` and return `ubic`, `vbic`, `wbic`,
  and the Gram-Schmidt step stops with a message when a basis vector has zero
  norm. The orthogonality measure printed for one OTD mode was undefined; it is
  zero. `cylinder_re100/500_otd` runs without sponge, which the 310 stage also
  does, and reproduces its growth rate 0.12479 to 1e-4 (seeded with the leading
  eigenmode), or to 1e-3 after about 100 time units from white noise.
- `cylinder_re100/610_modal_dmd` and `620_modal_spod` analysed the wrong snapshot
  window when a restart file at another time was present in the folder, because
  the run script copied the restart of `600_modal_pod` only when none existed,
  and the old reference listed no quantity. The three modal stages now run the
  DNS leg (t = 150 to 200) and the modal analysis in one script from tracked
  files, and their references check the POD energy, the DMD Strouhal numbers and
  the SPOD eigenvalues.
- Streaming SPOD no longer crashes. `spod_stream_update` and three more
  routines outside the module used an implicit real for `spod_s_nfreq` and
  failed at the first snapshot (ifort, segmentation fault). They import
  it and declare `implicit none`. The streaming eigenvalues equal those of the
  batch routine to six digits in every bin.
- The error for a non-finite Krylov norm reports the squared norm, which is
  the quantity that overflowed.
- The Krylov-Schur restart could keep all `k_dim` Ritz values and report
  `k_dim - 2` converged eigenvalues that were not converged (an ellipse with
  `k_dim = 160`: "keeping 160 of 160", then 158 converged). The cap at
  `k_dim - schur_tgt` removed single values, and the step that completes the
  conjugate pairs added them back. With no Arnoldi step after the restart,
  the residual test read the Schur form, where only the last 2x2 block has a
  nonzero last component, so every other Ritz value passed. The cap now removes
  whole pairs (`src/ks_select.f90`), and the restart stops with a message if it
  would keep all `k_dim` values. `validation/tests/test_ks_select.py` fails on
  the old selection (22 of 25 cases) and passes on the new one. A selection
  below the cap is unchanged; the logs at hand keep 21 to 44 of 50 values.
  The pre-releases 2.0.0-alpha.1 to 2.0.0-rc4 have the defect. Version 1.0.1 has no
  cap on the number of kept values at all.
- Thirty-five example stages called `hpts` and aborted at the first step in a
  fresh clone with `Cannot open history file in subroutine hpts()`. `hpts`
  reads the probe list from `<session>.his` and writes the probe values after
  it in the same file, and `.gitignore` excludes `*.his`. Each `run.local.slurm`
  (and `run.sh` of `cylinder_re1m/120_urans`) now writes its probe list before
  the run.
- `hpts` wrote zeros for the last block of probes when the number of probes was
  a multiple of `lhis` (2 probes with `lhis = 2` on `cylinder_re100`; with
  `lhis = 3`, 3 probes with `lhis = 2`, or 1 probe with `lhis = 2` the values
  were correct). The routine is in Nek5000; no shipped stage is affected.
  `nekStab_check_setup` prints a warning when the first line of `<session>.his`
  is a multiple of `lhis`, and `DOC.md` lists the cause under troubleshooting.
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
- The Floquet energy budget (mode 4.11) closes. It had three errors, and the
  closure it prints is the check that finds them. (1) The balance left out the
  damping of the sponge force, so a case with a sponge read the sponge as
  growth: `flip_flop_re62` gave `sigma_r = 8.478e-3` for the eigenvalue
  7.896e-3 (7.4 %). The routine now adds the sponge loss `S_bar`, and the
  closure is 1.2 % (0.99 % at half the time step). `S_bar` is printed and is
  line 15 of `PKE_floquet_*`, the energy-ratio exponent below is line 16; the earlier
  lines keep their position. The flux
  through an open outflow is still not counted: keep the sponge in front of it.
  (2) The period was `endTime` of the `.par` file, while the Floquet eigensolve
  takes it from the time stamp of `BF_*0.f00001`. The `flip_flop_re62` stage
  kept the initial guess of its Newton run (8.73356) and marched 0.9 % short of
  the converged period (8.815177). The routine now takes the period from the time
  stamp and prints a warning when the `.par` value differs. (3) Every mode used
  row 1 of `Spectre_NSd_conv.dat` for its growth rate; mode m now reads row m.
  On the closed `lid_driven_re3600` cavity (a steady base flow, no sponge, two
  modes with different exponents) the closures are 5.1e-4 for the real mode and
  3.5e-4 for the complex mode; mode 2 gave an error of 324 % before. The log also
  prints `sigma_r (energy)`, the exponent from the energy of the marched mode
  after one period, which has no space derivative in it.
- `scripts/inspect_nek_field.py`, the field reader behind the `field:` sources of
  `check_against_ref.py`, read a Nek5000 file wrongly when pymech was not
  installed. It skipped neither the byte-order test float nor the element map,
  read the velocity component by component instead of element by element, and
  expected three coordinates and three velocity components in 2D. A 2D file
  raised "field data ended unexpectedly"; a 3D file gave minima, maxima and
  norms of mixed components without a message (the same `vx`, `vy` and `vz`
  range for the cavity). The plain reader now follows the file layout and gives
  the same minimum, maximum and norm as pymech, to 4e-14, on ten field files of
  the examples (2D and 3D, single and double precision, `XU`, `XUP`, `XUT`,
  `XUPT`, `XT`). Its test wrote a file in the layout the old reader expected; it
  now writes the Nek5000 layout, 2D in single and 3D in double precision.
- The DMT mode (`userParam01 = 1.3`) had three faults. (1) After every step it
  copied the residual components into `t(:,:,:,:,1)` and `t(:,:,:,:,2)`: with
  `ldimt = 1` the second copy lay outside `t`, and with a temperature or scalar
  field both overwrote it. The copies are removed. (2) The converged `BF_*` file
  was written in single precision, which rounds the base flow to 6e-8 before the
  stability stages read it; it is written in double precision now, as SFD,
  BoostConv and Newton write theirs. (3) Nek5000 shortens the last step to land
  on `endTime`, and the BDF2 stencil of DMT assumes one time step: the residual
  of that step (divided by the short `dt`) jumped by a factor 10 to 100 in the
  last row of `residu_dmt.dat`. The step is skipped now, and any other change of
  the time step stops the run with a message. On the Re = 100 cylinder
  (`omega_c = 0.7337`, `beta = 1`, `chi = 0.5`, 8 ranks) DMT reaches 1e-7 at
  t = 766.27, and its base flow agrees with the Newton solution to 4e-5 (minimum
  of `vx`) and 8e-6 (norm), as the SFD base flows do.
- BoostConv (`userParam01 = 1.2`) does not sample at step 0 any more. No step has
  been taken there, so the "residual" compared the seed with an unset previous
  step (1.67 on the cylinder), its rate divided by `dt = 0`, and the first line of
  the log and of `residu.dat` read `Infinity`; `check_against_ref.py` rightly
  fails a run with that word in its log. The first snapshot of the subspace was
  that difference as well. The cylinder stage (`cylinder_re100/120_baseflow_boostconv`,
  Re = 50, the defaults `bst_snp = bst_skp = 10`) meets its gate 1e-6 at t = 159.4, and
  its base flow differs from the Newton solution at Re = 50 by 2e-5 at most in the extrema
  and 2.5e-4 in the norm. The stage ran with `bst_snp = 30` and `bst_skp = 5`, which needs
  t = 168; 30 and 10 does not reach the gate by t = 500. At Re = 100 neither the defaults
  nor 30 and 5 reach it by t = 500 (the stage README has the table).
- `scripts/check_against_ref.py` no longer fails a clean run because of one log
  line. With a fixed time step Nek5000 prints `DTCFL` and `DTFS` unset in the line
  `DT/DTCFL/DTFS/DTINIT`, and some builds show `NaN` there (gfortran at `-O3` did).
  The run gate skips that line; a `NaN` anywhere else still fails the run.

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
