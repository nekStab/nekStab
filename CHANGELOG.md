# Changelog

Differences between tagged versions. `RELEASE.md` is the complete change
description for the 2.0 series; entries here are the per-tag deltas.

## [Unreleased]

## [2.0.0-rc5] - 2026-10-17

### Added
- DMT filter bank, mode 6.4 (`userParam01 = 6.4`, string `'dmt_bank'`, flag
  `ifdmtbank`). Second-order filters of the Dynamic Mode Tracking kind, one per chosen
  frequency, run on the DNS without a force and return the spatial structure of the flow at
  their frequency as a complex field (`bRe`, `bIm`); a frequency 0 gives the mean flow. The
  frequencies are registered in `nekStab_usrchk` with `dmt_bank_add(omega, beta)`, at most 8.
  Each filter is a damped oscillator driven by the flow and integrated exactly for an input
  that varies linearly between the steps (`dmt_filter.f90`). The amplitude comes from the
  filtered signal and its derivative, which passes a constant mean flow with no gain; the
  combination `beta*(q + i*omega*p)` has the gain `i*beta/omega` for it. `dmt_bank_options`
  sets the length of the history `dmt_bank.dat` (norms of the flow, of what no filter holds,
  and of each amplitude), the output of the filtered fields at each output step, and a window
  over which the amplitude is averaged: a filter also passes about `beta/omega` of the
  neighbouring frequencies, and the average over whole periods removes that ripple (a test
  with three harmonics on a mean: 5 %, 15 % and 45 % at the last step, 4e-5, 5e-5 and 1e-4 on
  the average of 10 periods). Mode 1.3 is unchanged.
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
- `example/cylinder_re100/140_baseflow_meanflow` and
  `315_stability_direct_meanflow`: the mean flow of the shedding Re = 100 wake
  and the direct stability about it. The DNS averages velocity and pressure over
  25 whole periods, found from the upward zero crossings of the integral of v over
  the domain, writes the mean as a base-flow file and stops. The mean flow has a
  recirculation length of 1.42 diameters; the steady flow has 6.15. The leading
  eigenvalue about the mean flow is sigma = 0.0022, omega = 1.0290, within 1.2 %
  of the DNS frequency 1.0417; the steady flow gives omega = 0.7337. Marquet &
  Carini (2020, arXiv:2008.10002) give 0.0020 and 1.0322 for the mean flow at
  Re = 100.
- `example/cubic_cavity_re1920`: the 3D cubic cavity at Re = 1920, above the first
  Hopf bifurcation. A DNS at Re = 1900 gives the start field, Newton-GMRES finds
  the unstable steady flow in 12 iterations, and the direct stability stage finds
  one unstable pair, sigma = +5.84e-4, omega = 0.5868. At Re = 1914 the same pair
  has sigma = -1.83e-4; linear interpolation of the two puts the zero at
  Re = 1915.4 on this mesh (1916.63 in the thesis).
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
- `example/cubic_cavity_re1914/330_transient_growth` and `411_postproc_wavemaker`.
  The flow at Re = 1914 is stable (sigma = -1.8e-4), and a perturbation still gains a
  factor 229 in energy over the horizon tau = 10 (G = 229.19, the second gain is 215.96).
  The energy ratio of the stored response to the stored perturbation is 0.3 % below G;
  this measures how far the adjoint integration departs from the transpose of the forward
  one on this mesh. The wavemaker of the leading pair has its maximum, 45.0, at
  (x, y, z) = (-0.365, -0.244, 0) and a second maximum of 44.8 at (-0.307, 0.333, 0), so the
  position of the global maximum is not stable.
- `example/cubic_cavity_re1950`: DNS, POD, DMD and energy-budget stages on the limit cycle of
  the 3D cavity at Re = 1950 (period 10.78909). `000_dns` runs 100 time units on the cycle and
  writes 200 snapshots; `600_modal_pod` and `610_modal_dmd` read 100 of them, 4.6 periods. Two
  POD modes hold 94.3 % of the energy and turn at 1.0005 times 2 pi / T; the next four pairs
  turn at 1.9998, 2.9989, 3.9981 and 4.9975 times. The DMD modes lie on the harmonics 1 to 6 of
  the cycle, within 4.2e-4 in units of 2 pi / T, and none grows (largest modulus 0.999978).
  `411_postproc_energy_budget` runs the Floquet energy budget (mode 4.11) on the closed cavity,
  where P - D = sigma E holds exactly. The energy of the marched mode gives sigma to 0.17 %
  (2.90305e-3 against 2.907859e-3), and production and dissipation balance to 0.34 % of the
  dissipation. P/E and D/E are 175 times larger than sigma, so the exponent from the budget, the
  difference of the two, is 58.5 % high. The README says so, and the reference does not check
  that exponent.
- `example/cubic_cavity_re1950/630_modal_dmt` and `example/moving_cylinder_re100/630_modal_dmt`
  run the DMT filter bank (mode 6.4) with four filters: the frequency of the cycle or of the
  forcing, its second and third harmonic, and the mean flow. On the cavity cycle (15 periods,
  40 ranks, 1292 s) the amplitudes rebuild the last field to 0.227, 0.059 and 0.019 of its
  fluctuation after one, two and three harmonics; the first POD pair of stage 600 leaves
  0.239 after the fundamental. On the forced orbit of the oscillating cylinder (30 periods,
  40 ranks, 643 s) the amplitudes rebuild the orbit of stage 210 to 0.401, 0.236 and 0.145 in
  the near field. On the cylinder wall the amplitude of the first filter is 1.1148 against the
  wall speed A 2 pi f = 1.1149, and the other filters return less than 2e-6.
- `example/moving_cylinder_re100` (010 to 311): the cylinder in forced transverse oscillation at
  Re = 100, A = 1.08, f = 0.1643 on a moving mesh. The DNS locks to the forcing; Newton-GMRES
  finds the orbit in 11 iterations (squared residual 3.7e-3 to 9.7e-10, 3972 s on 40 ranks);
  the direct Floquet stage finds one converged multiplier, mu = 0.92988, real and inside the unit
  circle, so the lock-in is stable, and the DNS relaxes onto the orbit by 0.9252 per period,
  0.5 % below |mu|.
- `example/tpjet_re2000` (140 to 321): the harmonically forced axisymmetric jet at Re = 2000
  (St = 0.6, T = 1.66667, 4800 elements, `lx1 = 6`, 16 ranks). Time-delayed feedback (mode 1.4,
  stage 140, 4880 s) starts from the Re = 1900 orbit; its residual rises to 3.2e-4 at t = 21.5 and
  falls as t^-2.49 to 6.5e-8 at t = 600, and the field lies at 1.9e-3 of the Newton orbit (the seed
  at 1.05e-2). Stage 210 finds the orbit in 5 iterates (squared residual 3.5e-6 to 9.7e-9, 66
  Arnoldi vectors, 2663 s). The direct Floquet stage 311 finds one real multiplier outside the unit
  circle, mu = -1.1577 (a period doubling, sigma = 0.0878; the run stops at step 57 of 128 with
  13 Ritz values converged, orbit closure 4.6e-5, 2408 s). The other converged values lie on an
  arc inside the unit circle that moves with the step at which the run stops (-0.9241 +- 0.0458i
  after 57 steps, -0.8599 after all 128). The adjoint stage 321 gives mu = -1.1573 (step 56 of 128,
  22 converged, 1890 s) and an adjoint mode at the nozzle lip, z = 0.64 and r = 0.48, where the
  direct mode peaks at z = 20.6. The damped multipliers of the adjoint map are not those of the
  direct map (the second is -0.9063 after 56 steps, -0.9050 after 128 and -0.9654 at
  `k_dim` = 40), so the references check the leading multiplier only.

### Changed
- `example/tpjet_re2005` is `example/tpjet_re2000`. The stages ran the 4800-element mesh of the
  AMR_Krylov paper at `lx1 = 5` from a seed stored at `lx1 = 6`, and at `lx1 = 5` the Floquet
  spectrum at Re = 1900 had a complex pair of modulus 2.35; at `lx1 = 6` the leading multiplier is
  real, -1.139, a period doubling as in the paper. Every stage uses `lx1 = 6` and Re = 2000 now,
  and the references, the shipped orbit and the shipped modes of the old stages are replaced by
  those of the runs above. The DMT stage `130_baseflow_dmt` is removed: the band-pass filter has
  zero gain at zero frequency only, so it holds a steady state and not a periodic orbit.
- The stage folders of `example/` no longer track what a run writes. 150 copies of figures at
  the root of the folders (`plot_bf.png`, `plot_spectrum.png` and the like; the plot scripts
  write them at each run and `ref/` holds the maintained copy that the gallery draws), four
  `plot.out` logs (they held the paths of the machine that ran the plot script) and the twelve
  output fields `pr_`, `pi_`, `sr_`, `si_`, `tr_` and `ti_` of the wavemaker stages of
  `cylinder_re100` and `naca0012_re2000` (23 MB) are removed. The folders of
  `cylinder_re1m` and `plot_compare.png` of `110_baseflow_sfd`, a figure over six variants,
  stay.
- `bin/neksubmit` reads the remote host, the work directory and the Slurm account from
  `NEKSUBMIT_HOST`, `NEKSUBMIT_WORKDIR` and `NEKSUBMIT_ACCOUNT`, and stops with a message when one
  is missing. Before, a login, a work path and an account were written in the script. The queue
  checks of `bin/nekStab_tools.py` use `$USER`.
- `README.md` has a Disclaimer section: no warranty, and research software that is not validated
  for engineering design, certification or safety-critical use.
- `validation/validate.py` (with `validation/tests/test_validate.py`) and
  `example/collect_plots.py` are removed. Both listed example directories of an
  earlier layout that no longer exist. A stage is checked with
  `scripts/check_against_ref.py` against its `ref/reference.json`; the page of
  references and figures is built by `validation/build_gallery.py`, and
  `validation/tests/test_gallery_complete.py` fails for a stage without a
  reference, figure or README. The Examples and Validation sections of `DOC.md`
  describe the current layout.
- The validation gallery is built from the stage folders of `example/`
  (`python validation/build_gallery.py`). A card shows the figures of the
  `ref/` folder of a stage and the quantities of its `ref/reference.json`,
  with the new keys `summary` (one line) and `description` (what the values
  check). A stage without a reference, a summary, a figure or a README is not
  drawn, and `validation/tests/test_gallery_complete.py` fails for it. The
  test also fails for a deferred, blocked or pending card, a link that does not
  resolve, and a `validation/index.html` that differs from the generator
  output. The page no longer carries a date or a machine name, so the file in
  git is the output of the script. Before, the cards came from a hand-written
  case catalog. Almost every entry of it named a figure file that no longer
  existed, so most cards read "pending re-run". The catalog
  (`validation/catalog.py`, `catalog_check.py`, `scripts/evidence_report.py`) and its
  tests are removed. `validation/serve.py` allows the stage folders for the
  recompile and resubmit buttons; before, its list named folders that were
  gone.
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

- The Newton stages `example/naca0012_re2000/210_baseflow_newton/fp` and
  `example/tpjet_re2000/210_baseflow_newton` start from a field that is not
  converged. Before, each started from its own converged output: the first
  residual was already below the tolerance and the stage did not test Newton.
  The airfoil stage starts from the last DNS snapshot (`BF_seed_naca00120.f00001`,
  t = 70) and needs 5 iterations (squared residual 6.1e-7 to 9.5e-10, 956 s on 8
  ranks). The jet stage starts from the Re = 1900 field
  (`BF_Re1900_tpjet0.f00001`) and needs 5 iterates, the seed and four updates
  (squared residual 3.5e-6 to 9.7e-9, 66 Arnoldi vectors, 2663 s on 16 ranks).
  Both references check the first residual and the iteration count, so a run that
  starts from a converged field fails. The sponge forces the velocity toward the
  field that the run holds at its start, so the airfoil base flow depends on its
  seed inside the sponge zone (`DOC.md`, "Sponge reference field"). It differs by
  1.1e-2 (relative l2 norm of the velocity) from the base flow that the stability
  stages ship, and the leading growth rate of the direct stage is -0.0643 on it
  instead of -0.0665.

- `example/cylinder_re1m` (2D RANS) is marked experimental in `RELEASE.md`. It has the
  finite-difference operator only (no adjoint, wavemaker or transient growth), no reference
  data and no card in the gallery, and its stages are not part of the rerun of the examples.

- The Krylov-Schur eigensolver stops on the leading eigenvalues, and inside the Arnoldi
  pass. The run stopped when any `schur_tgt` Ritz values had converged, tested after all
  `k_dim` steps. A run could then stop with the leading values unconverged and damped values
  converged below them (cylinder at Re = 42, reported by a user: rows 1 and 2 unconverged,
  the converged shedding pair in rows 3 and 4). And a long basis was always paid in full: in
  the rerun of the examples 20 of 24 eigenvalue runs did exactly `k_dim` matvecs and ended
  with 13 to 83 converged values where 2 to 6 were wanted. The run now tests after every
  Arnoldi step and stops at the first step where the first `schur_tgt` entries of the list
  sorted by modulus have converged (`ks_leading_converged`; a conjugate pair counts as two).
  `k_dim` is the largest basis, not a fixed cost; `schur_tgt <= 0` still runs one full pass.
  `lid_driven_re3600/310_stability_direct` stops at step 114 of 120 with the reference rows
  unchanged. Had only the old count moved inside the pass, the lid run
  would have stopped at step 75 with its leading mode at residual 6e-4. One function,
  `ks_converged` (`residual < eigen_tol*max(abs(lambda), 1)`), now makes the test for the
  stop, the flag column, the export and the Arnoldi checkpoint, and rank 0 makes the stop
  decision for all ranks.
- The stop also waits until the leading Ritz values stand still. A residual below `eigen_tol`
  bounds the error of the eigenvalue only through its condition number, which is large for a
  non-normal operator. The adjoint Floquet run of `tpjet_re2000/321` stopped at step 37 of 128
  with the residual 1.1e-6 and the multiplier -1.16222, against -1.15792 of the direct stage
  and -1.15726 of the run before the stop moved into the pass (whole basis of 128); the
  condition number of that multiplier is 4e5. A leading value now has to move by less than
  `eigen_tol*max(abs(mu), 1)` in each of the last 5 steps (`ks_stationary`,
  `ks_leading_settled`); a conjugate pair may change places between steps. At the last step
  of a pass the residual test decides alone, as before, and the run prints a warning when
  the value was still moving. `validation/tests/test_ritz_stop_rule.py` runs the rule on 900
  random non-normal operators and counts the stops with an error above 3 times the error of
  the whole basis plus 1e-6: the residual test alone in 62 % of the operators, one step of
  rest in 26 %, three steps in 1.2 %, five steps in none. The mean stop moves from step 17 to
  step 26 of 128. The unit cases are in `validation/tests/test_ks_select.py`.
  `cylinder_re100/310_stability_direct` stops at step 50 of 50 (138 matvecs) with the
  residual 3.1e-7, and `lid_driven_re3600/310_stability_direct` still at step 114 of 120;
  both references are unchanged.
- Every exported mode is checked with one more application of the linearized map (one
  matvec for a real mode, two for a pair): `||M x - mu x|| / ||x||`. The Arnoldi residual
  equals this number only while `M V = V H + f e_k^T` holds for every column. Below
  `eigen_tol` the mode is written. Above it, but less than `eigen_tol` above the Arnoldi
  residual, it is written with a warning: the run stops at the first step where the Arnoldi
  residual falls below the tolerance, and a relation error of about 1e-9 can put the vector
  just above it (`cylinder_re100/310`: Arnoldi 9.761e-7, vector 9.772e-7). A larger
  difference means the vector is not an eigenvector: it is not written, and if it is one of
  the `schur_tgt` leading modes the run ends with an error once the files are written. A
  test build that moved one mesh node of the leading mode by 1 gave 5.6e-4 against 9.5e-7;
  the mode was rejected and the run ended with exit code 1. Row k of `Spectre_*_conv.dat` and
  mode file k always hold the same mode. `Spectre_?.info` gets `k_used=` (size of the last
  basis), `check matvecs=` and a `[vector residuals]` section with both residuals. After a
  restart with locking the flag column printed 2.2e-16 for the leading mode of
  `lid_driven_re3600/310` with `schur_tgt = 4`; its vector residual is 5.3e-7, the Arnoldi
  residual it had when it was locked.
- `lid_driven_re3600` 310 and 320 ask for 3 values, the real mode and the pair that the
  references read, not 4. With 4 the run needs a restart, 187 matvecs instead of 120, for a
  fourth value (real, sigma = -0.0503) that no reference reads. `tpjet_re2000` 311 and 321
  ask for 1: the leading multiplier is real (-1.1577), and the next entry, 0.9996, belongs to
  a cluster of slow multipliers between 0.988 and 0.9996 that had not converged after 128
  steps. The 311 reference no longer reads rows 2 and 3 of `Spectre_Hd_conv.dat` (-0.860 and
  -0.819): they were converged values below those three unconverged multipliers, not the
  second and third multipliers.

### Fixed
- Eigenvalue runs could converge to mu = 1.0000000 with a mode that is not a
  flow field: the copies that neighbouring elements hold for one grid point
  differ by O(1). This holds for every earlier version, v1.0.1 included. Copies
  on different MPI ranks pick up different round-off: the gather-scatter adds
  in another order on each rank where 3 or more ranks share a point, and MKL
  dgemm (the Krylov-Schur rotation) rounds the last `mod(nv, 12)` rows of each
  rank's array in another kernel. In the ellipse case below, all 206 bad points
  shared by 2 ranks lie in the last element of a rank. The perturbation
  solver adds a continuous increment to the old field, so it keeps such a
  difference with multiplier exactly 1.
  Below the critical Reynolds number the Arnoldi method then brings it to the
  top of the spectrum. Many ranks with few elements each make it likely: an
  AR 2 ellipse at Re = 70 (5920 elements, 64 ranks, 76 points shared by 3 or 4
  ranks) exported two such modes and not the shedding pair. `matvec` now
  averages the copies of each shared point on every vector that goes into
  the solver (`average_copies`). A pure jump field of 3.1e-3 put into the solver came out
  unchanged after 324 steps; with the averaging it comes out at 4e-19, and
  the output maximum agrees with that of the run without the jump to 5 digits. New tool `bin/c0check` finds such jumps
  in field files of older runs (see DOC.md, Troubleshooting).
- Adjoint runs on an axisymmetric mesh (`axiSymmetry = yes`) diverged. Nek5000 computes
  the adjoint convection term `(grad U)^T u` in `convop_adj` (`core/perturb.f`) with
  integration weights that have no radius (`set_dealias_rx`, `core/convect2.f`) and divides
  by the mass matrix `bm1`, which has it, so the term is too large by the factor 1/r. In
  the adjoint run of the forced jet (`tpjet_re2000/321_stability_adjoint_floquet`) the
  pressure residual grew by a factor of about 1.4 per time step after the first ten steps
  and reached 1e39 within one forcing period, and the run returned multipliers of modulus
  1e52. `nekStab_forcing` now adds the difference between the term of Nek5000 and the term with
  the radius in the weights to the forcing of every adjoint perturbation on an axisymmetric
  mesh (`axisymmetric_adjoint_term` in `forcing.f90`, the pattern of the mesh-velocity term).
  A Nek5000 that has the radius in its weights gives a difference of zero. The direct maps
  were not affected: Nek5000 selects the pointwise routine `conv1d` for an axisymmetric
  mesh. On the jet at Re = 2000 the adjoint Floquet analysis gives the leading multiplier
  -1.1573 (`k_dim` = 128, stage 321) against -1.1577 of the direct analysis. With `k_dim` = 40
  it gave -1.1579, and the same value to seven digits with a Nek5000 whose `set_dealias_rx`
  has the radius.
- `bin/makeneks` (`mks`) no longer deletes the shared `3rd_party/gslib` library when its
  compatibility test cannot run. The test compiled a small program with `mktemp -d` and the
  C compiler, and called the library incompatible when that failed for any reason. On a
  compute node whose `/tmp` was full, the library was deleted under the other builds of the
  node and every stage that started there ended in seconds. The check now stops without a
  verdict when `mktemp` fails or when the compiler cannot build a program that does not use
  gslib.
- The reference figures of four stages showed nothing or a wrong series, and the reference
  check did not notice because it reads numbers. The snapshot of `slot_fst_re495/000_dns` was the
  placeholder "No field files found", and its signal, phase and spectrum drew the rows of
  three probes as one series; `read_his` in its `plot.py` now takes every third row. The
  residual figure of `back_fstep_re500/210_baseflow_newton` had no data, and that of
  `cylinder_re30_thermal/210_baseflow_newton` showed the first Arnoldi trace only. The
  snapshot of `naca0012_re2000/000_dns` was a field at t = 9, almost empty; `find_snapshot`
  now takes the last field file of the run. All four are drawn from the rerun of the
  campaign.
- The export of the eigenmodes wrote the first N Ritz values of the list sorted by
  modulus, with N the number of converged values, and not the converged ones. A Ritz
  value that has not converged and has a larger modulus than a converged one took its
  place: the mode files (`dRe`, `dIm`, `aRe`, `aIm`, `pRe`, `pIm`) and the rows of
  `Spectre_H?_conv.dat` and `Spectre_NS?_conv.dat` then belonged to the unconverged
  value. `Spectre_H?.dat` and `Spectre_NS?.dat`, which list every value with the
  convergence flag, were right. In a full run of the examples 6 of 19 spectra were
  affected. The Floquet stage of `cylinder_re180` exported the multipliers 0.880 and
  0.872 (Arnoldi residual 9e-3) in place of the converged pair 0.233 +- 0.601i.
  `outpost_ks` now visits the converged values, selected by `ks_converged_indices`
  (`ks_select.f90`) with the test of the flag column,
  `residual < eigen_tol*max(abs(lambda), 1)`. The leading rows of a spectrum, which
  the references read, did not change.
- After a lucky breakdown (an invariant subspace found at step m < `k_dim`) the eigensolver
  computed the spectrum of all `k_dim` columns of the Hessenberg matrix, of which only m were
  filled in that pass. It now uses the m steps it has, and the export uses the size of the
  last basis.
- The direct linearized solver has the mesh-velocity term of a moving mesh. With
  `motion = user` Nek5000 solves the base flow with the term B div(w u) (`admeshv`, w
  the mesh velocity), and its perturbation solver does not add the term for the
  perturbation, so the linearized map was not the derivative of the nonlinear map.
  `nekStab_forcing` now adds B div(w u') to the forcing of every perturbation
  (`ale_mesh_velocity_term` in `forcing.f90`; `divws`, divided by B because the forcing is
  multiplied by B again). On the cylinder in forced transverse oscillation at Re = 100
  (`moving_cylinder_re100`, A = 1.08, f = 0.1643), Newton-GMRES for the forced orbit
  took the squared residual from 3.7e-3 to 6.5e-4, 2.6e-4, 9.9e-5, 4.3e-5, 3.8e-5 and
  1.9e-5 in seven iterations (factors 0.18 to 0.88) and was stopped there, with the
  GMRES solves needing more and more vectors. It now reaches 9e-10 in 11 iterations,
  with a factor 0.12 to 0.25 per iteration; the dynamic solver tolerance stops each linear
  solve at one fifth of the current residual. A check with 30 Arnoldi steps on the field of
  the forced DNS gives the same leading Floquet multipliers for the linearized map and for
  finite differences of the nonlinear map: 0.931905 and 0.931920 for the leading one, and
  a distance of at most 4.4e-4 over the eight leading ones. The solver without the term
  gives a pair of modulus 0.78 in place of 0.93, and its eight leading estimates lie 0.06
  to 0.18 from the correct ones. The adjoint, transient-growth and OTD modes have no
  mesh-velocity term; they stop at start-up when the mesh moves
  (`nekStab_validate_mode`).
- OTD modes after the first were forced with a shifted field whenever a rank held
  fewer elements than `lelt` in `SIZE`, which is the normal case (`lelt` is set
  for the smallest number of ranks of a case: 1000 elements on 16 ranks leave 62
  or 63 per rank against `lelt = 65`). `otd_generate_forces` formed the product
  of the basis with the reduced operator, `sum_j V(:,j) Lr(j,k)`, with `mxm`
  and the number of active points as column spacing, but the basis arrays are
  spaced by the compile-time number of points. Mode k received the terms of the
  modes before it shifted by several elements, and no forcing on its last
  elements. In `cubic_cavity_re1914/500_otd` the second mode left the eigenplane
  of the 310 stage by 20 degrees (sine 0.35) within 3 time units; the same run
  with `lelt` equal to the elements per rank, or with this fix, keeps the basis
  on the plane to 1e-6. Mode 1, and runs with one mode, were not affected. The
  cylinder wake stages lost little: the shifted terms, `Lr(j,k) + Lr(k,j)` for
  j < k, vanish when the real and the imaginary part of the seed are
  orthogonal. The same call pattern in `otd_compute_OTD_modes` (it fills only
  the fields that the disabled `otd_outpost_OTD_modes` writes) is fixed too.
  Both products are `dgemm` calls with the storage spacing.
- The viscous term of the OTD reduced operator is the weak form that the time
  stepper uses. `otd_construct_diffusive_terms` computed `nu grad^2 u` as two
  applications of the gradient inside each element, without a sum over the
  element faces, while the stepper works with the stiffness matrix; the two
  differ by the jump of the normal derivative across the faces. On the 1000
  elements of `cubic_cavity_re1914` at `lx1 = 6` the mean of the two diagonal
  entries of `Lr` was -4.0e-3 at t = 0, where the growth rate of the 310
  eigenpair on which the basis starts is -1.8e-4 (the mean is half the trace of
  `Lr` on the plane; each entry alone rotates with the pair). The term is
  `-nu A u / B` with `axhelm` now, and the mean is -1.832e-4 at t = 0. In 2D the
  change is small: the two entries of `cylinder_re100/500_otd` at t = 0 are
  0.124756 and 0.124753, against 0.124779 and 0.124776 before (sigma of 310:
  0.124788). `otd_compute_laplacian` stays in the module and is no longer
  called.
- The reduced operator of the OTD has no pressure term any more. The stepper keeps the
  pressure in the polynomials of degree N-2 inside each element, without continuity,
  and its incompressibility constraint makes every basis vector orthogonal to that
  space, so the pressure gradient adds `(p, div u) = 0` to every entry of `Lr`.
  `otd_construct_linear_operator` added the gradient of the pressure that `mappr`
  maps to the velocity mesh and averages across the element faces, which is outside
  that space. On `cubic_cavity_re1914/500_otd` the stepper's weak pressure term of
  the basis was 1e-13 of the entries, and the averaged gradient gave -1.9e-4 and
  -2.1e-4 for the two modes: the mean of the two entries moved from -1.83e-4 (the
  growth rate of the 310 eigenpair) to -3.86e-4 after the first step. Without the
  term the mean stays within 1e-6 of the growth rate over the first four steps. The
  basis stayed on the eigenplane of the 310 stage with the term as well (sine 8e-7):
  the term changes the mixing of the basis inside its subspace, the entries and the
  exponents of `otd_ftle.dat`. In the 2D wake of `cylinder_re100/500_otd` the term
  was 1e-7. `otd_construct_pressure_gradient_terms` stays in the module and is no
  longer called.
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
