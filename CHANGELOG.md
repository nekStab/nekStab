# Changelog

Differences between tagged versions. `RELEASE.md` is the complete change
description for the 2.0 series; entries here are the per-tag deltas.

## [2.0.0] - 2026-11-02

The 2.0 rewrite. `RELEASE.md` is the complete description. This section lists
the changes against 1.0.1.

### Known issues

- Two runs of one binary on the same input and rank count can differ in the last
  digits. gslib times its exchange methods at start-up and keeps the fastest, and
  each method adds the copies of a shared point in another order. In
  `cylinder_re100/310` (8 ranks) the leading growth rate moved by at most 5e-9.
  With the method forced, repeat runs are byte-identical. See Troubleshooting in
  `DOC.md`.

### Added

Code structure and modes
- The sources are free-form Fortran modules, except the entry points that Nek5000
  calls (`src/main.f90`, `src/usr_wrappers.f90`). `nekstab_nek_bridge` isolates
  the Nek5000 common-block imports.
- Mode selection has three routes with this priority: the string `nekstab_mode`
  (for example `'direct'`), then the mode flags, then `userParam01`. A
  `userParam01` that is not a mode code stops the run with an error that names
  the value; it must lie within 1e-4 of a code in the mode table. Solver
  dispatch uses integer mode-code constants. The seed loader of the
  eigensolver and the steady force sensitivity choose their files from the mode
  flags (`isDirect`, `isAdjoint`, `ifForceSensReal`). A seed load in a mode with
  no seed file, and a force sensitivity without either flag, stop the run with
  a message.
  `validation/tests/test_mode_selector.py` checks every legacy `userParam01`
  decimal, every mode string, the flag path, the Floquet modifier, the conflicts
  and the invalid inputs.
- Multiple passive scalars, RANS and conjugate heat transfer (`lt != lv`) work
  through the whole stability and modal pipeline. SFD, TDF and the sponge act on
  velocity, temperature and all passive scalars. BoostConv stays velocity-only
  and warns when scalars are active. `nekStab_avg` averages every scalar slot.
- RANS stability has two operators, both finite-difference Frechet operators
  with no adjoint, wavemaker or transient growth: `coupled` (velocity, k and tau
  perturbed, eddy viscosity responding; selected by `iffindiff`) and
  `quasilaminar` (velocity-only perturbation, eddy viscosity prescribed from the
  base flow; selected by `ifquasilaminar` together with `iffindiff`, paired
  with `ifKEnorm`). The prescribed viscosity is held in `src/base_mut.f90`
  (module `nekstab_base_mut`). The resolved Reynolds stress can be applied as a
  frozen force: `reynolds_commit` writes the moments and `reynolds_enable`
  arms the load (`src/reynolds_force.f90`).
- `thermal_norm_weight` scales the temperature term of the Krylov and Newton
  inner product. `thermal_norm_mode` selects `TN_MANUAL` (default), `TN_AUTO`
  (weight from the energy ratio and `thermal_buoyancy_coeff`) or `TN_CLIP`
  (`TN_AUTO` clipped to `[thermal_norm_min, thermal_norm_max]`).
- General Boussinesq buoyancy API (`nekStab_set_buoyancy`, `nekStab_qvol`) with
  the adjoint transpose terms. Thermal adjoint stability needs only the `.usr`
  file. Thermosyphon Ra = 500: adjoint 0.1022272, direct 0.1022191.
- Free-stream turbulence inflow as a reusable `src/` module (`nekstab_fst`) with
  a self-describing binary mode file.
- Modal analysis: POD, projected DMD and SPOD (batch and streaming), sharing the
  batch inner-product code. The modal outposts write the scalar field.
- Dynamic Mode Tracking (DMT) base flow, mode 1.3 (`dmt`). The band-pass filter
  has zero gain at zero frequency only, so mode 1.3 holds steady states, not
  periodic orbits. It skips the shortened last step and stops on any other
  change of the time step. The converged `BF_*` file is written in double
  precision.
- DMT filter bank, mode 6.4 (`userParam01 = 6.4`, string `'dmt_bank'`, flag
  `ifdmtbank`). Second-order filters, one per chosen frequency, run on the DNS
  without a force and return the spatial structure of the flow at their
  frequency as a complex field (`bRe`, `bIm`). Frequency 0 gives the mean flow.
  Register the frequencies in `nekStab_usrchk` with `dmt_bank_add(omega, beta)`,
  at most 8. Each filter is a damped oscillator driven by the flow and
  integrated exactly for an input that varies linearly between the steps
  (`src/dmt_filter.f90`). The amplitude comes from the filtered signal and its
  derivative, which passes a constant mean flow with no gain; the combination
  `beta*(q + i*omega*p)` has the gain `i*beta/omega` for it. `dmt_bank_options`
  sets the length of the history `dmt_bank.dat` (norms of the flow, of what no
  filter holds, and of each amplitude), the output of the filtered fields at
  each output step, and a window over which the amplitude is averaged. A filter
  also passes about `beta/omega` of the neighbouring frequencies; the average
  over whole periods removes that ripple (a test with three harmonics on a mean:
  5 %, 15 % and 45 % at the last step, 4e-5, 5e-5 and 1e-4 on the average of 10
  periods).
- OTD (mode 5) about a frozen base flow, and OTD about a base flow that the DNS
  advances: `userParam01 = 5.1`, string `'otd_dns'`, flag `ifotdDNS`. With
  `otd_FTLEPeriod` set, `otd_ftle_period.dat` holds one line per completed
  period, the exponents averaged over exactly that period. A basis vector of
  zero norm stops the run with a message.
- Floquet energy budget (mode 4.11) for time-periodic base flows, and the
  perturbation kinetic energy budget (mode 4.1) for steady flows. See Changed
  for what the budget counts.
- Fourier representation of a periodic orbit, off by default (`iffourierorbit`
  in `nekStab_usrchk`; `fourier_kmax` fixes the number of modes). A converged
  Newton orbit writes `FOU<session>0.f*` and `<session>.fourier`; Floquet, the
  energy budget and OTD read `U(t)` from those modes instead of storing the
  orbit. `K` follows a doubling rule on the missing energy, at most 128 modes.
- Zero-crossing period probe (`src/probes.f90`) as a DNS option:
  `ifzerocross = .true.` in `nekStab_usrchk`, probe position `(xck, yck, zck)`,
  component `zc_component` (1 = vx, 2 = vy [default], 3 = vz). The crossing time
  is interpolated linearly between steps. The mean is evaluated over the
  previous period. A grid point shared on a rank boundary is owned by the lowest
  rank. The reported defect is `||u(t)-u(t-T)||/||u(t)||` with the mass-weighted
  norm. When a period closes with a new level above the signal, the next
  crossing of that level starts the next period and closes none.
- `iforthocheck` (default `.false.`, set in `nekStab_usrchk`) checks the
  orthogonality of the Krylov vectors. The eigensolver writes the norm of each
  vector and the inner products of each pair to `orthonormality.dat`; GMRES
  prints the mean |<Q(k),Q(j)>| at each iteration. On `cylinder_re100` 310 the
  largest product is 3.5e-15 and on the 210 Newton run 3e-17.
- `schur_stop` (default 0) selects the stop rule of the eigensolver. `0` stops
  on the first `schur_tgt` entries of the list sorted by modulus. `1` stops when
  any `schur_tgt` values have converged. Below the critical Reynolds number,
  slow non-oscillating modes can lead the list: on an AR 2 ellipse at Re = 70
  the shedding pair converged after about 540 calls, while the leading real mode
  still had a residual of 4e-3. Set it in `nekStab_usrchk`; `Spectre_d.info`
  lists it.
- `ts_gmres` returns the squared residual of its solution in the optional
  argument `final_res`. The steady force sensitivity (modes 4.41 and 4.42) uses
  it to print a warning when GMRES ends above its target of 1e-6; the field is
  still written.
- Newton residual backtracking (a non-monotone Armijo line search) with a
  damped period correction for periodic orbits. It is on by default; set
  `ifnewton_backtrack = .false.` in the `.usr` to disable it. It rejects
  overshooting and NaN steps, evaluates each trial at its damped period, and
  grows the orbit storage when a trial orbit is longer than the allocation.
  When no damped step reduces the residual, the current state is held, and six
  iterations without progress end the loop and write the best state.
- Mesh-velocity term of a moving mesh in the direct linearized solver. With
  `motion = user` Nek5000 solves the base flow with the term B div(w u) (`w` is
  the mesh velocity) and does not add it for the perturbation.
  `nekStab_forcing` now adds B div(w u') to the forcing of every perturbation
  (`ale_mesh_velocity_term` in `forcing.f90`). On the cylinder in forced
  transverse oscillation at Re = 100 (A = 1.08, f = 0.1643), Newton-GMRES for
  the forced orbit reaches the squared residual 9e-10 in 11 iterations. Without
  the term it stalled at 1.9e-5 after seven. A check with 30 Arnoldi steps gives
  the same leading Floquet multipliers for the linearized map and for finite
  differences of the nonlinear map: 0.931905 and 0.931920 for the leading one,
  at most 4.4e-4 apart over the eight leading ones. The adjoint, transient
  growth and OTD modes have no mesh-velocity term and stop at start-up when the
  mesh moves (`nekStab_validate_mode`).
- Start-up checks (`nekStab_check_setup`). A Floquet run (direct, adjoint or
  transient growth) stops when the restart file (the orbit) was written at
  another polynomial order than `lx1` in `SIZE`; Nek5000 interpolates such a
  file, and the multiplier 1 of the phase mode is lost. The steady stability
  modes print a warning for the same mismatch. Newton runs are not checked: they
  may start from a seed of another order. A stability mode with an active
  sponge prints a note that the sponge also forces the base flow. A Floquet run
  prints the orbit closure `||U(T)-U(0)|| / ||U(0)||` of the base flow
  integrated over the period of the base-flow time stamp, and warns above 1e-3.
  The orbit of `flip_flop_re62` closes to 8e-6; with a period 9 % too long it
  closes to 8e-3. A Floquet run also warns when no converged multiplier lies
  within 1e-3 of 1; this test is weaker (with the wrong period the converged
  list still held two multipliers within 1e-11 of 1). A Newton or Krylov mode
  with `residualProj = yes` in `[VELOCITY]` prints a note that the velocity
  residual projection can stop the run with a non-finite norm. A warning
  appears when the first line of `<session>.his` is a multiple of `lhis`.
- Newton for periodic orbits (modes 2.1 and 2.2) prints a warning when
  `filtering = explicit` is set. The filter acts only while the base flow
  advances, so the Newton iterations after the first are unfiltered with
  `ifstorebase = .true.` and filtered with `.false.`. The two settings converge
  to different orbits: `flip_flop_re62` gave the periods 8.815177 and 8.815528.
  Newton for periodic orbits keeps the time-step count of its first iteration.
  A warning appears when the kept count gives a CFL number above 1.5 times
  `targetCFL`. The runs with and without a stored orbit agree to 4.2e-6 in the
  period and 1.1e-7 in the Floquet exponent.
- Tools and scripts: `bin/c0check` finds element-face jumps in field files
  (DOC.md, Troubleshooting). `scripts/check_against_ref.py` checks a stage
  against its `ref/reference.json` and fails a run whose log holds a line that
  starts with `ERROR:`; its run gate skips the `DT/DTCFL/DTFS/DTINIT` line,
  which some builds print with `NaN` for a fixed time step.
  `scripts/inspect_nek_field.py` prints the word size, the record code and the
  file size of a field file and reads Nek5000 field files in 2D and 3D, single
  and double precision, without pymech; it keeps the pymech statistics when its
  own header parser rejects the header. `scripts/compiler_battery.py` exports
  the tracked files of an example stage, builds it under nine compiler profiles
  (gfortran at -O0, -O2 and -O3, a portable build, floating-point traps and
  signalling NaN; ifx at -O0, the default level and a portable build), runs it,
  and lists every reference quantity under each profile with the spread. Seven
  stages (cylinder DNS, thermosyphon Newton and direct, BoostConv, flip-flop
  wavemaker and Floquet, dynamic Newton) pass the run gate and the reference
  check under every profile; the relative spread of the field statistics is
  1.7e-4 or less, and the cylinder DNS agrees to 5e-11.
- Validation and reproducibility. Each stage holds a `ref/` folder with
  `reference.json` and figures. The gallery
  (`https://nekstab.github.io/nekStab/validation/`) is built from the stage
  folders of `example/` with `python validation/build_gallery.py`. A card shows
  the figures of `ref/` and the quantities of `ref/reference.json`, with the keys
  `summary` (one line) and `description`. The page has no date or machine name,
  so the file in git is the output of the script. The title and blurb of each
  family come from the first heading and paragraph of
  `example/<family>/README.md`. A stage number with no lane (first digit 7, 8
  or 9) stops the build with an error that names the folder.
  `validation/tests/test_gallery_complete.py` fails for a stage without
  reference, summary, description, figure or README, for a link that does not
  resolve, for a deferred, blocked or pending card, and for a page that differs
  from the generator output. `validation/tests/test_references_not_vacuous.py`
  exports `example/` of HEAD without running anything and requires that no
  quantity of any `ref/reference.json` passes; a second test requires every
  quantity to have a numeric value, a source and exactly one tolerance.
  `validation/serve.py` starts a job through one route, `/api/run`, and allows
  the stage folders for the recompile and resubmit buttons.
- Example layout: numbered stages (`000_dns` to `6x0_modal_*`) for all families
  and `example/_templates`, with no runnable files at the root of a family.
  Families: `cylinder_re100`, `cylinder_re180`, `cylinder_re30_thermal`,
  `cylinder_re1m`, `flip_flop_re62`, `thermosyphon_ra500`, `lid_driven_re3600`,
  `back_fstep_re500`, `naca0012_re2000`, `slot_fst_re495`, `poiseuille_re5k`,
  `cubic_cavity_re1914`, `cubic_cavity_re1920`, `cubic_cavity_re1950`,
  `moving_cylinder_re100` and `tpjet_re2000`. The unified isothermal cylinder
  mesh has 2128 elements. Stages that call `hpts` write their probe list
  `<session>.his` in the run script, because `hpts` reads the probe list from
  that file and `.gitignore` excludes `*.his`. Run scripts use the default Slurm
  partition. `example/cylinder_re1m` (2D RANS) is experimental: finite-difference
  operator only, no reference data and no gallery card.
- `example/cylinder_re100/130_baseflow_dmt`: the Re = 100 base flow by DMT. With
  omega_c = 0.7337 (the leading mode of stage 310), beta = 1 and chi = 0.5 the
  residual falls below 1e-7 at t = 766, and the field agrees with the Newton
  solution to 4e-5 (minimum of `vx`) and 8e-6 (norm). The README states that DMT
  filters one frequency band, so a steady flow is the only state without force.
- `cylinder_re100/120_baseflow_boostconv` (Re = 50, defaults `bst_snp = bst_skp =
  10`) meets its gate 1e-6 at t = 159.4. Its base flow differs from the Newton
  solution at Re = 50 by 2e-5 at most in the extrema and 2.5e-4 in the norm. At
  Re = 100 neither the defaults nor 30 and 5 reach the gate by t = 500 (the stage
  README has the table).
- `cylinder_re100/140_baseflow_meanflow` and `315_stability_direct_meanflow`:
  the mean flow of the shedding Re = 100 wake and the direct stability about it.
  The DNS averages velocity and pressure over 25 whole periods, found from the
  upward zero crossings of the integral of v over the domain, writes the mean as
  a base-flow file and stops. The mean flow has a recirculation length of 1.42
  diameters; the steady flow has 6.15. The leading eigenvalue about the mean
  flow is sigma = 0.0022, omega = 1.0290, within 1.2 % of the DNS frequency
  1.0417; the steady flow gives omega = 0.7337. Marquet & Carini (2020,
  arXiv:2008.10002) give 0.0020 and 1.0322 for the mean flow at Re = 100.
- `cylinder_re100/412_postproc_steady_force_sensitivity` (mode 4.41) converges:
  the GMRES solve ends at 9.2e-7 after 2 restarts against a target of 1e-6. The
  case file changes only the outflow `O` into `v`, as stage 320 does. The
  reference holds the residual as a gate and the norms of the field. Mode 4.42
  has no example and is not validated in 2.0.0.
- `cylinder_re100/500_otd` runs without sponge, as the 310 stage does, and
  reproduces its growth rate 0.12479 to 1e-4 (seeded with the leading
  eigenmode), or to 1e-3 after about 100 time units from white noise.
- `cylinder_re100/610_modal_dmd` and `620_modal_spod`: the three modal stages
  run the DNS leg (t = 150 to 200) and the modal analysis in one script from
  tracked files, and their references check the POD energy, the DMD Strouhal
  numbers and the SPOD eigenvalues.
- `cylinder_re180`: `000_dns_seed` is one DNS at Re = 180 from rest. The DNS
  lands on the limit cycle, so its last field and wake-probe period start the
  Newton stage directly. The Newton stage runs without the explicit filter. The
  Floquet stages run at the polynomial order of the Newton stage (the multiplier
  of the phase mode then misses 1 by 1.8e-6) and with the sponge off.
- `cubic_cavity_re1914` runs with the uniform lid `u = 1`, as in the thesis. The
  Hopf pair is present: sigma = -1.8e-4, omega = 0.5868 (St = 0.0934); direct and
  adjoint agree to 2e-4 in sigma. The flow is stable, and a perturbation still
  gains a factor 229 in energy over the horizon tau = 10 (stages
  `330_transient_growth` and `411_postproc_wavemaker`; G = 229.19, the second
  gain is 215.96). The energy ratio of the stored response to the stored
  perturbation is 0.3 % below G; this measures how far the adjoint integration
  departs from the transpose of the forward one on this mesh. The wavemaker of
  the leading pair has its maximum, 45.0, at (x, y, z) = (-0.365, -0.244, 0) and
  a second maximum of 44.8 at (-0.307, 0.333, 0), so the position of the global
  maximum is not stable. The 500_otd stage keeps the OTD basis on the eigenplane
  of the 310 stage to 1e-6 (sine 8e-7) and the mean of the two diagonal entries
  of the reduced operator at -1.832e-4 at t = 0, the growth rate of the 310
  eigenpair.
- `cubic_cavity_re1920`: the 3D cubic cavity at Re = 1920, above the first Hopf
  bifurcation. A DNS at Re = 1900 gives the start field, Newton-GMRES finds the
  unstable steady flow in 12 iterations, and the direct stability stage finds one
  unstable pair, sigma = +5.84e-4, omega = 0.5868. At Re = 1914 the same pair has
  sigma = -1.83e-4; linear interpolation of the two puts the zero at Re = 1915.4
  on this mesh (1916.63 in the thesis).
- `cubic_cavity_re1950`: DNS, POD, DMD and energy-budget stages on the limit
  cycle of the 3D cavity at Re = 1950 (period 10.78909). `000_dns` runs 100 time
  units on the cycle and writes 200 snapshots; `600_modal_pod` and
  `610_modal_dmd` read 100 of them, 4.6 periods. Two POD modes hold 94.3 % of the
  energy and turn at 1.0005 times 2 pi / T; the next four pairs turn at 1.9998,
  2.9989, 3.9981 and 4.9975 times. The DMD modes lie on the harmonics 1 to 6 of
  the cycle, within 4.2e-4 in units of 2 pi / T, and none grows (largest modulus
  0.999978). `411_postproc_energy_budget` runs the Floquet energy budget on the
  closed cavity, where P - D = sigma E holds exactly. The energy of the marched
  mode gives sigma to 0.17 % (2.90305e-3 against 2.907859e-3), and production
  and dissipation balance to 0.012 % of the dissipation. P/E and D/E are 175
  times larger than sigma, so the exponent from the budget, the difference of the
  two, is 58.5 % high. The README says so, and the reference does not check that
  exponent.
- `cubic_cavity_re1950/630_modal_dmt` and `moving_cylinder_re100/630_modal_dmt`
  run the DMT filter bank with four filters: the frequency of the cycle or of the
  forcing, its second and third harmonic, and the mean flow. On the cavity cycle
  (15 periods, 40 ranks, 1292 s) the amplitudes rebuild the last field to 0.227,
  0.059 and 0.019 of its fluctuation after one, two and three harmonics; the
  first POD pair of stage 600 leaves 0.239 after the fundamental. On the forced
  orbit of the oscillating cylinder (30 periods, 40 ranks, 643 s) the amplitudes
  rebuild the orbit of stage 210 to 0.401, 0.236 and 0.145 in the near field. On
  the cylinder wall the amplitude of the first filter is 1.1148 against the wall
  speed A 2 pi f = 1.1149, and the other filters return less than 2e-6.
- `moving_cylinder_re100` (010 to 311): the cylinder in forced transverse
  oscillation at Re = 100, A = 1.08, f = 0.1643 on a moving mesh. The DNS locks
  to the forcing; Newton-GMRES finds the orbit in 11 iterations (squared residual
  3.7e-3 to 9.7e-10, 3972 s on 40 ranks); the direct Floquet stage finds one
  converged multiplier, mu = 0.92988, real and inside the unit circle, so the
  lock-in is stable, and the DNS relaxes onto the orbit by 0.9252 per period,
  0.5 % below |mu|.
- `tpjet_re2000` (140 to 321): the harmonically forced axisymmetric jet at
  Re = 2000 (St = 0.6, T = 1.66667, 4800 elements, `lx1 = 6`, 16 ranks).
  Time-delayed feedback (mode 1.4, stage 140, 4880 s) starts from the Re = 1900
  orbit; its residual rises to 3.2e-4 at t = 21.5 and falls as t^-2.49 to 6.5e-8
  at t = 600, and the field lies at 1.9e-3 of the Newton orbit (the seed at
  1.05e-2). Stage 210 starts from the Re = 1900 field
  (`BF_Re1900_tpjet0.f00001`), not from a converged field, and finds the orbit in
  5 iterates (squared residual 3.5e-6 to 9.7e-9, 66 Arnoldi vectors, 2663 s). The
  direct Floquet stage 311 finds one real multiplier outside the unit circle,
  mu = -1.1577 (a period doubling, sigma = 0.0878; the run stops at step 57 of
  128 with 13 Ritz values converged, orbit closure 4.6e-5, 2408 s). The other
  converged values lie on an arc inside the unit circle that moves with the step
  at which the run stops (-0.9241 +- 0.0458i after 57 steps, -0.8599 after all
  128). The adjoint stage 321 gives mu = -1.1573 (step 56 of 128, 22 converged,
  1890 s) and an adjoint mode at the nozzle lip, z = 0.64 and r = 0.48, where the
  direct mode peaks at z = 20.6. The damped multipliers of the adjoint map are
  not those of the direct map (the second is -0.9063 after 56 steps, -0.9050
  after 128 and -0.9654 at `k_dim` = 40), so the references check the leading
  multiplier only. At `lx1 = 5` the Floquet spectrum at Re = 1900 had a complex
  pair of modulus 2.35; at `lx1 = 6` the leading multiplier is real, -1.139, a
  period doubling as in the AMR_Krylov paper, so every stage uses `lx1 = 6`.
- `naca0012_re2000/210_baseflow_newton/fp` starts from the last DNS snapshot
  (`BF_seed_naca00120.f00001`, t = 70) and needs 5 iterations (squared residual
  6.1e-7 to 9.5e-10, 956 s on 8 ranks). The Newton references of the airfoil and
  the jet check the first residual and the iteration count, so a run that starts
  from a converged field fails. The sponge forces the velocity toward the field
  that the run holds at its start, so the airfoil base flow depends on its seed
  inside the sponge zone (`DOC.md`, "Sponge reference field"). It differs by
  1.1e-2 (relative l2 norm of the velocity) from the base flow that the
  stability stages ship, and the leading growth rate of the direct stage is
  -0.0643 on it instead of -0.0665. The Newton stage runs at any rank count up
  to 40 with `[VELOCITY] residualProj = no`; the residuals are identical to
  seven digits at 8, 16, 24 and 32 ranks.
- `lid_driven_re3600`: the closed lid cavity. Its 310 and 320 stages ask for 3
  values (the real mode and the pair that the references read). With 4 the run
  needs a restart, 187 matvecs instead of 120, for a fourth value that no
  reference reads.
- Floquet stages for `flip_flop_re62` (210, 311, 321, 411): T = 8.815129,
  sigma = 7.900166e-3, and the 210 reference checks the step count (1199). The
  210 stage sets no explicit filter, so its Newton and Floquet stages use the
  same unfiltered equations. Adjoint Floquet spectra match the direct ones:
  0.0078864 +/- 0.1407585i against 0.0078930 +/- 0.1407626i (the adjoint run
  replays the stored base-flow orbit time-reversed).
- Also with references: backward-facing step Re = 500 transient growth, validated
  against Barkley, Blackburn & Sherwin (2008) at six time horizons; moving
  cylinder DNS; NACA0012 DNS ringdown and wavemaker stages; thermosyphon and
  flip-flop wavemaker stages.

### Changed

- Krylov-Schur and Arnoldi (against 1.0.1): heap allocation, `dgemv`, eigenvalue
  locking, adaptive restart, conjugate-pair handling, combined absolute and
  relative convergence, and a stop at a lucky breakdown (an invariant subspace
  found at step m < `k_dim`). After a lucky breakdown the eigensolver uses the
  m steps it has, and the export uses the size of the last basis.
- The eigensolver stops on the leading eigenvalues, inside the Arnoldi pass. It
  tests after every Arnoldi step and stops at the first step where the first
  `schur_tgt` entries of the list sorted by modulus have converged
  (`ks_leading_converged`; a conjugate pair counts as two). `k_dim` is the
  largest basis, not a fixed cost; `schur_tgt <= 0` still runs one full pass.
  Before, the run stopped when any `schur_tgt` Ritz values had converged, tested
  after all `k_dim` steps, and could stop with the leading values unconverged
  and damped values converged below them (cylinder at Re = 42: rows 1 and 2
  unconverged, the converged shedding pair in rows 3 and 4). In a rerun of the
  examples 20 of 24 eigenvalue runs did exactly `k_dim` matvecs and ended with 13
  to 83 converged values where 2 to 6 were wanted. One function, `ks_converged`
  (`residual < eigen_tol*max(abs(lambda), 1)`), makes the test for the stop,
  the flag column, the export and the Arnoldi checkpoint, and rank 0 makes the
  stop decision for all ranks.
- The stop also waits until the leading Ritz values stand still. A residual
  below `eigen_tol` bounds the error of the eigenvalue only through its
  condition number, which is large for a non-normal operator. The adjoint
  Floquet run of `tpjet_re2000/321` stopped at step 37 of 128 with the residual
  1.1e-6 and the multiplier -1.16222, against -1.15792 of the direct stage; the
  condition number of that multiplier is 4e5. A leading value now has to move by
  less than `eigen_tol*max(abs(mu), 1)` in each of the last 5 steps
  (`ks_stationary`, `ks_leading_settled`); a conjugate pair may change places
  between steps. At the last step of a pass the residual test decides alone, and
  the run prints a warning when the value was still moving.
  `validation/tests/test_ritz_stop_rule.py` runs the rule on 900 random
  non-normal operators and counts the stops with an error above 3 times the
  error of the whole basis plus 1e-6: the residual test alone in 62 % of the
  operators, one step of rest in 26 %, three steps in 1.2 %, five steps in none.
  The mean stop moves from step 17 to step 26 of 128. The unit cases are in
  `validation/tests/test_ks_select.py`. `cylinder_re100/310_stability_direct`
  stops at step 50 of 50 (138 matvecs) with the residual 3.1e-7, and
  `lid_driven_re3600/310_stability_direct` at step 114 of 120. The jet stages
  `tpjet_re2000` 311 and 321 ask for 1 value: the leading multiplier is real
  (-1.1577), and the next entry, 0.9996, belongs to a cluster of slow
  multipliers between 0.988 and 0.9996 that had not converged after 128 steps.
- Every exported mode is checked with one more application of the linearized map
  (one matvec for a real mode, two for a pair): `||M x - mu x|| / ||x||`. The
  Arnoldi residual equals this number only while `M V = V H + f e_k^T` holds for
  every column. Below `eigen_tol` the mode is written. Above it, but less than
  `eigen_tol` above the Arnoldi residual, it is written with a warning
  (`cylinder_re100/310`: Arnoldi 9.761e-7, vector 9.772e-7). A larger difference
  means the vector is not an eigenvector: it is not written, and if it is one of
  the `schur_tgt` leading modes the run ends with an error once the files are
  written. Row k of `Spectre_*_conv.dat` and mode file k always hold the same
  mode. `Spectre_?.info` gets `k_used=` (size of the last basis), `check
  matvecs=` and a `[vector residuals]` section with both residuals.
- Newton, GMRES and Arnoldi guard against failures. The Newton residual log
  starts `k_sum` at 0. The closing message reports success only when the
  residual is below the target; a run that hits the iteration cap, stagnates or
  diverges is not reported as a success. Newton aborts when the first residual
  is exactly zero or not finite, at `Computed residual:`. Krylov normalization
  rejects a non-finite norm (NaN or Inf) before scaling and reports the squared
  norm. GMRES does not read one past `yvec` when Arnoldi uses every column
  without reaching its tolerance. Arnoldi keeps the Hessenberg columns already
  built across a Krylov-Schur restart. Newton uses an adaptive linear tolerance
  proportional to the residual (Eisenstat-Walker) with stagnation and
  divergence guards.
- The first residual of a steady Newton run is computed like the later ones. The
  steady linearized solver runs with `ifbase = .false.`, `q_filter` skips the
  velocity filter in that state, and the first linearized sweep left the flag
  false for the later nonlinear evaluations. The nonlinear map of a steady
  Newton run now runs unfiltered from its first evaluation, like the Jacobian.
  Converged base flows are unchanged, because they were fixed points of the
  unfiltered map. Floquet and periodic-orbit runs with an explicit filter still
  filter only the first sweep (see the warning above and `DOC.md`).
- Energy budget (modes 4.1 and 4.11). The dissipation is the weak form
  `nu int |grad q|^2` with element-local gradients, the term that the time
  stepper applies. The budget sigma of `flip_flop_re62` 411 closes to 0.40 %;
  on `cubic_cavity_re1950` 411 the gap of the balance is 0.012 % of the
  dissipation. The rest of the closure error is first order in `dt`. The
  balance adds the damping of the sponge force, `S_bar`, which is printed and is
  line 15 of `PKE_floquet_*`; the energy-ratio exponent is line 16. The flux
  through an open outflow is not counted: keep the sponge in front of it. The
  period is the time stamp of `BF_*0.f00001`, as in the Floquet eigensolve, and
  a warning appears when the `.par` value differs. Mode m reads row m of
  `Spectre_NSd_conv.dat`. On the closed `lid_driven_re3600` cavity the closures
  are 5.1e-4 for the real mode and 3.5e-4 for the complex mode. The log prints
  `sigma_r (energy)`, the exponent from the energy of the marched mode after one
  period, which has no space derivative in it.
- Transient growth writes the true optimal response (`ore`): the response call
  switches to the direct map for one step and restores the flags and the `evop`
  tag. The mass-weighted energy ratio of response to perturbation equals G
  (7.69523 against 7.695277 on the backward-facing step).
- The Ritz residual written to the spectrum is floored at
  `eps*max(|lambda|, 1)`, so a fully locked mode reports about 1e-16 and not
  values below float64 precision.
- Drag is computed directly on wall (`W`) faces with the correct scaling. The
  Floquet mode animation writes the correct deformed output. BoostConv QR
  re-orthogonalization uses modified Gram-Schmidt. Sponge strength is applied to
  every forcing branch (perturbation, temperature and velocity).
- MPI reduction batching: Gram matrices, projections and cross-spectral
  densities use BLAS with a single allreduce; CGS2 orthogonalization needs 2
  allreduces per Arnoldi step, independent of `k`. Temperton GPFA FFT is
  replaced by an FFTW3 interface.
- A run that stops on an error ends on all ranks: `nekStab_error` prints the
  message and calls `MPI_Abort`, so the run ends within seconds with exit code 1.
  Before, it called the collective `nek_end` from inside `if (nid == 0)` blocks
  and the job hung until its time limit. Orbit allocation failures report the
  requested size. Time-step and energy-budget divisions are guarded, and
  per-step solver logging is throttled (orbit-replay logfiles shrink by two
  orders of magnitude).
- Builds survive configuration changes in non-interactive shells (Slurm, CI);
  include dependency edges rebuild stale objects after `NEKSTAB.inc` or `SIZE`
  changes, and `makeneks` does not rewrite `NEKSTAB.inc` on every build. The
  long common blocks in `src/NEKSTAB` are split for the 132-column limit of
  fixed-form ifort. `bin/makeneks` (`mks`) does not delete the shared
  `3rd_party/gslib` library when its compatibility test cannot run (for example
  with a full `/tmp`): the check stops without a verdict when `mktemp` fails or
  when the compiler cannot build a program that does not use gslib.
- `bin/neksubmit` reads the remote host, the work directory and the Slurm account
  from `NEKSUBMIT_HOST`, `NEKSUBMIT_WORKDIR` and `NEKSUBMIT_ACCOUNT`, and stops
  with a message when one is missing. The queue checks of `bin/nekStab_tools.py`
  use `$USER`.
- `exports.sh` works from any directory and shell, keeps an existing
  `NEK_SOURCE_ROOT`, refuses `bash exports.sh` (it must be sourced) and sets
  `ulimit -s`. The README explains each line and has a Disclaimer section: no
  warranty, and research software that is not validated for engineering design,
  certification or safety-critical use.
- The Newton and Arnoldi loop tools in `tools/` use the current mode codes
  (`userParam01` = 2 and 3.1). `example/cylinder_re100/110_baseflow_sfd/plot_compare.py`
  also draws Newton, BoostConv (Re = 50) and DMT, and takes the folder pattern
  as an argument.
- The stage folders of `example/` track no run output: figures at the root of a
  folder (`plot_bf.png` and the like; `ref/` holds the maintained copy), `plot.out`
  logs and the output fields `pr_`, `pi_`, `sr_`, `si_`, `tr_` and `ti_` of the
  wavemaker stages are not tracked.

### Fixed

- Eigenvalue runs with Intel MKL could freeze: the log stops and all ranks spin in
  `mpi_allreduce` after a Krylov-Schur restart. Each rank solved the small Hessenberg
  problem by itself and took its own decisions, and MKL does not promise the same
  bits on each rank. On `thermosyphon_ra500/310` (4 ranks, `ifx`) the Schur vectors
  differed by 2.0 between ranks and the count of kept Ritz values by 2; 7 of 8 runs
  froze, and in one, two ranks wrote the converged output while two went on. All
  LAPACK wrappers now return the result of rank 0. No freeze in 87 restarts of the
  three test runs, and 174 of 174 stored modes show no jump between ranks (largest
  4.4e-11). gfortran results are byte-identical.
- In a closed domain the Arnoldi method returned a constant pressure as a converged
  eigenvalue at sigma = 0: the solver leaves a constant unchanged and the norm does
  not see the pressure. `matvec` now removes the mean pressure through Nek5000
  `ortho`, which acts only when all velocity boundaries are Dirichlet. See
  Troubleshooting in `DOC.md`.
- `solver = helm` removed from the thermosyphon stages and
  `cylinder_re30_thermal/210`: with this line Nek5000 ignores the temperature
  `residualTol` and keeps 1e-7.
- Eigenvalue runs could converge to mu = 1.0000000 with a mode that is not a flow
  field: the copies that neighbouring elements hold for one grid point differ by
  O(1). This holds for every earlier version, 1.0.1 included. Copies on different
  MPI ranks pick up different round-off: the gather-scatter adds in another order
  on each rank where 3 or more ranks share a point, and MKL dgemm (the
  Krylov-Schur rotation) rounds the last `mod(nv, 12)` rows of each rank's array
  in another kernel. The perturbation solver adds a continuous increment to the
  old field, so it keeps such a difference with multiplier exactly 1. Below the
  critical Reynolds number the Arnoldi method then brings it to the top of the
  spectrum. Many ranks with few elements each make it likely: an AR 2 ellipse at
  Re = 70 (5920 elements, 64 ranks, 76 points shared by 3 or 4 ranks) exported
  two such modes and not the shedding pair. `matvec` now averages the copies of
  each shared point on every vector that goes into the solver
  (`average_copies`). A pure jump field of 3.1e-3 put into the solver came out
  unchanged after 324 steps before, and comes out at 4e-19 now. `bin/c0check`
  finds such jumps in field files of older runs.
- Adjoint runs on an axisymmetric mesh (`axiSymmetry = yes`) diverged. Nek5000
  computes the adjoint convection term `(grad U)^T u` in `convop_adj`
  (`core/perturb.f`) with integration weights that have no radius
  (`set_dealias_rx`, `core/convect2.f`) and divides by the mass matrix `bm1`,
  which has it, so the term is too large by the factor 1/r. In the adjoint run of
  the forced jet the pressure residual grew by a factor of about 1.4 per time
  step and reached 1e39 within one forcing period, and the run returned
  multipliers of modulus 1e52. `nekStab_forcing` now adds the difference between
  the term of Nek5000 and the term with the radius in the weights to the forcing
  of every adjoint perturbation on an axisymmetric mesh
  (`axisymmetric_adjoint_term` in `forcing.f90`). A Nek5000 that has the radius
  in its weights gives a difference of zero. The direct maps were not affected.
  On the jet at Re = 2000 the adjoint Floquet analysis gives the leading
  multiplier -1.1573 (`k_dim` = 128, stage 321) against -1.1577 of the direct
  analysis.
- The Krylov-Schur restart could keep all `k_dim` Ritz values and report
  `k_dim - 2` converged eigenvalues that were not converged (an ellipse with
  `k_dim = 160`: "keeping 160 of 160", then 158 converged). The cap at
  `k_dim - schur_tgt` removed single values, and the step that completes the
  conjugate pairs added them back; with no Arnoldi step after the restart, the
  residual test read the Schur form and every Ritz value but the last block
  passed. The cap now removes whole pairs (`src/ks_select.f90`), and the restart
  stops with a message if it would keep all `k_dim` values.
  `validation/tests/test_ks_select.py` fails on the old selection (22 of 25
  cases) and passes on the new one. Version 1.0.1 has no cap on the number of
  kept values at all.
- The export of the eigenmodes wrote the first N Ritz values of the list sorted by
  modulus, with N the number of converged values, and not the converged ones. A
  Ritz value that has not converged and has a larger modulus than a converged one
  took its place: the mode files (`dRe`, `dIm`, `aRe`, `aIm`, `pRe`, `pIm`) and
  the rows of `Spectre_H?_conv.dat` and `Spectre_NS?_conv.dat` then belonged to
  the unconverged value. `Spectre_H?.dat` and `Spectre_NS?.dat`, which list every
  value with the convergence flag, were right. `outpost_ks` now visits the
  converged values, selected by `ks_converged_indices` (`ks_select.f90`). The
  leading rows of a spectrum, which the references read, did not change.
- `hpts` wrote zeros for the last block of probes when the number of probes was a
  multiple of `lhis` (2 probes with `lhis = 2`). The routine is in Nek5000;
  `nekStab_check_setup` prints a warning, and `DOC.md` lists the cause under
  troubleshooting.

### Removed

- The legacy flat-layout example directories, replaced by the numbered stages.
- The hand-written validation catalog and its tools: `validation/validate.py`,
  `example/collect_plots.py`, `validation/catalog.py`, `catalog_check.py`,
  `scripts/evidence_report.py` and their tests. `scripts/check_against_ref.py`
  and `validation/build_gallery.py` replace them.
- The "School A"/"School B" RANS names. The directories `findiff/` and
  `findiff_frozen/` of `cylinder_re1m/310_stability_direct` are `coupled/` and
  `quasilaminar/`, the flag `iffrozenEV` is `ifquasilaminar`, and
  `src/frozen_mut.f90` is `src/base_mut.f90`. No aliases remain.
- Procedures with no caller: `norm_grad` (diagnostics), `k_dot_complex`
  (modal_spod, replaced by `k_gram_complex`), `base_mut_reset` (base_mut) and
  `nek_fourier_decomposition` (fourier). `compute_laplacian` and
  `compute_gradients` (energy budget) are removed with the strong-form
  dissipation. `otd_compute_laplacian` and `otd_construct_pressure_gradient_terms`
  stay in the OTD module and are not called.
- Scripts: `bin/ident`, `validation/parsers/`, `write_manifest.py`,
  `generate_thumbnail.py`, `inspect_nek_ic.py`, `nstab_compare.py`,
  `matrix_audit.py`, `audit_screenshot.py`, `plot_cross_method_tts.py`,
  `plot_cross_method_convergence.py`, and from `tools/` `x_auto_cyl.py`,
  `test_fftw.f90`, `benchmark_fftw.f90` and `install_oneAPI.sh`.

## [1.0.1] - 2021-03-20

- Fixed an extra outpost when using Krylov-Schur.
- CFL limiter threshold raised from 1 to 10.
- Sensitivity subroutine with normalization; adjoint reference modes.
- New example cases (porous, MFM); GitHub Actions compiler check.

## [1.0] - 2021-01-22

First public release: steady-state computation (SFD, BoostConv, Newton),
direct and adjoint eigenvalue problems, transient growth, and
post-processing tools for Nek5000.
