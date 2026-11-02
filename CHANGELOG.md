# Changelog

All notable changes to nekStab are listed here, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the versions follow
[Semantic Versioning](https://semver.org/). `RELEASE.md` describes the 2.0 series
in more detail and `DOC.md` is the user documentation.

## [Unreleased]

### Fixed

- Eigensolver: `eigen_tol` now binds the tolerance of the solves. Nek5000 stops its
  solves on a volume-averaged residual and a Krylov vector has norm 1 over the
  domain, so the error of a matvec grows with `sqrt(volume)`: 3.6e-7 to 6.0e-7 at
  `residualTol = 1e-9` on a domain of volume 60000, close to `eigen_tol = 1e-6`.
  The eigensolver uses the smaller of `residualTol` and
  `eigen_tol/(25*sqrt(volume))`, with a volume below 1 counted as 1, and prints it
  at start-up. Domains of volume above 1600 at `residualTol = 1e-9` change. The
  bound also applies to the mesh, temperature and passive-scalar tolerances
  (`restol`). The vector check of the exported modes runs at the tolerance of
  `eigen_tol/100`, so a good mode no longer fails it because of the solver. On that
  domain the leading sigma moves from -0.0144116 to -0.0144271. See DOC.md,
  Troubleshooting.

## [2.0.0] - 2026-11-03

### Highlights

The first major release since 1.0.1 (March 2021): more than five years of work, and a
rewrite of the toolbox. This section lists the changes against 1.0.1. The previous
release is 1.1 (February 2023), whose own changes are in the 1.1 entry below; the
comparison link at the end of this file shows the difference from 1.1.

**A new code base**
- 43 `.f90` files in `src/`, about 19 500 lines, in place of the 18 files of `core/`.
- A mode is chosen by name (`nekstab_mode`) or by the `userParam01` code of 1.x, which
  still works.
- 519 automatic tests, CI on Linux and macOS with GNU and Intel compilers.

**More correct results**
- The Krylov-Schur eigensolver stops on the leading eigenvalues and checks every
  exported mode against the true operator.
- False modes are gone: mu = 1 modes from element copies, the constant pressure in a
  closed domain.
- Intel MKL runs no longer freeze: all MPI ranks use one LAPACK result.
- Newton-GMRES has a line search. The Floquet, thermal and axisymmetric adjoints are
  corrected.

**New analysis modes**
- Dynamic Mode Tracking and a DMT filter bank.
- OTD with a base flow advanced by the DNS, POD, DMD, SPOD, and the Floquet energy budget.
- Several passive scalars, conjugate heat transfer, a moving mesh in the direct solver,
  and a finite-difference RANS operator (experimental).

**80 validated examples and a website**
- 80 example stages in 15 flow families (1.0.1 had 8 example cases), from the cylinder
  and the lid-driven cavity to the cubic cavity, the flip-flop wake, the forced jet, the
  NACA 0012 and a thermosyphon. Each stage has reference data and a check script, and
  each family has a numbered workflow from DNS to stability and modal analysis.
- The website https://nekstab.github.io/nekStab/validation/ shows the figures, reference
  values and description of every stage.

**Long runs**
- `auto_restart/` takes a long DNS to a target time through a chain of jobs: each job
  restarts from the last, the histories are joined, and a failed job stops the chain.
  Its README teaches the use on your own case.

### Upgrading from 1.x

- The `userParam01` codes of 1.x keep their meaning. A value that is not a
  mode code now stops the run with a message. The mode table is in `README.md`.
- A `.usr` file that calls nekStab routines directly needs a `use` statement for
  the module of the routine (`krylov_subspace` for `inner_product` and `norm`,
  `nekstab_vectors` for the `nop*` routines, `nekstab_io`, `nekstab_diagnostics`).
  `nekStab_torque`, `nekStab_forcing` and `nekStab_define_obj` need none.
  See "Upgrading to v2.x" in `DOC.md`.
- The sources moved from `core/*.f` to `src/*.f90`. `Nek5000clone.sh` is now
  `Nek5000setup.sh`, and `exports.sh` sets the environment.
- Cases are built with `mks` (`bin/makeneks`) in place of `cmpile.sh`.
- The examples moved from `examples/` to `example/`, in numbered stages.
- The eigensolver stops on the leading eigenvalues (see Changed), and mu = 1
  modes that came from copies of shared points are removed (see Fixed). Spectra of earlier versions can hold such modes;
  `bin/c0check` finds them in field files.

### Added

Files and folders new in the repository root
- `src/`: the Fortran sources, in place of `core/`.
- `bin/`: the commands. `mks` builds a case, `c0check` finds jumps between element
  copies in a field file, `neksubmit` submits a job, `pwtore2` converts a mesh file
  through Gmsh to `.re2`, `zero` sets the time of a field file to 0 or empties a `.his`
  file, `fprettify_src` formats the sources.
- `example/`: the numbered stages, in place of `examples/`. `_templates/` holds a
  start for a new stage of each mode.
- `auto_restart/`: long DNS through a chain of restart jobs.
- `tools/`: `newton_loop/autorun.py` and `arnoldi_loop/autorun.py` copy a base case
  and run it over a grid of Krylov dimension and integration time, to choose these
  two numbers for a new case.
- `tutorials/cylinder/`: a worked cylinder case with its PDF guide and LaTeX source.
- `scripts/`: stage helpers. `check_against_ref.py` checks a run against its
  reference, `make_case_ref.py` writes a reference, `mkstage.sh` starts a stage,
  `plot_runner.py` runs the plots, `compiler_battery.py` tests compilers.
- `validation/`: the source of the gallery website and the tests.
- `python/`: the Python package `nekstab`, which generates free-stream turbulence
  inflow.
- `assets/`: the README banners.
- `Nek5000setup.sh` clones and builds Nek5000; `exports.sh` sets the environment.
- `DOC.md` is the user documentation, `EXAMPLES.md` the catalog of stages, `RELEASE.md`
  the full description of 2.0, `FORTRAN_STYLE.md` the coding rules, and `CITATIONS.md`
  and `CITATION.cff` say what to cite. `AUTHORS` and `LICENSE` are the usual files.
- `.coderabbit.yaml` configures the automatic review of pull requests.

Tools
- `auto_restart/`: a DNS to a target time through a chain of restart jobs, on Slurm
  or on a local machine, with two ready cylinder cases (Re = 50 and 150). On a test of
  4 jobs up to t = 20 the probe stays within 4e-4 of its peak-to-peak value of an
  unbroken run.
- Validation gallery at https://nekstab.github.io/nekStab/validation/, one card per
  stage, built by `python validation/build_gallery.py`.

Code structure and modes
- The sources are free-form Fortran modules. The entry points that Nek5000 calls
  (`main.f90`, `usr_wrappers.f90`) and the Nek5000 bridge module are fixed-form,
  and the common blocks stay in an include file. The bridge `nekstab_nek_bridge`
  isolates the Nek5000 common-block imports.
- Mode selection has three routes with this priority: the string `nekstab_mode`
  (for example `'direct'`), then the mode flags, then `userParam01`. A
  `userParam01` that is not a mode code stops the run with an error that names
  the value; it must lie within 1e-4 of a code in the mode table. The seed loader
  of the eigensolver and the steady force sensitivity choose their files from the
  mode flags. A seed load in a mode with no seed file, and a force sensitivity
  without either flag, stop the run with a message.
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
  with `ifKEnorm`). The resolved Reynolds stress can be applied as a frozen
  force: `reynolds_commit` writes the moments and `reynolds_enable` arms the load.
- `thermal_norm_weight` scales the temperature term of the Krylov and Newton
  inner product. `thermal_norm_mode` selects `TN_MANUAL` (default), `TN_AUTO`
  (weight from the energy ratio and `thermal_buoyancy_coeff`) or `TN_CLIP`
  (`TN_AUTO` clipped to `[thermal_norm_min, thermal_norm_max]`).
- General Boussinesq buoyancy API (`nekStab_set_buoyancy`, `nekStab_qvol`) with
  the adjoint transpose terms. Thermal adjoint stability needs only the `.usr`
  file. Thermosyphon Ra = 500: adjoint 0.1022272, direct 0.1022191.
- Free-stream turbulence inflow as a reusable module (`nekstab_fst`) with a
  self-describing binary mode file.
- Modal analysis: POD, projected DMD and SPOD (batch and streaming), sharing the
  batch inner-product code. The modal outposts write the scalar field.
- Dynamic Mode Tracking (DMT) base flow, mode 1.3 (`dmt`). The band-pass filter
  has zero gain at zero frequency only, so mode 1.3 holds steady states, not
  periodic orbits. It skips the shortened last step and stops on any other
  change of the time step. The converged `BF_*` file is written in double
  precision.
- DMT filter bank, mode 6.4 (`userParam01 = 6.4`, string `'dmt_bank'`, flag
  `ifdmtbank`). One second-order filter per chosen frequency runs on the DNS
  without a force and returns the spatial structure of the flow at that
  frequency as a complex field (`bRe`, `bIm`). Frequency 0 gives the mean flow.
  Register at most 8 frequencies in `nekStab_usrchk` with `dmt_bank_add(omega,
  beta)`; `dmt_bank_options` sets the history file `dmt_bank.dat`, the output of
  the filtered fields and an averaging window that removes the ripple that
  neighbouring frequencies leave (`DOC.md`, "Mode 6.4").
- OTD (mode 5) about a frozen base flow, and OTD about a base flow that the DNS
  advances: `userParam01 = 5.1`, string `'otd_dns'`, flag `ifotdDNS`. With
  `otd_FTLEPeriod` set, `otd_ftle_period.dat` holds one line per completed
  period, the exponents averaged over exactly that period. A basis vector of
  zero norm stops the run with a message.
- Floquet energy budget (mode 4.11) for time-periodic base flows, and the
  perturbation kinetic energy budget (mode 4.1) for steady flows. See Changed
  for what the budget counts.
- Fourier representation of a periodic orbit, off by default (`iffourierorbit`
  in `nekStab_usrchk`; `fourier_kmax` fixes the number of modes, at most 128).
  A converged Newton orbit writes `FOU<session>0.f*` and `<session>.fourier`;
  Floquet, the energy budget and OTD read `U(t)` from those modes instead of
  storing the orbit (`DOC.md`).
- Zero-crossing period probe as a DNS option: `ifzerocross = .true.` in
  `nekStab_usrchk`, probe position `(xck, yck, zck)`, component `zc_component`
  (1 = vx, 2 = vy [default], 3 = vz). The crossing time is interpolated linearly
  between steps and the mean is taken over the previous period. The reported
  defect is `||u(t)-u(t-T)||/||u(t)||` with the mass-weighted norm.
- `iforthocheck` (default `.false.`) checks the orthogonality of the Krylov
  vectors and writes `orthonormality.dat`; GMRES prints the mean
  |<Q(k),Q(j)>| at each iteration.
- `schur_stop` (default 0) selects the stop rule of the eigensolver. `0` stops
  on the first `schur_tgt` entries of the list sorted by modulus. `1` stops when
  any `schur_tgt` values have converged. Below the critical Reynolds number,
  slow non-oscillating modes can lead the list: on an AR 2 ellipse at Re = 70
  the shedding pair converged after about 540 calls, while the leading real mode
  still had a residual of 4e-3. Set it in `nekStab_usrchk`.
- `ts_gmres` returns the squared residual of its solution. The steady force
  sensitivity (modes 4.41 and 4.42) prints a warning when GMRES ends above its
  target of 1e-6; the field is still written.
- Newton residual backtracking (a non-monotone Armijo line search) with a
  damped period correction for periodic orbits. It is on by default; set
  `ifnewton_backtrack = .false.` in the `.usr` to disable it. It rejects
  overshooting and NaN steps and grows the orbit storage when a trial orbit is
  longer than the allocation. When no damped step reduces the residual, the
  current state is held, and six iterations without progress end the loop and
  write the best state.
- Mesh-velocity term of a moving mesh in the direct linearized solver. With
  `motion = user` Nek5000 solves the base flow with the term B div(w u) (`w` is
  the mesh velocity) and does not add it for the perturbation; nekStab now adds
  B div(w u') to the forcing of every perturbation (`DOC.md`, "Moving Mesh").
  On the oscillating cylinder at Re = 100, Newton-GMRES for the forced orbit
  reaches the squared residual 9e-10 in 11 iterations; without the term it
  stalled at 1.9e-5 after seven. A check with 30 Arnoldi steps on the field of
  the forced DNS gives the same leading Floquet multipliers for the linearized
  map and for finite differences of the nonlinear map: 0.931905 and 0.931920.
  The adjoint, transient growth and OTD modes have no mesh-velocity term and
  stop at start-up when the mesh moves.
- Start-up checks. A Floquet run (direct, adjoint or transient growth) stops
  when the restart file (the orbit) was written at another polynomial order than
  `lx1` in `SIZE`; Nek5000 interpolates such a file, and the multiplier 1 of the
  phase mode is lost. The steady stability modes print a warning for the same
  mismatch. Newton runs are not checked: they may start from a seed of another
  order. A stability mode with an active sponge prints a note that the sponge
  also forces the base flow. A Newton or Krylov mode with
  `residualProj = yes` in `[VELOCITY]` prints a note on the non-finite norm that it can cause. A Floquet run prints the orbit closure
  `||U(T)-U(0)|| / ||U(0)||` of the base flow integrated over the period of the
  base-flow time stamp, and warns above 1e-3. A Floquet run also warns when no converged
  multiplier lies within 1e-3 of 1; this test is weaker than the closure.
- Newton for periodic orbits (modes 2.1 and 2.2) prints a warning when
  `filtering = explicit` is set. The filter acts only while the base flow
  advances, so the Newton iterations after the first are unfiltered with
  `ifstorebase = .true.` and filtered with `.false.`. The two settings converge
  to different orbits (`flip_flop_re62`: periods 8.815177 and 8.815528).
  Newton for periodic orbits keeps the time-step count of its first iteration,
  and a warning appears when the kept count gives a CFL number above 1.5 times
  `targetCFL`.
- Tools and scripts: `bin/c0check` finds element-face jumps in field files
  (`DOC.md`, Troubleshooting). `scripts/check_against_ref.py` checks a stage
  against its `ref/reference.json` and fails a run whose log holds a line that
  starts with `ERROR:`.
  `scripts/inspect_nek_field.py` prints the word size, the record code and the
  file size of a field file and reads Nek5000 field files in 2D and 3D, single
  and double precision, without pymech. `scripts/compiler_battery.py` builds and
  runs a stage under nine compiler profiles and lists every reference quantity
  under each. Seven stages pass the run gate and the reference check under every
  profile; the relative spread of the field statistics is 1.7e-4 or less.
- Validation and reproducibility. Each stage holds a `ref/` folder with
  `reference.json` and figures. The gallery
  (`https://nekstab.github.io/nekStab/validation/`) is built from the stage
  folders of `example/` with `python validation/build_gallery.py`; the page has
  no date or machine name, so the file in git is the output of the script.
  `validation/tests/test_gallery_complete.py` fails for a stage without
  reference, summary, description, figure or README, for a link that does not
  resolve, and for a page that differs from the generator output. Another test
  requires that no quantity of any `ref/reference.json` passes on the tracked
  files alone, and that every quantity has a numeric value, a source and
  exactly one tolerance.
- Example layout: numbered stages (`000_dns` to `6x0_modal_*`) for all families
  and `example/_templates`, with no runnable files at the root of a family.
  There are 16 family folders in `example/`: the 15 families of the validated
  examples and `cylinder_re1m`, which has no reference data. The unified isothermal cylinder
  mesh has 2128 elements. Stages that call `hpts` track their probe list as
  `<session>.probes`, and `mks` copies it to `<session>.his` when no `.his` exists,
  because `hpts` reads the probe list from that file and `.gitignore` excludes `*.his`.
  45 run scripts also write a one-probe list to `<session>.his` before the run. Run scripts use the default Slurm
  partition. `example/cylinder_re1m` (2D RANS) is experimental: finite-difference
  operator only, no reference data and no gallery card.
- `example/cylinder_re100/130_baseflow_dmt`: the Re = 100 base flow by DMT. With
  omega_c = 0.7337 (the leading mode of stage 310), beta = 1 and chi = 0.5 the
  residual falls below 1e-7 at t = 766, and the field agrees with the Newton
  solution to 4e-5 (minimum of `vx`) and 8e-6 (norm).
- `cylinder_re100/120_baseflow_boostconv` (Re = 50, defaults `bst_snp = bst_skp =
  10`) meets its gate 1e-6 at t = 159.4. At Re = 100 neither the defaults nor 30 and 5
  reach the gate by t = 500.
- `cylinder_re100/140_baseflow_meanflow` and `315_stability_direct_meanflow`:
  the mean flow of the shedding Re = 100 wake and the direct stability about it.
  The DNS averages velocity and pressure over 25 whole periods and writes the
  mean as a base-flow file. The leading eigenvalue about the mean flow is
  sigma = 0.0022, omega = 1.0290, within 1.2 % of the DNS frequency 1.0417; the
  steady flow gives omega = 0.7337. Marquet & Carini (2020) give 0.0020 and
  1.0322.
- `cylinder_re100/412_postproc_steady_force_sensitivity` (mode 4.41) converges:
  the GMRES solve ends at 9.2e-7 after 2 restarts against a target of 1e-6.
  Mode 4.42 has no example and is not validated in 2.0.0.
- `cylinder_re100/500_otd` runs without sponge, as the 310 stage does, and
  reproduces its growth rate 0.12479 to 1e-4 (seeded with the leading
  eigenmode), or to 1e-3 after about 100 time units from white noise.
- `cylinder_re100/600_modal_pod`, `610_modal_dmd` and `620_modal_spod` run the
  DNS leg and the modal analysis in one script from tracked files; their
  references check the POD energy, the DMD Strouhal numbers and the SPOD
  eigenvalues.
- `cylinder_re180`: `000_dns_seed` is one DNS at Re = 180 from rest. The DNS
  lands on the limit cycle and starts the Newton stage directly. The Newton
  stage runs without the explicit filter, and the Floquet stages run at its
  polynomial order and with the sponge off.
- `cubic_cavity_re1914` runs with the uniform lid `u = 1`, as in the thesis. The
  Hopf pair is present: sigma = -1.8e-4, omega = 0.5868 (St = 0.0934); direct and
  adjoint agree to 2e-4 in sigma. The flow is stable, and a perturbation still
  gains a factor 229 in energy over the horizon tau = 10 (stages
  `330_transient_growth` and `411_postproc_wavemaker`). The wavemaker of the
  leading pair has two nearly equal maxima, so the position of the global
  maximum is not stable.
- `cubic_cavity_re1920`: the 3D cubic cavity at Re = 1920, above the first Hopf
  bifurcation. A DNS at Re = 1900 gives the start field, Newton-GMRES finds the
  unstable steady flow in 12 iterations, and the direct stability stage finds one
  unstable pair, sigma = +5.84e-4, omega = 0.5868. At Re = 1914 the same pair has
  sigma = -1.83e-4; linear interpolation of the two puts the zero at Re = 1915.4
  on this mesh (1916.63 in the thesis).
- `cubic_cavity_re1950`: DNS, POD, DMD and energy-budget stages on the limit
  cycle of the 3D cavity at Re = 1950 (period 10.78909). Two POD modes hold
  94.3 % of the energy, and the DMD modes lie on the harmonics 1 to 6 of the
  cycle. `411_postproc_energy_budget` runs the Floquet
  energy budget on the closed cavity, where P - D = sigma E holds exactly. The
  energy of the marched mode gives sigma to 0.17 % (2.90305e-3 against
  2.907859e-3), and production and dissipation balance to 0.012 % of the
  dissipation. P/E and D/E are 175 times larger than sigma, so the exponent from
  the budget, the difference of the two, is 2.0 % off; the reference does not
  check that exponent.
- `cubic_cavity_re1950/630_modal_dmt` and `moving_cylinder_re100/630_modal_dmt`
  run the DMT filter bank with four filters: the frequency of the cycle or of the
  forcing, its second and third harmonic, and the mean flow. On the cavity cycle
  the amplitudes rebuild the last field to 0.227, 0.059 and 0.019 of its
  fluctuation after one, two and three harmonics.
- `moving_cylinder_re100` (010 to 311): the cylinder in forced transverse
  oscillation at Re = 100, A = 1.08, f = 0.1643 on a moving mesh. The DNS locks
  to the forcing; Newton-GMRES finds the orbit in 11 iterations (squared residual
  3.7e-3 to 9.7e-10, 3972 s on 40 ranks); the direct Floquet stage finds one
  converged multiplier, mu = 0.92988, real and inside the unit circle, so the
  lock-in is stable.
- `tpjet_re2000` (140 to 321): the harmonically forced axisymmetric jet at
  Re = 2000 (St = 0.6, T = 1.66667, 4800 elements, `lx1 = 6`, 16 ranks).
  Time-delayed feedback (mode 1.4, stage 140, 4880 s) starts from the Re = 1900
  orbit and brings its residual to 6.5e-8 at t = 600. Stage 210 starts from the
  Re = 1900 field and finds the orbit in 5 iterates. The direct Floquet stage 311
  finds one real multiplier outside the unit circle, mu = -1.1577 (a period
  doubling, sigma = 0.0878). The adjoint stage 321 gives mu = -1.1573 and an
  adjoint mode at the nozzle lip, where the direct mode peaks at z = 20.6. The
  damped multipliers of the adjoint map are not those of
  the direct map, so the references check the leading multiplier only. At
  `lx1 = 5` the Floquet spectrum at Re = 1900 had a complex pair of modulus 2.35;
  at `lx1 = 6` the leading multiplier is real, so every stage uses `lx1 = 6`.
- `naca0012_re2000/210_baseflow_newton/fp` starts from the last DNS snapshot
  (`BF_seed_naca00120.f00001`, t = 70) and needs 5 iterations. The Newton references of the airfoil and
  the jet check the first residual and the iteration count, so a run that starts
  from a converged field fails. The sponge forces the velocity toward the field
  that the run holds at its start, so the airfoil base flow depends on its seed
  inside the sponge zone (`DOC.md`, "Sponge reference field"). It differs by
  1.1e-2 (relative l2 norm of the velocity) from the base flow that the
  stability stages ship. The Newton stage runs at any rank count up
  to 40 with `[VELOCITY] residualProj = no`.
- `lid_driven_re3600`: the closed lid cavity. Its 310 and 320 stages ask for 3
  values, the real mode and the pair that the references read.
- Floquet stages for `flip_flop_re62` (210, 311, 321, 411): T = 8.815129,
  sigma = 7.900166e-3, and the 210 reference checks the step count (1199). Its
  Newton and Floquet stages use the same unfiltered equations. The adjoint pair of stage 321, 0.0078935 +/-
  0.1407687i, matches the direct pair of stage 311, 0.0079002 +/- 0.1407729i, to
  7e-6 in sigma (the adjoint run replays the stored base-flow orbit
  time-reversed).
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
- The eigensolver stops on the leading eigenvalues, inside the Arnoldi pass
  (`schur_stop = 0`, the default). It tests after every Arnoldi step and stops at
  the first step where the first `schur_tgt` entries of the list sorted by
  modulus have converged (a conjugate pair counts as two). `k_dim` is the
  largest basis, not a fixed cost; `schur_tgt <= 0` still runs one full pass.
  Before, the run stopped when any `schur_tgt` Ritz values had converged, tested
  after all `k_dim` steps, and could stop with the leading values unconverged
  and damped values converged below them (cylinder at Re = 42). In a rerun of the
  examples 20 of 24 eigenvalue runs did exactly `k_dim` matvecs and ended with 13
  to 83 converged values where 2 to 6 were wanted. One convergence test,
  `residual < eigen_tol*max(abs(lambda), 1)`, serves the stop, the flag column,
  the export and the Arnoldi checkpoint, and rank 0 makes the stop decision for
  all ranks.
- The stop also waits until the leading Ritz values stand still. A residual
  below `eigen_tol` bounds the error of the eigenvalue only through its
  condition number, which is large for a non-normal operator. The adjoint
  Floquet run of `tpjet_re2000/321` stopped at step 37 of 128 with the residual
  1.1e-6 and the multiplier -1.16222, against -1.15792 of the direct stage; the
  condition number of that multiplier is 4e5. A leading value now has to move by
  less than `eigen_tol*max(abs(mu), 1)` in each of the last 5 steps; a conjugate
  pair may change places between steps. At the last step of a pass the residual
  test decides alone, and the run prints a warning when the value was still
  moving. `validation/tests/test_ritz_stop_rule.py` runs the rule on 900 random
  non-normal operators and counts the stops with an error above 3 times the
  error of the whole basis plus 1e-6: the residual test alone in 62 % of the
  operators, one step of rest in 26 %, three steps in 1.2 %, five steps in none.
  The mean stop moves from step 17 to step 26 of 128. The jet stages
  `tpjet_re2000` 311 and 321 ask for 1 value: the leading multiplier is real,
  and the next entry belongs to a cluster of slow multipliers that had not
  converged after 128 steps.
- Every exported mode is checked with one more application of the linearized map
  (one matvec for a real mode, two for a pair): `||M x - mu x|| / ||x||`. The
  Arnoldi residual equals this number only while `M V = V H + f e_k^T` holds for
  every column. Below `eigen_tol` the mode is written. Above it, but less than
  `eigen_tol` above the Arnoldi residual, it is written with a warning
  (`cylinder_re100/310`). A larger difference
  means the vector is not an eigenvector: it is not written, and if it is one of
  the `schur_tgt` leading modes the run ends with an error once the files are
  written. Row k of `Spectre_*_conv.dat` and mode file k always hold the same
  mode. `Spectre_?.info` gets `k_used=`, `check matvecs=` and a
  `[vector residuals]` section with both residuals.
- Newton, GMRES and Arnoldi guard against failures. The Newton residual log
  starts `k_sum` at 0. The closing message reports success only when the
  residual is below the target; a run that hits the iteration cap, stagnates or
  diverges is not reported as a success. Newton aborts when the first residual
  is exactly zero or not finite. Krylov normalization rejects a non-finite norm
  (NaN or Inf) before scaling and reports the squared norm. GMRES does not read
  one past its last column when Arnoldi uses every column without reaching its
  tolerance. Arnoldi keeps the Hessenberg columns already built across a
  Krylov-Schur restart. Newton uses an adaptive linear tolerance proportional to
  the residual (Eisenstat-Walker) with stagnation and divergence guards.
- The first residual of a steady Newton run is computed like the later ones. The
  steady linearized solver runs with `ifbase = .false.`, `q_filter` skips the
  velocity filter in that state, and the first linearized sweep left the flag
  false for the later nonlinear evaluations. The nonlinear map of a steady
  Newton run now runs unfiltered from its first evaluation, like the Jacobian.
  Converged base flows are unchanged, because they were fixed points of the
  unfiltered map. Floquet runs with an explicit filter filter only the first
  sweep, because they replay the stored orbit afterwards; Newton for periodic
  orbits does the same with `ifstorebase = .true.` (see Added and `DOC.md`).
- Energy budget (modes 4.1 and 4.11). The dissipation is the weak form
  `nu int |grad q|^2` with element-local gradients, the term that the time
  stepper applies. The budget sigma of `flip_flop_re62` 411 closes to 0.40 %.
  The rest of the closure error is first order in `dt`. The
  balance adds the damping of the sponge force, `S_bar`, which is printed and is
  line 15 of `PKE_floquet_*`; the energy-ratio exponent is line 16. The flux
  through an open outflow is not counted: keep the sponge in front of it. The
  period is the time stamp of `BF_*0.f00001`, as in the Floquet eigensolve, and
  a warning appears when the `.par` value differs. Mode m reads row m of
  `Spectre_NSd_conv.dat`. The log prints
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
- A run that stops on an error ends on all ranks: it prints the message and
  calls `MPI_Abort`, so the run ends within seconds with exit code 1. Before, it
  called the collective `nek_end` from inside `if (nid == 0)` blocks and the job
  hung until its time limit. Orbit allocation failures report the requested
  size. Time-step and energy-budget divisions are guarded, and per-step solver
  logging is throttled (orbit-replay logfiles shrink by two orders of
  magnitude).
- Builds survive configuration changes in non-interactive shells (Slurm, CI);
  include dependency edges rebuild stale objects after `NEKSTAB.inc` or `SIZE`
  changes, and `mks` does not rewrite `NEKSTAB.inc` on every build. The long
  common blocks are split for the 132-column limit of fixed-form ifort. `mks`
  does not delete the shared `3rd_party/gslib` library when its compatibility
  test cannot run (for example with a full `/tmp`).
- `bin/neksubmit` reads the remote host, the work directory and the Slurm account
  from `NEKSUBMIT_HOST`, `NEKSUBMIT_WORKDIR` and `NEKSUBMIT_ACCOUNT`, and stops
  with a message when one is missing. The queue checks of `bin/nekStab_tools.py`
  use `$USER`.
- `exports.sh` works from any directory and shell, keeps an existing
  `NEK_SOURCE_ROOT`, refuses `bash exports.sh` (it must be sourced) and sets
  `ulimit -s`. The README explains each line and has a Disclaimer section.
- The Newton and Arnoldi loop tools in `tools/` use the current mode codes.
  `plot_compare.py` of the SFD stage also draws Newton, BoostConv and DMT.
- Older output files left the stage folders of `example/`: old copies of figures at the
  root of a folder (`ref/` holds the maintained copy), `plot.out` logs and the output
  fields `pr_`, `pi_`, `sr_`, `si_`, `tr_` and `ti_` of the wavemaker stages.

### Fixed

- A run built with `mks` and started with `nekbmpi` or a bare `mpiexec` stopped with
  `Cannot open history file in subroutine hpts()` in the examples that use history points.
  Only `run.local.slurm` wrote the probe list `<case>.his`, and `.gitignore` excludes `*.his`.
  Every such example folder now tracks the list as `<case>.probes`, and `mks` copies it to
  `<case>.his` when that file does not exist. `mks` never overwrites an existing `.his`, and
  it warns when the `.usr` calls `hpts` and neither file exists.
- Eigenvalue runs with Intel MKL could freeze: the log stops and all ranks spin in
  `mpi_allreduce` after a Krylov-Schur restart. Each rank solved the small Hessenberg
  problem by itself and took its own decisions, and MKL does not promise the same
  bits on each rank. On `thermosyphon_ra500/310` (4 ranks, `ifx`) the Schur vectors
  differed by 2.0 between ranks and the count of kept Ritz values by 2; 7 of 8 runs
  froze, and in one, two ranks wrote the converged output while two went on. All
  LAPACK wrappers now return the result of rank 0. No freeze in 87 restarts of the
  three test runs, and 174 of 174 stored modes show no jump between ranks (largest
  4.4e-11). `cylinder_re100/310` (8 ranks, `ifx`) froze after its third restart and now
  finishes in 217 s at the reference value. gfortran results are byte-identical.
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
  MPI ranks pick up different round-off, and the perturbation solver keeps such a
  difference with multiplier exactly 1. Below the critical Reynolds number the
  Arnoldi method then brings it to the top of the spectrum. An AR 2 ellipse at Re = 70
  on 64 ranks exported two such modes and not the shedding pair. `matvec` now averages
  the copies of each shared point on every vector that goes into the solver. `bin/c0check` finds such jumps in
  field files of older runs (`DOC.md`, Troubleshooting).
- Adjoint runs on an axisymmetric mesh (`axiSymmetry = yes`) diverged. The
  adjoint convection term of Nek5000 is integrated with weights that have no
  radius and divided by the mass matrix, which has it, so the term is too large
  by the factor 1/r. In the adjoint run of the forced jet the pressure residual
  reached 1e39 within one forcing period, and the run returned multipliers of
  modulus 1e52. nekStab now adds the difference between the term of Nek5000 and
  the term with the radius in the weights to the forcing of every adjoint
  perturbation on an axisymmetric mesh. A Nek5000 that has the radius in its
  weights gives a difference of zero. The direct maps were not affected. On the
  jet at Re = 2000 the adjoint Floquet analysis gives the leading multiplier
  -1.1573 (`k_dim` = 128, stage 321) against -1.1577 of the direct analysis
  (`DOC.md`, "Adjoint on an Axisymmetric Mesh").
- The Krylov-Schur restart could keep all `k_dim` Ritz values and report
  `k_dim - 2` converged eigenvalues that were not converged (an ellipse with
  `k_dim = 160`: "keeping 160 of 160", then 158 converged). The cap at
  `k_dim - schur_tgt` removed single values, and the step that completes the
  conjugate pairs added them back; with no Arnoldi step after the restart, the
  residual test read the Schur form and every Ritz value but the last block
  passed. The cap now removes whole pairs, and the restart stops with a message
  if it would keep all `k_dim` values. `validation/tests/test_ks_select.py` fails
  on the old selection (22 of 25 cases) and passes on the new one. Version 1.0.1
  has no cap on the number of kept values at all.
- The export of the eigenmodes wrote the first N Ritz values of the list sorted by
  modulus, with N the number of converged values, and not the converged ones. A
  Ritz value that has not converged and has a larger modulus than a converged one
  took its place: the mode files (`dRe`, `dIm`, `aRe`, `aIm`, `pRe`, `pIm`) and
  the rows of `Spectre_H?_conv.dat` and `Spectre_NS?_conv.dat` then belonged to
  the unconverged value. `Spectre_H?.dat` and `Spectre_NS?.dat`, which list every
  value with the convergence flag, were right. The export now visits the
  converged values. The leading rows of a spectrum, which the references read,
  did not change.
- `hpts` wrote zeros for the last block of probes when the number of probes was a
  multiple of `lhis` (2 probes with `lhis = 2`). The routine is in Nek5000;
  nekStab prints a warning, and `DOC.md` lists the cause under troubleshooting.

### Removed

- The legacy flat-layout example directories, replaced by the numbered stages.
- The hand-written validation catalog and its tools: `validation/validate.py`,
  `example/collect_plots.py`, `validation/catalog.py`, `catalog_check.py`,
  `scripts/evidence_report.py` and their tests. `scripts/check_against_ref.py`
  and `validation/build_gallery.py` replace them.
- The "School A"/"School B" RANS names. The directories `findiff/` and
  `findiff_frozen/` of `cylinder_re1m/310_stability_direct` are `coupled/` and
  `quasilaminar/`, and the flag `iffrozenEV` is `ifquasilaminar`. No aliases
  remain.
- Routines with no caller: `norm_grad`, `k_dot_complex`, `base_mut_reset`,
  `nek_fourier_decomposition`, and `compute_laplacian` and `compute_gradients`
  of the energy budget, which go with the strong-form dissipation.
- Scripts: `bin/ident`, `validation/parsers/`, `write_manifest.py`,
  `generate_thumbnail.py`, `inspect_nek_ic.py`, `nstab_compare.py`,
  `matrix_audit.py`, `audit_screenshot.py`, `plot_cross_method_tts.py`,
  `plot_cross_method_convergence.py`, and from `tools/` `x_auto_cyl.py`,
  `test_fftw.f90`, `benchmark_fftw.f90` and `install_oneAPI.sh`.

### Known issues

- The solver error grows with domain size. The matvec error grows with the
  domain, so on a large domain a good mode can get the vector warning or be
  refused, and the eigenvalue can differ between runs at different
  `residualTol`: on an ellipse with V = 60000 the leading sigma is -0.0144116 at
  `1e-9` and -0.0144268 at `4e-12`, both with residuals below 1e-6. The values
  differ between these runs. Workaround: lower `residualTol` by hand. Version
  2.0.0 ships no change of the tolerance; a fix is planned for a later version.
- The norm has zero weight in the sponge zone: the inner product uses the mass
  matrix without the sponge, which is zero where the sponge function is not zero.
- Stage `cubic_cavity_re1950/350_dns_unstable_orbit` keeps the reference of an
  older run (13.7 h on 16 ranks) and did not run again. The final source and the
  source of that run give a byte-identical probe file over the first 3 of the
  2000 time units (2087 time steps, 16 ranks, `gfortran`).
- RANS stability (`example/cylinder_re1m`) is experimental, and mode 4.42 is not
  validated.
- The zero-crossing period probe (`ifzerocross`) can log one wrong period after
  its level changes. When the new level is below the signal at the crossing, the
  next period spans crossings of two levels: for a sine of period 1 and mean
  -0.25, with start level 0 and dt = 0.017, the logged periods are 1.0000, 0.9597,
  0.9992, 1.0008, 1.0000. Read the period from a later line of the log. A fix is
  planned for a later version.
- Two runs of one binary on the same input and rank count can differ in the last
  digits. gslib times its exchange methods at start-up and keeps the fastest, and
  each method adds the copies of a shared point in another order. In
  `cylinder_re100/310` (8 ranks) the leading growth rate moved by at most 5e-9.
  With the method forced, repeat runs are byte-identical. See Troubleshooting in
  `DOC.md`.

## [1.1] - 2023-02-08

- Support for several passive scalars: the base-flow and sponge arrays carry one
  field per scalar, and `nekStab_torque` computes the torque on a body.
- The Gram-Schmidt process of the Arnoldi method re-orthogonalizes, the starting
  vector of the eigensolver is corrected, and a run that diverges stops.
- New files `matvec.f`, `newton_krylov.f`, `sensitivity.f`, `krylov_subspace.f`,
  `krylov_decomposition.f` and `lapack_wrapper.f` split the former `linalg.f`
  and `reader_par_nekStab.f`.
- `core/compiler.sh` selects the compiler flags automatically in place of one
  script per compiler. `bin/makeneks` builds a case.
- The Nek5000 submodule is replaced by the script `Nek5000clone.sh`.
- New example cases: lid-driven cavity, forced jet, backward-facing step,
  thermosyphon and flip-flop; the cylinder examples are reorganized into DNS,
  base-flow, stability and post-processing folders. `lpmin` is greater than 1 in
  all examples, for compilation on macOS.
- `uparam1.md` lists the mode codes. `AUTHORS`, `LICENSE` and `CITATION.cff`
  are added. The Travis configuration is replaced by GitHub Actions workflows
  for gfortran and Intel.
- Removed: the flag `ifrans`, `legacy_Krylov_Schur.f`, `p_his.py` and `viz.nek3d`.

## [1.0.1] - 2021-03-20

- Fixed an extra outpost when using Krylov-Schur.
- CFL limiter threshold raised from 1 to 10.
- Sensitivity subroutine with normalization; adjoint reference modes.
- New example cases (porous, MFM); GitHub Actions compiler check.

## [1.0] - 2021-01-22

First public release: steady-state computation (SFD, BoostConv, Newton),
direct and adjoint eigenvalue problems, transient growth, and
post-processing tools for Nek5000.

[2.0.0]: https://github.com/nekStab/nekStab/compare/1.1...v2.0.0
[1.1]: https://github.com/nekStab/nekStab/compare/v1.0.1...1.1
[1.0.1]: https://github.com/nekStab/nekStab/compare/v1.0...v1.0.1
[1.0]: https://github.com/nekStab/nekStab/releases/tag/v1.0
