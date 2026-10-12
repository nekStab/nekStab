# Release 2.0 (2026)

Complete architectural rewrite.  User-facing API (uparam encoding) is
backward compatible with 1.x; internal code is restructured for
maintainability, numerical robustness, and HPC performance.


## What changed

### Architecture

- All source wrapped in Fortran modules (39 `.f90` files, was loose `.f`)
- `nekstab_nek_bridge` module isolates all Nek5000 common-block imports
- `mode_config` module: 3-priority mode resolution (string > flags > uparam)
- `NEKSTAB.inc`: single include for all nekStab common-block variables
- Version string auto-generated from `git describe --tag`

### Multiple passive scalars and RANS support

nekStab 1.x supported only velocity and (optionally) a single
temperature field.  Version 2.0 supports an arbitrary number of
passive scalar fields (`ldimt > 1`), each with its own element count
(`nelfld`).

All core operations — inner products, norms, Gram matrices, Krylov
vector copy/add/scale, eigensolver output, checkpointing — now loop
over active scalars via `ifpsco(m-1)` and use per-field lengths
`nelfld(m+1)` rather than a single shared length.  This is needed
for:

- **RANS stability analysis** — turbulence models (e.g., k-tau)
  introduce 2 extra transported scalars (tke, tau/omega) that must
  be included in the linearized state vector.  Since the analytical
  Jacobian for RANS source terms is not available, the linearized
  operator is approximated via **finite differences** (`iffindiff =
  .true.`).  Example:
  - `example/cylinder_re1m/` — high-Reynolds cylinder wake, finite-difference direct stability
    with both RANS operators as sibling stages: `310_stability_direct/coupled`
    (full RANS Frechet) and `310_stability_direct/quasilaminar` (`ifquasilaminar`,
    hydrodynamic operator with prescribed base eddy viscosity)
- **Conjugate heat transfer (CHT)** — temperature lives on `lelt`
  elements while velocity lives on `lelv` elements (`lelt > lelv`),
  so `lt != lv`.  The old code used `lv` for temperature arrays,
  causing incorrect strides in CHT cases.
- **Multi-species flows** — any problem with passive scalars
  (concentration, additional temperatures) now propagates those
  scalars through the full stability / modal analysis pipeline.

### Numerical improvements

#### MPI reduction batching in Gram-Schmidt and inner products

The old code called `glsc3` (which wraps `MPI_Allreduce`) once per
inner product.  In an Arnoldi step with k basis vectors, Modified
Gram-Schmidt (MGS) with reorthogonalization required 2k allreduces.
For k=100 on a Lustre filesystem this was the dominant cost.

The new code gathers field data into contiguous BLAS-ready arrays,
computes all inner products via `dgemm` or `dgemv`, and does a
**single `gop` (MPI_Allreduce)** over the entire result vector.

| Routine | What it computes | MPI calls (old) | MPI calls (new) |
|---------|-----------------|-----------------|-----------------|
| `k_gram_matrix` | G(i,j) = <P(i), Q(j)> | m*n | 1 |
| `k_project` | h(i) = <Q(i), f> | k | 1 |
| `k_gram_complex` | CSD(i,j) for complex vectors | 4*nblk^2 | 1 |

Arnoldi orthogonalization now has two modes:

- **CGS2** (default, `use_cgs=.true.`): Classical Gram-Schmidt with
  one reorthogonalization pass.  Two `k_project` calls = **2 allreduces
  per Arnoldi step**, independent of k.
- **MGS+reortho** (`use_cgs=.false.`): Modified Gram-Schmidt with
  reorthogonalization.  **2k allreduces per step** — preserved for
  comparison and edge cases where CGS2 is insufficient.

The same batch primitives (`k_gram_matrix`, `k_project`) are reused
by POD (correlation matrix), DMD (snapshot Gram + shifted Gram), and
SPOD (cross-spectral density matrix), so all modal methods also
benefit from the single-gop pattern.

#### Krylov-Schur eigensolver

- **Heap allocation** for Krylov basis copies — avoids stack overflow
  on large 3D meshes (lv * k_dim can exceed several GiB)
- **Eigenvalue locking** in Schur condensation — converged eigenvalues
  are decoupled by zeroing the coupling row, preventing them from
  being perturbed by subsequent restarts
- **Adaptive restart** with three selection criteria: unit-circle
  proximity, Schur residual, and an absolute floor to prevent
  discarding all eigenvalues
- **Conjugate-pair handling** via `ensure_conjugate_pairs` — keeps
  complex conjugate pairs together during Schur reordering
- **Combined absolute/relative convergence** — an eigenvalue is
  converged if residual < tol * max(|lambda|, 1), which prevents
  spurious acceptance of small eigenvalues

#### Newton-GMRES solver (inexact Newton with adaptive tolerance)

The Newton-Krylov solver now uses an **inexact Newton** strategy:
instead of solving the linear system J*dx = -F to full precision at
each Newton step, GMRES is terminated early with a tolerance that
adapts to the current nonlinear residual.  This avoids wasting
Krylov iterations when the Newton step is still far from convergence.

The adaptive tolerance is residual-proportional: when `ifdyntol` is
enabled, the inner GMRES/Nek tolerance is set to 20% of the current
squared nonlinear residual, bounded below by the user target tolerance.
This keeps the inner solve tighter than the current Newton residual
without forcing a full extra decade of accuracy at every iteration.

Other improvements:

- **Stagnation guard** — detects when 3 consecutive Newton steps
  fail to reduce the residual and triggers early termination
- **Open-once log files** — residu_newton.dat, residu_gmres.dat,
  residu_arnoldi.dat are opened at the start and closed at the end
  instead of open/write/close per iteration (~91k metadata operations
  saved on a typical Newton run on HPC filesystems)

#### DMD projected operator

The projected operator Atilde = Sinv * V^T * G_shift * V * Sinv is
now formed with 2 dgemm calls instead of a 4-nested scalar loop.
Cost goes from O(r^2 * n^2) scalar operations to O(n^2 * r) in BLAS,
which is substantially faster for large snapshot counts.

#### Persistent workspace buffers

Module-level `allocatable, save` arrays with grow-only capacity
tracking replace per-call allocate/deallocate in inner Krylov
routines.  Reallocation only happens when a caller requests more
columns than any previous call.  After the first Krylov cycle the
buffers are reused with zero allocation cost.  A single shared gop
scratch buffer serves all three batch inner product routines.

#### CHT dimension fix

Temperature arrays now use `lt = lx1*ly1*lz1*lelt` (not `lv`)
throughout.  Pressure uses `lp = lx2*ly2*lz2*lelv`.  Passive scalar
fields use per-field extents via `nelfld(m+1)`.  This fixes incorrect
array strides in conjugate heat transfer cases where `lelt > lelv`.

#### Adjoint correctness (rc3)

Two defects fixed.  In both cases the validation criterion is that the
adjoint spectrum must match the direct one (an operator and its adjoint
share eigenvalues).

- **Floquet adjoint base-flow replay** — the time-periodic base flow is
  built once in a dedicated forward sweep, and every adjoint matvec then
  replays it **time-reversed**: step k of the backward integration uses
  orbit snapshot nsteps-k+1.  The previous implementation replayed the
  orbit forward, so each Arnoldi matvec saw a different effective
  operator.  Storing the orbit remains a memory-for-compute option
  (`ifstorebase`); when the orbit does not fit in memory the base flow
  is recomputed instead.  Validated on the flip-flop UPO at Re=62:
  adjoint pair 0.0078864 +/- 0.1407585i vs direct
  0.0078930 +/- 0.1407626i.
- **Thermal (Boussinesq) adjoint** — the adjoint of the buoyancy
  coupling has two halves: the buoyancy transpose enters the adjoint
  scalar equation as a volumetric source, and the base-flow temperature
  gradients feed the adjoint momentum production term already present
  in the Nek5000 core (`advabp_adjoint`), which requires dTdx/dTdy/dTdz
  to be populated.  nekStab now computes the base temperature gradients
  before each adjoint matvec (and per orbit step in Floquet mode).
  Validated on the thermosyphon at Ra=500: adjoint eigenvalue 0.1022272
  vs direct 0.1022191.

A third defect was found in the final rerun of the examples.

- **Adjoint on an axisymmetric mesh** — the adjoint convection term of
  Nek5000 (`convop_adj`) is integrated with weights that have no radius
  and divided by the mass matrix, which has it, so on an axisymmetric
  mesh the term is too large by 1/r and the adjoint run diverges.
  `nekStab_forcing` adds the difference to the term with the radius
  in the weights (`axisymmetric_adjoint_term`, `forcing.f90`); it is
  zero for a Nek5000 that has the radius in its weights.  Validated on
  the forced jet at Re=2000: adjoint multiplier -1.1573 vs direct
  -1.1577, the adjoint mode at the nozzle lip and the direct mode at
  z = 20 (`example/tpjet_re2000/321_stability_adjoint_floquet`).

#### General buoyancy API (rc3)

Buoyant cases declare the buoyancy direction and coefficient once in
`nekStab_usrchk`:

```fortran
call nekStab_set_buoyancy(0.0d0, 1.0d0, 0.0d0, coeff)
```

and route the standard hooks through nekStab (`nekStab_forcing` in
`userf`, `nekStab_qvol` in `userq`).  The forcing hook applies the
direct buoyancy term for both nonlinear and perturbation solves; the
qvol hook applies the adjoint transpose.  A guard aborts a buoyant
adjoint run if `userq` is not wired, instead of returning silently
wrong eigenvalues.

### New features

| Mode | uparam(1) | String | Description |
|------|-----------|--------|-------------|
| Linearized DNS | 0.1 | `linear_dns` | Time-stepping of LNSE without eigensolver |
| Dynamic Mode Tracking | 1.3 | `dmt` | Steady base flow by band-pass filtering of the velocity (Queguineur et al., 2019) |
| Floquet energy budget | 4.11 | `energy_budget_floquet` | PKE budget for time-periodic base flows |
| Animate mode | 4.50 | `animate_mode` | Reconstruct eigenmode as time series |
| Animate BF deform | 4.51 | `animate_bf_deform` | Base flow + eigenmode deformation |
| Animate Floquet | 4.52 | `animate_floquet` | Floquet mode animation |
| OTD, evolving base flow | 5.1 | `otd_dns` | Optimally Time-Dependent modes with the base flow advanced by the DNS, for orbits and unsteady flows |
| POD | 6.1 | `pod` | Proper Orthogonal Decomposition |
| DMD | 6.2 | `dmd` | Dynamic Mode Decomposition |
| SPOD | 6.3 | `spod` | Spectral POD (Welch + CSD) |
| DMT filter bank | 6.4 | `dmt_bank` | Filters on the running DNS at chosen frequencies |
| All modal | 6.0 | — | POD + DMD + SPOD together |

#### Modal analysis framework (Mode 6)

A new top-level dispatcher (`modal_analysis`) loads a snapshot
sequence, subtracts the temporal mean, and runs the enabled modal
methods.  All methods operate on the same `krylov_vector` snapshot
array and share the batch inner product infrastructure.

**POD** (`modal_pod.f90`) — Proper Orthogonal Decomposition via the
method of snapshots (Sirovich, 1987).
- Correlation matrix C(i,j) = <snaps(i), snaps(j)> / n via
  `k_gram_matrix` (single gop)
- Symmetric eigensolve via LAPACK `dsyev`
- Mode reconstruction and output in descending energy order
- POD-FFT spectral analysis: projects snapshot time series onto POD
  modes and computes Welch power spectra of temporal coefficients

**DMD** (`modal_dmd.f90`) — Projected DMD (Schmid, 2010; Tu et al.,
2014).
- Gram matrix G and shifted Gram G_shift via `k_gram_matrix` and
  `k_project` (2 gop calls for the entire snapshot set)
- SVD of G via symmetric eigendecomp (G = V S^2 V^T)
- Automatic rank truncation: energy threshold or user-specified rank
- Projected operator Atilde via 2 dgemm calls
- Eigenvalues encode growth rate sigma = log|mu|/dt and frequency
  St = arg(mu)/(2*pi*dt)

**SPOD** — Spectral POD (Towne et al., 2018), two implementations (the modal
driver calls the streaming one; the batch routine gives the same eigenvalues):
- `modal_spod.f90` (batch): divide snapshots into overlapping blocks,
  window + FFT each block, form cross-spectral density (CSD) matrix
  at each frequency via `k_gram_complex` (single gop per frequency),
  eigensolve CSD for SPOD modes and eigenvalues
- `modal_spod_streaming.f90` (streaming): processes snapshots
  sequentially, accumulating DFT coefficients on the fly.  Same CSD
  eigensolve but with sequential access pattern for better cache
  locality on large datasets.  State management module tracks
  circular buffers, running mean, and DFT accumulators.

Both SPOD variants support Hamming windows with selectable
normalization (amplitude-preserving or energy-preserving / Parseval)
controlled by the `ifwinamp` flag.

#### Floquet energy budget

The perturbation kinetic energy budget for time-periodic base flows
(mode 4.11 / `energy_budget_floquet`) computes phase-averaged
production and dissipation over one orbit period.  Base flow is
re-loaded from disk at each orbit step to avoid storing the full
time-periodic base flow in memory.

#### Dynamic Mode Tracking (rc3)

Base-flow stabilization by band-pass filtering the velocity field
(Queguineur et al., Phys. Fluids 31, 034101, 2019), ported from the
legacy code into `dmt.f90`.  Target frequency, filter width, gain,
start time, and tolerance map to `uparam(11..15)`.  The force acts on the
band around the target frequency and vanishes for a steady flow, so the
run converges to a steady base flow: `cylinder_re100/130_baseflow_dmt`
(Re = 100) reaches a residual of 1e-7 and agrees with the Newton solution
to 4e-5.  A periodic orbit has content in the band at its own frequencies,
the force does not vanish on it, and the residual does not reach zero
(a trial on the forced jet orbit stalled near 0.7).  The filter bank of
mode 6.4 is the tool for the content of an orbit at given frequencies.

#### FST inflow module (rc3)

The free-stream turbulence inflow generation previously embedded in the
slot case is now a reusable module (`fst.f90`): any case enables it
with `use nekstab_fst` and a `call fst` in the inflow hook.  Mode data
is precomputed once and stored in a self-describing binary file, so the
expensive generation step is decoupled from the run.

#### Reproducibility framework (rc3)

Each validated example stage carries a `ref/` directory holding
`reference.json` — machine-checkable quantities with a source pointer
(file, column, row) and an absolute tolerance — plus the reference
figures produced by the stage's own `plot.py`.

```bash
python3 scripts/check_against_ref.py example/<family>/<stage>
```

re-extracts the quantities from the run outputs and reports PASS/FAIL
per quantity.  Supporting tooling under `scripts/`: residual and
spectrum parsers, a Slurm submit/check driver, field/IC inspection
utilities, and `make_case_ref.py` to classify case files and build
`ref/` skeletons.

#### Validation gallery (rc3)

The validation evidence is published as a static gallery
(`validation/index.html` + figures) deployed to GitHub Pages at
https://nekstab.github.io/nekStab/validation/. `validation/build_gallery.py`
builds it from the stage folders of `example/` (figures of `ref/`, quantities
of `ref/reference.json`); `validation/tests/test_gallery_complete.py` and a
Playwright test suite check it.

### String-based mode selection (new in 2.0)

In the `.usr` file, set `nekstab_mode` instead of encoding `uparam(1)`:
```fortran
subroutine nekStab_usrchk
   nekstab_mode = 'direct'          ! or 'adjoint', 'sfd', 'newton_fp', etc.
end subroutine
```
Priority: string > if-flags > uparam(1).  All three methods remain supported.

### Bug fixes

- Energy budget normalization: `alpha = 1/sqrt(...)` (was `sqrt(...)`)
- 2D energy budget: skip z-components (correctness — vz can contain
  uninitialized data in 2D, not guaranteed zero)
- SFD: velocity-only forcing, time-lag fix, file-open safety
- TDF: circular buffer for orbit storage, semicolon bug in single-line IF
- Module vs common-block variable collisions (OTD)
- Intel `ifx` compatibility (assumed-size arrays)
- BoostConv: QR re-orthogonalization restored to modified Gram-Schmidt
  (rc3; had regressed to classical GS, stalling convergence)
- POD/DMD/SPOD and mean-field outposts write the scalar field again
  (nfldt, rc3)
- Drag computed directly on wall (`W`) faces with the correct scaling
  (rc3)
- Floquet mode animation writes the correct deformed output (rc3)

### Robustness and diagnostics (rc3)

- The three orbit allocations (Floquet, TDF, Newton-UPO) — the largest
  arrays in the code — report the requested size in GB and exit cleanly
  on failure instead of hanging rank-asymmetrically on OOM
- Time-step guard in the sensitivity I/O setup and a division guard in
  the Floquet energy budget (same idiom as the existing matvec guard)
- Per-step solver logging inside the linearized maps throttled to ~10
  lines per matvec, and Nek core residual chatter silenced during
  matvecs; orbit-replay logfiles shrink by two orders of magnitude
  (flip-flop adjoint Floquet: 905k lines -> 2k)

### Build system

- `makefile_nekStab` with level-based dependency ordering
- `makeneks` archives all `.o` into `libnek5000.a`
- `nekstab_nek_bridge.f90`, `main.f90` and `usr_wrappers.f90` compile in fixed
  form (`-ffixed-form -ffixed-line-length-none`); every other `.f90` file
  compiles in free form
- Non-interactive builds survive configuration changes (rc3): `makeneks`
  pre-handles the makenek `.state` mismatch and answers the rebuild
  prompt itself, so Slurm/CI shells no longer die at an interactive
  `read` with a misleading "see build.log" message
- Include dependency edges (rc3): `NEKSTAB.inc` and `SIZE` are listed as
  prerequisites of the objects that include them, so incremental builds
  rebuild stale objects instead of relying on a full clean

### Examples (numbered stage layout, rc3)

The example tree was reorganized (rc3) as `example/<family>_<parameter>/`
with numbered stages encoding the workflow order, so each family reads
as a tutorial: produce a seed, converge a base flow, analyze it,
post-process it.

| Stage | Meaning |
|---|---|
| `000_dns`, `010_dns_forced` | DNS / seed generation, DNS with a forcing |
| `1x0_baseflow_*` | base flow by filtering (110 SFD, 120 BoostConv, 130 DMT, 140 TDF or mean flow) |
| `210_baseflow_newton` | base flow / UPO by Newton-GMRES |
| `3x0/3x1_stability_*` | 310 direct, 311 direct Floquet, 315 direct about the mean flow, 320 adjoint, 321 adjoint Floquet, 330 transient growth |
| `350_dns_unstable_orbit` | DNS that leaves an unstable orbit along its Floquet mode |
| `4xx_postproc_*` | 410 mode animation, 411 wavemaker or energy budget |
| `500_otd` | OTD modes |
| `6x0_modal_*` | 600 POD, 610 DMD, 620 SPOD, 630 DMT filter bank |

Families: `cylinder_re100` (flagship, full ladder), `cylinder_re180`
(Floquet chain; 2D proxy — Mode A/B physics needs 3D), `cylinder_re30_thermal`,
`cylinder_re1m` (2D RANS, experimental), `moving_cylinder_re100`, `flip_flop_re62`,
`thermosyphon_ra500`, `naca0012_re2000`, `lid_driven_re3600`,
`back_fstep_re500`, `cubic_cavity_re1914`/`_re1920`/`_re1950`, `tpjet_re2000`,
`slot_fst_re495`, `poiseuille_re5k` (OTD).

All cylinder isothermal cases share the same 2128-element mesh and a
shared seed initial condition for direct cross-case comparison.  The RANS
case (`cylinder_re1m`) uses the finite-difference Fréchet operator
only — no adjoint, wavemaker, or transient growth, since the analytical
RANS Jacobian is unavailable.  It is marked experimental: it has no
reference data and the gallery does not list it.

Validated stages carry `ref/` reference data checked by
`scripts/check_against_ref.py` (see Reproducibility framework above).
78 of the 79 stages outside `cylinder_re1m` ran again for 2.0.0 from a
tree that holds tracked files only: 52 stages on Jean Zay (4 to 40
ranks), 16 on c8 (8 to 16 ranks) and 10 on the Intel servers x1 and x3
(4 to 8 ranks).  The record of a run holds the digest of the source files
that built it.  Eight stages ran on the source files of commit `6c122882`.
The other 70 ran on earlier sources of the 2.0 series.  Each change of the
source after a run was read against the code path of the stages it can
reach (the OTD, DMT, BoostConv, Floquet energy budget and moving-mesh
changes, for example, reach only their own modes).  Commits `0ca2a0e1`
and `d65823c3` then changed the Krylov-Schur eigensolver (stop rule,
convergence test, check of each exported vector) and the targets of four
stages.  The 23 stages that run it (every stage numbered 310, 311, 315,
320, 321 or 330) ran before that change; their rerun is in the TODO list
below.  The 78 results passed their reference check when they ran and
took 1.4e3 core-hours (ranks times wall time).
One stage and the rerun of those 23 are open at this commit and listed in
the TODO list below.


## Full mode reference (uparam encoding)

### Mode 0 — DNS
- **0.0**: DNS (default if nothing set)
- **0.1**: Linearized DNS

### Mode 1 — Base flow computation
- **1.1**: Selective Frequency Damping (SFD)
- **1.2**: BoostConv acceleration
- **1.3**: Dynamic Mode Tracking (DMT)
- **1.4**: Time-Delayed Feedback (TDF)

### Mode 2 — Newton-Krylov
- **2.0**: Newton-GMRES for fixed points
- **2.1**: Newton-GMRES for periodic orbits (UPO)
- **2.2**: Newton-GMRES for forced UPOs

### Mode 3 — Stability analysis (Krylov-Schur)
- **3.1**: Direct LNSE
- **3.11**: Direct LNSE — Floquet
- **3.2**: Adjoint LNSE
- **3.21**: Adjoint LNSE — Floquet
- **3.3**: Transient growth (direct-adjoint)
- **3.31**: Transient growth — Floquet

### Mode 4 — Postprocessing
- **4.0**: Energy budget + wavemaker + BF sensitivity (all)
- **4.1**: Perturbation kinetic energy budget
- **4.11**: PKE budget — Floquet
- **4.2**: Wavemaker
- **4.3**: Sensitivity to base flow modifications
- **4.41**: Sensitivity to steady force (real part); not validated in 2.0.0, see `DOC.md`
- **4.42**: Sensitivity to steady force (imaginary part); not validated in 2.0.0, see `DOC.md`
- **4.43**: Delta forcing
- **4.50**: Animate eigenmode
- **4.51**: Animate base flow deformation
- **4.52**: Animate Floquet mode

### Mode 5 — OTD
- **5.0**: Optimally Time-Dependent modes, frozen base flow
- **5.1**: Optimally Time-Dependent modes, base flow advanced by the DNS

### Mode 6 — Modal analysis
- **6.0**: POD + DMD + SPOD (all)
- **6.1**: POD
- **6.2**: DMD
- **6.3**: SPOD
- **6.4**: DMT filter bank (online, on the DNS)


## TODO before tagging 2.0

- [x] Run validation suite across examples — 24 validated at rc2, 15 deferred
- [x] Tag `v2.0.0-rc2` and run rc cycle
- [x] Adjoint correctness: Floquet orbit replay + thermal adjoint — fixed
      and validated (rc3)
- [x] Tag `v2.0.0-rc3`
- [x] `ref/` reference data for every stage outside the experimental
      `cylinder_re1m` (cubic cavity, moving cylinder, thermal cylinder,
      wavemakers, forced jet)
- [x] `check_against_ref.py` sweep across the stages of the final rerun
      as the regression gate (78 of 79 stages had a passing result before
      the eigensolver change of `0ca2a0e1`)
- [x] Six stages ran again because their folders lost tracked log and
      output files (150 figure copies, four `plot.out` logs and twelve
      wavemaker output fields were removed): `back_fstep_re500/210_baseflow_newton`,
      `flip_flop_re62/210_baseflow_newton`,
      `thermosyphon_ra500/210_baseflow_newton`,
      `thermosyphon_ra500/310_stability_direct/direct`,
      `cylinder_re100/411_postproc_wavemaker` and
      `naca0012_re2000/411_postproc_wavemaker` (1 s to 20 min each, on
      Jean Zay; all six pass)
- [x] Results of the source of commit `6c122882` for
      `back_fstep_re500/330_transient_growth` (6.3 h on 8 ranks of c8) and
      `cylinder_re100/210_baseflow_newton/fp` (55 min on 8 ranks of c8)
- [ ] Rerun of the 23 stages that run the Krylov-Schur eigensolver (every
      stage numbered 310, 311, 315, 320, 321 or 330, about 790 core-hours)
      on the source of `d65823c3`; `lid_driven_re3600` 310 and 320 and
      `tpjet_re2000` 311 and 321 also have changed inputs
- [ ] `cubic_cavity_re1950/350_dns_unstable_orbit` (13.7 h on 16 ranks of c8,
      220 core-hours) keeps the reference of an earlier c8 run and runs
      again after the tag
- [x] Update `DOC.md` with string-mode documentation
- [x] Update copyright year (2020-2026) in `main.f90`
- [x] Review SPOD streaming module: it ran on `cylinder_re100/620_modal_spod`
      (8 ranks, Jean Zay), after a fix of the implicit real `spod_s_nfreq` in
      its external routines (they now import it and declare `implicit none`). Its eigenvalues equal those of the batch
      routine `spod_compute` to six digits in all 33 frequency bins, and the
      three eigenvalues of each bin sum to the total power of the POD-FFT
      analysis to four digits. The driver calls the streaming routine only;
      `spod_compute` is kept as the reference implementation.
- [ ] Tag `v2.0.0` and merge `dev` → `main`


---
