# Release 2.0 (2026)

Complete architectural rewrite.  User-facing API (uparam encoding) is
backward compatible with 1.x; internal code is restructured for
maintainability, numerical robustness, and HPC performance.


## What changed

### Architecture

- All source wrapped in Fortran modules (33 `.f90` files, was loose `.f`)
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
  .true.`).  Two RANS examples are included:
  - `example/poiseuille_RANS/` — turbulent channel at Re = 100,000
  - `example/cylinder/RANS/` — cylinder at Re = 40,000
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
| Dynamic Mode Tracking | 1.3 | `dmt` | UPO stabilization by band-pass filtering (Queguineur et al., 2019) |
| Floquet energy budget | 4.11 | `energy_budget_floquet` | PKE budget for time-periodic base flows |
| Animate mode | 4.50 | `animate_mode` | Reconstruct eigenmode as time series |
| Animate BF deform | 4.51 | `animate_bf_deform` | Base flow + eigenmode deformation |
| Animate Floquet | 4.52 | `animate_floquet` | Floquet mode animation |
| POD | 6.1 | `pod` | Proper Orthogonal Decomposition |
| DMD | 6.2 | `dmd` | Dynamic Mode Decomposition |
| SPOD | 6.3 | `spod` | Spectral POD (Welch + CSD) |
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

**SPOD** — Spectral POD (Towne et al., 2018), two implementations:
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

Base-flow stabilization of unstable periodic orbits by band-pass
filtering the velocity field (Queguineur et al., Phys. Fluids 31,
034101, 2019), ported from the legacy code into `dmt.f90`.  Target
frequency, filter width, gain, start time, and tolerance map to
`uparam(11..15)`.

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
https://nekstab.github.io/nekStab/validation/, catalog-driven and
checked by a Playwright test suite.

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
- `.f90` compiled with `-ffixed-form -ffixed-line-length-none`
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
| `000_dns` | DNS / seed generation |
| `1x0_baseflow_*` | base flow by filtering (110 SFD, 120 BoostConv, 130 DMT) |
| `210_baseflow_newton` | base flow / UPO by Newton-GMRES |
| `3x0/3x1_stability_*` | 310 direct, 311 direct Floquet, 320 adjoint, 321 adjoint Floquet, 330 transient growth |
| `4xx_postproc_*` | 410 mode animation, 411 wavemaker, 412 steady-force sensitivity |
| `500_otd` | OTD modes |
| `6x0_modal_*` | 600 POD, 610 DMD, 620 SPOD |

Families: `cylinder_re100` (flagship, full ladder), `cylinder_re180`
(Floquet chain; 2D proxy — Mode A/B physics needs 3D), `cylinder_re30_thermal`,
`cylinder_re1m` (2D RANS), `moving_cylinder_re100`, `flip_flop_re62`,
`thermosyphon_ra500`, `naca0012_re2000`, `lid_driven_re3600`,
`back_fstep_re500`, `cubic_cavity_re1914`/`_re1950`, `tpjet_re2005`,
`slot_fst_re495`, `poiseuille_re5k` (OTD), `poiseuille_re1e5` (RANS).

All cylinder isothermal cases share the same 2128-element mesh and a
shared seed initial condition for direct cross-case comparison.  RANS
cases (`cylinder_re1m`, `poiseuille_re1e5`) use the finite-difference
Fréchet operator only — no adjoint, wavemaker, or transient growth,
since the analytical RANS Jacobian is unavailable.

Validated stages carry `ref/` reference data checked by
`scripts/check_against_ref.py` (see Reproducibility framework above);
the rc2 campaign validated 24 cases on the previous flat layout
(8-rank Slurm runs, 2026-05-17), and those results carry over to the
renumbered stages.  Stages without `ref/` are tracked configurations
whose reference runs are still pending.


## Full mode reference (uparam encoding)

### Mode 0 — DNS
- **0.0**: DNS (default if nothing set)
- **0.1**: Linearized DNS

### Mode 1 — Base flow computation
- **1.1**: Selective Frequency Damping (SFD)
- **1.2**: BoostConv acceleration
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
- **4.41**: Sensitivity to steady force (real part)
- **4.42**: Sensitivity to steady force (imaginary part)
- **4.43**: Delta forcing
- **4.50**: Animate eigenmode
- **4.51**: Animate base flow deformation
- **4.52**: Animate Floquet mode

### Mode 5 — OTD
- **5.0**: Optimally Time-Dependent modes

### Mode 6 — Modal analysis
- **6.0**: POD + DMD + SPOD (all)
- **6.1**: POD
- **6.2**: DMD
- **6.3**: SPOD


## TODO before tagging 2.0

- [x] Run validation suite across examples — 24 validated at rc2, 15 deferred
- [x] Tag `v2.0.0-rc2` and run rc cycle
- [x] Adjoint correctness: Floquet orbit replay + thermal adjoint — fixed
      and validated (rc3)
- [x] Tag `v2.0.0-rc3`
- [ ] Finish the rerun campaign: `ref/` reference data for every remaining
      stage (RANS findiff, cubic cavity, moving cylinder, thermal cylinder,
      wavemakers)
- [ ] `check_against_ref.py` sweep across all migrated stages as the
      regression gate (replaces the old `validate.py` plan)
- [ ] Update `DOC.md` with string-mode documentation — partial
- [ ] Update copyright year (2020-2026) in `main.f90`
- [ ] Review SPOD streaming module status — deferred to v2.1
- [ ] Tag `v2.0.0` and merge `dev` → `main`


---
