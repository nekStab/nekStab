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

The adaptive tolerance follows the **Eisenstat-Walker type 2** rule
(Eisenstat & Walker, SIAM J. Sci. Comput. 17(1), 1996):

    eta_k = (||F_k|| / ||F_{k-1}||)^alpha,   alpha = 1.618

    GMRES tolerance = eta_k^2 * ||F_k||^2

where eta_k is the *forcing term* — large (loose) when the residual
is still high, small (tight) as Newton converges.  The exponent
alpha = 1.618 (golden ratio) gives a good balance between quadratic
convergence rate and computational cost.  The forcing term is capped
at eta_max = 0.9 to prevent the linear solve from being too loose.

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

### New features

| Mode | uparam(1) | String | Description |
|------|-----------|--------|-------------|
| Linearized DNS | 0.1 | `linear_dns` | Time-stepping of LNSE without eigensolver |
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

### Build system

- `makefile_nekStab` with level-based dependency ordering
- `makeneks` archives all `.o` into `libnek5000.a`
- `.f90` compiled with `-ffixed-form -ffixed-line-length-none`

### Examples

~45 documented examples covering: cylinder (Re 50), backward-facing
step, flip-flop, cubic cavity (UPO), thermosyphon, turbulent pulsed
jet, NACA 0012, Blasius, torus, lid-driven cavity, Poiseuille (RANS,
OTD), slot jet with FST.

Validation script: `validate.py` (automated regression tests).


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

- [ ] Run full validation suite across all examples
- [ ] Finalize `validate.py` regression tests
- [ ] Update `DOC.md` with string-mode documentation
- [ ] Update copyright year (2020-2026) in `main.f90`
- [ ] Review SPOD streaming module status
- [ ] Verify all example READMEs are current
- [ ] Tag and push: `git tag -a 2.0 -m "nekStab 2.0"`


---
