# nekStab Documentation

A toolbox for global stability and bifurcation analysis using the spectral element solver Nek5000.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Compilation](#compilation)
4. [Quick Start](#quick-start)
5. [Operating Modes](#operating-modes)
6. [Parameter Reference](#parameter-reference)
7. [Code Architecture](#code-architecture)
8. [Module Reference](#module-reference)
9. [Vector Operations](#vector-operations)
10. [Mesh Generation](#mesh-generation)
11. [Examples](#examples)
12. [Validation](#validation)
13. [Continuous Integration](#continuous-integration)
14. [Theoretical Background](#theoretical-background)
15. [Troubleshooting](#troubleshooting)
16. [Citation](#citation)
17. [License](#license)
18. [Contact](#contact)

---

## Overview

**nekStab** extends Nek5000 with capabilities for:

- **Direct Numerical Simulation (DNS)** - Standard and linearized time integration
- **Base Flow Computation** - Steady/periodic states via SFD, BoostConv, TDF, Newton-Krylov
- **Linear Stability Analysis** - Global eigenvalue problems using Krylov-Schur
- **Sensitivity Analysis** - Structural sensitivity, wavemaker, optimal forcing
- **Floquet Analysis** - Stability of time-periodic base flows
- **Transient Growth** - Optimal perturbation and non-modal analysis
- **OTD Modes** - Optimally Time-Dependent modes for chaotic flows

### Why nekStab?

Traditional stability analysis requires forming and storing large Jacobian matrices. For 3D spectral element discretizations with millions of degrees of freedom, this is impractical. nekStab uses **matrix-free time-stepping** methods:

1. The linearized Navier-Stokes operator **A** is never formed explicitly
2. Matrix-vector products **Aq** are computed by time-stepping the linearized equations
3. Krylov subspace methods extract eigenvalues from these products
4. Memory scales with the number of Krylov vectors, not the matrix size

---

## Installation

### Prerequisites

**Linux (Ubuntu/Debian) - GCC**
```bash
sudo apt install build-essential gfortran libopenmpi-dev liblapack-dev libblas-dev libfftw3-dev cmake
```

**macOS**
```bash
brew install gcc open-mpi fftw cmake
```

**HPC Systems**
```bash
module load gcc openmpi  # or intel-oneapi
```

### Intel oneAPI Installation (Recommended)

For optimal performance on Intel/AMD x86 CPUs, install Intel oneAPI. The `ifx` compiler with MKL typically provides 20-40% faster execution than GCC.

#### Repository Setup

```bash
# Install prerequisites
sudo apt update
sudo apt install -y gpg-agent wget

# Add Intel GPG key
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

# Add Intel repository
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt update
```

#### Installation Options

| Package | Command | Disk Space |
|---------|---------|------------|
| Fortran compiler only | `sudo apt install intel-oneapi-compiler-fortran` | ~3 GB |
| Fortran + MKL | `sudo apt install intel-oneapi-compiler-fortran intel-oneapi-mkl` | ~8 GB |
| Full HPC Toolkit | `sudo apt install intel-oneapi-hpc-toolkit` | ~15 GB |

The **Fortran + MKL** option is recommended for nekStab - it provides the compiler and optimized math libraries without unnecessary components.

#### Environment Setup

Add to `~/.bashrc` or `~/.zshrc`:
```bash
# Intel oneAPI environment
if [ -f /opt/intel/oneapi/setvars.sh ]; then
    source /opt/intel/oneapi/setvars.sh > /dev/null
fi
```

Verify installation:
```bash
source ~/.bashrc
ifx --version
# Intel(R) Fortran Compiler for oneAPI 2025.x.x
```

#### What's Included in oneAPI 2025

- **ifx** - Intel Fortran Compiler (LLVM-based, Fortran 2023 features)
- **MKL** - Math Kernel Library (optimized BLAS, LAPACK, FFT)
- **MPI** - Intel MPI Library (in HPC toolkit)
- **OpenMP 6.0** - Parallel programming support

> **Note:** Intel discontinued `ifort` in oneAPI 2025, but it remains available on many HPC systems with older installations. nekStab supports both `ifx` and `ifort`.

#### FFT Library Linking

nekStab uses FFTW3 for spectral analysis (temporal FFT, future SPOD support):

| Compiler | FFT Backend | Linking |
|----------|-------------|---------|
| **Intel (ifx/ifort)** | MKL FFTW wrappers | Automatic via `-qmkl` |
| **GCC (gfortran)** | FFTW3 | Requires `libfftw3-dev` |

Intel's MKL includes FFTW3-compatible wrappers, so no additional installation is needed. The `-qmkl` flag (already used for BLAS/LAPACK) provides optimized FFT routines. GCC builds link against the system FFTW3 library.

### Clone and Setup

```bash
git clone https://github.com/nekStab/nekStab.git
cd nekStab
./Nek5000setup.sh
```

The setup script will:
1. Clone Nek5000 into the `Nek5000/` subdirectory
2. Build the mesh tools (`genmap`, `genbox`)
3. Optionally configure environment variables

### Environment Variables

Add to `~/.bashrc` or `~/.zshrc`:

```bash
export NEKSTAB_SOURCE_ROOT=$HOME/nekStab
export NEK_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT/Nek5000
export PATH=$NEK_SOURCE_ROOT/bin:$NEKSTAB_SOURCE_ROOT/bin:$PATH
ulimit -s unlimited   # Required for large stack allocations
ulimit -c unlimited   # Enable core dumps for debugging
```

---

## Compilation

### Basic Usage

```bash
cd example/cylinder_re100/000_dns
mks 1cyl
```

The `mks` script:
1. Detects available compilers
2. Sets appropriate optimization flags
3. Copies `src/NEKSTAB` into the case-local `NEKSTAB.inc`
4. Generates `makefile_usr.inc` that includes `src/makefile_nekStab`
5. Invokes Nek5000's build system through `bin/makeneks`
6. Compiles with parallel make (`-j4`)

### Shell Completion

Enable tab completion for case names and common `mks` targets/options:

```bash
source "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash"
```

For zsh, enable Bash completion first:

```bash
autoload -Uz bashcompinit && bashcompinit
source "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash"
```

### Compiler Selection

Auto-detection order: `ifx` → `ifort` → `gfortran`

Force a specific compiler:
```bash
NEKSTAB_FC=ifx mks 1cyl    # Intel LLVM (recommended)
NEKSTAB_FC=ifort mks 1cyl  # Intel Classic (HPC systems)
NEKSTAB_FC=gcc mks 1cyl    # GCC/gfortran
```

### Compiler Comparison

| Compiler | BLAS/LAPACK | FFT | Vectorization | Notes |
|----------|-------------|-----|---------------|-------|
| **ifx** | MKL | MKL FFTW | AVX2/AVX-512 | Requires `source /opt/intel/oneapi/setvars.sh` |
| **ifort** | MKL | MKL FFTW | AVX2/AVX-512 | Discontinued in oneAPI 2025, but available on many HPC systems |
| **gfortran** | System libs | FFTW3 | Native | Uses `-framework Accelerate` on macOS, requires `libfftw3-dev` on Linux |

> **Note:** Intel discontinued `ifort` in oneAPI 2025, but it remains available on many HPC systems with older oneAPI or Intel Parallel Studio installations.

MKL provides runtime CPU dispatching - the same binary runs optimally on different Intel/AMD processors.

### Build Options

```bash
mks 1cyl --debug   # Enable: -g -fbacktrace -Wall (GCC) or -g3 -traceback (Intel)
mks 1cyl --fresh   # Clean build: removes obj/, *.mod, makefile, logs
```

### Custom Flags

```bash
NEKSTAB_EXTRA_FFLAGS="-O3 -march=native" mks 1cyl
```

---

## Quick Start

### 1. DNS Simulation

```bash
cd example/cylinder_re100/000_dns
mks 1cyl                    # Compile
nekbmpi 1cyl 4              # Run on 4 MPI ranks
tail -f logfile             # Monitor output
killall nek5000             # Stop if needed
```

### 2. Base Flow (Newton-Krylov)

```bash
cd example/cylinder_re100/210_baseflow_newton
# Ensure userParam01 = 2.0 in 1cyl.par
mks 1cyl && nekbmpi 1cyl 4
# Output: BF_1cyl0.f00001 (converged base flow)
```

### 3. Direct Stability Analysis

```bash
cd example/cylinder_re100/310_stability_direct
# Ensure userParam01 = 3.1 in 1cyl.par
# Ensure base flow file exists: BF_1cyl0.f00001
mks 1cyl && nekbmpi 1cyl 4
# Output: dRe*, dIm* (eigenmodes), eigenvalues in logfile
```

### 4. Adjoint Analysis

```bash
cd example/cylinder_re100/320_stability_adjoint
# Ensure userParam01 = 3.2 in 1cyl.par
mks 1cyl && nekbmpi 1cyl 4
# Output: aRe*, aIm* (adjoint eigenmodes)
```

### 5. Wavemaker Computation

```bash
cd example/cylinder_re100/410_postproc_animate_modes
# Ensure userParam01 = 4.2 in 1cyl.par
# Requires both direct and adjoint modes
mks 1cyl && nekbmpi 1cyl 4
```

---

## Operating Modes

nekStab supports **three equivalent methods** to select the operating mode, providing flexibility for different use cases.

### Mode Configuration Interface

#### Method 1: String Mode (Recommended)

The most readable approach. Set `nekstab_mode` in your `.usr` file's `nekStab_usrchk` subroutine:

```fortran
subroutine nekStab_usrchk
   nekstab_mode = 'floquet_adjoint'  ! Clear, self-documenting
   k_dim = 92                        ! Krylov dimension
end subroutine
```

**Available mode strings:**

| String | Description | Equivalent uparam |
|--------|-------------|-------------------|
| `'dns'` | Direct Numerical Simulation | 0 |
| `'linear_dns'` | Linearized DNS | 0.1 |
| `'sfd'` | Selective Frequency Damping | 1.1 |
| `'boostconv'` | BoostConv acceleration | 1.2 |
| `'tdf'` | Time-Delayed Feedback | 1.4 |
| `'newton_fp'` | Newton for fixed points | 2.0 |
| `'newton_po'` | Newton for periodic orbits | 2.1 |
| `'newton_po_t'` | Newton for forced periodic orbits | 2.2 |
| `'direct'` | Direct eigenmode analysis | 3.1 |
| `'floquet_direct'` | Direct Floquet analysis | 3.11 |
| `'adjoint'` | Adjoint eigenmode analysis | 3.2 |
| `'floquet_adjoint'` | Adjoint Floquet analysis | 3.21 |
| `'transient_growth'` | Transient growth analysis | 3.3 |
| `'floquet_tg'` | Floquet transient growth | 3.31 |
| `'energy_budget'` | Kinetic energy budget | 4.1 |
| `'energy_budget_floquet'` | Floquet energy budget | 4.11 |
| `'wavemaker'` | Wavemaker/structural sensitivity | 4.2 |
| `'bf_sensitivity'` | Base flow sensitivity | 4.3 |
| `'force_sensitivity_real'` | Real part forcing sensitivity | 4.41 |
| `'force_sensitivity_imag'` | Imaginary part forcing sensitivity | 4.42 |
| `'delta_forcing'` | Delta forcing response | 4.43 |
| `'animate_mode'` | Animate eigenmode | 4.50 |
| `'animate_bf_deform'` | Animate with base flow deformation | 4.51 |
| `'animate_floquet'` | Animate Floquet mode | 4.52 |
| `'otd'` | Optimally Time-Dependent modes | 5 |
| `'pod'` | Proper Orthogonal Decomposition | 6.1 |
| `'dmd'` | Dynamic Mode Decomposition | 6.2 |
| `'spod'` | Spectral POD | 6.3 |

Mode strings are **case-insensitive** and support aliases:

| Alias | Equivalent |
|-------|------------|
| `'lindns'`, `'linearized_dns'` | `'linear_dns'` |
| `'boost'` | `'boostconv'` |
| `'newton'` | `'newton_fp'` |
| `'upo'` | `'newton_po'` |
| `'forced_upo'` | `'newton_po_t'` |
| `'floquetdirect'` | `'floquet_direct'` |
| `'floquetadjoint'` | `'floquet_adjoint'` |
| `'tg'` | `'transient_growth'` |
| `'floquet_transient_growth'` | `'floquet_tg'` |
| `'baseflow_sensitivity'` | `'bf_sensitivity'` |
| `'force_sens_real'` | `'force_sensitivity_real'` |
| `'force_sens_imag'` | `'force_sensitivity_imag'` |
| `'animate'` | `'animate_mode'` |
| `'animate_deform'` | `'animate_bf_deform'` |

#### Method 2: Flag Mode (Flexible)

Set individual boolean flags for fine-grained control:

```fortran
subroutine nekStab_usrchk
   isAdjoint = .true.    ! Base mode
   ifFloquet = .true.    ! Modifier flag
   k_dim = 92
end subroutine
```

**Mode flags by category:**

| Category | Flags |
|----------|-------|
| DNS | `ifDNS`, `ifLinDNS` |
| Fixed Point | `ifSFD`, `ifBoostConv`, `ifTDF` |
| Newton-Krylov | `isNewtonFP`, `isNewtonPO`, `isNewtonPO_T` |
| Stability | `isDirect`, `isAdjoint`, `isTransientGrowth` |
| Floquet modifier | `ifFloquet` (transforms base stability modes) |
| Postprocessing | `ifEnergyBudget`, `ifWavemaker`, `ifBFSensitivity`, `ifForceSensReal`, `ifForceSensImag`, `ifDeltaForcing`, `ifAnimateMode`, `ifAnimateBFDeform`, `ifAnimateFloquet` |
| OTD | `ifotd` |
| Modal Analysis | `ifpod`, `ifdmd`, `ifspod` |

**The `ifFloquet` modifier** combines with base stability modes:
- `isDirect = .true.` + `ifFloquet = .true.` → Floquet direct analysis
- `isAdjoint = .true.` + `ifFloquet = .true.` → Floquet adjoint analysis
- `isTransientGrowth = .true.` + `ifFloquet = .true.` → Floquet transient growth

**Animation parameters:**
- `animate_mode_num = 3` — specifies which mode number to animate (default: 1)

#### Method 3: uparam (Automation-Friendly)

The original method, still fully supported. Set in your `.par` file:

```ini
[PROBLEMTYPE]
userParam01 = 3.21    # Floquet adjoint
userParam07 = 3       # Mode number for animation
```

This method is ideal for **parametric studies** and **automation scripts** where you need to change modes without recompiling.

#### Priority Order

When multiple methods are used, the priority is:

1. **String mode** (`nekstab_mode`) — highest priority
2. **Flag mode** (any `if*` flag set) — middle priority
3. **uparam(1)** — lowest priority (fallback)

This means user-set flags in `.usr` always override `.par` file settings, giving explicit user intent precedence.

#### Example: Three Ways to Run Floquet Adjoint

**Method 1: String (clearest)**
```fortran
! In .usr nekStab_usrchk
nekstab_mode = 'floquet_adjoint'
k_dim = 92
```

**Method 2: Flags (flexible)**
```fortran
! In .usr nekStab_usrchk
isAdjoint = .true.
ifFloquet = .true.
k_dim = 92
```

**Method 3: uparam (automation)**
```ini
# In .par file
userParam01 = 3.21
userParam07 = 92
```

All three produce **identical behavior**.

---

### Mode Reference

The following sections describe each mode in detail. The `userParam01` values are shown for reference but can be replaced with string or flag methods above.

### Mode 0: DNS

| Value | Description |
|-------|-------------|
| `0` | Standard DNS (nonlinear Navier-Stokes) |
| `0.1` | Linearized DNS (perturbation around base flow) |

### Mode 1: Fixed Point Methods

| Value | Method | Description |
|-------|--------|-------------|
| `1.1` | SFD | Selective Frequency Damping |
| `1.2` | BoostConv | Accelerating slow convergence |
| `1.3` | DMT | Dynamic Mode Tracking (not yet ported) |
| `1.4` | TDF | Time-Delayed Feedback |

**SFD Parameters** (`userParam01 = 1.1`):
- `userParam04`: Frequency `f` (positive: Åkervik formulation, negative: Casacuberta formulation)
- `userParam05`: Filter width `σ` (use `0` for continuation from previous run)

### Mode 2: Newton-Krylov

| Value | Description |
|-------|-------------|
| `2.0` | Newton for fixed points (steady states) |
| `2.1` | Newton for periodic orbits (`endTime` = initial guess) |
| `2.2` | Newton for forced periodic orbits (`endTime` = 1/f) |

Uses GMRES to solve the Newton system. The Jacobian-vector product is computed via finite differences of the time-stepper.
Arnoldi orthogonalization uses CGS2 by default (`use_cgs = .true.`), with classic MGS + reorthogonalization still available for debugging or comparison.

### Mode 3: Eigenvalue Problems

| Value | Description | Output Prefix |
|-------|-------------|---------------|
| `3.1` | Direct LNSE | `dRe*`, `dIm*` |
| `3.11` | Direct Floquet | `dRe*`, `dIm*` |
| `3.2` | Adjoint LNSE | `aRe*`, `aIm*` |
| `3.21` | Adjoint Floquet | `aRe*`, `aIm*` |
| `3.3` | Transient growth (optimal perturbation) | `oRe*`, `oIm*` |
| `3.31` | Transient growth Floquet | `oRe*`, `oIm*` |

**Key Parameter**: `userParam02`
- `0`: Start from random noise
- `> 0`: Restart from `userParam02` existing Krylov vectors

### Mode 4: Post-Processing

| Value | Description |
|-------|-------------|
| `4.0` | All: budget + wavemaker + BF sensitivity |
| `4.1` | Kinetic energy budget |
| `4.11` | Floquet kinetic energy budget |
| `4.2` | Wavemaker (structural sensitivity) |
| `4.3` | Base flow sensitivity |
| `4.41` | Real part forcing sensitivity |
| `4.42` | Imaginary part forcing sensitivity |
| `4.43` | Delta forcing response |
| `4.50` | Animate direct mode |
| `4.51` | Animate mode with base flow deformation |
| `4.52` | Animate Floquet mode |

### Mode 5: OTD

Optimally Time-Dependent modes for systems with chaotic base flows. Tracks the most unstable directions in real-time.

### Mode 6: Modal Analysis (POD, DMD, SPOD)

Data-driven modal decomposition methods for extracting coherent structures from time-resolved snapshot data.

| Value | Method | Description |
|-------|--------|-------------|
| `6.0` | All enabled | Run all methods where `ifpod`/`ifdmd`/`ifspod = .true.` |
| `6.1` | POD | Proper Orthogonal Decomposition only |
| `6.2` | DMD | Dynamic Mode Decomposition only |
| `6.3` | SPOD | Spectral POD only |

#### Overview

| Method | Basis | Output | Best For |
|--------|-------|--------|----------|
| **POD** | Energy-optimal spatial modes | Ranked by turbulent kinetic energy | Low-dimensional models, energy analysis |
| **DMD** | Temporal growth rate + frequency | Complex eigenvalues λ = σ + iω | Identifying oscillatory structures, linear dynamics |
| **SPOD** | Frequency-resolved energy-optimal modes | Eigenvalues per Strouhal number | Turbulent flows, broadband spectra |

#### Configuration

Set these parameters in `nekStab_usrchk`:

```fortran
subroutine nekStab_usrchk
   ! Required parameters
   modal_prefix = '   '       ! Blank = use SESSION stem directly (e.g. 1cyl0.f*)
   modal_nsnap  = 100         ! Number of snapshots to load
   modal_dt     = 0.5         ! Time step between snapshots
   modal_nsave  = 10          ! Number of modes to save to disk

   ! Method selection
   ifpod  = .true.            ! Enable POD
   ifdmd  = .true.            ! Enable DMD
   ifspod = .false.           ! Enable SPOD

   ! Optional: DMD rank truncation (0 = auto)
   dmd_rank = 0

   ! Optional: SPOD windowing
   spod_nfft = 64             ! FFT block size
   spod_noverlap = 32         ! Block overlap

   ! Optional: Window normalization for spectral analysis
   ifwinamp = .true.          ! .true.=PySPOD (amplitude), .false.=Parseval (energy)
end subroutine
```

#### Input Files

Snapshots must be Nek5000 field files with sequential numbering:
```
<prefix>0.f00001, <prefix>0.f00002, ..., <prefix>0.f<nsnap>
```

Example: `1cyl0.f00001` through `1cyl0.f00100` for 100 snapshots.

#### Output Files

| File | Contents | Format |
|------|----------|--------|
| `pod_energy.dat` | POD eigenvalue spectrum | `mode, eigenvalue, energy_%, cumulative_%` |
| `dmd_spectrum.dat` | DMD eigenvalues | `mode, |μ|, σ, ω, St, Re(μ), Im(μ)` |
| `spod_spectrum.dat` | SPOD spectrum | `St, λ₁, λ₂, ..., λ_nblk` |
| `mea*.f00001` | Temporal mean field | Nek5000 field file |
| `pod*.f00001` | POD modes | Nek5000 field files |
| `dm1*.f00001` | DMD modes (real part) | Nek5000 field files |
| `dm2*.f00001` | DMD modes (imag part) | Nek5000 field files |
| `SPG*.f00001` | SPOD modes | Nek5000 field files |

#### Plotting Results

Use the provided Python script:

```bash
cd example/cylinder_re100/600_modal_pod
python plot_modal.py           # Auto-detect and plot all
python plot_modal.py --pod     # POD only
python plot_modal.py --dmd     # DMD only
python plot_modal.py --spod    # SPOD only
```

Generates: `pod_spectrum.png`, `dmd_spectrum.png`, `spod_spectrum.png`

#### Method Details

**POD (Proper Orthogonal Decomposition)**

Computes energy-optimal modes via the method of snapshots:
1. Subtract temporal mean from snapshots
2. Form correlation matrix C_ij = ⟨q_i, q_j⟩
3. Solve eigenvalue problem: C v = λ v
4. Reconstruct spatial modes: φ_k = Σ_i v_ik q_i

Eigenvalues represent the kinetic energy captured by each mode. The first few modes typically capture >90% of fluctuation energy for periodic flows.

**DMD (Dynamic Mode Decomposition)**

Extracts modes with exponential temporal behavior using projected DMD:
1. Form data matrices: X = [q₁,...,q_{n-1}], Y = [q₂,...,q_n]
2. SVD of X: X = U Σ V*
3. Project dynamics: Ã = U* Y V Σ⁻¹
4. Eigendecomposition of Ã gives DMD eigenvalues μ

Eigenvalues relate to continuous-time growth rate σ and frequency ω:
- σ = log|μ| / Δt (growth rate)
- ω = arg(μ) / Δt (angular frequency)
- St = ω / (2π) (Strouhal number)

**SPOD (Spectral POD)**

Frequency-resolved POD using Welch's method:
1. Divide snapshots into overlapping blocks
2. Apply windowing (Hamming) and FFT each block
3. At each frequency, form cross-spectral density matrix
4. Eigendecomposition gives SPOD modes ranked by energy at that frequency

SPOD reduces to POD for a single block and to DMD for broadband-limited signals.

#### Example Workflow

```bash
# 1. Generate DNS snapshots (run DNS, save every 0.5 time units)
cd example/cylinder_re100/000_dns
# Edit 1cyl.par: writeInterval = 0.5, endTime = 50
mks 1cyl && nekbmpi 1cyl 4

# 2. Run modal analysis
cd ../modal
ln -s ../dns/1cyl0.f* .     # Link snapshots
# Edit 1cyl.usr nekStab_usrchk:
#   modal_prefix = '   '
#   modal_nsnap = 100
#   modal_dt = 0.5
mks 1cyl && nekbmpi 1cyl 4

# 3. Plot spectra
python plot_modal.py
```

---

## Parameter Reference

### .par File Structure

```ini
[GENERAL]
startFrom = BF_case0.f00001    # Initial condition file
stopAt = endTime
endTime = 100.0

userParam01 = 3.1    # Operating mode
userParam02 = 0      # Krylov restart (0 = fresh start)
userParam03 = ...    # Method-specific
userParam04 = ...    # SFD/TDF frequency
userParam05 = ...    # Forced frequency flag

dt = 0               # 0 = variable dt
variableDt = yes
targetCFL = 0.5

[VELOCITY]
viscosity = -100.0   # Negative = Reynolds number
```

### nekStab Parameters

Set in `nekStab_usrchk` subroutine in your `.usr` file:

#### Krylov Solver

| Parameter | Default | Description |
|-----------|---------|-------------|
| `k_dim` | 100 | Krylov subspace dimension (larger = more eigenvalues but more memory) |
| `use_cgs` | .true. | Use CGS2 orthogonalization with batched BLAS/MPI reductions (`.false.` = classic MGS + reorthogonalization) |
| `eigen_tol` | 1e-6 | Convergence tolerance for eigenvalues |
| `schur_tgt` | 2 | Number of eigenvalues to lock per restart |
| `schur_del` | 0.10 | Deflation threshold |
| `maxmodes` | 20 | Maximum eigenmodes to save to disk |

**Optimal Parameter Selection** (from parametric studies):

The product of Krylov dimension `m` and sampling time `τ` should satisfy:

```
4 < m × τ / T < 20
```

where `T` is the characteristic timescale of the instability (period for oscillatory modes, doubling time for stationary modes).

| Instability Type | Estimate T | Recommendation |
|-----------------|------------|----------------|
| Hopf (oscillatory) | T = 1/St | Use τ = T/8, m = 100-150 |
| Pitchfork (steady) | T ≈ 1 (diffusive) | Use τ = 1, m = 100-150 |
| Floquet | T = orbit period | Use τ = T, m = 64-128 |

**Example**: For cylinder wake at Re = 50 with St = 0.125:
- T = 1/0.125 = 8
- Choose τ = 1 (= T/8) and m = 100
- Product: m × τ / T = 100 × 1 / 8 = 12.5 ✓

**Convergence rule of thumb**: For open shear flows, eigenvalues converge after total integration time exceeds one flow-through time:
```
m × τ > L_x / U_∞
```
where L_x is the streamwise domain extent.

#### Newton-Krylov

| Parameter | Default | Description |
|-----------|---------|-------------|
| `findiff_order` | 2 | Finite difference order (1, 2, or 4) |
| `epsilon_base` | 1e-6 | Perturbation scale ε for Jacobian approximation |
| `ifdyntol` | .false. | Enable residual-proportional inner GMRES tolerances |
| `ew_tol_cap` | 0.0 | Optional upper cap for dynamic tolerance (`0` = uncapped) |

Newton iteration limits (30 Newton steps, 30 GMRES restarts) are hardcoded in `newton_krylov.f90`.

**Adaptive GMRES forcing**: When `ifdyntol = .true.`, the Newton solver sets the inner GMRES/Nek tolerance to 20% of the current squared nonlinear residual, bounded below by the user target `dtol = max(param(21), param(22))`. This keeps the inner solve tighter than the current Newton residual without spending a full extra decade of accuracy at every iteration. The original `param(21:22)` values are restored after Newton exits.

**Newton safeguards**:
- Warn after 3 consecutive iterations without meaningful residual decrease
- Abort if the residual exceeds `1e8` times the initial residual
- Use `ew_tol_cap` in `.usr` only if a stiff case needs to limit overly loose inner solves

**Performance Guidelines** (from parametric studies):

Similar to eigenvalue problems, optimal Newton-Krylov performance requires:
```
4 < m × τ / T < 20
```

| Case | Optimal m×τ | Time to Solution |
|------|-------------|------------------|
| Cylinder 2D (Re=80) | m×τ ≈ 67 | ~1 minute |
| Open cavity 2D (Re=4700) | m×τ ≈ 10 | ~30 seconds |

**Comparison with SFD**: For thermosyphon at Ra = 16,100:
- Newton-Krylov: 147 seconds
- SFD: 1071 seconds (7× slower)

Newton-Krylov is particularly advantageous for:
- Saddle-node fixed points (SFD cannot compute these)
- Unstable periodic orbits
- Cases far from bifurcation points

#### BoostConv

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bst_skp` | 10 | Iterations between acceleration updates |
| `bst_snp` | 10 | Size of residual subspace for extrapolation |

#### Output Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifres` | .false. | Save Krylov restart files (KRY*, HES*) |
| `ifvor` | .false. | Output vorticity components (vor*) |
| `ifvox` | .false. | Output vortex criteria Q, λ₂, ω (vox*) |
| `glob_skip` | 10 | Steps between energy/enstrophy output |

#### Initial Seed

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifseed_nois` | .true. | Initialize with random noise |
| `ifseed_symm` | .false. | Use symmetric initial perturbation |
| `ifseed_load` | .false. | Load seed from file |

If all are `.false.`, the `useric` subroutine defines the initial condition.

#### Base Flow

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifldbf` | .true. | Load base flow for linearized computations |
| `ifbf2D` | .false. | Force 2D base flow (zero w-component) |
| `ifstorebase` | .true. | Store base flow in memory for Floquet |

#### Sponge Zone

| Parameter | Default | Description |
|-----------|---------|-------------|
| `xLspg`, `xRspg` | 0.0 | Left/right sponge boundaries in x |
| `yLspg`, `yRspg` | 0.0 | Bottom/top sponge boundaries in y |
| `zLspg`, `zRspg` | 0.0 | Front/back sponge boundaries in z |
| `spng_st` | 0.0 | Sponge strength (0 = disabled) |
| `acc_spg` | 0.333 | Acceleration ramp fraction |

#### OTD

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifotd` | .false. | Enable OTD computation |
| `otd_printStep` | 100 | Output interval for OTD modes |
| `otd_gsStep` | 10 | Gram-Schmidt orthogonalization interval |
| `otd_FTLEPeriod` | 0.0 | FTLE averaging window |
| `otd_convTol` | 1e-6 | FTLE convergence tolerance |
| `otd_minSteps` | 200 | Minimum number of steps before FTLE convergence checks |

#### Modal Analysis (POD/DMD/SPOD)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `modal_prefix` | `'dns'` | 3-character snapshot file prefix; blank string uses the `SESSION` stem directly |
| `modal_nsnap` | 100 | Number of snapshots to load |
| `modal_dt` | 0.1 | Time step between snapshots |
| `modal_nsave` | 10 | Number of modes to save to disk |
| `ifpod` | .false. | Enable POD computation |
| `ifdmd` | .false. | Enable DMD computation |
| `ifspod` | .false. | Enable SPOD computation |
| `ifwinamp` | .true. | Window normalization: .true.=amplitude (PySPOD), .false.=energy (Parseval) |
| `dmd_rank` | 0 | DMD truncation rank (0 = auto, use all) |
| `spod_nfft` | 64 | SPOD FFT block size |
| `spod_noverlap` | 32 | SPOD block overlap (typically nfft/2) |

---

## Code Architecture

### Directory Structure

```
nekStab/
├── src/                          # Fortran source files (modules unless noted)
│   │
│   │── main.f90                  # Entry point, mode dispatcher (bare subroutines)
│   │── mode_config.f90           # Mode selection logic
│   │── usr_wrappers.f90          # Bare subroutine wrappers for .usr compatibility
│   │── nekstab_nek_bridge.f90    # Only module that includes Nek5000 headers
│   │── NEKSTAB                   # Shared common blocks / runtime parameters
│   │── makefile_nekStab          # Nek5000 build dependency rules
│   │
│   │── krylov_subspace.f90       # krylov_vector type, k_* operations
│   │── krylov_inner_products.f90 # Batched Gram/projection kernels (BLAS + single gop)
│   │── krylov_decomposition.f90  # Arnoldi factorization
│   │── eigensolvers.f90          # Krylov-Schur algorithm
│   │── lapack_wrapper.f90        # LAPACK wrappers (Schur, eig, etc.)
│   │── argsort.f90               # Sorting utilities
│   │── nek_vectors.f90           # nop* vector operations
│   │
│   │── matvec.f90                # Linearized NS operator
│   │── newton_krylov.f90         # Newton-GMRES solver
│   │── fixedp.f90                # SFD, BoostConv, TDF
│   │── sensitivity.f90           # Wavemaker, forcing response
│   │── otd.f90                   # Optimally Time-Dependent modes
│   │── fst.f90                   # Free-stream turbulence
│   │
│   │── energy_budget.f90         # PKE budget analysis
│   │── vortex.f90                # Vortex identification criteria
│   │── statistics.f90            # Time averaging
│   │── diagnostics.f90           # Energy, enstrophy, output
│   │
│   │── probes.f90                # Point monitoring
│   │── noise.f90                 # Initial condition perturbations
│   │── torque.f90                # Force/torque computation
│   │── forcing.f90               # External forcing, sponge
│   │
│   │── modal_analysis.f90        # POD/DMD/SPOD dispatcher
│   │── modal_pod.f90             # Proper Orthogonal Decomposition
│   │── modal_dmd.f90             # Dynamic Mode Decomposition
│   │── modal_spod.f90            # Spectral POD
│   │── modal_spod_streaming.f90  # Streaming SPOD
│   │── fourier.f90               # Fourier decomposition wrappers
│   │── fourier_fftw.f90          # FFT interface (FFTW3/MKL)
│   │
│   └── IO.f90                    # File I/O routines
│
├── bin/
│   │── mks                       # User-facing build wrapper
│   │── makeneks                  # Build backend invoked by mks
│   └── mks-completion.bash       # Shell completion for mks/makeneks
├── example/                      # Test cases
└── Nek5000/                      # Nek5000 solver
```

### Data Layout

nekStab uses Nek5000's spectral element data layout:

```fortran
! Velocity/temperature points per element
lx1, ly1, lz1 = polynomial order + 1 (e.g., 8 for N=7)

! Total points
lv = lx1 * ly1 * lz1 * lelv   ! Velocity grid
lt = lx1 * ly1 * lz1 * lelt   ! Temperature grid
lp = lx2 * ly2 * lz2 * lelv   ! Pressure grid (lx2 = lx1 or lx1-2)
```

### The `krylov_vector` Type

All stability computations use this derived type:

```fortran
type :: krylov_vector
   real, dimension(lv) :: vx, vy, vz    ! Velocity components
   real, dimension(lp) :: pr            ! Pressure
   real, dimension(lt, ldimt) :: t      ! Temperature + passive scalars
   real :: time                         ! For UPO period optimization
end type
```

---

## Module Reference

Most computational source files are wrapped in Fortran modules, providing explicit interfaces and compile-time type checking. The exceptions are `main.f90` and `usr_wrappers.f90`, which remain bare entry-point wrappers for Nek5000, plus `src/NEKSTAB`, which remains the shared include/common-block definition consumed by the build.

### Module Map

| Module | File | Key Public Interface |
|--------|------|---------------------|
| `nekstab_nek_bridge` | `nekstab_nek_bridge.f90` | Nek5000 `SIZE`/`TOTAL`/`ADJOINT` bridge, `nekStab_dp` |
| `krylov_subspace` | `krylov_subspace.f90` | `krylov_vector` type, `inner_product`, `norm`, `dnekclock`, `k_dot`, `k_norm`, `k_normalize`, `k_cmult`, `k_add2`, `k_add2s2`, `k_axpby`, `k_sub2`, `k_sub3`, `k_zero`, `k_copy`, `k_matmul`, `allocate_orbit`, `orbit_store`, `orbit_restore` |
| `krylov_inner_products` | `krylov_inner_products.f90` | `k_gram_matrix`, `k_project`, `k_gram_complex` |
| `nekstab_krylov_decomposition` | `krylov_decomposition.f90` | `arnoldi_factorization`, `update_hessenberg_matrix`, `arnoldi_checkpoint`, `log_transform` |
| `nekstab_eigensolvers` | `eigensolvers.f90` | `krylov_schur`, `inner_product`, `norm`, `outpost_ks`, `schur_condensation`, `select_eigenvalues`, `ensure_conjugate_pairs` |
| `nekstab_vectors` | `nek_vectors.f90` | `nopcopy`, `nopadd2`, `nopadd2s2`, `nopsub2`, `nopsub3`, `nopcmult`, `noprzero`, `nopaxpby`, `axpby`, `opadd3`, `opaddcol3` |
| `nekstab_lapack` | `lapack_wrapper.f90` | `schur`, `ordschur`, `eig`, `sort_eigendecomp`, `select_eigvals`, `lstsq`, `eig_symmetric`, `eig_hermitian` |
| `nekstab_argsort` | `argsort.f90` | `argsort` |
| `nekstab_matvec` | `matvec.f90` | `prepare_linearized_solver`, `matvec`, `forward_linearized_map`, `forward_finite_difference_map`, `adjoint_linearized_map` |
| `nekstab_io` | `IO.f90` | `whereyouwant`, `read_eigenvalue`, `load_mode_pair`, `load_files`, `k_load` |
| `nekstab_newton` | `newton_krylov.f90` | `newton_krylov`, `ts_gmres`, `initialize_gmres_vector`, `nonlinear_forward_map`, `set_nek5000_tolerances`, `spec_tole` |
| `nekstab_fixedpoint` | `fixedp.f90` | `sfd`, `BoostConv`, `tdf`, `boostconv_core`, `qr_dec`, `linear_system` |
| `nekstab_sensitivity` | `sensitivity.f90` | `wave_maker`, `bf_sensitivity`, `ts_steady_force_sensitivity`, `biorthogonalize`, `delta_forcing`, `animate_mode_only`, `animate_mode`, `animate_mode_Floquet` |
| `nekstab_energy_budget` | `energy_budget.f90` | `stability_energy_budget`, `stability_energy_budget_floquet`, `compute_velocity_gradient_tensor`, `compute_dissipation`, `compute_production` |
| `nekstab_diagnostics` | `diagnostics.f90` | `nekStab_energy`, `nekStab_enstrophy`, `outpost_vort`, `norm_grad`, `smooth_field`, `nekStab_outpost`, `nekStab_comment`, `nekStab_printNEKParams` |
| `nekstab_vortex` | `vortex.f90` | `vortex_core`, `compute_omega_jc`, `compute_omega`, `compute_q`, `compute_delta`, `compute_swirling` |
| `nekstab_statistics` | `statistics.f90` | `nekStab_avg` |
| `nekstab_noise` | `noise.f90` | `add_noise_scal`, `op_add_noise`, `add_symmetric_seed`, `mth_rand` |
| `nekstab_probes` | `probes.f90` | `pointcheck`, `zero_crossing` |
| `nekstab_torque_mod` | `torque.f90` | `nekStab_torque`, `nekStab_define_obj` |
| `nekstab_forcing_mod` | `forcing.f90` | `nekStab_forcing`, `nekStab_forcing_temp`, `activate_sponge`, `spng_init`, `spng_set`, `mth_stepf` |
| `nekstab_mode_config` | `mode_config.f90` | `nekStab_resolve_mode`, `nekStab_mode_from_string`, `nekStab_mode_from_flags`, `nekStab_mode_from_uparam`, `nekStab_validate_mode`, `nekStab_sync_uparam` |
| `nekstab_otd` | `otd.f90` | `otd`, `otd_construct_linear_operator`, `otd_compute_OTD_modes`, `otd_white_noise`, `otd_generate_forces`, `otd_orthonormalize_basis`, `otd_compute_FTLE` |
| `nekstab_fst` | `fst.f90` | `fst`, `fst_uin`, `fst_vin`, `fst_win` (public); internals: `readFSTinflow`, `defineBC`, `interpolateModes`, `computeBC`, `computeTurbu`, `spline`, `splint` |
| `nekstab_modal_analysis` | `modal_analysis.f90` | `modal_analysis`, `modal_compute_mean`, `modal_subtract_mean` |
| `modal_pod` | `modal_pod.f90` | `pod_compute`, `pod_fft_spectrum`, `hamming_window` |
| `modal_dmd` | `modal_dmd.f90` | `dmd_compute` |
| `modal_spod` | `modal_spod.f90` | `spod_compute`, `k_dot_complex`, `spod_save_modes` |
| `spod_streaming_state` | `modal_spod_streaming.f90` | `spod_s_init`, `spod_s_cleanup`, streaming SPOD state arrays and metadata |
| `modal_spod_streaming` | `modal_spod_streaming.f90` | `spod_streaming_batch` |
| `fourier` | `fourier.f90` | `nek_fourier_decomposition`, `fourier_decomposition`, `fourier_reconstruction`, `fft_init`, `fft_cleanup`, `fft_r2c`, `fft_c2r`, `fft_frequencies` |
| `fourier_fftw` | `fourier_fftw.f90` | `fft_init`, `fft_cleanup`, `fft_r2c`, `fft_c2r`, `fft_frequencies`, `fourier_decomposition`, `fourier_reconstruction` |

### Upgrading to v2.x

In nekStab v2.x, nearly all computational source files are wrapped in Fortran modules. This provides **explicit interfaces** and **compile-time type checking**, catching argument mismatches that previously caused silent runtime errors. The main exceptions are `main.f90` and `usr_wrappers.f90`, which remain bare entry points, and `src/NEKSTAB`, which remains the shared include/common-block definition.

The bridge module `nekstab_nek_bridge.f90` is now the only nekStab source file that directly includes Nek5000 headers (`SIZE`, `TOTAL`, `ADJOINT`). Other modules use the bridge instead of repeating `include` statements.

#### Changes for `.usr` files

If your `.usr` file calls nekStab routines directly, you need to add `use` statements at the top of the relevant subroutine. Without them, the compiler will report "undefined reference" errors.

```fortran
c     BEFORE (bare subroutine, implicit interface):
      call inner_product(alpha, px, py, pz, pp, pt,
     $                          qx, qy, qz, qp, qt)

c     AFTER (module, explicit interface -- compiler checks types):
      use krylov_subspace, only: inner_product
      call inner_product(alpha, px, py, pz, pp, pt,
     $                          qx, qy, qz, qp, qt)
```

#### Commonly needed modules in `.usr` files

| Module | Provides |
|--------|----------|
| `use krylov_subspace` | `inner_product`, `norm`, `krylov_vector` |
| `use nekstab_vectors` | `nopcopy`, `noprzero`, `nopadd2`, `nopcmult`, `nopaxpby` |
| `use nekstab_io` | `whereyouwant`, `load_files`, `k_load` |
| `use nekstab_diagnostics` | `nekStab_energy`, `nekStab_enstrophy`, `outpost_vort` |

#### Auto-wrapped routines (no `use` needed)

The following routines are called from `.usr` callback hooks and are automatically wrapped in `usr_wrappers.f90` with bare subroutine interfaces. You do **not** need `use` statements for these:

- `nekStab_torque` -- called from `userbc` or `userchk`
- `nekStab_forcing` -- called from `userf`
- `nekStab_define_obj` -- called from `usrdat2`

---

## Vector Operations

nekStab provides two levels of vector operations following BLAS conventions.

### Level 1: Nek5000's `op*` Routines

These operate on **velocity only** (vx, vy, vz):

| Routine | Operation | Notes |
|---------|-----------|-------|
| `opcopy(a1,a2,a3, b1,b2,b3)` | a = b | Copy vectors |
| `opadd2(a1,a2,a3, b1,b2,b3)` | a = a + b | In-place add |
| `opsub2(a1,a2,a3, b1,b2,b3)` | a = a - b | In-place subtract |
| `opsub3(a,b,c, d,e,f, g,h,i)` | a = d - g | Three-operand subtract |
| `opcmult(a1,a2,a3, c)` | a = c*a | Scalar multiply |
| `oprzero(a1,a2,a3)` | a = 0 | Zero vector |
| `opcolv(a1,a2,a3, c)` | a = a.*c | Element-wise multiply |

### Level 2: nekStab's `nop*` Routines

These operate on **all fields** (velocity + pressure + scalars):

```fortran
! Signature: (vx, vy, vz, pr, t)
call nopcopy(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)   ! a = b
call nopadd2(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)   ! a = a + b
call nopsub2(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)   ! a = a - b
call nopcmult(a1,a2,a3,a4,a5, c)               ! a = c*a
call noprzero(a1,a2,a3,a4,a5)                  ! a = 0
call nopaxpby(a,alpha, b,beta)                 ! a = alpha*a + beta*b
```

**Scalar Loop Pattern**: All `nop*` routines handle multiple scalars:

```fortran
! Temperature (if ifto = .true.)
if (ifto) call copy(a5(1,1), b5(1,1), n)

! Passive scalars (if ifpsco(k) = .true.)
if (ldimt > 1) then
   do k = 1, npscal
      if (ifpsco(k)) call copy(a5(1,k+1), b5(1,k+1), n)
   end do
end if
```

### Level 3: nekStab's `k_*` Routines

High-level operations on `krylov_vector` type:

| Routine | Operation | Description |
|---------|-----------|-------------|
| `k_dot(alpha, p, q)` | α = ⟨p,q⟩ | Mass-weighted inner product |
| `k_norm(alpha, p)` | α = ‖p‖ | L2 norm |
| `k_normalize(p, alpha)` | p = p/‖p‖ | Normalize, return norm |
| `k_cmult(p, c)` | p = c*p | Scalar multiply |
| `k_add2(p, q)` | p = p + q | Add vectors |
| `k_sub2(p, q)` | p = p - q | Subtract vectors |
| `k_copy(p, q)` | p = q | Copy vector |
| `k_zero(p)` | p = 0 | Zero vector |

**Inner Product**: Uses mass matrix weighting for spectral elements:

```fortran
alpha = glsc3(p%vx, q%vx, bm1s, nv) + glsc3(p%vy, q%vy, bm1s, nv)
if (if3d) alpha = alpha + glsc3(p%vz, q%vz, bm1s, nv)
if (ifto) alpha = alpha + glsc3(p%t(:,1), q%t(:,1), bm1s, nt)
! ... loop over passive scalars
```

For dense Krylov algebra, `krylov_inner_products.f90` provides batched kernels:
- `k_gram_matrix` builds Gram matrices with BLAS `dgemm` plus one global reduction
- `k_project` computes Arnoldi/GMRES projections with BLAS `dgemv`
- `k_gram_complex` assembles SPOD cross-spectral density blocks with a single packed reduction

---

## Mesh Generation

### PointWise to Nek5000 Workflow

This workflow converts a PointWise mesh to Nek5000 format via Gmsh.

#### Step 1: Create Mesh in PointWise

Create your mesh in PointWise (v18+) and save the project as `.pw` file.

#### Step 2: Configure Solver Settings

1. **CAE → Select Solver → Gmsh**
2. **Set Dimension → 3D**
3. **Set Boundary Conditions** - Create names with arbitrary boundary IDs:
   ```
   inf  1
   out  2
   top  3
   bot  4
   fsp  5
   ```
   > Note: These names don't need to follow Nek's 3-character style yet.

   For **periodic** spanwise (instead of free-slip), tag each side separately:
   ```
   P1  5
   P2  6
   ```

4. **Solver Attribute → Q1** (second-order elements will be set in Gmsh)

#### Step 3: Export from PointWise

1. **File → Export → CAE** → `mesh_file1.msh`
2. **Data Precision → Double** → OK

#### Step 4: Optimize in Gmsh

Open Gmsh and load `mesh_file1.msh`:

1. Expand the **Mesh** menu on the left panel
2. Click **Optimize 3D**
3. Click **Set order 2**
4. **File → Export** → `mesh_file2.msh`
5. **MSH Options → Version 2 ASCII** (ensure other options are not selected)

#### Step 5: Convert with gmsh2nek

Run `gmsh2nek` and follow the prompts:

```
$ gmsh2nek

3                          # 3D mesh

mesh_file2                 # mesh filename (without .msh)

# Output shows BID mapping - SAVE THESE VALUES!

1                          # 1 if periodic, 0 if SYM

4 5                        # periodic pair IDs (P1, P2 - may have changed!)

0 0 12                     # spanwise: x, y, z extent (e.g., -6 to 6 → 12)
```

#### Step 6: Generate Connectivity Map

Run `genmap`:

```
$ genmap

mesh_file2                 # mesh filename

                           # press Enter for default tolerance, or specify custom
```

You should now have:
- `mesh_file2.re2` - mesh file
- `mesh_file2.ma2` - connectivity map

#### Step 7: Configure Boundary Conditions in .usr

In `usrdat2` in your `.usr` file, set boundary conditions.

> **Important:** Use the boundary IDs from `gmsh2nek` output, NOT the original PointWise IDs!

> **Important:** Nek5000 requires 3-character BC codes. See [Nek5000 documentation](https://nek5000.github.io/NekDoc/problem_setup/boundary_conditions.html) for options.

```fortran
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

      integer iel, ifc, id_face

      do iel = 1, nelv
         do ifc = 1, 2*ndim
            id_face = bc(5, ifc, iel, 1)
            if (id_face .eq. 2) then       ! outflow
               cbc(ifc, iel, 1) = 'O  '
            elseif (id_face .eq. 3) then   ! normal outflow
               cbc(ifc, iel, 1) = 'ON '
            elseif (id_face .eq. 6) then   ! wall
               cbc(ifc, iel, 1) = 'W  '
            elseif (id_face .eq. 7) then   ! prescribed inflow
               cbc(ifc, iel, 1) = 'v  '
            endif
         enddo
      enddo

      return
      end
```

### Alternative: Gmsh Native

For simpler geometries, create meshes directly in Gmsh:

```bash
gmsh -3 geometry.geo -order 2 -o mesh.msh
gmsh2nek
genmap
```

---

## Examples

### Available Test Cases

| Directory | Flow | Features Demonstrated |
|-----------|------|----------------------|
| `cylinder/` | Circular cylinder wake | DNS, Newton, stability, Floquet, wavemaker, OTD, POD/DMD/SPOD, RANS |
| `back_fstep/` | Backward-facing step | Newton base flow, transient growth |
| `blasius/` | Flat plate boundary layer | TS waves |
| `cubic_cavity/` | 3D cubic cavity | 3D stability |
| `cubic_cavity_upo/` | 3D cubic cavity | Newton-GMRES periodic orbit |
| `flip_flop/` | Side-by-side cylinders | Neimark-Sacker, Floquet |
| `lid_driven/` | Lid-driven cavity | Confined flow bifurcations |
| `naca0012/` | Airfoil at incidence | Newton base flow |
| `naca0012_Re2500/` | Airfoil at incidence | Direct stability |
| `parque/` | Wind farm / actuator disks | RANS wake simulation |
| `poiseuille/` | Channel flow | Canonical stability / parity reference |
| `poiseuille_OTD/` | Plane Poiseuille flow | OTD and FTLE tracking |
| `poiseuille_RANS/` | Turbulent channel flow | RANS base flow and finite-difference stability |
| `slot_FST/` | Flat plate + FST | Free-stream turbulence inflow |
| `thermosyphon/` | Buoyancy-driven convection | Pitchfork + Hopf bifurcations |
| `torus/` | Toroidal pipe | Direct stability on baseline torus mesh |
| `torus2/` | Toroidal pipe | Direct stability on variant mesh |
| `torus_full/` | Toroidal pipe | Direct stability on full torus geometry |
| `torus_pulsed/` | Toroidal pipe | Newton-GMRES forced periodic orbit |
| `tpjet/` | Forced jet | Floquet period-doubling |

### Cylinder Workflow

```
example/cylinder/
├── dns/           # 1. Verify mesh with DNS
├── baseflow/
│   ├── newton/    # 2. Compute steady base flow
│   └── sfd/       # (alternative method)
├── stability/
│   ├── direct/    # 3. Direct eigenmodes
│   └── adjoint/   # 4. Adjoint eigenmodes
└── postproc/
    ├── sensitivity_budget_wavemaker/
    └── steady_force_sensitivity/
```

### Running a Complete Analysis

```bash
# 1. DNS verification
cd example/cylinder_re100/000_dns
mks 1cyl && nekbmpi 1cyl 4
# Check: flow develops vortex shedding

# 2. Base flow
cd ../baseflow/newton
cp ../dns/rst_1cyl0.f00001 .  # Use DNS snapshot as initial guess
# Edit 1cyl.par: userParam01 = 2.0
mks 1cyl && nekbmpi 1cyl 4
# Check: residual → 0, output BF_1cyl0.f00001

# 3. Direct stability
cd ../stability/direct
ln -s ../../baseflow/newton/BF_1cyl0.f00001 .
# Edit 1cyl.par: userParam01 = 3.1
mks 1cyl && nekbmpi 1cyl 4
# Check: eigenvalues in logfile, modes in dRe*, dIm*

# 4. Adjoint stability
cd ../adjoint
ln -s ../../baseflow/newton/BF_1cyl0.f00001 .
# Edit 1cyl.par: userParam01 = 3.2
mks 1cyl && nekbmpi 1cyl 4

# 5. Wavemaker
cd ../../postproc/sensitivity_budget_wavemaker
ln -s ../../stability/direct/dRe* ../../stability/direct/dIm* .
ln -s ../../stability/adjoint/aRe* ../../stability/adjoint/aIm* .
# Edit 1cyl.par: userParam01 = 4.2
mks 1cyl && nekbmpi 1cyl 4
```

### Plot Utilities

Shared plotting helpers now live in `example/nekplot.py`. It provides reusable field I/O, interpolation, residual-history plots, spectra plots, and OTD diagnostics for the example scripts.

To gather per-case `plot.png` outputs into one validation folder:

```bash
python example/collect_plots.py
python example/collect_plots.py --generate
python example/collect_plots.py --list
```

The collector writes consolidated figures to `validation/figures/`.

---

## Validation

### Automated Validation Suite

The repository now ships with `validation/validate.py`, which combines literature-facing validation checks with broad example coverage.

```bash
./validation/validate.py
./validation/validate.py --short
./validation/validate.py --check-only
./validation/validate.py --dry-run
./validation/validate.py --list
./validation/validate.py --nprocs 8 cylinder
./validation/validate.py --compile-all
```

Validation is split into two tiers:
- **Short tier**: AMR-oriented literature checks for `cylinder`, `thermosyphon`, `flipflop_bf`, `flipflop_floquet`, `backstep`, `tpjet_bf`, and `tpjet_floquet`
- **Full tier**: compile/run coverage across the cases currently listed in `validation/validate.py`, spanning cylinder variants plus `back_fstep`, `thermosyphon`, `flip_flop`, `tpjet`, `lid_driven`, `cubic_cavity`, `naca0012`, `poiseuille_OTD`, and `slot_FST`

Operational behavior of `validation/validate.py`:
- Auto-detects CPU topology and prefers physical cores
- Reads `SIZE` (`lelg`, `lpmin`) to choose a mesh-aware MPI rank count
- Materializes missing prerequisites from `NEKSTAB_DATA_ROOT` when set, otherwise from the default `$HOME/.data_baptiste.nosync`, before deciding whether to skip a case
- Skips a case only when required restart/base-flow files are missing both locally and in the configured external data roots
- Treats missing post-run `copies_to` artifacts as a validation failure so downstream dependent cases are not left with silent gaps
- `--compile-all` scans `example/`, ignores hidden placeholder `.usr` files and non-case directories, runs clean builds, writes `build.log`, and surfaces the first compiler error when a build fails

The following validation cases are documented in [Frantz et al. (2023)](https://doi.org/10.1115/1.4056808). Each demonstrates a specific bifurcation type and has been verified against published literature.

### Annular Thermosyphon (Pitchfork + Hopf)

A two-dimensional flow in a concentric annular enclosure heated from below. The inner-to-outer radius ratio is R₂/R₁ = 2. The lower wall is heated (θ = 1), upper wall cooled (θ = 0).

**Governing equations**: Incompressible Navier-Stokes with Boussinesq approximation.

**Control parameters**:
- Rayleigh number: Ra = ρgβΔT(R₂-R₁)³/(μα)
- Prandtl number: Pr = 5 (fixed)

**Bifurcations**:

| Bifurcation | Type | Critical Ra | Strouhal | Mode Character |
|-------------|------|-------------|----------|----------------|
| Primary | Pitchfork | Ra_c1 ≈ 494 | 0 (steady) | Symmetry-breaking convection cell |
| Secondary | Hopf | Ra_c2 ≈ 16,081 | St ≈ 7 | Oscillating convection |

**Reference**: [Loiseau et al., TCFD (2020)](https://doi.org/10.1007/s00162-019-00518-1)

**Validation details**:
- Mesh: 32×8 spectral elements (256 total), N=7 polynomial order
- Pitchfork: Single real eigenvalue crosses into unstable half-plane at Ra ≈ 494
- Hopf: Complex-conjugate pair crosses at Ra ≈ 16,081, frequency matches DNS

**Computational performance** (Ra = 16,100 base flow):
- Newton-Krylov: 147 seconds
- SFD: 1071 seconds (7× slower)

---

### Harmonically Forced Jet (Period-Doubling)

An axisymmetric jet forced at the inflow with 5% amplitude oscillation at Strouhal number St = 0.6. The configuration promotes vortex pairing via subharmonic instability.

**Geometry**: Domain 40D × 5D (streamwise × radial), 160×30 spectral elements, N=5.

**Inflow condition**:
```
u(r,t) = ½[1 - tanh(1/(4θ₀(r - 4r⁻¹)))] × (1 + A cos(ωt))
```
where A = 0.05, θ₀ = 0.025, ω = 2π St.

**Bifurcation**:

| Bifurcation | Type | Critical Re | Floquet μ | Physical Mechanism |
|-------------|------|-------------|-----------|-------------------|
| Secondary | Period-doubling | Re_c ≈ 1371 | μ = -1 | Vortex pairing |

**Reference**: [Leopold et al., JFM (2019)](https://doi.org/10.1017/jfm.2019.607)

**Validation**:
- nekStab: Re_c = 1371.18
- Reference: Re_c ≈ 1371
- Agreement: < 0.02%

The Floquet multiplier exits the unit circle at μ = -1 (characteristic of period-doubling). Above Re_c, vortices spontaneously pair, and the subharmonic St = 0.3 appears in spectra.

---

### Flow Past a Circular Cylinder (Hopf + Floquet Pitchfork)

The canonical bluff body wake exhibiting both primary (2D) and secondary (3D) instabilities.

**Mesh**: 1464 elements (2D), N=5. For 3D Floquet: extruded to 10 spanwise elements.

#### Primary Instability (Hopf)

| Parameter | nekStab | Reference | Source |
|-----------|---------|-----------|--------|
| Re_c1 | ≈ 46.6 | 46.6 | Jackson (1987), Kumar & Mittal (2006) |
| St_c | 0.125 | 0.125 | — |

The leading eigenvalue is a complex-conjugate pair crossing into the unstable half-plane. The corresponding eigenvector is the von Kármán mode.

#### Secondary Instability (Floquet/Mode A)

Three-dimensional perturbations on the 2D periodic wake. The critical spanwise wavelength is λ_z ≈ 4D (β_c = 1.585).

| Parameter | nekStab | Reference | Source |
|-----------|---------|-----------|--------|
| Re_c2 | ≈ 189 | 188.5 ± 1 | Barkley & Henderson (1996) |
| μ at Re=190 | 1.012 | 1.034 | Barkley (2005) |
| μ at Re=190 | 1.012 | 1.002 | Giannetti et al. (2010) |
| St at Re=190 | 0.196 | 0.195 | Barkley & Henderson (1996) |

The Floquet multiplier exits the unit circle at μ = +1 (synchronous mode), characteristic of a pitchfork bifurcation. The flow three-dimensionalizes while retaining the shedding frequency.

**References**:
- Primary: [Zebib, J. Eng. Math. (1987)](https://doi.org/10.1007/BF00127673)
- Secondary: [Barkley & Henderson, JFM (1996)](https://doi.org/10.1017/S0022112096008750)
- Floquet comparison: [Barkley, PoF (2005)](https://doi.org/10.1063/1.1868171), [Giannetti et al., JFM (2010)](https://doi.org/10.1017/S0022112010003083)

---

### Side-by-Side Cylinders (Neimark-Sacker)

Two circular cylinders placed side-by-side with gap g = 0.7D. This configuration exhibits the "flip-flop" instability—a quasi-periodic transition.

**Mesh**: 5092 elements, N=7. Domain: -50D to 75D (streamwise), ±50D (cross-stream).

**Bifurcations**:

| Bifurcation | Type | Critical Re | Frequency | Mode |
|-------------|------|-------------|-----------|------|
| Primary | Hopf | Re_c1 ≈ 55 | St = 0.11 | Synchronized vortex shedding |
| Secondary | Neimark-Sacker | Re_c2 ≈ 61.17 | St = 0.02 | Flip-flop |

**Reference**: [Carini et al., JFM (2014)](https://doi.org/10.1017/jfm.2014.362)

**Validation**:
- nekStab: Re_c2 = 61.17
- Reference: Re_c2 ≈ 61.6
- Agreement: < 0.7%

The Neimark-Sacker bifurcation is characterized by a complex-conjugate Floquet pair exiting the unit circle at angle φ ≈ 71°. Above Re_c2, the system exhibits quasi-periodic dynamics (torus in phase space) with two incommensurate frequencies.

---

### Backward-Facing Step (Transient Growth)

A linearly stable flow that exhibits strong transient amplification due to non-normality of the linearized operator.

**Geometry**: Step height h = 1, expansion ratio 1:2. Domain: -10h to 50h (streamwise), -1h to 1h (wall-normal). Mesh: 1670 elements, N=5.

**Reference**: [Blackburn et al., JFM (2008)](https://doi.org/10.1017/S002211200800267X)

**Results at Re = 500**:

| Quantity | nekStab | Reference |
|----------|---------|-----------|
| τ_opt | 58 | 58 |
| G(τ_opt) | Matches | See Fig. 22 in Frantz et al. (2023) |

The optimal perturbation consists of streamwise vortices in the shear layer. The optimal response (at τ = 58) shows amplified streaks downstream of the step. The gain envelope shows excellent quantitative agreement with the reference.

**Method**: Mode 3.3 (transient growth). The direct-adjoint iteration converges to the leading singular triplet of exp(τL).

---

### Summary of Validation Cases

| Case | Bifurcation | Critical Parameter | Agreement |
|------|-------------|-------------------|-----------|
| Thermosyphon | Pitchfork | Ra_c ≈ 494 | Literature |
| Thermosyphon | Hopf | Ra_c ≈ 16,081 | Literature |
| Forced jet | Period-doubling | Re_c ≈ 1371 | < 0.02% |
| Cylinder (2D) | Hopf | Re_c ≈ 46.6 | Literature |
| Cylinder (3D) | Floquet pitchfork | Re_c ≈ 189 | < 0.5% |
| Side-by-side | Neimark-Sacker | Re_c ≈ 61.17 | < 0.7% |
| Back step | Transient growth | τ_opt = 58 | Literature |

All validation cases are available in the `example/` directory with ready-to-run configurations.

---

## Theoretical Background

### Linearized Navier-Stokes

For base flow **U**(x), perturbation **u'**(x,t) satisfies:

```
∂u'/∂t + (U·∇)u' + (u'·∇)U = -∇p' + (1/Re)∇²u'
∇·u' = 0
```

### Matrix-Free Eigenvalue Problem

Seeking **u'** = **q̂** exp(σt), the eigenvalue problem is:

```
σq̂ = Aq̂
```

where **A** is the linearized operator. nekStab never forms **A** explicitly. Instead:

1. Time-step the linearized equations: q(T) = exp(AT)q(0)
2. The propagator exp(AT) shares eigenvectors with **A**
3. Eigenvalues: σ = log(μ)/T where μ is propagator eigenvalue

### Krylov-Schur Algorithm

1. **Arnoldi iteration**: Build orthonormal basis V and Hessenberg matrix H
   ```
   AV_m = V_m H_m + h_{m+1,m} v_{m+1} e_m^T
   ```
   Current default: CGS2 with reorthogonalization (`use_cgs = .true.`), which reduces MPI reductions relative to the older MGS path while keeping the fallback available.

2. **Ritz extraction**: Eigenvalues of H_m approximate eigenvalues of A

3. **Implicit restart**: Schur decomposition filters unwanted eigenvalues

4. **Deflation**: Lock converged eigenvalues, continue with reduced subspace

### Adjoint Operator

The continuous adjoint satisfies:

```
-∂u†/∂t - (U·∇)u† + (∇U)^T·u† = -∇p† + (1/Re)∇²u†
```

Discrete adjoint uses the **transpose** of the linearized time-stepper.

### Structural Sensitivity

The wavemaker S(x) identifies where feedback most affects eigenvalues:

```
S(x) = |q̂(x)| · |q̂†(x)| / ∫ q̂† · q̂ dV
```

High S(x) regions are sensitive to local modifications (e.g., control devices).

---

## Continuous Integration

nekStab uses GitHub Actions to automatically test compilation and execution across multiple platforms, compilers, and MPI implementations.

### Build Matrix

| OS | Arch | Compiler | MPI | Status |
|:---|:----:|:---------|:----|:------:|
| Ubuntu | x86_64 | gfortran 14 | OpenMPI | [![CI](https://github.com/nekStab/nekStab/actions/workflows/ci.yml/badge.svg)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| Ubuntu | x86_64 | gfortran 14 | MPICH | — |
| macOS 26 | ARM64 | gfortran 14 | OpenMPI | — |
| Ubuntu 24.04 | x86_64 | ifx 2025.2 | Intel MPI | — |

> **Note:** All configurations share the same CI badge. Individual job status can be viewed on the [Actions page](https://github.com/nekStab/nekStab/actions/workflows/ci.yml).

### What the CI Tests

1. **Compilation** — Builds the `example/cylinder_re100/000_dns/ci_test` case with each compiler/MPI combination
2. **Smoke Test** — Runs 10 DNS timesteps with 2 MPI ranks to verify basic functionality
3. **Output Verification** — Confirms the logfile contains successful completion markers

### Running Tests Locally

To replicate CI tests locally:

```bash
cd example/cylinder_re100/000_dns/ci_test
mks 1cyl                           # Compile
echo "1cyl" > SESSION.NAME
pwd >> SESSION.NAME                # Create session file
mpirun -np 2 ./nek5000 > logfile 2>&1
grep "run successful" logfile      # Verify completion
```

### Triggering CI

CI runs automatically on:
- Push to `dev` branch
- Pull requests to `dev` branch
- Manual dispatch via Actions page

### Adding New Test Cases

To add a new CI test case:

1. Create a minimal configuration in `example/<case>/ci_test/`
2. Use `stopAt = numSteps` with `numSteps = 10` for fast execution
3. Use cold start (no restart file) for simplicity
4. Add the build/run steps to `.github/workflows/ci.yml`

---

## Troubleshooting

### Compilation Issues

**Module file not found**
```
Cannot open module file 'krylov_subspace.mod'
```
Solution: `mks case --fresh` (parallel make race condition)

**MKL libraries not found**
```
ld: cannot find -lmkl_intel_lp64
```
Solution: `source /opt/intel/oneapi/setvars.sh` before compilation

**ifx segmentation fault**
Ensure you have the latest nekStab with ifx compatibility fixes (trim() on assumed-length strings).

### Runtime Issues

**Eigenvalues not converging**
- Increase `k_dim` (more Krylov vectors)
- Decrease `eigen_tol` (looser tolerance initially)
- Check time step: too large can miss fast modes

**Newton not converging**
- Verify initial guess is reasonable (close to solution)
- Try SFD first to get closer to steady state
- Reduce `epsilon_base` for better Jacobian approximation
- If `ifdyntol = .true.`, add `ew_tol_cap` for stiff cases that need a tighter inner-solve cap
- Check whether the new stagnation warning or divergence guard is firing in the Newton log

**Floquet modes incorrect**
- Ensure `ifstorebase = .true.`
- Verify base flow period matches actual period
- Check sufficient time resolution over one period

### Performance

- **Optimal element count**: 1000-5000 elements per MPI rank
- **Memory**: ~8 bytes × DOFs × k_dim × 2 (for V and H)
- **Intel compilers**: 20-40% faster than GCC on x86

---

## Citation

```bibtex
@article{frantz2023krylov,
    author = {Frantz, R. A. S. and Loiseau, J.-Ch. and Robinet, J.-Ch.},
    title = "{Krylov Methods for Large-Scale Dynamical Systems:
              Application in Fluid Dynamics}",
    journal = {Applied Mechanics Reviews},
    volume = {75},
    number = {3},
    pages = {030802},
    year = {2023},
    doi = {10.1115/1.4056808},
}
```

---

## License

BSD-3-Clause. See [LICENSE](LICENSE).

## Contact

- Ricardo Frantz: rasfrantz@gmail.com
- GitHub Issues: https://github.com/nekStab/nekStab/issues
