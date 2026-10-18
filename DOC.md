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
| `'dmt'` | Dynamic Mode Tracking base flow | 1.3 |
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
| `'otd'` | Optimally Time-Dependent modes, frozen base flow | 5 |
| `'otd_dns'` | Optimally Time-Dependent modes, base flow advanced by the DNS | 5.1 |
| `'pod'` | Proper Orthogonal Decomposition | 6.1 |
| `'dmd'` | Dynamic Mode Decomposition | 6.2 |
| `'spod'` | Spectral POD | 6.3 |
| `'dmt_bank'` | DMT filter bank on the running DNS | 6.4 |

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
| OTD | `ifotd`; `ifotdDNS` (also sets `ifotd`) for the evolving base flow |
| Modal Analysis | `ifpod`, `ifdmd`, `ifspod`; `ifdmtbank` (filters on the running DNS) |

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

Only the codes in the table above are accepted (the value times 100 is an integer code, so `3.11` and `3.110` are the same mode). Any other value, for example `3.15` or `7.3`, stops the run with `Unknown mode, userParam01 = ...`. Earlier versions ran a DNS with a warning.

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
| `1.3` | DMT | Dynamic Mode Tracking (Queguineur et al. 2019): band-pass filter around a target frequency |
| `1.4` | TDF | Time-Delayed Feedback |

**SFD Parameters** (`userParam01 = 1.1`):
- `userParam04`: Frequency `f` of the oscillation to damp, in cycles per time unit (`ω = 2πf`). The sign selects the formulation: positive for Åkervik, negative for Casacuberta.
- `userParam05`: Growth rate `σ` of that oscillation (use `0` for continuation from previous run). Åkervik gives the filter the cutoff `ω/2` and the gain `−2σ`; Casacuberta gives it the cutoff `(√(ω² + σ²) − σ)/2` and the gain `−(√(ω² + σ²) + σ)/2`. The run prints both values (`Akervik  cutoff, gain:` or `Casacub. cutoff, gain:`).

**DMT Parameters** (`userParam01 = 1.3`):
- `userParam11`: target angular frequency `ω_c` (`2π St`)
- `userParam12`: filter damping `β`
- `userParam13`: feedback gain `χ`
- `userParam14`: time before the forcing starts
- `userParam15`: convergence tolerance of the residual

DMT adds the force `-χ (u - ū)`, where `ū` follows `u` through `ū'' + β ū' + ω_c² ū = u'' + ω_c² u`. The difference `u - ū` is the band-pass of `u` around `ω_c`: its gain is 1 at `ω_c` and 0 at `ω = 0`. A steady flow is therefore the only state without force, and DMT finds steady base flows; a periodic orbit has content in the band at its own frequencies, so the residual does not reach zero on it. The run ends, and writes `BF_*` in double precision, when the residual (`residu_dmt.dat`, column 2: the root mean square over the domain of `u - ū`, divided by the time step) falls below `userParam15`. It needs a constant time step (`variableDt = no`). Take `ω_c` from a stability computation of a nearby flow. On the Re = 100 cylinder, `ω_c = 0.7337` (the leading mode of `cylinder_re100/310_stability_direct`, growth rate `σ = 0.1248`), `β = 1` and `χ = 0.5` converge at t = 766 (`cylinder_re100/130_baseflow_dmt`); `χ = 1` and `β = 0.35` did not converge by t = 800.

SFD ends, and writes the converged base flow `BF_*`, when its residual falls below the larger of the two `residualTol` of the `.par`. The residual (`residu.dat`, column 2) is the root mean square over the domain of the change of the state during one time step, so it does not depend on the size of the domain. With a thermal or RANS state the temperature and every scalar enter it, with the weight of the Krylov inner product. TDF uses the same measure for the difference between the state and the state one period earlier, and BoostConv for the change of the velocity. Once the state has converged, the residual stops falling at the noise of the solves, so `residualTol` has to lie above that level; for `cylinder_re100/110_baseflow_sfd` the level is about 8e-12, against `residualTol = 1e-10`.

**TDF Parameters** (`userParam01 = 1.4`):
- `userParam05`: frequency `St` of the forcing that sets the orbit. The delay is its period `T = 1/St`.

TDF holds an unstable periodic orbit of a forced flow. It adds the force `g (q(t) - q(t - T))` to every field of the state, with `g = -0.04432 · 2π / T` fixed in the code; `g` is negative, so the force opposes the distance from the orbit, and it is zero on the orbit. The run keeps the last period of the state in memory, `norbit = T/dt` steps of velocity (and of every stored scalar), and adjusts `dt` so that one period holds a whole number of steps. The delay is then the period to machine accuracy, which the method needs: a delay that differs from the period by `δT` leaves a force of the order `|g| δT |dq/dt|` on the orbit, and the residual cannot fall below that level. TDF therefore suits orbits whose period is set by a forcing (an inflow or a mesh motion), where `T` is known exactly. The memory is `2 lv norbit` double-precision values in 2D and `3 lv norbit` in 3D (`lv` of `SIZE`), and the run stops with a message when it cannot allocate them. The residual (`residu.dat`: time, the root mean square over the domain of `q(t) - q(t - T)`, its change per step) is the measure of SFD, and the run writes `BF_*` in double precision when it falls below the larger `residualTol`. The gain is not a parameter, so TDF does not hold every orbit: the forced jet of `tpjet_re2000` at polynomial order 5 converges (residual 6.5e-8 at t = 600, falling as t^-2.5), while an earlier test of the same jet at polynomial order 4 diverged (residual 6.3e-3 at t = 15 and a CFL of 6).

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
| `4.41` | Real part forcing sensitivity (not validated in 2.0.0) |
| `4.42` | Imaginary part forcing sensitivity (not validated in 2.0.0) |
| `4.43` | Delta forcing response |
| `4.50` | Animate direct mode |
| `4.51` | Animate mode with base flow deformation |
| `4.52` | Animate Floquet mode |

#### Steady force sensitivity (4.41, 4.42)

Not validated in 2.0.0, and no example ships for it. The modes solve `(I - exp(T L†)) x = c` for the
sensitivity to a steady force, by restarted GMRES on the time-stepper form of the adjoint problem. The right-hand side
`c` comes from the files `sr_<session>0.f00001` or `si_<session>0.f00001`, which the base flow sensitivity (mode 4.3, also
run by 4.0) writes. On the cylinder at Re = 100 (the stage `412_postproc_steady_force_sensitivity` of
2.0.0-rc3 and rc4, now removed) the solve did not reach its target:

| T (`endTime`) | `k_dim` | GMRES restarts | squared relative residual at the end (target 1e-6) |
|---|---|---|---|
| 1 | 50 | 10 | 1.35e-2 |
| 2 | 50 | 10 | 1.46e-2 |
| 3 | 50 | 6 (run stopped) | 1.64e-2 |
| 1 | 400 | 2 (run stopped) | 7.4e-3 |

The adjoint map also grows quickly with T. For an input of unit norm, the squared norm of the output was 0.10 at T = 1,
0.20 at T = 2, 19 at T = 4 and 3.9e25 at T = 10. At T = 10 two applications to the same input gave the same output, and the
map was linear to 6e-20 in the squared relative error, so the growth belongs to the operator and is not noise. The cause of
the residual plateau and of this growth is not established. The mode prints a warning when GMRES ends above its target, and
it still writes the field. Treat that field as unverified.

### Mode 5: OTD

Optimally Time-Dependent modes (Babaee and Sapsis, 2016). A basis of `lpert` orthonormal perturbations (`lpert` in `SIZE`) evolves
with the linearized Navier-Stokes equations and follows the subspace of fastest growth. The reduced operator and the finite-time
exponents are written to `otd_eigenvalues.dat`, `otd_growth_rates.dat` and `otd_ftle.dat`.

| Value | String | Base flow |
|-------|--------|-----------|
| `5` | `'otd'` | Frozen. The field loaded from `BF_<session>0.f00001` is not advanced; the `useric` base-flow branch (`JP == 0`) must return `ubic`, `vbic`, `wbic`, which the driver filled from that file. |
| `5.1` | `'otd_dns'` | Evolving. The DNS advances the base flow next to the perturbations, so the basis follows the trajectory. With `otd_FTLEPeriod` set to the period of an orbit, `otd_ftle_period.dat` holds the exponents of each period. |

The OTD perturbations start from the files `r01<session>0.f<n>`, `r02<session>0.f<n>`, ... (the highest `n` found) and from white
noise when there are none. A run writes files with these names, so remove those of an earlier run before seeding. About a steady base
flow, the real and imaginary parts of the leading eigenmode of `310_stability_direct` are an exact start: in
`cylinder_re100/500_otd` the eigenvalues of the reduced operator at the start equal the growth rate 0.12479 of the 310 stage to 4e-5
(the difference is the time discretisation of the 310 map, which `Lr` does not have: it falls by a factor 4 when the time step of the
310 stage is halved). The two diagonal entries of a complex pair need not be equal: in `cubic_cavity_re1914/500_otd` they swing between +0.16 and
-0.16 as the basis turns in the eigenplane, and their mean is the growth rate (-1.8320e-4 against -1.8296e-4 of the 310 stage).
From white noise the basis needs 100 to 200 convective time units to settle: 100 for the steady wake of
`cylinder_re100/500_otd` and 190 for the orbit of `cylinder_re180/500_otd`. A `userParam` value between 5 and 5.1 is not a mode
and stops the run.

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
| `pod_spectrum.dat` | POD eigenvalue spectrum | `mode, eigenvalue, energy_%, cumulative_%`; the first `modal_nsave` rows add `‖Φ‖` |
| `pod_coefficients.dat` | POD temporal coefficients | `snapshot, a_1, a_2, ...` |
| `pod_fft_spectrum.dat` | Welch power spectrum of the POD coefficients (written by `ifpod` and by `ifspod`) | `St, P_mode 1, P_mode 2, ...` |
| `dmd_spectrum.dat` | DMD eigenvalues, ordered by \|μ\| | `mode, \|μ\|, σ, ω, St, Re(μ), Im(μ)`; the first `modal_nsave` rows add `‖Φ_re‖` |
| `dmd_svd.dat` | Singular values of the snapshot matrix and the rank used | `i, σ_i, energy_%, cumulative_%` |
| `spod_stream_spectrum.dat` | SPOD eigenvalues of every frequency bin | `St, λ₁, λ₂, ..., λ_nblk` |
| `mea*.f00001` | Temporal mean field | Nek5000 field file |
| `pod*.f00001` | POD modes | One field file per saved mode (`modal_nsave`) |
| `dm1*`, `dm2*` | DMD modes: real part, imaginary part | One field file per saved mode, numbered in the order of `dmd_spectrum.dat` (`dm11cyl0.f00003` is mode 3) |
| `sRe*`, `sIm*` | SPOD modes: real part, imaginary part | The leading `modal_nsave` modes of the frequency bins 1, 1 + `nfreq/8`, 1 + 2 `nfreq/8`, ... and of the last bin; the peak bin is not always among them |

`SPG*` files hold the sponge-zone function (`forcing.f90`), not SPOD modes.

#### Plotting Results

Use the provided Python script:

```bash
cd example/cylinder_re100/600_modal_pod
python plot_modal.py           # Auto-detect and plot all
python plot_modal.py --pod     # POD only
python plot_modal.py --dmd     # DMD only
python plot_modal.py --spod    # SPOD only
```

Generates: `pod_spectrum.png`, `dmd_spectrum.png`, `spod_spectrum.png`, `pod_fft_spectrum.png`, `spectral_comparison.png`

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

### Mode 6.4: DMT Filter Bank

`userParam01 = 6.4` (or `nekstab_mode = 'dmt_bank'`) runs the DNS without any force and passes the flow through a set of filters, one per chosen frequency. Each filter returns the spatial structure of the flow at its frequency as a complex field `û(x)`, with `u(x,t) ≈ Re[û(x) e^{iωt}]`. It is the filter of mode 1.3 used for analysis: the filters watch the flow and do not act on it. Unlike POD, DMD and SPOD it needs no snapshot files; it runs in the time loop and writes the result at the last step.

**Filter.** A filter of angular frequency `ω` and width `β` is a damped oscillator driven by the flow `u`, point by point:

```
p' = q,    q' = u − β q − ω² p,    b = β q
```

`b` is `u` filtered around `ω`: gain 1 and phase 0 at `ω`, gain 0 at frequency 0. The complex amplitude is `û = (b − i b'/ω) e^{−iωt}`, where `b' = β (u − β q − ω² p)` follows from the equations, so no data are differentiated. The simpler combination `β (q + iωp)` has the same value in the band but passes a constant mean flow with the gain `iβ/ω`; for a mean of 1 and a tone of 0.01 it gave an error of 5 in units of the tone amplitude, the form above gave 2e-3 (the start-up transient). `ω = 0` gives the mean flow: `b = β q` follows `u` with the time constant `1/β`. The oscillator is integrated exactly over each step for an input that varies linearly between the steps (matrix exponential, `dmt_filter.f90`), so a variable time step needs no special case. Every active field is filtered: velocity, and pressure and scalars when they are written.

**Registration**, in `nekStab_usrchk`, on every rank:

```fortran
use nekstab_dmt, only: dmt_bank_add, dmt_bank_options
...
call dmt_bank_add(omega1, beta1)          ! angular frequency 2*pi*f; beta is optional, default 0.05*omega
call dmt_bank_add(2.0d0*omega1)
call dmt_bank_add(0.0d0, 0.1d0*omega1)    ! the mean flow, averaged over 1/beta
call dmt_bank_options(hist_step=100, snapshots=.false., avg_from=t_avg)
```

At most 8 filters; `beta` must lie below `omega`. The filters start from the field at the start of the run, in the equilibrium of that field, so a mean flow starts no transient.

**Width and run length.** The start-up transient decays as `exp(−β t/2)`: after `2 ln(1000)/β = 13.8/β` time units it is below 1e-3 of the amplitude. A narrow filter separates close frequencies (a tone at distance `Δω` appears with about `β/Δω` of its amplitude) and settles slowly. The start-up prints the settling time of each filter and warns when the run is shorter.

**Leakage and averaging.** A filter also passes part of the other frequencies of the flow, about `β/ω` of their amplitude. In a flow with harmonics of one frequency this appears as a ripple of `û` at the difference of the frequencies. With `avg_from` the files hold the average of `û` over `[avg_from, end of run]`, and over a whole number of periods of the fundamental the ripple cancels. On a flow made of three harmonics (amplitudes 0.10, 0.03, 0.01) on a mean, at `β = 0.1 ω`, the amplitude at the last step differs from the exact one by 5 %, 15 % and 45 % of the amplitude of each harmonic, and the average over 10 periods by 4e-5, 5e-5 and 1e-4. Without `avg_from`, expect a ripple of this size.

**Files.** The number `k` follows the order of the calls of `dmt_bank_add`.

| File | Contents |
|------|----------|
| `bRe<session>0.f0000k`, `bIm<session>0.f0000k` | `Re û` and `Im û` of filter `k`, at the last step or averaged over the window; the time stamp is the end of the run |
| `dmt_bank.dat` | every `hist_step` steps: `time, ‖u‖, ‖u − Σ b_k‖, ‖û_1‖, ‖û_2‖, ...`; the norm is the Krylov-vector norm, and `‖û_k‖` does not depend on the phase |
| `b01<session>0.f<n>`, `b02...` | `b_k` at every output step of Nek5000, with `snapshots=.true.` |

`‖u − Σ b_k‖` is the part of the flow that no filter holds. The filter states are not written to disk: a restarted run starts its filters again.

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
| `k_dim` | 100 | Largest Krylov basis. The run stops at the first Arnoldi step where the `schur_tgt` leading values have converged and have stood still for 5 steps, so a larger `k_dim` costs memory, not matvecs; the basis fills and restarts only when they need more steps |
| `use_cgs` | .true. | Use CGS2 orthogonalization with batched BLAS/MPI reductions (`.false.` = classic MGS + reorthogonalization) |
| `eigen_tol` | 1e-6 | Convergence tolerance: a Ritz value has converged when its residual is below `eigen_tol*max(abs(mu), 1)`. The same test applies to the residual of each exported mode vector. A leading value also has to stand still: it must move by less than the same bound in each of the last 5 steps before the run stops on it |
| `schur_tgt` | 2 | Number of wanted eigenvalues, counted from the top of the list sorted by `abs(mu)` (a conjugate pair counts as two). The run stops when the first `schur_tgt` entries have converged, not when any `schur_tgt` values have. `schur_tgt <= 0` runs one full pass of `k_dim` steps. A restart keeps at least `schur_tgt + 2` and at most `k_dim - schur_tgt` Ritz values, so the next Arnoldi run adds at least `schur_tgt` steps |
| `schur_stop` | 0 | Stop rule. `0`: the run stops when the first `schur_tgt` entries of the sorted list have converged. `1`: it stops when any `schur_tgt` values have converged. Use `1` below the critical Reynolds number, where slow non-oscillating modes can lead the list and converge much later than the modes of interest. Set it in `nekStab_usrchk` with `call bcast(schur_stop, isize)` |
| `schur_del` | 0.10 | Selection radius of the restart: it keeps the Ritz values with `abs(mu) >= 1 - schur_del`, and those with a Schur residual below `sqrt(eigen_tol)` |
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
| `ifquasilaminar` | .false. | Quasilaminar RANS stability: velocity-only perturbation, base eddy viscosity prescribed (requires `iffindiff`; pair with `ifKEnorm`) |
| `ifKEnorm` | .false. | Velocity-only kinetic-energy inner product (turbulence scalars excluded) |


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

BoostConv accelerates the velocity only. With a temperature or another scalar it prints a warning and leaves the scalar unaccelerated; use SFD for those cases. On the 2D cylinder (`example/cylinder_re100/120_baseflow_boostconv`) the defaults reach the gate 1e-6 at Re = 50 and not at Re = 100 within 500 time units. A subspace of 30 is not faster than the default there. The stage README has the table.

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

**Period of a Floquet run.** The period is the time stamp of the base-flow file `BF_<session>0.f00001`, not `endTime`. The first sweep integrates the base flow over this period from the start field and prints

```
 Orbit closure: ||U(T)-U(0)||/||U(0)|| =  8.356E-06  for T =  8.81517222497400
```

A Newton-converged orbit closes to 1e-5 or below. Above 1e-3 the run prints a warning, because the multipliers then belong to a path that does not close: the period of the time stamp is wrong, the polynomial order differs from that of the orbit, or the base flow is not on the orbit. With a period 9 % too long, the orbit of `flip_flop_re62` closed to 7.7e-3. The run also prints the converged multiplier closest to 1 and warns when none lies within 1e-3, but a multiplier near 1 does not prove the period: with that wrong period the converged list still held two multipliers within 1e-11 of 1 whose eigenvectors sat on a few mesh vertices.

**Explicit velocity filter.** `filtering = explicit` acts while the base flow advances: in a DNS, in SFD and in the first sweep of a Floquet run. The steady Newton run, the steady direct and adjoint runs and the later sweeps of a Floquet run (they replay the stored orbit) apply no velocity filter. A steady state from a filtered DNS or SFD run is therefore not an exact fixed point of the unfiltered map that Newton and the stability runs use: its first Newton residual is larger than the DNS suggests (`lid_driven_re3600`: 3e-5 for the state that the DNS had left at 7e-7 under the filter), and Newton moves it to the fixed point of the unfiltered map. With `filtering = explicit` the first sweep of a Floquet run therefore differs from the later ones. The Floquet stages of the examples run without the filter, except `cubic_cavity_re1950`, whose orbit comes from a filtered DNS.

#### Sponge Zone

| Parameter | Default | Description |
|-----------|---------|-------------|
| `xLspg`, `xRspg` | 0.0 | Left/right sponge boundaries in x |
| `yLspg`, `yRspg` | 0.0 | Bottom/top sponge boundaries in y |
| `zLspg`, `zRspg` | 0.0 | Front/back sponge boundaries in z |
| `spng_st` | 0.0 | Sponge strength (0 = disabled) |
| `acc_spg` | 0.333 | Acceleration ramp fraction |

**Sponge reference field.** The sponge of a nonlinear run (DNS, SFD, Newton) forces the velocity toward `spng_vr`, the field that the run holds when the sponge starts: the restart file, or the initial condition of a cold start. A perturbation run forces toward zero. A steady Newton run with the sponge on therefore solves the Navier-Stokes equations plus a forcing toward its start field, so the base flow depends on that field inside the sponge zone, and the forcing vanishes only if the start field is the solution. For `naca0012_re2000` (8 ranks, Jean Zay), Newton from the DNS snapshot differs by 1.1e-2 (relative l2 norm of the velocity) from the base flow that the stability stages ship: by 7e-4 to 8e-4 for x < 1 and by 3e-2 for 5 < x < 6, inside the sponge. A second Newton run seeded with the output of the first moves the field by 3.5e-3 and moves the leading growth rate of the direct problem from -0.0643 to -0.0661; the shipped base flow gives -0.0665. Periodic-orbit and Floquet runs use no sponge, because the sponge would force the orbit toward one snapshot.

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
| `ifdmtbank` | .false. | Enable the DMT filter bank (mode 6.4), see Mode 6.4 |
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
│   │── ks_select.f90             # Ritz value selection of the Krylov-Schur restart
│   │── dmt_filter.f90            # Oscillator-pair filter of the DMT bank (coefficients)
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
| `nekstab_eigensolvers` | `eigensolvers.f90` | `krylov_schur`, `inner_product`, `norm`, `outpost_ks`, `schur_condensation`, `select_eigenvalues` |
| `nekstab_vectors` | `nek_vectors.f90` | `nopcopy`, `nopadd2`, `nopadd2s2`, `nopsub2`, `nopsub3`, `nopcmult`, `noprzero`, `nopaxpby`, `axpby`, `opadd3`, `opaddcol3` |
| `nekstab_lapack` | `lapack_wrapper.f90` | `schur`, `ordschur`, `eig`, `sort_eigendecomp`, `select_eigvals`, `lstsq`, `eig_symmetric`, `eig_hermitian` |
| `nekstab_argsort` | `argsort.f90` | `argsort` |
| `nekstab_ks_select` | `ks_select.f90` | `ks_select_eigenvalues`, `ks_ensure_conjugate_pairs`, `ks_converged_indices` |
| `nekstab_dmt_filter` | `dmt_filter.f90` | `dmt_filter_step_coeffs`, `dmt_filter_amp_coeffs` |
| `nekstab_dmt` | `dmt.f90` | `dmt`, `dmt_init`, `dmt_bank_add`, `dmt_bank_options`, `dmt_bank` |
| `nekstab_matvec` | `matvec.f90` | `prepare_linearized_solver`, `matvec`, `forward_linearized_map`, `forward_finite_difference_map`, `adjoint_linearized_map` |
| `nekstab_io` | `IO.f90` | `whereyouwant`, `read_eigenvalue`, `load_mode_pair`, `load_files`, `k_load` |
| `nekstab_newton` | `newton_krylov.f90` | `newton_krylov`, `ts_gmres`, `initialize_gmres_vector`, `nonlinear_forward_map`, `set_nek5000_tolerances`, `spec_tole` |
| `nekstab_fixedpoint` | `fixedp.f90` | `sfd`, `BoostConv`, `tdf`, `boostconv_core`, `qr_dec`, `linear_system` |
| `nekstab_sensitivity` | `sensitivity.f90` | `wave_maker`, `bf_sensitivity`, `ts_steady_force_sensitivity`, `biorthogonalize`, `delta_forcing`, `animate_mode_only`, `animate_mode`, `animate_mode_Floquet` |
| `nekstab_energy_budget` | `energy_budget.f90` | `stability_energy_budget`, `stability_energy_budget_floquet`, `compute_velocity_gradient_tensor`, `compute_dissipation`, `compute_production` |
| `nekstab_diagnostics` | `diagnostics.f90` | `nekStab_energy`, `nekStab_enstrophy`, `outpost_vort`, `norm_grad`, `smooth_field`, `nekStab_outpost`, `nekStab_comment`, `nekStab_printNEKParams` |
| `nekstab_vortex` | `vortex.f90` | `vortex_core`, `compute_omega_jc`, `compute_omega`, `compute_q`, `compute_delta`, `compute_swirling` |
| `nekstab_statistics` | `statistics.f90` | `nekStab_avg` |
| `nekstab_reynolds` | `reynolds_force.f90` | `reynolds_commit`, `reynolds_load`, `reynolds_add`, `reynolds_enable`, `reynolds_enabled`, `reynolds_armed` |
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

### Layout

`example/` holds one folder per flow, such as `cylinder_re100` or `cubic_cavity_re1950`. Each flow folder holds numbered stages, for example `cylinder_re100/310_stability_direct`. The first digit of the number names the kind of stage:

| First digit | Kind of stage | Examples |
|---|---|---|
| 0 | DNS | `000_dns` |
| 1 | base flow by time stepping | SFD, BoostConv, DMT, TDF, mean flow |
| 2 | base flow by Newton | `210_baseflow_newton` |
| 3 | stability | direct, adjoint, Floquet, transient growth |
| 4 | post-processing | wavemaker, energy budget, mode animation |
| 5 | OTD | `500_otd` |
| 6 | modal decomposition | POD, DMD, SPOD, DMT filter bank |

A stage folder is a complete case: `.par`, `.usr`, `SIZE`, the mesh (`.re2`, `.ma2`), `run.local.slurm`, a `README.md`, a `plot.py` and a `ref/` folder with `reference.json` and the figures. The field files (base flows, restart files, modes) are not tracked. A stage that needs the field of an earlier stage says so in its README, and its run script copies the file from the sibling folder. [`EXAMPLES.md`](EXAMPLES.md) gives the rationale of each kind of stage, the case catalog and the verified results.

### Available Cases

| Directory | Flow |
|-----------|------|
| `cylinder_re100/` | 2D cylinder wake, Re = 100: the full pipeline from DNS to modal decomposition |
| `cylinder_re180/` | Periodic cylinder wake, Re = 180: Newton for the orbit, Floquet, OTD |
| `cylinder_re30_thermal/` | Heated cylinder wake, Re = 30 |
| `cylinder_re1m/` | RANS cylinder wake, Re = 1e6 (experimental, not validated) |
| `moving_cylinder_re100/` | Cylinder in forced oscillation on a moving mesh, Re = 100 |
| `flip_flop_re62/` | Side-by-side cylinders, Re = 62 |
| `tpjet_re2000/` | Harmonically forced jet, Re = 2000 |
| `thermosyphon_ra500/` | Buoyancy-driven annular thermosyphon, Ra = 500 |
| `back_fstep_re500/` | Backward-facing step, Re = 500 |
| `lid_driven_re3600/` | 2D lid-driven cavity, Re = 3600 |
| `cubic_cavity_re1914/`, `cubic_cavity_re1920/`, `cubic_cavity_re1950/` | 3D lid-driven cubic cavity below, just above and above the first Hopf bifurcation |
| `naca0012_re2000/` | NACA 0012 airfoil, Re = 2000 |
| `poiseuille_re5k/` | Channel flow, Re = 5000 |
| `slot_fst_re495/` | Slot jet with free-stream turbulence |

### Running a Stage

```bash
cd example/cylinder_re100/000_dns
mks 1cyl                                    # compile; the argument is the name of the .par and .usr
sbatch run.local.slurm                      # or: SLURM_NTASKS=8 bash run.local.slurm
python3 plot.py                             # figures of the run
python3 ../../../scripts/check_against_ref.py .   # compare the run with ref/reference.json
```

The README of the flow gives the order of the stages. The stage `210_baseflow_newton` of a flow starts from a DNS field, the stability stages from the base flow that it writes, the wavemaker from the direct and the adjoint modes, and so on. Each stage needs at least `lpmin` ranks (the value in its `SIZE`).

### Plot Utilities

Shared plotting helpers live in `example/nekplot.py`: field input, interpolation, residual histories, spectra and OTD diagnostics. The `plot.py` of a stage writes `plot_*.png`. The figures that the validation page shows are the files of the `ref/` folder of the stage.

---

## Validation

### Checks

Each stage has its own check:

```bash
python3 scripts/check_against_ref.py example/cylinder_re100/310_stability_direct/direct
```

`check_against_ref.py` reads `ref/reference.json` of the stage, extracts each listed quantity from the files of the run (a column of a data file, or a statistic of a field file) and compares it with the stored value within the stored relative or absolute tolerance. It fails when a quantity is outside its tolerance, when the log of the run shows a NaN or an error, and when the reference lists no quantity. The stored values are regression values of this repository unless the reference names a literature or analytical source.

Further checks:
- `validation/build_gallery.py` writes `validation/index.html`, the page of references and figures ([nekstab.github.io/nekStab/validation](https://nekstab.github.io/nekStab/validation/)). The test `validation/tests/test_gallery_complete.py` fails for a stage without a reference, summary, description, figure or README, and for a page that is out of date.
- `scripts/compiler_battery.py` builds and runs stages under nine compiler profiles and tabulates the spread of the reference quantities.
- `validation/tests/` holds the unit tests of the Ritz value selection, the DMT filters, the mode selector and the file readers (`pytest`).

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

### Moving Mesh

With `motion = user` in `[MESH]` the `.usr` sets the mesh velocity **w** (`wx`, `wy`, `wz`) in `userchk` and in
`nekStab_usrchk`, and Nek5000 solves the base flow on the moving mesh in the arbitrary Lagrangian-Eulerian form

```
d(B u)/dt = B div(w u) - B (u·∇)u - ∇p + ν A u + B f
```

where B is the mass matrix of the current mesh. The Nek5000 routine `admeshv` adds `B div(w u)` to the base flow. The
perturbation solver of Nek5000 has no matching term, and the mesh velocity does not depend on the flow, so nekStab adds
`B div(w u')` to the forcing of every perturbation (`ale_mesh_velocity_term` in `src/forcing.f90`, called from
`nekStab_forcing`). With it the linearized map is the derivative of the nonlinear map: the direct modes (0.1, 2.x,
3.1 and 3.11) run on a moving mesh. The linearized maps call `nekStab_usrchk` before each step, so it has to set **w**
with the same routine as `userchk`, and a periodic mesh motion has to close after one period for the Newton and Floquet
modes (`example/moving_cylinder_re100/`). The adjoint, transient-growth and OTD modes have no mesh-velocity term and
stop at start-up on a moving mesh.

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

2. **Ritz extraction**: Eigenvalues of H_m approximate eigenvalues of A. The residual of the Ritz pair
   (mu_i, V_m y_i) is |h_{m+1,m} y_i(m)|. It is tested after every step, and the run stops at the first step m where
   the `schur_tgt` leading values (largest `abs(mu)`) have converged and have stood still for 5 steps.

   The second condition is needed because the residual bounds the error of the eigenvalue only through its condition
   number kappa: |theta - mu| is up to kappa times the residual. For a strongly non-normal operator (a convectively
   unstable flow) kappa is large. For the forced jet `tpjet_re2000` kappa = 4e5 (from the direct and the adjoint mode,
   `<y,x>` = 2.5e-6 for unit norms); the adjoint run stopped at step 37 with the residual 1.1e-6 and its multiplier
   was off by 5e-3, the direct run at step 44 with 8e-7 and 2e-4. While the basis still resolves the mode the value
   moves from step to step by a part of the error that is left, so a value that moved by less than `eigen_tol` (relative
   for `abs(mu) > 1`) in each of the last 5 steps has settled. At the last step of a pass the residual test decides
   alone, as before, and a run whose value has not settled by then prints a warning.
   Check: `validation/tests/test_ritz_stop_rule.py` runs the rule on random non-normal operators (condition number of
   the leading value up to 1e6, clusters of slow values, a matvec exact to 1e-9). A stop counts as worse when its error
   is above 3 times the error of the whole basis plus `eigen_tol`. Of 900 operators the residual test alone was worse
   in 62 %, one step of rest in 26 %, three steps in 1.2 %, five steps in none; the mean stop moved from step 17 to
   step 26 of 128. A residual below the tolerance does not replace a comparison with a larger `k_dim` or with the
   adjoint value when `kappa` is large. The error of the matvec limits the eigenvalue too: it is of the order of `kappa`
   times the relative error of the matvec, whatever `k_dim` and `eigen_tol` are (5e-4 for `tpjet_re2000`, whose
   matvec repeats to 1e-9).

3. **Implicit restart**: Schur decomposition filters unwanted eigenvalues. It happens only when the basis reaches
   `k_dim` before the leading values have converged.

4. **Deflation**: Lock converged eigenvalues, continue with reduced subspace

5. **Check of the exported modes**: the Ritz residual equals ||M x - mu x|| only while the relation of step 1 holds
   for every column. Each mode to export is therefore given to the linearized solver once more (one matvec for a real
   mode, two for a complex pair) and its residual ||M x - mu x|| / ||x|| is measured. Below `eigen_tol` the mode is
   written. Above it, but less than `eigen_tol` above the Ritz residual, it is written with a warning: the run stops
   as soon as the Ritz residual is below the tolerance, and a small error of the relation can put the vector just
   above it. A larger difference means the vector is not an eigenvector: it is not written, and when it is one of
   the `schur_tgt` leading modes the run ends with an error after the files are written. Row k of
   `Spectre_*_conv.dat` and mode file k always hold the same mode. A locked mode shows a residual near 1e-16 in the
   spectrum files; its vector residual is the Ritz residual it had when it was locked. `Spectre_?.info` lists both
   numbers under `[vector residuals]`, with `k_used=` (size of the last basis) and `check matvecs=`.

### Adjoint Operator

The continuous adjoint satisfies:

```
-∂u†/∂t - (U·∇)u† + (∇U)^T·u† = -∇p† + (1/Re)∇²u†
```

Discrete adjoint uses the **transpose** of the linearized time-stepper.

### Adjoint on an Axisymmetric Mesh

With `axiSymmetry = yes` the mass matrix `bm1` of Nek5000 contains the radius. The adjoint term `(∇U)^T u†` is computed
by `convop_adj` (`core/perturb.f`) on the dealiased mesh: it integrates the product against the basis functions with
the weights of `set_dealias_rx` (`core/convect2.f`), which have no radius, and divides by `bm1`. The term is then too
large by the factor 1/r. On the forced jet (`example/tpjet_re2000/321_stability_adjoint_floquet`) the pressure residual
of the adjoint run grew by a factor of about 1.4 per time step after the first ten steps and reached 1e39 within one
forcing period, and the run returned multipliers of modulus 1e52. The direct perturbation
equation is not affected: for an axisymmetric mesh Nek5000 selects the pointwise convection routine `conv1d`
(`reader_par.f`: `if (ifaxis) param(99) = 3`).

nekStab corrects the term in the forcing of every adjoint perturbation on an axisymmetric mesh
(`axisymmetric_adjoint_term` in `src/forcing.f90`, called from `nekStab_forcing`, the same pattern as the mesh-velocity
term above). The routine computes `(∇U)^T u†` again with the radius in the integration weights and returns the
difference between the term of Nek5000 and this term. `makeufp` multiplies the forcing by B, and `advabp_adjoint`
subtracts B times the term of Nek5000, so what remains is B times the correct term. A Nek5000 that has the radius in
its weights gives a difference of zero, and the correction then does nothing. The weights are built in the routine, not
taken from `set_dealias_rx`, so that the result does not depend on that.

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
- A seed from a filtered DNS can start at a residual 40 times that of the DNS (see the explicit velocity filter under Base Flow), and the residual can stay flat or rise for three iterations in a row before it falls (`lid_driven_re3600`: 14 iterations in all). The loop ends only after six iterations without progress.

**Run stops with `k_normalize non-finite norm` in a Newton or Krylov mode**
- Set `residualProj = no` in `[VELOCITY]`. Each matvec starts from an unrelated
  vector, so the velocity residual projection has no history to use, and on
  `naca0012_re2000` its Helmholtz solve produced a non-finite norm at 24 or more
  ranks on Jean Zay (ifort). With `residualProj = no` the Newton run gave identical
  residuals at 8, 16, 24 and 32 ranks. `nekStab` prints a note at start-up when
  the projection is on in these modes.

**A converged eigenvalue at mu = 1.0000000 whose mode is a set of spikes**
- Seen in versions up to v2.0.0-rc4 (v1.0.1 included), most often below the
  critical Reynolds number and on many MPI ranks with few elements per rank.
  The Ritz residual is small and the mode can even pass the vector check, but
  the field is not a flow: the values that neighbouring elements hold for the
  same grid point differ by O(1). The physical modes sit further down the list
  and are often not exported.
- Cause. Nek5000 stores a point on an element face once for each element
  that holds it. Copies on different ranks pick up different round-off in
  two places. Where 3 or more ranks share a point, the gather-scatter adds
  the contributions in a different order on each rank (1e-16 in a base
  flow). And dgemm, which rotates the Krylov basis at each restart, rounds
  the last rows of a rank's array in another kernel than the rows before
  them: with MKL on an AMD Zen 5 CPU the last `mod(nv, 12)` rows, so on a
  rank of 92 elements of 64 points the last 8 points of its last element,
  and none on a rank of 93. The perturbation solver adds a
  continuous increment to the old field, so a difference between copies is
  never removed: the linearized map multiplies it by exactly 1. Below the
  critical Reynolds number every physical multiplier is below 1, and the
  Arnoldi method brings this round-off to the top of the spectrum after
  enough matvecs. More ranks give more such points: a 2D case of 5920
  elements on 64 ranks has 76, the cylinder of the examples on 8 ranks has 2.
- Fix. `matvec` now gives every copy of a shared point the mean of its copies
  on each vector that goes into the solver (`average_copies` in
  `src/matvec.f90`), as Nek5000 does for every field it
  reads. A difference between copies then has multiplier 0 instead of 1, and
  the vector check refuses any mode that still holds one.
- Results of older versions. `bin/c0check FIELD...` reports the largest jump
  between copies, relative to the field maximum, and the number of jumps at
  points shared by 1, 2, 3 or more ranks. Exit status 1 means a broken file.
  Base flows and continuous modes are correct. A broken mode cannot be
  repaired: Nek5000 averages a file when it reads it, but the average of a
  jump is not a mode, and the physical mode was never written. Recompute it.
  Do not restart a Krylov basis (`KRY*`, `HES*`) from such a run, because the
  averaged vectors no longer fit the stored Hessenberg matrix.

**Run stops with `Cannot open history file in subroutine hpts()`**
- `hpts` reads the probe list from `<session>.his` and writes the probe values
  after it, in the same file. The file name ends in `.his`, which `.gitignore`
  excludes, so a fresh clone does not have it. Each `run.local.slurm` of the
  examples writes the probe list before the run. In your own case, create the
  file: the first line holds the number of points, then one line per point with
  `x y z` (`x y` in 2D).

**History points hold zeros**
- Set `lhis` in `SIZE` to a value that does not divide the number of probes. When
  the number of probes is a multiple of `lhis`, `hpts` leaves the last block of
  probes on a rank that the output routine does not read, and the file receives
  zeros for them. On `cylinder_re100`, 2 probes with `lhis = 2` wrote zeros;
  with `lhis = 3` the same run wrote the values. `nekStab` prints a warning at
  start-up when it finds this case.

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
