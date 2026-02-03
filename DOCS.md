# nekStab Documentation

A comprehensive toolbox for bifurcation analysis using the spectral element CFD solver Nek5000.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Compilation](#compilation)
4. [Quick Start](#quick-start)
5. [Operating Modes](#operating-modes)
6. [Parameter Reference](#parameter-reference)
7. [Examples](#examples)
8. [Theoretical Background](#theoretical-background)
9. [Troubleshooting](#troubleshooting)
10. [Citation](#citation)

---

## Overview

**nekStab** provides tools for:

- **Direct Numerical Simulation (DNS)** - Time integration of Navier-Stokes equations
- **Base Flow Computation** - Finding steady states via SFD, BoostConv, TDF, Newton-Krylov
- **Linear Stability Analysis** - Eigenvalue problems using Krylov-Schur methods
- **Sensitivity Analysis** - Wavemaker, structural sensitivity, optimal forcing
- **Floquet Analysis** - Stability of time-periodic flows
- **Transient Growth** - Non-modal stability analysis

### Architecture

```
nekStab/
├── core/           # Main Fortran 90 modules
│   ├── main.f90           # Entry point and mode dispatcher
│   ├── eigensolvers.f90   # Krylov-Schur eigensolver
│   ├── newton_krylov.f90  # Newton-GMRES for fixed points/UPOs
│   ├── fixedp.f90         # SFD, BoostConv, TDF methods
│   ├── sensitivity.f90    # Sensitivity and wavemaker analysis
│   ├── postproc.f90       # Post-processing routines
│   ├── matvec.f90         # Matrix-vector products (linearized NS)
│   ├── krylov_subspace.f90    # Krylov vector type definitions
│   └── krylov_decomposition.f90  # Arnoldi/Lanczos decompositions
├── bin/            # Build scripts (mks, nekbmpi)
├── examples/       # Test cases
└── Nek5000/        # Nek5000 solver (submodule)
```

---

## Installation

### Prerequisites

**Linux (Ubuntu/Debian)**
```bash
sudo apt-get install build-essential gfortran libopenmpi-dev liblapack-dev libblas-dev cmake
```

**macOS**
```bash
brew install gcc open-mpi cmake
```

**Optional: Intel oneAPI (recommended for performance)**
```bash
# Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html
source /opt/intel/oneapi/setvars.sh
```

### Clone and Setup

```bash
git clone https://github.com/nekStab/nekStab.git
cd nekStab
./Nek5000setup.sh
```

### Environment Variables

Add to `~/.bashrc` or `~/.zshrc`:

```bash
export NEKSTAB_SOURCE_ROOT=$HOME/nekStab
export NEK_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT/Nek5000
export PATH=$NEK_SOURCE_ROOT/bin:$NEKSTAB_SOURCE_ROOT/bin:$PATH
ulimit -s unlimited
```

---

## Compilation

### Basic Usage

```bash
cd examples/cylinder/dns
mks 1cyl
```

### Compiler Selection

The `mks` script auto-detects compilers: `ifort` → `ifx` → `gfortran`

Force a specific compiler:
```bash
NEKSTAB_FC=ifx mks 1cyl    # Intel LLVM (fastest on x86)
NEKSTAB_FC=ifort mks 1cyl  # Intel Classic
NEKSTAB_FC=gcc mks 1cyl    # GCC/gfortran
```

### Performance Comparison

| Compiler | Performance | BLAS/LAPACK | Notes |
|----------|-------------|-------------|-------|
| **ifx** | Fastest | MKL | Intel oneAPI 2024+, AVX-512 optimized |
| **ifort** | Fast | MKL | Legacy, being deprecated |
| **gfortran** | Good | System BLAS | Universal compatibility |

### Debug Mode

```bash
mks 1cyl --debug
```

Enables: stack traces, bounds checking, floating-point traps

### Clean Build

```bash
mks 1cyl --fresh
```

---

## Quick Start

### 1. Run DNS

```bash
cd examples/cylinder/dns
mks 1cyl                    # Compile
nekbmpi 1cyl 4              # Run on 4 processors
tail -f logfile             # Monitor output
```

### 2. Compute Base Flow (Newton-Krylov)

```bash
cd examples/cylinder/baseflow/newton
# Edit 1cyl.par: userParam01 = 2.0
mks 1cyl
nekbmpi 1cyl 4
```

### 3. Linear Stability Analysis

```bash
cd examples/cylinder/stability/direct
# Edit 1cyl.par: userParam01 = 3.1
mks 1cyl
nekbmpi 1cyl 4
```

---

## Operating Modes

Set the mode via `userParam01` in the `.par` file:

### Mode 0: DNS

| Value | Description |
|-------|-------------|
| `0` | Standard DNS |
| `0.1` | Linearized DNS (perturbation evolution) |

### Mode 1: Fixed Point Methods

| Value | Method | Description |
|-------|--------|-------------|
| `1.1` | SFD | Selective Frequency Damping |
| `1.2` | BoostConv | Residual acceleration method |
| `1.3` | DMT | Dynamic Mode Tracking |
| `1.4` | TDF | Time-Delayed Feedback |

### Mode 2: Newton-Krylov

| Value | Description |
|-------|-------------|
| `2.0` | Newton-GMRES for fixed points |
| `2.1` | Newton-GMRES for UPOs (Unstable Periodic Orbits) |
| `2.2` | Newton-GMRES for forced UPOs |

### Mode 3: Eigenvalue Problems

| Value | Description |
|-------|-------------|
| `3.1` | Direct LNSE (Linearized Navier-Stokes) |
| `3.11` | Direct LNSE with Floquet |
| `3.2` | Adjoint LNSE |
| `3.21` | Adjoint LNSE with Floquet |
| `3.3` | Transient growth |
| `3.31` | Transient growth with Floquet |

### Mode 4: Post-Processing

| Value | Description |
|-------|-------------|
| `4.0` | All (budget + wavemaker + sensitivity) |
| `4.1` | Kinetic energy budget |
| `4.2` | Wavemaker |
| `4.3` | Base flow sensitivity |
| `4.41` | Sensitivity to steady force (type 1) |
| `4.42` | Sensitivity to steady force (type 2) |
| `4.43` | Delta forcing |
| `4.50` | Animate mode (direct) |
| `4.51` | Animate mode + base flow deformation |
| `4.52` | Animate Floquet mode |

### Mode 5: OTD

Optimally Time-Dependent modes for chaotic flows.

---

## Parameter Reference

### .par File Parameters

```ini
[GENERAL]
userParam01 = 3.1    # Operating mode (see above)
userParam02 = ...    # Method-specific
userParam03 = ...    # Method-specific
...
```

### nekStab Parameters (set in .usr file)

Set these in `nekStab_usrchk` subroutine:

#### Krylov Solver

| Parameter | Default | Description |
|-----------|---------|-------------|
| `k_dim` | 100 | Krylov subspace dimension |
| `eigen_tol` | 1e-6 | Eigenmode convergence tolerance |
| `schur_tgt` | 2 | Schur target for factorization |
| `schur_del` | 0.10 | Schur shift parameter |
| `maxmodes` | 20 | Max converged modes to save |

#### Newton-Krylov

| Parameter | Default | Description |
|-----------|---------|-------------|
| `findiff_order` | 1 | Finite difference order for Jacobian |
| `epsilon_base` | 1e-6 | Perturbation scale for finite differences |

#### BoostConv

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bst_skp` | 10 | Skip iterations between updates |
| `bst_snp` | 10 | Residual subspace matrix size |

#### Output Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifres` | .false. | Output restart files (KRY*, HES*) |
| `ifvor` | .false. | Output vorticity fields |
| `ifvox` | .false. | Output vortex criteria (Q, λ2, ω) |
| `glob_skip` | 10 | Energy computation frequency |

#### Initial Conditions

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifseed_nois` | .true. | Use noise as initial seed |
| `ifseed_symm` | .false. | Use symmetric initial seed |
| `ifseed_load` | .false. | Load initial seed from file |

#### Base Flow

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifldbf` | .true. | Load base flow for stability |
| `ifbf2D` | .false. | Force 2D base flow |
| `ifstorebase` | .true. | Store base flow for Floquet |

#### Sponge Zone

| Parameter | Default | Description |
|-----------|---------|-------------|
| `xLspg, xRspg` | 0.0 | Sponge zone x boundaries |
| `yLspg, yRspg` | 0.0 | Sponge zone y boundaries |
| `zLspg, zRspg` | 0.0 | Sponge zone z boundaries |
| `spng_st` | 0.0 | Sponge strength (0 = disabled) |
| `acc_spg` | 0.333 | Acceleration phase fraction |

#### OTD Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifotd` | .false. | Enable OTD computation |
| `otd_printStep` | 100 | OTD output frequency |
| `otd_gsStep` | 10 | Gram-Schmidt frequency |
| `otd_FTLEPeriod` | 0.0 | FTLE averaging period |

---

## Examples

### Available Test Cases

| Directory | Description |
|-----------|-------------|
| `cylinder/` | Flow around circular cylinder (Re=50-100) |
| `back_fstep/` | Backward-facing step |
| `lid_driven/` | Lid-driven cavity |
| `blasius/` | Blasius boundary layer |
| `cubic_cavity/` | 3D cubic cavity |
| `naca0012/` | NACA 0012 airfoil |
| `torus/` | Toroidal pipe flow |
| `poiseuille/` | Plane Poiseuille flow |

### Cylinder Example Structure

```
cylinder/
├── dns/                  # Direct numerical simulation
├── baseflow/
│   ├── sfd/             # SFD base flow
│   ├── newton/          # Newton-Krylov base flow
│   └── boostconv/       # BoostConv base flow
├── stability/
│   ├── direct/          # Direct stability
│   └── adjoint/         # Adjoint stability
└── postprocessing/
    ├── wavemaker/       # Wavemaker analysis
    └── sensitivity/     # Sensitivity analysis
```

### Typical Workflow

1. **Generate mesh** (using Nek5000 tools)
   ```bash
   genbox < box.txt
   genmap
   ```

2. **Run DNS** to verify setup
   ```bash
   # userParam01 = 0
   nekbmpi case 4
   ```

3. **Compute base flow**
   ```bash
   # userParam01 = 2.0 (Newton) or 1.1 (SFD)
   nekbmpi case 4
   ```

4. **Run stability analysis**
   ```bash
   # userParam01 = 3.1 (direct) or 3.2 (adjoint)
   nekbmpi case 4
   ```

5. **Post-process**
   ```bash
   # userParam01 = 4.2 (wavemaker)
   nekbmpi case 4
   ```

---

## Theoretical Background

### Linearized Navier-Stokes Equations (LNSE)

For a base flow **U**, the perturbation **u'** evolves according to:

```
∂u'/∂t = -U·∇u' - u'·∇U - ∇p' + ν∇²u'
∇·u' = 0
```

### Eigenvalue Problem

Seeking solutions **u'** = **û** exp(λt), we solve:

```
λû = Aû
```

where **A** is the linearized Navier-Stokes operator.

### Krylov-Schur Method

nekStab uses time-stepping to build a Krylov subspace:

```
Km = span{q, Aq, A²q, ..., A^(m-1)q}
```

The Arnoldi decomposition `AV = VH + residual` gives Ritz values approximating eigenvalues of **A**.

### Adjoint Equations

The adjoint eigenproblem provides sensitivity information:

```
λ*û† = A†û†
```

### Wavemaker

The structural sensitivity (wavemaker) identifies regions where perturbations most affect eigenvalues:

```
S(x) = |û(x)| · |û†(x)|
```

---

## Troubleshooting

### Common Issues

**Compilation fails with module errors**
```
Cannot open module file 'krylov_subspace.mod'
```
→ Clean and rebuild: `mks case --fresh`

**MKL not found (Intel compiler)**
```
cannot find -lmkl_intel_lp64
```
→ Source Intel environment: `source /opt/intel/oneapi/setvars.sh`

**Segmentation fault with ifx**
→ Ensure you're using the latest nekStab with ifx compatibility fixes

**Eigenvalues don't converge**
→ Increase `k_dim` or decrease `eigen_tol`

**Base flow doesn't converge**
→ Try different method (SFD for oscillating flows, Newton for steady)
→ Check Reynolds number is in stable regime
→ Adjust SFD parameters (χ, Δ)

### Performance Tips

1. **Use Intel compilers** on x86 for 20-30% speedup
2. **Increase polynomial order** (`lx1`) for accuracy, not elements
3. **Use appropriate time step** - CFL ~0.5 for stability analysis
4. **Parallelize wisely** - ~1000-5000 elements per MPI rank is optimal

---

## Citation

If you use nekStab, please cite:

```bibtex
@article{frantz2023krylov,
    author = {Frantz, R. A. S. and Loiseau, J.-Ch. and Robinet, J.-Ch.},
    title = "{Krylov Methods for Large-Scale Dynamical Systems: Application in Fluid Dynamics}",
    journal = {Applied Mechanics Reviews},
    volume = {75},
    number = {3},
    year = {2023},
    doi = {10.1115/1.4056808},
}
```

Additional references:

- Loiseau et al. (2019) - Time-stepping and Krylov methods, Springer
- Loiseau et al. (2014) - First Arnoldi implementation in Nek5000, J. Fluid Mech.

---

## License

BSD-3-Clause License. See [LICENSE](LICENSE) file.

## Contact

- Ricardo Frantz: rasfrantz@gmail.com
- Jean-Christophe Loiseau: loiseau.jc@gmail.com
- Website: https://nekstab.github.io/
