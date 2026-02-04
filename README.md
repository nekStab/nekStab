```
                 __   _____  __          __
   ____   ___   / /__/ ___/ / /_ ____ _ / /_
  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \
 / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /
/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/
```

| OS | Arch | Compiler | MPI | Status |
|:---|:----:|:---------|:----|:------:|
| Ubuntu 24.04 | x86_64 | gfortran 14 | OpenMPI / MPICH | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=dev&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| Ubuntu 24.04 | x86_64 | ifort 2024 / ifx 2025 | Intel MPI | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=dev&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| macOS 26 | ARM64 | gfortran 14 | OpenMPI | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=dev&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |

**nekStab** is a toolbox for global stability and bifurcation analysis using the spectral element solver [Nek5000](https://github.com/Nek5000/Nek5000). Released under BSD-3-Clause license.

## Features

### Steady-State & Periodic Orbit Computation
| Method | Description |
|--------|-------------|
| **SFD** | Selective Frequency Damping for unstable steady states |
| **BoostConv** | Residual acceleration for slow convergence |
| **TDF** | Time-Delayed Feedback for periodic orbits |
| **Newton-Krylov** | Quadratic convergence for fixed points and UPOs |

### Global Stability Analysis
| Analysis | Steady Flows | Time-Periodic (Floquet) |
|----------|:------------:|:-----------------------:|
| **Direct eigenmodes** | ✓ | ✓ |
| **Adjoint eigenmodes** | ✓ | ✓ |
| **Transient growth** | ✓ | ✓ |

### Sensitivity & Receptivity
- **Wavemaker** — structural sensitivity to feedback
- **Base flow sensitivity** — response to mean flow modifications
- **Forcing response** — optimal and localized forcing analysis

### Modal Decomposition
| Method | Description |
|--------|-------------|
| **POD** | Proper Orthogonal Decomposition (energy-ranked modes) |
| **DMD** | Dynamic Mode Decomposition (frequency-ranked modes) |
| **SPOD** | Spectral POD (frequency-resolved coherent structures) |

### Advanced Capabilities
- **OTD modes** — real-time Lyapunov vectors for chaotic flows
- **Linearized DNS** — perturbation evolution around base flows
- **Scalar transport** — temperature and passive scalar stability
- **Sponge zones** — non-reflecting boundaries for open flows
- **Vortex identification** — λ₂, Q-criterion, Ω-criterion output

### Why nekStab?
- **Matrix-free**: No Jacobian storage — scales to millions of DoFs
- **Spectral accuracy**: Leverages Nek5000's high-order elements
- **MPI parallel**: Efficient on laptops to supercomputers
- **Validated**: Benchmarked against canonical flows (cylinder, cavity, jets)

## Quick Reference

| Mode | String | Description |
|:----:|--------|-------------|
| 0 | `'dns'` | Direct Numerical Simulation |
| 0.1 | `'linear_dns'` | Linearized DNS (perturbation) |
| 1.1 | `'sfd'` | Selective Frequency Damping |
| 1.2 | `'boostconv'` | BoostConv acceleration |
| 1.4 | `'tdf'` | Time-Delayed Feedback |
| 2.0 | `'newton_fp'` | Newton for fixed points |
| 2.1 | `'newton_po'` | Newton for periodic orbits |
| 3.1 | `'direct'` | Direct stability eigenmodes |
| 3.11 | `'floquet_direct'` | Floquet direct analysis |
| 3.2 | `'adjoint'` | Adjoint stability eigenmodes |
| 3.21 | `'floquet_adjoint'` | Floquet adjoint analysis |
| 3.3 | `'transient_growth'` | Optimal perturbations |
| 3.31 | `'floquet_tg'` | Floquet transient growth |
| 4.1 | `'energy_budget'` | Kinetic energy budget |
| 4.2 | `'wavemaker'` | Structural sensitivity |
| 4.3 | `'bf_sensitivity'` | Base flow sensitivity |
| 5 | `'otd'` | Optimally Time-Dependent modes |
| 6.1 | `'pod'` | Proper Orthogonal Decomposition |
| 6.2 | `'dmd'` | Dynamic Mode Decomposition |
| 6.3 | `'spod'` | Spectral POD |

See [DOC.md](DOC.md) for full parameter reference (includes animation and forcing modes).

## Examples

Ready-to-run cases in `example/`:

| Case | Physics | Demonstrated Features |
|------|---------|----------------------|
| `cylinder/` | 2D/3D wake | DNS, Newton, stability, Floquet, wavemaker |
| `back_fstep/` | Separation | Convective instability, transient growth |
| `lid_driven/` | Confined | Steady bifurcations |
| `blasius/` | Boundary layer | Tollmien-Schlichting waves |
| `torus/` | Curved pipe | Dean instability, 3D modes |
| `thersyphon/` | Convection | Pitchfork + Hopf bifurcations |

## First Steps

**Prerequisites**

Linux (GCC)
```bash
sudo apt -y install build-essential gfortran libmpich-dev libopenblas-dev libfftw3-dev cmake m4 htop
```

macOS
```bash
brew install mpich gfortran fftw wget git cmake htop
```

**Cloning the repository and Nek5000**

```bash
git clone --depth=50 https://github.com/nekStab/nekStab.git
cd nekStab
./Nek5000setup.sh
```

Run **vim $HOME/.bashrc** and add the following :
```bash
export NEKSTAB_SOURCE_ROOT=$HOME/nekStab
export NEK_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT/Nek5000
export PATH=$NEK_SOURCE_ROOT/bin:$PATH
export PATH=$NEKSTAB_SOURCE_ROOT/bin:$PATH
ulimit -s unlimited
ulimit -c unlimited
```

## Compilation

Go to a given example folder and compile the code:
```bash
cd ~/nekStab/example/cylinder/baseflow/newton
mks 1cyl
```

### Compiler Selection

nekStab supports multiple Fortran compilers. The build script auto-detects available compilers in this order: `ifx` → `ifort` → `gfortran`. To force a specific compiler:

```bash
NEKSTAB_FC=ifx mks 1cyl    # Intel LLVM (recommended for Intel/AMD CPUs)
NEKSTAB_FC=ifort mks 1cyl  # Intel Classic (available on many HPC systems)
NEKSTAB_FC=gcc mks 1cyl    # GCC/gfortran
```

### Performance Considerations

| Compiler | Performance | Notes |
|----------|-------------|-------|
| **ifx** (Intel LLVM) | Fastest | Requires Intel oneAPI 2024+, uses MKL |
| **ifort** (Intel Classic) | Fast | Available on HPC systems, discontinued in oneAPI 2025 |
| **gfortran** (GCC) | Good | Universal, uses system BLAS/LAPACK |

> **Note:** Intel discontinued `ifort` in oneAPI 2025, but it remains available on many HPC systems. New installations should use `ifx`.

The Intel compilers use MKL (Math Kernel Library) which provides highly optimized BLAS/LAPACK routines with runtime CPU dispatching.

### Debug Mode

For debugging crashes or numerical issues:
```bash
mks 1cyl --debug
```
This enables stack traces, bounds checking, and additional warnings.

## Running

After successful compilation, run on 4 processors:
```bash
nekbmpi 1cyl 4
```
to follow the code output in the _logfile_ try:
```bash
tail -f logfile
```
To stop the code just:
```bash
killall nek5000
```

For more information, see the [Documentation](DOC.md).

## Mode Selection

Three equivalent ways to select operating mode:

| Method | Where | Example |
|--------|-------|---------|
| String | `.usr` | `nekstab_mode = 'floquet_adjoint'` |
| Flags | `.usr` | `isAdjoint = .true.` + `ifFloquet = .true.` |
| uparam | `.par` | `userParam01 = 3.21` |

## Development

**nekStab** is maintained by [Ricardo Frantz](https://github.com/ricardofrantz).

> **See also:** [LightKrylov](https://github.com/nekStab/LightKrylov) is a modern Fortran library providing abstract Krylov methods. [neklab](https://github.com/nekStab/neklab) is the next-generation bifurcation and stability analysis toolbox for Nek5000, built on LightKrylov. Both are actively developed by Jean-Christophe Loiseau and Simon Kern.

## Communication

- Mail: [Ricardo Frantz](mailto:rasfrantz@gmail.com?subject=[GitHub]%20Information%20about%20nekStab)
- Website: [https://nekstab.github.io/](https://nekstab.github.io/)

### Citation

When using **nekStab**, please cite [Frantz et al. (2023)](https://doi.org/10.1115/1.4056808). See [CITATIONS.md](CITATIONS.md) for BibTeX entries and additional references.
