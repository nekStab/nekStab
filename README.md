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

**Features:**
- Base flow computation (SFD, BoostConv, Newton-Krylov)
- Linear stability analysis (direct and adjoint eigenmodes)
- Floquet analysis for time-periodic flows
- Transient growth and optimal perturbations
- Sensitivity analysis (wavemaker, forcing response)
- Matrix-free Krylov methods, scales to millions of DoFs

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

## Mode Configuration

nekStab supports three equivalent ways to select the operating mode:

**Method 1: String mode** (human-readable, recommended)
```fortran
! In nekStab_usrchk subroutine in your .usr file:
nekstab_mode = 'floquet_adjoint'
```

**Method 2: Flag mode** (flexible)
```fortran
! In nekStab_usrchk subroutine:
isAdjoint = .true.
ifFloquet = .true.   ! modifier flag
```

**Method 3: uparam** (automation-friendly, backward compatible)
```
# In .par file [PROBLEMTYPE] section:
userParam01 = 3.21
```

All three methods produce identical behavior. See [DOC.md](DOC.md) for the complete mode reference table.

## Development

**nekStab** is maintained by [Ricardo Frantz](https://github.com/ricardofrantz).

> **See also:** [LightKrylov](https://github.com/nekStab/LightKrylov) is a modern Fortran library providing abstract Krylov methods. [neklab](https://github.com/nekStab/neklab) is the next-generation bifurcation and stability analysis toolbox for Nek5000, built on LightKrylov. Both are actively developed by Jean-Christophe Loiseau and Simon Kern.

## Communication

- Mail: [Jean-Christophe Loiseau](mailto:loiseau.jc@gmail.com?subject=[GitHub]%20Information%20about%20nekStab) or [Ricardo Frantz](mailto:rasfrantz@gmail.com?subject=[GitHub]%20Information%20about%20nekStab)
- Website: [https://nekstab.github.io/](https://nekstab.github.io/)

### Citation

When using **nekStab**, please cite [Frantz et al. (2023)](https://doi.org/10.1115/1.4056808). See [CITATIONS.md](CITATIONS.md) for BibTeX entries and additional references.
