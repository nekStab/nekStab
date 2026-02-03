```
                 __   _____  __          __
   ____   ___   / /__/ ___/ / /_ ____ _ / /_
  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \
 / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /
/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/
COPYRIGHT (c) 2020-2025 DynFluid Laboratoire Paris
```

[![CI Build](https://github.com/nekStab/nekStab/actions/workflows/ci.yml/badge.svg?branch=RANS_update)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml)

**nekStab** is a toolbox for performing bifurcation analysis using the spectral element CFD solver [Nek5000](https://github.com/Nek5000/Nek5000).
It is released under the BSD-3-Clause license.

The project started in 2010 with the PhD thesis of [Jean-Christophe Loiseau](https://loiseaujc.github.io/) and builds on the work of former PhD students of our group such as Frédéric Alizard, Stefania Cherubini, Alessandro Bucci, Mirko Farano, Francesco Picella and [Ricardo Frantz](https://github.com/ricardofrantz).
Ricardo is the one who brought all the previous contributions together in a single toolbox.

It is actively maintained mainly by [Ricardo Frantz](https://github.com/ricardofrantz) and [Jean-Christophe Loiseau](https://loiseaujc.github.io/).
Both of them work at [DynFluid](http://dynfluid.ensam.eu/), a fluid dynamics laboratory part of [Arts et Métiers Institute of Technology](https://artsetmetiers.fr/en).

## First Steps

**nekStab** is a toolbox written in 'Fortran 90' for the spectral element solver Nek500.
If you already have C and Fortran compilers, you can install both on Ubuntu/Debian distributions with the following commands.

**Prerequisites**

Linux (GCC)
```bash
sudo apt-get -y install build-essential gfortran libmpich-dev libopenblas-dev cmake m4 htop
```

macOS
```bash
brew install mpich gfortran wget git cmake htop
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

## Development

**nekStab** is maintained by [Ricardo Frantz](https://github.com/ricardofrantz).

> **See also:** [LightKrylov](https://github.com/nekStab/LightKrylov) is a modern Fortran library providing abstract Krylov methods. [neklab](https://github.com/nekStab/neklab) is the next-generation bifurcation and stability analysis toolbox for Nek5000, built on LightKrylov. Both are actively developed by Jean-Christophe Loiseau and Simon Kern.

## Communication

- Mail: [Jean-Christophe Loiseau](mailto:loiseau.jc@gmail.com?subject=[GitHub]%20Information%20about%20nekStab) or [Ricardo Frantz](mailto:rasfrantz@gmail.com?subject=[GitHub]%20Information%20about%20nekStab)
- Website: [https://nekstab.github.io/](https://nekstab.github.io/)

### Citation

When using **nekStab**, please cite [Frantz et al. (2023)](https://doi.org/10.1115/1.4056808). See [CITATIONS.md](CITATIONS.md) for BibTeX entries and additional references.
