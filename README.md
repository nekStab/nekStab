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

### News

- **February 2022:** The first official version of 'nekStab' is available online!
It is released under the BSD-3-Clause license.
Do not hesitate to contact us, have a look at the [documentation](https://nekstab.github.io/nekStabDoc/en/master/) or read the [release notes](https://github.com/nekStab/nekStab/blob/master/RELEASE.md) if you want to know more.

## First Steps

**nekStab** is a toolbox written in 'Fortran 90' for the spectral element solver Nek500.
If you already have C and Fortran compilers, you can install both on Ubuntu/Debian distributions with the following commands.

**Prerequisites**

Linux
```bash
sudo apt-get -y install libmpich-dev libopenblas-dev cmake m4 htop
```
MacOS
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

nekStab supports multiple Fortran compilers. The build script auto-detects available compilers in this order: `ifort` → `ifx` → `gfortran`. To force a specific compiler:

```bash
NEKSTAB_FC=ifx mks 1cyl    # Intel LLVM (recommended for Intel CPUs)
NEKSTAB_FC=ifort mks 1cyl  # Intel Classic
NEKSTAB_FC=gcc mks 1cyl    # GCC/gfortran
```

### Performance Considerations

| Compiler | Performance | Notes |
|----------|-------------|-------|
| **ifx** (Intel LLVM) | Fastest | Requires Intel oneAPI 2024+, uses MKL |
| **ifort** (Intel Classic) | Fast | Legacy, being phased out by Intel |
| **gfortran** (GCC) | Good | Universal, uses system BLAS/LAPACK |

**For best performance on Intel/AMD x86 CPUs:**
1. Install [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) (free)
2. Source the environment: `source /opt/intel/oneapi/setvars.sh`
3. Build with ifx: `NEKSTAB_FC=ifx mks 1cyl`

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

For more information on compiling the code on Mac OS or optional packages, see the [**Documentation**](https://nekstab.github.io/nekStabDoc/en/master/).

## Development

**nekStab** is mainly developed by Jean-Christophe Loiseau and Ricardo Frantz.
However, we welcome contributors of all levels of experience.
If you are planning a larger contribution, we encourage you to discuss the concept here on GitHub and to interact with us regularly to ensure your efforts are targeted.

## Help and support

Although **nekStab** relies on Nek5000, none of us are active developers of Nek5000.
If you have questions about Nek5000, please contact the dedicated [GitHub repo](https://github.com/Nek5000/Nek5000)
 and [documentation](https://nek5000.github.io/NekDoc/)

### Communication

- Mail : [Jean-Christophe Loiseau](mailto:loiseau.jc@gmail.com?subject=[GitHub]%20Information%20about%20nekStab) or [Ricardo Frantz](mailto:rasfrantz@gmail.com?subject=[GitHub]%20Information%20about%20nekStab)
- Website : [https://nekstab.github.io/](https://nekstab.github.io/)

### Citation

When using **nekStab**, please consider citing one of the following papers :

- [Frantz et al. (2023)](https://arxiv.org/abs/2301.12940) presents the complete theoretical framework underlying **nekStab**.
```bibtex
@article{frantz2023krylov,
    author = {Frantz, R. A. S. and Loiseau, J.-Ch. and Robinet, J.-Ch.},
    title = "{Krylov Methods for Large-Scale Dynamical Systems: Application in Fluid Dynamics}",
    journal = {Applied Mechanics Reviews},
    volume = {75},
    number = {3},
    year = {2023},
    month = {03},
    issn = {0003-6900},
    doi = {10.1115/1.4056808},
    url = {https://doi.org/10.1115/1.4056808},
    note = {030802},
    eprint = {https://asmedigitalcollection.asme.org/appliedmechanicsreviews/article-pdf/75/3/030802/6996354/amr\_075\_03\_030802.pdf},
}
```

- [Loiseau et al. (2019)](https://arxiv.org/pdf/1804.03859.pdf) presents most of the theoretical framework underlying **nekStab**.
```bibtex
@incollection{chapter:loiseau:2019,
  title={Time-stepping and Krylov methods for large-scale instability problems},
  author={Loiseau, J.-Ch. and Bucci, M. A. and Cherubini, S. and Robinet, J.-Ch.},
  booktitle={Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics},
  pages={33--73},
  year={2019},
  publisher={Springer}
}
```
- [Loiseau et al. (*J. Fluid Mech.*, 2014)](https://sam.ensam.eu/bitstream/handle/10985/8974/DYNFLUID-JFM-LOISEAU-2014.pdf?sequence=1&isAllowed=y) describes the first implementation of the Arnoldi solver in Nek5000.
```bibtex
@article{jfm:loiseau:2014,
    title={Investigation of the roughness-induced transition: global stability analyses and direct numerical simulations},
    author={Loiseau, J.-Ch. and Robinet, J.-Ch. and Cherubini, S. and Leriche, E.},
    journal={J. Fluid Mech.},
    volume={760},
    pages={175--211},
    year={2014},
}
```
