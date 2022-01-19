[![GCC](https://github.com/ricardofrantz/nekStab/actions/workflows/gcc.yml/badge.svg?branch=master)](https://github.com/ricardofrantz/nekStab/actions/workflows/gcc.yml)
[![Intel](https://github.com/ricardofrantz/nekStab/actions/workflows/intel.yml/badge.svg?branch=master)](https://github.com/ricardofrantz/nekStab/actions/workflows/intel.yml)
[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://ricardofrantz.github.io/nekStabDoc/en/master/)

**nekStab** is a toolbox to conduct bifurcation analyses using the spectral element CFD solver [Nek5000](https://github.com/Nek5000/Nek5000).
It is released under the ??? licence.

The project started in 2010 with the PhD thesis of [Jean-Christophe Loiseau](https://loiseaujc.github.io/) and is build on top of the work of former PhD students of our group such as Frédéric Alizard, Stefania Cherubini, Alessandro Bucci, Mirko Farano, Francesco Picella and [Ricardo Frantz](https://github.com/ricardofrantz).
Ricardo is the one responsible for having brought together all the previous contributions into a single toolbox.

It is actively maintained primarily by [Ricardo Frantz](https://github.com/ricardofrantz) and [Jean-Christophe Loiseau](https://loiseaujc.github.io/).
Both of them work at [DynFluid](http://dynfluid.ensam.eu/), a fluid dynamics laboratory part of [Arts et Métiers Institute of Technology](https://artsetmetiers.fr/en).

### News

- **February 2022:** The first release of `nekStab` is available online!
It is released under the [???]() licence.
Do not hesitate to get in touch with us or check the [documentation]() if you want to know more.

## Getting started

**nekStab** is a toolbox written in `Fortran 90` for the spectral element solver Nek500.
Assuming you already have C and Fortran compilers, both of them can be installed on Ubuntu/Debian distros using the following commands.

**Prerequisites**

```bash
sudo apt-get -y install libmpich-dev libopenblas-dev cmake m4
```

**Cloning the repository**

```bash
git clone --depth=50 --branch=master https://github.com/ricardofrantz/nekStab.git
cd nekStab
git submodule update --init --recursive
```

Run `vim $HOME/.bashrc` and add the following :

```bash
export NEKSTAB_SOURCE_ROOT=$HOME/nekStab
export NEK_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT/Nek5000
export PATH=$NEK_SOURCE_ROOT/bin:$PATH
```

Computing the fixed point for the cylinder flow example using the Newton-Krylov solver on 4 processors is as simple as

```bash
cd ~/nekStab/examples/cylinder/1_2baseflow/newton
./cmpile.sh all
nekbmpi 1cyl 4
```

More information about compiling the code on Mac OS or optional packages is available in the [**Documentation**](https://ricardofrantz.github.io/nekStabDoc/en/master/index.html).

## Development

TBA

## Help and Support

TBA

### Documentation

- HTML documentation (stable realse): https://ricardofrantz.github.io/nekStabDoc/en/master/

### Communication

- Mail : [Jean-Christophe Loiseau](mailto:loiseau.jc@gmail.com?subject=[GitHub]%20Information%20about%20nekStab) or [Ricardo Frantz](mailto:rasfrantz@gmail.com?subject=[GitHub]%20Information%20about%20nekStab)
- Website : [https://nekstab.github.io/](https://nekstab.github.io/)

### Citation

If you use **nekStab**, please consider citing one of the following papers :
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
