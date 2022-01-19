[![GCC](https://github.com/ricardofrantz/nekStab/actions/workflows/gcc.yml/badge.svg?branch=master)](https://github.com/ricardofrantz/nekStab/actions/workflows/gcc.yml)
[![Intel](https://github.com/ricardofrantz/nekStab/actions/workflows/intel.yml/badge.svg?branch=master)](https://github.com/ricardofrantz/nekStab/actions/workflows/intel.yml)
[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://ricardofrantz.github.io/nekStabDoc/en/master/)

**nekStab** is a toolbox to conduct bifurcation analyses using the spectral element CFD solver [Nek5000](https://github.com/Nek5000/Nek5000).
It is released under the ??? licence.

The project started in 2010 with the PhD thesis of [Jean-Christophe Loiseau](https://loiseaujc.github.io/) and is build on top of the work of former PhD students of our group such as Frédéric Alizard, Stefania Cherubini, Alessandro Bucci, Mirko Farano, Francesco Picella and [Ricardo Frantz](https://github.com/ricardofrantz).
Ricardo is the one responsible for having brought together all the previous contributions into a single toolbox.

It is actively maintained primarily by [Ricardo Frantz](https://github.com/ricardofrantz) and [Jean-Christophe Loiseau](https://loiseaujc.github.io/).
Both of them work at [DynFluid](http://dynfluid.ensam.eu/), a fluid dynamics laboratory part of [Arts et Métiers Institute of Technology](https://artsetmetiers.fr/en).

# Getting started

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
