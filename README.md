# nekStab
[![Build Status](https://travis-ci.com/ricardofrantz/nekStab.svg?token=DpocmcBgXShNTZ9nAQ5y&branch=master)](https://travis-ci.com/ricardofrantz/nekStab) Global stability framework for Nek5000

- uparam(1) = mode (0: DNS, 1:stabilization, 3:stability)

- uparam(2) = optional restart vector for stability

- uparam(3) = stabilzation technique (1:SFD)

- uparam(4) = frequency cuttor of forcing frequency

- uparam(5) = gain or forcing amplitude 

- uparam(6) = 

- uparam(7) = 

- uparam(8) = 

- uparam(9) = 

- uparam(10): sponge strenght (>0 to activate)

  

### Main features 



## OUR STABILITY FRAMEWORK


uparam(1) = 0 - > Direct numerical simulation

uparam(1) = 3.- > Direct mode computation

uparam(1) = 3.2 - > Adjoint mode computation

uparam(1) = 3.3 - > Direct-adjoint mode computation for transient growth analysis


# INSTALLING AND RUNNING

Fist install code dependencies if running with GCC

```bash
sudo apt-get -y install libmpich-dev mpich libopenblas-dev build-essential cmake m4
```

Clone the repository and Nek5000 as a submodule:

```bash
git clone --depth=50 --branch=master https://github.com/ricardofrantz/nekStab.git
cd nekStab
git submodule update --init --recursive
```

Run _*(sudo) vim $HOME/.bashrc*_ and add the following:

```bash
export NEKSTAB_SOURCE_ROOT=$HOME/nekStab
export NEK_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT/Nek5000
export PATH=$NEK_SOURCE_ROOT/bin:$PATH
```

Compile Nek5000 tools:

```bash
cd ~/nekStab/Nek5000/tools
./maketools all
```

<details>
  <summary>Outputs</summary>

  ```javascript
building genmap ... done
building gencon ... done
building genbox ... done
building n2to3 ... done
building reatore2 ... done
building nekmerge ... done
building prenek ... done
building postnek ... done
building nekamg_setup ... done
building gmsh2nek ... done
building exo2nek ... done
building cgns2nek ... done
  ```
</details>



Move to a specific example, for instance: _*cd examples/1cyl*_

```bash
./cmpile.sh all
nekbmpi 1cyl 4
```



# ADAPTING YOUR EXISTING CASE TO NEKSTAB

we recommend converting old _.rea_ files to __.par__ + __.re2__ and __.ma2__

in *SIZE* be sure to modify and include:

```fortran
parameter (lpelt=lelt)
include 'NEKSTAB.inc'
```

in your *.usr*, add to __userchk__ and __userf __respectively:

```
call nekStab
call nekStab_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
```



# PROGRAMMING BEST PRACTICES

> The best style is consistent style!

1. Always use **implicit none** with exception from _.usr_
2. Comments with ! instead of c
3. Use tab length with 3 space characters
4. small caps for fortran syntax : while, for, do, etc. 
5. proper define real numbers: 0.0d0, 0.50d0, 1.0d0
6. subroutine -> return, end 
7. to print output we use write(6,*)
8. stop the code with _call **exitt**_
9. use *nelv* instead of _*nelt*_ for consistency as we must respect velocity mesh
10. use *n* and *n2* for loops instead of *ntot1* 
11. local field declarations with *lt* and *lt2*

```fortran
real, save, intent(inout), dimension(lt) :: field1,field2,field3
```

11. shared variables in our custom x_SIZE.f  respecting the standard Nek style like in /core

```fortran
c     Solution data
      real vx     (lx1,ly1,lz1,lelv)
     $    ,vy     (lx1,ly1,lz1,lelv)
     $    ,vz     (lx1,ly1,lz1,lelv)
      common /vptsol/ vx, vy, vz
```

12. Limit line length increased to 132 to avoid excessive use of continuation tabs, in most cases we are slightly over the 72 characters standard
