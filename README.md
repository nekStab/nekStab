# nekStab
[![Build Status](https://travis-ci.com/ricardofrantz/nekStab.svg?token=DpocmcBgXShNTZ9nAQ5y&branch=master)](https://travis-ci.com/ricardofrantz/nekStab) Global stability framework for Nek5000

- uparam(1) = mode (0: DNS, 1:SD, 2:SA, 3:D-A, 4:Res)

- uparam(2) = stability restart starting vector 

- uparam(3) = forcing mode (1:SFD, 2:boosconv, 3:TDF)

- uparam(4) = frequency input

- uparam(5) = gain input or forcing amplitude 

- uparam(6) = tripping position 

- uparam(7) = 

- uparam(8) = 

- uparam(9) = 

- uparam(10): sponge intensity (>0 to activate)

  

### Main features 



## OUR STABILITY FRAMEWORK



uparam(1) = 0 - > Direct numerical simulation

uparam(1) = 1 - > Direct mode computation

uparam(1) = 2 - > Adjoint mode computation

uparam(1) = 3 - > Direct-adjoint mode computation for transient growth analysis

uparam(1) = 4 - > Resolvent analysis for optimal perturbation





0.1 turbulent statistics (full 3d (jet-cross-flow), 2 averaged (cylinder and boundary layer), 3: channel flow
0.2 DMT (space-time Fourier) - not sure how can be useful

1. when steady state we can dynamically extract: which norm to ensure steady-state?
   0.1 POD proper orthogonal decomposition (space-only)
   0.2 SPOD spectral proper orthogonal decomposition (space-time, single frequency modes)
   0.3 DMD dynamic mode decomposition  (https://arxiv.org/pdf/1406.7187.pdf)

   

### Forcing 

1.1 SFD 
1.2 BOOSTCONV
1.3 TDF

stratification: which parameter to activate coupling between T



### Stability


#### STEADY FLOWS: MEAN OR BASE

a) prescribed in useric (e.g. Blasius solution)
b) load steady state as:  BF_ + casename + 0.f00001
c) compute with 1.1 or 1.2, dumpt to disk and switch to 3 -> requires the frequency of the mode to T

 ####  LIMIT CYCLES WITH FORCED SYMMETRY

 ####   (1 cylinder and boundary layer) force vz=0 in userchk

a) run DNS alongside 3,  
b) run DNS, export PER and Fourier, reload and 3

####  LIMIT CYCLES OPTIONS FOR FORCED CASES

####   -> uparam(4) is the forced frequency

a) computed with TDF, export sequence and Fourier decomposition, switch to 3
b) load PER or Fourier modes and run 3

####  LIMIT CYCLES FOR UNFORCED CASES (2 cyl and Jet-in-Crossflow)



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
13. 