# nekStab
Global stability framework for Nek5000



- uparam(1) =
- uparam(2) =
- uparam(3) =
- uparam(4) =
- uparam(5) =
- uparam(6) =
- uparam(7) =
- uparam(8) =
- uparam(9) =
- uparam(10) =



### Main features 

0:  dns only
0.1 turbulent statistics (full 3d (jet-cross-flow), 2 averaged (cylinder and boundary layer), 3: channel flow
0.2 DMT (space-time Fourier) - not sure how can be useful

1. when steady state we can dynamically extract: which norm to ensure steady-state?
   0.1 POD proper orthogonal decomposition (space-only)
   0.2 SPOD spectral proper orthogonal decomposition (space-time, single frequency modes)
   0.3 DMD dynamic mode decomposition  (https://arxiv.org/pdf/1406.7187.pdf)

   

## STABILITY FRAMEWORK

### Forcing 

1.1 SFD 
1.2 BOOSTCONV
1.3 TDF



sponge: which parameter to activate

stratification: which parameter to activate coupling between T



### Stability
3.1 direct
3.2 adjoint
3.3 direct-adjoint (for transient growth)
3.4 resolvent (for optimal perturbation)

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









# PROGRAMMING BEST PRACTICES

> The best style is consistent style!
> Function and subroutines when possible...

1. Always use **implicit none** with exception from _.usr_

2. Comments with ! instead of c
3. Tav lenght 
4. small caps for fortran syntax : for, do, return
5. proper define real numbers: 0.0d0, 0.50d0, 1.0d0
6. subroutine -> return, end 
7. use *nelv* or *nelt*? we should pick only one
8. to print output we use write(6,*)
9. stop the code with _call **exitt**_

local field declarations with

```fortran
real, save, dimension(lt) :: field1,field2,field3
```

shared variables in SIZE we follow the standard nek style like in /core

```fortran
c     Solution data
      real vx     (lx1,ly1,lz1,lelv)
     $    ,vy     (lx1,ly1,lz1,lelv)
     $    ,vz     (lx1,ly1,lz1,lelv)
      common /vptsol/ vx, vy, vz
```

to discuss:
11) Limit line length increased to 132 to avoid excessive use of continuation tabs, modern editors can do soft wrap, bust most cases we are slightly over the 72 characters standard

5. Naming conventions
