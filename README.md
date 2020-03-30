# nekStab
Global stability framework for Nek5000

uparam(1)=
uparam(2)=
uparam(3)=
uparam(4)=
uparam(5)=
uparam(6)=
uparam(7)=
uparam(8)=
uparam(9)=
uparam(10)=

0: dns only
0.1 turbulent statistics (full 3d (jet-crossflow), 2 averged (cylinder and boundary layer), 3: channel flow
0.2 DMT (space-time Fourier) - not sure how can be useful

when steady state we can dynamically extract: which norm to ensure steadystate?
0.1 POD proper orthogonal decomposition (space-only)
0.2 SPOD spectral proper orthogonal decomposition (space-time, single frequency modes)
0.3 DMD dynamic mode decomposition  (https://arxiv.org/pdf/1406.7187.pdf) 

STABILITY FRAMEWORK

1: forcing
1.1 SFD
1.2 BOOSTCONV
1.3 TDF


3: stability
3.1 direct
3.2 adjoint
3.3 direct-adoint (for transient growth)
3.4 resolvent (for optimal perturbation)

BASE FLOW OPTIONS
a) prescribed in useric (e.g. Blasius solution)
b) load fixed point BF_ + casename + 0.f00001 
c) compute with 1.1 or 1.2, dumpt to disk and switch to 3 -> requires the frequency of the mode to T


PERIODIC ORBIT WITH FORCED SYMMETRY (1 cylinder and boundary layer) force vz=0 in userchk
a) run DNS alongside 3,  
b) run DNS, export PER and Fourier, reload and 3


PERIODIC ORBIT OPTIONS FOR FORCED CASES -> uparam(4) is the forced frequency
a) computed with TDF, export sequence and Fourier decomposition, switch to 3
b) load PER or Fourier modes and run 3

PERIODIC ORBIT FOR UNFORCED CASES (2 cyl and Jet-in-Crossflow)
























!!!!!!!!!!!!!!!!!!!!!! 
PROGRAMMING BEST PRACTICES: The best style is consistent style!
function and subroutine when possible
1) always use implicit none with exception from .usr file
2) comments with ! instead of c
3) one tab is 3 spaces
4) use small caps
5) proper define real numbers: 1 -> 1. -> 1.0d0, 0.0d0
6) subroutine -> return, end
7) use nelv or nelt? we should pick only one
8) to print output we use write(6,*)
9) stop the code with call exitt


local field declarations with
real, save, dimension(lt) :: field1,field2,field3,

shared variables:

in SIZE we follow the standard from vx,vy,vz

c     Solution data
      real vx     (lx1,ly1,lz1,lelv)
     $    ,vy     (lx1,ly1,lz1,lelv)
     $    ,vz     (lx1,ly1,lz1,lelv)
      common /vptsol/ vx, vy, vz

to discuss:
11) Limit line length increased to 132 to avoid excissive use of continuation tabs, modern editors can do soft wrap, bust most cases we are slightly over the 72 characters standard
5. File and folder structur: clustering subroutine in files with a clean structure 
6. Naming conventions





