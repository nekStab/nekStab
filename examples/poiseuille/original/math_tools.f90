!> @file math_tools.f
!! @ingroup math
!! @brief Set of math related tools for KTH modules
!! @author Adam Peplinski
!! @date Jan 31, 2017
!=======================================================================
!> @brief Step function
!! @ingroup math
!! @details Continuous step function:
!!  \f{eqnarray*}{
!!    stepf(x) = \left\{ \begin{array}{ll}
!!  0 &\mbox{ if $x \leq x_{min}$} \\
!!  \left(1+e^{\left((x-1)^{-1} + x^{-1}\right)}\right)^{-1} &\mbox{ if $x \leq x_{max}$} \\
!!  1 &\mbox{ if $x >  x_{max}$}
!!       \end{array} \right.
!!  \f}
!!  with \f$ x_{min} = 0.02\f$ and \f$ x_{max}=0.98\f$
!! @param[in] x       function argument
!! @return mth_stepf
      real function mth_stepf(x)
      implicit none

      ! argument list
      real x

      ! local variables
      real xdmin, xdmax
      parameter (xdmin = 0.001, xdmax = 0.999)
!-----------------------------------------------------------------------
      ! get function vale
      if (x.le.xdmin) then
         mth_stepf = 0.0
      else if (x.le.xdmax) then
         mth_stepf = 1./( 1. + exp(1./(x - 1.) + 1./x) )
      else
         mth_stepf = 1.
      end if

      end function mth_stepf
!=======================================================================
!> @brief Give random distribution depending on position
!! @ingroup math
!! @details The original Nek5000 rundom number generator is implementted
!!  in @ref ran1. This totally ad-hoc random number generator below
!!  could be preferable to the origina one for the simple reason that it
!!  gives the same initial cindition independent of the number of
!!  processors, which is important for code verification.
!! @param[in] ix,iy,iz     GLL point index
!! @param[in] ieg          global element number
!! @param[in] xl           physical point coordinates
!! @param[in] fcoeff       function coefficients
!! @return  random distribution
      real function mth_rand(ix,iy,iz,ieg,xl,fcoeff)
      implicit none

      include 'SIZE'
      include 'INPUT'       ! IF3D

      ! argument list
      integer ix,iy,iz,ieg
      real xl(LDIM)
      real fcoeff(3)
!-----------------------------------------------------------------------
      mth_rand = fcoeff(1)*(ieg+xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy +
     $     fcoeff(3)*ix
      if (IF3D) mth_rand = fcoeff(1)*(ieg +xl(NDIM)*sin(mth_rand)) +
     $     fcoeff(2)*iz*ix + fcoeff(3)*iz*iy
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = cos(mth_rand)

      return
      end function mth_rand

!=======================================================================
!     A simple portable random number generator
!
!     Requires 32-bit integer arithmetic
!     Taken from Numerical Recipes, William Press et al.
!     gives correlation free random numbers but does not have a very large
!     dynamic range, i.e only generates 714025 different numbers
!     for other use consult the above
!     Set idum negative for initialization

      real function ran2(idum)
      implicit none

      integer idum,ir(97),m,ia,ic,iff,iy,j
      real rm
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      save iff,ir,iy
      data iff /0/

      if (idum.lt.0.or.iff.eq.0) then

!     Initialize
!
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum
      end if
!
!     Generate random number
!
      j=1+(97*iy)/m
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum

      return
      end function ran2
    
