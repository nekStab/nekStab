!-----------------------------------------------------------------------





      subroutine newton_krylov()

!     Implementation of a simple Newton-Krylov solver for fixed point
!     computation using a time-stepper formulation. The resolution of
!     the linear system for each Newton iteration is sovled by means
!     of GMRES.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      return
      end




!-----------------------------------------------------------------------




      subroutine gmres()

!     Implementation of simple GMRES to be part of the Newton-Krylov solver
!     for fixed point computation.

      return
      end




!-----------------------------------------------------------------------
