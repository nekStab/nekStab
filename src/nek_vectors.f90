!-----------------------------------------------------------------------
! nek_vectors.f90 — Multi-field vector operations for Nek5000 arrays
!
! Purpose:
!   Extends Nek5000's op-style vector operations (opcopy, opadd2, etc.)
!   to the full (velocity, pressure, temperature, passive scalars)
!   field tuple. Each nop* routine operates on all active fields
!   respecting the if3D, ifpo, ifto, ifpsco flags.
!
! Public interface:
!   noprzero   — zero all fields
!   nopcmult   — scalar multiply all fields
!   axpby      — x = alpha*x + beta*y (single array)
!   nopaxpby   — x = alpha*x + beta*y (all fields)
!   nopcopy    — copy all fields
!   nopsub2    — a = a - b (all fields)
!   nopsub3    — c = a - b (all fields)
!   nopadd2    — a = a + b (all fields)
!   nopadd2s2  — a = a + c*b (all fields)
!   opadd3     — a = b + c (velocity only)
!   opaddcol3  — a = a + b*c (velocity only)
!
! Dependencies:
!   nekstab_nek_bridge
!-----------------------------------------------------------------------

module nekstab_vectors
   use nekstab_nek_bridge, only: lx1, ly1, lz1, lx2, ly2, lz2, lelt, nelv, &
                                 nelt, nelfld, ndim, if3D, ifpo, ifto, ldimt, &
                                 npscal, ifpsco, fcx, fcy, fcz, fct
   implicit none
   private
   public :: noprzero, nopcmult, axpby, nopaxpby, &
             nopcopy, nopsub2, nopsub3, nopadd2, &
             nopadd2s2, opadd3, opaddcol3, zero_forcing
contains

!-----------------------------------------------------------------------
! noprzero — Zero all active fields (velocity, pressure, temperature)
!-----------------------------------------------------------------------
   subroutine noprzero(a1, a2, a3, a4, a5)
      integer n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)

      n = lx1*ly1*lz1*nelv
      call rzero(a1, n)
      call rzero(a2, n)
      if (if3D) call rzero(a3, n)
      if (ifpo) call rzero(a4, lx2*ly2*lz2*nelv)
      if (ifto) call rzero(a5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call rzero(a5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine noprzero

!-----------------------------------------------------------------------
! nopcmult — Scalar multiply all active fields: a = c * a
!-----------------------------------------------------------------------
   subroutine nopcmult(a1, a2, a3, a4, a5, c)
      integer n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: c

      n = lx1*ly1*lz1*nelv
      call cmult(a1, c, n)
      call cmult(a2, c, n)
      if (if3D) call cmult(a3, c, n)
      if (ifpo) call cmult(a4, c, lx2*ly2*lz2*nelv)
      if (ifto) call cmult(a5(1, 1), c, lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call cmult(a5(1, k + 1), c, lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopcmult

!-----------------------------------------------------------------------
! axpby — x = alpha*x + beta*y for a single array of length n
!-----------------------------------------------------------------------
   subroutine axpby(x, alpha, y, beta, n)
      real, dimension(*), intent(inout) :: x
      real, dimension(*), intent(in) :: y
      real, intent(in) :: alpha, beta
      integer, intent(in) :: n
      integer :: i
      do i = 1, n
         x(i) = x(i)*alpha + y(i)*beta
      end do

   end subroutine axpby

!-----------------------------------------------------------------------
! nopaxpby — a = alpha*a + beta*b for all active fields
!-----------------------------------------------------------------------
   subroutine nopaxpby(a1, a2, a3, a4, a5, alpha, b1, b2, b3, b4, b5, beta)
      integer n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b4(lx2*ly2*lz2*nelv)
      real, intent(in) :: b5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: alpha, beta

      n = lx1*ly1*lz1*nelv
      call axpby(a1, alpha, b1, beta, n)
      call axpby(a2, alpha, b2, beta, n)
      if (if3D) call axpby(a3, alpha, b3, beta, n)
      if (ifpo) call axpby(a4, alpha, b4, beta, lx2*ly2*lz2*nelv)
      if (ifto) call axpby(a5(1, 1), alpha, b5(1, 1), beta, lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call axpby(a5(1, k + 1), alpha, b5(1, k + 1), beta, lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopaxpby

!-----------------------------------------------------------------------
! nopcopy — Copy all active fields: a = b
!-----------------------------------------------------------------------
   subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      integer :: n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b4(lx2*ly2*lz2*nelv)
      real, intent(in) :: b5(lx1*ly1*lz1*lelt, ldimt)

      n = lx1*ly1*lz1*nelv
      call copy(a1, b1, n)
      call copy(a2, b2, n)
      if (if3D) call copy(a3, b3, n)
      if (ifpo) call copy(a4, b4, lx2*ly2*lz2*nelv)
      if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopcopy

!-----------------------------------------------------------------------
! nopsub2 — Subtract all active fields in-place: a = a - b
!-----------------------------------------------------------------------
   subroutine nopsub2(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      integer n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b4(lx2*ly2*lz2*nelv)
      real, intent(in) :: b5(lx1*ly1*lz1*lelt, ldimt)

      n = lx1*ly1*lz1*nelv
      call sub2(a1, b1, n)
      call sub2(a2, b2, n)
      if (if3D) call sub2(a3, b3, n)
      if (ifpo) call sub2(a4, b4, lx2*ly2*lz2*nelv)
      if (ifto) call sub2(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call sub2(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopsub2

!-----------------------------------------------------------------------
! nopsub3 — Subtract fields into output: c = a - b
!-----------------------------------------------------------------------
   subroutine nopsub3(c1, c2, c3, c4, c5, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      integer n, k
      real, intent(inout) :: c1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: c2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: c3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: c4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: c5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: a1(lx1*ly1*lz1*nelv)
      real, intent(in) :: a2(lx1*ly1*lz1*nelv)
      real, intent(in) :: a3(lx1*ly1*lz1*nelv)
      real, intent(in) :: a4(lx2*ly2*lz2*nelv)
      real, intent(in) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b4(lx2*ly2*lz2*nelv)
      real, intent(in) :: b5(lx1*ly1*lz1*lelt, ldimt)

      n = lx1*ly1*lz1*nelv
      call sub3(c1, a1, b1, n)
      call sub3(c2, a2, b2, n)
      if (if3D) call sub3(c3, a3, b3, n)
      if (ifpo) call sub3(c4, a4, b4, lx2*ly2*lz2*nelv)
      if (ifto) call sub3(c5(1, 1), a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call sub3(c5(1, k + 1), a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopsub3

!-----------------------------------------------------------------------
! nopadd2 — Add all active fields in-place: a = a + b
!-----------------------------------------------------------------------
   subroutine nopadd2(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      integer n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b4(lx2*ly2*lz2*nelv)
      real, intent(in) :: b5(lx1*ly1*lz1*lelt, ldimt)

      n = lx1*ly1*lz1*nelv
      call add2(a1, b1, n)
      call add2(a2, b2, n)
      if (if3D) call add2(a3, b3, n)
      if (ifpo) call add2(a4, b4, lx2*ly2*lz2*nelv)
      if (ifto) call add2(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call add2(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopadd2

!-----------------------------------------------------------------------
! nopadd2s2 — a = a + c*b for all active fields (BLAS-style AXPY)
!-----------------------------------------------------------------------
   subroutine nopadd2s2(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, c)
      integer n, k
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a4(lx2*ly2*lz2*nelv)
      real, intent(inout) :: a5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b4(lx2*ly2*lz2*nelv)
      real, intent(in) :: b5(lx1*ly1*lz1*lelt, ldimt)
      real, intent(in) :: c

      n = lx1*ly1*lz1*nelv
      call add2s2(a1, b1, c, n)
      call add2s2(a2, b2, c, n)
      if (if3D) call add2s2(a3, b3, c, n)
      if (ifpo) call add2s2(a4, b4, c, lx2*ly2*lz2*nelv)
      if (ifto) call add2s2(a5(1, 1), b5(1, 1), c, lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call add2s2(a5(1, k + 1), b5(1, k + 1), c, lx1*ly1*lz1*nelfld(k + 2))
         end do
      end if

   end subroutine nopadd2s2

!-----------------------------------------------------------------------
! opadd3 — a = b + c for velocity fields
!-----------------------------------------------------------------------
   subroutine opadd3(a1, a2, a3, b1, b2, b3, c1, c2, c3)
      integer n
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: c1(lx1*ly1*lz1*nelv)
      real, intent(in) :: c2(lx1*ly1*lz1*nelv)
      real, intent(in) :: c3(lx1*ly1*lz1*nelv)

      n = lx1*ly1*lz1*nelv
      call add3(a1, b1, c1, n)
      call add3(a2, b2, c2, n)
      if (ndim == 3) call add3(a3, b3, c3, n)

   end subroutine opadd3

!-----------------------------------------------------------------------
! opaddcol3 — a = a + b*c for velocity fields (column-wise multiply)
!-----------------------------------------------------------------------
   subroutine opaddcol3(a1, a2, a3, b1, b2, b3, c1, c2, c3)
      integer n
      real, intent(inout) :: a1(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a2(lx1*ly1*lz1*nelv)
      real, intent(inout) :: a3(lx1*ly1*lz1*nelv)
      real, intent(in) :: b1(lx1*ly1*lz1*nelv)
      real, intent(in) :: b2(lx1*ly1*lz1*nelv)
      real, intent(in) :: b3(lx1*ly1*lz1*nelv)
      real, intent(in) :: c1(lx1*ly1*lz1*nelv)
      real, intent(in) :: c2(lx1*ly1*lz1*nelv)
      real, intent(in) :: c3(lx1*ly1*lz1*nelv)

      n = lx1*ly1*lz1*nelv
      call addcol3(a1, b1, c1, n)
      call addcol3(a2, b2, c2, n)
      if (ndim == 3) call addcol3(a3, b3, c3, n)

   end subroutine opaddcol3

!-----------------------------------------------------------------------
! zero_forcing — Zero the forcing arrays (fcx, fcy, fcz, fct)
!
!  Called at the top of each mode dispatch to prevent forcing from
!  accumulating across timesteps.  SFD/TDF/BoostConv set fcx/fcy/fcz
!  each call; without zeroing, old values would be re-applied.
!-----------------------------------------------------------------------
   subroutine zero_forcing
      implicit none
      call rzero(fcx, lx1*ly1*lz1*nelv)
      call rzero(fcy, lx1*ly1*lz1*nelv)
      call rzero(fcz, lx1*ly1*lz1*nelv)
      call rzero(fct, lx1*ly1*lz1*nelt*ldimt)

   end subroutine zero_forcing

end module nekstab_vectors
