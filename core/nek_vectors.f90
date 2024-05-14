      
      !-----------------------------------------------------------------------
      
      subroutine noprzero(a1, a2, a3, a4, a5)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n, k
         real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         n = nx1*ny1*nz1*nelv
         call rzero(a1, n)
         call rzero(a2, n)
         if (if3D) call rzero(a3, n)
         if (ifpo) call rzero(a4, nx2*ny2*nz2*nelv)
         if (ifto) call rzero(a5(1, 1), lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call rzero(a5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine noprzero
      !-----------------------------------------------------------------------
      subroutine nopcmult(a1, a2, a3, a4, a5, c)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n, k
         real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: c
         n = nx1*ny1*nz1*nelv
         call cmult(a1, c, n)
         call cmult(a2, c, n)
         if (if3D) call cmult(a3, c, n)
         if (ifpo) call cmult(a4, c, nx2*ny2*nz2*nelv)
         if (ifto) call cmult(a5(1, 1), c, lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call cmult(a5(1, k + 1), c, lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine nopcmult
      !-----------------------------------------------------------------------
      subroutine axpby(x, alpha, y, beta, n)
         real x(1), y(1), alpha, beta
         do i = 1, n
            x(i) = x(i)*alpha + y(i)*beta
         end do
         return
      end subroutine axpby
      subroutine nopaxpby(a1, a2, a3, a4, a5, alpha, b1, b2, b3, b4, b5, beta)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n, k
         real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: alpha, beta
         n = nx1*ny1*nz1*nelv
         call axpby(a1, alpha, b1, beta, n)
         call axpby(a2, alpha, b2, beta, n)
         if (if3D) call axpby(a3, alpha, b3, beta, n)
         if (ifpo) call axpby(a4, alpha, b4, beta, nx2*ny2*nz2*nelv)
         if (ifto) call axpby(a5(1, 1), alpha, b5(1, 1), beta, lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call axpby(a5(1, k + 1), alpha, b5(1, k + 1), beta, lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine nopaxpby
      !-----------------------------------------------------------------------
      subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer :: n, k
         real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
         n = nx1*ny1*nz1*nelv
         call copy(a1, b1, n)
         call copy(a2, b2, n)
         if (if3D) call copy(a3, b3, n)
         if (ifpo) call copy(a4, b4, nx2*ny2*nz2*nelv)
         if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine nopcopy
      !-----------------------------------------------------------------------
      subroutine nopsub2(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n, k
         real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
         n = nx1*ny1*nz1*nelv
         call sub2(a1, b1, n)
         call sub2(a2, b2, n)
         if (if3D) call sub2(a3, b3, n)
         if (ifpo) call sub2(a4, b4, nx2*ny2*nz2*nelv)
         if (ifto) call sub2(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call sub2(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine nopsub2
      !-----------------------------------------------------------------------
      subroutine nopsub3(c1, c2, c3, c4, c5, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n, k
         real, intent(inout) :: c1(1), c2(1), c3(1), c4(1), c5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
         n = nx1*ny1*nz1*nelv
         call sub3(c1, a1, b1, n)
         call sub3(c2, a2, b2, n)
         if (if3D) call sub3(c3, a3, b3, n)
         if (ifpo) call sub3(c4, a4, b4, nx2*ny2*nz2*nelv)
         if (ifto) call sub3(c5(1, 1), a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call sub3(c5(1, k + 1), a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine nopsub3
      !-----------------------------------------------------------------------
      subroutine nopadd2(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n, k
         real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
         real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
         n = nx1*ny1*nz1*nelv
         call add2(a1, b1, n)
         call add2(a2, b2, n)
         if (if3D) call add2(a3, b3, n)
         if (ifpo) call add2(a4, b4, nx2*ny2*nz2*nelv)
         if (ifto) call add2(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
         if (ldimt > 1) then
         do k = 1, npscal
            if (ifpsco(k)) call add2(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
         end do
         end if
         return
      end subroutine nopadd2
      !-----------------------------------------------------------------------
      subroutine opadd3(a1, a2, a3, b1, b2, b3, c1, c2, c3)
         implicit none
         include 'SIZE'
         integer n
         real a1(1), a2(1), a3(1), b1(1), b2(1), b3(1), c1(1), c2(1), c3(1)
         n = nx1*ny1*nz1*nelv
         call add3(a1, b1, c1, n)
         call add3(a2, b2, c2, n)
         if (ndim == 3) call add3(a3, b3, c3, n)
         return
      end subroutine opadd3
      !-----------------------------------------------------------------------
      subroutine opaddcol3(a1, a2, a3, b1, b2, b3, c1, c2, c3)
         implicit none
         include 'SIZE'
         integer n
         real a1(1), a2(1), a3(1), b1(1), b2(1), b3(1), c1(1), c2(1), c3(1)
         n = nx1*ny1*nz1*nelv
         call addcol3(a1, b1, c1, n)
         call addcol3(a2, b2, c2, n)
         if (ndim == 3) call addcol3(a3, b3, c3, n)
         return
      end subroutine opaddcol3
