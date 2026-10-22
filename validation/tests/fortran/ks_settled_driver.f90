! Feeds the Ritz values of successive Arnoldi steps to ks_leading_settled and prints its count per step.
! Line 1: nev nstill tol.  Then one block per step: a line with n, then n lines: re(lambda) im(lambda) residual.
! Output, one line per step: the settled count, then the count of the residual test alone.
program ks_settled_driver
   use nekstab_ks_select, only: ks_leading_settled
   implicit none
   integer :: nev, nstill, n, i, ios, nprev, nset, nconv
   real :: tol, re, im
   complex(kind(0.0d0)), allocatable :: vals(:), lead_prev(:)
   real, allocatable :: res(:)
   integer, allocatable :: still(:)

   read (*, *) nev, nstill, tol
   allocate (lead_prev(nev), still(nev))
   lead_prev = (0.0d0, 0.0d0)
   still = 0
   nprev = 0
   do
      read (*, *, iostat=ios) n
      if (ios /= 0) exit
      allocate (vals(n), res(n))
      do i = 1, n
         read (*, *) re, im, res(i)
         vals(i) = cmplx(re, im, kind(0.0d0))
      end do
      nset = ks_leading_settled(vals, res, n, nev, tol, nstill, lead_prev, nprev, still, nconv)
      write (6, '(2i6)') nset, nconv
      deallocate (vals, res)
   end do
end program ks_settled_driver
