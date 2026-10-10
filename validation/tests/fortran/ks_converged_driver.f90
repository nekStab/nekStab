! Reads one case from standard input and prints the converged positions of ks_converged_indices.
! Line 1: n tol.  Then n lines: re(lambda) im(lambda) residual.
! Output line 1: nconv.  Line 2: the positions, 1-based.
program ks_converged_driver
   use nekstab_ks_select, only: ks_converged_indices
   implicit none
   integer :: n, nconv, i
   real :: tol, re, im
   complex(kind(0.0d0)), allocatable :: vals(:)
   real, allocatable :: res(:)
   integer, allocatable :: idx(:)

   read (*, *) n, tol
   allocate (vals(n), res(n), idx(n))
   do i = 1, n
      read (*, *) re, im, res(i)
      vals(i) = cmplx(re, im, kind(0.0d0))
   end do
   call ks_converged_indices(idx, nconv, vals, res, n, tol)
   write (6, '(i6)') nconv
   write (6, '(*(i6))') (idx(i), i=1, nconv)
end program ks_converged_driver
