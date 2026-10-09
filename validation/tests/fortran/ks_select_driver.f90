! Reads one case from standard input and prints the selection of ks_select_eigenvalues.
! Line 1: n nev delta tol.  Then n lines: re(lambda) im(lambda) residual.
! Output line 1: nsel n_circle n_resid.  Line 2: one 0/1 per entry.
program ks_select_driver
   use nekstab_ks_select, only: ks_select_eigenvalues
   implicit none
   integer :: n, nev, nsel, n_circle, n_resid, i
   real :: delta, tol, re, im
   complex(kind(0.0d0)), allocatable :: vals(:)
   real, allocatable :: res(:)
   logical, allocatable :: selected(:)

   read (*, *) n, nev, delta, tol
   allocate (vals(n), res(n), selected(n))
   do i = 1, n
      read (*, *) re, im, res(i)
      vals(i) = cmplx(re, im, kind(0.0d0))
   end do
   call ks_select_eigenvalues(selected, nsel, vals, res, delta, nev, n, tol, n_circle, n_resid)
   write (6, '(3i6)') nsel, n_circle, n_resid
   write (6, '(*(i1))') (merge(1, 0, selected(i)), i=1, n)
end program ks_select_driver
