! Prints bound_solver_tol for the cases of test_ks_select.py.
program ks_bound_tol_driver
   use nekstab_ks_select, only: bound_solver_tol
   implicit none
   print *, bound_solver_tol([1.0d-9, 1.0d-9, 1.0d-9, 1.0d0, 0.0d0, -1.0d-6], &
                             [1.0d-7, 1.0d-8, 1.0d-7, 1.0d-7, 1.0d-7, 1.0d-7], &
                             [6.0d4, 6.0d4, 1.2d0, 1.0d-2, 6.0d4, 6.0d4])
end program ks_bound_tol_driver
