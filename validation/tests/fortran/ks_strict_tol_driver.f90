! Prints strict_solver_tol for the cases of test_ks_select.py.
program ks_strict_tol_driver
   use nekstab_ks_select, only: strict_solver_tol
   implicit none
   print *, strict_solver_tol([1.0d-9, 1.0d-9, 1.0d-9, 0.0d0], [6.0d4, 1.2d0, 0.5d0, 6.0d4])
end program ks_strict_tol_driver
