      !-----------------------------------------------------------------------
      ! newton_krylov: Main solver routine implementing Newton's method with GMRES for finding steady-state solutions or periodic orbits.
      !
      ! The algorithm follows these steps:
      ! 1. Initialize solver parameters and tolerances
      ! 2. For each Newton iteration:
      !    a. Compute nonlinear residual f(q)
      !    b. Check convergence against target tolerance (dtol)
      !    c. If using dynamic tolerances (ifdyntol), adjust GMRES tolerance and Nek5000 solver tolerances
      !    d. Solve linear system using GMRES to get Newton step
      !    e. Update solution: q = q - dq
      ! 3. Output final solution if converged
      !
      ! Key variables:
      ! - q: Current solution estimate
      ! - f: Nonlinear residual f(q)
      ! - dq: Newton correction step
      ! - dtol: Target tolerance for Newton convergence
      ! - tol: Current GMRES solver tolerance (may be relaxed if ifdyntol=true)
      !-----------------------------------------------------------------------
      subroutine newton_krylov()
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real :: spec_tole ! Function declaration
      
      !     ----- Krylov vectors
         type(krylov_vector) :: f ! Right-hand side vector for Newton solver
         type(krylov_vector) :: q ! Current estimate of the solution
         type(krylov_vector) :: dq ! Newton correction obtained from GMRES
      
      !     ----- Iteration parameters
         integer :: i, j, maxiter_newton, maxiter_gmres, calls
         real :: residual, tol, gmres_target, tottime = 0.0d0
         real :: prev_residual = 0.0d0 ! Track previous iteration's residual
         real :: newton_start_time, newton_iter_time ! Timing for Newton iterations
         integer :: total_gmres_calls ! Track GMRES calls per Newton iteration
         integer :: total_calls = 0, nonlin_calls = 0, lin_calls = 0 ! Call counters
         integer :: prev_total_calls = 0 ! Track calls from previous iteration
         integer, save :: k_out ! Store k from GMRES
         integer, save :: k_sum = 0 ! Accumulator for total k values across Newton iterations
         real, external :: dnekclock
      
      !     ----- Call Counting Logic -----
      ! Each nonlinear solve costs nsteps calls (in nonlinear_forward_map)
      ! Each Arnoldi iteration costs nsteps calls (in ts_gmres)
      ! Total linear calls = sum of all GMRES/Arnoldi calls + 1 per GMRES restart
      ! The +1 accounts for the matvec operation in initialize_gmres_vector (r = b - Ax)
      ! Total nonlinear calls = sum of all nonlinear solver calls
      ! Total calls = nonlin_calls + lin_calls
      ! Calls per iteration = total_calls - prev_total_calls
      ! Each GMRES iteration costs 1 call to initialize_gmres_vector
      ! Each GMRES restart costs 1 call to initialize_gmres_vector
      ! Total calls per GMRES iteration = nsteps + 1
      ! Total calls per GMRES restart = nsteps + 1
      ! Total calls per Newton iteration = nsteps + (nsteps + 1) * (GMRES iterations + GMRES restarts)
      
         real, save :: dtol = 0.0d0 ! Target residual for Newton convergence
         if (dtol == 0.0d0) then
            ! Initialize counters at first call
            total_calls = 0
            nonlin_calls = 0
            lin_calls = 0
            
            dtol = max(param(21), param(22))
            if (nid == 0) write (6, '(A,1PE15.6)') 'dtol saved as:', dtol
            tol = dtol ! Initialize time stepper tolerance
      !     no need to call set_nek5000_tolerances as they are already set
         end if
      
      !     ----- Set iteration limits
         maxiter_newton = 20
         maxiter_gmres = 20
      
      !     ----- Open output files
         if (istep == 0 .and. nid == 0) then
            write (6, *) 'Opening output files for residuals'
            open (unit=887, file='residu_newton.dat', status='replace'); close (887)
            open (unit=888, file='residu_gmres.dat', status='replace'); close (888)
            open (unit=889, file='residu_arnoldi.dat', status='replace'); close (889)
         end if
      
      !     ----- Initialize arrays
         if (nid == 0) write (6, *) 'Initializing Krylov vectors'
         call k_zero(f)
         call k_zero(dq)
      
      !     ----- Copy initial condition
         if (nid == 0) write (6, *) 'Copying initial condition'
         call nopcopy(q%vx, q%vy, q%vz, q%pr, q%t, vx, vy, vz, pr, t)
      
      !     ----- Newton iteration
         newton: do i = 1, maxiter_newton
            if (nid == 0) write (6, *) '------------------------------------------------'
            newton_start_time = dnekclock()
            total_gmres_calls = 0
      
            if (nid == 0) write (6, "('NEWTON   - Starting iteration ',I3,'/',I3,
     $   ' residual:', 1PE15.6, ' solver:', 1PE15.6, ' (target:', 1PE15.6, ')') ")
     $   i-1, maxiter_newton, prev_residual, max(param(21), param(22)), dtol
            if (i == 1) then ! Set time
               q%time = param(10) ! First guess from endTime in .par
               if (nid == 0) write (6, "('  Time    - Initial guess from param(10):', 1PE12.5) ") q%time
            else
               param(10) = q%time ! Other guesses included in nwt field
               if (nid == 0) write (6, "('  Time    - Updated from previous iteration:', 1PE12.5) ") param(10)
            end if
      
            call prepare_linearized_solver ! Compute nsteps and Nek parameters
      
      !     Allocate nonlinear solution variable for natural or forced UPO
            if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
               if (nid == 0) then
                  write (6, "('  Allocating orbit for GMRES')")
                  write (6, "('    Number of steps:',I6)") nsteps
               end if
               allocate (uor(lv, nsteps), vor(lv, nsteps))
               if (if3d) then
                  allocate (wor(lv, nsteps))
               else
                  allocate (wor(1, 1))
               end if
               if (ifto .or. ldimt > 1) allocate (tor(lt, nsteps, ldimt))
            end if
      
            call nonlinear_forward_map(f, q) ! rhs of Newton iteration f(q)
            nonlin_calls = nonlin_calls + nsteps ! nsteps of the nonlinear solver
            total_calls = total_calls + nonlin_calls
            tottime = tottime + time ! use time instead of nsteps*dt
      
      !     Check residual || f(q) ||
            call k_norm(residual, f) ! L2 norm: square root of dot product with weighted norms by bm1s
            residual = residual**2 ! squared L2 norm for consistent convergence check with GMRES
            if (nid == 0) write (6, *) '  Computed residual:', residual
      
      !     --> Outpost residual fields (optional)
            time = q%time ! adjust
            if (uparam(1) == 2.0) time = real(i - 1) ! to ease visu in paraview
            if (uparam(1) > 2.0) time = real(i - 1)*q%time ! to ease visu in paraview
            call outpost2(f%vx, f%vy, f%vz, f%pr, f%t, nof, 'res')
            time = q%time ! restore
      
      !     Output iteration information
            if (nid == 0) then
               newton_iter_time = dnekclock() - newton_start_time
               write (6, "('NEWTON   - Finished iteration ',I3,'/',I3,
     $   ' residual:', 1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter_newton, residual, dtol
               if (i > 1) then ! Only show rate after first iteration
                  write (6, "('          Change: ',A,1PE15.6, ' Rate:', 1PE15.6) ")
     $   merge('↑', '↓', residual > prev_residual), abs(residual - prev_residual), 
     $   residual/prev_residual
               end if
               write (6, "('          Time: ',1PE15.6,'s  GMRES calls: ', I4) ") newton_iter_time, total_gmres_calls
               open (887, file='residu_newton.dat', action='write', position='append')
               write (887, "(4I9,4(1PE15.6))") i-1, total_calls, total_calls - prev_total_calls, k_sum, tottime, 
     $          max(param(21), param(22)), residual, dtol
               close (887)
               k_sum = k_sum + k_out  ! Update k_sum with k from GMRES (moved after writing)
               prev_total_calls = total_calls ! Save for next iteration
               prev_residual = residual ! Save for next iteration
               write (6, *) '------------------------------------------------'
            end if
      
      !     Check for convergence
            if (residual < dtol) then
               if (nid == 0) write (6, *) '  Converged. Exiting Newton loop...'
               exit newton
            end if
      
      !     If ifdyntol=true, two steps occur:
      !     1. Compute relaxed tolerance: tol = spec_tole(residual, dtol, 0.1d0)
      !        - tol will be between dtol and min_tol
      !        - typically around 0.1*residual to avoid over-solving
      !     2. Update Nek5000 solver tolerances via set_nek5000_tolerances(tol)
            if (ifdyntol) then
               tol = spec_tole(residual, dtol, 0.1d0) ! compute the relaxed tolerance
               call set_nek5000_tolerances(tol) ! set the tolerance to the time-stepper
            end if
      
      !     Solve the linear system
            if (nid == 0) write (6, *) '  Solving linear system with GMRES'
            call ts_gmres(f, dq, maxiter_gmres, k_dim, tol, calls, k_out, i, dtol)
            lin_calls = lin_calls + calls + 1  ! Add GMRES calls + 1 for matvec in initialize_gmres_vector
            total_calls = total_calls + lin_calls
            tottime = tottime + calls*dt
            total_gmres_calls = total_gmres_calls + calls
      
      !     Update Newton solution
            if (nid == 0) write (6, *) '  Updating Newton solution'
            call k_sub2(q, dq)
      
      !     Outpost current estimate of the solution ! NOT SURE ABOUT THIS POSITION
            if (nid == 0) write (6, *) '  Outposting current solution estimate'
            time = q%time
            if (uparam(1) == 2.0) time = real(i - 1) ! Ease visualization in paraview: file 1 is t = 0
            if (uparam(1) == 2.2) time = real(i - 1)*q%time
            call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, 'nwt')
            time = q%time ! Restore
      
      !     Deallocate nonlinear solution variable
            if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
               if (nid == 0) write (6, *) '  Deallocating orbit storage'
               deallocate (uor, vor, wor)
               if (ifto .or. ldimt > 1) deallocate (tor)
            end if
      
         end do newton
      
      !     Output final results
         if (nid == 0) then
            if (i == maxiter_newton) then
               write (6, *) 'Reached maxiter_newto. STOPPING! (verify convergence)'
            else
               if (uparam(1) == 2) then
                  write (6, *) 'NEWTON finished successfully after', i, 'iterations.'
               elseif (uparam(1) == 2.1) then
                  write (6, *) 'NEWTON UPO finished successfully', i, 'iterations.'
                  write (6, *) ' period found:', time, 1.0d0/time
               elseif (uparam(1) == 2.2) then
                  write (6, *) 'NEWTON for forced UPO finished successfully', i, 'iterations.'
                  write (6, *) ' period found:', time, 1.0d0/time
               end if
               write (6, *) 'Calls to the linearized solver: ', total_calls
               write (6, *) 'Total nondimensional time:', tottime
               if (ifdyntol) write (6, *) 'ifdyntol active!'
            end if
         end if
      
      !     Output solution if converged
         if (residual < dtol) then
            if (nid == 0) write (6, *) 'Outputting converged solution'
            param(63) = 1.0d0 ! Enforce 64-bit output
            call bcast(param(63), wdsize)
            call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, "BF_")
            param(63) = 0.0d0 ! Enforce 32-bit output
            call bcast(param(63), wdsize)
            call outpost_vort(vx, vy, vz, 'BFV')
         end if
      
         if (nid == 0) write (6, *) 'newton_krylov subroutine completed'
         return
      end subroutine
      
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
      ! ts_gmres: Time-stepper GMRES solver for the Newton correction equation
      !
      ! This implements a restarted GMRES method to solve the linear system:
      !    J(q)dq = -f(q)
      ! where J is the Jacobian matrix (approximated via finite differences)
      !
      ! Algorithm steps:
      ! 1. Initialize Krylov subspace and Hessenberg matrix
      ! 2. For each GMRES iteration:
      !    a. Perform Arnoldi process to build orthonormal basis
      !    b. Solve least squares problem for coefficients
      !    c. Check convergence of residual against tol
      !    d. If not converged, continue building Krylov subspace
      ! 3. Update solution using computed coefficients
      !
      ! Key parameters:
      ! - ksize: Size of Krylov subspace before restart
      ! - maxiter: Maximum number of GMRES iterations
      ! - tol: Convergence tolerance for GMRES
      ! - dtol: Target Newton tolerance (for output only)
      !-----------------------------------------------------------------------
      subroutine ts_gmres(rhs, sol, maxiter, ksize, tol, calls, k_out, newton_iter, dtol)
      
      !     Implementation of simple time-stepper GMRES to be part of the Newton-Krylov solver
      !     for fixed point computation. The rank of the Krylov subspace is set as the user parameter k_dim.
      !
      
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         type(krylov_vector), intent(in) :: rhs
         type(krylov_vector), intent(inout) :: sol
         integer, intent(in) :: maxiter, ksize
         real, intent(in) :: tol, dtol
         integer, intent(out) :: calls, k_out
         integer, intent(in) :: newton_iter
      
         type(krylov_vector) :: dq
         type(krylov_vector), dimension(:), allocatable :: Q
         real, dimension(:, :), allocatable :: H
         real, dimension(:), allocatable :: yvec, evec
         real :: beta, beta2
         integer :: i, j, k ! Added m for orthogonalization loop
         real :: gmres_start_time, gmres_iter_time
         real :: prev_beta2 = 0.0d0
         real :: arnoldi_start_time, arnoldi_iter_time
         real :: ortho_metric ! Measure of orthogonalization quality
         real :: dot_val ! For weighted inner product calculation
         real, external :: dnekclock
         integer, save :: k_sum = 0  ! Accumulator for total k values in GMRES
      
         calls = 0
      
      !     ----- Allocate arrays for GMRES -----
         allocate (Q(ksize + 1), H(ksize + 1, ksize), yvec(ksize), evec(ksize + 1))
      
      !     Initialize Krylov basis vectors and Hessenberg matrix
         H = 0.0d0
         yvec = 0.0d0
         evec = 0.0d0
         call k_zero(sol)
         do i = 1, ksize + 1
            call k_zero(Q(i))
         end do
      
      !     Set first Krylov vector as normalized residual
         call k_copy(Q(1), rhs)
         call k_normalize(Q(1), beta)
         beta2 = beta**2 ! L2 norm squared for convergence check
      
         gmres: do i = 1, maxiter
            if (nid == 0) write (6, *) '  ------------------------------------'
            gmres_start_time = dnekclock()
      
            if (nid == 0) write (6, "('  GMRES   - Starting iteration ',I3,'/',I3,
     $   ' residual:', 1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter, beta2, tol
      
      !     Reset arrays for new restart
            H = 0.0d+00
            yvec = 0.0d+00
            evec = 0.0d+00
            evec(1) = beta
            do j = 2, ksize + 1
               call k_zero(Q(j))
            end do
      
            arnoldi: do k = 1, k_dim
      
               if (nid == 0) write (6, "('    ARNOLDI [GMRES ',I3,'/',I3,'] Starting iteration ', I3, '/', I3) ")
     $   i, maxiter, k, ksize
               arnoldi_start_time = dnekclock()
      
      !     Save residual for convergence rate calculation
               prev_beta2 = beta2
      
      !     Build k-th column of Hessenberg matrix and orthogonalize
               call arnoldi_factorization(Q, H, k, k, ksize)
               calls = calls + nsteps  ! Each Arnoldi step costs nsteps calls
      
      !     Solve least squares problem for coefficients
               call lstsq(H(1:k + 1, 1:k), evec(1:k + 1), yvec(1:k), k + 1, k)
      
      !     Compute residual norm of least squares solution
               beta = norm2(evec(1:k + 1) - matmul(H(1:k + 1, 1:k), yvec(1:k)))
               beta2 = beta**2
      
      !     Compute orthogonalization quality metric using proper weighted inner products
               ortho_metric = 0.0d0
               if (k > 1) then
                  do j = 1, k - 1
                     call k_dot(dot_val, Q(k), Q(j))  ! Use weighted inner product
                     ortho_metric = ortho_metric + abs(dot_val)
                  end do
                  ortho_metric = ortho_metric/(k - 1)  ! Average orthogonalization error
               end if
      
               if (nid == 0) then
                  arnoldi_iter_time = dnekclock() - arnoldi_start_time
                  write (6, "('    ARNOLDI [GMRES ',I3,'/',I3,'] residual:',
     $   1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter, beta2, tol
                  if (k > 1) then ! Only show rate and ortho after first iteration
                     write (6, "('              Rate:',1PE15.6,' Time:',1PE15.6,
     $   's  Ortho:', 1PE15.6) ") beta2/prev_beta2, arnoldi_iter_time, ortho_metric
                  else
                     write (6, "('              Time:',1PE15.6,'s') ") arnoldi_iter_time
                  end if

                  open (889, file='residu_arnoldi.dat', action='write', position='append')
                  write (889, "(I9,5(1PE15.6))") k, arnoldi_iter_time, ortho_metric, tol, beta2, dtol
                  close (889)
                  write (6, *) '    ...........................'
               end if
      
      !     Store current residual for next iteration's rate calculation
               prev_beta2 = beta2
      
               if (beta2 < tol) then ! count of calls to linearized solver
                  exit arnoldi
               end if
      
      ! --> Relaxed exit condition if finite-difference approximation of the operator is considered.
               if ((iffindiff) .and. (beta2 < 1e-8)) then ! count of calls to linearized solver
                  exit arnoldi
               end if
            end do arnoldi
      
      !     Update solution using Krylov basis and coefficients
            call k_matmul(dq, Q(1:k), yvec(1:k), k)
            call k_add2(sol, dq)
      
      !     Verify residual and prepare for next restart if needed
            call k_copy(Q(1), sol)
            call initialize_gmres_vector(beta, Q(1), rhs)
            beta2 = beta**2
      
            if (nid == 0) then
               gmres_iter_time = dnekclock() - gmres_start_time
               write (6, "('  GMRES   - Finished iteration:',I3,'/',I3, ' residual:', 1PE15.6) ") i, maxiter, beta2
               if (i > 1) then ! Only show rate after first iteration
                  write (6, "('           Rate:',1PE15.6,' Time:',1PE15.6,
     $   's  Arnoldi steps:', I4) ") beta2/prev_beta2, gmres_iter_time, k
               else
                  write (6, "('           Time:',1PE15.6,'s  Arnoldi steps:', I4) ") 
     $   gmres_iter_time, k
               end if
               k_sum = k_sum + k  ! Update accumulator after writing
               open (888, file='residu_gmres.dat', action='write', position='append')
               write (888, "(4I9,3(1PE15.6))") newton_iter, i, k, k_sum, tol, beta2, dtol
               close (888)
               
               prev_beta2 = beta2
               write (6, *) '  ------------------------------------'
            end if
      
            if (beta2 < tol) then
               exit gmres
            end if
            if ((iffindiff) .and. (beta2 < 1e-6)) then
               exit gmres
            end if
         end do gmres
      
      !     ----- Deallocate arrays -----
         deallocate (Q, H, yvec, evec)
      
         k_out = k  ! Save final k value
      end subroutine ts_gmres
      
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
      ! initialize_gmres_vector: Prepares initial vector for GMRES iterations
      !
      ! This routine:
      ! 1. Computes initial residual r = b - Ax
      ! 2. Normalizes the residual to create first Krylov vector
      ! 3. Returns the norm (beta) for use in GMRES
      !
      ! The normalized residual becomes the first basis vector of the
      ! Krylov subspace used in the Arnoldi process.
      !-----------------------------------------------------------------------
      subroutine initialize_gmres_vector(beta, q, rhs)
         use krylov_subspace
      
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Right-hand side vector of A x = b (input) -----
         type(krylov_vector), intent(in) :: rhs
      
      !     ----- Seed for the GMRES Krylov subspace (output) -----
         type(krylov_vector), intent(inout) :: q
         type(krylov_vector) :: f
         real, intent(out) :: beta
      
      !     --> Compute initial residual: r = b - Ax
         call matvec(f, q)
         call k_sub2(f, rhs)
         call k_cmult(f, -1.0d0)
      
      !     --> Normalize the starting vector.
         call k_normalize(f, beta)
         call k_copy(q, f)
      
         return
      end subroutine initialize_gmres_vector
      
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
      ! nonlinear_forward_map: Computes nonlinear residual f(q) = F(q) - q
      !
      ! This routine:
      ! 1. Sets initial condition from current solution estimate q
      ! 2. Advances solution forward in time for nsteps
      ! 3. Optionally stores trajectory for UPO computations
      ! 4. Computes residual as difference between final and initial states
      ! 5. Updates base flow for linearized calculations
      !
      ! For periodic orbits (UPOs):
      ! - Stores full trajectory in uor, vor, wor arrays
      ! - Handles both natural (2.1) and forced (2.2) UPO cases
      !-----------------------------------------------------------------------
      subroutine nonlinear_forward_map(f, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Right-hand side of the Newton.
         type(krylov_vector), intent(out) :: f
      !     ----- Initial condition for the forward simulation.
         type(krylov_vector), intent(in) :: q
      
         integer :: m
         real :: gmres_target
         nt = nx1*ny1*nz1*nelt
      
      !     --> Copy the initial condition to Nek.
         call k_copy(ic_nwt, q) ! for Newton UPO newton_linearized_map
         call nopcopy(vx, vy, vz, pr, t, q%vx, q%vy, q%vz, q%pr, q%t)
      
      !     --> Turn-off the linearized solver.
         ifpert = .false.
         call bcast(ifpert, lsize)
      
      !     --> Run the simulation forward.
         time = 0.0d0
         do istep = 1, nsteps
            call nekStab_usrchk()
            call nek_advance()
            if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
               if (nid == 0) then
                  write (6, *) 'Storing nonlinear solution for GMRES:', istep, '/', nsteps
               end if
               call opcopy(uor(:, istep), vor(:, istep), wor(:, istep), vx, vy, vz)
               if (ifto) call copy(tor(:, istep, 1), t(:, :, :, :, 1), nt)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) call copy(tor(:, istep, m), t(:, :, :, :, m), nt)
                  end do
               end if
            end if
         end do
      
      !     --> Compute the right hand side of the time-stepper Newton.
         call nopcopy(f%vx, f%vy, f%vz, f%pr, f%t, vx, vy, vz, pr, t)
         call k_copy(fc_nwt, f) ! for Newton UPO newton_linearized_map
         call k_sub2(f, q)
         f%time = 0.0d0
      
      !     --> Pass current guess as base flow for the linearized calculation.
         call nopcopy(ubase, vbase, wbase, pbase, tbase, q%vx, q%vy, q%vz, q%pr, q%t)
      
         return
      end subroutine nonlinear_forward_map
      
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
      ! set_nek5000_tolerances: Updates all Nek5000 solver tolerances
      !
      ! This routine ensures consistency across all solvers by:
      ! 1. Setting pressure and velocity solver tolerances (param 21,22)
      !
      ! All changes are broadcast to maintain consistency across MPI ranks
      !-----------------------------------------------------------------------
      subroutine set_nek5000_tolerances(solver_tol) ! set solver tolerances
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, intent(in) :: solver_tol ! New tolerance value to be set
         if (nid == 0) write (6, "('  TOLERANCE set from:',1PE15.6,' to:', 1PE15.6) ") param(21), abs(solver_tol)
      
      ! Set both param(21) and param(22) to the absolute value of the new tolerance.
      ! Broadcast these changes to all nodes.
         param(21) = abs(solver_tol)
         call bcast(param(21), wdsize)
         param(22) = abs(solver_tol)
         call bcast(param(22), wdsize)
      
      ! Update TOLPDF and TOLHDF with the new tolerance and broadcast the changes.
         TOLPDF = param(21)
         call bcast(TOLPDF, wdsize)
         TOLHDF = param(22)
         call bcast(TOLHDF, wdsize)
      
      ! Update restol and atol with the new tolerance and broadcast the changes.
         restol(:) = param(22)
         call bcast(restol, (ldimt1 + 1)*wdsize)
         atol(:) = param(22)
         call bcast(atol, (ldimt1 + 1)*wdsize)
      
      end subroutine set_nek5000_tolerances
      
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
      ! spec_tole: Computes relaxed tolerance for GMRES solver
      !
      ! This implements an adaptive tolerance strategy:
      ! - Early iterations: Use relaxed tolerance ≈ relaxation_factor * residual
      !   to avoid over-solving when far from solution
      ! - Later iterations: Tolerance approaches dtol as residual decreases
      ! - Bounded between dtol (min) and min_tol (max) for stability
      !
      ! Parameters:
      ! - residual: Current Newton residual norm
      ! - dtol: Target tolerance for Newton convergence
      ! - relaxation_factor: Controls how much to relax tolerance (typically 0.1)
      !-----------------------------------------------------------------------
      function spec_tole(residual, dtol, relaxation_factor) result(nwtol)
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, intent(in) :: residual ! Current residual norm
         real, intent(in) :: dtol ! Target tolerance
         real, intent(in) :: relaxation_factor ! Relaxation factor for solver tolerances
         real :: nwtol ! Returned new tolerance
         real, parameter :: min_tol = 1.0d-5 ! Minimum allowed tolerance
         
         ! 1.51485E+02 e-5
         ! 1.89187E+02 e-6
      
      ! Compute new time stepper tolerances based on Newton residual
      ! Adjusts how accurately we solve the time stepping problem:
      ! - Early Newton iterations: Relaxed tolerances (≈ relaxation_factor * residual)
      ! - Later iterations: Stricter tolerances approaching dtol
      ! - Bounded by dtol (min) and min_tol (max) for stability
         nwtol = max(min(residual*relaxation_factor, min_tol), dtol)
         if (nid == 0) then
            if (nwtol == min_tol) then
               write (6, "('  [TOLERANCE at max limit:',1PE15.6,']')") min_tol
            else if (nwtol == dtol) then
               write (6, "('  [TOLERANCE at min limit:',1PE15.6,']')") dtol
            end if
         end if
      
      end function spec_tole
