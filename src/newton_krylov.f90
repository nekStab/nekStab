      !-----------------------------------------------------------------------
      ! newton_krylov.f90 — Newton-Krylov solver for fixed points and UPOs
      !
      ! Purpose:
      !   Implements Newton's method with GMRES for finding steady-state
      !   solutions or periodic orbits. Includes adaptive tolerance and
      !   residual monitoring.
      !
      ! Public interface:
      !   newton_krylov, ts_gmres, initialize_gmres_vector,
      !   nonlinear_forward_map, set_nek5000_tolerances, spec_tole
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL
      !-----------------------------------------------------------------------

      module nekstab_newton
         use krylov_subspace
         use nekstab_vectors
         use nekstab_matvec
         use nekstab_krylov_decomposition
         use nekstab_lapack
         use nekstab_diagnostics
         implicit none
         private
         public :: newton_krylov, ts_gmres,
     $             initialize_gmres_vector,
     $             nonlinear_forward_map,
     $             set_nek5000_tolerances, spec_tole
      contains

      !-----------------------------------------------------------------------
      ! newton_krylov — Main Newton iteration loop
      !
      ! The algorithm solves F(q) = q using Newton iteration:
      ! 1. Compute residual: f = F(q) - q
      ! 2. Solve linear system: J(q)dq = -f using GMRES
      ! 3. Update solution: q = q - dq
      ! 4. Repeat until ||f|| < dtol
      !
      ! Sign convention:
      ! - Residual f = F(q) - q (computed in nonlinear_forward_map)
      ! - GMRES solves J(q)dq = -f for the Newton step dq
      ! - Update uses q = q - dq to move toward the fixed point
      !
      ! Key variables:
      ! - q: Current solution estimate
      ! - f: Nonlinear residual f(q) = F(q) - q
      ! - dq: Newton correction step (dq)
      ! - dtol: Target tolerance for Newton convergence
      ! - tol: Current GMRES solver tolerance (may be relaxed if ifdyntol=true)
      !-----------------------------------------------------------------------
      subroutine newton_krylov()
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

      !     ----- Krylov vectors
         type(krylov_vector) :: f ! Right-hand side vector for Newton solver
         type(krylov_vector) :: q ! Current estimate of the solution
         type(krylov_vector) :: dq ! Newton correction obtained from GMRES

      !     ----- Iteration parameters
         integer :: i, j, maxiter_newton, maxiter_gmres, calls
         real :: residual, tol, tottime = 0.0d0
         real :: prev_residual = 0.0d0 ! Previous iteration residual (EW + stagnation)
         real :: initial_residual = 0.0d0 ! For divergence guard
         integer :: stagnation_count = 0 ! Stagnation detection counter
         real :: newton_start_time, newton_iter_time ! Timing for Newton iterations
         integer :: total_gmres_calls ! Track GMRES calls per Newton iteration
         integer :: total_calls = 0, nonlin_calls = 0, lin_calls = 0 ! Call counters
         integer :: prev_total_calls = 0 ! Track calls from previous iteration
         integer, save :: k_out ! Store k from GMRES
         integer, save :: k_sum = 0 ! Accumulator for total k values across Newton iterations
         real, external :: dnekclock
         real :: saved_tol21, saved_tol22 ! Save/restore param(21:22)

      !     ----- Call Counting -----
! total_calls  = cumulative time-steps (nonlinear + linear matvecs)
! nonlin_calls = cumulative nonlinear time-steps
! lin_calls    = cumulative linear matvecs (+1 per GMRES for init vector)

         real, save :: dtol = 0.0d0 ! Target residual for Newton convergence
         if (dtol == 0.0d0) then ! Initialize counters at first call
            total_calls = 0
            nonlin_calls = 0
            lin_calls = 0

            dtol = max(param(21), param(22))
            if (nid == 0) write (6, '(A,1PE15.6)') 'dtol saved as:', dtol
            tol = dtol ! Initialize time stepper tolerance: no need to call set_nek5000_tolerances as they are already set
         end if

         maxiter_newton = 30; maxiter_gmres = 30

         if (istep == 0 .and. nid == 0) then
            write (6, *) 'Opening output files for residuals'
            open (unit=887, file='residu_newton.dat', status='replace'); close (887)
            open (unit=888, file='residu_gmres.dat', status='replace'); close (888)
            open (unit=889, file='residu_arnoldi.dat', status='replace'); close (889)
         end if

         if (nid == 0) write (6, *) 'Initializing Krylov vectors'
         call k_zero(f); call k_zero(dq)

         if (nid == 0) write (6, *) 'Copying initial condition'
         call nopcopy(q%vx, q%vy, q%vz, q%pr, q%t, vx, vy, vz, pr, t)

c        Save original tolerances (restored after Newton exits).
c        Without this, EW leaves param(21:22) at whatever the last
c        Newton iteration set — corrupting any subsequent DNS or
c        stability run that reads param(21:22) for its own tolerances.
         saved_tol21 = param(21)
         saved_tol22 = param(22)

         newton: do i = 1, maxiter_newton
            if (nid == 0) write (6, *) '------------------------------------------------'
            newton_start_time = dnekclock()
            total_gmres_calls = 0

            if (nid == 0) write (6, "('NEWTON   - Starting iteration ',I3,'/',I3,
     $   ' residual:', 1PE15.6, ' solver:', 1PE15.6, ' (target:', 1PE15.6,')')")
     $   i, maxiter_newton, prev_residual, max(param(21), param(22)), dtol
            if (i == 1) then ! Set time
               q%time = param(10) ! First guess from endTime in .par
               if (nid == 0) write (6, "('  Time    - Initial guess from param(10):', 1PE15.6) ") q%time
            else
               param(10) = q%time ! Other guesses included in nwt field
               if (nid == 0) write (6, "('  Time    - Updated from previous iteration:', 1PE15.6) ") param(10)
            end if

            call prepare_linearized_solver ! Compute nsteps and Nek parameters

            if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
               if (nid == 0) then ! Allocate nonlinear solution variable for natural or forced UPO
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
            nonlin_calls = nonlin_calls + nsteps ! cumulative nonlinear time-steps
            total_calls = total_calls + nsteps
            tottime = tottime + time ! use time instead of nsteps*dt

c           Cache bvec/btvec for UPO (avoids recomputing per matvec)
            if (isNewtonPO)
     $         call cache_newton_bvec(ic_nwt, fc_nwt)

      !     Check residual || f(q) ||! L2 norm: square root of dot product with weighted norms by bm1s
            call k_norm(residual, f) ! Computes ||f||
            residual = residual**2 ! squared L2 norm for consistent convergence check with GMRES
            if (nid == 0) write (6, *) '  Computed residual:', residual

c           Save initial residual for divergence guard
            if (i == 1) initial_residual = residual

      !     --> Outpost residual fields (optional)
            time = q%time ! adjust
            if (uparam(1) == 2.0) time = real(i - 1) ! to ease visu in paraview
            if (uparam(1) > 2.0) time = real(i - 1)*q%time ! to ease visu in paraview
            call outpost2(f%vx, f%vy, f%vz, f%pr, f%t, nof, 'res')
            time = q%time ! restore


c           Stagnation detection (before prev_residual is overwritten)
            if (i > 1) then
               if (residual > 0.9d0*prev_residual) then
                  stagnation_count = stagnation_count + 1
               else
                  stagnation_count = 0
               end if
            end if

            if (nid == 0) then ! Output iteration information
               newton_iter_time = dnekclock() - newton_start_time
               k_sum = k_sum + k_out  ! Update k_sum with k from GMRES
               write (6, "('NEWTON   - Finished iteration ',I3,'/',I3,
     $   ' residual:', 1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter_newton, residual, dtol
               if (i > 1) then ! Only show rate after first iteration
                  write (6, "('          Change: ',A,1PE15.6, ' Rate:', 1PE15.6) ")
     $   merge('↑', '↓', residual > prev_residual), abs(residual - prev_residual),
     $   residual/prev_residual
               end if
               if (stagnation_count >= 3) then
                  write (6, *) 'WARNING: Newton stagnation',
     $               ' (3 iterations without progress)'
               end if
               write (6, "('          Time: ',1PE15.6,'s  GMRES calls: ', I4) ") newton_iter_time, total_gmres_calls
               open (887, file='residu_newton.dat', action='write', position='append')
               write (887, "(4I9,4(1PE15.6))") i, total_calls, total_calls - prev_total_calls, k_sum, tottime,
     $          max(param(21), param(22)), residual, dtol
               close (887)
               prev_total_calls = total_calls ! Save for next iteration
               write (6, *) '------------------------------------------------'
            end if
            if (residual < dtol) then
               if (nid == 0) write (6, *) '  Converged. Exiting Newton loop...'
               exit newton
            end if

c           Divergence guard: abort if residual grows excessively
            if (i > 1 .and. residual > 1.0d8*initial_residual) then
               if (nid == 0) write (6, *)
     $            'NEWTON: DIVERGENCE — residual exceeds',
     $            ' 1e8 × initial. Aborting.'
               exit newton
            end if

c           EW adaptive tolerance; sqrt converts ||r||^2 -> ||r|| for param(21:22)
            if (ifdyntol) then
               tol = spec_tole(residual, prev_residual, dtol)
               call set_nek5000_tolerances(sqrt(tol))
            end if
            prev_residual = residual ! after stagnation check + EW use the old value

            if (nid == 0) write (6, *) '  Solving linear system with GMRES for rhs = f = F(q) - q'
            call ts_gmres(f, dq, maxiter_gmres, k_dim, tol, calls, k_out, i, dtol)
            ! J(q)dq=rhs=F(q)-q, with dq being the solution (denoted sol in ts_gmres).
            lin_calls = lin_calls + calls + 1  ! cumulative linear matvecs (+1 for initialize_gmres_vector)
            total_calls = total_calls + calls + 1
            tottime = tottime + calls*dt
            total_gmres_calls = total_gmres_calls + calls

            call k_sub2(q, dq) ! accepting the full step of the Newton update q = q - dq

            if (nid == 0) write (6, *) '  Outposting current solution estimate'
            time = q%time
            if (uparam(1) == 2.0) time = real(i - 1) ! Ease visualization in paraview: file 1 is t = 0
            if (uparam(1) == 2.2) time = real(i - 1)*q%time
            call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, 'nwt')
            time = q%time ! Restore

            if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
               if (nid == 0) write (6, *) '  Deallocating orbit storage'
               deallocate (uor, vor, wor)
               if (ifto .or. ldimt > 1) deallocate (tor)
            end if

         end do newton

c        Restore original tolerances from .par so subsequent DNS or
c        stability runs are not corrupted by EW's last setting.
         param(21) = saved_tol21
         param(22) = saved_tol22
         call bcast(param(21:22), 2*wdsize)

         if (nid == 0) then
            if (i == maxiter_newton) then
               write (6, *) 'Reached maxiter_newton. STOPPING! (verify convergence)'
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

      !     Output converged solution (time = orbit period)
         if (residual < dtol) then
            if (nid == 0) write (6, *) 'Outputting converged solution'
            param(63) = 1.0d0; call bcast(param(63), wdsize)  ! Double precision
            call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, "BF_")
            param(63) = 0.0d0; call bcast(param(63), wdsize)  ! Single precision
            call outpost_vort(vx, vy, vz, 'BFV')
         end if

         return
      end subroutine

      !-----------------------------------------------------------------------
      ! ts_gmres — Time-stepper GMRES solver for the Newton correction equation
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
         real :: dot_val ! For weighted inner product calculation
         real, external :: dnekclock
         integer, save :: gmres_k_sum = 0  ! Accumulator for total k values in GMRES

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
            H = 0.0d0
            yvec = 0.0d0
            evec = 0.0d0
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
      !     Note: This check is disabled by default as most problems maintain good orthogonality.
      !     Enable if you suspect loss of orthogonality in your specific problem.
      !     Warning signs:
      !     1. Slow GMRES convergence
      !     2. Unexpected increase in residuals
      !     3. Large number of GMRES restarts
!              ortho_metric = 0.0d0
!              if (k > 1) then
!                 do j = 1, k - 1
!                    call k_dot(dot_val, Q(k), Q(j))  ! Use weighted inner product
!                    ortho_metric = ortho_metric + abs(dot_val)
!                 end do
!                 ortho_metric = ortho_metric/(k - 1)  ! Average orthogonalization error
!              end if

               if (nid == 0) then
                  arnoldi_iter_time = dnekclock() - arnoldi_start_time
                  write (6, "('    ARNOLDI [GMRES ',I3,'/',I3,'] residual:',
     $   1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter, beta2, tol
                  if (k > 1) then ! Only show rate after first iteration
                     write (6, "('              Rate:',1PE15.6,' Time:',1PE15.6,
     $   's') ") beta2/prev_beta2, arnoldi_iter_time
                  else
                     write (6, "('              Time:',1PE15.6,'s') ") arnoldi_iter_time
                  end if

                  open (889, file='residu_arnoldi.dat', action='write', position='append')
                  write (889, "(I9,4(1PE15.6))") k, arnoldi_iter_time, tol, beta2, dtol
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
               gmres_k_sum = gmres_k_sum + k  ! Update accumulator after writing
               open (888, file='residu_gmres.dat', action='write', position='append')
               write (888, "(4I9,3(1PE15.6))") newton_iter, i, k, gmres_k_sum, tol, beta2, dtol
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
      ! initialize_gmres_vector — Prepare initial vector for GMRES iterations
      !
      ! This routine:
      ! 1. Computes initial residual r = b - Ax for the system Ax = b
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

         call matvec(f, q)         ! f = A*q
         call k_sub2(f, rhs)       ! f = A*q - rhs (residual calculation)
         call k_cmult(f, -1.0d0)   ! f = -(A*q - rhs) = rhs - A*q (standard residual form)
         call k_normalize(f, beta) ! Normalizes f
         call k_copy(q, f)         ! Sets q as normalized residual

         return
      end subroutine initialize_gmres_vector

      !-----------------------------------------------------------------------
      ! nonlinear_forward_map — Compute nonlinear residual f(q) = F(q) - q
      !
      ! Newton iteration requires the residual r = F(q) - q where:
      ! - F(q) is computed by time-stepping from q
      ! - The residual f = F(q) - q measures how far q is from a fixed point
      ! - For convergence, we want ||f|| -> 0
      !
      ! Steps:
      ! 1. Set initial condition from q
      ! 2. Time-step nsteps to get F(q)
      ! 3. Compute residual f = F(q) - q
      ! 4. Store base flow for linearized calculations
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
         call nopcopy(f%vx, f%vy, f%vz, f%pr, f%t, vx, vy, vz, pr, t) ! F(q)
         call k_copy(fc_nwt, f) ! for Newton UPO newton_linearized_map
         call k_sub2(f, q) ! f = F(q) - q
         f%time = 0.0d0

      !     --> Pass current guess as base flow for the linearized calculation.
         call nopcopy(ubase, vbase, wbase, pbase, tbase, q%vx, q%vy, q%vz, q%pr, q%t)

         return
      end subroutine nonlinear_forward_map

      !-----------------------------------------------------------------------
      ! set_nek5000_tolerances — Update solver tolerances
      !-----------------------------------------------------------------------
      subroutine set_nek5000_tolerances(solver_tol)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         real, intent(in) :: solver_tol ! New tolerance value to be set

         if (nid == 0) write (6, "('  TOLERANCE set from:',1PE15.6,' to:', 1PE15.6) ") param(21), abs(solver_tol)
         param(21:22) = abs(solver_tol)  ! Set both tolerances at once
         call bcast(param(21:22), 2*wdsize)  ! Broadcast both values in one call

      end subroutine set_nek5000_tolerances

      !-----------------------------------------------------------------------
      ! spec_tole — Eisenstat-Walker type 2 adaptive GMRES tolerance
      !
      ! Ref: Eisenstat & Walker, SIAM J. Sci. Comput. 17(1), 1996.
      !
      ! eta = (||F_k||/||F_{k-1}||)^alpha, alpha = 1.618 (golden ratio).
      ! Solve GMRES loosely when far from convergence, tighten as Newton
      ! converges. eta_max = 0.9 caps forcing (tol <= 0.81*||f||^2).
      !
      ! All residuals are squared norms (||f||^2). Caller must pass
      ! sqrt(tol) to set_nek5000_tolerances for norm-based param(21:22).
      !
      ! Bounds: tol >= dtol (lower); ew_tol_cap > 0 adds upper bound
      ! (set in .usr for stiff/high-Re problems, default 0 = uncapped).
      !-----------------------------------------------------------------------
      function spec_tole(residual, prev_res, dtol) result(nwtol)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         real, intent(in) :: residual ! Current ||f||^2
         real, intent(in) :: prev_res ! Previous ||f||^2
         real, intent(in) :: dtol ! Target tolerance
         real :: nwtol ! Returned new tolerance
         real :: eta ! Forcing term
         real, parameter :: eta_max = 0.9d0 ! Cap on forcing term
         real, parameter :: eta_initial = 0.5d0 ! First-iteration forcing
         real, parameter :: alpha_ew = 1.618d0 ! Golden ratio exponent

         if (prev_res > 1.0d-100 .and. residual > 0.0d0) then
c           Eisenstat-Walker type 2 forcing
c           eta = (||F_k||/||F_{k-1}||)^alpha = (res/prev_res)^(alpha/2)
            eta = (residual / prev_res) ** (alpha_ew / 2.0d0)
            eta = min(eta, eta_max)
c           GMRES tolerance: ||r||^2 < eta^2 * ||F_k||^2
            nwtol = eta**2 * residual
         else
c           First iteration: moderate relaxation
            eta = eta_initial
            nwtol = eta_initial**2 * residual
         end if

c        Lower bound: never go below Newton target
         nwtol = max(nwtol, dtol)

c        Optional upper bound (ew_tol_cap > 0 activates the cap)
         if (ew_tol_cap > 0.0d0) nwtol = min(nwtol, ew_tol_cap)

         if (nid == 0) write (6,
     $      "('  [EW: eta=',1PE10.3,' tol=',1PE10.3,
     $        ' solver=',1PE10.3,']')")
     $      eta, nwtol, sqrt(nwtol)

      end function spec_tole

      end module nekstab_newton
