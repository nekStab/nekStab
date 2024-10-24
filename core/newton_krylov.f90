      subroutine newton_krylov()
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

!     ----- Krylov vectors
      type(krylov_vector) :: f  ! Right-hand side vector for Newton solver
      type(krylov_vector) :: q  ! Current estimate of the solution
      type(krylov_vector) :: dq ! Newton correction obtained from GMRES

!     ----- Iteration parameters
      integer :: i, j, k, maxiter_newton, maxiter_gmres, calls
      real :: tol, residual, tottime

!     ----- Counters and flags
      integer, save :: calls_counter = 0
      logical, save :: dyntolinit = .false.
      real, save :: dtol = 0.0d0

!     ----- Initialize tolerance
      if (dtol == 0.0d0) then
         dtol = max(param(21), param(22))
         if (nid == 0) then
            write(6, *) 'Initializing tolerance:'
            write(6, *) '  User specified tolerance:', dtol
         end if
      end if
      tol = max(param(21), param(22))

!     ----- Set iteration limits
      maxiter_newton = 100
      maxiter_gmres = 100

!     ----- Open output files
      if (istep == 0 .and. nid == 0) then
         write(6, *) 'Opening output files for residuals'
         
!     File for overall residual data
!     Columns: calls_counter, tottime, residual
         open(unit=886, file='residu.dat', status='replace')
         
!     File for Newton iteration residuals
!     Columns: i (iteration number), residual
         open(unit=887, file='residu_newton.dat', status='replace')
         
!     File for GMRES iteration residuals (opened and closed in ts_gmres subroutine)
!     Columns: iteration, residual (written in ts_gmres)
         open(unit=888, file='residu_gmres.dat', status='replace')
         close(888)
         
!     File for Arnoldi iteration data (opened and closed in arnoldi subroutine)
!     Columns: iteration, data (written in arnoldi)
         open(unit=889, file='residu_arnoldi.dat', status='replace')
         close(889)
      end if

!     ----- Initialize arrays
      if (nid == 0) write(6, *) 'Initializing Krylov vectors'
      call k_zero(f)
      call k_zero(dq)

!     ----- Copy initial condition
      if (nid == 0) write(6, *) 'Copying initial condition'
      call nopcopy(q%vx, q%vy, q%vz, q%pr, q%t, vx, vy, vz, pr, t)

!     ----- Newton iteration
      if (nid == 0) write(6, *) 'Starting Newton iterations'
      newton: do i = 1, maxiter_newton
      if (nid == 0) write(6, *) 'Newton iteration:', i

!     Set time
      if (i == 1) then
         q%time = param(10)     ! First guess from endTime in .par
         if (nid == 0) write(6, *) '  Setting initial time:', q%time
      else
         param(10) = q%time     ! Other guesses included in nwt field
         if (nid == 0) write(6, *) '  Updating time:', q%time
      end if

!     Setup/Update the nek-parameters for the Newton solver
      if (nid == 0) write(6, *) '  Preparing linearized solver'
      call prepare_linearized_solver ! Compute nsteps

!     Outpost current estimate of the solution
      if (nid == 0) write(6, *) '  Outposting current solution estimate'
      time = q%time
      if (uparam(1) == 2.0) time = real(i) ! Ease visualization in paraview
      if (uparam(1) == 2.2) time = real(i) * q%time
      call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, 'nwt')
      time = q%time             ! Restore

!     Allocate nonlinear solution variable for natural or forced UPO
      if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
      if (nid == 0) then
         write(6, *) '  Allocating orbit for GMRES'
         write(6, *) '    Number of steps:', nsteps
      end if
      allocate(uor(lv, nsteps), vor(lv, nsteps))
      if (if3d) then
         allocate(wor(lv, nsteps))
      else
         allocate(wor(1, 1))
      end if
      if (ifto .or. ldimt > 1) allocate(tor(lt, nsteps, ldimt))
      end if

!     Variable tolerances for speed up
      if (ifdyntol .and. dyntolinit) then
         if (nid == 0) write(6, *) '  Adjusting tolerance dynamically'
         call spec_tole(residual, dtol)
      end if
      tol = max(param(21), param(22))
      if (nid == 0) write(6, *) '  Current tolerance:', tol

!     Compute rhs of Newton iteration f(q)
      if (nid == 0) write(6, *) '  Computing nonlinear forward map'
      call nonlinear_forward_map(f, q)
      calls_counter = calls_counter + nsteps
      tottime = tottime + nsteps * dt

!     Check residual || f(q) ||
      if (nid == 0) write(6, *) '  Calculating residual'
      call k_norm(residual, f)
      residual = residual ** 2

!     --> Outpost residual fields (optional)
      time = q%time             ! adjust
      if (uparam(1) == 2.0) time = real(i) ! to ease visu in paraview
      if (uparam(1) == 2.2) time = real(i)*q%time ! to ease visu in paraview
      call outpost2(f%vx, f%vy, f%vz, f%pr, f%t, nof, 'res')
      time = q%time             ! restore
      
!     Output iteration information
      if (nid == 0) then
         write(6, "('  NEWTON  - Iteration:',I3,'/',I3,' residual:',E15.7)") i, maxiter_newton, residual
         write(6, *) '             Tolerance target:', tol
         write(886, "(I6,2E15.7)") calls_counter, tottime, residual
         write(887, "(I6,1E15.7)") i, residual
      end if

!     Check for convergence
      if (residual < dtol) then
         if (nid == 0) write(6, *) '  Convergence achieved, exiting Newton loop'
         exit newton
      end if

!     Adjust tolerance if needed
      if (ifdyntol .and. .not. dyntolinit) then
         if (nid == 0) write(6, *) '  Initializing dynamic tolerance'
         call spec_tole(residual, dtol)
         dyntolinit = .true.
         tol = max(param(21), param(22))
      end if

!     Solve the linear system
      if (nid == 0) write(6, *) '  Solving linear system with GMRES'
      call ts_gmres(f, dq, maxiter_gmres, k_dim, calls)
      calls_counter = calls_counter + calls
      tottime = tottime + calls * dt

!     Update Newton solution
      if (nid == 0) write(6, *) '  Updating Newton solution'
      call k_sub2(q, dq)

!     Deallocate nonlinear solution variable
      if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
      if (nid == 0) write(6, *) '  Deallocating orbit storage'
      deallocate(uor, vor, wor)
      if (ifto .or. ldimt > 1) deallocate(tor)
      end if
      end do newton

!     Output final results
      if (nid == 0 .and. i <= maxiter_newton) then
         close(886)
         close(887)
         if (i == maxiter_newton) then
            write(6, *) 'Reached maxiter_newto. STOPPING! (verify convergence)'
         else
            if (uparam(1) == 2) then
               write(6, *) 'NEWTON finished successfully after', i, 'iterations.'
            elseif (uparam(1) == 2.1) then
               write(6, *) 'NEWTON UPO finished successfully', i, 'iterations.'
               write(6, *) ' period found:', time, 1.0d0/time
            elseif (uparam(1) == 2.2) then
               write(6, *) 'NEWTON for forced UPO finished successfully', i, 'iterations.'
               write(6, *) ' period found:', time, 1.0d0/time
            end if
            write(6, *) 'Calls to the linearized solver: ', calls_counter
            write(6, *) 'Total nondimensional time:', tottime
            if (ifdyntol) write(6, *) 'ifdyntol active!'
         end if
      end if

!     Output solution if converged
      if (residual < dtol) then
         if (nid == 0) write(6, *) 'Outputting converged solution'
         param(63) = 1          ! Enforce 64-bit output
         call bcast(param, 200*wdsize)
         call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, "BF_")
         param(63) = 0          ! Enforce 32-bit output
         call bcast(param, 200*wdsize)
         call outpost_vort(vx, vy, vz, 'BFV')
      end if

      if (nid == 0) write(6, *) 'newton_krylov subroutine completed'
      return
      end subroutine
      
      !-----------------------------------------------------------------------
      
      subroutine ts_gmres(rhs, sol, maxiter, ksize, calls)
      
      !     Implementation of simple time-stepper GMRES to be part of the Newton-Krylov solver
      !     for fixed point computation. The rank of the Krylov subspace is set as the user parameter k_dim.
      !
      !     INPUT
      !     -----
      !
      !     rhs_x, rhs_y, rhs_z, rhs_t : nek arrays of size (lv).
      !     Arrays containing the right-hand side of the linear problem to be solved.
      !
      !     rhs_p : nek array of size (lp)
      !     Array containing the right-hand side of the linear problem to be solved (pressure component).
      !
      !     maxiter : integer
      !     Maximum number of restarts for the GMRES computation.
      !
      !     ksize : integer
      !     Dimension of the Krylov subspace.
      !
      !     RETURNS
      !     -------
      !
      !     sol_x, sol_y, sol_z, sol_t : nek arrays of size (lv).
      !     Arrays containing the solution of the linear problem.
      !
      !     sol_p : nek array of size (lp).
      !     Array containing the solution of the linear problem (pressure component).
      !
      !     calls : total number of calls to the linearized solver.
      !     Integer containing the total sum of nsteps*k necessary to converge the solution.
      !
      !     NOTE : This is a plain implementation of GMRES following the algorithm given in
      !     Y. Saad. Iterative methods for sparse linear systems. Section 6.5 GMRES, alg. 6.9
      !
      !     Last Edit : March 26th 2021 by JC Loiseau.
      !
      
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         integer :: ksize
      
      !     ----- Right-hand side vector of A x = b -----
         type(krylov_vector) :: rhs
      
      !     ----- GMRES solution vector.
         type(krylov_vector) :: sol, dq
      
      !     ----- Krylov basis for the Arnoldi factorization.
         type(krylov_vector), dimension(:), allocatable :: Q
      
      !     ----- Upper Hessenberg matrix.
         real, allocatable, dimension(:, :) :: H
         real, allocatable, dimension(:) :: yvec, evec
      
      !     ----- Miscellaneous.
         integer :: i, k, maxiter, calls
         real :: beta, tol
      
         tol = max(param(21), param(22))
      
      !     ----- Allocate arrays -----
         allocate (Q(ksize + 1), H(ksize + 1, ksize), yvec(ksize), evec(ksize + 1))
      
         call k_zero(sol)
         call k_zero(Q(1:ksize + 1))
         H = 0.0d+00; yvec = 0.0d+00; evec = 0.0d+00
      
         call k_copy(Q(1), rhs)
         call k_normalize(Q(1), beta)
      
         calls = 0
         gmres: do i = 1, maxiter
      
      !     --> Zero-out stuff.
            H = 0.0d+00; yvec = 0.0d+00; evec = 0.0d+00; evec(1) = beta; call k_zero(Q(2:ksize + 1))
      
            arnoldi: do k = 1, k_dim
      
               call arnoldi_factorization(Q, H, k, k, ksize)
      
      !     --> Least-squares problem.
               call lstsq(H(1:k + 1, 1:k), evec(1:k + 1), yvec(1:k), k + 1, k)
      
      !     --> Compute residual.
               beta = norm2(evec(1:k + 1) - matmul(H(1:k + 1, 1:k), yvec(1:k)))
      
               if (nid == 0) then
                  open (889, file='residu_arnoldi.dat', action='write', position='append')
                  write (6, "(' ARNOLDI --- Iteration:',I5,'/',I5,' residual:',E15.7)") k, ksize, beta**2
                  write (889, "(I6,1E15.7)") k, beta**2; close (889)
               end if
      
               if (beta**2 < tol) then ! count of calls to linearized solver
                  calls = calls + k*nsteps
                  exit arnoldi
               end if
      
      ! --> Relaxed exit condition if finite-difference approximation of the operator is considered.
               if ((iffindiff) .and. (beta**2 < 1e-8)) then ! count of calls to linearized solver
                  calls = calls + k*nsteps
                  exit arnoldi
               end if
            end do arnoldi
      
      !     --> Update solution.
            call k_matmul(dq, Q(1:k), yvec(1:k), k)
            call k_add2(sol, dq)
      
      !     --> Recompute residual for sanity check and initialize new Krylov seed if needed.
            call k_copy(Q(1), sol)
            call initialize_gmres_vector(beta, Q(1), rhs)
      
            if (nid == 0) then
               open (888, file='residu_gmres.dat', action='write', position='append')
               write (6, "(' GMRES   -- Iteration:',I4,'/',I4,' residual:',E15.7)") i, maxiter, beta**2
               write (888, "(I6,1E15.7)") i, beta**2; close (888)
            end if
      
            if (beta**2 < tol) exit gmres
            if ((iffindiff) .and. (beta**2 < 1e-6)) exit gmres
         end do gmres
      
      !     ----- Deallocate arrays -----
         deallocate (Q, H, yvec, evec)
      
      end subroutine ts_gmres
      
      !-----------------------------------------------------------------------
      
      subroutine initialize_gmres_vector(beta, q, rhs)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Right-hand side vector of A x = b (input) -----
         type(krylov_vector) :: rhs
      
      !     ----- Seed for the GMRES Krylov subspace (output) -----
         type(krylov_vector) :: q, f
         real :: beta
      
      !     --> Initial Krylov vector.
         call matvec(f, q)
         call k_sub2(f, rhs)
         call k_cmult(f, -1.0d+00)
      
      !     --> Normalize the starting vector.
         call k_normalize(f, beta)
         call k_copy(q, f)
      
         return
      end subroutine initialize_gmres_vector
      
      !-----------------------------------------------------------------------
      
      subroutine nonlinear_forward_map(f, q)
         use krylov_subspace
      
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Right-hand side of the Newton.
         type(krylov_vector) :: f
      !     ----- Initial condition for the forward simulation.
         type(krylov_vector) :: q
      
         integer m
         nt = nx1*ny1*nz1*nelt
      
      !     --> Copy the initial condition to Nek.
         call k_copy(ic_nwt, q) ! for Newton UPO newton_linearized_map
         call nopcopy(vx, vy, vz, pr, t, q%vx, q%vy, q%vz, q%pr, q%t)
      
      !     --> Turn-off the linearized solver.
         ifpert = .false.; call bcast(ifpert, lsize)
      
      !     --> Run the simulation forward.
         time = 0.0d0
         do istep = 1, nsteps
            call nekStab_usrchk()
            call nek_advance()
            if (ifstorebase .and. (uparam(1) == 2.1 .or. uparam(1) == 2.2)) then
               if (nid == 0) write (6, *) 'Storing nonlinear solution for GMRES:', istep, '/', nsteps
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
         f%time = 0.0d+00
      
      !     --> Pass current guess as base flow for the linearized calculation.
         call nopcopy(ubase, vbase, wbase, pbase, tbase, q%vx, q%vy, q%vz, q%pr, q%t)
      
         return
      end subroutine nonlinear_forward_map
      
      !-----------------------------------------------------------------------
      
      subroutine set_solv_tole(new_tol) ! set solver tolerances
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, intent(in) :: new_tol  ! New tolerance value to be set
         if (nid == 0) write (6, *) 'ifdyntol Changing tol from', param(21), 'to', abs(new_tol)
      
      ! Set both param(21) and param(22) to the absolute value of the new tolerance.
      ! Broadcast these changes to all nodes.
         param(21) = abs(new_tol); call bcast(param(21), wdsize)
         param(22) = abs(new_tol); call bcast(param(22), wdsize)
      
      ! Update TOLPDF and TOLHDF with the new tolerance and broadcast the changes.
         TOLPDF = param(21); call bcast(TOLPDF, wdsize)
         TOLHDF = param(22); call bcast(TOLHDF, wdsize)
      
      ! Update restol and atol with the new tolerance and broadcast the changes.
         restol(:) = param(22); call bcast(restol, (ldimt1 + 1)*wdsize)
         atol(:) = param(22); call bcast(atol, (ldimt1 + 1)*wdsize)
      
      end subroutine set_solv_tole
      
      !-----------------------------------------------------------------------
      
      subroutine spec_tole(residual, dtol) ! Tolerance scheduling
      ! Progressively tighten tolerances to minimize computational time
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, intent(in) :: residual, dtol
         real :: nwtol
      
         nwtol = 10**(log10(residual) - 1) ! Always two decades smaller
      
         if (nid == 0) then
            write (6, *) 'Current residual:', residual
            write (6, *) 'New tolerance:', nwtol
            write (6, *) 'Target tolerance:', dtol
         end if
      
         if (nwtol <= dtol) then
            call set_solv_tole(dtol)
            if (nid == 0) write (6, *) 'Forcing user specified tolerance:', dtol
         else
            nwtol = min(nwtol, 1e-4) ! Never exceed a tolerance of 1e-4
            call set_solv_tole(nwtol)
            if (nwtol == 1e-4 .and. nid == 0) then
               write (6, *) 'Forcing minimal tolerances:', nwtol
            end if
         end if
      
      end subroutine spec_tole
