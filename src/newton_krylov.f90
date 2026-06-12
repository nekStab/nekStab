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
      !   nonlinear_forward_map, set_nek5000_tolerances,
      !   set_nek5000_velocity_tolerance, spec_tole
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL
      !-----------------------------------------------------------------------

   module nekstab_newton
   use krylov_subspace
   use nekstab_nek_bridge
   use nekstab_vectors
   use nekstab_matvec
   use nekstab_krylov_decomposition
   use nekstab_lapack
   use nekstab_diagnostics
   implicit none
   private
   public :: newton_krylov, ts_gmres, &
      initialize_gmres_vector, &
      nonlinear_forward_map, &
      set_nek5000_tolerances, set_nek5000_velocity_tolerance, spec_tole
   ! ── Log file management ──────────────────────────────────────────
   !
   ! The old code did open/write/close on every Newton, GMRES, and
   ! Arnoldi iteration — 3 open+close cycles per inner Arnoldi step.
   ! On networked HPC filesystems the metadata round-trips add up.
   !
   ! Now each log file is opened once at the start of newton_krylov
   ! and closed once at the end.  ts_gmres uses ensure_*_open so it
   ! also works when called standalone (without newton_krylov).
   integer, parameter :: newton_log_unit = 887
   integer, parameter :: gmres_log_unit = 888
   integer, parameter :: arnoldi_log_unit = 889
   contains

   !  Open (or reopen) a log file.  reset_file=.true. truncates;
   !  .false. appends (restart).
   subroutine open_log_file(unit_no, filename, reset_file)
   integer, intent(in) :: unit_no
   character(len=*), intent(in) :: filename
   logical, intent(in) :: reset_file
   logical :: opened

   inquire(unit=unit_no, opened=opened)
   if (opened) close(unit_no)

   if (reset_file) then
   open(unit=unit_no, file=filename, status='replace', action='write')
   else
   open(unit=unit_no, file=filename, action='write', position='append')
   end if

   end subroutine open_log_file

   subroutine open_newton_log_units(reset_files)
   logical, intent(in) :: reset_files

   if (nid /= 0) return

   call open_log_file(newton_log_unit, 'residu_newton.dat', reset_files)
   call open_log_file(gmres_log_unit, 'residu_gmres.dat', reset_files)
   call open_log_file(arnoldi_log_unit, 'residu_arnoldi.dat', reset_files)

   end subroutine open_newton_log_units

   subroutine close_newton_log_units()

   if (nid /= 0) return

   close(newton_log_unit)
   close(gmres_log_unit)
   close(arnoldi_log_unit)

   end subroutine close_newton_log_units

   !  ts_gmres may be called from newton_krylov (units already open)
   !  or standalone.  opened_here is an ownership flag: .true. means
   !  this call opened the log files and the caller must close them
   !  on exit (standalone case); .false. means newton_krylov owns them.
   subroutine ensure_gmres_log_units_open(opened_here)
   logical, intent(out) :: opened_here
   logical :: opened

   opened_here = .false.
   if (nid /= 0) return

   inquire(unit=gmres_log_unit, opened=opened)
   if (.not. opened) then
   call open_log_file(gmres_log_unit, 'residu_gmres.dat', .false.)
   opened_here = .true.
   end if

   inquire(unit=arnoldi_log_unit, opened=opened)
   if (.not. opened) then
   call open_log_file(arnoldi_log_unit, 'residu_arnoldi.dat', .false.)
   opened_here = .true.
   end if

   end subroutine ensure_gmres_log_units_open

   subroutine close_gmres_log_units()

   if (nid /= 0) return

   close(gmres_log_unit)
   close(arnoldi_log_unit)

   end subroutine close_gmres_log_units

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

   !-----------------------------------------------------------------------
   !  Why this routine reads isNewtonFP / isNewtonPO / isNewtonPO_T
   !  instead of `uparam(1) == 2.0/2.1/2.2`:
   !  the legacy float equality was fragile (2.1 and 2.2 are not exactly
   !  representable in IEEE-754) and relied on caller-specific tolerance
   !  conventions. The booleans are resolved once at startup by
   !  nekStab_mode_from_uparam in src/mode_config.f90 from the
   !  MODE_NEWTON_FP / MODE_NEWTON_PO / MODE_NEWTON_POT integer codes
   !  in src/mode_codes.f90.
   !-----------------------------------------------------------------------

      !     ----- Krylov vectors
   type(krylov_vector) :: f ! Right-hand side vector for Newton solver
   type(krylov_vector) :: q ! Current estimate of the solution
   type(krylov_vector) :: dq ! Newton correction obtained from GMRES

      !     ----- Iteration parameters
   integer :: i, j, maxiter_newton, maxiter_gmres, calls
   real :: residual, tol, tottime = 0.0d0
   real :: prev_residual = 0.0d0 ! Previous iteration residual for scheduling + stagnation
   real :: initial_residual = 0.0d0 ! For divergence guard
   integer :: stagnation_count = 0 ! Stagnation detection counter
   real :: newton_start_time, newton_iter_time ! Timing for Newton iterations
   integer :: total_gmres_calls ! Track GMRES calls per Newton iteration
   integer :: total_calls = 0, nonlin_calls = 0, lin_calls = 0 ! Call counters
   integer :: prev_total_calls = 0 ! Track calls from previous iteration
   integer, save :: k_out ! Store k from GMRES
   integer, save :: k_sum = 0 ! Accumulator for total k values across Newton iterations
   integer :: alloc_stat
   real :: saved_tol21, saved_tol22 ! Save/restore param(21:22)
   real :: request_gb
   real :: accepted_residuals(3) ! Non-monotone Armijo window (last 3 accepted)
   integer :: accepted_residual_count
   integer :: accepted_residual_next
   real :: armijo_reference

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
   end if
   ! Initialize from the current Nek tolerances. Dynamic scheduling below will
   ! overwrite this with a residual-proportional target when ifdyntol is active.
   tol = max(param(21), param(22))

   maxiter_newton = 30; maxiter_gmres = 30

   ! Open all 3 log files once.  istep==0 truncates (fresh run);
   ! istep>0 appends (restart).  They stay open until close_newton_log_units().
   if (nid == 0) then
   write (6, *) 'Opening output files for residuals'
   call open_newton_log_units(istep == 0)
   end if

   if (nid == 0) write (6, *) 'Initializing Krylov vectors'
   call k_zero(f); call k_zero(dq)

   if (nid == 0) write (6, *) 'Copying initial condition'
   call nopcopy(q%vx, q%vy, q%vz, q%pr, q%t, vx, vy, vz, pr, t)

!        Save original tolerances (restored after Newton exits).
!        Without this, dynamic scheduling leaves param(21:22) at whatever the last
!        Newton iteration set — corrupting any subsequent DNS or
!        stability run that reads param(21:22) for its own tolerances.
   saved_tol21 = param(21)
   saved_tol22 = param(22)
   accepted_residuals = 0.0d0
   accepted_residual_count = 0
   accepted_residual_next = 1

   newton: do i = 1, maxiter_newton
   if (nid == 0) write (6, *) '------------------------------------------------'
   newton_start_time = dnekclock()
   total_gmres_calls = 0

   if (nid == 0) write (6, "('NEWTON   - Starting iteration ',I3,'/',I3, &
      ' residual:', 1PE15.6, ' solver:', 1PE15.6, ' (target:', 1PE15.6,')') ") &
      i, maxiter_newton, prev_residual, max(param(21), param(22)), dtol
   if (i == 1) then ! Set time
   q%time = param(10) ! First guess from endTime in .par
   if (nid == 0) write (6, "('  Time    - Initial guess from param(10):', 1PE15.6) ") q%time
   else
   param(10) = q%time ! Other guesses included in nwt field
   if (nid == 0) write (6, "('  Time    - Updated from previous iteration:', 1PE15.6) ") param(10)
   end if

   call prepare_linearized_solver ! Compute nsteps and Nek parameters

   if (ifstorebase .and. (isNewtonPO .or. isNewtonPO_T)) then
   if (nid == 0) then ! Allocate nonlinear solution variable for natural or forced UPO
   write (6, "('  Allocating orbit for GMRES')")
   write (6, "('    Number of steps:',I6)") nsteps
   end if
   request_gb = 2.0d0*real(lv)*real(nsteps)*8.0d0/1.0d9
   allocate (uor(lv, nsteps), vor(lv, nsteps), stat=alloc_stat)
   if (alloc_stat /= 0) then
   if (nid == 0) write (6, *) &
      'ERROR: orbit allocation failed for uor/vor; requested GB:', request_gb
   call exitti('orbit allocation failed$', alloc_stat)
   end if
   if (if3d) then
   request_gb = real(lv)*real(nsteps)*8.0d0/1.0d9
   allocate (wor(lv, nsteps), stat=alloc_stat)
   else
   request_gb = 8.0d0/1.0d9
   allocate (wor(1, 1), stat=alloc_stat)
   end if
   if (alloc_stat /= 0) then
   if (nid == 0) write (6, *) &
      'ERROR: orbit allocation failed for wor; requested GB:', request_gb
   call exitti('orbit allocation failed$', alloc_stat)
   end if
   if (ifto .or. ldimt > 1) then
   request_gb = real(lt)*real(nsteps)*real(ldimt)*8.0d0/1.0d9
   allocate (tor(lt, nsteps, ldimt), stat=alloc_stat)
   if (alloc_stat /= 0) then
   if (nid == 0) write (6, *) &
      'ERROR: orbit allocation failed for tor; requested GB:', request_gb
   call exitti('orbit allocation failed$', alloc_stat)
   end if
   end if
   end if

   call nonlinear_forward_map(f, q) ! rhs of Newton iteration f(q)
   nonlin_calls = nonlin_calls + nsteps ! cumulative nonlinear time-steps
   total_calls = total_calls + nsteps
   tottime = tottime + time ! use time instead of nsteps*dt

!           Cache bvec/btvec for UPO (avoids recomputing per matvec)
   if (isNewtonPO) &
      call cache_newton_bvec(ic_nwt, fc_nwt)

      !     Check residual || f(q) ||! L2 norm: square root of dot product with weighted norms by bm1s
   call k_norm(residual, f) ! Computes ||f||
   residual = residual**2 ! squared L2 norm for consistent convergence check with GMRES
   if (nid == 0) write (6, *) '  Computed residual:', residual

!           Save initial residual for divergence guard
   if (i == 1) initial_residual = residual

!           A true machine-zero first residual means F(q)==q before Newton has
!           taken a step.  In practice this marks a degenerate forward map
!           (for example nsteps<=0 after dt/CFL setup), which would otherwise
!           report fake convergence for any nontrivial seed field.
   if (i == 1 .and. residual == 0.0d0) then
   if (nid == 0) write (6, *) &
      'ERROR: first Newton residual is exactly zero; degenerate forward map likely used nsteps<=0'
   call exitti('newton degenerate forward map nsteps<=0$', 1)
   end if

   if (ifnewton_backtrack) then
   accepted_residuals(accepted_residual_next) = residual
   if (accepted_residual_count < 3) accepted_residual_count = accepted_residual_count + 1
   accepted_residual_next = accepted_residual_next + 1
   if (accepted_residual_next > 3) accepted_residual_next = 1
   armijo_reference = maxval(accepted_residuals(1:accepted_residual_count))
   end if

      !     --> Outpost residual fields (optional)
   time = q%time ! adjust
   if (isNewtonFP) time = real(i - 1) ! to ease visu in paraview
   if (isNewtonPO .or. isNewtonPO_T) time = real(i - 1)*q%time ! to ease visu in paraview
   call outpost2(f%vx, f%vy, f%vz, f%pr, f%t, nof, 'res')
   time = q%time ! restore


!           Stagnation detection (before prev_residual is overwritten)
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
   write (6, "('NEWTON   - Finished iteration ',I3,'/',I3, &
      ' residual:', 1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter_newton, residual, dtol
   if (i > 1) then ! Only show rate after first iteration
   write (6, "('          Change: ',A,1PE15.6, ' Rate:', 1PE15.6) ") &
      merge('↑', '↓', residual > prev_residual), abs(residual - prev_residual), &
      residual/prev_residual
   end if
   if (stagnation_count >= 3) then
   write (6, *) 'WARNING: Newton stagnation', &
      ' (3 iterations without progress)'
   end if
   write (6, "('          Time: ',1PE15.6,'s  GMRES calls: ', I4) ") newton_iter_time, total_gmres_calls
   write (newton_log_unit, "(4I9,4(1PE15.6))") i, total_calls, total_calls - prev_total_calls, k_sum, tottime, &
      max(param(21), param(22)), residual, dtol
   prev_total_calls = total_calls ! Save for next iteration
   write (6, *) '------------------------------------------------'
   end if
   if (residual < dtol) then
   if (nid == 0) write (6, *) '  Converged. Exiting Newton loop...'
   exit newton
   end if

!           Divergence guard: abort if residual grows excessively
   if (i > 1 .and. residual > 1.0d8*initial_residual) then
   if (nid == 0) write (6, *) &
      'NEWTON: DIVERGENCE — residual exceeds', &
      ' 1e8 × initial. Aborting.'
   exit newton
   end if

!           Dynamic tolerance scheduler.  residual is ||F(q)-q||^2 and
!           ts_gmres also checks squared residuals, so keep the same units for
!           both the GMRES target and the Nek solver tolerance.
   if (ifdyntol) then
   tol = spec_tole(residual, prev_residual, dtol)
   call set_nek5000_tolerances(tol)
   end if
   prev_residual = residual ! after stagnation check + scheduler use the old value

   if (nid == 0) write (6, *) '  Solving linear system with GMRES for rhs = f = F(q) - q'
   call ts_gmres(f, dq, maxiter_gmres, k_dim, tol, calls, k_out, i, dtol)
            ! J(q)dq=rhs=F(q)-q, with dq being the solution (denoted sol in ts_gmres).
   lin_calls = lin_calls + calls + 1  ! cumulative linear matvecs (+1 for initialize_gmres_vector)
   total_calls = total_calls + calls + 1
   tottime = tottime + calls*dt
   total_gmres_calls = total_gmres_calls + calls

   if (ifnewton_backtrack) then
   call newton_backtrack(q, dq, f, residual, armijo_reference, i, &
      nonlin_calls, total_calls, tottime)
   else
   call k_sub2(q, dq) ! accepting the full step of the Newton update q = q - dq
   end if

   if (nid == 0) write (6, *) '  Outposting current solution estimate'
   time = q%time
   if (isNewtonFP) time = real(i - 1) ! Ease visualization in paraview: file 1 is t = 0
   if (isNewtonPO_T) time = real(i - 1)*q%time
   call outpost2(q%vx, q%vy, q%vz, q%pr, q%t, nof, 'nwt')
   time = q%time ! Restore

   if (ifstorebase .and. (isNewtonPO .or. isNewtonPO_T)) then
   if (nid == 0) write (6, *) '  Deallocating orbit storage'
   deallocate (uor, vor, wor)
   if (ifto .or. ldimt > 1) deallocate (tor)
   end if

   end do newton

!        Restore original tolerances from .par so subsequent DNS or
!        stability runs are not corrupted by the scheduler's last setting.
   param(21) = saved_tol21
   param(22) = saved_tol22
   call bcast(param(21:22), 2*wdsize)

   if (nid == 0) then
   if (i > maxiter_newton) then
   write (6, *) 'Reached maxiter_newton. STOPPING! (verify convergence)'
   else
   if (isNewtonFP) then
   write (6, *) 'NEWTON finished successfully after', i, 'iterations.'
   elseif (isNewtonPO) then
   write (6, *) 'NEWTON UPO finished successfully', i, 'iterations.'
   write (6, *) ' period found:', time, 1.0d0/time
   elseif (isNewtonPO_T) then
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

   call close_newton_log_units()

   end subroutine

      !-----------------------------------------------------------------------
      ! newton_backtrack — Globalized acceptance of the Newton update
      !
      ! Replaces the unconditional full step q = q - dq with damped trials
      ! q - alpha*dq, alpha = 1, 1/2, 1/4, ...  A trial is accepted by a
      ! NON-MONOTONE Armijo test against the worst of the last 3 accepted
      ! residuals (armijo_reference) rather than the previous residual:
      ! with a rough initial guess the Newton residual can legitimately
      ! rise for an iteration before falling, and a strict monotone test
      ! would strangle exactly those runs.  alpha scales the WHOLE
      ! correction, including the UPO period component dq%time.
      !
      ! Each trial costs one nonlinear forward map (nsteps time-steps);
      ! the counters are passed inout so the caller's accounting stays
      ! exact.  On acceptance, q/f/residual are updated in place and the
      ! accepted trial is the last forward map run, so base-flow and
      ! ifstorebase orbit data left behind belong to the accepted state.
      !-----------------------------------------------------------------------
   subroutine newton_backtrack(q, dq, f, residual, armijo_reference, &
      newton_iter, nonlin_calls, total_calls, tottime)

   type(krylov_vector), intent(inout) :: q ! Accepted state (updated in place)
   type(krylov_vector), intent(in) :: dq ! Newton correction from GMRES
   type(krylov_vector), intent(inout) :: f ! Residual field of the accepted state
   real, intent(inout) :: residual ! Squared L2 residual of the accepted state
   real, intent(in) :: armijo_reference ! Non-monotone acceptance reference
   integer, intent(in) :: newton_iter ! Outer iteration (logging only)
   integer, intent(inout) :: nonlin_calls, total_calls
   real, intent(inout) :: tottime

   integer, parameter :: max_backtracks = 4
   real, parameter :: c1 = 1.0d-4 ! Armijo sufficient-decrease constant
   real, parameter :: shrink = 0.5d0 ! Geometric step reduction

      ! Trial state on the heap: a krylov_vector holds full velocity,
      ! pressure and temperature fields, far too large for the stack.
   type(krylov_vector), allocatable :: q_trial, f_trial
   integer :: trial, best_trial, alloc_stat
   integer :: base_nsteps, trial_nsteps
   real :: alpha, trial_residual, best_residual, best_alpha
   real :: base_period, trial_period
   logical :: accepted, trial_finite

   allocate (q_trial, f_trial, stat=alloc_stat)
   if (alloc_stat /= 0) &
      call exitti('newton backtracking allocation failed$', alloc_stat)

   accepted = .false.
   best_residual = huge(best_residual)
   best_alpha = 0.0d0
   best_trial = -1
   base_period = q%time
   base_nsteps = nsteps
   alpha = 1.0d0

   do trial = 0, max_backtracks
   call evaluate_trial(alpha, trial_residual, trial_finite)

   if (trial_finite .and. trial_residual < best_residual) then
   best_residual = trial_residual
   best_alpha = alpha
   best_trial = trial
   end if

      ! Trial diagnostics go to stdout ONLY.  residu_newton.dat must
      ! contain one row per ACCEPTED Newton iterate: its consumers
      ! (validation/parsers/residual.py) read every row as an iterate,
      ! so rejected trials in that file would corrupt convergence plots.
   if (nid == 0) then
   if (isNewtonPO .or. isNewtonPO_T) then
   write (6, "('  BACKTRACK iter=',I3,' trial=',I2, &
      ' alpha=',1PE12.4,' period=',1PE15.6,' nsteps=',I6, &
      ' base_period=',1PE15.6,' base_nsteps=',I6, &
      ' residual=',1PE15.6,' ref=',1PE15.6)") &
      newton_iter, trial, alpha, trial_period, trial_nsteps, &
      base_period, base_nsteps, trial_residual, armijo_reference
   else
   write (6, "('  BACKTRACK iter=',I3,' trial=',I2, &
      ' alpha=',1PE12.4,' residual=',1PE15.6,' ref=',1PE15.6)") &
      newton_iter, trial, alpha, trial_residual, armijo_reference
   end if
   end if

   if (trial_finite .and. &
      trial_residual <= (1.0d0 - c1*alpha)*armijo_reference) then
   call accept_trial(trial_residual)
   exit
   end if

   alpha = alpha*shrink
   end do

      ! All trials failed the Armijo test.  Policy: take the best finite
      ! trial that still decreased the residual (hard cases keep moving);
      ! abort loudly only when every trial made things worse — silently
      ! accepting growth is how the old full-step divergences happened.
   if (.not. accepted .and. best_residual < residual) then
      ! Re-run the forward map at best_alpha: the arrays currently hold
      ! the LAST trial, and base-flow/orbit data stored during the map
      ! must belong to the state actually accepted.
   call evaluate_trial(best_alpha, trial_residual, trial_finite)
   if (nid == 0) write (6, "('  BACKTRACK fallback accepted trial=',I2, &
      ' alpha=',1PE12.4,' residual=',1PE15.6)") &
      best_trial, best_alpha, trial_residual
   call accept_trial(trial_residual)
   end if

   if (.not. accepted) then
   if (nid == 0) write (6, "('ERROR: Newton backtracking found no ', &
      'decreasing trial; base=',1PE15.6,' best=',1PE15.6, &
      ' trials=',I2)") residual, best_residual, max_backtracks + 1
   call exitti('newton backtracking no decreasing trial$', 1)
   end if

   deallocate (q_trial, f_trial)

   contains

      ! Build q - alpha_in*dq and measure its true nonlinear residual.
      ! Shared by the search loop and the fallback so the two paths can
      ! never drift apart.
   subroutine evaluate_trial(alpha_in, res_out, finite_out)
   real, intent(in) :: alpha_in
   real, intent(out) :: res_out
   logical, intent(out) :: finite_out

   call k_copy(q_trial, q)
   call k_add2s2(q_trial, dq, -alpha_in) ! damps dq%time too (UPO period)
   trial_period = q_trial%time
   trial_nsteps = nsteps

   if (isNewtonPO .or. isNewtonPO_T) then
      ! A non-finite or non-positive trial period must never reach Nek.
   if (.not. (is_finite(q_trial%time) .and. q_trial%time > 0.0d0)) then
   res_out = huge(res_out)
   finite_out = .false.
   trial_nsteps = -1
   return
   end if
   param(10) = q_trial%time
   call prepare_linearized_solver
   trial_nsteps = nsteps
   call ensure_trial_orbit_storage(nsteps)
   end if

   call nonlinear_forward_map(f_trial, q_trial)
   nonlin_calls = nonlin_calls + nsteps
   total_calls = total_calls + nsteps
   tottime = tottime + time
   call k_norm(res_out, f_trial)
   res_out = res_out**2 ! squared L2, same units as every Newton check
   finite_out = is_finite(res_out)
   end subroutine evaluate_trial

   subroutine ensure_trial_orbit_storage(required_nsteps)
   integer, intent(in) :: required_nsteps
   integer :: allocated_nsteps

   if (.not. ifstorebase) return
   if (.not. (isNewtonPO .or. isNewtonPO_T)) return

   allocated_nsteps = 0
   if (allocated(uor)) allocated_nsteps = size(uor, 2)
   if (allocated_nsteps >= required_nsteps) return

      ! q_trial%time = q%time - alpha*dq%time is linear in alpha.
      ! Depending on the sign of dq%time, a damped UPO trial can be longer
      ! than the orbit arrays allocated for the iteration's original period,
      ! so grow the storage before nonlinear_forward_map writes snapshots.
   if (nid == 0) write (6, "('  BACKTRACK reallocating orbit storage from ',I6, &
      ' to ',I6,' steps')") allocated_nsteps, required_nsteps
   call allocate_orbit(required_nsteps)

   end subroutine ensure_trial_orbit_storage

   subroutine accept_trial(res_in)
   real, intent(in) :: res_in
   call k_copy(q, q_trial)
   call k_copy(f, f_trial)
   residual = res_in
   accepted = .true.
   end subroutine accept_trial

   logical function is_finite(x)
   real, intent(in) :: x
      ! x == x is false only for NaN (no IEEE intrinsics needed);
      ! the huge() bound rejects +/-Inf.
   is_finite = x == x .and. abs(x) < huge(x)
   end function is_finite

   end subroutine newton_backtrack

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
   integer, save :: gmres_k_sum = 0  ! Accumulator for total k values in GMRES
   logical :: opened_local_logs  ! true if WE opened the logs (standalone call)

   calls = 0
   opened_local_logs = .false.

   ! If called from newton_krylov, logs are already open (opened_local_logs stays false).
   ! If called standalone, this opens them and sets opened_local_logs=true
   ! so we close them at the end of this subroutine.
   call ensure_gmres_log_units_open(opened_local_logs)

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

   if (nid == 0) write (6, "('  GMRES   - Starting iteration ',I3,'/',I3, &
      ' residual:', 1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter, beta2, tol

      !     Reset arrays for new restart
   H = 0.0d0
   yvec = 0.0d0
   evec = 0.0d0
   evec(1) = beta
   do j = 2, ksize + 1
   call k_zero(Q(j))
   end do

   arnoldi: do k = 1, k_dim

   if (nid == 0) write (6, "('    ARNOLDI [GMRES ',I3,'/',I3,'] Starting iteration ', I3, '/', I3) ") &
      i, maxiter, k, ksize
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
   write (6, "('    ARNOLDI [GMRES ',I3,'/',I3,'] residual:', &
      1PE15.6, ' (target:', 1PE15.6, ')') ") i, maxiter, beta2, tol
   if (k > 1) then ! Only show rate after first iteration
   write (6, "('              Rate:',1PE15.6,' Time:',1PE15.6, &
      's') ") beta2/prev_beta2, arnoldi_iter_time
   else
   write (6, "('              Time:',1PE15.6,'s') ") arnoldi_iter_time
   end if

   write (arnoldi_log_unit, "(I9,4(1PE15.6))") k, arnoldi_iter_time, tol, beta2, dtol
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
   write (6, "('           Rate:',1PE15.6,' Time:',1PE15.6, &
      's  Arnoldi steps:', I4) ") beta2/prev_beta2, gmres_iter_time, k
   else
   write (6, "('           Time:',1PE15.6,'s  Arnoldi steps:', I4) ") &
      gmres_iter_time, k
   end if
   gmres_k_sum = gmres_k_sum + k  ! Update accumulator after writing
   write (gmres_log_unit, "(4I9,3(1PE15.6))") newton_iter, i, k, gmres_k_sum, tol, beta2, dtol

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

   if (opened_local_logs) call close_gmres_log_units()

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
   if (ifstorebase .and. (isNewtonPO .or. isNewtonPO_T)) then
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
   f%time = q%time
   call k_copy(fc_nwt, f) ! for Newton UPO newton_linearized_map
   call k_sub2(f, q) ! f = F(q) - q
   f%time = 0.0d0

      !     --> Pass current guess as base flow for the linearized calculation.
   call nopcopy(ubase, vbase, wbase, pbase, tbase, q%vx, q%vy, q%vz, q%pr, q%t)

   end subroutine nonlinear_forward_map

      !-----------------------------------------------------------------------
      ! set_nek5000_tolerances — Update solver tolerances
      !-----------------------------------------------------------------------
   subroutine set_nek5000_tolerances(solver_tol)

   real, intent(in) :: solver_tol ! New tolerance value to be set

   if (nid == 0) write (6, "('  TOLERANCE set from:',1PE15.6,' to:', 1PE15.6) ") param(21), abs(solver_tol)
   param(21:22) = abs(solver_tol)  ! Set pressure and velocity tolerances

   ! Nek5000's velocity Helmholtz path does not read param(22) directly in all
   ! cases.  It first consults restol(ifield), which is initialized from the
   ! input file and otherwise stays stale unless we refresh it here.  Without
   ! this update, the scheduler appears to change the tolerance on paper while
   ! the actual velocity solves keep running with the old fixed value.
   restol = abs(solver_tol)

   call bcast(param(21:22), 2*wdsize)
   call bcast(restol, (ldimt1 + 1)*wdsize)

   end subroutine set_nek5000_tolerances

      !-----------------------------------------------------------------------
      ! set_nek5000_velocity_tolerance — Update velocity Helmholtz tolerance
      !
      ! SFD only needs to relax the velocity Helmholtz solves. Keeping the
      ! pressure tolerance fixed avoids perturbing the pressure/projection path.
      !-----------------------------------------------------------------------
   subroutine set_nek5000_velocity_tolerance(solver_tol)

   real, intent(in) :: solver_tol

   if (nid == 0) write (6, &
      "('  VELOCITY TOLERANCE set from:',1PE15.6,' to:',1PE15.6)") &
      param(22), abs(solver_tol)

   param(22) = abs(solver_tol)
   restol(1) = abs(solver_tol)

   call bcast(param(22), wdsize)
   call bcast(restol(1), wdsize)

   end subroutine set_nek5000_velocity_tolerance

      !-----------------------------------------------------------------------
      ! spec_tole — residual-proportional adaptive GMRES tolerance
      !
      ! Keep the inner linear/flow solves below the current nonlinear
      ! residual without forcing them a full decade tighter:
      !
      !     nwtol = 0.2 * ||F(q)-q||^2
      !
      ! The Newton and GMRES convergence checks in this module use squared
      ! norms, so nwtol intentionally remains in squared-residual units.
      !
      ! Bounds: tol >= dtol (lower); ew_tol_cap > 0 adds upper bound
      ! (legacy variable name, set in .usr for stiff/high-Re problems,
      ! default 0 = uncapped).
      !-----------------------------------------------------------------------
   function spec_tole(residual, prev_res, dtol) result(nwtol)

   real, intent(in) :: residual ! Current ||f||^2
   real, intent(in) :: prev_res ! Retained for interface stability; unused
   real, intent(in) :: dtol ! Target tolerance
   real :: nwtol ! Returned new tolerance
   real, parameter :: residual_fraction = 0.2d0

   if (residual > 0.0d0 .and. residual == residual) then
   nwtol = residual_fraction*residual
   else
   nwtol = dtol
   end if

!        Lower bound: never go below Newton target
   nwtol = max(nwtol, dtol)

!        Optional upper bound (ew_tol_cap > 0 activates the cap)
   if (ew_tol_cap > 0.0d0) nwtol = min(nwtol, ew_tol_cap)

   if (nid == 0) write (6, &
      "('  [DYN: residual=',1PE10.3,' tol=',1PE10.3, &
      ' cap=',1PE10.3,']')") &
      residual, nwtol, ew_tol_cap

   end function spec_tole

   end module nekstab_newton
