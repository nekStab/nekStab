!-----------------------------------------------------------------------
! matvec.f90 — Matrix-vector products and linearized forward maps
!
! Purpose:
!   Provides the matrix-vector product dispatcher and all linearized
!   map implementations (direct, adjoint, finite-difference, transient
!   growth, force sensitivity, and Newton linearized maps).
!
! Public interface:
!   prepare_linearized_solver, matvec, forward_linearized_map,
!   forward_finite_difference_map, adjoint_linearized_map,
!   transient_growth_map, ts_force_sensitivity_map,
!   newton_linearized_map, compute_bvec
!
! Dependencies:
!   krylov_subspace, SIZE, TOTAL, ADJOINT
!-----------------------------------------------------------------------

module nekstab_matvec
    use krylov_subspace
    use nekstab_vectors
    use nekstab_nek_bridge
    implicit none
    private

    !        Cached bvec/btvec for Newton UPO (avoid recomputing per matvec)
    type(krylov_vector), save :: bvec_cache, btvec_cache
    logical, save :: bvec_cached = .false.
    ! Initialization flags for one-time setup in various subroutines
    logical, save :: matvec_init = .false.
    logical, save :: forward_linearized_map_init = .false.
    logical, save :: forward_finite_difference_map_init = .false.
    logical, save :: adjoint_linearized_map_init = .false.

   public :: prepare_linearized_solver, matvec, &
             forward_linearized_map, &
             forward_finite_difference_map, &
             adjoint_linearized_map, &
             transient_growth_map, &
             ts_force_sensitivity_map, &
             newton_linearized_map, compute_bvec, &
             cache_newton_bvec
contains

!-----------------------------------------------------------------------
! prepare_linearized_solver — Set up Nek parameters for linearized solver
!-----------------------------------------------------------------------
subroutine prepare_linearized_solver

   if (nid == 0) write (6, *) 'Preparing linearized solver...'

!  Force single perturbation mode
   if (param(31) > 1) then
      if (nid == 0) then
         write (6, *) 'ERROR: nekStab not ready for multiple perturbation modes.'
         write (6, *) 'Setting number of perturbations to 1.'
      end if
   end if
   param(31) = int(1); npert = int(param(31)) ! param is real !
   if (nid == 0) write (6, *) 'Number of perturbations set to:', npert

!  Adjust time step and number of steps if end time is specified
   if (param(10) > 0) then
      if (nid == 0) then
         write (6, *) 'End time specified:', param(10)
         write (6, *) 'Current time:', time
         write (6, *) 'Recomputing dt and nsteps to match end time...'
      end if

!     Use par-file dt if set; otherwise compute from CFL
!     param(12) < 0 means constant dt was set in the par file
      if (abs(param(12)) > 0) then
         dt = abs(param(12))
         if (nid == 0) write (6, *) 'Using dt from par file:', dt
      else
!        Compute maximum allowable time step based on CFL condition
         call compute_cfl(ctarg, vx, vy, vz, 1.0d0)
         if (nid == 0) write (6, *) 'Maximum spatial restriction:', ctarg
!        Calculate time step based on CFL target
         dt = param(26)/ctarg
      end if

!     Calculate number of steps needed to reach end time
      nsteps = ceiling(param(10)/dt)

!     Adjust time step to exactly reach end time
      dt = param(10)/nsteps

      if (nid == 0) then
         write (6, *) 'Adjusted time step dt =', dt
         write (6, *) 'Number of steps nsteps =', nsteps
         write (6, *) 'Total simulation time =', nsteps*dt
      end if

!     Update parameters
      param(12) = dt
      lastep = 0
      fintim = nsteps*dt

!     Calculate actual CFL
      call compute_cfl(ctarg, vx, vy, vz, dt)
      if (nid == 0) write (6, *) 'Actual CFL:', ctarg
   end if

!  Force constant time step
   param(12) = -abs(param(12))
   if (nid == 0) write (6, *) 'Constant time step enforced, dt =', -param(12)

!  Broadcast updated parameters to all processes
   call bcast(param, 200*wdsize) ! broadcast all params

   call nekgsync ! ensures that all processes reach before any can proceed further

   if (nid == 0) write (6, *) 'Linearized solver preparation complete.'

end subroutine prepare_linearized_solver

!-----------------------------------------------------------------------
! matvec — Dispatch matrix-vector product to the Arnoldi factorization
!
!     All subroutines need to have the same interface.
!
!     NOTE : The baseflow needs to be pass to (ubase, vbase, wbase, tbase)
!     before this function is called.
!
!     INPUTS
!     ------
!
!     qx, qy, qz, qt : nek-arrays of size lv
!     Initial velocity and temperature components.
!
!     qp : nek-array of size lp
!     Initial pressure component.
!
!     OUTPUTS
!     -------
!
!     fx, fy, fz, ft : nek-arrays of size lv
!     Final velocity and temperature components.
!
!     fp : nek-array of size lp
!     Final pressure component.
!-----------------------------------------------------------------------
subroutine matvec(f, q)

    use krylov_subspace

    type(krylov_vector), intent(out) :: f
    type(krylov_vector), intent(inout) :: q

    !     --> Pass the baseflow to vx, vy, vz
    call nopcopy(vx, vy, vz, pr, t, ubase, vbase, wbase, pbase, tbase)

    if (ifbf2d .and. if3d) then
       call rzero(vz, nx1*ny1*nz1*nelv)
       if (nid == 0) write (6, *) 'Forcing vz=0'
    end if

    !     --> Standard setup for the linearized solver.
    !     NOTE: matvec_init is reset together with forward_linearized_map_init
    !     at the end of newton_linearized_map when bvec is not cached, so that
    !     prepare_linearized_solver re-runs each Newton iteration (dt/nsteps
    !     may change). Resetting only forward_linearized_map_init would freeze
    !     the operator parameters from the first iteration.
    if (.not. matvec_init) then
       call prepare_linearized_solver
       matvec_init = .true.
    end if

   lastep = 0
   fintim = param(10)

!     --> Direct solver only steady and periodic!
   if (uparam(01) >= 3.0 .and. uparam(01) < 3.2) then
      evop = 'd'
      if (iffindiff) then
         if (nid == 0) write (*, *) "Using the finite-difference approximation of the Fréchet derivative."
         call forward_finite_difference_map(f, q)
      else
         call forward_linearized_map(f, q)
      end if
   end if

!     --> Adjoint solver only steady and periodic!
   if (uparam(01) >= 3.2 .and. uparam(01) < 3.3) then
      evop = 'a'
      call adjoint_linearized_map(f, q)
   end if

!     --> Direct-Adjoint for optimal transient growth.
   if (uparam(01) >= 3.3 .and. uparam(01) < 3.4) then
      evop = "p"
      call transient_growth_map(f, q)
   end if

!     --> Adjoint solver for the steady force sensitivity analysis.
   if (floor(uparam(01)) == 4) then
      call ts_force_sensitivity_map(f, q)
   end if

!     --> Linearized forward map for the Newton-Krylov solver.
   if (floor(uparam(01)) == 2) then
      evop = 'n'
      call newton_linearized_map(f, q)
      !  When bvec is not cached (fixed-point Newton, or before the UPO
      !  cache is populated), reset BOTH guards so that the next Newton
      !  iteration re-runs prepare_linearized_solver (fresh dt/nsteps)
      !  and forward_linearized_map recomputes the orbit.
      if (.not. bvec_cached) then
         matvec_init = .false.
         forward_linearized_map_init = .false.
      end if
   end if

end subroutine matvec

!-----------------------------------------------------------------------
! forward_linearized_map — Integrate linearized Navier-Stokes forward
!
!     Denoting by L the Jacobian of the Navier-Stokes equations, the
!     corresponding matrix vector product is thus
!
!     x(t) = exp(t * L) * x(0)
!
!     where x(0) is the initial condition (qx, qy, qz, qp, qt) and x(t)
!     the final one (fx, fy, fz, fp, ft).
!-----------------------------------------------------------------------
subroutine forward_linearized_map(f, q)

    use krylov_subspace

    type(krylov_vector), intent(out) :: f
    type(krylov_vector), intent(inout) :: q

    integer m

    nt = nx1*ny1*nz1*nelt

    !     --> Setup the parameters for the linearized solver.
    ifpert = .true.; ifadj = .false.
    call bcast(ifpert, lsize); call bcast(ifadj, lsize)

    !     --> Base flow computation control.
    !     Two-phase logic when ifstorebase is enabled:
    !       1st call  (forward_linearized_map_init=F): keep ifbase=T so
    !                 nek_advance computes the base flow; orbit_store saves it.
    !       Later calls (forward_linearized_map_init=T): set ifbase=F so
    !                 nek_advance skips it; orbit_restore replays the stored data.
    !     CAUTION: do NOT clear ifbase unconditionally before the conditional
    !     below — that would make the allocate_orbit / first-store path dead.
    ifbase = .false.
    if (uparam(01) == 3.11) ifbase = .true. ! activate Floquet
    if (uparam(01) == 3.31) ifbase = .true. ! activate Floquet for intracycle transient growth
    !  uparam==2.1/2.2 selects the UPO baseflow-evolution path specifically.
    !  Do NOT replace with isNewtonPO: that flag only signals the period-
    !  augmented inner product, while 2.1/2.2 identifies orbit storage needs.
    if (uparam(01) == 2.1 .or. uparam(01) == 2.2) then
       ifbase = .true. ! activate baseflow evolution for UPO
    end if
    if (ifstorebase) then
       if (forward_linearized_map_init) then
          ifbase = .false. ! subsequent calls: use stored orbit
       elseif (ifbase) then
          call allocate_orbit(nsteps) ! first call: allocate storage
       end if
    end if

!     --> Pass the initial condition for the perturbation.
   call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1), &
   q%vx, q%vy, q%vz, q%pr, q%t)

   time = 0.0d+00
   do istep = 1, nsteps
!     --> Output current info to logfile.
      if (nid == 0) write (6, "(' DIRECT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

!     Integrate forward in time.
      call nekstab_usrchk()
      call nek_advance()

      if (ifstorebase .and. ifbase .and. .not. forward_linearized_map_init) then !storing first time
         if (nid == 0) write (6, *) 'storing first series:', istep, '/', nsteps
         call orbit_store(istep)
      elseif (ifstorebase .and. forward_linearized_map_init .and. .not. ifbase) then !just moving in memory
         if (nid == 0) write (6, *) 'using stored baseflow'
         call orbit_restore(istep)
      end if
   end do
   if (ifstorebase .and. .not. forward_linearized_map_init .and. ifbase) then
      ifbase = .false.; forward_linearized_map_init = .true.
   end if

!     --> Copy the solution.
   call nopcopy(f%vx, f%vy, f%vz, f%pr, f%t, &
   vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1))

end subroutine forward_linearized_map

!-----------------------------------------------------------------------
! forward_finite_difference_map — Approximate linearized map via FD
!-----------------------------------------------------------------------
subroutine forward_finite_difference_map(f, q)

    use krylov_subspace

    type(krylov_vector), intent(out) :: f
    type(krylov_vector), intent(inout) :: q
    type(krylov_vector) :: pert
    type(krylov_vector) :: work

    real :: epsilon0, dummy
    integer :: i, m

    logical, save :: init = .false.

   nv = nx1*ny1*nz1*nelv

!  Setup the Nek parameters for the finite-differences approximation.
   ifpert = .false.; ifadj = .false.
   call bcast(ifpert, lsize); call bcast(ifadj, lsize)

   call k_zero(f)
   call k_zero(work)
   call nopcopy(work%vx, work%vy, work%vz, work%pr, work%t, ubase, vbase, wbase, pbase, tbase)

   call k_norm(dummy, work)
   epsilon0 = 1e-6*dummy

   if (findiff_order == 2) then
      ampls(1) = 1; ampls(2) = -1
      coefs(1) = 1; coefs(2) = -1
      coefs = coefs/2.0d+00
   else if (findiff_order == 4) then
      ampls(1) = 1; ampls(2) = -1
      ampls(3) = 2; ampls(4) = -2
      coefs(1) = 8; coefs(2) = -8
      coefs(3) = -1; coefs(4) = 1
      coefs = coefs/12.0d+00
   end if
   ampls = ampls*epsilon0

!  Turning-off the base flow side-by-side computation.
!  uparam==2.1/2.2 identifies the UPO orbit-storage path (not isNewtonPO,
!  which only controls the period-augmented inner product).
   ifbase = .false.
   if (uparam(01) == 3.11) ifbase = .true. ! activate Floquet
   if (uparam(01) == 3.31) ifbase = .true. ! activate Floquet for intracycle transient growth
   if (uparam(01) == 2.1 .or. uparam(01) == 2.2) then
      init = .true. ! Use stored baseflow if ifstorebase.
      ifbase = .true. ! activate baseflow evolution for UPO.
   end if
   if (ifstorebase .and. init) ifbase = .false. ! deactivate ifbase if baseflow stored.

   if (ifstorebase .and. ifbase .and. .not. init) then
      call allocate_orbit(nsteps)
   end if

!-----------------------------------------------------------------------
!-----                                                             -----
!-----     FINITE-DIFFERENCE APPROX. OF THE FRECHET DERIVATIVE     -----
!-----                                                             -----
!-----------------------------------------------------------------------

   do i = 1, findiff_order

!     Scale the perturbation.
      call k_copy(pert, q)
      call k_cmult(pert, ampls(i))

!     Initial condition for the each evaluation.
      call nopcopy(vx, vy, vz, pr, t, ubase, vbase, wbase, pbase, tbase)
      call nopadd2(vx, vy, vz, pr, t, pert%vx, pert%vy, pert%vz, pert%pr, pert%t)
      if (ifbf2d .and. if3d) then
         call rzero(vz, nv); if (nid == 0) write (6, *) 'Forcing vz=0'
      end if

!     Time-integration of the nonlinear Nek5000 equations.
      time = 0.0d+00
      do istep = 1, nsteps
!     --> Output current info to logfile.
         if (nid == 0) write (6, "(' DIRECT FD [',I1,'/',I1,']:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") i, &
   findiff_order, istep, nsteps, mstep, k_dim, schur_cnt

!        Nek5000 computational core.
         call nekstab_usrchk()
         call nek_advance()

         if (i == 1 .and. ifstorebase .and. ifbase .and. .not. init) then !storing first time
            if (nid == 0) write (6, *) 'storing first series:', istep, '/', nsteps
            call orbit_store(istep)
         elseif (i > 1 .and. ifstorebase .and. init .and. .not. ifbase) then !just moving in memory
            if (nid == 0) write (6, *) 'using stored baseflow'
            call orbit_restore(istep)
         end if
      end do
      if (ifstorebase .and. .not. init .and. ifbase) then
         ifbase = .false.; init = .true.
      end if

!     --> Copy the solution and compute the approximation of the Frechet derivative.
      call nopcopy(work%vx, work%vy, work%vz, work%pr, work%t, vx, vy, vz, pr, t)
      call k_cmult(work, coefs(i))
      call k_add2(f, work)

   end do

!  Rescale the approximate Frechet derivative with the step size.
   call k_cmult(f, 1.0d+00/epsilon0)

end subroutine forward_finite_difference_map

!-----------------------------------------------------------------------
! adjoint_linearized_map — Integrate adjoint Navier-Stokes forward
!
!     Denoting by L the adjoint Navier-Stokes operator, the corresponding
!     matrix vector product is thus
!
!     x(t) = exp(t * L) * x(0)
!
!     where x(0) is the initial condition (qx, qy, qz, qp, qt) and x(t)
!     the final one (fx, fy, fz, fp, ft).
!-----------------------------------------------------------------------
subroutine adjoint_linearized_map(f, q)

   use krylov_subspace

   type(krylov_vector), intent(out) :: f
   type(krylov_vector), intent(inout) :: q

   logical, save :: init
   data init/.false./

   integer m

   nt = nx1*ny1*nz1*nelt

!     --> Setup the parameters for the linearized solver.
   ifpert = .true.; ifadj = .true.
   call bcast(ifpert, lsize); call bcast(ifadj, lsize)

!     --> Turning-off the base flow side-by-side computation. Need to change for Floquet.
   ifbase = .false.
   if (uparam(01) == 3.21) ifbase = .true. ! activate floquet

   if (ifstorebase .and. init) ifbase = .false.

   if (ifstorebase .and. ifbase .and. .not. init) then
      call allocate_orbit(nsteps)
   end if
   if (uparam(01) == 3.31) init = .true. ! activate Floquet for intracycle transient growth (base flow already computed)

!     --> Pass the initial condition for the perturbation.
   call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1), &
   q%vx, q%vy, q%vz, q%pr, q%t)

   time = 0.0d+00
   do istep = 1, nsteps
!     --> Output current info to logfile.
      if (nid == 0) write (6, "(' ADJOINT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") &
   istep, nsteps, mstep, k_dim, schur_cnt

!     Integrate backward in time.
      call nekstab_usrchk()
      call nek_advance()

      if (ifstorebase .and. ifbase .and. .not. init) then !storing first time
         if (nid == 0) write (6, *) 'storing first series:', istep, '/', nsteps
         call orbit_store(istep)
      elseif (ifstorebase .and. init .and. .not. ifbase) then !just moving in memory
         if (nid == 0) write (6, *) 'using stored baseflow'
         call orbit_restore(istep)
      end if
   end do
   if (ifstorebase .and. .not. init .and. ifbase) then
      ifbase = .false.; init = .true.
   end if

!     --> Copy the solution.
   call nopcopy(f%vx, f%vy, f%vz, f%pr, f%t, &
   vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1))

end subroutine adjoint_linearized_map

!-----------------------------------------------------------------------
! transient_growth_map — Direct-adjoint composition for optimal growth
!-----------------------------------------------------------------------
subroutine transient_growth_map(f, q)

   use krylov_subspace

   type(krylov_vector), intent(out) :: f
   type(krylov_vector), intent(inout) :: q
   type(krylov_vector) :: wrk

!     --> Evaluate the forward map.
   call forward_linearized_map(wrk, q)

!     --> Evaluate the adjoint map.
   call adjoint_linearized_map(f, wrk)

end subroutine transient_growth_map

!-----------------------------------------------------------------------
! ts_force_sensitivity_map — Adjoint map for steady force sensitivity
!-----------------------------------------------------------------------
subroutine ts_force_sensitivity_map(f, q)

   use krylov_subspace

   type(krylov_vector), intent(out) :: f
   type(krylov_vector), intent(inout) :: q

!     --> Evaluate exp(t*L) * q0.
   call adjoint_linearized_map(f, q)

!     --> Evaluate (I - exp(t*L)) * q0.
   call k_sub2(f, q)
   call k_cmult(f, -1.0d+00)

end subroutine ts_force_sensitivity_map

!-----------------------------------------------------------------------
! newton_linearized_map — Linearized forward map for Newton-Krylov
!-----------------------------------------------------------------------
subroutine newton_linearized_map(f, q)

   use krylov_subspace

   type(krylov_vector), intent(out) :: f
   type(krylov_vector), intent(inout) :: q
   type(krylov_vector) :: bvec, btvec

!     ----------------------------------
!     -----     REGULAR NEWTON     -----
!     ----------------------------------

!     --> Evaluate exp(t*L) * q0.
   if (iffindiff) then
      if (nid == 0) write (*, *) "Using the finite-difference approximation of the Fréchet derivative."
      call forward_finite_difference_map(f, q)
   else
      call forward_linearized_map(f, q)
   end if

!     --> Evaluate (exp(t*L) - I) * q0.
   call k_sub2(f, q)

!     ----------------------------------
!     -----     NEWTON FOR UPO     -----
!     ----------------------------------

   if (isNewtonPO) then

      if (bvec_cached) then
!              Use cached bvec/btvec (saves 2 time-steps per matvec)
         call k_add2s2(f, bvec_cache, q%time)
         call k_dot(f%time, btvec_cache, q)
      else
!              Fallback: compute on the fly (first call or no caching)
         call k_zero(bvec)
         call k_zero(btvec)
         call compute_bvec(bvec, fc_nwt)
         call k_cmult(bvec, q%time)
         call k_add2(f, bvec)
         call compute_bvec(btvec, ic_nwt)
         call k_dot(f%time, btvec, q)
      end if

      if (nid == 0) write (6, *) &
         'Newton period correction:', f%time

   else

      f%time = 0.0d+00

   end if

end subroutine newton_linearized_map

!-----------------------------------------------------------------------
! compute_bvec — Approximate time derivative for Newton UPO
!-----------------------------------------------------------------------
subroutine compute_bvec(bvec, qbase)

   use krylov_subspace

   type(krylov_vector), intent(out) :: bvec
   type(krylov_vector), intent(in) :: qbase
   type(krylov_vector) :: wrk1, wrk2

!     --> Setup the paramtemers for the solver.
   ifpert = .false.; ifadj = .false.
   call bcast(ifpert, lsize); call bcast(ifadj, lsize)

   ifbase = .false.
   call bcast(ifbase, lsize)

!     --> Pass the initial condition.
   call nopcopy(vx, vy, vz, pr, t, qbase%vx, qbase%vy, qbase%vz, qbase%pr, qbase%t)
   param(10) = qbase%time
   call prepare_linearized_solver
   call k_copy(wrk1, qbase)

!     --> Single time-step to approximate the time-derivative.
   time = 0.0d+00
   do istep = 1, 1
      call nekStab_usrchk()
      call nek_advance()
   end do
   call nopcopy(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%t, vx, vy, vz, pr, t)

!     --> Approximate the time-derivative.
   call k_sub2(wrk2, wrk1)
   call k_copy(bvec, wrk2)
   call k_cmult(bvec, 1.0/dt)
   bvec%time = 0.0d+00

end subroutine compute_bvec

!-----------------------------------------------------------------------
! cache_newton_bvec — Precompute and cache bvec/btvec for UPO Newton
!
! Called once per Newton iteration (after nonlinear_forward_map).
! Avoids recomputing bvec/btvec at every matvec call within GMRES,
! saving 2 nonlinear time-steps per Arnoldi iteration.
!-----------------------------------------------------------------------
subroutine cache_newton_bvec(ic, fc)

   use krylov_subspace

   type(krylov_vector), intent(in) :: ic, fc
   real :: saved_param10

!        Save param(10) — compute_bvec overwrites it with qbase%time
   saved_param10 = param(10)

   call compute_bvec(bvec_cache, fc)
   call compute_bvec(btvec_cache, ic)

!        Restore param(10) and re-prepare solver state
!        (compute_bvec overwrites param(10), sets ifpert=.false.)
   param(10) = saved_param10
   call prepare_linearized_solver
   bvec_cached = .true.

   if (nid == 0) write (6, *) &
      'Newton bvec/btvec cached for GMRES'

end subroutine cache_newton_bvec

end module nekstab_matvec
