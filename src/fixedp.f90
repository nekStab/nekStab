!-----------------------------------------------------------------------
! fixedp.f90 — Fixed-point solvers for base flow computation
!
! Purpose:
!   Implements iterative methods to compute unstable steady states
!   (base flows) of the Navier-Stokes equations: Time-Delayed
!   Feedback (TDF), Selective Frequency Damping (SFD), and
!   BoostConv acceleration.
!
! Public interface:
!   tdf           — Time-Delayed Feedback stabilization
!   SFD           — Selective Frequency Damping
!   BoostConv     — BoostConv convergence acceleration
!   boostconv_core — BoostConv inner iteration
!   qr_dec        — QR decomposition for BoostConv
!   linear_system — Back-substitution for upper triangular system
!
! Dependencies:
!   krylov_subspace, SIZE, TOTAL
!-----------------------------------------------------------------------

module nekstab_fixedpoint
   use nekstab_nek_bridge, only: nekStab_error, nid, uparam, param, ctarg, &
                                 vx, vy, vz, pr, t, dt, time, istep, nsteps, &
                                 if3d, ifto, ifpsco, ldimt, bm1, fcx, fcy, fcz, &
                                 fct, nx1, ny1, nz1, nelv, nelt, nfield, &
                                 vxlag, vylag, vzlag, ifbfcv, ifdyntol, lsize, &
                                 wdsize, nof, ew_tol_cap, bst_snp, bst_skp, &
                                 NEKSTAB_UNIT_RESIDU, NEKSTAB_UNIT_DYNTOL
   use nekstab_vectors, only: nopcopy, nopsub2, opadd3
   use nekstab_diagnostics, only: outpost_vort
   use nekstab_newton, only: set_nek5000_velocity_tolerance
   implicit none
   private
   logical, save :: tdf_init = .false.
   logical, save :: sfd_init = .false. ! prevents double-open of unit 10 in SFD
   logical, save :: boostconv_init = .false.
   logical, save :: boostconv_core_init = .false.
   real, save :: SFD_oldRes, SFD_dtol
   real, parameter :: sfd_tol_residual_factor_default = 0.9d0
   real, parameter :: sfd_tol_update_rtol = 1.0d-12
   integer, parameter :: sfd_tol_update_stride_default = 100
   real, save :: sfd_tol_residual_factor = sfd_tol_residual_factor_default
   integer, save :: sfd_tol_update_stride = sfd_tol_update_stride_default
   real :: sfd_solver_tol

   public :: tdf, SFD, BoostConv, boostconv_core, &
             qr_dec, linear_system
contains

   subroutine configure_sfd_dynamic_tolerance()

      sfd_tol_residual_factor = sfd_tol_residual_factor_default
      if (uparam(8) > 0.0d0) sfd_tol_residual_factor = uparam(8)

      sfd_tol_update_stride = sfd_tol_update_stride_default
      if (uparam(9) > 0.0d0) sfd_tol_update_stride = max(1, nint(uparam(9)))

      if (nid == 0 .and. ifdyntol) write (6, &
                                          "('  SFD dynamic tolerance: factor=',1PE10.3,' stride=',I0)") &
         sfd_tol_residual_factor, sfd_tol_update_stride

   end subroutine configure_sfd_dynamic_tolerance

   !-----------------------------------------------------------------------
   ! apply_feedback_forcing — add gain*d to the global forcing arrays for
   ! ALL fields: velocity (fcx/fcy/fcz) and every scalar (fct(:,m),
   ! m=1..nfield-1: temperature + passive scalars such as RANS tke/tau).
   ! Shared by SFD/TDF/BoostConv so every feedback method drives the full
   ! state, not just velocity.  No-op on scalars for velocity-only cases
   ! (nfield==1).  fct is zeroed each step by zero_forcing, so this adds
   ! cleanly without accumulating.
   !-----------------------------------------------------------------------
   subroutine apply_feedback_forcing(d, gain)
      use nekstab_krylov_subspace, only: krylov_vector
      type(krylov_vector), intent(in) :: d
      real, intent(in) :: gain
      integer :: m, ntotv, ntott
      ntotv = nx1*ny1*nz1*nelv
      ntott = nx1*ny1*nz1*nelt
      call add2s2(fcx, d%vx, gain, ntotv)
      call add2s2(fcy, d%vy, gain, ntotv)
      if (if3d) call add2s2(fcz, d%vz, gain, ntotv)
      do m = 1, nfield - 1
         call add2s2(fct(1, 1, 1, 1, m), d%t(1, m), gain, ntott)
      end do
   end subroutine apply_feedback_forcing

   !  Safely close a Fortran unit if it is open.
   !  Called from residual/scheduler cleanup paths so repeated close
   !  attempts remain harmless.
   !  Note: historical path from v1.0-style unit hygiene in dispatchers.
   subroutine close_if_open(unit_id)
      integer, intent(in) :: unit_id
      logical :: opened

      if (nid /= 0) return

      inquire (unit=unit_id, opened=opened)
      if (opened) close (unit_id)

   end subroutine close_if_open

   !  Safely close residu.dat if it is open.
   subroutine close_residu_unit()
      call close_if_open(NEKSTAB_UNIT_RESIDU)

   end subroutine close_residu_unit

   !  Safely close the dynamic-tolerance scheduler log if open.
   subroutine close_dyn_tol_unit()
      call close_if_open(NEKSTAB_UNIT_DYNTOL)

   end subroutine close_dyn_tol_unit

   !-----------------------------------------------------------------------
   ! tdf — Time-Delayed Feedback stabilization
   !
   ! Purpose:
   !   Stabilizes unstable periodic orbits by applying a feedback
   !   force proportional to the difference between the current
   !   state and the state one period earlier.
   !-----------------------------------------------------------------------
   subroutine tdf
      use nekstab_krylov_subspace, only: lv, lt, nv, nt, uor, vor, wor, tor, NEKSTAB_PI, &
                                         krylov_vector, k_zero, k_norm

      real, allocatable, save :: do1(:), do2(:), do3(:)
      real :: h1, semi, l2, linf, rate, tol
      real :: request_gb
      real, save :: residu0, gain, porbit
      integer, save :: i, norbit, m, ibuf
      integer :: alloc_stat
      type(krylov_vector) :: tdfD

      if (.not. tdf_init) then

         if (nid == 0) write (6, *) '   Target frequency specified in uparam(5)=', uparam(5)
         porbit = 1.0d0/uparam(5) ! user input from .usr -> also used on the inflow BC

         call compute_cfl(ctarg, vx, vy, vz, 1.0d0) ! ctarg contains the sum ( ux_i / dx_i )
         dt = param(26)/ctarg ! dt given target CFL
         norbit = ceiling(porbit/dt) ! computing a safe value of norbit
         if (nid == 0) write (6, *) ' Computing norbit=', norbit

         dt = porbit/norbit ! reducing dt to match forced orbit to machine accuracy
         param(12) = dt

         if (nid == 0) write (6, *) ' Adjusting timeStep dt=', dt
         call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
         if (nid == 0) write (6, *) ' veryfing current CFL and target=', ctarg, param(26)
         param(12) = -abs(param(12))

         gain = -0.04432d0*2.0d0*NEKSTAB_PI/porbit ! Theoretical optimal feedback parameter see reference

         if (nid == 0) write (6, *) 'Allocating TDF orbit with nsteps:', norbit, norbit*dt
         ! Allocation guard pattern (allocate once at init for requested norbit size;
         ! model from energy_budget.f90 ensure_* + krylov_inner_products workspace helpers;
         ! dealloc on re-init or end). WHY: orbit storage grows with user norbit; fixed after init.
         request_gb = 2.0d0*real(lv)*real(norbit)*8.0d0/1.0d9
         allocate (uor(lv, norbit), vor(lv, norbit), stat=alloc_stat)
         if (alloc_stat /= 0) then
            call nekStab_error('orbit allocation failed for uor/vor')
            call exitti('orbit allocation failed$', alloc_stat)
         end if
         if (if3d) then
            request_gb = real(lv)*real(norbit)*8.0d0/1.0d9
            allocate (wor(lv, norbit), stat=alloc_stat)
         else
            request_gb = 8.0d0/1.0d9
            allocate (wor(1, 1), stat=alloc_stat)
         end if
         if (alloc_stat /= 0) then
            call nekStab_error('orbit allocation failed for wor')
            call exitti('orbit allocation failed$', alloc_stat)
         end if
         call oprzero(uor(:, :), vor(:, :), wor(:, :))

         if (ifto .or. ldimt > 1) then
            request_gb = real(lt)*real(norbit)*real(ldimt)*8.0d0/1.0d9
            allocate (tor(lt, norbit, ldimt), stat=alloc_stat)
            if (alloc_stat /= 0) then
               call nekStab_error('orbit allocation failed for tor')
               call exitti('orbit allocation failed$', alloc_stat)
            end if
            tor(:, :, :) = 0.0d0
         end if

         rate = 0.0d0; residu0 = 0.0d0; ibuf = 0
         open (unit=NEKSTAB_UNIT_RESIDU, file='residu.dat')

         tdf_init = .true.

      else

         if (istep <= norbit) then !t<T->save solutions

            if (nid == 0) write (6, *) ' Storing initial solution in memory:', istep, '/', norbit
            call opcopy(uor(:, istep), vor(:, istep), wor(:, istep), vx, vy, vz)
            if (ifto) call copy(tor(1, istep, 1), t(:, :, :, :, 1), nt)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(tor(1, istep, m), t(:, :, :, :, m), nt)
               end do
            end if

         else !t>T->compute forcing !f(t)= - \Lambda * 2*pi*St * ( u(t) - u(t-T) )

            !  Circular buffer: ibuf cycles 1..norbit so old snapshots are overwritten
            !  in-place.  This replaces an O(norbit) array shift that was the original
            !  implementation.  ibuf advances BEFORE reading so that uor(:,ibuf) holds
            !  the snapshot from exactly T seconds ago (the oldest in the ring).
            ibuf = mod(ibuf, norbit) + 1
            if (.not. allocated(do1)) allocate (do1(lv), do2(lv), do3(lv))
            !  TDF feedback on the FULL state  f = gain*(q(t) - q(t-T)):
            !  build the current-minus-orbit difference for velocity AND every
            !  stored scalar (temperature + passive scalars), take the
            !  multi-field residual, and force all fields via the shared helper.
            !  k_zero leaves non-stored fields at 0 -> helper adds nothing there.
            call opsub3(do1, do2, do3, vx, vy, vz, uor(:, ibuf), vor(:, ibuf), wor(:, ibuf))
            call k_zero(tdfD)
            call opcopy(tdfD%vx, tdfD%vy, tdfD%vz, do1, do2, do3)
            if (ifto) call sub3(tdfD%t(1, 1), t(1, 1, 1, 1, 1), tor(1, ibuf, 1), nt)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call sub3(tdfD%t(1, m), t(1, 1, 1, 1, m), tor(1, ibuf, m), nt)
               end do
            end if
            call k_norm(l2, tdfD)
            rate = (l2 - residu0)
            residu0 = l2
            call apply_feedback_forcing(tdfD, gain)

            call opcopy(uor(1, ibuf), vor(1, ibuf), wor(1, ibuf), vx, vy, vz)
            if (ifto) call copy(tor(1, ibuf, 1), t(:, :, :, :, 1), nt)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(tor(1, ibuf, m), t(:, :, :, :, m), nt)
               end do
            end if

            if (nid == 0) then
               write (NEKSTAB_UNIT_RESIDU, "(3E15.7)") time, l2, rate
               write (6, "(' TDF residu=',1pE11.4,' rate of change= ',1pE11.4)") l2, rate
               write (6, *) ' '
            end if

            tol = max(param(21), param(22))
            if (l2 > 0.0d0 .and. l2 < tol) then
               if (nid == 0) write (6, *) ' Converged base flow to:', tol
               call close_residu_unit()
               ifbfcv = .true.
               call bcast(ifbfcv, lsize)
               param(63) = 1.0d0 ! Enforce 64-bit output
               call bcast(param(63), wdsize)
               call outpost2(vx, vy, vz, pr, t, nof, 'BF_')
               param(63) = 0.0d0 ! Reset to 32-bit output
               call bcast(param(63), wdsize)
               call outpost_vort(vx, vy, vz, 'BFV')
            end if

         end if ! else
         if (istep == nsteps) call close_residu_unit()
      end if ! not tdf_init

   end subroutine tdf

   !-----------------------------------------------------------------------
   ! SFD — Selective Frequency Damping
   !
   ! Purpose:
   !   Damps unsteady oscillations to converge to a steady base flow
   !   using either the Akervik or Casacuberta formulation.
   !-----------------------------------------------------------------------
   subroutine SFD
      use nekstab_krylov_subspace, only: krylov_vector, k_zero, k_copy, k_cmult, k_add2, &
                                         k_sub3, k_norm, NEKSTAB_PI
      type(krylov_vector), save :: oldQ, oldV, qa, qb, qc
      type(krylov_vector) :: tempD, tempM
      real adt, bdt, cdt, cutoff, gain, res, h1, l2, semi, linf, frq, sig, rate
      real :: current_solver_tol, tol_delta, tol_scale

      frq = abs(uparam(04))*2.0d0*NEKSTAB_PI ! St to omega
      sig = abs(uparam(05))

      if (uparam(5) > 0) then

         if (uparam(4) > 0) then ! Akervik
            cutoff = 0.5d0*frq
            gain = -2.0d0*sig
         else ! Casacuberta
            cutoff = 0.5d0*(sqrt(frq**2 + sig**2) - sig)
            gain = -0.5d0*(sqrt(frq**2 + sig**2) + sig)
         end if

         if (istep == 0) then
            if (.not. sfd_init) then
               call configure_sfd_dynamic_tolerance()
               if (nid == 0) open (unit=NEKSTAB_UNIT_RESIDU, file='residu.dat')
               if (nid == 0 .and. ifdyntol) then
                  open (unit=NEKSTAB_UNIT_DYNTOL, file='dyn_tol.dat')
                  write (NEKSTAB_UNIT_DYNTOL, '(A)') '# time residual current_tol requested_tol used_tol cap'
               end if
               sfd_init = .true.
            end if
            SFD_dtol = max(param(21), param(22)); SFD_oldRes = 0.0d0
            call k_zero(tempM)
            call k_zero(qa)
            call k_zero(qb)
            call k_zero(qc)
            call nopcopy(oldQ%vx, oldQ%vy, oldQ%vz, oldQ%pr, oldQ%t, vx, vy, vz, pr, t)
            call nopcopy(oldV%vx, oldV%vy, oldV%vz, oldV%pr, oldV%t, vx, vy, vz, pr, t)
         else
            call setab3(adt, bdt, cdt)
            call k_copy(qc, qb)
            call k_copy(qb, qa)
            call k_sub3(qa, oldV, oldQ)
            call k_copy(tempD, qa)
            call k_cmult(tempD, adt)
            call k_copy(tempM, tempD)
            call k_copy(tempD, qb)
            call k_cmult(tempD, bdt)
            call k_add2(tempM, tempD)
            call k_copy(tempD, qc)
            call k_cmult(tempD, cdt)
            call k_add2(tempM, tempD)
            call k_cmult(tempM, cutoff*dt)
            call k_add2(oldQ, tempM)
         end if

         !  SFD feedback forcing on the FULL state  f = gain*(q - q_filtered):
         !  velocity -> fcx/fcy/fcz and every scalar (temperature + passive
         !  scalars such as the RANS tke/tau) -> fct, via apply_feedback_forcing.
         !  A RANS fixed point needs the scalars damped too; otherwise tke/tau
         !  keep drifting while only the velocity has converged and the result
         !  is not a true fixed point of the full operator.
         call nopcopy(tempM%vx, tempM%vy, tempM%vz, tempM%pr, tempM%t, vx, vy, vz, pr, t)
         call k_sub3(tempD, tempM, oldQ)        ! tempD = current - filtered state
         call apply_feedback_forcing(tempD, gain)

      elseif (.not. sfd_init) then

         call configure_sfd_dynamic_tolerance()
         if (nid == 0) then
            open (unit=NEKSTAB_UNIT_RESIDU, file='residu.dat')
            write (6, *) ' SFD in continuation mode'
            if (ifdyntol) then
               open (unit=NEKSTAB_UNIT_DYNTOL, file='dyn_tol.dat')
               write (NEKSTAB_UNIT_DYNTOL, '(A)') '# time residual current_tol requested_tol used_tol cap'
            end if
         end if
         sfd_init = .true.

      end if

      if (istep >= 1) then
         !  Convergence residual on the FULL state (multi-field, thermal-weighted
         !  k_norm): SFD only declares convergence once velocity AND every scalar
         !  have stopped changing.  Velocity-only normvc would call a still-
         !  drifting RANS scalar field "converged".
         call nopcopy(tempM%vx, tempM%vy, tempM%vz, tempM%pr, tempM%t, vx, vy, vz, pr, t)
         call k_sub3(tempD, tempM, oldV)        ! tempD = current - previous step
         call k_norm(res, tempD)
         rate = (res - SFD_oldRes)/dt; SFD_oldRes = res
         call k_copy(oldV, tempM)               ! oldV = current state (next step)
         if (nid == 0) then
            write (NEKSTAB_UNIT_RESIDU, "(4E15.7)") time, res, rate, param(21)
            write (6, "(A,3E15.7)") '  SFD residual, rate:', res, rate
            if (uparam(4) > 0) then
               write (6, *) ' Akervik  cutoff, gain:', cutoff, gain
            else
               write (6, *) ' Casacub. cutoff, gain:', cutoff, gain
            end if
         end if
         if (ifdyntol .and. mod(istep, sfd_tol_update_stride) == 0 &
             .and. res > 0.0d0) then
            ! Dynamic SFD scheduling uses the current fixed-point residual as a
            ! cheap proxy for how tightly the inner Nek5000 solves need to be
            ! converged.  Keep the target just below the current residual, while
            ! never tightening below the requested fixed-point convergence gate.
            ! The legacy ew_tol_cap variable remains an optional upper bound for
            ! stiff cases.
            sfd_solver_tol = max(sfd_tol_residual_factor*res, SFD_dtol)
            if (ew_tol_cap > 0.0d0) sfd_solver_tol = min(sfd_solver_tol, ew_tol_cap)
            if (nid == 0) then
               write (NEKSTAB_UNIT_DYNTOL, '(6E15.7)') time, res, param(22), &
                  sfd_tol_residual_factor*res, &
                  sfd_solver_tol, ew_tol_cap
            end if
            current_solver_tol = param(22)
            tol_delta = abs(sfd_solver_tol - current_solver_tol)
            tol_scale = max(abs(current_solver_tol), abs(sfd_solver_tol), 1.0d0)
            if (tol_delta > sfd_tol_update_rtol*tol_scale) then
               call set_nek5000_velocity_tolerance(sfd_solver_tol)
            end if
         end if

         if (istep > 100 .and. res < SFD_dtol) then
            if (nid == 0) write (6, *) ' Converged base flow to:', res
            call close_residu_unit()
            call close_dyn_tol_unit()
            ifbfcv = .true.
            call bcast(ifbfcv, lsize)
            param(63) = 1.0d0 ! Enforce 64-bit output
            call bcast(param(63), wdsize)
            call outpost2(vx, vy, vz, pr, t, nof, 'BF_')
            param(63) = 0.0d0 ! Reset to 32-bit output
            call bcast(param(63), wdsize)
            call outpost_vort(vx, vy, vz, 'BFV')
         end if

         if (istep == nsteps) then
            call close_residu_unit()
            call close_dyn_tol_unit()
         end if

      end if

   end subroutine SFD

   !-----------------------------------------------------------------------
   ! BoostConv — BoostConv convergence acceleration
   !
   ! Purpose:
   !   Accelerates convergence to steady state using Anderson-like
   !   extrapolation from residual snapshots.
   !-----------------------------------------------------------------------
   subroutine BoostConv
      use nekstab_krylov_subspace, only: lv

      real, allocatable, save, dimension(:) :: dvx, dvy, dvz
      real :: residu, h1, semi, linf, rate, tol
      real, save :: residu0
      if (mod(istep, bst_skp) == 0) then

         if (.not. boostconv_init) then
            allocate (dvx(lv), dvy(lv), dvz(lv))
            residu = 0.0d0; rate = 0.0d0; residu0 = 0.0d0
            open (unit=NEKSTAB_UNIT_RESIDU, file='residu.dat')
            !  BoostConv accelerates the VELOCITY update only: its QR residual
            !  subspace (q_x/q_y/q_z) carries no scalar components, so for a
            !  multi-field state (RANS tke/tau, thermal) the scalars are NOT
            !  accelerated and the converged result is not a simultaneous fixed
            !  point of the full operator.  Warn loudly rather than silently
            !  return a velocity-only base flow; use SFD for multi-field cases.
            if (nfield > 1 .and. nid == 0) then
               write (6, *) 'nekStab WARNING: BoostConv is velocity-only;'
               write (6, *) '  scalar fields (temperature / tke / tau) are NOT'
               write (6, *) '  accelerated. Use SFD for multi-field base flows.'
            end if
            boostconv_init = .true.
         end if

         call opsub3(dvx, dvy, dvz, vx, vy, vz, vxlag(1, 1, 1, 1, 1), vylag(1, 1, 1, 1, 1), vzlag(1, 1, 1, 1, 1)) !dv=v-vold
         call normvc(h1, semi, residu, linf, dvx, dvy, dvz); rate = (residu - residu0)/dt; residu0 = residu
         call boostconv_core(dvx, dvy, dvz)
         call opadd3(vx, vy, vz, vxlag(1, 1, 1, 1, 1), vylag(1, 1, 1, 1, 1), vzlag(1, 1, 1, 1, 1), dvx, dvy, dvz) !v=vold+dv

         if (nid == 0) then
            write (NEKSTAB_UNIT_RESIDU, "(3E15.7)") time, residu, rate
            write (6, "(' BoostConv residu=',1pE11.4,' delta= ',1pE11.4)") residu, rate
            write (6, *) ' '
         end if

         tol = max(param(21), param(22))
         if (residu < tol) then
            if (nid == 0) write (6, *) ' Converged base flow to:', tol
            call close_residu_unit()
            ifbfcv = .true.
            call bcast(ifbfcv, lsize)
            param(63) = 1.0d0 ! Enforce 64-bit output
            call bcast(param(63), wdsize)
            call outpost2(vx, vy, vz, pr, t, nof, 'BF_')
            param(63) = 0.0d0 ! Reset to 32-bit output
            call bcast(param(63), wdsize)
         end if

      end if
      if (istep == nsteps) call close_residu_unit()

   end subroutine BoostConv

   !-----------------------------------------------------------------------
   ! boostconv_core — Inner iteration of BoostConv acceleration
   !
   ! Purpose:
   !   Performs QR-based least-squares extrapolation on residual
   !   snapshots to accelerate convergence.
   !
   ! Arguments:
   !   rbx, rby, rbz [inout] -- velocity residual, updated in-place
   !-----------------------------------------------------------------------
   subroutine boostconv_core(rbx, rby, rbz)
      use nekstab_krylov_subspace, only: lv, nv

      integer, save :: rot
      integer :: j

      real, allocatable, save, dimension(:) :: cc, ccb
      real, allocatable, save, dimension(:, :) :: dd
      real, allocatable, save, dimension(:, :) :: q_x, q_y, q_z
      real, allocatable, save, dimension(:, :) :: x_x, x_y, x_z, y_x, y_y, y_z

      real, dimension(lv), intent(inout) :: rbx, rby, rbz
      real, allocatable, save, dimension(:) :: dumx, dumy, dumz
      real, allocatable, save, dimension(:) :: fwrk
      real, dimension(bst_snp) :: wk_gop

      real :: glsc3

      if (.not. boostconv_core_init) then

         allocate (cc(bst_snp), ccb(bst_snp), dd(bst_snp, bst_snp))
         allocate (q_x(lv, bst_snp), q_y(lv, bst_snp), q_z(lv, bst_snp))
         allocate (x_x(lv, bst_snp), x_y(lv, bst_snp), x_z(lv, bst_snp))
         allocate (y_x(lv, bst_snp), y_y(lv, bst_snp), y_z(lv, bst_snp))
         allocate (dumx(lv), dumy(lv), dumz(lv), fwrk(lv))

         if (nid == 0) write (6, *) 'Allocating BoostConv variables with:', bst_snp
         if (nid == 0) write (6, *) '                     skipping every:', bst_skp

         call oprzero(x_x(:, :), x_y(:, :), x_z(:, :))
         call oprzero(y_x(:, :), y_y(:, :), y_z(:, :))
         call oprzero(q_x(:, :), q_y(:, :), q_z(:, :))
         call opcopy(y_x(:, 1), y_y(:, 1), y_z(:, 1), rbx, rby, rbz)
         call opcopy(x_x(:, 1), x_y(:, 1), x_z(:, 1), rbx, rby, rbz)
         dd(:, :) = 1.0d0; rot = 1; boostconv_core_init = .true.

      else

         call opsub2(y_x(:, rot), y_y(:, rot), y_z(:, rot), rbx, rby, rbz)
         call opsub2(x_x(:, rot), x_y(:, rot), x_z(:, rot), y_x(:, rot), y_y(:, rot), y_z(:, rot))
         call qr_dec(dd, q_x, q_y, q_z, y_x, y_y, y_z)

!           Batch projection: cc = Q^T * (bm1 * rb) via dgemv.
!           LDA = lv (compile-time leading dimension of q_x/q_y/q_z),
!           M = nv (runtime active rows). BLAS requires LDA >= M and
!           the stride must match the array layout — using nv as LDA
!           would be wrong since q_x is allocated as (lv, bst_snp).
         call copy(fwrk, rbx, nv)
         call col2(fwrk, bm1, nv)
         call dgemv('T', nv, bst_snp, 1.0d0, &
                    q_x, lv, fwrk, 1, 0.0d0, cc, 1)
         call copy(fwrk, rby, nv)
         call col2(fwrk, bm1, nv)
         call dgemv('T', nv, bst_snp, 1.0d0, &
                    q_y, lv, fwrk, 1, 1.0d0, cc, 1)
         if (if3d) then
            call copy(fwrk, rbz, nv)
            call col2(fwrk, bm1, nv)
            call dgemv('T', nv, bst_snp, 1.0d0, &
                       q_z, lv, fwrk, 1, 1.0d0, cc, 1)
         end if
         call gop(cc, wk_gop, '+  ', bst_snp)

         call linear_system(ccb, cc, dd, bst_snp); rot = mod(rot, bst_snp) + 1
         call opcopy(y_x(:, rot), y_y(:, rot), y_z(:, rot), rbx, rby, rbz)

         do j = 1, bst_snp
            call opcopy(dumx, dumy, dumz, x_x(:, j), x_y(:, j), x_z(:, j))
            call opcmult(dumx, dumy, dumz, ccb(j))
            call opadd2(rbx, rby, rbz, dumx, dumy, dumz)
         end do
         call opcopy(x_x(:, rot), x_y(:, rot), x_z(:, rot), rbx, rby, rbz)

      end if

   end subroutine boostconv_core

   !-----------------------------------------------------------------------
   ! qr_dec — QR decomposition for BoostConv residual subspace
   !
   ! Purpose:
   !   Computes modified Gram-Schmidt QR factorization of the
   !   residual difference subspace used by BoostConv.
   !
   ! Arguments:
   !   rr                    [out]   -- upper triangular R factor
   !   q_x, q_y, q_z        [out]   -- orthonormal Q basis
   !   x_x, x_y, x_z        [in]    -- input residual vectors
   !-----------------------------------------------------------------------
   subroutine qr_dec(rr, q_x, q_y, q_z, x_x, x_y, x_z)
      use nekstab_krylov_subspace, only: lv, nv
      integer i, j
      real, dimension(lv, bst_snp), intent(in) :: x_x, x_y, x_z
      real, dimension(lv, bst_snp), intent(out) :: q_x, q_y, q_z
      real, dimension(lv) :: dum_x, dum_y, dum_z
      real, dimension(bst_snp, bst_snp), intent(out) :: rr
      real norma, norma_loc, glsc3, proj

      ! nv from globals
      rr = 0.0d0; norma = 0.0d0
      call oprzero(Q_x(:, :), Q_y(:, :), Q_z(:, :))

      call opcopy(dum_x, dum_y, dum_z, &
                  x_x(:, 1), x_y(:, 1), x_z(:, 1))

!        First column norm (scalar, no batch needed)
      norma = glsc3(dum_x, bm1, dum_x, nv) &
              + glsc3(dum_y, bm1, dum_y, nv)
      if (if3d) norma = norma &
                        + glsc3(dum_z, bm1, dum_z, nv)
      norma = sqrt(norma)
      if (norma < 1.0d-60) then
         if (nid == 0) write (6, *) &
            'WARNING: qr_dec near-zero first column norm'
         return
      end if

      ! norma >= 1e-60 guaranteed by the early return above.
      call opcmult(dum_x, dum_y, dum_z, 1.0d0/norma)
      call opcopy(q_x(:, 1), q_y(:, 1), q_z(:, 1), &
                  dum_x, dum_y, dum_z)
      rr(1, 1) = norma

      do j = 2, bst_snp
         call opcopy(dum_x, dum_y, dum_z, &
                     x_x(:, j), x_y(:, j), x_z(:, j))

!           Modified Gram-Schmidt: project onto and subtract each q one at a
!           time, vs the old classical GS that projected onto all previous q
!           at once. Sequential subtraction keeps Q orthonormal even when the
!           residual-difference subspace becomes near-linearly-dependent close
!           to convergence; classical GS loses orthogonality there (error grows
!           O(eps*kappa^2) vs O(eps*kappa)), which corrupted R and degraded the
!           BoostConv correction to noise -> stall at the limit-cycle amplitude.
!           glsc3 carries its own global reduction, so no separate gop is needed.
         do i = 1, j - 1
            proj = glsc3(dum_x, bm1, q_x(:, i), nv) &
                   + glsc3(dum_y, bm1, q_y(:, i), nv)
            if (if3d) proj = proj + glsc3(dum_z, bm1, q_z(:, i), nv)
            rr(i, j) = proj
            call add2s2(dum_x, q_x(:, i), -proj, nv)
            call add2s2(dum_y, q_y(:, i), -proj, nv)
            if (if3d) call add2s2(dum_z, q_z(:, i), -proj, nv)
         end do

!           Column norm
         norma = glsc3(dum_x, bm1, dum_x, nv) &
                 + glsc3(dum_y, bm1, dum_y, nv)
         if (if3d) norma = norma &
                           + glsc3(dum_z, bm1, dum_z, nv)

         if (norma < 1.0d-60) then
            norma = 1.0d0
            q_x(:, j) = 0.0d0; q_y(:, j) = 0.0d0
            if (if3d) q_z(:, j) = 0.0d0
         else
            call opcmult(dum_x, dum_y, dum_z, &
                         1.0d0/sqrt(norma))
            call opcopy(q_x(:, j), q_y(:, j), q_z(:, j), &
                        dum_x, dum_y, dum_z)
         end if

         rr(j, j) = sqrt(norma)

      end do

   end subroutine qr_dec

   !-----------------------------------------------------------------------
   ! linear_system — Back-substitution for upper triangular system
   !
   ! Arguments:
   !   outp   [out] -- solution vector
   !   inp    [in]  -- right-hand side vector
   !   m      [in]  -- upper triangular matrix
   !   size_m [in]  -- system dimension
   !-----------------------------------------------------------------------
   subroutine linear_system(outp, inp, m, size_m)
      integer, intent(in) :: size_m
      integer :: j, k
      real, intent(in) :: m(size_m, size_m), inp(size_m)
      real, intent(out) :: outp(size_m)
      outp = 0.0d0
      do j = size_m, 1, -1
         outp(j) = inp(j)
         do k = j + 1, size_m
            outp(j) = outp(j) - m(j, k)*outp(k)
         end do
         if (abs(m(j, j)) < 1.0d-60) then
            if (nid == 0) write (6, *) 'WARNING: linear_system singular diagonal at j=', j
            outp(j) = 0.0d0
         else
            outp(j) = outp(j)/m(j, j)
         end if
      end do

   end subroutine linear_system
   !-----------------------------------------------------------------------

end module nekstab_fixedpoint
