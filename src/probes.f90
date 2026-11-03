!-----------------------------------------------------------------------
! probes.f90 — Point probes and zero-crossing detection
!
! Purpose:
!   Locates grid points nearest to user-specified coordinates and
!   monitors velocity at those points for zero-crossing period
!   detection during DNS.
!
! Public interface:
!   pointcheck    — find nearest grid point to probe location
!   zero_crossing — monitor velocity zero-crossings for period
!
! Dependencies:
!   krylov_subspace, SIZE, TOTAL
!-----------------------------------------------------------------------

module nekstab_probes
   use nekstab_nek_bridge, only: nekStab_log, nid, lx1, ly1, lz1, lelt, lelv, &
                                 xck, yck, zck, nx1, ny1, nz1, nelt, if3D, &
                                 xm1, ym1, zm1, istep, time, dt, vx, vy, vz, &
                                 NEKSTAB_UNIT_ZC1, NEKSTAB_UNIT_ZC2, &
                                 zc_component
   use nekstab_zero_cross, only: crossing_reset, crossing_state, crossing_update
   implicit none
   private
   public :: pointcheck, zero_crossing
contains

!-----------------------------------------------------------------------
! pointcheck — Locate nearest grid point to probe coordinates
!
! Purpose:
!   Finds the grid point closest to (xck, yck, zck) and returns
!   its local index and owning processor flag.
!
! Arguments:
!   posiz   [out] — local index of closest grid point
!   procmin [out] — 1 if this rank owns the point, 0 otherwise
!-----------------------------------------------------------------------
   subroutine pointcheck(posiz, procmin)
      integer, intent(out) :: posiz, procmin
      integer :: m, n, owner_nid
      real :: chk(lx1*ly1*lz1*lelt), chkmin, glmin, cand_nid

      n = nx1*ny1*nz1*nelt
      procmin = 0
      posiz = 0
      call nekStab_log('Evaluating probe at:')
      if (nid == 0) write (6, *) xck, yck, zck
      do m = 1, n
         chk(m) = (xm1(m, 1, 1, 1) - xck)**2 + (ym1(m, 1, 1, 1) - yck)**2
         if (if3D) chk(m) = chk(m) + (zm1(m, 1, 1, 1) - zck)**2
      end do
      chkmin = glmin(chk, n)
      do m = 1, n
         if (abs(chk(m) - chkmin) <= 1.0d-12*max(1.0d0, chkmin)) then
            posiz = m
            exit
         end if
      end do

      ! Fix defect 3: pointcheck marked all ranks at the minimum distance as owner;
      ! on element boundaries shared by ranks, glsum doubled the velocity.
      ! Now only the lowest nid among ranks at chkmin is the owner.
      if (posiz > 0) then
         cand_nid = real(nid, 8)
      else
         cand_nid = 1.0d9
      end if
      owner_nid = nint(glmin(cand_nid, 1))

      if (nid == owner_nid) then
         procmin = 1
         call nekStab_log('Point found:')
         write (6, *) posiz ! owner rank prints the index (not necessarily nid==0)
      else
         procmin = 0
         posiz = 0
      end if

   end subroutine pointcheck
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! zero_crossing — Monitor velocity zero-crossings for period detection
!
! Purpose:
!   Tracks a probe velocity signal and detects upward zero-crossings
!   relative to a running mean. Outputs detected period and Poincare
!   section data to files.
!
! Arguments:
!   v_mean_init [in] — initial estimate of mean velocity at probe
!-----------------------------------------------------------------------
   subroutine zero_crossing(v_mean_init)
      use nekstab_krylov_subspace, only: lv
      real, intent(in), optional :: v_mean_init
      real, save :: T_delayed(lv, 3), do1(lv), do2(lv), do3(lv)
      integer, parameter :: plor = 4
      real :: h1, l2, semi, linf, l2_u, rel_norm
      save l2
      real :: glsum, vdot, vddot
      real, save :: velp(plor)
      type(crossing_state), save :: zc
      integer, save :: probe_nel, probe_nid
      integer :: i, phase
      real :: vel_now, v_mean_start
      real :: theta, t_hit, period, v_logged, period_change

      if (present(v_mean_init)) then
         v_mean_start = v_mean_init
      else
         v_mean_start = 0.0d0
      end if

      if (istep == 0) then
         call nekStab_log('Initializing zero-crossing routine...')
         probe_nel = 0; probe_nid = 0; vdot = 0.0d0; vddot = 0.0d0
         l2 = 0.0d0
         call crossing_reset(zc, v_mean_start)
         call pointcheck(probe_nel, probe_nid) !alter xck, yck, zck in usrchk
         if (nid == 0) open (unit=NEKSTAB_UNIT_ZC1, file='zc_period.dat')
         if (nid == 0) open (unit=NEKSTAB_UNIT_ZC2, file='zc_poincare.dat')
         call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
      end if

      vel_now = 0.0d0
      if (probe_nid == 1) then
         ! Defect 5 / component selection: 1=vx, 2=vy (default), 3=vz
         if (zc_component == 1) then
            vel_now = vx(probe_nel, 1, 1, 1)
         elseif (zc_component == 3) then
            vel_now = vz(probe_nel, 1, 1, 1)
         else
            vel_now = vy(probe_nel, 1, 1, 1)
         end if
      end if
      vel_now = glsum(vel_now, 1) ! broadcast from the single owner rank

      if (istep == 0) then
         velp(:) = vel_now
      else
         velp(plor) = vel_now
      end if

      ! Period, crossing time and the mean over the last period.
      call crossing_update(zc, velp(plor - 1), velp(plor), time, dt, istep > 0, &
                           phase, theta, t_hit, period, v_logged, period_change)

      if (phase == 1) then
         call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
         if (nid == 0) write (6, '(A,ES15.7)') ' Zero-crossing: first crossing at time = ', t_hit
      else if (phase == 2) then
         ! Relative mass-weighted norm ||u(t)-u(t-T)|| / ||u(t)||
         call opsub3(do1, do2, do3, vx, vy, vz, T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3))
         call normvc(h1, semi, l2, linf, do1, do2, do3)
         call normvc(h1, semi, l2_u, linf, vx, vy, vz)
         if (l2_u > 0.0d0) then
            rel_norm = l2/l2_u
         else
            rel_norm = 0.0d0
         end if
         call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
         if (nid == 0) then
            write (6, '(A,ES14.7,A,ES14.7,A,ES14.7,A,ES14.7)') &
               ' Zero-crossing time = ', t_hit, '  T = ', period, &
               '  |T-T_old| = ', period_change, '  rel_norm = ', rel_norm
            write (NEKSTAB_UNIT_ZC1, "(5E15.7)") t_hit, period, period_change, v_logged, rel_norm
            flush (6)
            flush (NEKSTAB_UNIT_ZC1)
         end if
      else if (phase == 3) then
         ! The period restarts at a later crossing of the new level: keep the field of that point.
         call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
      end if
!     https://en.wikipedia.org/wiki/Finite_difference_coefficient
      if (istep > 3) then
         vdot = ((11.0d0/6.0d0)*velp(plor) - 3.0d0*velp(plor - 1) + 1.5d0*velp(plor - 2) - (1.0d0/3.0d0)*velp(plor - 3))*dt**(-1)
         vddot = (2.0d0*velp(plor) - 5.0d0*velp(plor - 1) + 4.0d0*velp(plor - 2) - 1.0d0*velp(plor - 3))*dt**(-2)
      else
         vdot = 0.0d0; vddot = 0.0d0
      end if
      if (nid == 0) write (NEKSTAB_UNIT_ZC2, "(4E15.7)") time, velp(plor), vdot, vddot

      do i = 1, plor - 1
         velp(i) = velp(i + 1)
      end do

   end subroutine zero_crossing
!-----------------------------------------------------------------------

end module nekstab_probes
