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
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer, intent(out) :: posiz, procmin
         integer :: m, n
         real chk(lx1*ly1*lz1*lelt), chkmin, glmin

         n = nx1*ny1*nz1*nelt
         procmin = 0
         if (nid == 0) write (6, *) 'Evaluating probe at:', xck, yck, zck
         do m = 1, n
            chk(m) = (xm1(m, 1, 1, 1) - xck)**2 + (ym1(m, 1, 1, 1) - yck)**2
            if (if3D) chk(m) = chk(m) + (zm1(m, 1, 1, 1) - zck)**2
         end do
         chkmin = glmin(chk, n)
         do m = 1, n
            if (chkmin == chk(m)) then
               procmin = 1
               posiz = m
               write(6,*) 'Point found:', m ! do not use nid = 0 as it could be any rank with this value !
            end if
         end do

         return
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
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, intent(in) :: v_mean_init
         real, save :: T_delayed(lv, 3), do1(lv), do2(lv), do3(lv), v_mean
         integer, parameter :: plor = 2
         real h1, l2, semi, linf
         save l2
         real :: glsum, dtime, vdot, vddot
         real, save :: velp(plor), v_sum, time0
         real, save :: t_cross, t_cross_old
         real, save :: v_cross, v_cross_old
         real, save :: p_now, p_sum, p_mean, p_old
         integer, save :: probe_nel, probe_nid, t_cross_count
         integer :: i

         if (istep == 0) then
            if (nid == 0) write (6, *) 'Initializing zero-crossing routine...'
            probe_nel = 0; probe_nid = 0; vdot = 0.0d0; vddot = 0.0d0
            velp(:) = 0.0d0; v_sum = 0.0d0; v_mean = v_mean_init
            p_now = 0.0d0; p_sum = 0.0d0; p_mean = 0.0d0
            t_cross = 0.0d0; t_cross_old = 0.0d0; t_cross_count = 0
            time0 = time; p_old = 0.0d0; l2 = 0.0d0
            call pointcheck(probe_nel, probe_nid) !alter xck, yck, zck in usrchck
            if (nid == 0) open (unit=17, file='zc_period.dat')
            if (nid == 0) open (unit=19, file='zc_poincare.dat')
            call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
         end if

         velp(plor) = 0.0d0
         if (probe_nid == 1) velp(plor) = vy(probe_nel, 1, 1, 1)
         velp(plor) = glsum(velp(plor), 1) !broadcast
         dtime = time - time0

         if (istep > 1) then
            v_sum = v_sum + velp(plor)*dt !(v_mean*(dtime-dt)+v_now*dt)/dtime
            v_mean = v_mean_init + v_sum/dtime
         end if
         if (velp(plor - 1) <= v_mean .and. velp(plor) >= v_mean) then !period found

            p_old = p_now !save old value
            t_cross_old = t_cross !save old value
            t_cross = dtime !update new value
            p_now = t_cross - t_cross_old !compute period

            call opsub3(do1, do2, do3, vx, vy, vz, T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3)) !ub=v-vold
            call normvc(h1, semi, l2, linf, do1, do2, do3)
            call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
            if (nid == 0) write (6, *) ' Zero-crossing T=', p_now, abs(p_now - p_old), l2
            v_cross_old = v_cross; v_cross = velp(plor)
            if (nid == 0) write (17, "(5E15.7)") time, p_now, abs(p_now - p_old), v_mean, l2

         end if
      !     https://en.wikipedia.org/wiki/Finite_difference_coefficient
         if (istep > 3) then
            vdot = ((11./6.)*velp(plor) - 3*velp(plor - 1) + 1.5*velp(plor - 2) - (1./3.)*velp(plor - 3))*dt**(-1)
            vddot = (2*velp(plor) - 5*velp(plor - 1) + 4.0*velp(plor - 2) - 1.0*velp(plor - 3))*dt**(-2)
         else
            vdot = 0.0d0; vddot = 0.0d0
         end if
         if (nid == 0) write (19, "(4E15.7)") time, velp(plor), vdot, vddot

         do i = 1, plor - 1
            velp(i) = velp(i + 1)
         end do

         return
      end subroutine zero_crossing
      !-----------------------------------------------------------------------

      end module nekstab_probes
