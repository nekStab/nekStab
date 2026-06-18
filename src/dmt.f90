!-----------------------------------------------------------------------
! dmt.f90 — Dynamic Mode Tracking (DMT) for base flow stabilization
!
! Purpose:
!   Implements the Dynamic Mode Tracking technique of Queguineur et al.
!   (Phys. Fluids 31(3), 034101, 2019) for stabilizing unstable periodic
!   orbits by band-pass filtering the velocity field.
!
! Algorithm:
!   PDE:   d2ubar/dt2 - d2u/dt2 + beta*dubar/dt = omc2*(u - ubar)
!   Force: f = -chi*(u - ubar)
!
!   Discrete form uses a BDF2 stencil for d2u/dt2 (3 history levels)
!   and a forward Euler step for the coupled (Y, ubar) system.
!
!   Change of variables: Y = dubar/dt (auxiliary variable)
!   BDF2:  d2u/dt2 ~ (2u - 5u(-1) + 4u(-2) - u(-3)) / dt2
!   Y:     Y = Y(-1) + dt*[(ubar(-1)-u)*(-omc2) + Y(-1)*(-beta) + d2u/dt2]
!   ubar:  ubar = ubar(-1) + dt*Y(-1)
!
! Parameter mapping (uparam indices):
!   uparam(11) = omc (target angular frequency; legacy used dmt_fki=freq, omc2=(2*pi*f)^2)
!   uparam(12) = beta (filter width / damping)
!   uparam(13) = chi (feedback gain)
!   uparam(14) = t_start (time before forcing kicks in)
!   uparam(15) = tol (L2 residual convergence tolerance)
!
! Note on omc2: legacy computed omc2=(2*pi*f)^2; new module takes omega_c directly,
!   so omc2 = uparam(11)**2.
!
! Public interface:
!   dmt      — per-step driver; call from main.f90 fixed-point branch
!   dmt_init — allocator + first-call setup; also called lazily from dmt
!
! Dependencies:
!   nekstab_nek_bridge (wraps SIZE, TOTAL, ADJOINT)
!
! Port history:
!   Legacy source: nekStab_old_bkp/core_dev_floquet_old/x_periodic_orbit.f:266-393
!   Ported to free-form F90 module by Ricardo Frantz, 2026.
!-----------------------------------------------------------------------

module nekstab_dmt
   use nekstab_nek_bridge, only: uparam, lx1, ly1, lz1, lelt, nelt, nid, &
                                 nekStab_error, nekStab_log, &
                                 NEKSTAB_UNIT_DYNTOL, istep, time, dt, &
                                 nsteps, vx, vy, vz, pr, t, fcx, fcy, fcz, &
                                 ifbfcv, lastep
   implicit none
   private

   logical, save :: dmt_initialized = .false.
   logical, save :: dmt_log_open = .false.
   real, allocatable, save :: vxo(:, :), vyo(:, :), vzo(:, :)
   real, allocatable, save :: xxo(:), xyo(:), xzo(:)
   real, allocatable, save :: yxo(:), yyo(:), yzo(:)
   real, allocatable, save :: do1(:), do2(:), do3(:)
   real, allocatable, save :: vxdd(:), vydd(:), vzdd(:)
   real, allocatable, save :: xx(:), xy(:), xz(:)
   real, allocatable, save :: yx(:), yy(:), yz(:)
   real, save :: residu0, omc2
   integer, save :: n

   public :: dmt_init, dmt

contains

   subroutine dmt_init
      integer, parameter :: lt = lx1*ly1*lz1*lelt
      real :: dmt_omc, dmt_ban, dmt_gan, dmt_skp, dmt_tol
      logical :: file_exists
      integer :: open_stat

      if (dmt_initialized) return

      dmt_omc = uparam(11)
      dmt_ban = uparam(12)
      dmt_gan = uparam(13)
      dmt_skp = uparam(14)
      dmt_tol = uparam(15)

      allocate (vxo(lt, 3), vyo(lt, 3), vzo(lt, 3))
      allocate (xxo(lt), xyo(lt), xzo(lt))
      allocate (yxo(lt), yyo(lt), yzo(lt))
      allocate (do1(lt), do2(lt), do3(lt))
      allocate (vxdd(lt), vydd(lt), vzdd(lt))
      allocate (xx(lt), xy(lt), xz(lt))
      allocate (yx(lt), yy(lt), yz(lt))

      vxo = 0.0d0; vyo = 0.0d0; vzo = 0.0d0
      xxo = 0.0d0; xyo = 0.0d0; xzo = 0.0d0
      yxo = 0.0d0; yyo = 0.0d0; yzo = 0.0d0
      do1 = 0.0d0; do2 = 0.0d0; do3 = 0.0d0
      vxdd = 0.0d0; vydd = 0.0d0; vzdd = 0.0d0
      xx = 0.0d0; xy = 0.0d0; xz = 0.0d0
      yx = 0.0d0; yy = 0.0d0; yz = 0.0d0

      n = lx1*ly1*lz1*nelt
      residu0 = 0.0d0
      omc2 = dmt_omc**2

      if (nid == 0) then
         inquire (file='residu_dmt.dat', exist=file_exists)
         if (file_exists) then
            open (unit=NEKSTAB_UNIT_DYNTOL, file='residu_dmt.dat', status='old', &
                  position='append', iostat=open_stat)
         else
            open (unit=NEKSTAB_UNIT_DYNTOL, file='residu_dmt.dat', status='new', iostat=open_stat)
         end if
         if (open_stat /= 0) then
            call nekStab_error('unable to open residu_dmt.dat')
         else
            dmt_log_open = .true.
         end if

         call nekStab_log('DMT: Dynamic Mode Tracking (Queguineur et al. 2019)')
         write (6, *) '   uparam(11) dmt_omc (omega_c) =', dmt_omc
         write (6, *) '   uparam(12) dmt_ban (beta)    =', dmt_ban
         write (6, *) '   uparam(13) dmt_gan (chi)      =', dmt_gan
         write (6, *) '   uparam(14) dmt_skp (t_start)  =', dmt_skp
         write (6, *) '   uparam(15) dmt_tol (tol)      =', dmt_tol
      end if

      dmt_initialized = .true.

   end subroutine dmt_init

   subroutine dmt
      integer, parameter :: lt = lx1*ly1*lz1*lelt
      real :: dmt_ban, dmt_gan, dmt_skp, dmt_tol
      real :: h1, semi, l2, linf, residu, rate, dt2
      integer :: i

      if (.not. dmt_initialized) call dmt_init

      ! Skip before first time step (istep==0 pre-step call has dt=0)
      if (istep <= 0) return

      dmt_ban = uparam(12)
      dmt_gan = uparam(13)
      dmt_skp = uparam(14)
      dmt_tol = uparam(15)

      n = lx1*ly1*lz1*nelt

      if (istep == 1) then
         call opcopy(vxo(1, 3), vyo(1, 3), vzo(1, 3), vx, vy, vz)
         return
      elseif (istep == 2) then
         call opcopy(vxo(1, 2), vyo(1, 2), vzo(1, 2), vx, vy, vz)
         return
      elseif (istep == 3) then
         call opcopy(vxo(1, 1), vyo(1, 1), vzo(1, 1), vx, vy, vz)
         call opcopy(xxo, xyo, xzo, vx, vy, vz)
         return
      end if

      call opcopy(do1, do2, do3, vx, vy, vz)

      dt2 = dt*dt
      do i = 1, n
         vxdd(i) = (2.0d0*do1(i) - 5.0d0*vxo(i, 1) &
                    + 4.0d0*vxo(i, 2) - vxo(i, 3))/dt2
         vydd(i) = (2.0d0*do2(i) - 5.0d0*vyo(i, 1) &
                    + 4.0d0*vyo(i, 2) - vyo(i, 3))/dt2
         vzdd(i) = (2.0d0*do3(i) - 5.0d0*vzo(i, 1) &
                    + 4.0d0*vzo(i, 2) - vzo(i, 3))/dt2

         yx(i) = ((xxo(i) - do1(i))*(-omc2) + yxo(i)*(-dmt_ban) &
                  + vxdd(i))*dt + yxo(i)
         yy(i) = ((xyo(i) - do2(i))*(-omc2) + yyo(i)*(-dmt_ban) &
                  + vydd(i))*dt + yyo(i)
         yz(i) = ((xzo(i) - do3(i))*(-omc2) + yzo(i)*(-dmt_ban) &
                  + vzdd(i))*dt + yzo(i)

         xx(i) = xxo(i) + yxo(i)*dt
         xy(i) = xyo(i) + yyo(i)*dt
         xz(i) = xzo(i) + yzo(i)*dt
      end do

      call opcopy(xxo, xyo, xzo, xx, xy, xz)
      call opcopy(yxo, yyo, yzo, yx, yy, yz)

      call opcopy(vxo(1, 3), vyo(1, 3), vzo(1, 3), &
                  vxo(1, 2), vyo(1, 2), vzo(1, 2))
      call opcopy(vxo(1, 2), vyo(1, 2), vzo(1, 2), &
                  vxo(1, 1), vyo(1, 1), vzo(1, 1))
      call opcopy(vxo(1, 1), vyo(1, 1), vzo(1, 1), vx, vy, vz)

      call opsub3(do1, do2, do3, vx, vy, vz, xx, xy, xz)
      call normvc(h1, semi, l2, linf, do1, do2, do3)
      residu = l2/dt
      rate = residu - residu0
      residu0 = residu

      call copy(t(1, 1, 1, 1, 1), do1, n)
      call copy(t(1, 1, 1, 1, 2), do2, n)

      call opcmult(do1, do2, do3, -dmt_gan)
      if (time > dmt_skp) call opadd2(fcx, fcy, fcz, do1, do2, do3)

      if (nid == 0) then
         if (dmt_log_open) write (NEKSTAB_UNIT_DYNTOL, "(3E15.7)") time, residu, rate
         write (6, "(' DMT residu=',1pE11.4,' rate of change= ',1pE11.4)") residu, rate
         write (6, *) ' '
      end if

      if (residu < dmt_tol) then
         if (nid == 0) write (6, *) ' Converged DMT base flow to:', residu
         call outpost(vx, vy, vz, pr, t, 'BF_')
         ifbfcv = .true.
         lastep = 1
      end if

      if (istep == nsteps .and. nid == 0 .and. dmt_log_open) then
         close (NEKSTAB_UNIT_DYNTOL)
         dmt_log_open = .false.
      end if

   end subroutine dmt

end module nekstab_dmt
