!-----------------------------------------------------------------------
! dmt.f90 — Dynamic Mode Tracking (DMT) for base flow stabilization
!
! Purpose:
!   Implements the Dynamic Mode Tracking technique of Queguineur et al.
!   (Phys. Fluids 31(3), 034101, 2019) for computing unstable steady states
!   by band-pass filtering the velocity field. The band-pass has zero gain
!   only at zero frequency, so a steady state is the only state without
!   force. A periodic orbit has content at its own frequencies and is not held.
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
! Filter bank (mode 6.4, 'dmt_bank'):
!   The same second-order filter, run without a force, isolates the spatial
!   content of the flow at chosen frequencies. Each filter is a damped
!   oscillator of natural frequency omega and damping beta driven by the
!   flow field u (dmt_filter.f90):
!     p' = q,   q' = u - beta*q - omega**2*p,   b = beta*q
!   b is u filtered around omega. The complex amplitude uhat follows from b
!   and its derivative, with u ~ Re[uhat exp(i*omega*t)]. The DNS is not
!   forced: the filters watch it. A frequency omega = 0 gives the mean.
!   The case registers the filters in nekStab_usrchk, on every rank:
!     use nekstab_dmt, only: dmt_bank_add, dmt_bank_options
!     call dmt_bank_add(omega1, beta1)   ! beta optional, default 0.05*omega
!     call dmt_bank_add(2.0d0*omega1)
!     call dmt_bank_options(hist_step=100, snapshots=.false., avg_from=t0)
!   Output, with k the number of the filter in the order of the calls:
!     bRe<session>0.f<k>, bIm<session>0.f<k> — Re and Im of uhat (all fields
!       of the state: velocity, and pressure and scalars if they are active),
!       at the end of the run, or averaged over [avg_from, end] if avg_from is set
!     dmt_bank.dat — time, |u|, |u - sum b_k|, |uhat_k| for each filter
!     b0k<session>0.f<n> — b_k at every output step of Nek5000, if snapshots=.true.
!
! Public interface:
!   dmt      — per-step driver; call from main.f90 fixed-point branch
!   dmt_init — allocator + first-call setup; also called lazily from dmt
!   dmt_bank_add     — register a filter of the bank (nekStab_usrchk)
!   dmt_bank_options — history interval, snapshots, averaging window (nekStab_usrchk)
!   dmt_bank         — per-step driver of mode 6.4; call from main.f90
!
! Dependencies:
!   nekstab_nek_bridge (wraps SIZE, TOTAL, ADJOINT), nekstab_krylov_subspace,
!   nekstab_vectors, nekstab_io, nekstab_dmt_filter
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
                                 ifbfcv, lastep, param, wdsize, ifoutfld, nof
   use nekstab_krylov_subspace, only: krylov_vector, k_zero, k_copy, k_cmult, &
                                      k_add2s2, k_norm
   use nekstab_vectors, only: nopcopy
   use nekstab_io, only: whereyouwant
   use nekstab_dmt_filter, only: dmt_filter_step_coeffs, dmt_filter_amp_coeffs
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
   real, save :: residu0, omc2, dt_ref
   integer, save :: n

   !  Filter bank (mode 6.4)
   integer, parameter :: dp = kind(0.0d0)
   integer, parameter :: BANK_MAX = 8 ! filters of the bank
   real, parameter :: BANK_BETA_REL = 0.05d0 ! default beta/omega
   !  Time for the start-up transient exp(-beta*t/2) to fall to 1e-3
   real, parameter :: BANK_SETTLE = 2.0d0*log(1.0d3)

   type :: bank_filter
      real :: omega, beta
      type(krylov_vector) :: p, q ! oscillator state
      type(krylov_vector) :: are, aim ! integral of Re uhat, Im uhat over the window
   end type bank_filter

   real, save :: bank_omega(BANK_MAX), bank_beta(BANK_MAX)
   integer, save :: bank_n = 0 ! filters registered
   integer, save :: bank_hist = 50 ! steps between lines of dmt_bank.dat
   integer, save :: bank_nstep = 0 ! steps taken by the bank
   logical, save :: bank_ready = .false.
   logical, save :: bank_snap = .false. ! write b_k at the output steps
   logical, save :: bank_avg = .false. ! average uhat from bank_avg_from
   real, save :: bank_avg_from = 0.0d0
   real, save :: bank_wsum = 0.0d0 ! length of the averaging window so far
   real, save :: bank_tprev = 0.0d0 ! time of the previous sample of u
   type(bank_filter), allocatable, save :: bank(:)
   !  Previous and current sample of u, and work vectors
   type(krylov_vector), allocatable, save :: bk_uprev, bk_ucur, bk_wp, bk_wq, bk_wr, bk_wi

   public :: dmt_init, dmt, dmt_bank_add, dmt_bank_options, dmt_bank

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
         dt_ref = dt
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

      !  The BDF2 stencil and the Euler step below take one time step for the three steps back.
      !  Nek5000 shortens the last step to land on endTime (0.0011 against 0.0100 on the cylinder):
      !  the residual l2/dt then jumped by a factor 10 to 100 in the last row of residu_dmt.dat.
      !  That step is skipped. Any other change of dt (variableDt = yes) gives a wrong stencil
      !  and stops the run.
      if (abs(dt - dt_ref) > 1.0d-6*dt_ref) then
         if (abs(time - param(10)) < 1.0d-6*dt_ref) return
         call nekStab_error('DMT needs a constant time step: set variableDt = no in the .par file')
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

      call opcmult(do1, do2, do3, -dmt_gan)
      if (time > dmt_skp) call opadd2(fcx, fcy, fcz, do1, do2, do3)

      if (nid == 0) then
         if (dmt_log_open) write (NEKSTAB_UNIT_DYNTOL, "(3E15.7)") time, residu, rate
         write (6, "(' DMT residu=',1pE11.4,' rate of change= ',1pE11.4)") residu, rate
         write (6, *) ' '
      end if

      if (residu < dmt_tol) then
         if (nid == 0) write (6, *) ' Converged DMT base flow to:', residu
         !  64-bit output, as SFD, BoostConv and Newton write it: a base flow read back
         !  in single precision carries 6e-8 of error into every stability stage after it.
         param(63) = 1.0d0
         call bcast(param(63), wdsize)
         call outpost(vx, vy, vz, pr, t, 'BF_')
         param(63) = 0.0d0
         call bcast(param(63), wdsize)
         ifbfcv = .true.
         lastep = 1
      end if

      if (istep == nsteps .and. nid == 0 .and. dmt_log_open) then
         close (NEKSTAB_UNIT_DYNTOL)
         dmt_log_open = .false.
      end if

   end subroutine dmt

!-----------------------------------------------------------------------
! dmt_bank_add — Register one filter of the bank
!
!   Call in nekStab_usrchk on every rank, once per frequency, at most
!   BANK_MAX times. The filters are numbered in the order of the calls.
!
!   omega — angular frequency in 1 per time unit (2*pi*f); 0 gives the mean
!   beta  — damping, 1 per time unit; the bandwidth of the filter. Default
!           0.05*omega. It must lie below omega. A filter settles in about
!           2*ln(1000)/beta time units (start-up transient exp(-beta*t/2));
!           a narrow filter separates close frequencies and settles slowly.
!           For omega = 0 it is required: the mean is taken over 1/beta.
!-----------------------------------------------------------------------
   subroutine dmt_bank_add(omega, beta)
      real, intent(in) :: omega
      real, intent(in), optional :: beta
      real :: b

      if (bank_ready) call nekStab_error( &
         'dmt_bank_add: call it in nekStab_usrchk, before the first time step')
      if (bank_n >= BANK_MAX) call nekStab_error( &
         'dmt_bank_add: the bank holds at most 8 filters')
      if (omega < 0.0d0) call nekStab_error('dmt_bank_add: omega must not be negative')
      if (present(beta)) then
         b = beta
      elseif (omega > 0.0d0) then
         b = BANK_BETA_REL*omega
      else
         call nekStab_error('dmt_bank_add: the mean (omega = 0) needs beta, the inverse of its averaging time')
      end if
      if (b <= 0.0d0) call nekStab_error('dmt_bank_add: beta must be positive')
      if (omega > 0.0d0 .and. b >= omega) call nekStab_error( &
         'dmt_bank_add: beta must be below omega, a filter wider than its centre frequency isolates nothing')

      bank_n = bank_n + 1
      bank_omega(bank_n) = omega
      bank_beta(bank_n) = b

   end subroutine dmt_bank_add

!-----------------------------------------------------------------------
! dmt_bank_options — Output options of the bank (nekStab_usrchk)
!
!   hist_step — steps between lines of dmt_bank.dat (default 50)
!   snapshots — write b_k, the field filtered around omega_k, at every
!               output step of Nek5000 (files b01, b02, ...; default .false.)
!   avg_from  — average uhat over [avg_from, end of run] and write the
!               average. A filter also passes part of the other frequencies
!               of the flow, about beta/omega of their amplitude, which
!               shows as a ripple in uhat at the difference of the
!               frequencies. For a flow with harmonics of one frequency the
!               average over a whole number of its periods removes the
!               ripple: in a test with harmonics of amplitude 0.10, 0.03 and
!               0.01 at beta = 0.1*omega, the last step differed from the
!               exact amplitude by 5 %, 15 % and 45 % of each harmonic, the
!               average over 10 periods by 4e-5, 5e-5 and 1e-4. Without
!               avg_from the files hold uhat at the last step.
!-----------------------------------------------------------------------
   subroutine dmt_bank_options(hist_step, snapshots, avg_from)
      integer, intent(in), optional :: hist_step
      logical, intent(in), optional :: snapshots
      real, intent(in), optional :: avg_from

      if (present(hist_step)) bank_hist = max(hist_step, 1)
      if (present(snapshots)) bank_snap = snapshots
      if (present(avg_from)) then
         bank_avg = .true.
         bank_avg_from = avg_from
      end if

   end subroutine dmt_bank_options

!-----------------------------------------------------------------------
! dmt_bank — Per-step driver of the filter bank (mode 6.4)
!
!   Called from nekStab at every step after the Navier-Stokes solve. At
!   istep = 0 the filters start from the equilibrium of the initial field.
!   At the last step (lastep = 1) it writes the amplitude fields.
!-----------------------------------------------------------------------
   subroutine dmt_bank

      if (bank_n == 0) call nekStab_error( &
         'mode dmt_bank: no filter. Call dmt_bank_add(omega, beta) in nekStab_usrchk')
      if (.not. bank_ready) call bank_init
      if (istep == 0) return

      call bank_step
      if (mod(istep, bank_hist) == 0 .or. lastep == 1) call bank_history
      if (bank_snap .and. ifoutfld) call bank_snapshots
      if (lastep == 1) call bank_finish

   end subroutine dmt_bank

!-----------------------------------------------------------------------
! bank_init — Allocate the filters and start them from the initial field
!
!   A constant input u0 holds the oscillator at p = u0/omega**2, q = 0, so
!   the mean of the flow starts no transient. The mean filter (omega = 0)
!   starts at q = u0/beta, which is b = u0.
!-----------------------------------------------------------------------
   subroutine bank_init
      integer :: k
      real :: tend, tsettle

      allocate (bank(bank_n))
      allocate (bk_uprev, bk_ucur, bk_wp, bk_wq, bk_wr, bk_wi)
      call k_zero(bk_uprev)
      call k_zero(bk_ucur)
      call k_zero(bk_wp)
      call k_zero(bk_wq)
      call k_zero(bk_wr)
      call k_zero(bk_wi)

      call nopcopy(bk_uprev%vx, bk_uprev%vy, bk_uprev%vz, bk_uprev%pr, bk_uprev%t, &
                   vx, vy, vz, pr, t)

      do k = 1, bank_n
         bank(k)%omega = bank_omega(k)
         bank(k)%beta = bank_beta(k)
         call k_zero(bank(k)%p)
         call k_zero(bank(k)%q)
         call k_zero(bank(k)%are)
         call k_zero(bank(k)%aim)
         if (bank(k)%omega > 0.0d0) then
            call k_copy(bank(k)%p, bk_uprev)
            call k_cmult(bank(k)%p, 1.0d0/bank(k)%omega**2)
         else
            call k_copy(bank(k)%q, bk_uprev)
            call k_cmult(bank(k)%q, 1.0d0/bank(k)%beta)
         end if
      end do

      bank_tprev = time
      bank_nstep = 0
      bank_wsum = 0.0d0

      !  Length of the run, to warn about filters that cannot settle
      tend = param(10) - time
      if (tend <= 0.0d0) tend = nsteps*dt

      if (nid == 0) then
         call nekStab_log('DMT filter bank (mode 6.4)')
         write (6, '(a,i0,a,es10.3,a)') '   filters:', bank_n, '   run length:', tend, ' time units'
         write (6, '(a)') '   k      omega       beta     f=omega/2pi   settling time'
         do k = 1, bank_n
            tsettle = BANK_SETTLE/bank(k)%beta
            write (6, '(i4,4es12.4)') k, bank(k)%omega, bank(k)%beta, &
               bank(k)%omega/(8.0d0*atan(1.0d0)), tsettle
            if (tsettle > tend) then
               write (6, '(a,i0,a)') '   WARNING: the run is shorter than the settling time of filter ', k, &
                  ': its amplitude has not reached 1e-3 of the final value'
            end if
         end do
         if (bank_avg) write (6, '(a,es10.3)') '   amplitude averaged from t =', bank_avg_from
         open (unit=NEKSTAB_UNIT_DYNTOL, file='dmt_bank.dat', status='replace')
         write (NEKSTAB_UNIT_DYNTOL, '(a)', advance='no') '# time  |u|  |u-sum(b_k)|  |uhat_k| for k = 1..'
         write (NEKSTAB_UNIT_DYNTOL, '(i0)') bank_n
      end if

      bank_ready = .true.

   end subroutine bank_init

!-----------------------------------------------------------------------
! bank_step — Advance every filter over the step just taken
!
!   z(n+1) = phi*z(n) + a*u(n) + b*u(n+1), z = (p, q). The step length is
!   taken from the time, so a variable time step and the shortened last step
!   need no special case.
!-----------------------------------------------------------------------
   subroutine bank_step
      real(dp) :: h, phi(2, 2), a(2), b(2)
      integer :: k

      h = time - bank_tprev
      if (h <= 0.0d0) return

      call nopcopy(bk_ucur%vx, bk_ucur%vy, bk_ucur%vz, bk_ucur%pr, bk_ucur%t, &
                   vx, vy, vz, pr, t)

      do k = 1, bank_n
         call dmt_filter_step_coeffs(real(bank(k)%omega, dp), real(bank(k)%beta, dp), h, phi, a, b)

         call k_copy(bk_wp, bank(k)%p)
         call k_cmult(bk_wp, real(phi(1, 1)))
         call k_add2s2(bk_wp, bank(k)%q, real(phi(1, 2)))
         call k_add2s2(bk_wp, bk_uprev, real(a(1)))
         call k_add2s2(bk_wp, bk_ucur, real(b(1)))

         call k_copy(bk_wq, bank(k)%p)
         call k_cmult(bk_wq, real(phi(2, 1)))
         call k_add2s2(bk_wq, bank(k)%q, real(phi(2, 2)))
         call k_add2s2(bk_wq, bk_uprev, real(a(2)))
         call k_add2s2(bk_wq, bk_ucur, real(b(2)))

         call k_copy(bank(k)%p, bk_wp)
         call k_copy(bank(k)%q, bk_wq)
      end do

      if (bank_avg .and. time > bank_avg_from) then
         do k = 1, bank_n
            call bank_amplitude(k, bk_wr, bk_wi)
            call k_add2s2(bank(k)%are, bk_wr, real(h))
            call k_add2s2(bank(k)%aim, bk_wi, real(h))
         end do
         bank_wsum = bank_wsum + h
      end if

      call k_copy(bk_uprev, bk_ucur)
      bank_tprev = time
      bank_nstep = bank_nstep + 1

   end subroutine bank_step

!-----------------------------------------------------------------------
! bank_amplitude — Re and Im of uhat for filter k at the present time
!
!   Uses the state of the filter and the present field u (bk_ucur).
!-----------------------------------------------------------------------
   subroutine bank_amplitude(k, re, im)
      integer, intent(in) :: k
      type(krylov_vector), intent(inout) :: re, im
      real(dp) :: cre(3), cim(3)

      call dmt_filter_amp_coeffs(real(bank(k)%omega, dp), real(bank(k)%beta, dp), real(time, dp), cre, cim)

      call k_copy(re, bank(k)%q)
      call k_cmult(re, real(cre(1)))
      call k_add2s2(re, bank(k)%p, real(cre(2)))
      call k_add2s2(re, bk_ucur, real(cre(3)))

      call k_copy(im, bank(k)%q)
      call k_cmult(im, real(cim(1)))
      call k_add2s2(im, bank(k)%p, real(cim(2)))
      call k_add2s2(im, bk_ucur, real(cim(3)))

   end subroutine bank_amplitude

!-----------------------------------------------------------------------
! bank_history — One line of dmt_bank.dat
!
!   |uhat_k| = |Z_k| = beta*sqrt(|q|**2 + |q'/omega|**2), the modulus of the
!   complex amplitude, which does not depend on the phase. |u - sum b_k| is
!   the part of the flow that no filter holds.
!-----------------------------------------------------------------------
   subroutine bank_history
      real :: nu, nres, nq, ndq, amp(BANK_MAX)
      integer :: k

      call k_norm(nu, bk_ucur)
      call k_copy(bk_wr, bk_ucur)
      do k = 1, bank_n
         call k_add2s2(bk_wr, bank(k)%q, -bank(k)%beta)
      end do
      call k_norm(nres, bk_wr)

      do k = 1, bank_n
         call k_norm(nq, bank(k)%q)
         if (bank(k)%omega > 0.0d0) then
            call k_copy(bk_wq, bk_ucur)
            call k_add2s2(bk_wq, bank(k)%q, -bank(k)%beta)
            call k_add2s2(bk_wq, bank(k)%p, -bank(k)%omega**2)
            call k_norm(ndq, bk_wq)
            amp(k) = bank(k)%beta*sqrt(nq**2 + (ndq/bank(k)%omega)**2)
         else
            amp(k) = bank(k)%beta*nq
         end if
      end do

      if (nid == 0) then
         write (NEKSTAB_UNIT_DYNTOL, '(*(es15.7))') time, nu, nres, amp(1:bank_n)
         flush (NEKSTAB_UNIT_DYNTOL)
      end if

   end subroutine bank_history

!-----------------------------------------------------------------------
! bank_snapshots — b_k = beta_k*q_k at an output step of Nek5000
!-----------------------------------------------------------------------
   subroutine bank_snapshots
      character(len=3) :: pre
      integer :: k

      do k = 1, bank_n
         call k_copy(bk_wr, bank(k)%q)
         call k_cmult(bk_wr, bank(k)%beta)
         write (pre, '(a,i2.2)') 'b', k
         call outpost2(bk_wr%vx, bk_wr%vy, bk_wr%vz, bk_wr%pr, bk_wr%t, nof, pre)
      end do

   end subroutine bank_snapshots

!-----------------------------------------------------------------------
! bank_finish — Write the amplitude fields at the end of the run
!
!   File k of bRe and bIm belongs to the k-th call of dmt_bank_add.
!-----------------------------------------------------------------------
   subroutine bank_finish
      integer :: k

      if (bank_nstep == 0) return

      do k = 1, bank_n
         if (bank_avg .and. bank_wsum > 0.0d0) then
            call k_copy(bk_wr, bank(k)%are)
            call k_cmult(bk_wr, 1.0d0/bank_wsum)
            call k_copy(bk_wi, bank(k)%aim)
            call k_cmult(bk_wi, 1.0d0/bank_wsum)
         else
            call bank_amplitude(k, bk_wr, bk_wi)
         end if
         call whereyouwant('bRe', k)
         call outpost2(bk_wr%vx, bk_wr%vy, bk_wr%vz, bk_wr%pr, bk_wr%t, nof, 'bRe')
         call whereyouwant('bIm', k)
         call outpost2(bk_wi%vx, bk_wi%vy, bk_wi%vz, bk_wi%pr, bk_wi%t, nof, 'bIm')
      end do

      if (nid == 0) then
         close (NEKSTAB_UNIT_DYNTOL)
         if (bank_avg .and. bank_wsum > 0.0d0) then
            write (6, '(a,es10.3,a)') ' DMT bank: amplitude averaged over', bank_wsum, ' time units'
         else
            write (6, '(a)') ' DMT bank: amplitude at the last step'
         end if
         write (6, '(a,i0,a)') ' DMT bank: wrote bRe and bIm for ', bank_n, ' filters'
      end if

   end subroutine bank_finish

end module nekstab_dmt
