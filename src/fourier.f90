!-----------------------------------------------------------------------
! fourier.f90 — Fourier representation of a closed periodic orbit
!
! Purpose:
!   Holds a periodic base flow as K Fourier modes instead of N stored
!   snapshots (iffourierorbit). The modes are accumulated during one
!   integration of the orbit (a DFT that adds one sample per time step),
!   written by a converged Newton run, and read back by Floquet, energy
!   budget and OTD runs:
!
!     U(t) = c_0 + 2 Re sum_{k=1}^{K} c_k exp(2 pi i k t / T)
!
!   Only a closed orbit is held this way. A Newton guess does not close.
!   Its periodic extension has a jump, so |c_k| falls as 1/k.
!   Newton keeps its snapshots and writes the modes once it converges.
!   Also re-exports the FFTW primitives used by POD and SPOD.
!
! Public interface:
!   fourier_orbit_available          — modes are ready for evaluation
!   fourier_orbit_init_streaming     — start the DFT over one period
!   fourier_orbit_stream_step        — add the current fields as sample istep
!   fourier_orbit_finalize_streaming — select K, keep modes 0..K
!   fourier_orbit_eval               — set vx, vy, vz, pr, t to U(t)
!   fourier_orbit_write              — FOU<session>0.f* and <session>.fourier
!   fourier_orbit_load               — read them back
!
! Dependencies:
!   nekstab_nek_bridge, nekstab_vectors, nekstab_fourier_fftw,
!   nekstab_fourier_series
!-----------------------------------------------------------------------
module nekstab_fourier
   use nekstab_nek_bridge, only: lx1, ly1, lz1, lx2, ly2, lz2, lelv, lelt, &
                                 ldimt, nx1, ny1, nz1, nelv, nelt, if3D, &
                                 ifto, ifpsco, nid, nekStab_error, bm1s, &
                                 vx, vy, vz, pr, t, param, SESSION, wdsize, &
                                 lsize, isize, thermal_norm_weight, ifKEnorm, &
                                 iffourierorbit, fourier_kmax, nof, time
   use nekstab_vectors, only: nopcopy, nopadd2s2, noprzero
   use nekstab_fourier_fftw, only: fourier_decomposition, fourier_reconstruction, &
                                   fft_init, fft_cleanup, fft_r2c, fft_c2r, &
                                   fft_frequencies
   use nekstab_fourier_series, only: fourier_k_accum, fourier_dft_weight, &
                                     fourier_eval_weight, fourier_select_k
   implicit none
   private

   ! Re-exported procedures from low-level FFTW module
   public :: fourier_decomposition, fourier_reconstruction
   public :: fft_init, fft_cleanup, fft_r2c, fft_c2r, fft_frequencies

   public :: fourier_orbit_available, fourier_orbit_init_streaming, &
             fourier_orbit_stream_step, fourier_orbit_finalize_streaming, &
             fourier_orbit_eval, fourier_orbit_write, fourier_orbit_load

   integer, parameter :: lv = lx1*ly1*lz1*lelv
   integer, parameter :: lp = lx2*ly2*lz2*lelv
   integer, parameter :: lt = lx1*ly1*lz1*lelt

   logical, save :: loaded = .false., streaming = .false.
   real, save :: period = 0.0d0, sample_dt = 0.0d0, missing_energy = 0.0d0
   integer, save :: nsamples = 0, K = 0
   real, allocatable, save :: Ek(:)

   ! Mode m of a field: m = 0 is c_0, m = 2k-1 is Re(c_k), m = 2k is Im(c_k).
   ! Each mode holds the full state, so the field operations of nekstab_vectors
   ! (and outpost2 / load_fld) work on it directly.
   real, allocatable, save :: fou_u(:, :), fou_v(:, :), fou_w(:, :), fou_p(:, :)
   real, allocatable, save :: fou_t(:, :, :)

contains

   logical function fourier_orbit_available()
      fourier_orbit_available = iffourierorbit .and. loaded
   end function fourier_orbit_available

   !-----------------------------------------------------------------------
   ! allocate_modes — zeroed storage for modes 0..2*kmodes
   !-----------------------------------------------------------------------
   subroutine allocate_modes(kmodes)
      integer, intent(in) :: kmodes
      integer :: m

      if (allocated(fou_u)) deallocate (fou_u, fou_v, fou_w, fou_p, fou_t)
      allocate (fou_u(lv, 0:2*kmodes), fou_v(lv, 0:2*kmodes), fou_w(lv, 0:2*kmodes), &
                fou_p(lp, 0:2*kmodes), fou_t(lt, ldimt, 0:2*kmodes))
      do m = 0, 2*kmodes
         call noprzero(fou_u(:, m), fou_v(:, m), fou_w(:, m), fou_p(:, m), fou_t(:, :, m))
      end do
      loaded = .false.
   end subroutine allocate_modes

   !-----------------------------------------------------------------------
   ! mode_energy — <c, c> with the mass matrix, as in the Krylov norm
   !-----------------------------------------------------------------------
   real function mode_energy(m)
      integer, intent(in) :: m
      real :: glsc3
      integer :: n, nt, i

      n = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt
      mode_energy = glsc3(fou_u(1, m), fou_u(1, m), bm1s, n) + glsc3(fou_v(1, m), fou_v(1, m), bm1s, n)
      if (if3D) mode_energy = mode_energy + glsc3(fou_w(1, m), fou_w(1, m), bm1s, n)
      if (ifKEnorm) return
      if (ifto) mode_energy = mode_energy + thermal_norm_weight*glsc3(fou_t(1, 1, m), fou_t(1, 1, m), bm1s, nt)
      do i = 2, ldimt
         if (ifpsco(i - 1)) mode_energy = mode_energy + glsc3(fou_t(1, i, m), fou_t(1, i, m), bm1s, nt)
      end do
   end function mode_energy

   !-----------------------------------------------------------------------
   ! fourier_orbit_init_streaming — start the DFT of one period
   !   N samples U(dt), ..., U(N dt) with N dt = T; U(N dt) stands for U(0).
   !   At most (N-1)/2 modes: below the Nyquist mode every c_k has a
   !   distinct conjugate c_{N-k}, so the factor 2 in U(t) holds for all k.
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_init_streaming(nsteps_in, T_in)
      integer, intent(in) :: nsteps_in
      real, intent(in) :: T_in

      nsamples = nsteps_in
      period = T_in
      sample_dt = T_in/real(nsteps_in)
      K = fourier_k_accum(nsamples, fourier_kmax)
      call allocate_modes(K)
      streaming = .true.
      if (nid == 0) write (6, "('Fourier orbit: accumulating',I5,' modes over',I7,' steps')") K, nsamples
   end subroutine fourier_orbit_init_streaming

   !-----------------------------------------------------------------------
   ! fourier_orbit_stream_step — c_k += U(istep dt) exp(-2 pi i k istep / N) / N
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_stream_step(istep_in)
      integer, intent(in) :: istep_in
      real :: w_re, w_im
      integer :: kk

      if (.not. streaming) return
      call fourier_dft_weight(istep_in, nsamples, 0, w_re, w_im)
      call nopadd2s2(fou_u(:, 0), fou_v(:, 0), fou_w(:, 0), fou_p(:, 0), fou_t(:, :, 0), &
                     vx, vy, vz, pr, t, w_re)
      do kk = 1, K
         call fourier_dft_weight(istep_in, nsamples, kk, w_re, w_im)
         call nopadd2s2(fou_u(:, 2*kk - 1), fou_v(:, 2*kk - 1), fou_w(:, 2*kk - 1), &
                        fou_p(:, 2*kk - 1), fou_t(:, :, 2*kk - 1), &
                        vx, vy, vz, pr, t, w_re)
         call nopadd2s2(fou_u(:, 2*kk), fou_v(:, 2*kk), fou_w(:, 2*kk), &
                        fou_p(:, 2*kk), fou_t(:, :, 2*kk), &
                        vx, vy, vz, pr, t, w_im)
      end do
   end subroutine fourier_orbit_stream_step

   !-----------------------------------------------------------------------
   ! fourier_orbit_finalize_streaming — choose K and keep modes 0..K
   !
   !   K starts at 8 and doubles while a doubling cuts the missing
   !   fluctuation energy by 10 or more, down to 1e-10. A fixed
   !   fourier_kmax keeps every accumulated mode.
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_finalize_streaming()
      integer :: kk, Kacc, ntrace, i
      integer :: trace_k(32), trace_next(32)
      real :: trace_miss(32)

      if (.not. streaming) return
      streaming = .false.
      Kacc = K

      if (allocated(Ek)) deallocate (Ek)
      allocate (Ek(0:Kacc))
      Ek(0) = mode_energy(0)
      do kk = 1, Kacc
         Ek(kk) = mode_energy(2*kk - 1) + mode_energy(2*kk)
      end do
      call fourier_select_k(Ek, Kacc, fourier_kmax, K, missing_energy, &
                            ntrace, trace_k, trace_next, trace_miss)
      if (nid == 0) then
         do i = 1, ntrace
            write (6, "('Fourier orbit doubling: K =',I4,' -> ',I4,' missing energy =',1PE12.4)") &
               trace_k(i), trace_next(i), trace_miss(i)
         end do
      end if

      ! Free the modes above K.
      if (K < Kacc) then
         call shrink(fou_u, 2*K)
         call shrink(fou_v, 2*K)
         call shrink(fou_w, 2*K)
         call shrink(fou_p, 2*K)
         call shrink_t(2*K)
      end if
      loaded = .true.
      if (nid == 0) write (6, "('Fourier orbit selected K =',I4,' missing energy =',1PE12.4)") K, missing_energy
   end subroutine fourier_orbit_finalize_streaming

   subroutine shrink(a, mlast)
      real, allocatable, intent(inout) :: a(:, :)
      integer, intent(in) :: mlast
      real, allocatable :: keep(:, :)
      allocate (keep(size(a, 1), 0:mlast))
      keep = a(:, 0:mlast)
      call move_alloc(keep, a)
   end subroutine shrink

   subroutine shrink_t(mlast)
      integer, intent(in) :: mlast
      real, allocatable :: keep(:, :, :)
      allocate (keep(lt, ldimt, 0:mlast))
      keep = fou_t(:, :, 0:mlast)
      call move_alloc(keep, fou_t)
   end subroutine shrink_t

   !-----------------------------------------------------------------------
   ! fourier_orbit_eval — set the Nek fields to U(t_in), t_in taken modulo T
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_eval(t_in)
      real, intent(in) :: t_in
      real :: w_re, w_im
      integer :: kk

      if (.not. loaded) return
      call nopcopy(vx, vy, vz, pr, t, fou_u(:, 0), fou_v(:, 0), fou_w(:, 0), fou_p(:, 0), fou_t(:, :, 0))
      do kk = 1, K
         call fourier_eval_weight(kk, t_in, period, w_re, w_im)
         call nopadd2s2(vx, vy, vz, pr, t, fou_u(:, 2*kk - 1), fou_v(:, 2*kk - 1), fou_w(:, 2*kk - 1), &
                        fou_p(:, 2*kk - 1), fou_t(:, :, 2*kk - 1), w_re)
         call nopadd2s2(vx, vy, vz, pr, t, fou_u(:, 2*kk), fou_v(:, 2*kk), fou_w(:, 2*kk), &
                        fou_p(:, 2*kk), fou_t(:, :, 2*kk), w_im)
      end do
   end subroutine fourier_orbit_eval

   !-----------------------------------------------------------------------
   ! fourier_orbit_write — header <session>.fourier, then FOU<session>0.f00001
   !   (c_0), f00002/f00003 (Re, Im of c_1), ... in double precision
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_write()
      integer :: io_unit, kk, m
      real :: p63_sav

      if (.not. loaded) return
      if (nid == 0) then
         open (newunit=io_unit, file=trim(SESSION)//'.fourier', status='replace', action='write')
         write (io_unit, "(A,1PE22.14)") "period: ", period
         write (io_unit, "(A,I10)") "N:      ", nsamples
         write (io_unit, "(A,1PE22.14)") "dt:     ", sample_dt
         write (io_unit, "(A,I10)") "K:      ", K
         write (io_unit, "(A,1PE22.14)") "missing_energy: ", missing_energy
         write (io_unit, "(A)") "# mode  energy"
         do kk = 0, K
            write (io_unit, "(I6,1PE22.14)") kk, Ek(kk)
         end do
         close (io_unit)
      end if

      p63_sav = param(63)
      param(63) = 1.0d0
      do m = 0, 2*K
         call outpost2(fou_u(1, m), fou_v(1, m), fou_w(1, m), fou_p(1, m), fou_t(1, 1, m), nof, 'FOU')
      end do
      param(63) = p63_sav
      if (nid == 0) write (6, "('Wrote ',A,'.fourier and ',I4,' mode files FOU',A,'0.f*')") &
         trim(SESSION), 2*K + 1, trim(SESSION)
   end subroutine fourier_orbit_write

   !-----------------------------------------------------------------------
   ! fourier_orbit_load — read the files of fourier_orbit_write
   !   loaded_out = .false. when they are not there. The Nek fields are
   !   overwritten; the caller sets its initial condition afterwards.
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_load(loaded_out)
      logical, intent(out) :: loaded_out
      logical :: exists
      integer :: io_unit, ierr, kk, idx
      character(len=132) :: line, fname
      real :: time_sav, param10_sav

      exists = .false.
      if (nid == 0) then
         inquire (file='FOU'//trim(SESSION)//'0.f00001', exist=exists)
         if (exists) inquire (file=trim(SESSION)//'.fourier', exist=exists)
      end if
      call bcast(exists, lsize)
      loaded_out = exists
      if (.not. exists) return

      if (nid == 0) then
         open (newunit=io_unit, file=trim(SESSION)//'.fourier', status='old', action='read', iostat=ierr)
         if (ierr /= 0) call nekStab_error('could not open '//trim(SESSION)//'.fourier')
         read (io_unit, "(A)") line; read (line(index(line, ':') + 1:), *) period
         read (io_unit, "(A)") line; read (line(index(line, ':') + 1:), *) nsamples
         read (io_unit, "(A)") line; read (line(index(line, ':') + 1:), *) sample_dt
         read (io_unit, "(A)") line; read (line(index(line, ':') + 1:), *) K
         read (io_unit, "(A)") line; read (line(index(line, ':') + 1:), *) missing_energy
      end if
      call bcast(period, wdsize)
      call bcast(nsamples, isize)
      call bcast(sample_dt, wdsize)
      call bcast(K, isize)
      call bcast(missing_energy, wdsize)

      if (allocated(Ek)) deallocate (Ek)
      allocate (Ek(0:K))
      if (nid == 0) then
         read (io_unit, "(A)") line ! comment line
         do kk = 0, K
            read (io_unit, *) idx, Ek(kk)
         end do
         close (io_unit)
      end if
      call bcast(Ek, (K + 1)*wdsize)

      if (abs(abs(param(12)) - sample_dt) > 1.0d-7*sample_dt .and. nid == 0) then
         write (6, *) 'nekStab: NOTE: dt', abs(param(12)), ' differs from the dt of the modes', sample_dt
         write (6, *) '  U(t) is continuous in t; the neutral multiplier moves by the time-step error.'
      end if

      call allocate_modes(K)
      ! load_fld also sets the time and endTime of the run from each file.
      time_sav = time
      param10_sav = param(10)
      do idx = 1, 2*K + 1
         write (fname, '(A,A,"0.f",I5.5)') 'FOU', trim(SESSION), idx
         call load_fld(fname)
         call nopcopy(fou_u(:, idx - 1), fou_v(:, idx - 1), fou_w(:, idx - 1), fou_p(:, idx - 1), &
                      fou_t(:, :, idx - 1), vx, vy, vz, pr, t)
      end do
      time = time_sav
      param(10) = param10_sav
      loaded = .true.
      if (nid == 0) write (6, "('Fourier orbit loaded: K =',I4,' N =',I6,' T =',1PE15.7,' missing energy =',1PE12.4)") &
         K, nsamples, period, missing_energy
   end subroutine fourier_orbit_load

end module nekstab_fourier
