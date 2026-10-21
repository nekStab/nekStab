!-----------------------------------------------------------------------
! fourier.f90 — Fourier orbit representation and modal FFT interface
!
! Purpose:
!   Provides Fourier mode representation for periodic orbits in nekStab:
!   - Storage and evaluation of periodic orbits in Fourier modes (U(t))
!   - Building modes via FFTW in time on stored orbit (all N/2 modes)
!   - Building modes via streaming DFT when orbit is not stored in memory
!   - Automatic K selection by doubling rule (missing energy < 1e-10)
!   - Loading/saving FOU<session>0.f* field files and <session>.fourier header
!   - Re-exports low-level FFTW primitives for modal POD/SPOD
!
! Public interface:
!   fourier_orbit_available        — check if Fourier modes are ready
!   fourier_orbit_load             — load modes from FOU* and .fourier
!   fourier_orbit_write            — write modes to FOU* and .fourier
!   fourier_orbit_init_streaming   — start streaming DFT accumulation
!   fourier_orbit_stream_step      — accumulate one step of streaming DFT
!   fourier_orbit_finalize_streaming — finish streaming DFT, select K
!   fourier_orbit_build_from_stored — build modes from uor via FFTW
!   fourier_orbit_eval             — evaluate U(t) from Fourier series
!   fourier_orbit_cleanup          — deallocate module mode storage
!   fourier_orbit_reset_trial      — reset trial orbit state for Newton PO
!   fourier_orbit_is_streaming     — check if streaming DFT is active
!   fourier_orbit_get_K            — get number of active Fourier modes
!   fourier_orbit_get_missing_energy — get missing fluctuation energy
!   fourier_orbit_get_period       — get orbit period
!
!   Re-exports from nekstab_fourier_fftw:
!   fft_init, fft_cleanup, fft_r2c, fft_c2r, fft_frequencies,
!   fourier_decomposition, fourier_reconstruction
!
! Dependencies:
!   nekstab_nek_bridge, nekstab_fourier_fftw, iso_c_binding
!-----------------------------------------------------------------------
module nekstab_fourier
   use nekstab_nek_bridge, only: lx1, ly1, lz1, lx2, ly2, lz2, lelv, lelt, &
                                 ldimt, nx1, ny1, nz1, nelv, nelt, if3D, if3d, &
                                 ifto, ifpsco, nid, nekStab_error, bm1s, &
                                 vx, vy, vz, pr, t, param, SESSION, wdsize, &
                                 lsize, isize, thermal_norm_weight, ifKEnorm, &
                                 iffourierorbit, fourier_kmax, fourier_mem_gb, &
                                 nof, ndim
   use nekstab_fourier_fftw, only: fourier_decomposition, fourier_reconstruction, &
                                   fft_init, fft_cleanup, fft_r2c, fft_c2r, &
                                   fft_frequencies
   use, intrinsic :: iso_c_binding, only: c_double, c_double_complex
   implicit none
   private

   ! Re-exported procedures from low-level FFTW module
   public :: fourier_decomposition, fourier_reconstruction
   public :: fft_init, fft_cleanup, fft_r2c, fft_c2r, fft_frequencies

   ! Public procedures for Fourier orbit representation
   public :: fourier_orbit_available, fourier_orbit_load, fourier_orbit_write, &
             fourier_orbit_init_streaming, fourier_orbit_stream_step, &
             fourier_orbit_finalize_streaming, fourier_orbit_build_from_stored, &
             fourier_orbit_eval, fourier_orbit_cleanup, fourier_orbit_reset_trial, &
             fourier_orbit_is_streaming, fourier_orbit_get_K, &
             fourier_orbit_get_missing_energy, fourier_orbit_get_period

   ! Grid dimensions (fixed at compile time)
   integer, parameter :: lv = lx1*ly1*lz1*lelv
   integer, parameter :: lt = lx1*ly1*lz1*lelt
   real(c_double), parameter :: TWOPI = 6.283185307179586476925286766559d0

   ! Module state for Fourier orbit
   logical, save :: fourier_orbit_loaded = .false.
   real, save :: fourier_T = 0.0d0
   integer, save :: fourier_N = 0
   real, save :: fourier_dt = 0.0d0
   integer, save :: fourier_K = 0
   integer, save :: fourier_replay_K = 0
   real, save :: fourier_missing_energy = 0.0d0
   real, allocatable, save :: fourier_Ek(:)

   ! Mode coefficient storage:
   ! Slice 0: c_0 (DC mean component)
   ! Slice 2*k-1: Re(c_k)
   ! Slice 2*k:   Im(c_k)
   real, allocatable, save :: fou_u(:, :), fou_v(:, :), fou_w(:, :)
   real, allocatable, save :: fou_t(:, :, :)

   ! Streaming DFT state
   logical, save :: streaming_active = .false.
   integer, save :: stream_step = 0
   integer, save :: stream_N = 0
   real, save :: stream_T = 0.0d0
   integer, save :: stream_Kcap = 128

contains

   !-----------------------------------------------------------------------
   ! fourier_orbit_available — check if Fourier orbit is ready
   !-----------------------------------------------------------------------
   logical function fourier_orbit_available()
      fourier_orbit_available = iffourierorbit .and. fourier_orbit_loaded
   end function fourier_orbit_available

   !-----------------------------------------------------------------------
   ! fourier_orbit_is_streaming — check if streaming DFT is active
   !-----------------------------------------------------------------------
   logical function fourier_orbit_is_streaming()
      fourier_orbit_is_streaming = iffourierorbit .and. streaming_active
   end function fourier_orbit_is_streaming

   !-----------------------------------------------------------------------
   ! Getters
   !-----------------------------------------------------------------------
   integer function fourier_orbit_get_K()
      fourier_orbit_get_K = fourier_K
   end function fourier_orbit_get_K

   real function fourier_orbit_get_missing_energy()
      fourier_orbit_get_missing_energy = fourier_missing_energy
   end function fourier_orbit_get_missing_energy

   real function fourier_orbit_get_period()
      fourier_orbit_get_period = fourier_T
   end function fourier_orbit_get_period

   !-----------------------------------------------------------------------
   ! fourier_orbit_reset_trial — reset trial orbit state for Newton PO
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_reset_trial()
      fourier_orbit_loaded = .false.
   end subroutine fourier_orbit_reset_trial

   !-----------------------------------------------------------------------
   ! fourier_orbit_cleanup — deallocate module arrays
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_cleanup()
      if (allocated(fou_u)) deallocate (fou_u)
      if (allocated(fou_v)) deallocate (fou_v)
      if (allocated(fou_w)) deallocate (fou_w)
      if (allocated(fou_t)) deallocate (fou_t)
      if (allocated(fourier_Ek)) deallocate (fourier_Ek)
      fourier_orbit_loaded = .false.
      streaming_active = .false.
      fourier_K = 0
      fourier_replay_K = 0
      fourier_N = 0
      fourier_T = 0.0d0
      fourier_dt = 0.0d0
      fourier_missing_energy = 0.0d0
   end subroutine fourier_orbit_cleanup

   !-----------------------------------------------------------------------
   ! set_replay_depth — file K stays the doubling-rule cutoff
   !
   ! Replay uses every coefficient already computed when that cutoff still
   ! leaves a jump tail. K = N/2 interpolates the samples; the file K does not.
   !-----------------------------------------------------------------------
   subroutine set_replay_depth(k_have)
      integer, intent(in) :: k_have

      if (fourier_missing_energy >= 1.0d-10 .and. k_have > fourier_K) then
         fourier_replay_K = k_have
         if (nid == 0) then
            write (6, "('Fourier orbit replay keeps K =',I5,' so sample times match; file K =',I5)") &
               fourier_replay_K, fourier_K
         end if
      else
         fourier_replay_K = fourier_K
      end if
   end subroutine set_replay_depth


   !-----------------------------------------------------------------------
   ! compute_mode_energy — mass-weighted norm squared: <q, q>_bm1s
   !-----------------------------------------------------------------------
   subroutine compute_mode_energy(energy, ux, uy, uz, ut)
      real, dimension(lv), intent(in) :: ux, uy, uz
      real, dimension(lt, ldimt), intent(in) :: ut
      real, intent(out) :: energy

      real :: glsc3, temp_inner
      integer :: nv, nt, m

      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

      energy = glsc3(ux, ux, bm1s, nv) + glsc3(uy, uy, bm1s, nv)
      if (if3D) energy = energy + glsc3(uz, uz, bm1s, nv)
      if (ifto .and. .not. ifKEnorm) then
         temp_inner = glsc3(ut(:, 1), ut(:, 1), bm1s, nt)
         energy = energy + thermal_norm_weight*temp_inner
      end if
      if (ldimt > 1 .and. .not. ifKEnorm) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) then
               energy = energy + glsc3(ut(:, m), ut(:, m), bm1s, nt)
            end if
         end do
      end if
   end subroutine compute_mode_energy

   !-----------------------------------------------------------------------
   ! select_K_doubling — automatic K by doubling rule of plan
   !
   ! Rule:
   ! - Never more than N/2
   ! - Start K = 8
   ! - Double K while doubling cuts missing fluctuation energy by 10x or more
   ! - Stop when missing energy < 1e-10 or a doubling gains less than 10x
   ! - Log energy curve per doubling
   ! - User overrides: fourier_kmax (fixed K) and fourier_mem_gb (cap)
   !-----------------------------------------------------------------------
   subroutine select_K_doubling(Ek, Kmax, N, K_out, Emiss_out)
      integer, intent(in) :: Kmax, N
      real, intent(in) :: Ek(0:Kmax)
      integer, intent(out) :: K_out
      real, intent(out) :: Emiss_out

      real :: Etot, Emiss, Enext, gain, bytes_per_mode
      integer :: K, Knext, K_mem_cap, K_bound

      K_bound = min(Kmax, N/2)

      if (K_bound >= 1) then
         Etot = sum(Ek(1:K_bound))
      else
         Etot = 0.0d0
      end if

      ! User override: fixed K
      if (fourier_kmax > 0) then
         K_out = min(fourier_kmax, K_bound)
         if (Etot > 0.0d0) then
            Emiss_out = max(0.0d0, 1.0d0 - sum(Ek(1:K_out))/Etot)
         else
            Emiss_out = 0.0d0
         end if
         if (nid == 0) write (6, "('Fourier orbit: override fourier_kmax =',I5)") K_out
         return
      end if

      ! Check memory cap override
      K_mem_cap = K_bound
      if (fourier_mem_gb > 0.0d0) then
         bytes_per_mode = 2.0d0*8.0d0*real(ndim*lv + nof*lt)
         if (bytes_per_mode > 0.0d0) then
            K_mem_cap = int((fourier_mem_gb*1.0d9 - 8.0d0*real(ndim*lv + nof*lt))/bytes_per_mode)
            K_mem_cap = max(1, min(K_mem_cap, K_bound))
         end if
         if (nid == 0) write (6, "('Fourier orbit: fourier_mem_gb cap K <=',I5)") K_mem_cap
      end if

      if (Etot <= 0.0d0) then
         K_out = min(8, K_bound)
         Emiss_out = 0.0d0
         return
      end if

      ! Start doubling rule: K = 8 (or min(8, K_bound))
      K = min(8, K_bound)
      Emiss = max(0.0d0, 1.0d0 - sum(Ek(1:K))/Etot)
      if (nid == 0) then
         write (6, "('Fourier orbit doubling: K =',I4,' missing energy =',1PE12.4)") K, Emiss
      end if

      do while (K < K_bound)
         if (Emiss < 1.0d-10) exit
         Knext = min(2*K, K_bound)
         if (Knext == K) exit
         Enext = max(0.0d0, 1.0d0 - sum(Ek(1:Knext))/Etot)
         if (Enext > 0.0d0) then
            gain = Emiss/Enext
         else
            gain = 1.0d30
         end if
         if (nid == 0) then
            write (6, "('Fourier orbit doubling: K =',I4,' -> ',I4,' missing energy =',1PE12.4,' gain =',1PE12.4)") &
               K, Knext, Enext, gain
         end if
         if (gain < 10.0d0 .and. Enext >= 1.0d-10) exit
         K = Knext
         Emiss = Enext
         if (K >= K_mem_cap) then
            K = K_mem_cap
            Emiss = max(0.0d0, 1.0d0 - sum(Ek(1:K))/Etot)
            exit
         end if
      end do

      if (K > K_mem_cap) then
         K = K_mem_cap
         Emiss = max(0.0d0, 1.0d0 - sum(Ek(1:K))/Etot)
      end if

      K_out = K
      Emiss_out = Emiss
      if (nid == 0) then
         write (6, "('Fourier orbit selected K =',I4,' final missing energy =',1PE12.4)") K_out, Emiss_out
      end if
   end subroutine select_K_doubling

   !-----------------------------------------------------------------------
   ! fourier_orbit_init_streaming — start streaming DFT accumulation
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_init_streaming(nsteps_in, T_in)
      integer, intent(in) :: nsteps_in
      real, intent(in) :: T_in

      call fourier_orbit_cleanup()

      stream_N = nsteps_in
      stream_T = T_in
      stream_Kcap = min(128, stream_N/2)
      if (fourier_kmax > 0) stream_Kcap = min(fourier_kmax, stream_N/2)

      allocate (fou_u(lv, 0:2*stream_Kcap), fou_v(lv, 0:2*stream_Kcap))
      call rzero(fou_u, lv*(2*stream_Kcap + 1))
      call rzero(fou_v, lv*(2*stream_Kcap + 1))

      if (if3d) then
         allocate (fou_w(lv, 0:2*stream_Kcap))
         call rzero(fou_w, lv*(2*stream_Kcap + 1))
      else
         allocate (fou_w(1, 0:2*stream_Kcap))
         call rzero(fou_w, 2*stream_Kcap + 1)
      end if

      if (ifto .or. ldimt > 1) then
         allocate (fou_t(lt, ldimt, 0:2*stream_Kcap))
         call rzero(fou_t, lt*ldimt*(2*stream_Kcap + 1))
      else
         allocate (fou_t(1, 1, 0:2*stream_Kcap))
         call rzero(fou_t, 2*stream_Kcap + 1)
      end if

      streaming_active = .true.
      stream_step = 0
   end subroutine fourier_orbit_init_streaming

   !-----------------------------------------------------------------------
   ! fourier_orbit_stream_step — accumulate one step of streaming DFT
   !   c_k += U_n * exp(-2 pi i k n / N) / N
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_stream_step(istep_in)
      integer, intent(in) :: istep_in
      real(c_double) :: theta, C, S, factor0
      integer :: k, nv, nt, m

      if (.not. streaming_active) return

      stream_step = istep_in
      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

      ! DC component (k=0)
      factor0 = 1.0d0/dble(stream_N)
      call add2s2(fou_u(1, 0), vx, real(factor0), nv)
      call add2s2(fou_v(1, 0), vy, real(factor0), nv)
      if (if3d) call add2s2(fou_w(1, 0), vz, real(factor0), nv)
      if (ifto) call add2s2(fou_t(1, 1, 0), t(1, 1, 1, 1, 1), real(factor0), nt)
      if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) then
               call add2s2(fou_t(1, m, 0), t(1, 1, 1, 1, m), real(factor0), nt)
            end if
         end do
      end if

      ! Harmonic components (k=1..Kcap)
      do k = 1, stream_Kcap
         theta = TWOPI*dble(k)*dble(istep_in)/dble(stream_N)
         C = cos(theta)/dble(stream_N)
         S = -sin(theta)/dble(stream_N)

         call add2s2(fou_u(1, 2*k - 1), vx, real(C), nv)
         call add2s2(fou_u(1, 2*k),     vx, real(S), nv)

         call add2s2(fou_v(1, 2*k - 1), vy, real(C), nv)
         call add2s2(fou_v(1, 2*k),     vy, real(S), nv)

         if (if3d) then
            call add2s2(fou_w(1, 2*k - 1), vz, real(C), nv)
            call add2s2(fou_w(1, 2*k),     vz, real(S), nv)
         end if

         if (ifto) then
            call add2s2(fou_t(1, 1, 2*k - 1), t(1, 1, 1, 1, 1), real(C), nt)
            call add2s2(fou_t(1, 1, 2*k),     t(1, 1, 1, 1, 1), real(S), nt)
         end if

         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) then
                  call add2s2(fou_t(1, m, 2*k - 1), t(1, 1, 1, 1, m), real(C), nt)
                  call add2s2(fou_t(1, m, 2*k),     t(1, 1, 1, 1, m), real(S), nt)
               end if
            end do
         end if
      end do
   end subroutine fourier_orbit_stream_step

   !-----------------------------------------------------------------------
   ! fourier_orbit_finalize_streaming — complete streaming DFT, select K
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_finalize_streaming()
      real, allocatable :: Ek(:)
      real :: Ea, Eb
      integer :: k, final_K

      if (.not. streaming_active) return

      allocate (Ek(0:stream_Kcap))

      call compute_mode_energy(Ek(0), fou_u(:, 0), fou_v(:, 0), fou_w(:, 0), fou_t(:, :, 0))

      do k = 1, stream_Kcap
         call compute_mode_energy(Ea, fou_u(:, 2*k - 1), fou_v(:, 2*k - 1), fou_w(:, 2*k - 1), fou_t(:, :, 2*k - 1))
         call compute_mode_energy(Eb, fou_u(:, 2*k), fou_v(:, 2*k), fou_w(:, 2*k), fou_t(:, :, 2*k))
         Ek(k) = Ea + Eb
      end do

      call select_K_doubling(Ek, stream_Kcap, stream_N, final_K, fourier_missing_energy)

      fourier_K = final_K
      call set_replay_depth(stream_Kcap)
      allocate (fourier_Ek(0:fourier_K))
      fourier_Ek(0:fourier_K) = Ek(0:fourier_K)
      deallocate (Ek)

      fourier_T = stream_T
      fourier_N = stream_N
      fourier_dt = stream_T/real(stream_N)
      streaming_active = .false.
      fourier_orbit_loaded = .true.
   end subroutine fourier_orbit_finalize_streaming

   !-----------------------------------------------------------------------
   ! fourier_orbit_build_from_stored — build modes via temporal FFTW
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_build_from_stored(nsteps_in, T_in, uor, vor, wor, tor)
      integer, intent(in) :: nsteps_in
      real, intent(in) :: T_in
      real, dimension(lv, nsteps_in), intent(in) :: uor, vor, wor
      real, dimension(lt, nsteps_in, ldimt), intent(in) :: tor

      real(c_double), allocatable :: signal(:)
      complex(c_double_complex), allocatable :: spectrum(:)
      real, allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_w(:, :)
      real, allocatable :: tmp_t(:, :, :)
      real, allocatable :: Ek(:)
      real :: Ea, Eb, Emiss_file
      integer :: i, k, m, N, Kmax, K_file, nv, nt

      N = nsteps_in
      Kmax = N/2
      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

      allocate (tmp_u(lv, 0:2*Kmax), tmp_v(lv, 0:2*Kmax))
      if (if3d) then
         allocate (tmp_w(lv, 0:2*Kmax))
      else
         allocate (tmp_w(1, 0:2*Kmax))
         tmp_w = 0.0d0
      end if

      if (ifto .or. ldimt > 1) then
         allocate (tmp_t(lt, ldimt, 0:2*Kmax))
         tmp_t = 0.0d0
      else
         allocate (tmp_t(1, 1, 0:2*Kmax))
         tmp_t = 0.0d0
      end if

      allocate (signal(N), spectrum(N/2 + 1))
      call fft_init(N)

      ! FFT in time for each velocity grid point.
      ! Samples are U(T), U(dt), ..., U((N-1)*dt): one period, U(T) standing
      ! in for U(0). That is exact only when the orbit closes. If it does not,
      ! the periodic extension has a jump and the tail falls as 1/k^2.
      do i = 1, nv
         ! u component
         signal(1) = dble(uor(i, N))
         signal(2:N) = dble(uor(i, 1:N - 1))
         call fft_r2c(N, signal, spectrum)
         tmp_u(i, 0) = real(real(spectrum(1), c_double)/dble(N))
         do k = 1, Kmax
            tmp_u(i, 2*k - 1) = real(real(spectrum(k + 1), c_double)/dble(N))
            tmp_u(i, 2*k)     = real(aimag(spectrum(k + 1))/dble(N))
         end do

         ! v component
         signal(1) = dble(vor(i, N))
         signal(2:N) = dble(vor(i, 1:N - 1))
         call fft_r2c(N, signal, spectrum)
         tmp_v(i, 0) = real(real(spectrum(1), c_double)/dble(N))
         do k = 1, Kmax
            tmp_v(i, 2*k - 1) = real(real(spectrum(k + 1), c_double)/dble(N))
            tmp_v(i, 2*k)     = real(aimag(spectrum(k + 1))/dble(N))
         end do

         ! w component (3D)
         if (if3d) then
            signal(1) = dble(wor(i, N))
            signal(2:N) = dble(wor(i, 1:N - 1))
            call fft_r2c(N, signal, spectrum)
            tmp_w(i, 0) = real(real(spectrum(1), c_double)/dble(N))
            do k = 1, Kmax
               tmp_w(i, 2*k - 1) = real(real(spectrum(k + 1), c_double)/dble(N))
               tmp_w(i, 2*k)     = real(aimag(spectrum(k + 1))/dble(N))
            end do
         end if
      end do

      ! Passive scalars / temperature
      if (ifto .or. ldimt > 1) then
         do m = 1, ldimt
            if ((m == 1 .and. ifto) .or. (m > 1 .and. ifpsco(m - 1))) then
               do i = 1, nt
                  signal(1) = dble(tor(i, N, m))
                  signal(2:N) = dble(tor(i, 1:N - 1, m))
                  call fft_r2c(N, signal, spectrum)
                  tmp_t(i, m, 0) = real(real(spectrum(1), c_double)/dble(N))
                  do k = 1, Kmax
                     tmp_t(i, m, 2*k - 1) = real(real(spectrum(k + 1), c_double)/dble(N))
                     tmp_t(i, m, 2*k)     = real(aimag(spectrum(k + 1))/dble(N))
                  end do
               end do
            end if
         end do
      end if

      deallocate (signal, spectrum)

      ! Mode energy computation
      allocate (Ek(0:Kmax))
      call compute_mode_energy(Ek(0), tmp_u(:, 0), tmp_v(:, 0), tmp_w(:, 0), tmp_t(:, :, 0))

      do k = 1, Kmax
         call compute_mode_energy(Ea, tmp_u(:, 2*k - 1), tmp_v(:, 2*k - 1), tmp_w(:, 2*k - 1), tmp_t(:, :, 2*k - 1))
         call compute_mode_energy(Eb, tmp_u(:, 2*k), tmp_v(:, 2*k), tmp_w(:, 2*k), tmp_t(:, :, 2*k))
         Ek(k) = Ea + Eb
      end do

      ! Doubling rule selects the file K. cleanup zeros module state, so the
      ! selected K is held in locals and restored afterwards. Replay keeps
      ! every computed mode when that K still leaves a jump tail: K = N/2
      ! interpolates the samples, K = 8 does not.
      call select_K_doubling(Ek, Kmax, N, K_file, Emiss_file)

      call fourier_orbit_cleanup()

      fourier_K = K_file
      fourier_missing_energy = Emiss_file
      call set_replay_depth(Kmax)

      allocate (fou_u(lv, 0:2*fourier_replay_K), fou_v(lv, 0:2*fourier_replay_K))
      fou_u(:, 0:2*fourier_replay_K) = tmp_u(:, 0:2*fourier_replay_K)
      fou_v(:, 0:2*fourier_replay_K) = tmp_v(:, 0:2*fourier_replay_K)

      if (if3d) then
         allocate (fou_w(lv, 0:2*fourier_replay_K))
         fou_w(:, 0:2*fourier_replay_K) = tmp_w(:, 0:2*fourier_replay_K)
      else
         allocate (fou_w(1, 0:2*fourier_replay_K))
         fou_w = 0.0d0
      end if

      if (ifto .or. ldimt > 1) then
         allocate (fou_t(lt, ldimt, 0:2*fourier_replay_K))
         fou_t(:, :, 0:2*fourier_replay_K) = tmp_t(:, :, 0:2*fourier_replay_K)
      else
         allocate (fou_t(1, 1, 0:2*fourier_replay_K))
         fou_t = 0.0d0
      end if

      allocate (fourier_Ek(0:fourier_K))
      fourier_Ek(0:fourier_K) = Ek(0:fourier_K)

      deallocate (tmp_u, tmp_v, tmp_w, tmp_t, Ek)

      fourier_T = T_in
      fourier_N = N
      fourier_dt = T_in/real(N)
      fourier_orbit_loaded = .true.
   end subroutine fourier_orbit_build_from_stored

   !-----------------------------------------------------------------------
   ! fourier_orbit_eval — evaluate U(t) from Fourier series at time t_in
   !   U(t) = c_0 + 2 Re sum_{k=1..K} c_k exp(2 pi i k t / T)
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_eval(t_in, vx_out, vy_out, vz_out, t_out)
      real, intent(in) :: t_in
      real, dimension(lv), intent(out) :: vx_out, vy_out, vz_out
      real, dimension(lt, ldimt), intent(out) :: t_out

      real(c_double) :: t_phase, theta, C, S, fac
      integer :: k, nv, nt, m

      if (.not. fourier_orbit_loaded) return

      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

      ! Periodic phase in [0, T)
      if (fourier_T > 0.0d0) then
         t_phase = mod(dble(t_in), dble(fourier_T))
         if (t_phase < 0.0d0) t_phase = t_phase + dble(fourier_T)
      else
         t_phase = 0.0d0
      end if

      ! DC component (k=0)
      vx_out(1:nv) = fou_u(1:nv, 0)
      vy_out(1:nv) = fou_v(1:nv, 0)
      if (if3d) vz_out(1:nv) = fou_w(1:nv, 0)
      if (ifto) t_out(1:nt, 1) = fou_t(1:nt, 1, 0)
      if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) t_out(1:nt, m) = fou_t(1:nt, m, 0)
         end do
      end if

      ! Harmonic reconstruction (k=1..K)
      do k = 1, fourier_replay_K
         theta = TWOPI*dble(k)*t_phase/dble(fourier_T)

         ! Nyquist mode at k = N/2 for even N has factor 1.0 (no conjugate pair)
         if (k == fourier_N/2 .and. mod(fourier_N, 2) == 0) then
            fac = 1.0d0
         else
            fac = 2.0d0
         end if

         C = fac*cos(theta)
         S = -fac*sin(theta)

         vx_out(1:nv) = vx_out(1:nv) + real(C)*fou_u(1:nv, 2*k - 1) &
                                     + real(S)*fou_u(1:nv, 2*k)

         vy_out(1:nv) = vy_out(1:nv) + real(C)*fou_v(1:nv, 2*k - 1) &
                                     + real(S)*fou_v(1:nv, 2*k)

         if (if3d) then
            vz_out(1:nv) = vz_out(1:nv) + real(C)*fou_w(1:nv, 2*k - 1) &
                                         + real(S)*fou_w(1:nv, 2*k)
         end if

         if (ifto) then
            t_out(1:nt, 1) = t_out(1:nt, 1) + real(C)*fou_t(1:nt, 1, 2*k - 1) &
                                            + real(S)*fou_t(1:nt, 1, 2*k)
         end if

         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) then
                  t_out(1:nt, m) = t_out(1:nt, m) + real(C)*fou_t(1:nt, m, 2*k - 1) &
                                                  + real(S)*fou_t(1:nt, m, 2*k)
               end if
            end do
         end if
      end do
   end subroutine fourier_orbit_eval

   !-----------------------------------------------------------------------
   ! fourier_orbit_write — write modes to FOU<session>0.f* and <session>.fourier
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_write()
      integer :: io_unit, k, m
      real :: p63_sav

      if (.not. fourier_orbit_loaded) return

      ! 1. Rank 0 writes text header <session>.fourier
      if (nid == 0) then
         open (newunit=io_unit, file=trim(SESSION)//'.fourier', status='replace', action='write')
         write (io_unit, "(A,1PE22.14)") "period: ", fourier_T
         write (io_unit, "(A,I10)")      "N:      ", fourier_N
         write (io_unit, "(A,1PE22.14)") "dt:     ", fourier_dt
         write (io_unit, "(A,I10)")      "K:      ", fourier_K
         write (io_unit, "(A,1PE22.14)") "missing_energy: ", fourier_missing_energy
         write (io_unit, "(A)")          "# mode  energy"
         do k = 0, fourier_K
            write (io_unit, "(I6,1PE22.14)") k, fourier_Ek(k)
         end do
         close (io_unit)
         write (6, *) 'Wrote Fourier orbit header: ', trim(SESSION)//'.fourier'
      end if

      ! 2. All ranks output FOU<session>0.f00001 ... FOU<session>0.f(2K+1)
      p63_sav = param(63)
      param(63) = 1.0d0
      call bcast(param(63), wdsize)

      do m = 0, 2*fourier_K
         call outpost2(fou_u(1, m), fou_v(1, m), fou_w(1, m), pr, fou_t(1, 1, m), nof, 'FOU')
      end do

      param(63) = p63_sav
      call bcast(param(63), wdsize)

      if (nid == 0) then
         write (6, "('Wrote ',I4,' Fourier mode checkpoint files (FOU',A,'0.f*)')") &
            2*fourier_K + 1, trim(SESSION)
      end if
   end subroutine fourier_orbit_write

   !-----------------------------------------------------------------------
   ! fourier_orbit_load — load modes from FOU<session>0.f* and <session>.fourier
   !-----------------------------------------------------------------------
   subroutine fourier_orbit_load(loaded)
      logical, intent(out) :: loaded

      logical :: f_exist1, f_exist2, f_exist
      integer :: io_unit, ierr, k, idx, m, mm, col
      character(len=132) :: line, fname
      real :: header_period, header_dt, header_missing_energy
      integer :: header_N, header_K
      real, allocatable :: header_Ek(:)
      real :: cur_dt, time_sav
      integer :: nv, nt

      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

      ! Check file existence
      f_exist = .false.
      if (nid == 0) then
         inquire (file='FOU'//trim(SESSION)//'0.f00001', exist=f_exist1)
         inquire (file=trim(SESSION)//'.fourier', exist=f_exist2)
         f_exist = f_exist1 .and. f_exist2
      end if
      call bcast(f_exist, lsize)

      if (.not. f_exist) then
         loaded = .false.
         return
      end if

      ! Rank 0 reads header
      if (nid == 0) then
         open (newunit=io_unit, file=trim(SESSION)//'.fourier', status='old', action='read', iostat=ierr)
         if (ierr /= 0) then
            call nekStab_error('could not open '//trim(SESSION)//'.fourier')
         end if

         read (io_unit, "(A)") line; col = index(line, ':'); read (line(col + 1:), *) header_period
         read (io_unit, "(A)") line; col = index(line, ':'); read (line(col + 1:), *) header_N
         read (io_unit, "(A)") line; col = index(line, ':'); read (line(col + 1:), *) header_dt
         read (io_unit, "(A)") line; col = index(line, ':'); read (line(col + 1:), *) header_K
         read (io_unit, "(A)") line; col = index(line, ':'); read (line(col + 1:), *) header_missing_energy
         read (io_unit, "(A)") line ! comment line

         allocate (header_Ek(0:header_K))
         do k = 0, header_K
            read (io_unit, *) idx, header_Ek(k)
         end do
         close (io_unit)
      end if

      ! Broadcast header metadata
      call bcast(header_period, wdsize)
      call bcast(header_N, isize)
      call bcast(header_dt, wdsize)
      call bcast(header_K, isize)
      call bcast(header_missing_energy, wdsize)

      call fourier_orbit_cleanup()

      allocate (fourier_Ek(0:header_K))
      if (nid == 0) fourier_Ek = header_Ek
      call bcast(fourier_Ek, (header_K + 1)*wdsize)
      if (nid == 0) deallocate (header_Ek)

      ! Warn if current dt differs from header dt
      cur_dt = abs(param(12))
      if (cur_dt > 0.0d0 .and. header_dt > 0.0d0) then
         if (abs(cur_dt - header_dt) > 1.0d-7*header_dt) then
            if (nid == 0) then
               write (6, *) 'nekStab: NOTE: current dt differs from Fourier header dt:'
               write (6, *) '  current dt =', cur_dt, ' header dt =', header_dt
               write (6, *) '  Reconstructed orbit is smooth; neutral multiplier will move by time discretisation error.'
            end if
         end if
      end if

      ! Allocate module mode arrays
      allocate (fou_u(lv, 0:2*header_K), fou_v(lv, 0:2*header_K))
      if (if3d) then
         allocate (fou_w(lv, 0:2*header_K))
      else
         allocate (fou_w(1, 0:2*header_K))
         fou_w = 0.0d0
      end if

      if (ifto .or. ldimt > 1) then
         allocate (fou_t(lt, ldimt, 0:2*header_K))
         fou_t = 0.0d0
      else
         allocate (fou_t(1, 1, 0:2*header_K))
         fou_t = 0.0d0
      end if

      ! Load each mode file
      time_sav = param(10)
      do idx = 1, 2*header_K + 1
         m = idx - 1
         write (fname, '(A,A,"0.f",I5.5)') 'FOU', trim(SESSION), idx
         call load_fld(fname)
         call copy(fou_u(1, m), vx, nv)
         call copy(fou_v(1, m), vy, nv)
         if (if3d) call copy(fou_w(1, m), vz, nv)
         if (ifto) call copy(fou_t(1, 1, m), t(1, 1, 1, 1, 1), nt)
         if (ldimt > 1) then
            do mm = 2, ldimt
               if (ifpsco(mm - 1)) call copy(fou_t(1, mm, m), t(1, 1, 1, 1, mm), nt)
            end do
         end if
      end do
      param(10) = time_sav

      fourier_T = header_period
      fourier_N = header_N
      fourier_dt = header_dt
      fourier_K = header_K
      fourier_replay_K = header_K
      fourier_missing_energy = header_missing_energy
      fourier_orbit_loaded = .true.
      loaded = .true.

      if (nid == 0) then
         write (6, "('Fourier orbit loaded: K =',I4,' N =',I6,' T =',1PE15.7,' missing energy =',1PE12.4)") &
            fourier_K, fourier_N, fourier_T, fourier_missing_energy
      end if
   end subroutine fourier_orbit_load

end module nekstab_fourier
