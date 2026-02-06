      !-----------------------------------------------------------------------
      ! fourier_fftw.f90 — Low-level FFT module using FFTW3 interface
      !
      ! Purpose:
      !   Provides temporal FFT decomposition/reconstruction for spectral
      !   analysis. Compatible with both FFTW3 (GCC) and MKL (Intel)
      !   via the FFTW3 wrapper interface.
      !
      ! Public interface:
      !   fft_init, fft_cleanup       — plan management
      !   fft_r2c, fft_c2r            — forward/inverse FFT
      !   fft_frequencies              — compute frequency array
      !   fourier_decomposition        — multi-point temporal FFT
      !   fourier_reconstruction       — inverse from Fourier modes
      !
      ! Dependencies:
      !   iso_c_binding (FFTW3 C bindings)
      !-----------------------------------------------------------------------
module fourier_fftw
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! Public procedures
  public :: fft_init, fft_cleanup
  public :: fft_r2c, fft_c2r
  public :: fft_frequencies
  public :: fourier_decomposition, fourier_reconstruction

  ! FFTW3 constants (from fftw3.f03)
  integer(C_INT), parameter :: FFTW_ESTIMATE = 64
  integer(C_INT), parameter :: FFTW_MEASURE = 0
  integer(C_INT), parameter :: FFTW_PATIENT = 32
  integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8

  ! FFTW3 C interface bindings
  interface
    type(C_PTR) function fftw_plan_dft_r2c_1d(n, in, out, flags) bind(C, name='fftw_plan_dft_r2c_1d')
      import :: C_PTR, C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
      integer(C_INT), value :: n
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function

    type(C_PTR) function fftw_plan_dft_c2r_1d(n, in, out, flags) bind(C, name='fftw_plan_dft_c2r_1d')
      import :: C_PTR, C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
      integer(C_INT), value :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function

    subroutine fftw_execute_dft_r2c(plan, in, out) bind(C, name='fftw_execute_dft_r2c')
      import :: C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX
      type(C_PTR), value :: plan
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine

    subroutine fftw_execute_dft_c2r(plan, in, out) bind(C, name='fftw_execute_dft_c2r')
      import :: C_PTR, C_DOUBLE, C_DOUBLE_COMPLEX
      type(C_PTR), value :: plan
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine

    subroutine fftw_destroy_plan(plan) bind(C, name='fftw_destroy_plan')
      import :: C_PTR
      type(C_PTR), value :: plan
    end subroutine
  end interface

  ! Module state
  integer, save :: plan_size = 0
  type(C_PTR), save :: plan_r2c = C_NULL_PTR
  type(C_PTR), save :: plan_c2r = C_NULL_PTR
  logical, save :: initialized = .false.

  ! Mathematical constants
  real(C_DOUBLE), parameter :: PI = 3.14159265358979323846d0
  real(C_DOUBLE), parameter :: TWOPI = 2.0d0 * PI

contains

  !-----------------------------------------------------------------------
  ! fft_init — initialize FFTW plans for a given transform size
  !
  ! Purpose:
  !   Create forward/inverse FFT plans. Plans are cached; calling
  !   with the same n is a no-op. Calling with different n destroys
  !   old plans and creates new ones.
  !
  ! Arguments:
  !   n [in] — number of real samples (time points)
  !-----------------------------------------------------------------------
  subroutine fft_init(n)
    integer, intent(in) :: n
    real(C_DOUBLE), allocatable :: tmp_in(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: tmp_out(:)
    integer(C_INT) :: nc

    ! Skip if already initialized with same size
    if (initialized .and. plan_size == n) return

    ! Cleanup old plans if reinitializing
    if (initialized) call fft_cleanup()

    ! Allocate temporary arrays for planning
    nc = n / 2 + 1
    allocate(tmp_in(n), tmp_out(nc))

    ! Create plans (FFTW_ESTIMATE for fast planning)
    plan_r2c = fftw_plan_dft_r2c_1d(int(n, C_INT), tmp_in, tmp_out, FFTW_ESTIMATE)
    plan_c2r = fftw_plan_dft_c2r_1d(int(n, C_INT), tmp_out, tmp_in, FFTW_ESTIMATE)

    deallocate(tmp_in, tmp_out)

    plan_size = n
    initialized = .true.
  end subroutine fft_init

  !-----------------------------------------------------------------------
  ! fft_cleanup — destroy FFTW plans and free resources
  !-----------------------------------------------------------------------
  subroutine fft_cleanup()
    if (c_associated(plan_r2c)) call fftw_destroy_plan(plan_r2c)
    if (c_associated(plan_c2r)) call fftw_destroy_plan(plan_c2r)
    plan_r2c = C_NULL_PTR
    plan_c2r = C_NULL_PTR
    plan_size = 0
    initialized = .false.
  end subroutine fft_cleanup

  !-----------------------------------------------------------------------
  ! fft_r2c — real-to-complex forward FFT
  !
  ! Purpose:
  !   Compute forward FFT. Output is NOT normalized (multiply by 1/n
  !   for standard normalization). Only positive frequencies stored.
  !
  ! Arguments:
  !   n      [in]    — number of real samples
  !   input  [inout] — real input array (n elements)
  !   output [out]   — complex output array (n/2+1 elements)
  !-----------------------------------------------------------------------
  subroutine fft_r2c(n, input, output)
    integer, intent(in) :: n
    real(C_DOUBLE), intent(inout) :: input(n)
    complex(C_DOUBLE_COMPLEX), intent(out) :: output(n/2+1)

    if (.not. initialized .or. plan_size /= n) call fft_init(n)
    call fftw_execute_dft_r2c(plan_r2c, input, output)
  end subroutine fft_r2c

  !-----------------------------------------------------------------------
  ! fft_c2r — complex-to-real inverse FFT
  !
  ! Purpose:
  !   Compute inverse FFT. Input array IS MODIFIED by FFTW (c2r
  !   destroys input). Output is scaled by n.
  !
  ! Arguments:
  !   n      [in]    — number of real samples (output size)
  !   input  [inout] — complex input array (n/2+1 elements, modified!)
  !   output [out]   — real output array (n elements)
  !-----------------------------------------------------------------------
  subroutine fft_c2r(n, input, output)
    integer, intent(in) :: n
    complex(C_DOUBLE_COMPLEX), intent(inout) :: input(n/2+1)
    real(C_DOUBLE), intent(out) :: output(n)

    if (.not. initialized .or. plan_size /= n) call fft_init(n)
    call fftw_execute_dft_c2r(plan_c2r, input, output)
  end subroutine fft_c2r

  !-----------------------------------------------------------------------
  ! fft_frequencies — compute frequency array for FFT output
  !
  ! Purpose:
  !   Returns frequencies in Hz: f(k) = k / (n * dt) for k = 0..n/2.
  !   For angular frequency (rad/s): omega = 2*pi*freq.
  !
  ! Arguments:
  !   n    [in]  — number of samples
  !   dt   [in]  — time step (sampling interval)
  !   freq [out] — frequency array (n/2+1 elements)
  !-----------------------------------------------------------------------
  subroutine fft_frequencies(n, dt, freq)
    integer, intent(in) :: n
    real(C_DOUBLE), intent(in) :: dt
    real(C_DOUBLE), intent(out) :: freq(n/2+1)
    integer :: k
    real(C_DOUBLE) :: df

    df = 1.0d0 / (n * dt)
    do k = 0, n/2
      freq(k+1) = k * df
    end do
  end subroutine fft_frequencies

  !-----------------------------------------------------------------------
  ! fourier_decomposition — decompose velocity into Fourier modes
  !
  ! Purpose:
  !   Performs temporal FFT at each spatial point. Returns normalized
  !   complex coefficients and frequency array.
  !
  ! Arguments:
  !   npts     [in]    — number of spatial points (lv)
  !   nsnap    [in]    — number of time snapshots
  !   vx,vy,vz [inout] — velocity components (npts x nsnap)
  !   time     [in]    — time array (nsnap)
  !   bm1      [in]    — mass matrix for energy weighting (npts)
  !   nfreq    [out]   — number of positive frequencies (nsnap/2+1)
  !   freq     [out]   — frequency array (nfreq)
  !   vx_hat   [out]   — FFT of vx (npts x nfreq), complex
  !   vy_hat   [out]   — FFT of vy (npts x nfreq), complex
  !   vz_hat   [out]   — FFT of vz (npts x nfreq), complex
  !-----------------------------------------------------------------------
  subroutine fourier_decomposition(npts, nsnap, vx, vy, vz, time, bm1, &
                                    nfreq, freq, vx_hat, vy_hat, vz_hat)
    integer, intent(in) :: npts, nsnap
    real(C_DOUBLE), intent(inout) :: vx(npts, nsnap)
    real(C_DOUBLE), intent(inout) :: vy(npts, nsnap)
    real(C_DOUBLE), intent(inout) :: vz(npts, nsnap)
    real(C_DOUBLE), intent(in) :: time(nsnap)
    real(C_DOUBLE), intent(in) :: bm1(npts)
    integer, intent(out) :: nfreq
    real(C_DOUBLE), intent(out) :: freq(nsnap/2+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: vx_hat(npts, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: vy_hat(npts, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: vz_hat(npts, nsnap/2+1)

    integer :: i, k
    real(C_DOUBLE) :: dt, norm_factor
    real(C_DOUBLE) :: signal(nsnap)
    complex(C_DOUBLE_COMPLEX) :: spectrum(nsnap/2+1)

    ! Initialize FFT
    call fft_init(nsnap)
    nfreq = nsnap / 2 + 1

    ! Compute sampling interval (assuming uniform sampling)
    dt = (time(nsnap) - time(1)) / (nsnap - 1)

    ! Compute frequency array
    call fft_frequencies(nsnap, dt, freq)

    ! Normalization factor for standard FFT normalization
    norm_factor = 1.0d0 / nsnap

    ! FFT each spatial point
    do i = 1, npts
      ! vx component
      signal(:) = vx(i, :)
      call fft_r2c(nsnap, signal, spectrum)
      vx_hat(i, :) = spectrum * norm_factor

      ! vy component
      signal(:) = vy(i, :)
      call fft_r2c(nsnap, signal, spectrum)
      vy_hat(i, :) = spectrum * norm_factor

      ! vz component
      signal(:) = vz(i, :)
      call fft_r2c(nsnap, signal, spectrum)
      vz_hat(i, :) = spectrum * norm_factor
    end do
  end subroutine fourier_decomposition

  !-----------------------------------------------------------------------
  ! fourier_reconstruction — reconstruct velocity at time t from modes
  !
  ! Purpose:
  !   Inverse Fourier synthesis from complex coefficients at a single
  !   time instant. Accounts for Hermitian symmetry (factor 2) and
  !   Nyquist frequency for even-length transforms.
  !
  ! Arguments:
  !   npts     [in]  — number of spatial points
  !   nfreq    [in]  — number of frequencies
  !   nsnap    [in]  — original number of snapshots (even/odd check)
  !   freq     [in]  — frequency array (Hz)
  !   vx_hat   [in]  — FFT coefficients for vx (npts x nfreq)
  !   vy_hat   [in]  — FFT coefficients for vy (npts x nfreq)
  !   vz_hat   [in]  — FFT coefficients for vz (npts x nfreq)
  !   t        [in]  — time at which to reconstruct
  !   vx,vy,vz [out] — reconstructed velocity (npts)
  !-----------------------------------------------------------------------
  subroutine fourier_reconstruction(npts, nfreq, nsnap, freq, vx_hat, vy_hat, vz_hat, &
                                     t, vx, vy, vz)
    integer, intent(in) :: npts, nfreq, nsnap
    real(C_DOUBLE), intent(in) :: freq(nfreq)
    complex(C_DOUBLE_COMPLEX), intent(in) :: vx_hat(npts, nfreq)
    complex(C_DOUBLE_COMPLEX), intent(in) :: vy_hat(npts, nfreq)
    complex(C_DOUBLE_COMPLEX), intent(in) :: vz_hat(npts, nfreq)
    real(C_DOUBLE), intent(in) :: t
    real(C_DOUBLE), intent(out) :: vx(npts), vy(npts), vz(npts)

    integer :: i, k
    real(C_DOUBLE) :: omega, cos_wt, sin_wt, factor
    logical :: has_nyquist

    ! Nyquist frequency only exists for even-length FFT
    has_nyquist = (mod(nsnap, 2) == 0)

    vx = 0.0d0
    vy = 0.0d0
    vz = 0.0d0

    ! DC component (k=0)
    vx(:) = real(vx_hat(:, 1))
    vy(:) = real(vy_hat(:, 1))
    vz(:) = real(vz_hat(:, 1))

    ! Positive frequencies (k=1 to nfreq-1)
    ! Factor of 2 accounts for negative frequencies (Hermitian symmetry)
    do k = 2, nfreq
      omega = TWOPI * freq(k)
      cos_wt = cos(omega * t)
      sin_wt = sin(omega * t)

      ! Nyquist frequency (last bin for even nsnap) has no conjugate pair
      if (k == nfreq .and. has_nyquist) then
        factor = 1.0d0
      else
        factor = 2.0d0
      end if

      do i = 1, npts
        vx(i) = vx(i) + factor * (real(vx_hat(i,k)) * cos_wt - aimag(vx_hat(i,k)) * sin_wt)
        vy(i) = vy(i) + factor * (real(vy_hat(i,k)) * cos_wt - aimag(vy_hat(i,k)) * sin_wt)
        vz(i) = vz(i) + factor * (real(vz_hat(i,k)) * cos_wt - aimag(vz_hat(i,k)) * sin_wt)
      end do
    end do
  end subroutine fourier_reconstruction

end module fourier_fftw
