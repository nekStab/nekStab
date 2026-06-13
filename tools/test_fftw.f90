!-----------------------------------------------------------------------
! test_fftw.f90: Comprehensive test suite for fourier_fftw module
!
! Tests:
!   1. Basic FFT initialization and cleanup
!   2. Forward/inverse FFT round-trip accuracy
!   3. Known signal reconstruction (sine, cosine, mixed)
!   4. Parseval's theorem (energy conservation)
!   5. Different FFT sizes (power of 2, non-power of 2)
!   6. Edge cases (DC signal, Nyquist frequency)
!   7. High-level API (fourier_decomposition/reconstruction)
!
! Build:
!   GCC:   gfortran -o test_fftw fourier_fftw.f90 test_fftw.f90 -lfftw3
!   Intel: ifort/ifx -o test_fftw fourier_fftw.f90 test_fftw.f90 -qmkl
!-----------------------------------------------------------------------
program test_fftw
  use fourier_fftw
  use, intrinsic :: iso_c_binding
  implicit none

  integer :: n_passed, n_failed
  real(C_DOUBLE), parameter :: PI = 3.14159265358979323846d0
  real(C_DOUBLE), parameter :: TOL = 1.0d-10

  n_passed = 0
  n_failed = 0

  print '(A)', '========================================================'
  print '(A)', '  FFTW3 Interface Test Suite for nekStab'
  print '(A)', '========================================================'
  print '(A)', ''

  call test_init_cleanup()
  call test_roundtrip_accuracy()
  call test_known_signals()
  call test_parseval_theorem()
  call test_different_sizes()
  call test_edge_cases()
  call test_high_level_api()

  print '(A)', ''
  print '(A)', '========================================================'
  print '(A,I3,A,I3,A)', '  Results: ', n_passed, ' passed, ', n_failed, ' failed'
  print '(A)', '========================================================'

  if (n_failed > 0) stop 1

contains

  !-------------------------------------------------------------------
  subroutine report(test_name, passed)
    character(*), intent(in) :: test_name
    logical, intent(in) :: passed
    if (passed) then
      print '(A,A)', '  [PASS] ', test_name
      n_passed = n_passed + 1
    else
      print '(A,A)', '  [FAIL] ', test_name
      n_failed = n_failed + 1
    end if
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_init_cleanup()
    print '(A)', 'Test 1: Initialization and Cleanup'
    print '(A)', '-----------------------------------'

    ! Test basic init
    call fft_init(64)
    call report('Init with n=64', .true.)

    ! Test re-init with same size (should be no-op)
    call fft_init(64)
    call report('Re-init with same size', .true.)

    ! Test re-init with different size
    call fft_init(128)
    call report('Re-init with n=128', .true.)

    ! Cleanup
    call fft_cleanup()
    call report('Cleanup', .true.)
    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_roundtrip_accuracy()
    integer, parameter :: n = 64
    real(C_DOUBLE) :: signal(n), reconstructed(n), original(n)
    complex(C_DOUBLE_COMPLEX) :: spectrum(n/2+1)
    real(C_DOUBLE) :: max_error
    integer :: i

    print '(A)', 'Test 2: Round-trip Accuracy'
    print '(A)', '---------------------------'

    ! Random-ish signal
    do i = 1, n
      original(i) = sin(0.1d0*i) + 0.5d0*cos(0.3d0*i) + 0.1d0*i
    end do
    signal = original

    call fft_init(n)
    call fft_r2c(n, signal, spectrum)
    call fft_c2r(n, spectrum, reconstructed)
    reconstructed = reconstructed / n  ! Normalize

    max_error = maxval(abs(original - reconstructed))
    call report('Round-trip error < 1e-10: ' // trim(fmt_exp(max_error)), max_error < TOL)

    call fft_cleanup()
    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_known_signals()
    ! Use frequencies that exactly align with FFT bins to avoid spectral leakage
    integer, parameter :: n = 128
    real(C_DOUBLE) :: signal(n), dt, df
    complex(C_DOUBLE_COMPLEX) :: spectrum(n/2+1)
    real(C_DOUBLE) :: freq(n/2+1)
    real(C_DOUBLE) :: amp_f1, amp_f2
    real(C_DOUBLE) :: f1, f2  ! Exact bin frequencies
    integer :: i, idx_f1, idx_f2

    print '(A)', 'Test 3: Known Signal Reconstruction'
    print '(A)', '------------------------------------'

    dt = 0.01d0  ! 100 Hz sampling, 1.28s total
    df = 1.0d0 / (n * dt)  ! Frequency resolution = 0.78125 Hz

    ! Use frequencies that align exactly with FFT bins
    idx_f1 = 8   ! bin 8 -> f1 = 7 * df = 5.46875 Hz
    idx_f2 = 16  ! bin 16 -> f2 = 15 * df = 11.71875 Hz
    f1 = (idx_f1 - 1) * df
    f2 = (idx_f2 - 1) * df

    call fft_init(n)
    call fft_frequencies(n, dt, freq)

    ! Test pure cosine: 3*cos(2*pi*f1*t)
    do i = 1, n
      signal(i) = 3.0d0 * cos(2.0d0 * PI * f1 * (i-1) * dt)
    end do
    call fft_r2c(n, signal, spectrum)
    amp_f1 = 2.0d0 * abs(spectrum(idx_f1)) / n  ! Factor 2 for one-sided
    call report('Cosine amplitude (expect 3.0): ' // trim(fmt_real(amp_f1)), &
                abs(amp_f1 - 3.0d0) < 1.0d-10)

    ! Test pure sine: 2*sin(2*pi*f2*t)
    do i = 1, n
      signal(i) = 2.0d0 * sin(2.0d0 * PI * f2 * (i-1) * dt)
    end do
    call fft_r2c(n, signal, spectrum)
    amp_f2 = 2.0d0 * abs(spectrum(idx_f2)) / n
    call report('Sine amplitude (expect 2.0): ' // trim(fmt_real(amp_f2)), &
                abs(amp_f2 - 2.0d0) < 1.0d-10)

    ! Test mixed signal
    do i = 1, n
      signal(i) = 3.0d0 * cos(2.0d0 * PI * f1 * (i-1) * dt) + &
                  2.0d0 * sin(2.0d0 * PI * f2 * (i-1) * dt)
    end do
    call fft_r2c(n, signal, spectrum)
    amp_f1 = 2.0d0 * abs(spectrum(idx_f1)) / n
    amp_f2 = 2.0d0 * abs(spectrum(idx_f2)) / n
    call report('Mixed signal: f1 component (expect 3.0): ' // trim(fmt_real(amp_f1)), &
                abs(amp_f1 - 3.0d0) < 1.0d-10)
    call report('Mixed signal: f2 component (expect 2.0): ' // trim(fmt_real(amp_f2)), &
                abs(amp_f2 - 2.0d0) < 1.0d-10)

    call fft_cleanup()
    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_parseval_theorem()
    ! Parseval: sum(|x|^2) = (1/N) * sum(|X|^2)
    integer, parameter :: n = 64
    real(C_DOUBLE) :: signal(n), original(n)
    complex(C_DOUBLE_COMPLEX) :: spectrum(n/2+1)
    real(C_DOUBLE) :: energy_time, energy_freq, ratio
    integer :: i

    print '(A)', 'Test 4: Parseval''s Theorem (Energy Conservation)'
    print '(A)', '-------------------------------------------------'

    ! Create test signal
    do i = 1, n
      original(i) = sin(0.2d0*i) + 0.3d0*cos(0.5d0*i)
    end do
    signal = original

    ! Time-domain energy
    energy_time = sum(original**2)

    ! Frequency-domain energy
    call fft_init(n)
    call fft_r2c(n, signal, spectrum)

    ! For real FFT: E = |X(0)|^2 + 2*sum(|X(k)|^2, k=1..N/2-1) + |X(N/2)|^2
    energy_freq = abs(spectrum(1))**2
    do i = 2, n/2
      energy_freq = energy_freq + 2.0d0 * abs(spectrum(i))**2
    end do
    energy_freq = energy_freq + abs(spectrum(n/2+1))**2
    energy_freq = energy_freq / n  ! Normalize

    ratio = energy_freq / energy_time
    call report('Energy ratio (expect 1.0): ' // trim(fmt_real(ratio)), &
                abs(ratio - 1.0d0) < TOL)

    call fft_cleanup()
    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_different_sizes()
    integer :: sizes(6), n, i, j
    real(C_DOUBLE), allocatable :: signal(:), reconstructed(:), original(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: spectrum(:)
    real(C_DOUBLE) :: max_error
    logical :: passed
    character(len=64) :: msg

    print '(A)', 'Test 5: Different FFT Sizes'
    print '(A)', '----------------------------'

    sizes = [32, 64, 128, 256, 100, 144]  ! Mix of power-of-2 and others

    do j = 1, size(sizes)
      n = sizes(j)
      allocate(signal(n), reconstructed(n), original(n), spectrum(n/2+1))

      do i = 1, n
        original(i) = sin(0.1d0*i) + cos(0.2d0*i)
      end do
      signal = original

      call fft_init(n)
      call fft_r2c(n, signal, spectrum)
      call fft_c2r(n, spectrum, reconstructed)
      reconstructed = reconstructed / n

      max_error = maxval(abs(original - reconstructed))
      passed = max_error < TOL
      write(msg, '(A,I4,A,ES9.2)') 'n=', n, ', error=', max_error
      call report(trim(msg), passed)

      call fft_cleanup()
      deallocate(signal, reconstructed, original, spectrum)
    end do

    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_edge_cases()
    integer, parameter :: n = 64
    real(C_DOUBLE) :: signal(n), reconstructed(n)
    complex(C_DOUBLE_COMPLEX) :: spectrum(n/2+1)
    real(C_DOUBLE) :: freq(n/2+1), dt, max_error

    print '(A)', 'Test 6: Edge Cases'
    print '(A)', '------------------'

    dt = 0.01d0
    call fft_init(n)
    call fft_frequencies(n, dt, freq)

    ! DC signal (constant)
    signal = 5.0d0
    call fft_r2c(n, signal, spectrum)
    call report('DC signal: only bin 0 nonzero', &
                abs(real(spectrum(1))/n - 5.0d0) < TOL .and. &
                maxval(abs(spectrum(2:))) < TOL * n)

    ! Nyquist frequency signal
    signal(1:n:2) = 1.0d0
    signal(2:n:2) = -1.0d0
    call fft_r2c(n, signal, spectrum)
    call report('Nyquist signal: only last bin nonzero', &
                abs(spectrum(n/2+1)) > 0.1d0 * n .and. &
                maxval(abs(spectrum(2:n/2))) < TOL * n)

    ! Zero signal
    signal = 0.0d0
    call fft_r2c(n, signal, spectrum)
    call report('Zero signal: all bins zero', maxval(abs(spectrum)) < TOL)

    call fft_cleanup()
    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  subroutine test_high_level_api()
    integer, parameter :: nv = 10, nsnap = 64
    real(C_DOUBLE) :: vx(nv, nsnap), vy(nv, nsnap), vz(nv, nsnap)
    real(C_DOUBLE) :: time(nsnap), bm1(nv)
    real(C_DOUBLE) :: freq(nsnap/2+1)
    complex(C_DOUBLE_COMPLEX) :: vx_hat(nv, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX) :: vy_hat(nv, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX) :: vz_hat(nv, nsnap/2+1)
    real(C_DOUBLE) :: vx_rec(nv), vy_rec(nv), vz_rec(nv)
    real(C_DOUBLE) :: dt, df, t_test, max_error
    real(C_DOUBLE) :: f1, f2, f3  ! Bin-aligned frequencies
    integer :: nfreq, i, j

    print '(A)', 'Test 7: High-Level API (fourier_decomposition/reconstruction)'
    print '(A)', '-------------------------------------------------------------'

    ! Setup with bin-aligned frequencies
    dt = 0.01d0
    df = 1.0d0 / (nsnap * dt)  ! = 1.5625 Hz
    f1 = 5 * df   ! 7.8125 Hz (bin 6)
    f2 = 8 * df   ! 12.5 Hz (bin 9)
    f3 = 3 * df   ! 4.6875 Hz (bin 4)

    do i = 1, nsnap
      time(i) = (i-1) * dt
    end do
    bm1 = 1.0d0  ! Uniform mass matrix

    ! Create velocity field: bin-aligned frequencies at each point
    do i = 1, nv
      do j = 1, nsnap
        vx(i, j) = 2.0d0 * cos(2.0d0 * PI * f1 * time(j) + 0.1d0*i)
        vy(i, j) = 1.5d0 * sin(2.0d0 * PI * f2 * time(j) + 0.2d0*i)
        vz(i, j) = 1.0d0 * cos(2.0d0 * PI * f3 * time(j) + 0.3d0*i)
      end do
    end do

    ! Decompose
    call fourier_decomposition(nv, nsnap, vx, vy, vz, time, bm1, &
                               nfreq, freq, vx_hat, vy_hat, vz_hat)
    call report('Decomposition: nfreq = ' // trim(itoa(nfreq)), nfreq == nsnap/2+1)

    ! Reconstruct at a sample point in time array (exact match)
    t_test = time(17)  ! Use exact sample point
    call fourier_reconstruction(nv, nfreq, nsnap, freq, vx_hat, vy_hat, vz_hat, &
                                 t_test, vx_rec, vy_rec, vz_rec)

    ! Compare with expected values
    max_error = 0.0d0
    do i = 1, nv
      max_error = max(max_error, abs(vx_rec(i) - 2.0d0 * cos(2.0d0*PI*f1*t_test + 0.1d0*i)))
      max_error = max(max_error, abs(vy_rec(i) - 1.5d0 * sin(2.0d0*PI*f2*t_test + 0.2d0*i)))
      max_error = max(max_error, abs(vz_rec(i) - 1.0d0 * cos(2.0d0*PI*f3*t_test + 0.3d0*i)))
    end do
    call report('Reconstruction error: ' // trim(fmt_exp(max_error)), &
                max_error < 1.0d-10)

    call fft_cleanup()
    print '(A)', ''
  end subroutine

  !-------------------------------------------------------------------
  ! Helper functions for formatting
  !-------------------------------------------------------------------
  function fmt_real(x) result(str)
    real(C_DOUBLE), intent(in) :: x
    character(len=16) :: str
    write(str, '(F8.4)') x
    str = adjustl(str)
  end function

  function fmt_exp(x) result(str)
    real(C_DOUBLE), intent(in) :: x
    character(len=16) :: str
    write(str, '(ES10.2)') x
    str = adjustl(str)
  end function

  function itoa(i) result(str)
    integer, intent(in) :: i
    character(len=16) :: str
    write(str, '(I0)') i
  end function

end program test_fftw
