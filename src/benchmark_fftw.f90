!-----------------------------------------------------------------------
! benchmark_fftw.f90: Performance benchmark for FFTW3 FFT operations
!
! Tests:
!   - Single FFT performance for various sizes
!   - Batch FFT (simulating spatial points) performance
!   - Comparison of power-of-2 vs non-power-of-2 sizes
!
! Build:
!   GCC:   gfortran -O3 -o benchmark_fftw fourier_fftw.f90 benchmark_fftw.f90 -lfftw3
!   Intel: ifort/ifx -O3 -o benchmark_fftw fourier_fftw.f90 benchmark_fftw.f90 -qmkl
!-----------------------------------------------------------------------
program benchmark_fftw
  use fourier_fftw
  use, intrinsic :: iso_c_binding
  implicit none

  integer :: i
  integer, parameter :: n_warmup = 3
  integer, parameter :: n_repeat = 10

  print '(A)', '========================================================'
  print '(A)', '  FFTW3 Performance Benchmark for nekStab'
  print '(A)', '========================================================'
  print '(A)', ''

  ! Test 1: Single FFT at various sizes
  call benchmark_single_fft()

  ! Test 2: Batch FFT (simulating npts spatial points)
  call benchmark_batch_fft()

  ! Test 3: Full decomposition/reconstruction cycle
  call benchmark_decomposition()

  call fft_cleanup()

contains

  !-------------------------------------------------------------------
  subroutine benchmark_single_fft()
    integer :: sizes(8), n, j, k
    real(C_DOUBLE), allocatable :: signal(:), reconstructed(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: spectrum(:)
    real(C_DOUBLE) :: t_start, t_end, t_fwd, t_inv, t_total
    real(C_DOUBLE) :: mflops

    print '(A)', 'Test 1: Single FFT Performance'
    print '(A)', '------------------------------'
    print '(A)', ''
    print '(A)', '    N       Forward(us)  Inverse(us)  Total(us)    MFLOPS'
    print '(A)', '  -----     ----------   ----------   ---------    ------'

    ! Mix of power-of-2 and other sizes
    sizes = [64, 128, 256, 512, 1024, 100, 144, 360]

    do j = 1, size(sizes)
      n = sizes(j)
      allocate(signal(n), reconstructed(n), spectrum(n/2+1))

      ! Initialize with test signal
      do i = 1, n
        signal(i) = sin(0.1d0*i) + 0.5d0*cos(0.3d0*i)
      end do

      ! Warmup
      do k = 1, n_warmup
        call fft_r2c(n, signal, spectrum)
        call fft_c2r(n, spectrum, reconstructed)
      end do

      ! Benchmark forward FFT
      call cpu_time(t_start)
      do k = 1, n_repeat
        call fft_r2c(n, signal, spectrum)
      end do
      call cpu_time(t_end)
      t_fwd = (t_end - t_start) / n_repeat * 1.0d6  ! microseconds

      ! Benchmark inverse FFT
      call cpu_time(t_start)
      do k = 1, n_repeat
        call fft_c2r(n, spectrum, reconstructed)
      end do
      call cpu_time(t_end)
      t_inv = (t_end - t_start) / n_repeat * 1.0d6  ! microseconds

      t_total = t_fwd + t_inv

      ! MFLOPS estimate: FFT is O(N log N), ~5N log2(N) operations
      mflops = 5.0d0 * n * log(real(n,8)) / log(2.0d0) / (t_total * 1.0d-6) / 1.0d6

      print '(I7, 4X, F10.2, 4X, F10.2, 4X, F9.2, 4X, F8.1)', &
            n, t_fwd, t_inv, t_total, mflops

      call fft_cleanup()
      deallocate(signal, reconstructed, spectrum)
    end do
    print '(A)', ''
  end subroutine benchmark_single_fft

  !-------------------------------------------------------------------
  subroutine benchmark_batch_fft()
    integer :: npts_list(4), nsnap_list(3)
    integer :: npts, nsnap, j, k, m
    real(C_DOUBLE), allocatable :: vx(:,:), vy(:,:), vz(:,:)
    real(C_DOUBLE), allocatable :: time_arr(:), bm1(:), freq(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: vx_hat(:,:), vy_hat(:,:), vz_hat(:,:)
    real(C_DOUBLE) :: t_start, t_end, t_decomp
    real(C_DOUBLE) :: dt
    integer :: nfreq

    print '(A)', 'Test 2: Batch FFT Performance (simulating nekStab workflow)'
    print '(A)', '-----------------------------------------------------------'
    print '(A)', ''
    print '(A)', '  npts      nsnap     Time(ms)    Rate(pts*snaps/s)'
    print '(A)', '  -----     -----     --------    ------------------'

    npts_list = [1000, 10000, 50000, 100000]
    nsnap_list = [64, 128, 256]

    do j = 1, size(npts_list)
      do m = 1, size(nsnap_list)
        npts = npts_list(j)
        nsnap = nsnap_list(m)

        allocate(vx(npts, nsnap), vy(npts, nsnap), vz(npts, nsnap))
        allocate(time_arr(nsnap), bm1(npts), freq(nsnap/2+1))
        allocate(vx_hat(npts, nsnap/2+1))
        allocate(vy_hat(npts, nsnap/2+1))
        allocate(vz_hat(npts, nsnap/2+1))

        ! Initialize
        dt = 0.01d0
        do i = 1, nsnap
          time_arr(i) = (i-1) * dt
        end do
        bm1 = 1.0d0

        ! Random-ish velocity field
        do i = 1, npts
          do k = 1, nsnap
            vx(i,k) = sin(0.1d0*i + 0.2d0*k)
            vy(i,k) = cos(0.15d0*i + 0.25d0*k)
            vz(i,k) = sin(0.05d0*i + 0.1d0*k)
          end do
        end do

        ! Warmup
        call fourier_decomposition(npts, nsnap, vx, vy, vz, time_arr, bm1, &
                                    nfreq, freq, vx_hat, vy_hat, vz_hat)

        ! Re-initialize (decomposition modifies input)
        do i = 1, npts
          do k = 1, nsnap
            vx(i,k) = sin(0.1d0*i + 0.2d0*k)
            vy(i,k) = cos(0.15d0*i + 0.25d0*k)
            vz(i,k) = sin(0.05d0*i + 0.1d0*k)
          end do
        end do

        ! Benchmark
        call cpu_time(t_start)
        call fourier_decomposition(npts, nsnap, vx, vy, vz, time_arr, bm1, &
                                    nfreq, freq, vx_hat, vy_hat, vz_hat)
        call cpu_time(t_end)
        t_decomp = (t_end - t_start) * 1.0d3  ! milliseconds

        print '(I7, 4X, I5, 5X, F8.2, 6X, ES12.2)', &
              npts, nsnap, t_decomp, real(npts,8)*nsnap/((t_end-t_start))

        call fft_cleanup()
        deallocate(vx, vy, vz, time_arr, bm1, freq, vx_hat, vy_hat, vz_hat)
      end do
    end do
    print '(A)', ''
  end subroutine benchmark_batch_fft

  !-------------------------------------------------------------------
  subroutine benchmark_decomposition()
    integer, parameter :: npts = 10000
    integer, parameter :: nsnap = 128
    real(C_DOUBLE) :: vx(npts, nsnap), vy(npts, nsnap), vz(npts, nsnap)
    real(C_DOUBLE) :: time_arr(nsnap), bm1(npts), freq(nsnap/2+1)
    complex(C_DOUBLE_COMPLEX) :: vx_hat(npts, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX) :: vy_hat(npts, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX) :: vz_hat(npts, nsnap/2+1)
    real(C_DOUBLE) :: vx_rec(npts), vy_rec(npts), vz_rec(npts)
    real(C_DOUBLE) :: t_start, t_end, t_decomp, t_recon
    real(C_DOUBLE) :: dt, t_test
    integer :: nfreq, i, k

    print '(A)', 'Test 3: Full Decomposition/Reconstruction Cycle'
    print '(A)', '------------------------------------------------'
    print '(A,I6,A,I4,A)', '  Configuration: npts=', npts, ', nsnap=', nsnap, ''
    print '(A)', ''

    ! Initialize
    dt = 0.01d0
    do i = 1, nsnap
      time_arr(i) = (i-1) * dt
    end do
    bm1 = 1.0d0

    do i = 1, npts
      do k = 1, nsnap
        vx(i,k) = sin(0.1d0*i + 0.2d0*k)
        vy(i,k) = cos(0.15d0*i + 0.25d0*k)
        vz(i,k) = sin(0.05d0*i + 0.1d0*k)
      end do
    end do

    ! Benchmark decomposition
    call cpu_time(t_start)
    call fourier_decomposition(npts, nsnap, vx, vy, vz, time_arr, bm1, &
                                nfreq, freq, vx_hat, vy_hat, vz_hat)
    call cpu_time(t_end)
    t_decomp = (t_end - t_start) * 1.0d3

    ! Benchmark reconstruction at multiple time points
    call cpu_time(t_start)
    do k = 1, 100
      t_test = time_arr(1) + (k-1) * dt * 0.5d0
      call fourier_reconstruction(npts, nfreq, freq, vx_hat, vy_hat, vz_hat, &
                                   t_test, vx_rec, vy_rec, vz_rec)
    end do
    call cpu_time(t_end)
    t_recon = (t_end - t_start) * 1.0d3 / 100.0d0  ! per reconstruction

    print '(A,F10.2,A)', '  Decomposition time:      ', t_decomp, ' ms'
    print '(A,F10.2,A)', '  Reconstruction time:     ', t_recon, ' ms (per time point)'
    print '(A,F10.2,A)', '  Memory for coefficients: ', &
          real(npts,8) * nfreq * 3 * 16 / 1024.0d0 / 1024.0d0, ' MB'
    print '(A)', ''

    call fft_cleanup()
  end subroutine benchmark_decomposition

end program benchmark_fftw
