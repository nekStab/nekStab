!-----------------------------------------------------------------------
! fourier.f90: Fourier decomposition/reconstruction for nekStab
!
! Purpose:
!   - Temporal FFT decomposition of velocity snapshots
!   - Reconstruction of velocity fields from Fourier modes
!   - Foundation for SPOD (Spectral Proper Orthogonal Decomposition)
!
! This module provides nekStab-compatible wrappers around the
! fourier_fftw low-level FFT interface. Works with both FFTW3 (GCC)
! and Intel MKL via the FFTW3 wrapper interface.
!
! Author: nekStab team
! Date: 2026
!-----------------------------------------------------------------------
module fourier
  use krylov_subspace, only: lv
  use fourier_fftw
  use, intrinsic :: iso_c_binding, only: C_DOUBLE, C_DOUBLE_COMPLEX
  implicit none
  private

  ! Public procedures - re-export low-level routines
  public :: fourier_decomposition, fourier_reconstruction
  public :: fft_init, fft_cleanup, fft_r2c, fft_c2r, fft_frequencies

  ! High-level routines with energy sorting
  public :: nek_fourier_decomposition

contains

  !---------------------------------------------------------------------
  ! nek_fourier_decomposition: Fourier decomposition with energy ranking
  !
  ! Performs FFT, computes energy per mode, sorts by energy descending.
  !
  ! Arguments:
  !   npts     - Number of spatial points
  !   nsnap    - Number of time snapshots
  !   vx,vy,vz - Velocity components (npts x nsnap), MODIFIED
  !   time     - Time array (nsnap)
  !   bm1      - Mass matrix for energy weighting (npts)
  !   nmodes   - Output: number of modes (nsnap/2+1)
  !   freq_out - Output: frequencies sorted by energy (nmodes)
  !   energy   - Output: energy of each mode, sorted (nmodes)
  !   vx_hat   - Output: FFT coefficients sorted by energy (npts x nmodes)
  !   vy_hat   - Output: FFT coefficients sorted by energy (npts x nmodes)
  !   vz_hat   - Output: FFT coefficients sorted by energy (npts x nmodes)
  !---------------------------------------------------------------------
  subroutine nek_fourier_decomposition(npts, nsnap, vx, vy, vz, time, bm1, &
                                        nmodes, freq_out, energy, &
                                        vx_hat, vy_hat, vz_hat)
    integer, intent(in) :: npts, nsnap
    real(C_DOUBLE), intent(inout) :: vx(npts, nsnap)
    real(C_DOUBLE), intent(inout) :: vy(npts, nsnap)
    real(C_DOUBLE), intent(inout) :: vz(npts, nsnap)
    real(C_DOUBLE), intent(in) :: time(nsnap)
    real(C_DOUBLE), intent(in) :: bm1(npts)
    integer, intent(out) :: nmodes
    real(C_DOUBLE), intent(out) :: freq_out(nsnap/2+1)
    real(C_DOUBLE), intent(out) :: energy(nsnap/2+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: vx_hat(npts, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: vy_hat(npts, nsnap/2+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: vz_hat(npts, nsnap/2+1)

    integer :: i, k, nfreq
    real(C_DOUBLE) :: freq(nsnap/2+1)
    real(C_DOUBLE) :: mode_energy(nsnap/2+1)
    integer :: order(nsnap/2+1)
    ! Use allocatable for large arrays to avoid stack overflow
    complex(C_DOUBLE_COMPLEX), allocatable :: vx_tmp(:,:), vy_tmp(:,:), vz_tmp(:,:)

    allocate(vx_tmp(npts, nsnap/2+1))
    allocate(vy_tmp(npts, nsnap/2+1))
    allocate(vz_tmp(npts, nsnap/2+1))

    ! Perform FFT decomposition
    call fourier_decomposition(npts, nsnap, vx, vy, vz, time, bm1, &
                                nfreq, freq, vx_tmp, vy_tmp, vz_tmp)

    ! Compute energy per mode: E_k = sum_i bm1(i) * |v_hat(i,k)|^2
    do k = 1, nfreq
      mode_energy(k) = 0.0d0
      do i = 1, npts
        mode_energy(k) = mode_energy(k) + bm1(i) * ( &
          abs(vx_tmp(i,k))**2 + abs(vy_tmp(i,k))**2 + abs(vz_tmp(i,k))**2)
      end do
    end do

    ! Sort modes by energy (descending)
    call sort_by_energy(nfreq, mode_energy, order)

    ! Copy sorted results to output
    nmodes = nfreq
    do k = 1, nfreq
      freq_out(k) = freq(order(k))
      energy(k) = mode_energy(k)
      vx_hat(:, k) = vx_tmp(:, order(k))
      vy_hat(:, k) = vy_tmp(:, order(k))
      vz_hat(:, k) = vz_tmp(:, order(k))
    end do

    deallocate(vx_tmp, vy_tmp, vz_tmp)
  end subroutine nek_fourier_decomposition

  !---------------------------------------------------------------------
  ! sort_by_energy: Sort indices by energy in descending order
  ! Note: energy array is modified (sorted in place)
  !---------------------------------------------------------------------
  subroutine sort_by_energy(n, energy, order)
    integer, intent(in) :: n
    real(C_DOUBLE), intent(inout) :: energy(n)
    integer, intent(out) :: order(n)

    integer :: i, j, max_idx
    real(C_DOUBLE) :: max_val, temp_e
    integer :: temp_o

    ! Initialize order
    do i = 1, n
      order(i) = i
    end do

    ! Selection sort (descending) - O(n^2) but n is typically small (<256)
    do i = 1, n - 1
      max_idx = i
      max_val = energy(i)
      do j = i + 1, n
        if (energy(j) > max_val) then
          max_val = energy(j)
          max_idx = j
        end if
      end do
      if (max_idx /= i) then
        temp_e = energy(i)
        energy(i) = energy(max_idx)
        energy(max_idx) = temp_e
        temp_o = order(i)
        order(i) = order(max_idx)
        order(max_idx) = temp_o
      end if
    end do
  end subroutine sort_by_energy

end module fourier
