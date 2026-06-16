!-----------------------------------------------------------------------
! fourier.f90 — Fourier decomposition with energy ranking
!
! Purpose:
!   Provides nekStab-compatible wrappers around the nekstab_fourier_fftw
!   low-level FFT interface with energy-based mode sorting.
!   Works with both FFTW3 (GCC) and Intel MKL.
!
! Public interface:
!   nek_fourier_decomposition — FFT with energy-ranked modes
!   (re-exports: fourier_decomposition, fourier_reconstruction,
!    fft_init, fft_cleanup, fft_r2c, fft_c2r, fft_frequencies)
!
! Dependencies:
!   nekstab_krylov_subspace, nekstab_fourier_fftw, iso_c_binding
!-----------------------------------------------------------------------
module nekstab_fourier
   use nekstab_krylov_subspace, only: lv
   use nekstab_fourier_fftw, only: fourier_decomposition, fourier_reconstruction, &
                                   fft_init, fft_cleanup, fft_r2c, fft_c2r, &
                                   fft_frequencies
   use, intrinsic :: iso_c_binding, only: c_double, c_double_complex
   implicit none
   private

   ! Public procedures - re-export low-level routines
   public :: fourier_decomposition, fourier_reconstruction
   public :: fft_init, fft_cleanup, fft_r2c, fft_c2r, fft_frequencies

   ! High-level routines with energy sorting
   public :: nek_fourier_decomposition

contains

   !-----------------------------------------------------------------------
   ! nek_fourier_decomposition — Fourier decomposition with energy ranking
   !
   ! Purpose:
   !   Performs FFT, computes energy per mode, sorts by energy descending.
   !
   ! Arguments:
   !   nv       [in]    — number of velocity spatial points
   !   nsnap    [in]    — number of time snapshots
   !   vx,vy,vz [inout] — velocity components (nv x nsnap)
   !   time     [in]    — time array (nsnap)
   !   bm1      [in]    — mass matrix for energy weighting (nv)
   !   nmodes   [out]   — number of modes (nsnap/2+1)
   !   freq_out [out]   — frequencies sorted by energy (nmodes)
   !   energy   [out]   — energy of each mode, sorted (nmodes)
   !   vx_hat   [out]   — FFT coefficients sorted by energy (nv x nmodes)
   !   vy_hat   [out]   — FFT coefficients sorted by energy (nv x nmodes)
   !   vz_hat   [out]   — FFT coefficients sorted by energy (nv x nmodes)
   !-----------------------------------------------------------------------
   subroutine nek_fourier_decomposition(nv, nsnap, vx, vy, vz, time, bm1, &
                                        nmodes, freq_out, energy, &
                                        vx_hat, vy_hat, vz_hat)
      integer, intent(in) :: nv, nsnap
      real(c_double), intent(inout) :: vx(nv, nsnap)
      real(c_double), intent(inout) :: vy(nv, nsnap)
      real(c_double), intent(inout) :: vz(nv, nsnap)
      real(c_double), intent(in) :: time(nsnap)
      real(c_double), intent(in) :: bm1(nv)
      integer, intent(out) :: nmodes
      real(c_double), intent(out) :: freq_out(nsnap/2 + 1)
      real(c_double), intent(out) :: energy(nsnap/2 + 1)
      complex(c_double_complex), intent(out) :: vx_hat(nv, nsnap/2 + 1)
      complex(c_double_complex), intent(out) :: vy_hat(nv, nsnap/2 + 1)
      complex(c_double_complex), intent(out) :: vz_hat(nv, nsnap/2 + 1)

      integer :: i, k, nfreq
      real(c_double) :: freq(nsnap/2 + 1)
      real(c_double) :: mode_energy(nsnap/2 + 1)
      integer :: order(nsnap/2 + 1)
      ! Use allocatable for large arrays to avoid stack overflow
      complex(c_double_complex), allocatable :: vx_tmp(:, :), vy_tmp(:, :), vz_tmp(:, :)

      allocate (vx_tmp(nv, nsnap/2 + 1))
      allocate (vy_tmp(nv, nsnap/2 + 1))
      allocate (vz_tmp(nv, nsnap/2 + 1))

      ! Perform FFT decomposition
      call fourier_decomposition(nv, nsnap, vx, vy, vz, time, bm1, &
                                 nfreq, freq, vx_tmp, vy_tmp, vz_tmp)

      ! Compute energy per mode: E_k = sum_i bm1(i) * |v_hat(i,k)|^2
      do k = 1, nfreq
         mode_energy(k) = 0.0d0
         do i = 1, nv
            mode_energy(k) = mode_energy(k) + bm1(i)*( &
                             abs(vx_tmp(i, k))**2 + abs(vy_tmp(i, k))**2 + abs(vz_tmp(i, k))**2)
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

      deallocate (vx_tmp, vy_tmp, vz_tmp)

   end subroutine nek_fourier_decomposition

   !-----------------------------------------------------------------------
   ! sort_by_energy — sort indices by energy (descending, in-place)
   !-----------------------------------------------------------------------
   subroutine sort_by_energy(n, energy, order)
      integer, intent(in) :: n
      real(c_double), intent(inout) :: energy(n)
      integer, intent(out) :: order(n)

      integer :: i, j, max_idx
      real(c_double) :: max_val, temp_e
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

end module nekstab_fourier
