!-----------------------------------------------------------------------
! modal_pod.f90 — POD via method of snapshots (Sirovich, 1987)
!
! Purpose:
!   Computes Proper Orthogonal Decomposition using the snapshot
!   method. Includes POD-FFT spectral analysis of temporal
!   coefficients via Welch's method.
!
! Public interface:
!   pod_compute      — Full POD: correlation, eigensolve, modes
!   pod_fft_spectrum — Welch power spectrum of POD temporal coeffs
!   hamming_window   — Hamming window with amplitude/energy norm
!
! Dependencies:
!   krylov_subspace, fourier
!-----------------------------------------------------------------------
module modal_pod

   use krylov_subspace
   use krylov_inner_products
   use nekstab_lapack
   use nekstab_vectors
   use nekstab_nek_bridge
   use fourier

   implicit none
   private

   public :: pod_compute
   public :: pod_fft_spectrum
   public :: hamming_window

contains

!-----------------------------------------------------------------------
! pod_compute — POD via method of snapshots
!
!   C(i,j) = <snaps(i), snaps(j)>_E / n
!   C v = lambda v
!   Phi_k = Sum_i v_k(i) snaps(i) / sqrt(lambda_k * n)
!-----------------------------------------------------------------------
subroutine pod_compute(snaps, nsnap, nsave)

   integer, intent(in) :: nsnap, nsave
   type(krylov_vector), intent(in) :: snaps(nsnap)

   real, allocatable :: C(:,:), eigvals(:), eigvecs(:,:)
   real, allocatable :: norms(:)
   integer :: i, j


!  Form symmetric correlation matrix
   if (nid == 0) write(6,*) '  Forming correlation matrix...'

   allocate(C(nsnap, nsnap))
   allocate(eigvals(nsnap))
   allocate(eigvecs(nsnap, nsnap))

   call k_gram_matrix(C, snaps, snaps, nsnap, nsnap, nsnap)
   C = C / dble(nsnap)

!  Eigensolve (descending order)
   if (nid == 0) write(6,*) '  Solving eigenvalue problem...'
   call eig_symmetric(C, eigvals, eigvecs, nsnap)

!  Reconstruct and save modes
   if (nid == 0) write(6,*) '  Reconstructing POD modes...'

   allocate(norms(min(nsave, nsnap)))
   call pod_reconstruct_modes(snaps, eigvecs, eigvals, &
                              nsnap, nsave, norms)

!  Write eigenvalue spectrum
   call pod_write_spectrum(eigvals, nsnap, norms, &
                           min(nsave, nsnap))

   deallocate(C, eigvals, eigvecs, norms)

   if (nid == 0) write(6,*) '  POD complete'

end subroutine pod_compute

!-----------------------------------------------------------------------
! pod_reconstruct_modes — Reconstruct and output POD spatial modes
!
!   Phi_m = Sum_i v_m(i) * snaps(i) / sqrt(lambda_m * nsnap)
!-----------------------------------------------------------------------
subroutine pod_reconstruct_modes(snaps, eigvecs, eigvals, &
                                  nsnap, nsave, norms)

   integer, intent(in) :: nsnap, nsave
   type(krylov_vector), intent(in) :: snaps(nsnap)
   real, intent(in) :: eigvecs(nsnap, nsnap), eigvals(nsnap)
   real, intent(out) :: norms(min(nsave, nsnap))

   type(krylov_vector) :: mode
   real :: scale, mode_norm
   integer :: m, i, nmodes
   character(len=3) :: prefix

   prefix = 'pod'
   nmodes = min(nsave, nsnap)

   do m = 1, nmodes
      call k_zero(mode)

      do i = 1, nsnap
         call k_add2s2(mode, snaps(i), eigvecs(i, m))
      end do

      if (eigvals(m) > 0.0d0) then
         scale = 1.0d0 / sqrt(eigvals(m) * dble(nsnap))
         call k_cmult(mode, scale)
      end if

      call nopcopy(vx, vy, vz, pr, t, &
           mode%vx, mode%vy, mode%vz, mode%pr, mode%t)
      call outpost2(vx, vy, vz, pr, t, 1, prefix)

      call k_norm(mode_norm, mode)
      norms(m) = mode_norm

      if (nid == 0) then
         write(6,'(A,I4,A,E12.4,A,F8.4)') &
              '    Mode', m, ': lambda =', eigvals(m), &
              ', ||Phi|| =', mode_norm
      end if
   end do

end subroutine pod_reconstruct_modes

!-----------------------------------------------------------------------
! pod_write_spectrum — Write POD eigenvalue spectrum to file
!-----------------------------------------------------------------------
subroutine pod_write_spectrum(eigvals, n, norms, nnorms)

   integer, intent(in) :: n, nnorms
   real, intent(in) :: eigvals(n), norms(nnorms)

   real :: total_energy, cumsum
   integer :: m

   if (nid /= 0) return

   total_energy = sum(eigvals)
   if (total_energy <= 0.0d0) total_energy = 1.0d0

   open(unit=77, file='pod_spectrum.dat', status='replace')
   write(77, '(A)') '# POD Eigenvalue Spectrum'
   write(77, '(A)') &
        '# mode   eigenvalue      energy_%    cumulative_%    ||Phi||'

   cumsum = 0.0d0
   do m = 1, n
      cumsum = cumsum + eigvals(m)
      if (m <= nnorms) then
         write(77, '(I6, 4E16.8)') m, eigvals(m), &
              100.0d0 * eigvals(m) / total_energy, &
              100.0d0 * cumsum / total_energy, &
              norms(m)
      else
         write(77, '(I6, 3E16.8)') m, eigvals(m), &
              100.0d0 * eigvals(m) / total_energy, &
              100.0d0 * cumsum / total_energy
      end if
   end do

   close(77)
   write(6,*) '  Wrote pod_spectrum.dat'

end subroutine pod_write_spectrum

!-----------------------------------------------------------------------
! pod_fft_spectrum — Welch power spectrum of POD temporal coefficients
!
!   Fast alternative to field-based SPOD: FFT the POD temporal
!   coefficients a_k(t) = sqrt(lambda_k * n) * v_k(t) instead of
!   every spatial point.
!-----------------------------------------------------------------------
subroutine pod_fft_spectrum(snaps, nsnap, delta_t, nfft, noverlap)

   integer, intent(in) :: nsnap, nfft, noverlap
   real, intent(in) :: delta_t
   type(krylov_vector), intent(in) :: snaps(nsnap)

   real, allocatable :: C(:,:), eigvals(:), eigvecs(:,:)
   real, allocatable :: coeffs(:,:)
   real, allocatable :: window(:), power(:,:), freq(:)
   real, allocatable :: time_series(:)
   complex(nekStab_dp), allocatable :: spectrum(:)
   real :: win_weight, df, norm_factor
   integer :: i, j, k, m, nblk, nfreq, blk_start, r
   integer :: npeak
   integer, allocatable :: peak_idx(:)


!  Compute POD correlation matrix and eigenvectors
   if (nid == 0) write(6,*) '  Computing POD correlation matrix...'

   allocate(C(nsnap, nsnap))
   allocate(eigvals(nsnap), eigvecs(nsnap, nsnap))

   call k_gram_matrix(C, snaps, snaps, nsnap, nsnap, nsnap)
   C = C / dble(nsnap)

   call eig_symmetric(C, eigvals, eigvecs, nsnap)

!  Determine number of modes to analyze (up to 99.99% energy, capped at 30)
   r = nsnap
   do i = 1, nsnap
      if (sum(eigvals(1:i)) / sum(eigvals) > 0.9999d0) then
         r = i
         exit
      end if
   end do
   r = max(r, min(10, nsnap))
   r = min(r, 30)

   if (nid == 0) write(6,'(A,I4,A)') &
        '  Analyzing', r, ' POD modes'

!  Extract temporal coefficients: a_k(t) = sqrt(lambda_k * n) * v_k(t)
   allocate(coeffs(nsnap, r))
   do k = 1, r
      if (eigvals(k) > 0) then
         do i = 1, nsnap
            coeffs(i, k) = eigvecs(i, k) * sqrt(eigvals(k) * &
                 nsnap)
         end do
      else
         coeffs(:, k) = 0.0d0
      end if
   end do

!  Write temporal coefficients to file
   if (nid == 0) then
      open(unit=77, file='pod_coefficients.dat', &
           status='replace')
      write(77, '(A)') '# POD Temporal Coefficients'
      write(77, '(A,I6,A,I4)') &
           '# nsnap = ', nsnap, ',  nmodes = ', r
      write(77, '(A)', advance='no') '# snapshot'
      do k = 1, r
         write(77, '(A,I3)', advance='no') '         a_', k
      end do
      write(77, *)

      do i = 1, nsnap
         write(77, '(I6)', advance='no') i
         do k = 1, r
            write(77, '(E14.6)', advance='no') coeffs(i, k)
         end do
         write(77, *)
      end do

      close(77)
      write(6,*) '  Wrote pod_coefficients.dat'
   end if

   deallocate(C, eigvals, eigvecs)

!  Setup FFT parameters
   nblk = (nsnap - nfft) / (nfft - noverlap) + 1
   nfreq = nfft / 2 + 1
   df = 1.0d0 / (nfft * delta_t)

   if (nblk < 1) then
      if (nid == 0) then
         write(6,*) '  ERROR: Not enough snapshots for FFT'
         write(6,*) '    nsnap =', nsnap, ', nfft =', nfft
      end if
      deallocate(coeffs)
      return
   end if

   if (nid == 0) then
      write(6,'(A,I6)') '    nfft:     ', nfft
      write(6,'(A,I6)') '    noverlap: ', noverlap
      write(6,'(A,I6)') '    nblk:     ', nblk
      write(6,'(A,I6)') '    nfreq:    ', nfreq
      write(6,'(A,E12.4)') '    dSt:      ', df
   end if

!  Compute Welch power spectrum for each mode
   if (nid == 0) write(6,*) '  Computing power spectra...'

   allocate(window(nfft), freq(nfreq))
   allocate(power(nfreq, r))
   allocate(time_series(nfft), spectrum(nfreq))

   call hamming_window(nfft, window, win_weight)

   do i = 1, nfreq
      freq(i) = dble(i - 1) * df
   end do

   power = 0.0d0

   do k = 1, r
      do m = 1, nblk
         blk_start = (m - 1) * (nfft - noverlap) + 1

         do i = 1, nfft
            time_series(i) = window(i) * coeffs(blk_start+i-1, k)
         end do

         call fft_r2c(nfft, time_series, spectrum)

         do i = 1, nfreq
            norm_factor = 2.0d0 / (dble(nfft) * win_weight)
            if (i == 1 .or. i == nfreq) &
                 norm_factor = norm_factor / 2.0d0
            power(i, k) = power(i, k) + &
                 (real(spectrum(i))**2 + aimag(spectrum(i))**2) &
                 * norm_factor**2
         end do
      end do

      power(:, k) = power(:, k) / dble(nblk)
   end do

!  Find peaks and write spectrum file
   npeak = min(r, 10)
   allocate(peak_idx(npeak))

   if (nid == 0) then
      do k = 1, npeak
         peak_idx(k) = 2  ! Skip DC
         do i = 3, nfreq
            if (power(i, k) > power(peak_idx(k), k)) &
                 peak_idx(k) = i
         end do
      end do

      open(unit=79, file='pod_fft_spectrum.dat', &
           status='replace')
      write(79, '(A)') '# POD-FFT Spectral Analysis'
      write(79, '(A,I6,A,I6,A,I6,A,I6,A,E12.4)') &
           '# nfft = ', nfft, ', noverlap = ', noverlap, &
           ', nblk = ', nblk, ', nfreq = ', nfreq, &
           ', dSt = ', df
      write(79, '(A)') &
           '# Peak frequencies (excluding DC):'
      do k = 1, npeak
         write(79, '(A,I3,A,F8.4,A,E12.4)') &
              '#   Mode', k, ': St = ', freq(peak_idx(k)), &
              ', power = ', power(peak_idx(k), k)
      end do

      write(79, '(A)', advance='no') '# St'
      do k = 1, r
         write(79, '(A,I3)', advance='no') '  P_mode', k
      end do
      write(79, *)

      do i = 1, nfreq
         write(79, '(E14.6)', advance='no') freq(i)
         do k = 1, r
            write(79, '(E14.6)', advance='no') power(i, k)
         end do
         write(79, *)
      end do

      close(79)
      write(6,*) '  Wrote pod_fft_spectrum.dat'

      write(6,*) ''
      write(6,*) '  Peak frequencies (excluding DC):'
      do k = 1, npeak
         write(6,'(A,I3,A,F8.4,A,E12.4)') &
              '    Mode', k, ': St =', freq(peak_idx(k)), &
              ', power =', power(peak_idx(k), k)
      end do
   end if

   deallocate(coeffs, window, freq, power, time_series, spectrum)
   deallocate(peak_idx)

   if (nid == 0) write(6,*) '  POD-FFT analysis complete'

end subroutine pod_fft_spectrum

!-----------------------------------------------------------------------
! hamming_window — Compute Hamming window and normalization weight
!
!   ifwinamp = .true.  -> Amplitude: win_weight = mean(w)  [PySPOD]
!   ifwinamp = .false. -> Energy: win_weight = sqrt(sum(w^2)/n)
!-----------------------------------------------------------------------
subroutine hamming_window(n, window, win_weight)

   integer, intent(in) :: n
   real, intent(out) :: window(n), win_weight
   real :: win_mean, win_energy
   integer :: i

   do i = 1, n
      window(i) = 0.54d0 - 0.46d0 * cos(2.0d0 * NEKSTAB_PI * &
           dble(i-1) / dble(n-1))
   end do

   win_mean = sum(window) / dble(n)
   win_energy = sum(window**2) / dble(n)

   if (ifwinamp) then
      win_weight = win_mean
   else
      win_weight = sqrt(win_energy)
   end if

end subroutine hamming_window

end module modal_pod
