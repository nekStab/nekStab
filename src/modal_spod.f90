      !-----------------------------------------------------------------------
      ! modal_spod.f90 — Spectral POD (Towne et al., 2018) batch version
      !
      ! Purpose:
      !   1. Divide snapshots into overlapping blocks
      !   2. Apply window and FFT each block
      !   3. For each frequency, form CSD matrix and eigensolve
      !
      ! Public interface:
      !   spod_compute    — Full batch SPOD computation
      !   k_dot_complex   — Hermitian inner product for complex fields
      !   spod_save_modes — Save SPOD modes at selected frequencies
      !
      ! Dependencies:
      !   krylov_subspace, fourier, modal_pod, SIZE, TOTAL
      !-----------------------------------------------------------------------
      module modal_spod

         use krylov_subspace
         use fourier
         use modal_pod, only: hamming_window

         implicit none
         private

         public :: spod_compute
         public :: k_dot_complex
         public :: spod_save_modes

      contains

      !-----------------------------------------------------------------------
      ! spod_compute — Spectral POD via batch processing
      !-----------------------------------------------------------------------
      subroutine spod_compute(snaps, nsnap, delta_t, nfft, noverlap,
     $                        nsave)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, nfft, noverlap, nsave
         real, intent(in) :: delta_t
         type(krylov_vector), intent(in) :: snaps(nsnap)

         integer :: nblk, nfreq, iblk, ifreq, blk_start, i, j
         real, allocatable :: window(:), freq(:)
         real :: win_weight, df
         type(krylov_vector), allocatable :: blk_re(:), blk_im(:)
         complex(kind=kind(0.0d0)), allocatable :: CSD(:,:)
         real, allocatable :: spod_evals(:)
         complex(kind=kind(0.0d0)), allocatable :: spod_evecs(:,:)
         complex(kind=kind(0.0d0)) :: cval

         external :: eig_hermitian

         if (nid == 0) then
            write(6,*) '  SPOD parameters:'
            write(6,'(A,I6)') '    nfft:     ', nfft
            write(6,'(A,I6)') '    noverlap: ', noverlap
         end if

         nblk = (nsnap - nfft) / (nfft - noverlap) + 1
         nfreq = nfft / 2 + 1
         df = 1.0d0 / (nfft * delta_t)

         if (nblk < 2) then
            if (nid == 0) then
               write(6,*) '  ERROR: Not enough snapshots for SPOD'
               write(6,*) '    Need at least 2 blocks'
               write(6,*) '    Current: nblk =', nblk
            end if
            return
         end if

         if (nid == 0) then
            write(6,'(A,I6)') '    nblk:     ', nblk
            write(6,'(A,I6)') '    nfreq:    ', nfreq
            write(6,'(A,E12.4)') '    dSt:      ', df
         end if

         allocate(window(nfft), freq(nfreq))
         allocate(blk_re(nblk), blk_im(nblk))
         allocate(CSD(nblk, nblk))
         allocate(spod_evals(nblk), spod_evecs(nblk, nblk))

         call hamming_window(nfft, window, win_weight)
         call fft_frequencies(nfft, delta_t, freq)

!        Open spectrum file
         if (nid == 0) then
            open(unit=79, file='spod_spectrum.dat',
     $           status='replace')
            write(79, '(A)')
     $         '# SPOD Spectrum: St, eigenvalues (nblk columns)'
         end if

!        Loop over frequencies
         if (nid == 0) write(6,'(A,I4,A)')
     $      '   Processing ', nfreq, ' frequencies...'

         do ifreq = 1, nfreq

            if (nid == 0 .and. mod(ifreq,5) == 1) then
               write(6,'(A,I4,A,I4)') '     freq ', ifreq, ' / ',nfreq
            end if

!           FFT each block at this frequency
            do iblk = 1, nblk
               blk_start = (iblk - 1) * (nfft - noverlap) + 1
               call fft_block_at_freq(snaps, blk_start, nfft,
     $              window, win_weight, delta_t, ifreq,
     $              blk_re(iblk), blk_im(iblk))
            end do

!           Form Hermitian CSD matrix
            do j = 1, nblk
               do i = 1, j
                  call k_dot_complex(cval,
     $                 blk_re(i), blk_im(i),
     $                 blk_re(j), blk_im(j))
                  CSD(i,j) = cval
                  CSD(j,i) = conjg(cval)
               end do
            end do
            CSD = CSD / dble(nblk)

            call eig_hermitian(CSD, spod_evals, spod_evecs, nblk)

            if (nid == 0) then
               write(79, '(E14.6)', advance='no') freq(ifreq)
               do i = 1, nblk
                  write(79, '(E14.6)', advance='no') spod_evals(i)
               end do
               write(79, *)
            end if

            call spod_save_modes(blk_re, blk_im, spod_evecs,
     $           spod_evals, nblk, ifreq, nfreq, freq, nsave)

         end do

         if (nid == 0) then
            close(79)
            write(6,*) '  Wrote spod_spectrum.dat'
         end if

         deallocate(window, freq, blk_re, blk_im)
         deallocate(CSD, spod_evals, spod_evecs)

         if (nid == 0) write(6,*) '  SPOD complete'

      end subroutine spod_compute

      !-----------------------------------------------------------------------
      ! fft_block_at_freq — FFT a block of snapshots, extract one frequency
      !-----------------------------------------------------------------------
      subroutine fft_block_at_freq(snaps, blk_start, nfft,
     $     window, win_weight, delta_t, ifreq, out_re, out_im)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: blk_start, nfft, ifreq
         real, intent(in) :: delta_t, win_weight
         real, intent(in) :: window(nfft)
         type(krylov_vector), intent(in) :: snaps(*)
         type(krylov_vector), intent(out) :: out_re, out_im

         real :: time_series(nfft)
         complex(kind=kind(0.0d0)) :: spectrum(nfft/2+1)
         real :: norm_factor
         integer :: ipt, j, nfreq

         nfreq = nfft / 2 + 1

         norm_factor = 2.0d0 / (dble(nfft) * win_weight)
         if (ifreq == 1 .or. ifreq == nfreq) then
            norm_factor = norm_factor / 2.0d0
         end if

         call k_zero(out_re)
         call k_zero(out_im)

         nv = nx1*ny1*nz1*nelv
         nt = nx1*ny1*nz1*nelt

!        Velocity x
         do ipt = 1, nv
            do j = 1, nfft
               time_series(j) = window(j) *
     $              snaps(blk_start + j - 1)%vx(ipt)
            end do
            call fft_r2c(nfft, time_series, spectrum)
            out_re%vx(ipt) = real(spectrum(ifreq)) * norm_factor
            out_im%vx(ipt) = aimag(spectrum(ifreq)) * norm_factor
         end do

!        Velocity y
         do ipt = 1, nv
            do j = 1, nfft
               time_series(j) = window(j) *
     $              snaps(blk_start + j - 1)%vy(ipt)
            end do
            call fft_r2c(nfft, time_series, spectrum)
            out_re%vy(ipt) = real(spectrum(ifreq)) * norm_factor
            out_im%vy(ipt) = aimag(spectrum(ifreq)) * norm_factor
         end do

!        Velocity z (3D only)
         if (if3d) then
            do ipt = 1, nv
               do j = 1, nfft
                  time_series(j) = window(j) *
     $                 snaps(blk_start + j - 1)%vz(ipt)
               end do
               call fft_r2c(nfft, time_series, spectrum)
               out_re%vz(ipt) = real(spectrum(ifreq)) * norm_factor
               out_im%vz(ipt) = aimag(spectrum(ifreq)) * norm_factor
            end do
         end if

!        Temperature (if thermal)
         if (ifto) then
            do ipt = 1, nt
               do j = 1, nfft
                  time_series(j) = window(j) *
     $                 snaps(blk_start + j - 1)%t(ipt, 1)
               end do
               call fft_r2c(nfft, time_series, spectrum)
               out_re%t(ipt, 1) = real(spectrum(ifreq)) * norm_factor
               out_im%t(ipt, 1) = aimag(spectrum(ifreq)) * norm_factor
            end do
         end if

      end subroutine fft_block_at_freq

      !-----------------------------------------------------------------------
      ! k_dot_complex — Hermitian inner product: <p, q> = p^H * W * q
      !
      !   <p,q> = <p_re,q_re> + <p_im,q_im>
      !         + i(<p_re,q_im> - <p_im,q_re>)
      !-----------------------------------------------------------------------
      subroutine k_dot_complex(alpha, p_re, p_im, q_re, q_im)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         type(krylov_vector), intent(in) :: p_re, p_im, q_re, q_im
         complex(kind=kind(0.0d0)), intent(out) :: alpha

         real :: rr, ii, ri, ir

         call k_dot(rr, p_re, q_re)
         call k_dot(ii, p_im, q_im)
         call k_dot(ri, p_re, q_im)
         call k_dot(ir, p_im, q_re)

         alpha = dcmplx(rr + ii, ri - ir)

      end subroutine k_dot_complex

      !-----------------------------------------------------------------------
      ! spod_save_modes — Save SPOD modes at selected frequencies
      !
      !   Saves at DC, Nyquist, and every 8th frequency.
      !-----------------------------------------------------------------------
      subroutine spod_save_modes(blk_re, blk_im, evecs, evals,
     $     nblk, ifreq, nfreq, freq, nsave)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nblk, ifreq, nfreq, nsave
         type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)
         complex(kind=kind(0.0d0)), intent(in) :: evecs(nblk, nblk)
         real, intent(in) :: evals(nblk), freq(nfreq)

         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm
         integer :: m, i, nmodes
         logical :: save_this_freq

!        Decide which frequencies to save
         save_this_freq = .false.
         if (ifreq == 1) save_this_freq = .true.
         if (ifreq == nfreq) save_this_freq = .true.
         if (mod(ifreq-1, max(1, nfreq/8)) == 0) save_this_freq = .true.

         if (.not. save_this_freq) return

         nmodes = min(nsave, nblk)

         do m = 1, nmodes
            call k_zero(mode_re)
            call k_zero(mode_im)

!           Phi_m = Sum_i w_m(i) * Qhat(i) / sqrt(lambda_m * nblk)
!           Complex multiplication: (a+ib)(c+id) = (ac-bd) + i(ad+bc)
            do i = 1, nblk
               coef_re = real(evecs(i, m))
               coef_im = aimag(evecs(i, m))

               call k_add2s2(mode_re, blk_re(i), coef_re)
               call k_add2s2(mode_re, blk_im(i), -coef_im)
               call k_add2s2(mode_im, blk_re(i), coef_im)
               call k_add2s2(mode_im, blk_im(i), coef_re)
            end do

!           Normalize
            if (evals(m) > 0.0d0) then
               coef_re = 1.0d0 / sqrt(evals(m) * dble(nblk))
               call k_cmult(mode_re, coef_re)
               call k_cmult(mode_im, coef_re)
            end if

!           Output real part
            call nopcopy(vx, vy, vz, pr, t,
     $           mode_re%vx, mode_re%vy, mode_re%vz,
     $           mode_re%pr, mode_re%t)
            call outpost2(vx, vy, vz, pr, t, 0, 'sRe')

!           Output imaginary part
            call nopcopy(vx, vy, vz, pr, t,
     $           mode_im%vx, mode_im%vy, mode_im%vz,
     $           mode_im%pr, mode_im%t)
            call outpost2(vx, vy, vz, pr, t, 0, 'sIm')

!           Print info for leading mode only
            if (m == 1) then
               call k_norm(mode_norm, mode_re)
               if (nid == 0) then
                  write(6,'(A,F10.4,A,E12.4,A,E12.4)')
     $                 '    St =', freq(ifreq), ': lambda_1 =',
     $                 evals(1), ', ||Phi|| =', mode_norm
               end if
            end if
         end do

      end subroutine spod_save_modes

      end module modal_spod
