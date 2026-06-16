!-----------------------------------------------------------------------
! modal_spod_streaming.f90 — Streaming SPOD algorithm
!
! Purpose:
!   Implements a streaming variant of Spectral POD that processes
!   snapshots one at a time, accumulating DFT coefficients on the
!   fly. Faster than batch SPOD for large datasets due to sequential
!   access pattern and perfect cache locality.
!
! Public interface:
!   spod_streaming_batch  — process pre-loaded snapshots via streaming
!   spod_stream_init      — initialize streaming SPOD state
!   spod_stream_update    — accumulate one snapshot into DFT
!   spod_stream_finalize  — CSD eigensolve and mode output
!
! Dependencies:
!   nekstab_krylov_subspace, nekstab_spod_streaming_state, nekstab_modal_pod, nekstab_modal_spod
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! nekstab_spod_streaming_state — State management module for streaming SPOD
!-----------------------------------------------------------------------
module nekstab_spod_streaming_state
   use nekstab_krylov_subspace, only: krylov_vector, k_zero
   implicit none
   private

   public :: spod_s_init, spod_s_cleanup
   public :: spod_s_nfft, spod_s_noverlap, spod_s_nfreq, spod_s_nblk
   public :: spod_s_dt, spod_s_win_weight
   public :: spod_s_window, spod_s_t_idx, spod_s_n_complete
   public :: spod_s_n_snaps, spod_s_initialized
   public :: spod_s_mean, spod_s_center
   public :: spod_s_dft_re, spod_s_dft_im
   public :: spod_s_work, spod_s_freq

   ! Scalar parameters
   integer, save :: spod_s_nfft = 0
   integer, save :: spod_s_noverlap = 0
   integer, save :: spod_s_nfreq = 0
   integer, save :: spod_s_nblk = 0
   real, save :: spod_s_dt = 0.0d0
   real, save :: spod_s_win_weight = 0.0d0
   integer, save :: spod_s_n_snaps = 0
   logical, save :: spod_s_initialized = .false.
   logical, save :: spod_s_center = .true.

   ! Allocatable arrays
   real, allocatable, save :: spod_s_window(:)
   real, allocatable, save :: spod_s_freq(:)
   integer, allocatable, save :: spod_s_t_idx(:)
   integer, allocatable, save :: spod_s_n_complete(:)

   ! Mean snapshot (running average)
   type(krylov_vector), save :: spod_s_mean

   ! Work vector (reused to avoid stack allocation)
   type(krylov_vector), save :: spod_s_work

   ! DFT accumulators indexed as (ifreq-1)*nblk + iblk
   type(krylov_vector), allocatable, save :: spod_s_dft_re(:)
   type(krylov_vector), allocatable, save :: spod_s_dft_im(:)

contains

   ! spod_s_init — Allocate and initialize streaming SPOD state
   subroutine spod_s_init(nfft, noverlap, dt, nblk_in)
      integer, intent(in) :: nfft, noverlap, nblk_in
      real, intent(in) :: dt
      integer :: i, idx

      spod_s_nfft = nfft
      spod_s_noverlap = noverlap
      spod_s_dt = dt
      spod_s_nfreq = nfft/2 + 1
      spod_s_nblk = nblk_in
      spod_s_n_snaps = 0
      spod_s_center = .true.

      allocate (spod_s_window(nfft))
      allocate (spod_s_freq(spod_s_nfreq))
      allocate (spod_s_t_idx(spod_s_nblk))
      allocate (spod_s_n_complete(spod_s_nblk))

      do i = 1, spod_s_nfreq
         spod_s_freq(i) = dble(i - 1)/(nfft*dt)
      end do

      ! Staggered block starts
      do i = 1, spod_s_nblk
         spod_s_t_idx(i) = -(i - 1)*(nfft - noverlap)
         spod_s_n_complete(i) = 0
      end do

      allocate (spod_s_dft_re(spod_s_nfreq*spod_s_nblk))
      allocate (spod_s_dft_im(spod_s_nfreq*spod_s_nblk))

      call k_zero(spod_s_mean)
      do idx = 1, spod_s_nfreq*spod_s_nblk
         call k_zero(spod_s_dft_re(idx))
         call k_zero(spod_s_dft_im(idx))
      end do

      spod_s_initialized = .true.

   end subroutine spod_s_init

   ! spod_s_cleanup — Deallocate all streaming SPOD state arrays
   subroutine spod_s_cleanup()
      if (allocated(spod_s_window)) deallocate (spod_s_window)
      if (allocated(spod_s_freq)) deallocate (spod_s_freq)
      if (allocated(spod_s_t_idx)) deallocate (spod_s_t_idx)
      if (allocated(spod_s_n_complete)) &
         deallocate (spod_s_n_complete)
      if (allocated(spod_s_dft_re)) deallocate (spod_s_dft_re)
      if (allocated(spod_s_dft_im)) deallocate (spod_s_dft_im)
      spod_s_initialized = .false.
      spod_s_n_snaps = 0

   end subroutine spod_s_cleanup

end module nekstab_spod_streaming_state

!-----------------------------------------------------------------------
!     Main streaming SPOD module
!-----------------------------------------------------------------------
module nekstab_modal_spod_streaming

   use nekstab_krylov_subspace, only: krylov_vector
   use nekstab_nek_bridge, only: nid
   use nekstab_spod_streaming_state, only: spod_s_center

   implicit none
   private

   public :: spod_streaming_batch

contains

!-----------------------------------------------------------------------
! spod_streaming_batch — Process pre-loaded snapshots via streaming SPOD
!
! Arguments:
!   snaps    [in] -- array of centered snapshots
!   nsnap    [in] -- number of snapshots
!   delta_t  [in] -- sampling interval
!   nfft     [in] -- FFT block length
!   noverlap [in] -- block overlap
!   nsave    [in] -- number of modes to save per frequency
!-----------------------------------------------------------------------
   subroutine spod_streaming_batch(snaps, nsnap, delta_t, &
                                   nfft, noverlap, nsave)

      integer, intent(in) :: nsnap, nfft, noverlap, nsave
      real, intent(in) :: delta_t
      type(krylov_vector), intent(in) :: snaps(nsnap)

      integer :: i

      if (nid == 0) then
         write (6, *) ''
         write (6, *) '  Using STREAMING SPOD algorithm (batch mode)'
         write (6, *) ''
      end if

      call spod_stream_init(nfft, noverlap, delta_t)
      spod_s_center = .false.

      do i = 1, nsnap
         call spod_stream_update(snaps(i), i)
      end do

      call spod_stream_finalize(nsave)

   end subroutine spod_streaming_batch

end module nekstab_modal_spod_streaming

!-----------------------------------------------------------------------
! spod_stream_init — Initialize streaming SPOD computation
!
! Arguments:
!   nfft_in     [in] -- FFT block length
!   noverlap_in [in] -- block overlap
!   delta_t_in  [in] -- sampling interval
!-----------------------------------------------------------------------
subroutine spod_stream_init(nfft_in, noverlap_in, delta_t_in)

   use nekstab_nek_bridge, only: nid
   use nekstab_spod_streaming_state, only: spod_s_init, spod_s_window, spod_s_win_weight
   use nekstab_modal_pod, only: hamming_window

   integer, intent(in) :: nfft_in, noverlap_in
   real, intent(in) :: delta_t_in

   integer :: nblk

   nblk = 10

   if (nid == 0) then
      write (6, *) ''
      write (6, *) '==============================================='
      write (6, *) '  STREAMING SPOD INITIALIZATION'
      write (6, *) '==============================================='
      write (6, '(A,I6)') '    nfft:     ', nfft_in
      write (6, '(A,I6)') '    noverlap: ', noverlap_in
      write (6, '(A,I6)') '    nfreq:    ', nfft_in/2 + 1
      write (6, '(A,I6)') '    nblk:     ', nblk
      write (6, '(A,E12.4)') '    dt:       ', delta_t_in
      write (6, '(A,E12.4)') '    dSt:      ', &
         1.0d0/(nfft_in*delta_t_in)
   end if

   call spod_s_init(nfft_in, noverlap_in, delta_t_in, nblk)
   call hamming_window(nfft_in, spod_s_window, spod_s_win_weight)

   if (nid == 0) then
      write (6, *) '  Streaming SPOD initialized'
      write (6, *) ''
   end if

end subroutine spod_stream_init

!-----------------------------------------------------------------------
! spod_stream_update — Accumulate one snapshot into DFT coefficients
!
! Arguments:
!   snap  [in] -- current snapshot (krylov_vector)
!   isnap [in] -- snapshot index (for progress output)
!-----------------------------------------------------------------------
subroutine spod_stream_update(snap, isnap)

   use nekstab_nek_bridge, only: nid
   use nekstab_krylov_subspace, only: krylov_vector, k_cmult, k_add2s2, k_copy, k_sub2, &
                                      NEKSTAB_PI
   use nekstab_spod_streaming_state, only: spod_s_center, spod_s_dft_im, spod_s_dft_re, &
                                           spod_s_initialized, spod_s_mean, spod_s_n_snaps, &
                                           spod_s_nblk, spod_s_nfft, spod_s_n_complete, &
                                           spod_s_t_idx, spod_s_window, spod_s_win_weight, &
                                           spod_s_work

   type(krylov_vector), intent(in) :: snap
   integer, intent(in) :: isnap

   real :: scale, w_norm, cos_t, sin_t
   integer :: iblk, ifreq, idx, t_local
   real :: base_angle

   if (.not. spod_s_initialized) then
      if (nid == 0) write (6, *) 'nekStab: ERROR: spod_stream_update called before init'
      return
   end if

   spod_s_n_snaps = spod_s_n_snaps + 1

   if (spod_s_center) then
      ! Update running mean: mean = (n*mean + snap) / (n+1)
      scale = dble(spod_s_n_snaps - 1)/dble(spod_s_n_snaps)
      call k_cmult(spod_s_mean, scale)
      scale = 1.0d0/dble(spod_s_n_snaps)
      call k_add2s2(spod_s_mean, snap, scale)

      ! Center snapshot (subtract running mean)
      call k_copy(spod_s_work, snap)
      call k_sub2(spod_s_work, spod_s_mean)
   else
      ! Data already centered by caller
      call k_copy(spod_s_work, snap)
   end if

   ! Accumulate DFT for each active block
   do iblk = 1, spod_s_nblk

      ! Skip blocks not yet active
      if (spod_s_t_idx(iblk) < 0) then
         spod_s_t_idx(iblk) = spod_s_t_idx(iblk) + 1
         cycle
      end if

      t_local = spod_s_t_idx(iblk)

      ! Block complete: mark done and preserve DFT data
      if (t_local >= spod_s_nfft) then
         if (spod_s_n_complete(iblk) == 0) &
            spod_s_n_complete(iblk) = 1
         cycle
      end if

      ! Combined weight: window * normalization
      w_norm = spod_s_window(t_local + 1)*2.0d0 &
               /(dble(spod_s_nfft)*spod_s_win_weight)

      base_angle = 2.0d0*NEKSTAB_PI*dble(t_local) &
                   /dble(spod_s_nfft)

      ! DC component (half-weight for one-sided spectrum)
      idx = iblk
      call k_add2s2(spod_s_dft_re(idx), spod_s_work, &
                    w_norm*0.5d0)

      ! Interior frequencies
      do ifreq = 2, spod_s_nfreq - 1
         cos_t = cos(base_angle*dble(ifreq - 1))*w_norm
         sin_t = -sin(base_angle*dble(ifreq - 1))*w_norm

         idx = (ifreq - 1)*spod_s_nblk + iblk
         call k_add2s2(spod_s_dft_re(idx), spod_s_work, cos_t)
         call k_add2s2(spod_s_dft_im(idx), spod_s_work, sin_t)
      end do

      ! Nyquist component (half-weight)
      ifreq = spod_s_nfreq
      cos_t = cos(base_angle*dble(ifreq - 1)) &
              *w_norm*0.5d0
      sin_t = -sin(base_angle*dble(ifreq - 1)) &
              *w_norm*0.5d0
      idx = (ifreq - 1)*spod_s_nblk + iblk
      call k_add2s2(spod_s_dft_re(idx), spod_s_work, cos_t)
      call k_add2s2(spod_s_dft_im(idx), spod_s_work, sin_t)

      spod_s_t_idx(iblk) = t_local + 1

   end do

   ! Progress output
   if (nid == 0 .and. mod(spod_s_n_snaps, 50) == 0) then
      write (6, '(A,I6,A,I3,A)') '    Processed ', spod_s_n_snaps, &
         ' snapshots, ', sum(spod_s_n_complete), ' blocks done'
   end if

end subroutine spod_stream_update

!-----------------------------------------------------------------------
! spod_stream_finalize — Finalize streaming SPOD (CSD eigensolve)
!
! Purpose:
!   Forms the cross-spectral density matrix at each frequency from
!   the accumulated DFT coefficients, solves the Hermitian eigenproblem,
!   and saves SPOD modes and spectrum to disk.
!
! Arguments:
!   nsave [in] -- number of modes to save per frequency
!-----------------------------------------------------------------------
subroutine spod_stream_finalize(nsave)

   use nekstab_nek_bridge, only: nid, nekStab_dp, NEKSTAB_UNIT_SPOD
   use nekstab_krylov_subspace, only: krylov_vector, k_copy
   use nekstab_krylov_inner_products, only: k_gram_complex
   use nekstab_lapack, only: eig_hermitian
   use nekstab_spod_streaming_state, only: spod_s_cleanup, spod_s_dt, spod_s_dft_im, &
                                           spod_s_dft_re, spod_s_freq, spod_s_initialized, &
                                           spod_s_nblk, spod_s_nfft, spod_s_nfreq, &
                                           spod_s_noverlap, spod_s_n_complete, spod_s_n_snaps
   use nekstab_modal_spod, only: spod_save_modes

   integer, intent(in) :: nsave

   integer :: ifreq, iblk, jblk, idx_i, idx_j, nblk_done, nmodes
   real, allocatable :: spod_evals(:)
   complex(nekStab_dp), allocatable :: CSD(:, :)
   complex(nekStab_dp), allocatable :: spod_evecs(:, :)
   type(krylov_vector), allocatable :: tmp_re(:), tmp_im(:)
   integer :: i, idx

   if (.not. spod_s_initialized) then
      if (nid == 0) write (6, *) 'nekStab: ERROR: spod_stream_finalize called before init'
      return
   end if

   if (nid == 0) then
      write (6, *) ''
      write (6, *) '==============================================='
      write (6, *) '  STREAMING SPOD FINALIZATION'
      write (6, *) '==============================================='
      write (6, '(A,I6)') '    Total snapshots:   ', spod_s_n_snaps
      write (6, '(A,I6)') '    Completed blocks:  ', &
         sum(spod_s_n_complete)
   end if

   ! Count valid blocks
   nblk_done = 0
   do iblk = 1, spod_s_nblk
      if (spod_s_n_complete(iblk) > 0) nblk_done = nblk_done + 1
   end do

   if (nblk_done < 2) then
      if (nid == 0) write (6, *) 'nekStab: ERROR: Need at least 2 complete blocks for streaming SPOD'
      call spod_s_cleanup()
      return
   end if

   if (nid == 0) then
      write (6, '(A,I6)') '    Valid blocks:      ', nblk_done
   end if

   allocate (CSD(nblk_done, nblk_done))
   allocate (spod_evals(nblk_done))
   allocate (spod_evecs(nblk_done, nblk_done))
   allocate (tmp_re(nblk_done), tmp_im(nblk_done))

   if (nid == 0) then
      open (unit=NEKSTAB_UNIT_SPOD, file='spod_stream_spectrum.dat', &
            status='replace')
      write (NEKSTAB_UNIT_SPOD, '(A)') &
         '# Streaming SPOD Spectrum'
      write (NEKSTAB_UNIT_SPOD, '(A,I6,A,I6,A,E12.4)') &
         '# nfft = ', spod_s_nfft, &
         ', noverlap = ', spod_s_noverlap, &
         ', dt = ', spod_s_dt
      write (NEKSTAB_UNIT_SPOD, '(A,I6,A,I6,A,I6)') &
         '# total_snapshots = ', spod_s_n_snaps, &
         ', completed_blocks = ', sum(spod_s_n_complete), &
         ', valid_blocks = ', nblk_done
      write (NEKSTAB_UNIT_SPOD, '(A)', advance='no') '# St'
      do i = 1, nblk_done
         write (NEKSTAB_UNIT_SPOD, '(A,I3)', advance='no') '    lambda_', i
      end do
      write (NEKSTAB_UNIT_SPOD, *)
   end if

   ! Loop over frequencies
   if (nid == 0) write (6, '(A,I4,A)') &
      '    Processing ', spod_s_nfreq, ' frequencies...'

   nmodes = min(nsave, nblk_done)

   do ifreq = 1, spod_s_nfreq

      ! Gather DFT slices into contiguous arrays
      do iblk = 1, nblk_done
         idx = (ifreq - 1)*spod_s_nblk + iblk
         call k_copy(tmp_re(iblk), spod_s_dft_re(idx))
         call k_copy(tmp_im(iblk), spod_s_dft_im(idx))
      end do
      ! Form Hermitian CSD matrix (batch)
      call k_gram_complex(CSD, tmp_re, tmp_im, &
                          nblk_done, nblk_done)
      CSD = CSD/dble(nblk_done)

      call eig_hermitian(CSD, spod_evals, spod_evecs, nblk_done)

      if (nid == 0) then
         write (NEKSTAB_UNIT_SPOD, '(E14.6)', advance='no') &
            spod_s_freq(ifreq)
         do i = 1, nblk_done
            write (NEKSTAB_UNIT_SPOD, '(E14.6)', advance='no') &
               spod_evals(i)
         end do
         write (NEKSTAB_UNIT_SPOD, *)
      end if

      call spod_stream_save_modes(ifreq, spod_s_nfreq, &
                                  spod_s_freq, spod_evecs, spod_evals, &
                                  nblk_done, nmodes)

   end do

   if (nid == 0) then
      close (NEKSTAB_UNIT_SPOD)
      write (6, *) '    Wrote spod_stream_spectrum.dat'
   end if

   deallocate (CSD, spod_evals, spod_evecs)
   deallocate (tmp_re, tmp_im)
   call spod_s_cleanup()

   if (nid == 0) then
      write (6, *) '  Streaming SPOD complete'
      write (6, *) '==============================================='
   end if

end subroutine spod_stream_finalize

!-----------------------------------------------------------------------
! spod_stream_save_modes — Save SPOD modes at one frequency
!
! Purpose:
!   Extracts DFT slices from module state and delegates to
!   spod_save_modes for field reconstruction and output.
!-----------------------------------------------------------------------
subroutine spod_stream_save_modes(ifreq, nfreq, freq, &
                                  evecs, evals, nblk, nsave)

   use nekstab_nek_bridge, only: nekStab_dp
   use nekstab_spod_streaming_state, only: spod_s_dft_im, spod_s_dft_re, spod_s_nblk
   use nekstab_modal_spod, only: spod_save_modes

   integer, intent(in) :: ifreq, nfreq, nblk, nsave
   real, intent(in) :: freq(nfreq), evals(nblk)
   complex(nekStab_dp), intent(in) :: evecs(nblk, nblk)

   integer :: idx_base

   idx_base = (ifreq - 1)*spod_s_nblk + 1

   call spod_save_modes( &
      spod_s_dft_re(idx_base), spod_s_dft_im(idx_base), &
      evecs, evals, nblk, ifreq, nfreq, freq, nsave)

end subroutine spod_stream_save_modes
