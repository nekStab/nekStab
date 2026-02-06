!-----------------------------------------------------------------------
! modal_analysis.f90 -- Modal analysis dispatcher for nekStab (Mode 6)
!
! Purpose:
!   Top-level driver that loads snapshots, computes temporal mean,
!   and dispatches to the enabled modal decomposition methods
!   (POD, DMD, SPOD). Controlled by uparam(1):
!     6.0 = run all enabled, 6.1 = POD only,
!     6.2 = DMD only, 6.3 = SPOD only
!
! Public interface:
!   modal_analysis       -- main dispatcher
!   modal_compute_mean   -- compute ensemble temporal mean
!   modal_subtract_mean  -- subtract mean from all snapshots
!
! Dependencies:
!   krylov_subspace, nekstab_vectors, nekstab_io,
!   modal_pod, modal_dmd, modal_spod, modal_spod_streaming,
!   SIZE, TOTAL
!-----------------------------------------------------------------------

      module nekstab_modal_analysis
         use krylov_subspace
         use nekstab_vectors
         use nekstab_io
         use modal_pod
         use modal_dmd
         use modal_spod
         use modal_spod_streaming
         implicit none
         private
         public :: modal_analysis, modal_compute_mean,
     $             modal_subtract_mean
      contains

!-----------------------------------------------------------------------
! modal_analysis -- Main modal decomposition dispatcher
!
! Purpose:
!   Loads snapshot sequence, subtracts mean, then runs POD, DMD,
!   and/or SPOD depending on user flags (ifpod, ifdmd, ifspod).
!-----------------------------------------------------------------------
      subroutine modal_analysis

         implicit none
         include 'SIZE'
         include 'TOTAL'

         type(krylov_vector), allocatable :: snaps(:)
         type(krylov_vector) :: mean_snap
         integer :: i

!        Print header
         if (nid == 0) then
            write(6,*) ''
            write(6,*) '==============================================='
            write(6,*) '        MODAL ANALYSIS (Mode 6)'
            write(6,*) '==============================================='
            write(6,*) ''
            write(6,'(A,A)')    '  Prefix:     ', trim(modal_prefix)
            write(6,'(A,I6)')   '  Snapshots:  ', modal_nsnap
            write(6,'(A,E12.4)')'  dt:         ', modal_dt
            write(6,'(A,I6)')   '  Modes save: ', modal_nsave
            write(6,'(A,L1)')   '  POD:        ', ifpod
            write(6,'(A,L1)')   '  DMD:        ', ifdmd
            write(6,'(A,L1)')   '  SPOD:       ', ifspod
            if (ifwinamp) then
               write(6,'(A)')   '  Window:     Amplitude norm (PySPOD)'
            else
               write(6,'(A)')   '  Window:     Energy norm (Parseval)'
            end if
            write(6,*) ''
         end if

!        Validate parameters
         if (modal_nsnap < 2) then
            if (nid == 0) write(6,*) 'ERROR: modal_nsnap must be >= 2'
            call nek_end
         end if

         if (.not. ifpod .and. .not. ifdmd .and. .not. ifspod) then
            if (nid == 0) write(6,*) 'WARNING: No methods enabled'
            if (nid == 0) write(6,*)
     $         '  Set ifpod, ifdmd, or ifspod = .true.'
            return
         end if

!        Load snapshots
         if (nid == 0) write(6,*) 'Loading snapshots...'
         allocate(snaps(modal_nsnap))
         call load_files(snaps, modal_nsnap, modal_nsnap, modal_prefix)
         if (nid == 0) then
            write(6,*) '  Loaded', modal_nsnap, 'snapshots'
            call flush(6)
         end if

!        Compute and subtract mean
         if (nid == 0) write(6,*) 'Computing temporal mean...'
         call modal_compute_mean(snaps, modal_nsnap, mean_snap)
         call modal_subtract_mean(snaps, modal_nsnap, mean_snap)

!        Save mean field
         call nopcopy(vx, vy, vz, pr, t,
     $        mean_snap%vx, mean_snap%vy, mean_snap%vz,
     $        mean_snap%pr, mean_snap%t)
         call outpost2(vx, vy, vz, pr, t, 0, 'mea')
         if (nid == 0) write(6,*) '  Saved mean field as mea*'

!        POD (also needed for POD-FFT spectral analysis)
         if (ifpod .or. ifspod) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '-------------------------------------------'
               write(6,*) '  POD (Proper Orthogonal Decomposition)'
               write(6,*) '-------------------------------------------'
            end if
            call pod_compute(snaps, modal_nsnap, modal_nsave)
         end if

!        DMD
         if (ifdmd) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '-------------------------------------------'
               write(6,*) '  DMD (Dynamic Mode Decomposition)'
               write(6,*) '-------------------------------------------'
            end if
            call dmd_compute(snaps, modal_nsnap, modal_dt,
     $                       dmd_rank, modal_nsave)
         end if

!        POD-FFT Spectral Analysis
         if (ifpod .or. ifspod) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '-------------------------------------------'
               write(6,*) '  POD-FFT Spectral Analysis'
               write(6,*) '-------------------------------------------'
            end if
            call pod_fft_spectrum(snaps, modal_nsnap, modal_dt,
     $                            spod_nfft, spod_noverlap)
         end if

!        SPOD (Spectral POD) - Streaming Algorithm
         if (ifspod) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '-------------------------------------------'
               write(6,*) '  SPOD (Spectral POD) - Streaming'
               write(6,*) '-------------------------------------------'
            end if
            call spod_streaming_batch(snaps, modal_nsnap, modal_dt,
     $                        spod_nfft, spod_noverlap, modal_nsave)
         end if

!        Cleanup
         deallocate(snaps)

         if (nid == 0) then
            write(6,*) ''
            write(6,*) '==============================================='
            write(6,*) '  Modal analysis complete'
            write(6,*) '==============================================='
         end if

      end subroutine modal_analysis

!-----------------------------------------------------------------------
! modal_compute_mean -- Compute temporal mean of snapshot ensemble
!
! Arguments:
!   snaps     [in]  -- array of snapshots
!   nsnap     [in]  -- number of snapshots
!   mean_snap [out] -- temporal mean
!-----------------------------------------------------------------------
      subroutine modal_compute_mean(snaps, nsnap, mean_snap)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap
         type(krylov_vector), intent(in) :: snaps(nsnap)
         type(krylov_vector), intent(out) :: mean_snap

         integer :: i
         real :: scale

         call k_zero(mean_snap)
         do i = 1, nsnap
            call k_add2(mean_snap, snaps(i))
         end do

         scale = 1.0d0 / dble(nsnap)
         call k_cmult(mean_snap, scale)

      end subroutine modal_compute_mean

!-----------------------------------------------------------------------
! modal_subtract_mean -- Subtract temporal mean from all snapshots
!
! Arguments:
!   snaps     [inout] -- snapshots to center
!   nsnap     [in]    -- number of snapshots
!   mean_snap [in]    -- temporal mean to subtract
!-----------------------------------------------------------------------
      subroutine modal_subtract_mean(snaps, nsnap, mean_snap)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap
         type(krylov_vector), intent(inout) :: snaps(nsnap)
         type(krylov_vector), intent(in) :: mean_snap

         integer :: i

         do i = 1, nsnap
            call k_sub2(snaps(i), mean_snap)
         end do

      end subroutine modal_subtract_mean

      end module nekstab_modal_analysis
