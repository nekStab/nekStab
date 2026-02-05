!-----------------------------------------------------------------------
! modal_analysis.f90: POD, DMD, and SPOD for nekStab
!
! Purpose:
!   - Proper Orthogonal Decomposition (POD) via method of snapshots
!   - Dynamic Mode Decomposition (DMD) - projected algorithm
!   - Spectral POD (SPOD) - frequency-resolved modes
!
! Usage (Mode 6):
!   uparam(1) = 6.0   ! Run all enabled methods
!   uparam(1) = 6.1   ! POD only
!   uparam(1) = 6.2   ! DMD only
!   uparam(1) = 6.3   ! SPOD only
!
! Set in nekStab_usrchk:
!   modal_prefix  = 'DNS'    ! Snapshot file prefix
!   modal_nsnap   = 500      ! Number of snapshots
!   modal_dt      = 0.1      ! Time between snapshots
!   modal_nsave   = 10       ! Number of modes to save
!   ifpod/ifdmd/ifspod       ! Which methods to run
!
! Author: nekStab team
! Date: 2026
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine modal_analysis
!     Main dispatcher for modal analysis (Mode 6)

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         type(krylov_vector), allocatable :: snaps(:)
         type(krylov_vector) :: mean_snap
         integer :: i

!  -----------------------------------------------------------------
!  Print header
!  -----------------------------------------------------------------
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

!  -----------------------------------------------------------------
!  Validate parameters
!  -----------------------------------------------------------------
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

!  -----------------------------------------------------------------
!  PHASE 1: Load snapshots
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) 'Loading snapshots...'
         allocate(snaps(modal_nsnap))

         call load_files(snaps, modal_nsnap, modal_nsnap, modal_prefix)

         if (nid == 0) then
            write(6,*) '  Loaded', modal_nsnap, 'snapshots'
            call flush(6)
         end if

!  -----------------------------------------------------------------
!  PHASE 2: Compute and subtract mean
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) 'Computing temporal mean...'
         call modal_compute_mean(snaps, modal_nsnap, mean_snap)
         call modal_subtract_mean(snaps, modal_nsnap, mean_snap)

!        Save mean field
         call nopcopy(vx, vy, vz, pr, t,
     $        mean_snap%vx, mean_snap%vy, mean_snap%vz,
     $        mean_snap%pr, mean_snap%t)
         call outpost2(vx, vy, vz, pr, t, 0, 'mea')
         if (nid == 0) write(6,*) '  Saved mean field as mea*'

!  -----------------------------------------------------------------
!  PHASE 3: POD (also needed for POD-FFT spectral analysis)
!  -----------------------------------------------------------------
         if (ifpod .or. ifspod) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '-------------------------------------------'
               write(6,*) '  POD (Proper Orthogonal Decomposition)'
               write(6,*) '-------------------------------------------'
            end if
            call pod_compute(snaps, modal_nsnap, modal_nsave)
         end if

!  -----------------------------------------------------------------
!  PHASE 4: DMD
!  -----------------------------------------------------------------
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

!  -----------------------------------------------------------------
!  PHASE 5: POD-FFT Spectral Analysis
!  -----------------------------------------------------------------
!  Analyze frequency content of POD temporal coefficients.
!  Runs automatically with POD to identify dominant frequencies.
!  -----------------------------------------------------------------
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

!  -----------------------------------------------------------------
!  PHASE 6: SPOD (Spectral POD) - Streaming Algorithm
!  -----------------------------------------------------------------
!  Uses streaming SPOD for efficiency (perfect cache locality).
!  Batch SPOD (spod_compute) is retained for reference but disabled
!  because per-point FFT is ~100x slower than streaming DFT.
!  -----------------------------------------------------------------
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

!  -----------------------------------------------------------------
!  Cleanup
!  -----------------------------------------------------------------
         deallocate(snaps)

         if (nid == 0) then
            write(6,*) ''
            write(6,*) '==============================================='
            write(6,*) '  Modal analysis complete'
            write(6,*) '==============================================='
         end if

      end subroutine modal_analysis

!-----------------------------------------------------------------------
      subroutine modal_compute_mean(snaps, nsnap, mean_snap)
!     Compute temporal mean of snapshot ensemble

         use krylov_subspace
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
      subroutine modal_subtract_mean(snaps, nsnap, mean_snap)
!     Subtract temporal mean from all snapshots

         use krylov_subspace
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

!-----------------------------------------------------------------------
      subroutine pod_compute(snaps, nsnap, nsave)
!     POD via method of snapshots (Sirovich, 1987)
!
!     Forms correlation matrix C(i,j) = <snaps(i), snaps(j)>_E
!     Solves eigenvalue problem C v = lambda v
!     Reconstructs spatial modes Phi_k = Sum_i v_k(i) snaps(i) / sqrt(lambda_k n)

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, nsave
         type(krylov_vector), intent(in) :: snaps(nsnap)

         real, allocatable :: C(:,:), eigvals(:), eigvecs(:,:)
         integer :: i, j

!  -----------------------------------------------------------------
!  Step 1: Form correlation matrix using k_dot
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Forming correlation matrix...'

         allocate(C(nsnap, nsnap))
         allocate(eigvals(nsnap))
         allocate(eigvecs(nsnap, nsnap))

!        Exploit symmetry: only compute upper triangle
         do j = 1, nsnap
            do i = 1, j
               call k_dot(C(i,j), snaps(i), snaps(j))
               C(j,i) = C(i,j)
            end do
         end do

!        Normalize by number of snapshots
         C = C / dble(nsnap)

!  -----------------------------------------------------------------
!  Step 2: Solve eigenvalue problem
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Solving eigenvalue problem...'

         call eig_symmetric(C, eigvals, eigvecs, nsnap)

!  -----------------------------------------------------------------
!  Step 3: Reconstruct and save modes
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Reconstructing POD modes...'

         call pod_reconstruct_modes(snaps, eigvecs, eigvals,
     $                              nsnap, nsave)

!  -----------------------------------------------------------------
!  Step 4: Write eigenvalue spectrum
!  -----------------------------------------------------------------
         call pod_write_spectrum(eigvals, nsnap)

         deallocate(C, eigvals, eigvecs)

         if (nid == 0) write(6,*) '  POD complete'

      end subroutine pod_compute

!-----------------------------------------------------------------------
      subroutine eig_symmetric(A, eigvals, eigvecs, n)
!     Wrapper for LAPACK DSYEV (symmetric eigenvalue problem)
!     Returns eigenvalues in DESCENDING order (largest first)

         implicit none
         integer, intent(in) :: n
         real, intent(in) :: A(n, n)
         real, intent(out) :: eigvals(n), eigvecs(n, n)

         real, allocatable :: work(:), Acopy(:,:)
         real :: tmp_val
         real, allocatable :: tmp_vec(:)
         integer :: lwork, info, i, j

!        Copy input (DSYEV overwrites)
         allocate(Acopy(n,n), tmp_vec(n))
         Acopy = A

!        Query optimal workspace
         allocate(work(1))
         call dsyev('V', 'U', n, Acopy, n, eigvals, work, -1, info)
         lwork = int(work(1))
         deallocate(work)
         allocate(work(lwork))

!        Compute eigenvalues/vectors (returns ascending order)
         call dsyev('V', 'U', n, Acopy, n, eigvals, work, lwork, info)

         if (info /= 0) then
            write(6,*) 'ERROR: DSYEV failed with info =', info
            call nek_end
         end if

         eigvecs = Acopy

!        Reverse to descending order (largest eigenvalue first)
         do i = 1, n/2
            j = n - i + 1

            tmp_val = eigvals(i)
            eigvals(i) = eigvals(j)
            eigvals(j) = tmp_val

            tmp_vec = eigvecs(:, i)
            eigvecs(:, i) = eigvecs(:, j)
            eigvecs(:, j) = tmp_vec
         end do

         deallocate(work, Acopy, tmp_vec)

      end subroutine eig_symmetric

!-----------------------------------------------------------------------
      subroutine pod_reconstruct_modes(snaps, eigvecs, eigvals,
     $                                  nsnap, nsave)
!     Reconstruct POD spatial modes from temporal coefficients
!     Phi_m = Sum_i v_m(i) * snaps(i) / sqrt(lambda_m * nsnap)

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, nsave
         type(krylov_vector), intent(in) :: snaps(nsnap)
         real, intent(in) :: eigvecs(nsnap, nsnap), eigvals(nsnap)

         type(krylov_vector) :: mode
         real :: scale, mode_norm
         integer :: m, i, nmodes
         character(len=3) :: prefix

         prefix = 'pod'
         nmodes = min(nsave, nsnap)

!        Modes are sorted descending by eigenvalue
         do m = 1, nmodes
            call k_zero(mode)

!           Phi_m = Sum_i v_m(i) * snaps(i)
            do i = 1, nsnap
               call k_add2s2(mode, snaps(i), eigvecs(i, m))
            end do

!           Normalize: divide by sqrt(lambda_m * nsnap)
            if (eigvals(m) > 0.0d0) then
               scale = 1.0d0 / sqrt(eigvals(m) * dble(nsnap))
               call k_cmult(mode, scale)
            end if

!           Output mode
            call nopcopy(vx, vy, vz, pr, t,
     $           mode%vx, mode%vy, mode%vz, mode%pr, mode%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix)

!           Compute norm (ALL ranks must call - uses MPI_Allreduce)
            call k_norm(mode_norm, mode)

            if (nid == 0) then
               write(6,'(A,I4,A,E12.4,A,F8.4)')
     $              '    Mode', m, ': lambda =', eigvals(m),
     $              ', ||Phi|| =', mode_norm
            end if
         end do

      end subroutine pod_reconstruct_modes

!-----------------------------------------------------------------------
      subroutine pod_write_spectrum(eigvals, n)
!     Write POD eigenvalue spectrum to file

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: n
         real, intent(in) :: eigvals(n)

         real :: total_energy, cumsum
         integer :: m

         if (nid /= 0) return

         total_energy = sum(eigvals)
         if (total_energy <= 0.0d0) total_energy = 1.0d0

         open(unit=77, file='pod_energy.dat', status='replace')
         write(77, '(A)') '# POD Eigenvalue Spectrum'
         write(77, '(A)')
     $      '# mode   eigenvalue      energy_%    cumulative_%'

         cumsum = 0.0d0
         do m = 1, n
            cumsum = cumsum + eigvals(m)
            write(77, '(I6, 3E16.8)') m, eigvals(m),
     $           100.0d0 * eigvals(m) / total_energy,
     $           100.0d0 * cumsum / total_energy
         end do

         close(77)
         write(6,*) '  Wrote pod_energy.dat'

      end subroutine pod_write_spectrum

!-----------------------------------------------------------------------
      subroutine pod_fft_spectrum(snaps, nsnap, delta_t, nfft, noverlap)
!  =====================================================================
!  POD-FFT SPECTRAL ANALYSIS
!  =====================================================================
!
!  Fast alternative to field-based SPOD. Instead of FFT on every spatial
!  point (expensive!), we FFT the POD temporal coefficients (cheap!).
!
!  Method:
!  -------
!  1. Compute POD via correlation matrix: C(i,j) = <x_i, x_j>
!  2. Eigensolve: C v_k = λ_k v_k
!  3. Temporal coefficients: a_k(t) = sqrt(λ_k * n) * v_k(t)
!  4. FFT each coefficient: â_k(f)
!  5. Power spectrum: P_k(f) = |â_k(f)|^2
!
!  Output:
!  -------
!  pod_fft_spectrum.dat with columns:
!    St, P_1(St), P_2(St), ..., P_r(St)
!
!  This tells us which frequencies are present in each POD mode.
!  For vortex shedding, modes 1-2 should show a peak at St ≈ 0.166.
!
!  =====================================================================

         use krylov_subspace
         use fourier
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, nfft, noverlap
         real, intent(in) :: delta_t
         type(krylov_vector), intent(in) :: snaps(nsnap)

         real, allocatable :: C(:,:), eigvals(:), eigvecs(:,:)
         real, allocatable :: coeffs(:,:)
         real, allocatable :: window(:), power(:,:), freq(:)
         real, allocatable :: time_series(:)
         complex(kind=kind(0.0d0)), allocatable :: spectrum(:)
         real :: win_weight, df, norm_factor
         integer :: i, j, k, m, nblk, nfreq, blk_start, r
         integer :: unit_spec
         real, parameter :: PI_VAL = 3.14159265358979323846d0

!  -----------------------------------------------------------------
!  Step 1: Compute POD correlation matrix and eigenvectors
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Computing POD correlation matrix...'

         allocate(C(nsnap, nsnap))
         allocate(eigvals(nsnap), eigvecs(nsnap, nsnap))

!        Form correlation matrix C(i,j) = <x_i, x_j>
         do j = 1, nsnap
            do i = 1, j
               call k_dot(C(i,j), snaps(i), snaps(j))
               C(j,i) = C(i,j)
            end do
         end do
         C = C / dble(nsnap)

!        Eigensolve (returns descending order)
         call eig_symmetric(C, eigvals, eigvecs, nsnap)

!  -----------------------------------------------------------------
!  Step 2: Determine number of modes to analyze
!  -----------------------------------------------------------------
!        Analyze at least 10 modes, or up to 99.99% energy
         r = nsnap
         do i = 1, nsnap
            if (sum(eigvals(1:i)) / sum(eigvals) > 0.9999d0) then
               r = i
               exit
            end if
         end do
         r = max(r, min(10, nsnap))  ! At least 10 modes
         r = min(r, 30)              ! Cap at 30 for readability

         if (nid == 0) write(6,'(A,I4,A)')
     $        '  Analyzing', r, ' POD modes'

!  -----------------------------------------------------------------
!  Step 3: Extract temporal coefficients
!  -----------------------------------------------------------------
!        a_k(t) = sqrt(λ_k * n) * v_k(t)
!        eigvecs(t, k) is the temporal coefficient for mode k at time t
         allocate(coeffs(nsnap, r))
         do k = 1, r
            if (eigvals(k) > 0) then
               do i = 1, nsnap
                  coeffs(i, k) = eigvecs(i, k) * sqrt(eigvals(k) * nsnap)
               end do
            else
               coeffs(:, k) = 0.0d0
            end if
         end do

         deallocate(C, eigvals, eigvecs)

!  -----------------------------------------------------------------
!  Step 4: Setup FFT parameters
!  -----------------------------------------------------------------
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

!  -----------------------------------------------------------------
!  Step 5: Compute Welch power spectrum for each mode
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Computing power spectra...'

         allocate(window(nfft), freq(nfreq))
         allocate(power(nfreq, r))
         allocate(time_series(nfft), spectrum(nfreq))

!        Hamming window
         call hamming_window(nfft, window, win_weight)

!        Frequency array
         do i = 1, nfreq
            freq(i) = dble(i - 1) * df
         end do

!        Initialize power to zero
         power = 0.0d0

!        Welch method: average over blocks
         do k = 1, r
            do m = 1, nblk
               blk_start = (m - 1) * (nfft - noverlap) + 1

!              Window the time series
               do i = 1, nfft
                  time_series(i) = window(i) * coeffs(blk_start + i - 1, k)
               end do

!              FFT
               call fft_r2c(nfft, time_series, spectrum)

!              Accumulate power (one-sided spectrum)
               do i = 1, nfreq
                  norm_factor = 2.0d0 / (dble(nfft) * win_weight)
                  if (i == 1 .or. i == nfreq) norm_factor = norm_factor / 2.0d0
                  power(i, k) = power(i, k) +
     $                 (real(spectrum(i))**2 + aimag(spectrum(i))**2)
     $                 * norm_factor**2
               end do
            end do

!           Average over blocks
            power(:, k) = power(:, k) / dble(nblk)
         end do

!  -----------------------------------------------------------------
!  Step 6: Write spectrum file
!  -----------------------------------------------------------------
         unit_spec = 79
         if (nid == 0) then
            open(unit=unit_spec, file='pod_fft_spectrum.dat',
     $           status='replace')
            write(unit_spec, '(A)', advance='no') '# St'
            do k = 1, r
               write(unit_spec, '(A,I3)', advance='no') '  P_mode', k
            end do
            write(unit_spec, *)

!           Write data
            do i = 1, nfreq
               write(unit_spec, '(E14.6)', advance='no') freq(i)
               do k = 1, r
                  write(unit_spec, '(E14.6)', advance='no') power(i, k)
               end do
               write(unit_spec, *)
            end do

            close(unit_spec)
            write(6,*) '  Wrote pod_fft_spectrum.dat'

!           Report peak frequencies for first few modes
            write(6,*) ''
            write(6,*) '  Peak frequencies (excluding DC):'
            do k = 1, min(r, 10)
               j = 2  ! Start from index 2 to skip DC
               do i = 3, nfreq
                  if (power(i, k) > power(j, k)) j = i
               end do
               write(6,'(A,I3,A,F8.4,A,E12.4)')
     $              '    Mode', k, ': St =', freq(j),
     $              ', power =', power(j, k)
            end do
         end if

         deallocate(coeffs, window, freq, power, time_series, spectrum)

         if (nid == 0) write(6,*) '  POD-FFT analysis complete'

      end subroutine pod_fft_spectrum

!-----------------------------------------------------------------------
      subroutine dmd_compute(snaps, nsnap, delta_t, rank, nsave)
!  =====================================================================
!  PROJECTED DYNAMIC MODE DECOMPOSITION (DMD)
!  =====================================================================
!
!  References:
!    - Schmid, P.J. (2010) "Dynamic mode decomposition of numerical
!      and experimental data", JFM 656:5-28
!    - Tu et al. (2014) "On dynamic mode decomposition: Theory and
!      applications", J. Comp. Dyn. 1(2):391-421
!
!  Algorithm Overview:
!  -------------------
!  Given snapshots x_1, x_2, ..., x_n at times t, t+dt, ..., t+(n-1)*dt,
!  DMD finds the best-fit linear operator A such that x_{k+1} ≈ A x_k.
!
!  The eigenvalues μ of A (called Ritz values) encode:
!    - Growth/decay rate: σ = log|μ| / dt
!    - Frequency: ω = arg(μ) / dt, or St = ω / (2π)
!
!  For oscillatory dynamics (like vortex shedding), eigenvalues come
!  in complex conjugate pairs: μ = |μ| e^{±iω dt}
!
!  Projected DMD Algorithm:
!  ------------------------
!  1. Define data matrices:
!       X = [x_1, x_2, ..., x_{n-1}]   (input snapshots)
!       Y = [x_2, x_3, ..., x_n]       (time-shifted snapshots)
!     So Y ≈ A X represents the dynamics.
!
!  2. Compute SVD of X via Gram matrix (memory efficient):
!       G = X^T X  (Gram matrix, n×n instead of storing full X)
!       G = V Λ V^T where Λ = diag(σ_1^2, ..., σ_n^2)
!       This gives X ≈ U Σ V^T where U_i = X V_i / σ_i
!
!  3. Project operator onto POD basis:
!       Ã = U^H A U = U^H Y V Σ^{-1}
!
!     In terms of inner products (key formula!):
!       Ã(i,j) = Σ_k Σ_m ⟨x_{k+1}, x_m⟩ V(m,i) V(k,j) / (σ_i σ_j) × σ_i
!
!     IMPORTANT: We need the FULL shifted Gram matrix:
!       G_shift(k,m) = ⟨Y_k, X_m⟩ = ⟨x_{k+1}, x_m⟩  for ALL k,m pairs
!
!     This is NOT symmetric: G_shift(k,m) ≠ G_shift(m,k) in general.
!     The asymmetry is what produces complex eigenvalues!
!
!  4. Eigendecomposition of Ã:
!       Ã W = W Λ_dmd
!     where Λ_dmd contains the DMD eigenvalues (Ritz values).
!
!  5. DMD modes (optional):
!       Φ = X V Σ^{-1} W  (projected DMD modes)
!
!  =====================================================================

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, rank, nsave
         real, intent(in) :: delta_t
         type(krylov_vector), intent(in) :: snaps(nsnap)

         integer :: n, r, i, j, k, m
         real, allocatable :: G(:,:), S(:), Vt(:,:), Sinv(:)
         real, allocatable :: G_shift(:,:)
         real, allocatable :: Atilde(:,:), Atilde_work(:,:)
         complex(kind=kind(0.0d0)), allocatable :: dmd_evals(:)
         complex(kind=kind(0.0d0)), allocatable :: dmd_evecs(:,:)
         real :: alpha, total_energy, cumsum, tol
         real, parameter :: PI_VAL = 3.14159265358979323846d0

!        Number of snapshot pairs for DMD
!        X = [x_1, ..., x_{n-1}], Y = [x_2, ..., x_n]
         n = nsnap - 1

         if (nid == 0) write(6,*) '  Using', n, 'snapshot pairs'

!  -----------------------------------------------------------------
!  STEP 1: Form Gram matrix G = X^T X
!  -----------------------------------------------------------------
!  G(i,j) = ⟨x_i, x_j⟩ for i,j = 1,...,n (using first n snapshots)
!  This is symmetric: G(i,j) = G(j,i), so we only compute upper triangle.
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Computing Gram matrix...'

         allocate(G(n, n))

         do j = 1, n
            do i = 1, j
               call k_dot(G(i,j), snaps(i), snaps(j))
               G(j,i) = G(i,j)  ! Symmetry
            end do
         end do

!  -----------------------------------------------------------------
!  STEP 2: SVD of X via eigendecomposition of Gram matrix
!  -----------------------------------------------------------------
!  G = V Λ V^T where Λ = diag(σ_1^2, ..., σ_n^2)
!  This is equivalent to X = U Σ V^T where:
!    - V are the right singular vectors (temporal modes)
!    - Σ = diag(σ_1, ..., σ_n) are singular values
!    - U = X V Σ^{-1} are left singular vectors (spatial modes)
!
!  We truncate to rank r to reduce computational cost and noise.
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Computing SVD via eigendecomp...'

         allocate(S(n), Vt(n, n))
         call eig_symmetric(G, S, Vt, n)

!        S now contains σ^2 (eigenvalues of G), sorted descending
!        Determine rank r based on energy threshold or user input
         total_energy = sum(S)
         if (rank > 0) then
            r = min(rank, n)
         else
!           Auto-rank: capture 99% of total energy
            cumsum = 0.0d0
            r = n
            do i = 1, n
               cumsum = cumsum + S(i)
               if (cumsum / total_energy > 0.99d0) then
                  r = i
                  exit
               end if
            end do
         end if

         if (nid == 0) write(6,'(A,I4,A,F6.2,A)')
     $        '  Using rank r =', r, ' (',
     $        100.0d0 * sum(S(1:r)) / total_energy, '% energy)'

!        Convert eigenvalues (σ^2) to singular values (σ) and inverses
         allocate(Sinv(r))
         tol = 1.0d-12 * S(1)
         do i = 1, r
            if (S(i) > tol) then
               S(i) = sqrt(S(i))
               Sinv(i) = 1.0d0 / S(i)
            else
               S(i) = 0.0d0
               Sinv(i) = 0.0d0
            end if
         end do

!  -----------------------------------------------------------------
!  STEP 3: Form shifted Gram matrix and projected operator Ã
!  -----------------------------------------------------------------
!  The shifted Gram matrix captures the time-shift dynamics:
!    G_shift(k,m) = ⟨Y_k, X_m⟩ = ⟨x_{k+1}, x_m⟩
!
!  CRITICAL: This is NOT symmetric! G_shift(k,m) ≠ G_shift(m,k)
!  The asymmetry encodes the temporal dynamics and is essential
!  for producing complex eigenvalues in oscillatory systems.
!
!  The projected operator is:
!    Ã(i,j) = [Σ^{-1} V^T G_shift V Σ^{-1}](i,j) × σ_i
!           = Σ_k Σ_m G_shift(k,m) V(m,i) V(k,j) Sinv(i) Sinv(j) × S(i)
!           = Σ_k Σ_m ⟨x_{k+1}, x_m⟩ V(m,i) V(k,j) Sinv(j)
!  -----------------------------------------------------------------
!        Allocate shifted Gram matrix (NOT symmetric!)
         allocate(G_shift(n, n))

!        Compute G_shift(k,m) = ⟨x_{k+1}, x_m⟩ for all k,m = 1,...,n
!        Note: x_{k+1} is snaps(k+1), x_m is snaps(m)
!
!        Optimization: Many elements overlap with G (already computed)
!        G(i,j) = ⟨x_i, x_j⟩, so G_shift(k,m) = G(k+1,m) for k+1 <= n
!        Only need to compute elements involving snaps(n+1) = snaps(nsnap)
!
!        Reuse from G where possible (k+1 ranges 2..n, m ranges 1..n)
         do m = 1, n
            do k = 1, n-1
!              G_shift(k,m) = ⟨x_{k+1}, x_m⟩ = G(k+1, m)
               G_shift(k,m) = G(k+1, m)
            end do
!           k=n requires snaps(n+1) = snaps(nsnap), not in G
            call k_dot(G_shift(n,m), snaps(nsnap), snaps(m))
         end do

         if (nid == 0) write(6,*) '  Forming projected operator...'

!        Form Ã using the shifted Gram matrix
!        Ã = U^H Y V Σ^{-1} = Σ^{-1} V^T (X^T Y) V Σ^{-1}
!        Ã(i,j) = Σ_k Σ_m G_shift(k,m) × V(m,i) × V(k,j) × Sinv(i) × Sinv(j)
!
!        IMPORTANT: Both Sinv(i) and Sinv(j) are needed!
!        - Sinv(i) comes from U_i = X V_i / σ_i
!        - Sinv(j) comes from the explicit Σ^{-1} in U^H Y V Σ^{-1}
         allocate(Atilde(r, r))
         Atilde = 0.0d0

         do j = 1, r
            do i = 1, r
               do k = 1, n
                  do m = 1, n
                     Atilde(i,j) = Atilde(i,j) + G_shift(k,m) *
     $                    Vt(m,i) * Vt(k,j) * Sinv(i) * Sinv(j)
                  end do
               end do
            end do
         end do

!        Note: Ã is generally NOT symmetric, which is correct!
!        For oscillatory systems, Ã will have complex eigenvalue pairs.

!  -----------------------------------------------------------------
!  STEP 4: Eigendecomposition of Ã
!  -----------------------------------------------------------------
!  Solve Ã W = W Λ where Λ = diag(μ_1, ..., μ_r) are DMD eigenvalues.
!
!  For oscillatory dynamics:
!    - Eigenvalues come in complex conjugate pairs: μ, μ*
!    - |μ| < 1: decaying mode, |μ| > 1: growing mode, |μ| = 1: neutral
!    - arg(μ) / dt = angular frequency ω
!    - Strouhal number St = ω / (2π) = arg(μ) / (2π dt)
!  -----------------------------------------------------------------
         if (nid == 0) write(6,*) '  Computing DMD eigenvalues...'

         allocate(dmd_evals(r), dmd_evecs(r, r))
         allocate(Atilde_work(r, r))
         Atilde_work = Atilde

!        General eigenvalue solver (not symmetric - returns complex)
         call eig(Atilde_work, dmd_evecs, dmd_evals, r)

!  -----------------------------------------------------------------
!  STEP 5: Output spectrum and reconstruct modes
!  -----------------------------------------------------------------
         call dmd_write_spectrum(dmd_evals, r, delta_t)
         call dmd_reconstruct_modes(snaps, Vt, S, Sinv,
     $        dmd_evecs, dmd_evals, n, r, nsave)

         deallocate(G, G_shift, S, Vt, Sinv, Atilde, Atilde_work)
         deallocate(dmd_evals, dmd_evecs)

         if (nid == 0) write(6,*) '  DMD complete'

      end subroutine dmd_compute

!-----------------------------------------------------------------------
      subroutine dmd_write_spectrum(evals, r, delta_t)
!     Write DMD eigenvalue spectrum

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: r
         complex(kind=kind(0.0d0)), intent(in) :: evals(r)
         real, intent(in) :: delta_t

         real :: sigma, omega, mu_mag, growth, freq
         real, parameter :: PI_VAL = 3.14159265358979323846d0
         integer :: m

         if (nid /= 0) return

         open(unit=78, file='dmd_spectrum.dat', status='replace')
         write(78, '(A)') '# DMD Eigenvalue Spectrum'
         write(78, '(A)') '# mode   |mu|      sigma       omega       '
     $      // 'St          mu_real     mu_imag'

         do m = 1, r
            mu_mag = abs(evals(m))
            sigma = log(mu_mag) / delta_t
            omega = atan2(aimag(evals(m)), real(evals(m))) / delta_t
            freq = omega / (2.0d0 * PI_VAL)

            write(78, '(I6, 6E14.6)') m, mu_mag, sigma, omega, freq,
     $           real(evals(m)), aimag(evals(m))
         end do

         close(78)
         write(6,*) '  Wrote dmd_spectrum.dat'

      end subroutine dmd_write_spectrum

!-----------------------------------------------------------------------
      subroutine dmd_reconstruct_modes(snaps, V, S, Sinv,
     $     evecs, evals, n, r, nsave)
!     Reconstruct DMD modes: Phi = X V Sum^{-1} W (projected DMD modes)
!     Output real and imaginary parts separately

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: n, r, nsave
         type(krylov_vector), intent(in) :: snaps(n)
         real, intent(in) :: V(n, n), S(r), Sinv(r)
         complex(kind=kind(0.0d0)), intent(in) :: evecs(r, r), evals(r)

         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm, mu_mag, freq
         integer :: m, i, j, nmodes
         real, parameter :: PI_VAL = 3.14159265358979323846d0
         character(len=3) :: prefix_re, prefix_im

         prefix_re = 'dm1'
         prefix_im = 'dm2'
         nmodes = min(nsave, r)

!        DMD mode m: Phi_m = X V Sum^{-1} w_m
!        where w_m is eigenvector of Atilde
         if (nid == 0) write(6,'(A,I4,A)')
     $      '   Reconstructing ', nmodes, ' DMD modes...'

         do m = 1, nmodes
            call k_zero(mode_re)
            call k_zero(mode_im)

!           Phi_m = Sum_i snaps(i) * (Sum_j V(i,j) Sinv(j) w_m(j))
            do i = 1, n
               coef_re = 0.0d0
               coef_im = 0.0d0
               do j = 1, r
                  coef_re = coef_re + V(i,j) * Sinv(j) *
     $                 real(evecs(j,m))
                  coef_im = coef_im + V(i,j) * Sinv(j) *
     $                 aimag(evecs(j,m))
               end do
               call k_add2s2(mode_re, snaps(i), coef_re)
               call k_add2s2(mode_im, snaps(i), coef_im)
            end do

!           Compute mode norm (ALL ranks must call - uses MPI_Allreduce)
            call k_norm(mode_norm, mode_re)

!           Skip output if mode is invalid (NaN or very large)
            if (mode_norm /= mode_norm .or. mode_norm > 1.0d30) then
               if (nid == 0) write(6,*) '  WARNING: DMD mode', m,
     $              'has invalid values, skipping'
               cycle
            end if

!           Output real part
            call nopcopy(vx, vy, vz, pr, t,
     $           mode_re%vx, mode_re%vy, mode_re%vz,
     $           mode_re%pr, mode_re%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix_re)

!           Output imaginary part
            call nopcopy(vx, vy, vz, pr, t,
     $           mode_im%vx, mode_im%vy, mode_im%vz,
     $           mode_im%pr, mode_im%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix_im)

!           Print mode info
            if (nid == 0) then
               mu_mag = abs(evals(m))
               freq = atan2(aimag(evals(m)), real(evals(m))) /
     $              (2.0d0 * PI_VAL * modal_dt)
               write(6,'(A,I4,A,F8.4,A,F10.4,A,E12.4)')
     $              '    Mode', m, ': |mu| =', mu_mag,
     $              ', St =', freq, ', ||Phi_re|| =', mode_norm
            end if
         end do

      end subroutine dmd_reconstruct_modes

!-----------------------------------------------------------------------
      subroutine spod_streaming_batch(snaps, nsnap, delta_t,
     $                                nfft, noverlap, nsave)
!  =====================================================================
!  STREAMING SPOD (BATCH MODE)
!  =====================================================================
!
!  Process pre-loaded snapshots using streaming SPOD algorithm.
!  This is faster than the per-point FFT approach in spod_compute()
!  because it has perfect cache locality (O(nfreq × nblk) k_add2s2
!  operations per snapshot vs per-point FFT with terrible stride).
!
!  Usage: Can replace spod_compute() in modal_analysis() dispatcher.
!
!  =====================================================================

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, nfft, noverlap, nsave
         real, intent(in) :: delta_t
         type(krylov_vector), intent(in) :: snaps(nsnap)

         integer :: i

         if (nid == 0) then
            write(6,*) ''
            write(6,*) '  Using STREAMING SPOD algorithm (batch mode)'
            write(6,*) ''
         end if

!        Initialize streaming SPOD
         call spod_stream_init(nfft, noverlap, delta_t)

!        Process each snapshot
         do i = 1, nsnap
            call spod_stream_update(snaps(i), i)
         end do

!        Finalize and output results
         call spod_stream_finalize(nsave)

      end subroutine spod_streaming_batch

!-----------------------------------------------------------------------
      subroutine spod_compute(snaps, nsnap, delta_t, nfft, noverlap, nsave)
!     Spectral POD (Towne et al., 2018)
!
!     1. Divide snapshots into overlapping blocks
!     2. Apply window and FFT each block
!     3. For each frequency, form cross-spectral density matrix
!     4. Eigendecomposition gives SPOD modes at that frequency

         use krylov_subspace
         use fourier
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
         integer :: unit_spec
         real, parameter :: PI_VAL = 3.14159265358979323846d0

         if (nid == 0) then
            write(6,*) '  SPOD parameters:'
            write(6,'(A,I6)') '    nfft:     ', nfft
            write(6,'(A,I6)') '    noverlap: ', noverlap
         end if

!        Compute number of blocks
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

!        Allocate arrays
         allocate(window(nfft), freq(nfreq))
         allocate(blk_re(nblk), blk_im(nblk))
         allocate(CSD(nblk, nblk))
         allocate(spod_evals(nblk), spod_evecs(nblk, nblk))

!        Compute Hamming window
         call hamming_window(nfft, window, win_weight)

!        Compute frequency array
         call fft_frequencies(nfft, delta_t, freq)

!        Open spectrum file
         unit_spec = 79
         if (nid == 0) then
            open(unit=unit_spec, file='spod_spectrum.dat',
     $           status='replace')
            write(unit_spec, '(A)')
     $         '# SPOD Spectrum: St, eigenvalues (nblk columns)'
         end if

!  -----------------------------------------------------------------
!  Loop over frequencies
!  -----------------------------------------------------------------
         if (nid == 0) write(6,'(A,I4,A)')
     $      '   Processing ', nfreq, ' frequencies...'

         do ifreq = 1, nfreq

            if (nid == 0 .and. mod(ifreq,5) == 1) then
               write(6,'(A,I4,A,I4)') '     freq ', ifreq, ' / ', nfreq
            end if

!           FFT each block at this frequency
            do iblk = 1, nblk
               blk_start = (iblk - 1) * (nfft - noverlap) + 1
               call fft_block_at_freq(snaps, blk_start, nfft,
     $              window, win_weight, delta_t, ifreq,
     $              blk_re(iblk), blk_im(iblk))
            end do

!           Form CSD matrix (Hermitian)
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

!           Eigensolve
            call eig_hermitian(CSD, spod_evals, spod_evecs, nblk)

!           Write spectrum row
            if (nid == 0) then
               write(unit_spec, '(E14.6)', advance='no') freq(ifreq)
               do i = 1, nblk
                  write(unit_spec, '(E14.6)', advance='no')
     $                 spod_evals(i)
               end do
               write(unit_spec, *)
            end if

!           Save modes at selected frequencies
            call spod_save_modes(blk_re, blk_im, spod_evecs,
     $           spod_evals, nblk, ifreq, nfreq, freq, nsave)

         end do

         if (nid == 0) then
            close(unit_spec)
            write(6,*) '  Wrote spod_spectrum.dat'
         end if

         deallocate(window, freq, blk_re, blk_im)
         deallocate(CSD, spod_evals, spod_evecs)

         if (nid == 0) write(6,*) '  SPOD complete'

      end subroutine spod_compute

!-----------------------------------------------------------------------
      subroutine fft_block_at_freq(snaps, blk_start, nfft,
     $     window, win_weight, delta_t, ifreq, out_re, out_im)
!     FFT a block of snapshots and extract frequency component ifreq
!
!     Uses per-point FFT for efficiency (O(nfft log nfft) vs O(nfft²))
!     Point loop is internal implementation detail, not exposed to caller

         use krylov_subspace
         use fourier
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

!        Normalization for one-sided spectrum
         norm_factor = 2.0d0 / (dble(nfft) * win_weight)
         if (ifreq == 1 .or. ifreq == nfreq) then
            norm_factor = norm_factor / 2.0d0  ! DC and Nyquist
         end if

         call k_zero(out_re)
         call k_zero(out_im)

         nv = nx1*ny1*nz1*nelv
         nt = nx1*ny1*nz1*nelt

!        Process velocity components
         do ipt = 1, nv
            do j = 1, nfft
               time_series(j) = window(j) *
     $              snaps(blk_start + j - 1)%vx(ipt)
            end do
            call fft_r2c(nfft, time_series, spectrum)
            out_re%vx(ipt) = real(spectrum(ifreq)) * norm_factor
            out_im%vx(ipt) = aimag(spectrum(ifreq)) * norm_factor
         end do

         do ipt = 1, nv
            do j = 1, nfft
               time_series(j) = window(j) *
     $              snaps(blk_start + j - 1)%vy(ipt)
            end do
            call fft_r2c(nfft, time_series, spectrum)
            out_re%vy(ipt) = real(spectrum(ifreq)) * norm_factor
            out_im%vy(ipt) = aimag(spectrum(ifreq)) * norm_factor
         end do

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

!        Process temperature (if thermal)
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
      subroutine k_dot_complex(alpha, p_re, p_im, q_re, q_im)
!     Hermitian inner product: <p, q> = p^H * W * q
!     where p = p_re + i*p_im, q = q_re + i*q_im
!     <p,q> = <p_re,q_re> + <p_im,q_im> + i*(<p_re,q_im> - <p_im,q_re>)

         use krylov_subspace
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
      subroutine eig_hermitian(A, eigvals, eigvecs, n)
!     Wrapper for LAPACK ZHEEV (complex Hermitian eigenvalue problem)
!     Returns eigenvalues in DESCENDING order

         implicit none
         integer, intent(in) :: n
         complex(kind=kind(0.0d0)), intent(in) :: A(n, n)
         real, intent(out) :: eigvals(n)
         complex(kind=kind(0.0d0)), intent(out) :: eigvecs(n, n)

         complex(kind=kind(0.0d0)), allocatable :: work(:), Acopy(:,:)
         real, allocatable :: rwork(:)
         complex(kind=kind(0.0d0)), allocatable :: tmp_vec(:)
         real :: tmp_val
         integer :: lwork, info, i, j

         allocate(Acopy(n,n), rwork(max(1, 3*n-2)), tmp_vec(n))
         Acopy = A

!        Query optimal workspace
         allocate(work(1))
         call zheev('V', 'U', n, Acopy, n, eigvals,
     $        work, -1, rwork, info)
         lwork = int(real(work(1)))
         deallocate(work)
         allocate(work(max(1, lwork)))

!        Compute eigenvalues/vectors (ascending order)
         call zheev('V', 'U', n, Acopy, n, eigvals,
     $        work, lwork, rwork, info)

         if (info /= 0) then
            write(6,*) 'ERROR: ZHEEV failed with info =', info
            call nek_end
         end if

         eigvecs = Acopy

!        Reverse to descending order
         do i = 1, n/2
            j = n - i + 1

            tmp_val = eigvals(i)
            eigvals(i) = eigvals(j)
            eigvals(j) = tmp_val

            tmp_vec = eigvecs(:, i)
            eigvecs(:, i) = eigvecs(:, j)
            eigvecs(:, j) = tmp_vec
         end do

         deallocate(work, rwork, Acopy, tmp_vec)

      end subroutine eig_hermitian

!-----------------------------------------------------------------------
      subroutine spod_save_modes(blk_re, blk_im, evecs, evals,
     $     nblk, ifreq, nfreq, freq, nsave)
!     Save SPOD modes at selected frequencies
!     Only saves at frequencies with significant energy

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nblk, ifreq, nfreq, nsave
         type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)
         complex(kind=kind(0.0d0)), intent(in) :: evecs(nblk, nblk)
         real, intent(in) :: evals(nblk), freq(nfreq)

         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm
         integer :: m, i, nmodes, ifreq_save
         character(len=3) :: prefix_re, prefix_im
         logical :: save_this_freq

!        Decide which frequencies to save
!        Strategy: save every N-th frequency, plus peak frequencies
         save_this_freq = .false.

!        Save DC (ifreq=1), Nyquist, and every 8th frequency
         if (ifreq == 1) save_this_freq = .true.
         if (ifreq == nfreq) save_this_freq = .true.
         if (mod(ifreq-1, max(1, nfreq/8)) == 0) save_this_freq = .true.

         if (.not. save_this_freq) return

         prefix_re = 'sRe'
         prefix_im = 'sIm'
         nmodes = min(nsave, nblk)

!        Reconstruct and save modes
         do m = 1, nmodes
            call k_zero(mode_re)
            call k_zero(mode_im)

!           Phi_m = Sum_i w_m(i) * Qhat_k(i) / sqrt(lambda_m * nblk)
            do i = 1, nblk
               coef_re = real(evecs(i, m))
               coef_im = aimag(evecs(i, m))

!              mode += coef * blk (complex multiplication)
!              (a + ib)(c + id) = (ac - bd) + i(ad + bc)
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
            call outpost2(vx, vy, vz, pr, t, 0, prefix_re)

!           Output imaginary part
            call nopcopy(vx, vy, vz, pr, t,
     $           mode_im%vx, mode_im%vy, mode_im%vz,
     $           mode_im%pr, mode_im%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix_im)

            if (nid == 0 .and. m == 1) then
               call k_norm(mode_norm, mode_re)
               write(6,'(A,F10.4,A,E12.4,A,E12.4)')
     $              '    St =', freq(ifreq), ': lambda_1 =', evals(1),
     $              ', ||Phi|| =', mode_norm
            end if
         end do

      end subroutine spod_save_modes

!-----------------------------------------------------------------------
      subroutine hamming_window(n, window, win_weight)
!     Compute Hamming window and normalization weight
!
!     Usage: norm_factor = 2 / (nfft * win_weight)  for one-sided spectrum
!
!     Normalization convention controlled by ifwinamp:
!       ifwinamp = .true.   -> Amplitude: win_weight = mean(window) [PySPOD]
!       ifwinamp = .false.  -> Energy: win_weight = sqrt(sum(w²)/n) [Parseval]
!
!     For Hamming window:
!       - mean(w) ≈ 0.54 (amplitude-preserving normalization)
!       - sqrt(sum(w²)/n) ≈ 0.63 (energy-preserving normalization)
!       - Power spectrum ratio: (0.63/0.54)² ≈ 1.36
!
!     PySPOD compatible: `ifwinamp = .true.` (default)

         implicit none
         include 'SIZE'
!        Note: NEKSTAB variables (ifwinamp) available via SIZE -> NEKSTAB.inc

         integer, intent(in) :: n
         real, intent(out) :: window(n), win_weight
         real, parameter :: PI_VAL = 3.14159265358979323846d0
         real :: win_mean, win_energy
         integer :: i

         do i = 1, n
            window(i) = 0.54d0 - 0.46d0 * cos(2.0d0 * PI_VAL *
     $           dble(i-1) / dble(n-1))
         end do

!        Compute window statistics
         win_mean = sum(window) / dble(n)           ! ≈ 0.54 for Hamming
         win_energy = sum(window**2) / dble(n)      ! ≈ 0.397 for Hamming

!        Select win_weight based on normalization convention
!        Formula used: norm_factor = 2 / (nfft * win_weight)
         if (ifwinamp) then
!           Amplitude normalization (PySPOD compatible)
!           A pure sinusoid of amplitude A produces FFT peak of amplitude A
            win_weight = win_mean
         else
!           Energy normalization (Parseval's theorem)
!           Total energy in time equals total energy in frequency domain
            win_weight = sqrt(win_energy)
         end if

      end subroutine hamming_window

!-----------------------------------------------------------------------
! STREAMING SPOD MODULE
!-----------------------------------------------------------------------
! Persistent state for streaming (online) SPOD computation.
! Processes one snapshot at a time, accumulating DFT coefficients.
!
! Advantages over batch SPOD:
!   - No need to store all snapshots in memory
!   - Perfect cache locality (one snapshot at a time)
!   - Can run online during DNS simulation
!
! Usage:
!   call spod_stream_init(nfft, noverlap, dt)
!   do istep = 1, nsteps
!      ... DNS step produces snapshot ...
!      call spod_stream_update(snap, istep)
!   end do
!   call spod_stream_finalize(nsave)
!-----------------------------------------------------------------------
      module spod_streaming_state
         use krylov_subspace
         implicit none
         private

         public :: spod_s_init, spod_s_cleanup
         public :: spod_s_nfft, spod_s_noverlap, spod_s_nfreq, spod_s_nblk
         public :: spod_s_dt, spod_s_win_weight
         public :: spod_s_window, spod_s_t_idx, spod_s_n_complete
         public :: spod_s_n_snaps, spod_s_initialized
         public :: spod_s_mean
         public :: spod_s_dft_re, spod_s_dft_im

!        Scalar parameters (set during init)
         integer, save :: spod_s_nfft = 0
         integer, save :: spod_s_noverlap = 0
         integer, save :: spod_s_nfreq = 0
         integer, save :: spod_s_nblk = 0
         real, save :: spod_s_dt = 0.0d0
         real, save :: spod_s_win_weight = 0.0d0
         integer, save :: spod_s_n_snaps = 0
         logical, save :: spod_s_initialized = .false.

!        Allocatable arrays
         real, allocatable, save :: spod_s_window(:)
         integer, allocatable, save :: spod_s_t_idx(:)
         integer, allocatable, save :: spod_s_n_complete(:)

!        Mean snapshot (running average)
         type(krylov_vector), save :: spod_s_mean

!        DFT accumulators: dft_re(ifreq, iblk), dft_im(ifreq, iblk)
!        Stored as 1D array indexed as (ifreq-1)*nblk + iblk
         type(krylov_vector), allocatable, save :: spod_s_dft_re(:)
         type(krylov_vector), allocatable, save :: spod_s_dft_im(:)

      contains

         subroutine spod_s_init(nfft, noverlap, dt, nblk_in)
            integer, intent(in) :: nfft, noverlap, nblk_in
            real, intent(in) :: dt
            integer :: i, idx

            spod_s_nfft = nfft
            spod_s_noverlap = noverlap
            spod_s_dt = dt
            spod_s_nfreq = nfft / 2 + 1
            spod_s_nblk = nblk_in
            spod_s_n_snaps = 0

!           Allocate window
            allocate(spod_s_window(nfft))
            allocate(spod_s_t_idx(spod_s_nblk))
            allocate(spod_s_n_complete(spod_s_nblk))

!           Initialize block indices with staggered starts
!           Block 1 starts at 0, block 2 at -(nfft-noverlap), etc.
            do i = 1, spod_s_nblk
               spod_s_t_idx(i) = -(i-1) * (nfft - noverlap)
               spod_s_n_complete(i) = 0
            end do

!           Allocate DFT accumulators
            allocate(spod_s_dft_re(spod_s_nfreq * spod_s_nblk))
            allocate(spod_s_dft_im(spod_s_nfreq * spod_s_nblk))

!           Initialize to zero
            call k_zero(spod_s_mean)
            do idx = 1, spod_s_nfreq * spod_s_nblk
               call k_zero(spod_s_dft_re(idx))
               call k_zero(spod_s_dft_im(idx))
            end do

            spod_s_initialized = .true.

         end subroutine spod_s_init

         subroutine spod_s_cleanup()
            if (allocated(spod_s_window)) deallocate(spod_s_window)
            if (allocated(spod_s_t_idx)) deallocate(spod_s_t_idx)
            if (allocated(spod_s_n_complete)) deallocate(spod_s_n_complete)
            if (allocated(spod_s_dft_re)) deallocate(spod_s_dft_re)
            if (allocated(spod_s_dft_im)) deallocate(spod_s_dft_im)
            spod_s_initialized = .false.
            spod_s_n_snaps = 0
         end subroutine spod_s_cleanup

      end module spod_streaming_state

!-----------------------------------------------------------------------
      subroutine spod_stream_init(nfft_in, noverlap_in, delta_t_in)
!     Initialize streaming SPOD computation
!
!     Parameters:
!       nfft_in     - FFT block size (e.g., 64)
!       noverlap_in - Block overlap (e.g., 32 for 50%)
!       delta_t_in  - Time between snapshots
!
!     Call this once before starting to process snapshots.

         use spod_streaming_state
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nfft_in, noverlap_in
         real, intent(in) :: delta_t_in

         integer :: nblk, i
         real :: win_mean, win_energy
         real, parameter :: PI_VAL = 3.14159265358979323846d0

!        Estimate number of blocks (will grow if more snapshots arrive)
!        Start with a reasonable number
         nblk = 10

         if (nid == 0) then
            write(6,*) ''
            write(6,*) '==============================================='
            write(6,*) '  STREAMING SPOD INITIALIZATION'
            write(6,*) '==============================================='
            write(6,'(A,I6)')    '    nfft:     ', nfft_in
            write(6,'(A,I6)')    '    noverlap: ', noverlap_in
            write(6,'(A,I6)')    '    nfreq:    ', nfft_in/2+1
            write(6,'(A,I6)')    '    nblk:     ', nblk
            write(6,'(A,E12.4)') '    dt:       ', delta_t_in
            write(6,'(A,E12.4)') '    dSt:      ',
     $           1.0d0/(nfft_in*delta_t_in)
         end if

!        Initialize state module
         call spod_s_init(nfft_in, noverlap_in, delta_t_in, nblk)

!        Compute Hamming window
         do i = 1, nfft_in
            spod_s_window(i) = 0.54d0 - 0.46d0 * cos(2.0d0 * PI_VAL *
     $           dble(i-1) / dble(nfft_in-1))
         end do

!        Window normalization (amplitude-preserving for PySPOD compat)
         win_mean = sum(spod_s_window) / dble(nfft_in)
         win_energy = sum(spod_s_window**2) / dble(nfft_in)

         if (ifwinamp) then
            spod_s_win_weight = win_mean
         else
            spod_s_win_weight = sqrt(win_energy)
         end if

         if (nid == 0) then
            write(6,*) '  Streaming SPOD initialized'
            write(6,*) '  Call spod_stream_update() for each snapshot'
            write(6,*) ''
         end if

      end subroutine spod_stream_init

!-----------------------------------------------------------------------
      subroutine spod_stream_update(snap, isnap)
!     Process one snapshot for streaming SPOD
!
!     Parameters:
!       snap  - Current snapshot (krylov_vector)
!       isnap - Current snapshot index (1-based)
!
!     This accumulates DFT coefficients for all overlapping blocks.
!     Call this for each snapshot during DNS or post-processing.

         use krylov_subspace
         use spod_streaming_state
         implicit none
         include 'SIZE'
         include 'TOTAL'

         type(krylov_vector), intent(in) :: snap
         integer, intent(in) :: isnap

         type(krylov_vector) :: x_centered
         real :: scale, theta, w, cos_t, sin_t, norm_factor
         integer :: iblk, ifreq, idx, t_local
         real, parameter :: PI_VAL = 3.14159265358979323846d0

         if (.not. spod_s_initialized) then
            if (nid == 0) write(6,*)
     $         'ERROR: spod_stream_update called before init'
            return
         end if

!  -----------------------------------------------------------------
!  Step 1: Update running mean
!  -----------------------------------------------------------------
!        mean = (n * mean + snap) / (n + 1)
         spod_s_n_snaps = spod_s_n_snaps + 1
         scale = dble(spod_s_n_snaps - 1) / dble(spod_s_n_snaps)
         call k_cmult(spod_s_mean, scale)
         scale = 1.0d0 / dble(spod_s_n_snaps)
         call k_add2s2(spod_s_mean, snap, scale)

!  -----------------------------------------------------------------
!  Step 2: Center snapshot (subtract running mean)
!  -----------------------------------------------------------------
         call k_copy(x_centered, snap)
         call k_sub2(x_centered, spod_s_mean)

!  -----------------------------------------------------------------
!  Step 3: Accumulate DFT for each block
!  -----------------------------------------------------------------
!        Normalization for one-sided spectrum
         norm_factor = 2.0d0 / (dble(spod_s_nfft) * spod_s_win_weight)

         do iblk = 1, spod_s_nblk

!           Check if this block is active (t_idx >= 0)
            if (spod_s_t_idx(iblk) < 0) then
!              Block not yet started
               spod_s_t_idx(iblk) = spod_s_t_idx(iblk) + 1
               cycle
            end if

            t_local = spod_s_t_idx(iblk)

!           Check if block is complete
            if (t_local >= spod_s_nfft) then
!              Block complete - reset with overlap
               spod_s_n_complete(iblk) = spod_s_n_complete(iblk) + 1

!              Reset DFT accumulators for this block
               do ifreq = 1, spod_s_nfreq
                  idx = (ifreq - 1) * spod_s_nblk + iblk
                  call k_zero(spod_s_dft_re(idx))
                  call k_zero(spod_s_dft_im(idx))
               end do

!              Reset time index (keep overlap worth of samples)
               spod_s_t_idx(iblk) = 0
               t_local = 0
            end if

!           Get window coefficient
            w = spod_s_window(t_local + 1)

!           Accumulate DFT for each frequency
!           X̂(f) += w * x * exp(-2πi(f-1)*t/nfft)
!                 = w * x * (cos(θ) - i*sin(θ))
            do ifreq = 1, spod_s_nfreq
               theta = 2.0d0 * PI_VAL * dble(ifreq - 1) *
     $              dble(t_local) / dble(spod_s_nfft)
               cos_t = cos(theta) * w * norm_factor
               sin_t = -sin(theta) * w * norm_factor

!              DC and Nyquist corrections
               if (ifreq == 1 .or. ifreq == spod_s_nfreq) then
                  cos_t = cos_t / 2.0d0
                  sin_t = sin_t / 2.0d0
               end if

               idx = (ifreq - 1) * spod_s_nblk + iblk
               call k_add2s2(spod_s_dft_re(idx), x_centered, cos_t)
               call k_add2s2(spod_s_dft_im(idx), x_centered, sin_t)
            end do

!           Advance time index for this block
            spod_s_t_idx(iblk) = t_local + 1

         end do

!        Progress output
         if (nid == 0 .and. mod(spod_s_n_snaps, 50) == 0) then
            write(6,'(A,I6,A,I3,A)') '    Processed ', spod_s_n_snaps,
     $           ' snapshots, ', sum(spod_s_n_complete), ' blocks done'
         end if

      end subroutine spod_stream_update

!-----------------------------------------------------------------------
      subroutine spod_stream_finalize(nsave)
!     Finalize streaming SPOD: compute modes and write output
!
!     Parameters:
!       nsave - Number of modes to save per frequency
!
!     This computes the cross-spectral density matrix at each frequency,
!     performs eigendecomposition, and outputs SPOD modes.

         use krylov_subspace
         use spod_streaming_state
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsave

         integer :: ifreq, iblk, jblk, idx_i, idx_j, nblk_done, nmodes
         real :: df
         real, allocatable :: freq(:), spod_evals(:)
         complex(kind=kind(0.0d0)), allocatable :: CSD(:,:)
         complex(kind=kind(0.0d0)), allocatable :: spod_evecs(:,:)
         complex(kind=kind(0.0d0)) :: cval
         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm
         integer :: unit_spec, m, i
         character(len=3) :: prefix_re, prefix_im
         real, parameter :: PI_VAL = 3.14159265358979323846d0

         if (.not. spod_s_initialized) then
            if (nid == 0) write(6,*)
     $         'ERROR: spod_stream_finalize called before init'
            return
         end if

         if (nid == 0) then
            write(6,*) ''
            write(6,*) '==============================================='
            write(6,*) '  STREAMING SPOD FINALIZATION'
            write(6,*) '==============================================='
            write(6,'(A,I6)') '    Total snapshots:   ', spod_s_n_snaps
            write(6,'(A,I6)') '    Completed blocks:  ',
     $           sum(spod_s_n_complete)
         end if

!        Count valid blocks (those with at least one complete cycle)
         nblk_done = 0
         do iblk = 1, spod_s_nblk
            if (spod_s_n_complete(iblk) > 0) nblk_done = nblk_done + 1
         end do

         if (nblk_done < 2) then
            if (nid == 0) then
               write(6,*) '  ERROR: Need at least 2 complete blocks'
               write(6,*) '    Have:', nblk_done
               write(6,*) '    Need more snapshots or smaller nfft'
            end if
            call spod_s_cleanup()
            return
         end if

         if (nid == 0) then
            write(6,'(A,I6)') '    Valid blocks:      ', nblk_done
         end if

!        Allocate working arrays
         df = 1.0d0 / (spod_s_nfft * spod_s_dt)
         allocate(freq(spod_s_nfreq))
         allocate(CSD(nblk_done, nblk_done))
         allocate(spod_evals(nblk_done))
         allocate(spod_evecs(nblk_done, nblk_done))

!        Frequency array
         do i = 1, spod_s_nfreq
            freq(i) = dble(i - 1) * df
         end do

!        Open spectrum file
         unit_spec = 79
         prefix_re = 'sRe'
         prefix_im = 'sIm'

         if (nid == 0) then
            open(unit=unit_spec, file='spod_stream_spectrum.dat',
     $           status='replace')
            write(unit_spec, '(A)')
     $         '# Streaming SPOD Spectrum: St, eigenvalues'
         end if

!  -----------------------------------------------------------------
!  Loop over frequencies
!  -----------------------------------------------------------------
         if (nid == 0) write(6,'(A,I4,A)')
     $      '    Processing ', spod_s_nfreq, ' frequencies...'

         nmodes = min(nsave, nblk_done)

         do ifreq = 1, spod_s_nfreq

!           Form CSD matrix (Hermitian): CSD(i,j) = <Qhat_i, Qhat_j>
            do jblk = 1, nblk_done
               do iblk = 1, jblk
                  idx_i = (ifreq - 1) * spod_s_nblk + iblk
                  idx_j = (ifreq - 1) * spod_s_nblk + jblk

                  call k_dot_complex(cval,
     $                 spod_s_dft_re(idx_i), spod_s_dft_im(idx_i),
     $                 spod_s_dft_re(idx_j), spod_s_dft_im(idx_j))

                  CSD(iblk, jblk) = cval
                  CSD(jblk, iblk) = conjg(cval)
               end do
            end do
            CSD = CSD / dble(nblk_done)

!           Eigendecomposition
            call eig_hermitian(CSD, spod_evals, spod_evecs, nblk_done)

!           Write spectrum row
            if (nid == 0) then
               write(unit_spec, '(E14.6)', advance='no') freq(ifreq)
               do i = 1, nblk_done
                  write(unit_spec, '(E14.6)', advance='no')
     $                 spod_evals(i)
               end do
               write(unit_spec, *)
            end if

!           Save modes at selected frequencies
            call spod_stream_save_modes(ifreq, spod_s_nfreq, freq,
     $           spod_evecs, spod_evals, nblk_done, nmodes)

         end do

         if (nid == 0) then
            close(unit_spec)
            write(6,*) '    Wrote spod_stream_spectrum.dat'
         end if

!        Cleanup
         deallocate(freq, CSD, spod_evals, spod_evecs)
         call spod_s_cleanup()

         if (nid == 0) then
            write(6,*) '  Streaming SPOD complete'
            write(6,*) '==============================================='
         end if

      end subroutine spod_stream_finalize

!-----------------------------------------------------------------------
      subroutine spod_stream_save_modes(ifreq, nfreq, freq,
     $     evecs, evals, nblk, nsave)
!     Save SPOD modes at selected frequencies

         use krylov_subspace
         use spod_streaming_state
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: ifreq, nfreq, nblk, nsave
         real, intent(in) :: freq(nfreq), evals(nblk)
         complex(kind=kind(0.0d0)), intent(in) :: evecs(nblk, nblk)

         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm
         integer :: m, iblk, idx, nmodes
         character(len=3) :: prefix_re, prefix_im
         logical :: save_this_freq

!        Decide which frequencies to save
         save_this_freq = .false.
         if (ifreq == 1) save_this_freq = .true.
         if (ifreq == nfreq) save_this_freq = .true.
         if (mod(ifreq-1, max(1, nfreq/8)) == 0) save_this_freq = .true.

         if (.not. save_this_freq) return

         prefix_re = 'sRe'
         prefix_im = 'sIm'
         nmodes = min(nsave, nblk)

!        Reconstruct and save modes
         do m = 1, nmodes
            call k_zero(mode_re)
            call k_zero(mode_im)

!           Phi_m = Sum_i w_m(i) * Qhat(i) / sqrt(lambda_m * nblk)
            do iblk = 1, nblk
               coef_re = real(evecs(iblk, m))
               coef_im = aimag(evecs(iblk, m))

               idx = (ifreq - 1) * spod_s_nblk + iblk

!              Complex multiplication: mode += coef * blk
               call k_add2s2(mode_re, spod_s_dft_re(idx), coef_re)
               call k_add2s2(mode_re, spod_s_dft_im(idx), -coef_im)
               call k_add2s2(mode_im, spod_s_dft_re(idx), coef_im)
               call k_add2s2(mode_im, spod_s_dft_im(idx), coef_re)
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
            call outpost2(vx, vy, vz, pr, t, 0, prefix_re)

!           Output imaginary part
            call nopcopy(vx, vy, vz, pr, t,
     $           mode_im%vx, mode_im%vy, mode_im%vz,
     $           mode_im%pr, mode_im%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix_im)

!           Compute norm (ALL ranks must call - uses MPI_Allreduce)
            if (m == 1) then
               call k_norm(mode_norm, mode_re)
               if (nid == 0) then
                  write(6,'(A,F10.4,A,E12.4,A,E12.4)')
     $                 '      St =', freq(ifreq), ': lambda_1 =',
     $                 evals(1), ', ||Phi|| =', mode_norm
               end if
            end if
         end do

      end subroutine spod_stream_save_modes
