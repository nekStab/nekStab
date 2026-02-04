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

!  ─────────────────────────────────────────────────────────────────
!  Print header
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) then
            write(6,*) ''
            write(6,*) '═══════════════════════════════════════════════'
            write(6,*) '        MODAL ANALYSIS (Mode 6)'
            write(6,*) '═══════════════════════════════════════════════'
            write(6,*) ''
            write(6,'(A,A)')    '  Prefix:     ', trim(modal_prefix)
            write(6,'(A,I6)')   '  Snapshots:  ', modal_nsnap
            write(6,'(A,E12.4)')'  dt:         ', modal_dt
            write(6,'(A,I6)')   '  Modes save: ', modal_nsave
            write(6,'(A,L1)')   '  POD:        ', ifpod
            write(6,'(A,L1)')   '  DMD:        ', ifdmd
            write(6,'(A,L1)')   '  SPOD:       ', ifspod
            write(6,*) ''
         end if

!  ─────────────────────────────────────────────────────────────────
!  Validate parameters
!  ─────────────────────────────────────────────────────────────────
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

!  ─────────────────────────────────────────────────────────────────
!  PHASE 1: Load snapshots
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) 'Loading snapshots...'
         allocate(snaps(modal_nsnap))

         call load_files(snaps, modal_nsnap, modal_nsnap, modal_prefix)

         if (nid == 0) write(6,*) '  Loaded', modal_nsnap, 'snapshots'

!  ─────────────────────────────────────────────────────────────────
!  PHASE 2: Compute and subtract mean
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) 'Computing temporal mean...'
         call modal_compute_mean(snaps, modal_nsnap, mean_snap)
         call modal_subtract_mean(snaps, modal_nsnap, mean_snap)

!        Save mean field
         call nopcopy(vx, vy, vz, pr, t,
     $        mean_snap%vx, mean_snap%vy, mean_snap%vz,
     $        mean_snap%pr, mean_snap%t)
         call outpost2(vx, vy, vz, pr, t, 0, 'mea')
         if (nid == 0) write(6,*) '  Saved mean field as mea*'

!  ─────────────────────────────────────────────────────────────────
!  PHASE 3: POD
!  ─────────────────────────────────────────────────────────────────
         if (ifpod) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '───────────────────────────────────────────'
               write(6,*) '  POD (Proper Orthogonal Decomposition)'
               write(6,*) '───────────────────────────────────────────'
            end if
            call pod_compute(snaps, modal_nsnap, modal_nsave)
         end if

!  ─────────────────────────────────────────────────────────────────
!  PHASE 4: DMD
!  ─────────────────────────────────────────────────────────────────
         if (ifdmd) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '───────────────────────────────────────────'
               write(6,*) '  DMD (Dynamic Mode Decomposition)'
               write(6,*) '───────────────────────────────────────────'
            end if
            call dmd_compute(snaps, modal_nsnap, modal_dt,
     $                       dmd_rank, modal_nsave)
         end if

!  ─────────────────────────────────────────────────────────────────
!  PHASE 5: SPOD
!  ─────────────────────────────────────────────────────────────────
         if (ifspod) then
            if (nid == 0) then
               write(6,*) ''
               write(6,*) '───────────────────────────────────────────'
               write(6,*) '  SPOD (Spectral POD)'
               write(6,*) '───────────────────────────────────────────'
            end if
            call spod_compute(snaps, modal_nsnap, modal_dt,
     $                        spod_nfft, spod_noverlap, modal_nsave)
         end if

!  ─────────────────────────────────────────────────────────────────
!  Cleanup
!  ─────────────────────────────────────────────────────────────────
         deallocate(snaps)

         if (nid == 0) then
            write(6,*) ''
            write(6,*) '═══════════════════════════════════════════════'
            write(6,*) '  Modal analysis complete'
            write(6,*) '═══════════════════════════════════════════════'
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
!     Solves eigenvalue problem C v = λ v
!     Reconstructs spatial modes Φ_k = Σ_i v_k(i) snaps(i) / √(λ_k n)

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, nsave
         type(krylov_vector), intent(in) :: snaps(nsnap)

         real, allocatable :: C(:,:), eigvals(:), eigvecs(:,:)
         integer :: i, j

!  ─────────────────────────────────────────────────────────────────
!  Step 1: Form correlation matrix using k_dot
!  ─────────────────────────────────────────────────────────────────
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

!  ─────────────────────────────────────────────────────────────────
!  Step 2: Solve eigenvalue problem
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) '  Solving eigenvalue problem...'

         call eig_symmetric(C, eigvals, eigvecs, nsnap)

!  ─────────────────────────────────────────────────────────────────
!  Step 3: Reconstruct and save modes
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) '  Reconstructing POD modes...'

         call pod_reconstruct_modes(snaps, eigvecs, eigvals,
     $                              nsnap, nsave)

!  ─────────────────────────────────────────────────────────────────
!  Step 4: Write eigenvalue spectrum
!  ─────────────────────────────────────────────────────────────────
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
!     Φ_m = Σ_i v_m(i) * snaps(i) / sqrt(λ_m * nsnap)

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

!           Φ_m = Σ_i v_m(i) * snaps(i)
            do i = 1, nsnap
               call k_add2s2(mode, snaps(i), eigvecs(i, m))
            end do

!           Normalize: divide by sqrt(λ_m * nsnap)
            if (eigvals(m) > 0.0d0) then
               scale = 1.0d0 / sqrt(eigvals(m) * dble(nsnap))
               call k_cmult(mode, scale)
            end if

!           Output mode
            call nopcopy(vx, vy, vz, pr, t,
     $           mode%vx, mode%vy, mode%vz, mode%pr, mode%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix)

            if (nid == 0) then
               call k_norm(mode_norm, mode)
               write(6,'(A,I4,A,E12.4,A,F8.4)')
     $              '    Mode', m, ': λ =', eigvals(m),
     $              ', ||Φ|| =', mode_norm
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
      subroutine dmd_compute(snaps, nsnap, delta_t, rank, nsave)
!     Projected DMD (Schmid, 2010; Tu et al., 2014)
!
!     1. Form data matrices: X = [x1,...,x_{n-1}], Y = [x2,...,x_n]
!     2. SVD of X via eigendecomposition of X^T X
!     3. Project dynamics: Ã = U^H Y V Σ^{-1}
!     4. Eigendecomposition of Ã gives DMD eigenvalues
!     5. DMD modes: Φ = U W (or Y V Σ^{-1} W for exact DMD)

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: nsnap, rank, nsave
         real, intent(in) :: delta_t
         type(krylov_vector), intent(in) :: snaps(nsnap)

         integer :: n, r, i, j, m
         real, allocatable :: G(:,:), S(:), Vt(:,:), Sinv(:)
         real, allocatable :: Atilde(:,:), Atilde_work(:,:)
         real, allocatable :: Y_proj(:,:)
         complex(kind=kind(0.0d0)), allocatable :: dmd_evals(:)
         complex(kind=kind(0.0d0)), allocatable :: dmd_evecs(:,:)
         real :: alpha, total_energy, cumsum, tol
         real, parameter :: PI_VAL = 3.14159265358979323846d0

!        n-1 snapshot pairs for DMD
         n = nsnap - 1

         if (nid == 0) write(6,*) '  Using', n, 'snapshot pairs'

!  ─────────────────────────────────────────────────────────────────
!  Step 1: Form Gram matrix G = X^T X
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) '  Computing Gram matrix...'

         allocate(G(n, n))

         do j = 1, n
            do i = 1, j
               call k_dot(G(i,j), snaps(i), snaps(j))
               G(j,i) = G(i,j)
            end do
         end do

!  ─────────────────────────────────────────────────────────────────
!  Step 2: SVD via eigendecomposition of G
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) '  Computing SVD via eigendecomp...'

         allocate(S(n), Vt(n, n))
         call eig_symmetric(G, S, Vt, n)

!        S now contains σ², need sqrt for singular values
!        Also determine rank (truncation)
         total_energy = sum(S)
         if (rank > 0) then
            r = min(rank, n)
         else
!           Auto-rank: capture 99% energy
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

!        Convert eigenvalues to singular values and compute inverse
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

!  ─────────────────────────────────────────────────────────────────
!  Step 3: Project Y onto POD basis and form Ã
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) '  Forming projected operator...'

!        Ã(i,j) = Σ_k Σ_l <U_i, Y_k> V(k,j) Σ^{-1}_j
!        where U_i = Σ_m V(m,i) X_m / σ_i
!        So Ã(i,j) = Σ_k <snaps(k+1), snaps(m)> V(m,i) V(k,j) Σ^{-1}_i Σ^{-1}_j

         allocate(Y_proj(r, r))
         allocate(Atilde(r, r))

!        Y_proj(i,j) = <U_i, Y V Σ^{-1} e_j> = Σ_k <U_i, snaps(k+1)> V(k,j) Σ^{-1}_j
!        But U_i = Σ_m snaps(m) V(m,i) Σ^{-1}_i
!        So Y_proj(i,j) = Σ_k Σ_m <snaps(m), snaps(k+1)> V(m,i) V(k,j) Σ^{-1}_i Σ^{-1}_j

         Y_proj = 0.0d0
         do i = 1, r
            do j = 1, r
               do m = 1, n
!                 <snaps(m), snaps(m+1)> contributes to shifted correlation
                  call k_dot(alpha, snaps(m), snaps(m+1))
                  Y_proj(i,j) = Y_proj(i,j) + alpha *
     $                 Vt(m,i) * Vt(m,j) * Sinv(i) * Sinv(j)
               end do
            end do
         end do

!        Ã = Σ^{-1} U^H Y V Σ^{-1} = Y_proj (already scaled)
!        Actually need to rescale: Atilde(i,j) = σ_i * Y_proj(i,j)
         do i = 1, r
            do j = 1, r
               Atilde(i,j) = S(i) * Y_proj(i,j)
            end do
         end do

!  ─────────────────────────────────────────────────────────────────
!  Step 4: Eigendecomposition of Ã
!  ─────────────────────────────────────────────────────────────────
         if (nid == 0) write(6,*) '  Computing DMD eigenvalues...'

         allocate(dmd_evals(r), dmd_evecs(r, r))
         allocate(Atilde_work(r, r))
         Atilde_work = Atilde

         call eig(Atilde_work, dmd_evecs, dmd_evals, r)

!  ─────────────────────────────────────────────────────────────────
!  Step 5: Write spectrum and reconstruct modes
!  ─────────────────────────────────────────────────────────────────
         call dmd_write_spectrum(dmd_evals, r, delta_t)
         call dmd_reconstruct_modes(snaps, Vt, S, Sinv,
     $        dmd_evecs, dmd_evals, n, r, nsave)

         deallocate(G, S, Vt, Sinv, Y_proj, Atilde, Atilde_work)
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
!     Reconstruct DMD modes: Φ = X V Σ^{-1} W (projected DMD modes)
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

         prefix_re = 'dRe'
         prefix_im = 'dIm'
         nmodes = min(nsave, r)

!        DMD mode m: Φ_m = X V Σ^{-1} w_m
!        where w_m is eigenvector of Ã
         do m = 1, nmodes
            call k_zero(mode_re)
            call k_zero(mode_im)

!           Φ_m = Σ_i snaps(i) * (Σ_j V(i,j) Σ^{-1}_j w_m(j))
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

            if (nid == 0) then
               mu_mag = abs(evals(m))
               freq = atan2(aimag(evals(m)), real(evals(m))) /
     $              (2.0d0 * PI_VAL * modal_dt)
               call k_norm(mode_norm, mode_re)
               write(6,'(A,I4,A,F8.4,A,F10.4,A,E12.4)')
     $              '    Mode', m, ': |μ| =', mu_mag,
     $              ', St =', freq, ', ||Φ_re|| =', mode_norm
            end if
         end do

      end subroutine dmd_reconstruct_modes

!-----------------------------------------------------------------------
      subroutine spod_compute(snaps, nsnap, delta_t, nfft, noverlap, nsave)
!     Spectral POD (Towne et al., 2018)
!
!     1. Divide snapshots into overlapping blocks
!     2. Apply window and FFT each block
!     3. For each frequency, form cross-spectral density matrix
!     4. Eigendecomposition gives SPOD modes at that frequency

         use krylov_subspace
         use fourier_fftw
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

!  ─────────────────────────────────────────────────────────────────
!  Loop over frequencies
!  ─────────────────────────────────────────────────────────────────
         do ifreq = 1, nfreq

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
!     For each spatial point:
!       1. Apply window: x_windowed(j) = window(j) * snaps(j)%field
!       2. FFT the time series
!       3. Extract real/imag at frequency ifreq
!       4. Apply normalization

         use krylov_subspace
         use fourier_fftw
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
!        Factor of 2 for discarding negative frequencies
         norm_factor = 2.0d0 / (dble(nfft) * win_weight)
         if (ifreq == 1 .or. ifreq == nfreq) then
            norm_factor = norm_factor / 2.0d0  ! DC and Nyquist
         end if

         call k_zero(out_re)
         call k_zero(out_im)

         nv = nx1*ny1*nz1*nelv
         nt = nx1*ny1*nz1*nelt

!        Process vx component
         do ipt = 1, nv
            do j = 1, nfft
               time_series(j) = window(j) *
     $              snaps(blk_start + j - 1)%vx(ipt)
            end do
            call fft_r2c(nfft, time_series, spectrum)
            out_re%vx(ipt) = real(spectrum(ifreq)) * norm_factor
            out_im%vx(ipt) = aimag(spectrum(ifreq)) * norm_factor
         end do

!        Process vy component
         do ipt = 1, nv
            do j = 1, nfft
               time_series(j) = window(j) *
     $              snaps(blk_start + j - 1)%vy(ipt)
            end do
            call fft_r2c(nfft, time_series, spectrum)
            out_re%vy(ipt) = real(spectrum(ifreq)) * norm_factor
            out_im%vy(ipt) = aimag(spectrum(ifreq)) * norm_factor
         end do

!        Process vz component (if 3D)
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

!           Φ_m = Σ_i w_m(i) * Q̂_k(i) / sqrt(λ_m * nblk)
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
     $              '    St =', freq(ifreq), ': λ₁ =', evals(1),
     $              ', ||Φ|| =', mode_norm
            end if
         end do

      end subroutine spod_save_modes

!-----------------------------------------------------------------------
      subroutine hamming_window(n, window, win_weight)
!     Compute Hamming window and its energy normalization factor

         implicit none
         integer, intent(in) :: n
         real, intent(out) :: window(n), win_weight
         real, parameter :: PI_VAL = 3.14159265358979323846d0
         integer :: i

         do i = 1, n
            window(i) = 0.54d0 - 0.46d0 * cos(2.0d0 * PI_VAL *
     $           dble(i-1) / dble(n-1))
         end do

!        Window energy for normalization
         win_weight = sum(window**2) / dble(n)

      end subroutine hamming_window
