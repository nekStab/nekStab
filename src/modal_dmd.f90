      !-----------------------------------------------------------------------
      ! modal_dmd.f90 — Projected DMD (Schmid, 2010; Tu et al., 2014)
      !
      ! Purpose:
      !   Given snapshots x_1..x_n, finds best-fit linear operator A
      !   such that x_{k+1} ~ A x_k. Eigenvalues mu encode:
      !     growth rate sigma = log|mu| / dt
      !     frequency   St    = arg(mu) / (2*pi*dt)
      !
      ! Public interface:
      !   dmd_compute — Full DMD: Gram matrix, SVD, projected operator
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL
      !-----------------------------------------------------------------------
module modal_dmd

    use krylov_subspace
    use krylov_inner_products
    use nekstab_lapack
    use nekstab_vectors
    use nekstab_nek_bridge

    implicit none
         private

         public :: dmd_compute

      contains

      !-----------------------------------------------------------------------
      ! dmd_compute — Projected DMD algorithm
      !
      !   1. Form Gram matrix G = X^T X (symmetric)
      !   2. SVD via eigendecomp of G: G = V Lambda V^T
      !   3. Shifted Gram: G_shift(k,m) = <x_{k+1}, x_m> (NOT symmetric)
      !   4. Project operator: Atilde = Sinv V^T G_shift V Sinv
      !   5. Eigendecomp of Atilde -> DMD eigenvalues and modes
      !-----------------------------------------------------------------------
       subroutine dmd_compute(snaps, nsnap, delta_t, rank, nsave)

          integer, intent(in) :: nsnap, rank, nsave
         real, intent(in) :: delta_t
         type(krylov_vector), intent(in) :: snaps(nsnap)

         integer :: n, r, i, j, k, m
         real, allocatable :: G(:,:), S(:), Vt(:,:), Sinv(:)
         real, allocatable :: G_shift(:,:)
         real, allocatable :: Atilde(:,:), Atilde_work(:,:)
         real, allocatable :: Atilde_tmp(:,:), Atilde_proj(:,:)
         real, allocatable :: dmd_norms(:)
         complex(nekStab_dp), allocatable :: dmd_evals(:)
         complex(nekStab_dp), allocatable :: dmd_evecs(:,:)
         real :: total_energy, cumsum, tol

!        Number of snapshot pairs: X = [x_1..x_{n-1}], Y = [x_2..x_n]
         n = nsnap - 1

         if (nid == 0) write(6,*) '  Using', n, 'snapshot pairs'

!        Form symmetric Gram matrix G(i,j) = <x_i, x_j>
         if (nid == 0) write(6,*) '  Computing Gram matrix...'

         allocate(G(n, n))

         call k_gram_matrix(G, snaps, snaps, n, n, n)

!        SVD of X via eigendecomposition of G
         if (nid == 0) write(6,*) '  Computing SVD via eigendecomp...'

         allocate(S(n), Vt(n, n))
         call eig_symmetric(G, S, Vt, n)

!        Determine truncation rank
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

         if (nid == 0) write(6,'(A,I4,A,F6.2,A)') &
            '  Using rank r =', r, ' (', &
            100.0d0 * sum(S(1:r)) / total_energy, '% energy)'

!        Write SVD diagnostics (before sqrt overwrites S)
         if (nid == 0) then
            open(unit=77, file='dmd_svd.dat', status='replace')
            write(77, '(A)') '# DMD SVD Singular Values'
            write(77, '(A,I6,A,F6.2,A)') &
               '# rank_used = ', r, '  (', &
               100.0d0 * sum(S(1:r)) / total_energy, '% energy)'
            write(77, '(A,I6)') &
               '# total_snapshots = ', n
            write(77, '(A)') &
               '# i       sigma_i         energy_%    cumulative_%'

            cumsum = 0.0d0
            do i = 1, n
               cumsum = cumsum + S(i)
               write(77, '(I6, 3E16.8)') i, sqrt(S(i)), &
                  100.0d0 * S(i) / total_energy, &
                  100.0d0 * cumsum / total_energy
            end do

            close(77)
            write(6,*) '  Wrote dmd_svd.dat'
         end if

!        Convert eigenvalues (sigma^2) to singular values (sigma) and inverses
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

!        Form shifted Gram matrix G_shift(k,m) = <x_{k+1}, x_m>
!        NOT symmetric -- the asymmetry produces complex eigenvalues.
!        Reuse G where possible: G_shift(k,m) = G(k+1,m) for k < n
         allocate(G_shift(n, n))

         do m = 1, n
            do k = 1, n-1
               G_shift(k,m) = G(k+1, m)
            end do
         end do
!        k=n: single batch projection for last row
         block
            real, allocatable :: last_row(:)
            allocate(last_row(n))
            call k_project(last_row, snaps(1), snaps(nsnap), n)
            do m = 1, n
               G_shift(n, m) = last_row(m)
            end do
            deallocate(last_row)
         end block

!        Form projected operator via BLAS (replaces O(r^2 n^2) scalar loops).
!
!        NOTE: despite its name, Vt stores eigenvector COLUMNS (i.e. V,
!        not V^T).  eig_symmetric returns columns = eigenvectors.
!
!        Derivation: tracing the original loop indices gives
!          Atilde(i,j) = Sinv(i) * sum_{k,m} G_shift(k,m) Vt(m,i) Vt(k,j) * Sinv(j)
!                      = Sinv(i) * [V^T G_shift V](j,i) * Sinv(j)
!        i.e.  Atilde = diag(Sinv) * (V^T G_shift V)^T * diag(Sinv)
!        Note the transpose: G_shift is NOT symmetric, so the (j,i) index
!        order matters.  The two dgemm calls compute M = V^T G_shift V,
!        then the double loop reads M(j,i) (transposed) with Sinv scaling.
!
!        Cost: 2 dgemm of size (n x r) instead of 4-nested scalar loop.
         if (nid == 0) write(6,*) '  Forming projected operator...'

         allocate(Atilde(r, r))
         allocate(Atilde_tmp(n, r), Atilde_proj(r, r))

!        Step 1: Atilde_tmp(n,r) = G_shift(n,n) * V_r(n,r)
         call dgemm('N', 'N', n, r, n, 1.0d0, &
            G_shift, n, Vt, n, 0.0d0, Atilde_tmp, n)
!        Step 2: Atilde_proj(r,r) = V_r^T(r,n) * Atilde_tmp(n,r)
         call dgemm('T', 'N', r, r, n, 1.0d0, &
            Vt, n, Atilde_tmp, n, 0.0d0, Atilde_proj, r)
!        Step 3: transpose + diagonal scaling
         do j = 1, r
            do i = 1, r
               Atilde(i,j) = Sinv(i) * Atilde_proj(j,i) * Sinv(j)
            end do
         end do

!        Eigendecomposition of Atilde (non-symmetric -> complex eigenvalues)
         if (nid == 0) write(6,*) '  Computing DMD eigenvalues...'

         allocate(dmd_evals(r), dmd_evecs(r, r))
         allocate(Atilde_work(r, r))
         Atilde_work = Atilde

         call eig(Atilde_work, dmd_evecs, dmd_evals, r)

!        Reconstruct modes and write spectrum
         allocate(dmd_norms(min(nsave, r)))
          call dmd_reconstruct_modes(snaps, Vt, S, Sinv, &
     &   dmd_evecs, dmd_evals, n, r, nsave, dmd_norms)
          call dmd_write_spectrum(dmd_evals, r, delta_t, &
     &   dmd_norms, min(nsave, r))

         deallocate(G, G_shift, S, Vt, Sinv, Atilde, Atilde_work)
         deallocate(Atilde_tmp, Atilde_proj)
         deallocate(dmd_evals, dmd_evecs, dmd_norms)

         if (nid == 0) write(6,*) '  DMD complete'

      end subroutine dmd_compute

      !-----------------------------------------------------------------------
      ! dmd_write_spectrum — Write DMD eigenvalue spectrum to file
      !-----------------------------------------------------------------------
       subroutine dmd_write_spectrum(evals, r, delta_t, norms, nnorms)


          integer, intent(in) :: r, nnorms
         complex(nekStab_dp), intent(in) :: evals(r)
         real, intent(in) :: delta_t, norms(nnorms)

         real :: sigma, omega, mu_mag, freq
         integer :: m

         if (nid /= 0) return

         open(unit=78, file='dmd_spectrum.dat', status='replace')
         write(78, '(A)') '# DMD Eigenvalue Spectrum'
         write(78, '(A)') '# mode   |mu|      sigma       omega       ' &
            // 'St          mu_real     mu_imag     ||Phi_re||'

         do m = 1, r
            mu_mag = abs(evals(m))
            sigma = log(mu_mag) / delta_t
            omega = atan2(aimag(evals(m)), real(evals(m))) / delta_t
            freq = omega / (2.0d0 * NEKSTAB_PI)

             if (m <= nnorms) then
                write(78, '(I6, 7E14.6)') m, mu_mag, sigma, omega, &
     &              freq, real(evals(m)), aimag(evals(m)), norms(m)
            else
                write(78, '(I6, 6E14.6)') m, mu_mag, sigma, omega, &
     &              freq, real(evals(m)), aimag(evals(m))
            end if
         end do

         close(78)
         write(6,*) '  Wrote dmd_spectrum.dat'

      end subroutine dmd_write_spectrum

      !-----------------------------------------------------------------------
      ! dmd_reconstruct_modes — Reconstruct DMD spatial modes
      !
      !   DMD modes: Phi = X V Sinv W (projected DMD modes)
      !-----------------------------------------------------------------------
       subroutine dmd_reconstruct_modes(snaps, V, S, Sinv, evecs, evals, n, r, nsave, norms)


          integer, intent(in) :: n, r, nsave
         type(krylov_vector), intent(in) :: snaps(n)
         real, intent(in) :: V(n, n), S(r), Sinv(r)
         complex(nekStab_dp), intent(in) :: evecs(r, r), evals(r)
         real, intent(out) :: norms(min(nsave, r))

         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm, mu_mag, freq
         integer :: m, i, j, nmodes
         character(len=3) :: prefix_re, prefix_im

         prefix_re = 'dm1'
         prefix_im = 'dm2'
         nmodes = min(nsave, r)

         if (nid == 0) write(6,'(A,I4,A)') '   Reconstructing ', nmodes, ' DMD modes...'


         do m = 1, nmodes
            call k_zero(mode_re)
            call k_zero(mode_im)

!           Phi_m = Sum_i snaps(i) * (Sum_j V(i,j) Sinv(j) w_m(j))
            do i = 1, n
               coef_re = 0.0d0
               coef_im = 0.0d0
               do j = 1, r
                  coef_re = coef_re + V(i,j) * Sinv(j) * real(evecs(j,m))

                  coef_im = coef_im + V(i,j) * Sinv(j) * aimag(evecs(j,m))

               end do
               call k_add2s2(mode_re, snaps(i), coef_re)
               call k_add2s2(mode_im, snaps(i), coef_im)
            end do

            call k_norm(mode_norm, mode_re)
            norms(m) = mode_norm

!           Skip output if mode is invalid (NaN or very large)
            if (mode_norm /= mode_norm .or. mode_norm > 1.0d30) then
               if (nid == 0) write(6,*) '  WARNING: DMD mode', m, 'has invalid values, skipping'

               norms(m) = 0.0d0
               cycle
            end if

            call nopcopy(vx, vy, vz, pr, t, mode_re%vx, mode_re%vy, mode_re%vz, mode_re%pr, mode_re%t)

            call outpost2(vx, vy, vz, pr, t, 0, prefix_re)

            call nopcopy(vx, vy, vz, pr, t, mode_im%vx, mode_im%vy, mode_im%vz, mode_im%pr, mode_im%t)

            call outpost2(vx, vy, vz, pr, t, 0, prefix_im)

            if (nid == 0) then
               mu_mag = abs(evals(m))
               freq = atan2(aimag(evals(m)), real(evals(m))) / (2.0d0 * NEKSTAB_PI * modal_dt)

               write(6,'(A,I4,A,F8.4,A,F10.4,A,E12.4)') '    Mode', m, &
                  ': |mu| =', mu_mag, ', St =', freq, ', ||Phi_re|| =', mode_norm

            end if
         end do

      end subroutine dmd_reconstruct_modes

      end module modal_dmd
