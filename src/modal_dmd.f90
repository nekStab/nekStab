!-----------------------------------------------------------------------
!     modal_dmd.f90: Projected DMD (Schmid, 2010; Tu et al., 2014)
!
!     Given snapshots x_1..x_n, finds best-fit linear operator A
!     such that x_{k+1} ~ A x_k. Eigenvalues mu encode:
!       growth rate sigma = log|mu| / dt
!       frequency   St    = arg(mu) / (2*pi*dt)
!-----------------------------------------------------------------------
      module modal_dmd

         use krylov_subspace

         implicit none
         private

         public :: dmd_compute

      contains

!-----------------------------------------------------------------------
      subroutine dmd_compute(snaps, nsnap, delta_t, rank, nsave)
!     Projected DMD algorithm:
!     1. Form Gram matrix G = X^T X (symmetric)
!     2. SVD via eigendecomp of G: G = V Lambda V^T
!     3. Shifted Gram: G_shift(k,m) = <x_{k+1}, x_m> (NOT symmetric)
!     4. Project operator: Atilde = Sinv V^T G_shift V Sinv
!     5. Eigendecomp of Atilde -> DMD eigenvalues and modes

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
         real, allocatable :: dmd_norms(:)
         complex(kind=kind(0.0d0)), allocatable :: dmd_evals(:)
         complex(kind=kind(0.0d0)), allocatable :: dmd_evecs(:,:)
         real :: total_energy, cumsum, tol
         real, parameter :: PI_VAL = 3.14159265358979323846d0

         external :: eig_symmetric, eig

!        Number of snapshot pairs: X = [x_1..x_{n-1}], Y = [x_2..x_n]
         n = nsnap - 1

         if (nid == 0) write(6,*) '  Using', n, 'snapshot pairs'

!        Form symmetric Gram matrix G(i,j) = <x_i, x_j>
         if (nid == 0) write(6,*) '  Computing Gram matrix...'

         allocate(G(n, n))

         do j = 1, n
            do i = 1, j
               call k_dot(G(i,j), snaps(i), snaps(j))
               G(j,i) = G(i,j)
            end do
         end do

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

         if (nid == 0) write(6,'(A,I4,A,F6.2,A)')
     $        '  Using rank r =', r, ' (',
     $        100.0d0 * sum(S(1:r)) / total_energy, '% energy)'

!        Write SVD diagnostics (before sqrt overwrites S)
         if (nid == 0) then
            open(unit=77, file='dmd_svd.dat', status='replace')
            write(77, '(A)') '# DMD SVD Singular Values'
            write(77, '(A,I6,A,F6.2,A)')
     $         '# rank_used = ', r, '  (',
     $         100.0d0 * sum(S(1:r)) / total_energy, '% energy)'
            write(77, '(A,I6)')
     $         '# total_snapshots = ', n
            write(77, '(A)')
     $    '# i       sigma_i         energy_%    cumulative_%'

            cumsum = 0.0d0
            do i = 1, n
               cumsum = cumsum + S(i)
               write(77, '(I6, 3E16.8)') i, sqrt(S(i)),
     $              100.0d0 * S(i) / total_energy,
     $              100.0d0 * cumsum / total_energy
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
!           k=n requires snaps(nsnap), not available in G
            call k_dot(G_shift(n,m), snaps(nsnap), snaps(m))
         end do

!        Form projected operator Atilde = Sinv V^T G_shift V Sinv
         if (nid == 0) write(6,*) '  Forming projected operator...'

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

!        Eigendecomposition of Atilde (non-symmetric -> complex eigenvalues)
         if (nid == 0) write(6,*) '  Computing DMD eigenvalues...'

         allocate(dmd_evals(r), dmd_evecs(r, r))
         allocate(Atilde_work(r, r))
         Atilde_work = Atilde

         call eig(Atilde_work, dmd_evecs, dmd_evals, r)

!        Reconstruct modes and write spectrum
         allocate(dmd_norms(min(nsave, r)))
         call dmd_reconstruct_modes(snaps, Vt, S, Sinv,
     $        dmd_evecs, dmd_evals, n, r, nsave, dmd_norms)
         call dmd_write_spectrum(dmd_evals, r, delta_t,
     $        dmd_norms, min(nsave, r))

         deallocate(G, G_shift, S, Vt, Sinv, Atilde, Atilde_work)
         deallocate(dmd_evals, dmd_evecs, dmd_norms)

         if (nid == 0) write(6,*) '  DMD complete'

      end subroutine dmd_compute

!-----------------------------------------------------------------------
      subroutine dmd_write_spectrum(evals, r, delta_t,
     $                              norms, nnorms)
!     Write DMD eigenvalue spectrum

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: r, nnorms
         complex(kind=kind(0.0d0)), intent(in) :: evals(r)
         real, intent(in) :: delta_t, norms(nnorms)

         real :: sigma, omega, mu_mag, freq
         real, parameter :: PI_VAL = 3.14159265358979323846d0
         integer :: m

         if (nid /= 0) return

         open(unit=78, file='dmd_spectrum.dat', status='replace')
         write(78, '(A)') '# DMD Eigenvalue Spectrum'
         write(78, '(A)') '# mode   |mu|      sigma       omega       '
     $      // 'St          mu_real     mu_imag     ||Phi_re||'

         do m = 1, r
            mu_mag = abs(evals(m))
            sigma = log(mu_mag) / delta_t
            omega = atan2(aimag(evals(m)), real(evals(m))) / delta_t
            freq = omega / (2.0d0 * PI_VAL)

            if (m <= nnorms) then
               write(78, '(I6, 7E14.6)') m, mu_mag, sigma, omega,
     $              freq, real(evals(m)), aimag(evals(m)), norms(m)
            else
               write(78, '(I6, 6E14.6)') m, mu_mag, sigma, omega,
     $              freq, real(evals(m)), aimag(evals(m))
            end if
         end do

         close(78)
         write(6,*) '  Wrote dmd_spectrum.dat'

      end subroutine dmd_write_spectrum

!-----------------------------------------------------------------------
      subroutine dmd_reconstruct_modes(snaps, V, S, Sinv,
     $     evecs, evals, n, r, nsave, norms)
!     DMD modes: Phi = X V Sinv W (projected DMD modes)

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: n, r, nsave
         type(krylov_vector), intent(in) :: snaps(n)
         real, intent(in) :: V(n, n), S(r), Sinv(r)
         complex(kind=kind(0.0d0)), intent(in) :: evecs(r, r), evals(r)
         real, intent(out) :: norms(min(nsave, r))

         type(krylov_vector) :: mode_re, mode_im
         real :: coef_re, coef_im, mode_norm, mu_mag, freq
         integer :: m, i, j, nmodes
         real, parameter :: PI_VAL = 3.14159265358979323846d0
         character(len=3) :: prefix_re, prefix_im

         prefix_re = 'dm1'
         prefix_im = 'dm2'
         nmodes = min(nsave, r)

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

            call k_norm(mode_norm, mode_re)
            norms(m) = mode_norm

!           Skip output if mode is invalid (NaN or very large)
            if (mode_norm /= mode_norm .or. mode_norm > 1.0d30) then
               if (nid == 0) write(6,*) '  WARNING: DMD mode', m,
     $              'has invalid values, skipping'
               norms(m) = 0.0d0
               cycle
            end if

            call nopcopy(vx, vy, vz, pr, t,
     $           mode_re%vx, mode_re%vy, mode_re%vz,
     $           mode_re%pr, mode_re%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix_re)

            call nopcopy(vx, vy, vz, pr, t,
     $           mode_im%vx, mode_im%vy, mode_im%vz,
     $           mode_im%pr, mode_im%t)
            call outpost2(vx, vy, vz, pr, t, 0, prefix_im)

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

      end module modal_dmd
