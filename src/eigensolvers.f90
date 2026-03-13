      !-----------------------------------------------------------------------
      ! eigensolvers.f90 -- Eigensolver framework for nekStab
      !
      ! Purpose:
      !   Implements the Krylov-Schur eigensolver for computing leading
      !   eigenvalues/eigenmodes of the linearized Navier-Stokes operator,
      !   including Arnoldi factorization, Schur condensation with
      !   adaptive restart, eigenmode output, and checkpointing.
      !
      ! Public interface:
      !   inner_product         -- weighted L2 inner product
      !   norm                  -- vector norm from inner product
      !   krylov_schur          -- main Krylov-Schur eigensolver
      !   outpost_ks            -- eigenmode output and spectrum files
      !   schur_condensation    -- Krylov-Schur restart via Schur form
      !   select_eigenvalues    -- adaptive eigenvalue selection
      !   ensure_conjugate_pairs -- keep conjugate pairs together
      !
      ! Dependencies:
      !   krylov_subspace, nekstab_krylov_decomposition,
      !   nekstab_lapack, nekstab_argsort, nekstab_vectors,
      !   nekstab_matvec, nekstab_io, SIZE, TOTAL
      !-----------------------------------------------------------------------

   module nekstab_eigensolvers
         use krylov_subspace
         use nekstab_nek_bridge
   use nekstab_krylov_decomposition
   use nekstab_lapack
   use nekstab_argsort
   use nekstab_vectors
   use nekstab_matvec
   use nekstab_io
   use nekstab_diagnostics
   use nekstab_noise
   implicit none
   private
   public :: inner_product, norm, krylov_schur,&
      outpost_ks, schur_condensation,&
      select_eigenvalues, ensure_conjugate_pairs
   contains

      !-----------------------------------------------------------------------
      ! NOTE: inner_product and norm are now defined in krylov_subspace.f90.
      ! They are re-exported here via the public statement for backwards
      ! compatibility with callers that 'use nekstab_eigensolvers'.
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! krylov_schur -- Main Krylov-Schur eigensolver
      !
      ! Purpose:
      !   Builds a Krylov subspace via Arnoldi factorization, computes
      !   eigenvalues of the Hessenberg matrix, and applies Schur
      !   condensation restarts until convergence. Outputs converged
      !   eigenmodes to disk.
      !-----------------------------------------------------------------------
    subroutine krylov_schur
    use krylov_subspace

       !     -----Krylov basis V for the projection M*V = V*H -----
   type(krylov_vector), allocatable, dimension(:) :: Q

      !     ----- Upper Hessenberg matrix -----
   real, allocatable, dimension(:, :) :: H
   real, allocatable, dimension(:, :) :: b_vec

      !     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----
   complex(nekStab_dp), allocatable, dimension(:) :: vals
   complex(nekStab_dp), allocatable, dimension(:, :) :: vecs

   real, allocatable, dimension(:) :: residual

      !     ----- Miscellaneous -----
   type(krylov_vector) :: wrk, wrk2

   integer :: mstart, converged_eigenvalues, m, total_matvecs
   real :: alpha
   logical :: converged
   integer :: i, j
   character(len=30) :: filename

      !     ----- Allocate arrays -----
   allocate (Q(k_dim + 1))
   allocate (H(k_dim + 1, k_dim), b_vec(1, k_dim), vals(k_dim), vecs(k_dim, k_dim), residual(k_dim))

   time = 0.0d0
   H(:, :) = 0.0d0
   b_vec = 0.0d0
   residual = 0.0d0
   do i = 1, k_dim + 1
   call k_zero(Q(i))
   end do
      !call k_zero(Q(1:k_dim + 1))

      !     ----- Loading baseflow from disk (optional) -----

   if (ifldbf) then !skip loading if single run
   if (nid == 0) write (*, *) 'Loading base flow from disk:'
   write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)
   if (nid == 0) write (*, *) ' Number os scalars found (npscal): ', npscal
   if (nid == 0) write (*, *) ' ifldbf done.'
   else
   if (nid == 0) write (*, *) 'Baseflow prescribed by the useric function in the .usr'
   end if

      !     ----- Save baseflow to disk (recommended) -----
   call nopcopy(ubase, vbase, wbase, pbase, tbase, vx, vy, vz, pr, t)

      !     ----- Prepare stability parameters -----

   if (istep == 0 .and. (&
      uparam(1) == 3.11 .or. & ! Floquet direct
      uparam(1) == 3.21 .or. & ! Floquet adjoint
      uparam(1) == 3.31 & ! Floquet direct-adjoint
      )) then
   param(10) = time ! upo period in field
   if (nid == 0) write (6, *) 'Floquet mode !!!'
   if (nid == 0) write (6, *) ' getting endTime from file: endTime=', param(10)
   end if
   call bcast(param(10), wdsize)

      !     ----- First vector (new from noise or restart) -----

   if (uparam(2) == 0) then

   if (nid == 0) write (6, *) 'Starting first Arnoldi decomposition...'

      !     ----- Creates seed vector for the Krylov subspace -----

   if (ifseed_nois) then ! noise as initial seed

   if (nid == 0) write (6, *) 'Filling fields with noise...'
   call op_add_noise(wrk2%vx, wrk2%vy, wrk2%vz)
   if (ifto) call add_noise_scal(wrk2%t(:, 1), 9.0e4, 3.0e3, 4.0e5)
   if (ldimt > 1) then
   do m = 2, ldimt
   if (ifpsco(m - 1)) call add_noise_scal(wrk2%t(:, m), 9.0e1*m, 3.0e2*m, 4.0e1*m)
   end do
   end if
   call k_normalize(wrk2, alpha)
      !call outpost2(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%t, nof, 'NOS')
   call matvec(wrk, wrk2)
      !call outpost2(wrk%vx, wrk%vy, wrk%vz, wrk%pr, wrk%t, nof, 'NOS')

   elseif (ifseed_symm) then ! symmetry initial seed

   if (nid == 0) write (6, *) 'Enforcing symmetric seed perturb...'
   call add_symmetric_seed(wrk%vx, wrk%vy, wrk%vz, wrk%t(:, 1))

   elseif (ifseed_load) then ! loading initial seed (e.g. Re_ )

   if (uparam(01) >= 3.0 .and. uparam(01) < 3.2) then
   write (filename, '(a,a,a)') 'dRe', trim(SESSION), '0.f00001'

   elseif (uparam(01) >= 3.2 .and. uparam(01) < 3.3) then
   write (filename, '(a,a,a)') 'aRe', trim(SESSION), '0.f00001'
   end if

   if (nid == 0) write (*, *) 'Load real part of mode 1 as seed: ', filename
   call load_fld(filename)
   call nopcopy(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%t, vx, vy, vz, pr, t)
   call k_normalize(wrk2, alpha)
   call matvec(wrk, wrk2)

   else

   call nopcopy(wrk%vx, wrk%vy, wrk%vz, wrk%pr, wrk%t, ubase, vbase, wbase, pbase, tbase)
   call k_normalize(wrk, alpha)

   end if

      !     ----- Normalized to unit-norm -----
   mstart = 1; istep = 1; time = 0.0d0

   call k_copy(Q(1), wrk)

   if (ifres) then
   call whereyouwant('KRY', 1)
   call outpost2(Q(1)%vx, Q(1)%vy, Q(1)%vz, Q(1)%pr, Q(1)%t, nof, 'KRY')
   end if

   elseif (uparam(2) > 0) then

   mstart = int(uparam(2))

   if (nid == 0) then

   write (6, *) 'Restarting from:', mstart
   write (6, '(a,a,i4.4)') ' Loading Hessenberg matrix: HES', trim(SESSION), mstart
   write (filename, '(a,a,i4.4)') 'HES', trim(SESSION), mstart

   open (67, file=trim(filename), status='unknown', form='formatted')

   if (k_dim < mstart) then !subsampling
   do i = 1, k_dim + 1
   do j = 1, mstart
   if (j <= k_dim) read (67, "(1E15.7)") H(i, j)
   end do
   end do
   else

   read (67, *) ((H(i, j), j=1, mstart), i=1, mstart + 1)

   end if

   close (67)
   write (6, *) 'Broadcast H matrix to all procs...'

   end if !nid.eq.0
   call bcast(H, (k_dim + 1)*k_dim*wdsize) !broadcast H matrix to all procs

      !     if(nid.eq.2)then !-> debug only
      !     write(filename,'(a,a,i4.4)')'HESloaded',trim(SESSION),mstart
      !     write(6,*) ''
      !     open(67,file=trim(filename),status='unknown',form='formatted')
      !     write(67,*) beta
      !     do i = 1,mstart+1
      !     do j = 1,mstart
      !     write(67,*) H(i,j)
      !     enddo
      !     enddo
      !     close(67)
      !     endif

   mstart = mstart + 1 !careful here!
   call load_files(Q, mstart, k_dim + 1, 'KRY')
   if (nid == 0) write (6, *) 'Restart fields loaded to memory!'

   end if

      !     ======================================
      !     =====                            =====
      !     ===== Krylov-Schur decomposition =====
      !     =====                            =====
      !     ======================================

   schur_cnt = 0
   total_matvecs = 0
   converged = .false.

   do while (.not. converged)

      !     --> Arnoldi factorization.
   total_matvecs = total_matvecs + (k_dim - mstart + 1)
   call arnoldi_factorization(Q, H, mstart, k_dim, k_dim)

      !     if(nid.eq.0) then
      !        open(unit=12345, file="Hessenberg_matrix.dat")
      !        write(12345, *) H(1:k_dim, 1:k_dim)
      !        close(12345)
      !     endif

      !     --> Compute the eigenspectrum of the Hessenberg matrix.
   call eig(H(1:k_dim, 1:k_dim), vecs, vals, k_dim)

      !     --> Check the residual of the eigenvalues.
      !         The standard Arnoldi residual is: r_i = |beta * e_k^T * y_i|
      !         where beta = H(k+1,k), e_k is the k-th unit vector, and y_i is the
      !         i-th eigenvector of the Hessenberg matrix H.
      !
      !         IMPORTANT: Using only an absolute tolerance can be problematic:
      !         - For |lambda| << 1: tolerance may be too loose, accepting unconverged values
      !         - For |lambda| >> 1: tolerance may be too strict
      !         We use BOTH absolute and relative criteria: an eigenvalue is converged if
      !         residual < eigen_tol OR residual < eigen_tol * |lambda|
      !         This ensures small eigenvalues use absolute tolerance while large ones
      !         use relative tolerance, preventing spurious acceptance of unconverged modes.

   residual = abs(H(k_dim + 1, k_dim)*vecs(k_dim, :))

      !     --> Count converged eigenvalues using combined absolute/relative criterion.
      !         We use: converged if residual < eigen_tol * max(|lambda|, 1)
      !         This ensures:
      !         - For |lambda| < 1: we use absolute tolerance (don't over-tighten)
      !         - For |lambda| >= 1: we use relative tolerance (don't under-tighten)
      !         This is important for stability analysis where |lambda| ~ 1 matters most.
   converged_eigenvalues = 0
   do i = 1, k_dim
   if (residual(i) < eigen_tol * max(abs(vals(i)), 1.0d0)) then
   converged_eigenvalues = converged_eigenvalues + 1
   end if
   end do

   if (nid == 0) write (6, *) 'total eigenvalues converged:', converged_eigenvalues

      !     --> Select whether to stop or apply Schur condensation depending on schur_tgt.
   select case (schur_tgt)

      !     --> k-step Arnoldi factorization completed.
   case (:0)
   converged = .true.
   if (nid == 0) write (6, *) 'Arnoldi factorization completed.'

      !     --> Krylov-Schur factorization.
   case (1:)
   if (converged_eigenvalues >= schur_tgt) then ! Krylov-Schur factorization completed.
   converged = .true.
   else ! Apply Schur condensation before restarting the factorization.
   schur_cnt = schur_cnt + 1
   if (nid == 0) write (6, *) 'Starting Schur condensation phase.', schur_cnt
   call schur_condensation(mstart, H, Q, k_dim)
   end if

   end select

   end do

      !  if (nid == 0) open (unit=99, file="orthonormality.dat")
      !  do i = 1, k_dim
      !     call k_norm(alpha, Q(i))
      !     if (nid == 0) write (99, '("Norm of the ", I4, "th mode = ", F20.14)') i, alpha
      !     do j = i + 1, k_dim
      !        call k_dot(alpha, q(i), q(j))
      !        if (nid == 0) write (99, '("Orthogonality between mode ", I4, " and mode ", I4, " = ", E15.7)') i, j, alpha
      !     end do
      !     if (nid == 0) write (99, *)
      !  end do
      !  if (nid == 0) close (unit=99)

   if (nid == 0) then
   write (6, *) 'Converged eigenvalues: ',&
      converged_eigenvalues
   write (6, *) 'Total matrix-vector products: ',&
      total_matvecs
   end if

   if (converged_eigenvalues > 0) then
   if (nid == 0) then
   write (6, *) 'Exporting modes...'
   end if
   call outpost_ks(vals, vecs, Q, residual,&
      converged_eigenvalues, total_matvecs)
   end if

   if (nid == 0) write (6, *) 'Eigenproblem solver finished.'

      !--> Deallocation.
   if (allocated(Q)) deallocate (Q)
   if (allocated(H)) deallocate (H)
   if (allocated(b_vec)) deallocate (b_vec)
   if (allocated(vals)) deallocate (vals)
   if (allocated(vecs)) deallocate (vecs)
   if (allocated(residual)) deallocate (residual)

   end subroutine krylov_schur

      !-----------------------------------------------------------------------
      ! outpost_ks -- Output converged eigenmodes and spectrum files
      !
      ! Purpose:
      !   Reconstructs eigenmodes from Krylov basis via dgemv, normalizes
      !   them, and writes real/imaginary parts plus spectrum data files.
      !-----------------------------------------------------------------------
   subroutine outpost_ks(vals, vecs, Q, residual, converged, &
      total_matvecs)
    use krylov_subspace

       ! Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix
   complex(nekStab_dp), dimension(k_dim), intent(in) :: vals
   complex(nekStab_dp), dimension(k_dim, k_dim), intent(in) :: vecs

      ! Krylov basis V for the projection M*V = V*H
   type(krylov_vector), dimension(k_dim + 1), intent(in) :: Q
   real, dimension(k_dim), intent(in) :: residual
   integer, intent(in) :: converged, total_matvecs

      ! Krylov vectors
   type(krylov_vector) :: qq
   type(krylov_vector) :: ff
      ! Arrays for Krylov basis (heap-allocated to avoid stack overflow)
   real, allocatable :: qx(:,:), qy(:,:), qz(:,:)
   real, allocatable :: qp(:,:)
   real, allocatable :: qt(:,:,:)

      ! Arrays for the storage/output of a given eigenmode of the NS operator
   complex(nekStab_dp), allocatable :: fp_cx(:), fp_cy(:), fp_cz(:)
   complex(nekStab_dp), allocatable :: fp_cp(:)
   complex(nekStab_dp), allocatable :: fp_ct(:,:)

      ! Work arrays for dgemv (real/imag split)
   real, allocatable :: vecs_re(:), vecs_im(:)
   real, allocatable :: work_re(:), work_im(:)

      ! Miscellaneous variables
   integer :: i, m
   real :: speriod, trim, spurious_tol
   real :: alpha, alpha_r, alpha_i, beta, old_uparam1, omega
      ! File handling variables
   character(len=80) filename
   character(len=20) fich1, fich2, fich3, fich4, fmt2, fmt3, fmt4, fmt5, fmt6
   character(len=3) nRe, nIm, nRv
   character(len=2) ci
   integer :: outp

   nv = nx1*ny1*nz1*nelv

      ! Allocate on heap (avoids stack overflow for large 3D cases)
   allocate(qx(lv, k_dim), qy(lv, k_dim))
   allocate(qz(lv, k_dim))
   allocate(qp(lp, k_dim))
   allocate(qt(lt, ldimt, k_dim))
   allocate(fp_cx(lv), fp_cy(lv), fp_cz(lv))
   allocate(fp_cp(lp), fp_ct(lt, ldimt))
   allocate(vecs_re(k_dim), vecs_im(k_dim))
   allocate(work_re(lv), work_im(lv))
   speriod = dt*nsteps ! sampling period

      !  evop (evolution operator) defined in matvec.f90
   nRe = trim(evop)//'Re'
   nIm = trim(evop)//'Im'
   nRv = trim(evop)//'Rv'

   fich1 = 'Spectre_H'//trim(evop)//'.dat'
   fich2 = 'Spectre_NS'//trim(evop)//'.dat'
   fich3 = 'Spectre_NS'//trim(evop)//'_conv.dat'
   fich4 = 'Spectre_H'//trim(evop)//'_conv.dat'

   if (nid == 0) then
   open (unit=10, file=fich1, form='formatted', status='unknown')
   open (unit=20, file=fich2, form='formatted', status='unknown')
   open (unit=30, file=fich3, form='formatted', status='unknown')
   open (unit=40, file=fich4, form='formatted', status='unknown')
   end if

      ! Copy Krylov basis into contiguous arrays for BLAS
   do i = 1, k_dim
   qx(:, i) = Q(i)%vx(:)
   qy(:, i) = Q(i)%vy(:)
   if (if3D) qz(:, i) = Q(i)%vz(:)
   if (ifpo) qp(:, i) = Q(i)%pr(:)
   if (ifto) qt(:, 1, i) = Q(i)%t(:, 1)
   if (ldimt > 1) then
   do m = 2, ldimt
   if (ifpsco(m - 1)) qt(:, m, i) = Q(i)%t(:, m)
   end do
   end if
   end do

      ! outpost full spectrum with convergence flag (4th column: 1=converged, 0=not)
   do i = 1, k_dim
   if (nid == 0) then
   if (residual(i) < eigen_tol&
      * max(abs(vals(i)), 1.0d0)) then
   write (10, "(3E15.7,I2)") real(vals(i)),&
      aimag(vals(i)), residual(i), 1
   write (20, "(3E15.7,I2)")&
      real(log_transform(vals(i)))/speriod,&
      aimag(log_transform(vals(i)))/speriod,&
      residual(i), 1
   else
   write (10, "(3E15.7,I2)") real(vals(i)),&
      aimag(vals(i)), residual(i), 0
   write (20, "(3E15.7,I2)")&
      real(log_transform(vals(i)))/speriod,&
      aimag(log_transform(vals(i)))/speriod,&
      residual(i), 0
   end if
   end if
   end do

   spurious_tol = max(param(21), param(22))
   outp = 0 ! outposted modes counter
   do i = 1, converged ! loop over converged modes

   if (outp >= maxmodes) then

   if (nid == 0) write (6, *) 'maxmodes reached, skipping converged mode ', i
   cycle ! skip nonconverged modes or if too many modes have been outposted

   else ! converged modes

      !     ----- Computation of eigenmode via dgemv (real/imag split) -----
   vecs_re(:) = real(vecs(:, i))
   vecs_im(:) = aimag(vecs(:, i))

   call dgemv('N', nv, k_dim, 1.0d0, qx, lv, vecs_re, 1, 0.0d0, work_re, 1)
   call dgemv('N', nv, k_dim, 1.0d0, qx, lv, vecs_im, 1, 0.0d0, work_im, 1)
   fp_cx(:) = dcmplx(work_re, work_im)

   call dgemv('N', nv, k_dim, 1.0d0, qy, lv, vecs_re, 1, 0.0d0, work_re, 1)
   call dgemv('N', nv, k_dim, 1.0d0, qy, lv, vecs_im, 1, 0.0d0, work_im, 1)
   fp_cy(:) = dcmplx(work_re, work_im)

   if (if3D) then
   call dgemv('N', nv, k_dim, 1.0d0, qz, lv, vecs_re, 1, 0.0d0, work_re, 1)
   call dgemv('N', nv, k_dim, 1.0d0, qz, lv, vecs_im, 1, 0.0d0, work_im, 1)
   fp_cz(:) = dcmplx(work_re, work_im)
   end if

   if (ifpo) fp_cp(:) = matmul(qp(:, 1:k_dim), vecs(:, i))

   if (ifto) fp_ct(:, 1) = matmul(qt(:, 1, 1:k_dim), vecs(:, i))
   if (ldimt > 1) then
   do m = 2, ldimt
   if (ifpsco(m - 1)) fp_ct(:, m) = matmul(qt(:, m, 1:k_dim), vecs(:, i))
   end do
   end if

      !        normalization to unit-norm (volume integral of FP*conj(FP) = 1.)
   call norm(real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct), alpha_r)
   call norm(aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), aimag(fp_ct), alpha_i)

   if (nid == 0) write (6, *) 'Checking eigenvector', i
   if (nid == 0) write (6, *)
   if (nid == 0) write (6, *) '       norm Re/Im:', alpha_r, alpha_i

   alpha = alpha_r**2 + alpha_i**2
   if (alpha < 1e-60) then
   if (nid == 0) write (6, *) 'WARNING: zero-norm eigenvector', i
   cycle
   end if
   beta = 1.0d0/sqrt(alpha)

      !call norm_grad(real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct), norma_Re)
      !call norm_grad(aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), aimag(fp_ct), norma_Im)
      !if (nid == 0) write (6, *) '  grad norm Re/Im:', norma_Re, norma_Im
      !if (nid == 0) write (6, *)

      !    if (norma_Re > 1.1 .or. norma_Im > 1.1) then
      !       if (nid == 0) write (6, *) ' Skipping spurious (non-physical) eigenvector:', i, real(vals(i)), aimag(vals(i))
      !       cycle  ! skip this iteration if all real parts are zero
      !    end if !

   if (nid == 0) then
   write (6, '(A, I0, A, I0)') 'Outposting eigenvector: ', i, ' / ', maxmodes
   write (6, '(A, E15.7)') '  sigma = ', real(log_transform(vals(i)))/speriod
   omega = aimag(log_transform(vals(i)))/speriod
   write (6, '(A, E15.7)') '  omega = ', omega
   write (6, '(A, E15.7)') '      f = ', omega/(2.0d0*NEKSTAB_PI)
   end if
   outp = outp + 1

   if (nid == 0) then
   write (30, "(2E15.7)") real(log_transform(vals(i)))/speriod, aimag(log_transform(vals(i)))/speriod
   write (40, "(2E15.7)") real(vals(i)), aimag(vals(i))
   end if

   time = real(outp) !outp files are numbered from 1 to k_dim

      !     ----- Output the real part -----
   call nopcopy(vx, vy, vz, pr, t, real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct))
   call nopcmult(vx, vy, vz, pr, t, beta)
   call outpost2(vx, vy, vz, pr, t, nof, nRe)
   call outpost_vort(vx, vy, vz, nRv)

      !     ----- Output the imaginary part -----
   call nopcopy(vx, vy, vz, pr, t, aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), aimag(fp_ct))
   call nopcmult(vx, vy, vz, pr, t, beta)
   call outpost2(vx, vy, vz, pr, t, nof, nIm)

      !     computing and outposting optimal response from real part ! works with Floquet!
   if ((uparam(1) == 3.3 .or. uparam(1) == 3.31)) then
   old_uparam1 = uparam(1)
   if (uparam(1) == 3.3) uparam(1) = 3.1 ! changing to linearized solver !
   if (uparam(1) == 3.31) uparam(1) = 3.11 ! changing to linearized solver in Floquet
   call bcast(uparam(1), wdsize)
   call nopcopy(ff%vx, ff%vy, ff%vz, ff%pr, ff%t, real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct))
   call matvec(qq, ff) ! baseflow already in ubase
   call outpost2(qq%vx, qq%vy, qq%vz, qq%pr, qq%t, nof, 'ore')
   call outpost_vort(qq%vx, qq%vy, qq%vz, 'orv')
   uparam(1) = old_uparam1
   call bcast(uparam(1), wdsize)
   end if ! uparam(1).eq.3.3.or.uparam(1).eq.3.31
   end if

   end do ! i=1, k_dim

   if (nid == 0) then

   close (10); close (20); close (30); close (40)
      !
   fmt2 = '(A,I16)'
   fmt3 = '(A,F16.4)'
   fmt4 = '(A,F16.12)'
   fmt5 = '(A,E15.7)' ! max precision
   fmt6 = '(A,E13.4)' ! same as hmhlz
      !
   filename = 'Spectre_'//trim(evop)//'.info'
   open (844, file=filename, action='write', status='replace')

   write (844, '(A,A)') 'Nek5000 version:', NVERSION
   write (844, '(A,A)') 'nekStab version:', NSVERSION
   write (844, '(A)') '[mesh]'
   write (844, fmt2) 'lx1=             ', lx1
   write (844, fmt2) 'polyOrder N=     ', lx1 - 1
   write (844, fmt2) 'tot elemts=      ', nelgv
   write (844, fmt2) 'tot points=      ', nelgv*(lx1)**ldim
   write (844, fmt2) 'MPI ranks=       ', np
   write (844, fmt2) 'e/rank=          ', nelgv/np
   write (844, '(A)') '[userParams]'
   do i = 1, 10
   write (ci, '(I2.2)') i
   write (844, fmt4) 'uparam'//ci//'=        ',uparam(i)
   end do
   write (844, '(A)') '[solver]'
   write (844, fmt3) 'ctarg=           ', ctarg
   write (844, fmt2) 'nsteps=          ', nsteps
   write (844, fmt5) 'dt=              ', dt
   write (844, fmt3) 'Re=              ', 1.0/param(2)
   write (844, fmt6) 'residualTol PRE= ', param(21)
   write (844, fmt6) 'residualTol VEL= ', param(22)
   if (ifheat) then
   write (844, fmt6) 'residualTol TEM= ', param(22)
   write (844, fmt3) 'Pe=              ', 1.0/param(8)
   end if
   write (844, '(A)') '[eigensolver]'
   write (844, fmt4) 'sampling period =', speriod
   write (844, fmt2) 'k_dim=           ', k_dim
   write (844, fmt6) 'eigentol=        ', eigen_tol
   write (844, fmt2) 'schur_target=    ', schur_tgt
   write (844, fmt3) 'schur_del=       ', schur_del
   write (844, fmt2) 'schur iterations=', schur_cnt
   write (844, fmt2) 'total matvecs=   ', total_matvecs
   write (844, fmt2) 'converged=       ', converged
   write (844, fmt2) 'outp=            ', outp
   close (844)
   end if

   deallocate(qx, qy, qz, qp, qt)
   deallocate(fp_cx, fp_cy, fp_cz, fp_cp, fp_ct)
   deallocate(vecs_re, vecs_im, work_re, work_im)
   end subroutine outpost_ks

      !-----------------------------------------------------------------------

   subroutine schur_condensation(mstart, H, Q, ksize)

      !     Krylov-Schur restart (Stewart, SIMAX 2001).
      !
      !     Given a k-step Arnoldi factorization  A * Q = Q * H + f * e_k^T,
      !     this subroutine compresses it to a p-step Krylov-Schur factorization
      !     by keeping only the p "wanted" Ritz pairs and discarding the rest.
      !     Arnoldi then resumes from position p+1.
      !
      !     Algorithm outline:
      !     1. Compute real Schur decomposition:  H = V * T * V^T
      !     2. Select p wanted eigenvalues (near unit circle + nev+4 largest)
      !     3. Reorder Schur form so wanted eigenvalues are in T(1:p, 1:p)
      !     4. Rotate Krylov basis:  Q_new(:,1:p) = Q_old(:,1:k) * V(:,1:p)
      !     5. Set coupling row:  H(p+1, 1:p) = beta * e_k^T * V(:, 1:p)
      !     6. Copy Q(k+1) -> Q(p+1) as new starting vector
      !
      !     Memory optimization (this revision):
      !     - Fields are rotated one at a time (vx, then vy, then vz, ...)
      !       to reduce peak memory from O(nfields * nv * k) to O(nv * k).
      !     - Only p (= nsel) output columns are computed via dgemm, not
      !       all k columns.  FLOP savings: O(nv * k * (k - p)) per field.
      !     - Explicit out-of-place dgemm avoids hidden temporaries that
      !       Fortran's matmul would create for in-place A = matmul(A, B).
      !
      !     INPUTS / OUTPUTS
      !     ----------------
      !     mstart (inout) : on entry, previous restart index (from uparam(2))
      !                      on exit, set to nsel + 1 (restart position)
      !     H (inout)      : (ksize+1, ksize) Hessenberg matrix from Arnoldi;
      !                      overwritten with condensed Schur form
      !     Q (inout)      : (ksize+1) array of krylov_vectors;
      !                      Q(1:nsel) overwritten with rotated Schur basis,
      !                      Q(nsel+1) = former Q(ksize+1) (residual vector)
      !     ksize (in)     : dimension of the Krylov subspace (= k_dim)
      !
    use krylov_subspace

    integer, intent(inout) :: mstart
   integer, intent(in) :: ksize

      !     ----- Krylov basis V for the projection A*V = V*H + f*e_k^T -----
   type(krylov_vector), dimension(ksize + 1), intent(inout) :: Q

      !     ----- Upper Hessenberg matrix -----
   real, dimension(ksize + 1, ksize), intent(inout) :: H

      !     b_vec = H(k+1,k) * e_k  (residual coupling row before rotation)
   real, dimension(ksize) :: b_vec

      !     ----- Eigenvalues and Schur vectors -----
      !     vals  : eigenvalues of H(1:k, 1:k) from Schur decomposition
      !     vecs  : orthogonal Schur vectors (columns of V in H = V*T*V^T)
      !     selected : logical mask for wanted eigenvalues
   complex(nekStab_dp), dimension(ksize) :: vals
   real, dimension(ksize, ksize) :: vecs
   logical, dimension(ksize) :: selected

      !     ----- Workspace for field-by-field basis rotation -----
      !     basis(lv, ksize)  : contiguous copy of one field from Q(1:k)
      !     rotated(lv, nsel) : result of dgemm (out-of-place, no aliasing)
      !
      !     Processing one field at a time means the velocity workspace is
      !     reused for vx, vy, vz in sequence.  Peak memory is O(nv * k)
      !     + O(nv * nsel) regardless of the number of fields (vx/vy/vz/pr/t).
      !     The old code allocated ALL fields simultaneously: O(3*nv*k) just
      !     for velocity, plus Fortran's matmul created hidden temporaries.
   real, allocatable :: basis(:,:), rotated(:,:)

      !     ----- Miscellaneous -----
   integer :: i, j, m, nsel, n_locked
   real :: res
   real, dimension(ksize) :: time_old, schur_res

   nv = nx1*ny1*nz1*nelv

      !     =====================================================
      !     Step 1 : Extract residual coupling and Schur decompose
      !     =====================================================

      !     b_vec captures the only nonzero in row k+1 of the Arnoldi relation.
      !     After rotation: b_new = b_vec * V, giving the coupling between
      !     the residual vector and the new Schur basis.
   b_vec = 0.0d0; b_vec(ksize) = H(ksize + 1, ksize)
   vals = (0.0d0, 0.0d0); vecs = 0.0d0

      !     Compute H(1:k,1:k) = V * T * V^T  via LAPACK dgees.
      !     On return: H is overwritten with T (quasi-upper-triangular),
      !                vecs contains V (orthogonal Schur vectors),
      !                vals contains eigenvalues.
   call schur(H(1:ksize, 1:ksize), vecs, vals, ksize)

      !     =====================================================
      !     Step 2 : Select wanted eigenvalues
      !     =====================================================

      !     Compute per-Schur-vector residuals for adaptive selection.
      !     Since b_vec = [0,...,0, H(k+1,k)], the dot product with
      !     Schur vector j reduces to: beta * V(k, j).
   do i = 1, ksize
   schur_res(i) = abs(b_vec(ksize) * vecs(ksize, i))
   end do

      !     Adaptive selection: keep eigenvalues near unit circle,
      !     partially converged (residual < sqrt(tol)), and at least
      !     nev+2 by magnitude.  Conjugate pairs kept together.
   call select_eigenvalues(selected, mstart, vals,&
      schur_res, schur_del, schur_tgt, ksize)
   nsel = mstart
   if (nid == 0) write (6, *) nsel,&
      'Ritz eigenpairs have been selected.'

      !     =====================================================
      !     Step 3 : Reorder Schur form (LAPACK dtrsen)
      !     =====================================================

      !     After reordering, the quasi-upper-triangular T has the nsel wanted
      !     eigenvalues in T(1:nsel, 1:nsel) and the unwanted in the lower-right.
      !     The Schur vectors V are reordered correspondingly.
   call ordschur(H(1:ksize, 1:ksize), vecs, selected, ksize)

      !     =====================================================
      !     Step 4 : Truncate Hessenberg matrix
      !     =====================================================

      !     Zero the off-diagonal coupling between wanted and unwanted subspaces,
      !     and all rows below nsel (these will be filled by the next Arnoldi run).
   H(1:nsel, nsel + 1:ksize) = 0.0d0
   H(nsel + 1:ksize + 1, :) = 0.0d0

      !     =====================================================
      !     Step 5 : Residual coupling in the new basis
      !     =====================================================

      !     The Arnoldi residual transforms as:  f * e_k^T * V.
      !     Since b_vec = beta * e_k, the new coupling row is b_vec * V.
      !     Only the first nsel entries are needed (the rest couple to zeroed blocks).
      !     This row goes into H(nsel+1, 1:nsel); Arnoldi will overwrite
      !     H(nsel+1, nsel+1:k) when it resumes.
   H(nsel + 1, 1:nsel) = matmul(b_vec, vecs(:, 1:nsel))

      !     =====================================================
      !     Step 5b: Lock converged eigenvalues (zero coupling)
      !     =====================================================
      !
      !     Walk the quasi-upper-triangular T(1:nsel, 1:nsel) to
      !     identify converged eigenvalues.  A 2x2 block at (j,j)
      !     has H(j+1,j) /= 0 (complex conjugate pair); otherwise
      !     it's a 1x1 block (real eigenvalue).
      !
      !     If the coupling entry |H(nsel+1, j)| < eigen_tol (or
      !     norm2 for a 2x2 block), the eigenvalue is converged.
      !     Zeroing the coupling entry "locks" it: the eigenvalue
      !     remains in the Schur basis but Arnoldi no longer
      !     perturbs it through the residual vector.  This prevents
      !     converged eigenvalues from drifting after restarts.
      !
      !     This is "option 2" from Stewart (2001): zero coupling
      !     row entries only, without full block decoupling.
   n_locked = 0
   j = 1
   do while (j <= nsel)
   if (j < nsel .and.&
      H(j+1, j) /= 0.0d0) then
      !        2x2 block (complex conjugate pair).
   res = sqrt(H(nsel+1, j)**2 + H(nsel+1, j+1)**2)
   if (res < eigen_tol) then
   H(nsel+1, j)   = 0.0d0
   H(nsel+1, j+1) = 0.0d0
   n_locked = n_locked + 2
   end if
   j = j + 2
   else
      !        1x1 block (real eigenvalue).
   if (abs(H(nsel+1, j)) < eigen_tol) then
   H(nsel+1, j) = 0.0d0
   n_locked = n_locked + 1
   end if
   j = j + 1
   end if
   end do
   if (nid == 0 .and. n_locked > 0) then
   write(6, *) n_locked,&
      'eigenvalues locked (coupling zeroed).'
   end if

      !     =====================================================
      !     Step 6 : Basis rotation via dgemm (field-by-field)
      !     =====================================================
      !
      !     The key operation is:  Q_new(:, j) = sum_i Q(:, i) * V(i, j)
      !     for j = 1, ..., nsel.  In matrix form:
      !
      !        Q_new(:, 1:nsel) = Q(:, 1:ksize) * V(:, 1:nsel)
      !
      !     This is a matrix-matrix product computed via BLAS dgemm.
      !
      !     Why dgemm instead of matmul?
      !     - matmul(A, B) where A is (lv, k) creates a hidden temporary of
      !       size (lv, k) before assigning back.  For lv = 10M, k = 100,
      !       that's ~8 GB per field — invisible in the source but real.
      !     - dgemm writes directly into a separate output buffer (rotated),
      !       avoiding the hidden allocation entirely.
      !
      !     Why field-by-field?
      !     - The krylov_vector derived type stores vx, vy, vz, pr, t as
      !       separate arrays.  They cannot be passed as a single contiguous
      !       block to dgemm without copying.
      !     - By processing one field at a time, we need only ONE workspace
      !       of size (lv, ksize) instead of THREE (or more) simultaneously.
      !       Peak memory: ~(lv * ksize + lv * nsel) * 8 bytes total,
      !       reused for each field in turn.
      !
      !     Why only nsel columns?
      !     - After reordering, columns nsel+1:ksize of V correspond to the
      !       unwanted eigenvalues.  Those rotated Krylov vectors are never
      !       used — Arnoldi will overwrite positions nsel+1:ksize anyway.
      !     - Savings: for ksize=100, nsel=34, we compute 34 columns
      !       instead of 100 → ~3x fewer FLOPs in the rotation.
      !
      !     dgemm calling convention:
      !       call dgemm('N','N', M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
      !       C(M,N) = alpha * A(M,K) * B(K,N) + beta * C(M,N)
      !     Here: M = nv (active grid points), N = nsel, K = ksize,
      !           LDA = lv (compile-time leading dimension of basis),
      !           LDB = ksize (leading dimension of vecs),
      !           LDC = lv (compile-time leading dimension of rotated).

      !     --- Rotate velocity fields (vx, vy, vz) ---
      !     Allocate workspace: basis holds one field from ALL k vectors,
      !     rotated holds the nsel wanted columns of the product.
   allocate(basis(lv, ksize), rotated(lv, nsel))

      !     vx: gather from Q(1:k) into contiguous array, rotate, scatter back.
      !     After scatter, Q(1:nsel)%vx contains the rotated data.
      !     Q(nsel+1:ksize)%vx still has old data (irrelevant, will be overwritten).
   do i = 1, ksize
   basis(:, i) = Q(i)%vx(:)
   end do
   call dgemm('N', 'N', nv, nsel, ksize,&
      1.0d0, basis, lv, vecs, ksize,&
      0.0d0, rotated, lv)
   do i = 1, nsel
   Q(i)%vx(:) = rotated(:, i)
   end do

      !     vy: same pattern, reusing basis and rotated arrays.
      !     Safe because we only modify Q(i)%vy, not Q(i)%vx (already done).
   do i = 1, ksize
   basis(:, i) = Q(i)%vy(:)
   end do
   call dgemm('N', 'N', nv, nsel, ksize,&
      1.0d0, basis, lv, vecs, ksize,&
      0.0d0, rotated, lv)
   do i = 1, nsel
   Q(i)%vy(:) = rotated(:, i)
   end do

      !     vz: only for 3D problems.
   if (if3D) then
   do i = 1, ksize
   basis(:, i) = Q(i)%vz(:)
   end do
   call dgemm('N', 'N', nv, nsel, ksize,&
      1.0d0, basis, lv, vecs, ksize,&
      0.0d0, rotated, lv)
   do i = 1, nsel
   Q(i)%vz(:) = rotated(:, i)
   end do
   end if

      !     Free velocity workspace before allocating pressure (different size).
   deallocate(basis, rotated)

      !     --- Rotate pressure field ---
      !     Pressure uses PN-2/PN discretization: different grid size lp, n2.
   if (ifpo) then
   n2 = nx2*ny2*nz2*nelv
   allocate(basis(lp, ksize), rotated(lp, nsel))
   do i = 1, ksize
   basis(:, i) = Q(i)%pr(:)
   end do
   call dgemm('N', 'N', n2, nsel, ksize,&
      1.0d0, basis, lp, vecs, ksize,&
      0.0d0, rotated, lp)
   do i = 1, nsel
   Q(i)%pr(:) = rotated(:, i)
   end do
   deallocate(basis, rotated)
   end if

      !     --- Rotate temperature and passive scalars ---
      !     krylov_vector%t is (lt, ldimt) where lt = lx1*ly1*lz1*lelt.
      !     For CHT (lelt > lelv), lt > lv and nt > nv; workspace must
      !     use lt as the leading dimension to match the array stride.
   if (ifto .or. ldimt > 1) then
   nt = nx1*ny1*nz1*nelt
   allocate(basis(lt, ksize), rotated(lt, nsel))
   end if

      !     Temperature (first scalar field).
   if (ifto) then
   do i = 1, ksize
   basis(:, i) = Q(i)%t(:, 1)
   end do
   call dgemm('N', 'N', nt, nsel, ksize,&
      1.0d0, basis, lt, vecs, ksize,&
      0.0d0, rotated, lt)
   do i = 1, nsel
   Q(i)%t(:, 1) = rotated(:, i)
   end do
   end if

      !     Additional passive scalars (ifpsco flags which are active).
      !     Reuses the same basis/rotated workspace for each scalar.
   if (ldimt > 1) then
   do m = 2, ldimt
   if (ifpsco(m - 1)) then
   do i = 1, ksize
   basis(:, i) = Q(i)%t(:, m)
   end do
   call dgemm('N', 'N', nt, nsel, ksize,&
      1.0d0, basis, lt, vecs, ksize,&
      0.0d0, rotated, lt)
   do i = 1, nsel
   Q(i)%t(:, m) = rotated(:, i)
   end do
   end if
   end do
   end if

   if (allocated(basis)) deallocate(basis, rotated)

      !     --- Rotate time component (Newton periodic orbits) ---
      !     For Newton PO (uparam(1)==2.1), the krylov_vector includes a
      !     scalar time component in the inner product.  It must be rotated
      !     along with the spatial fields to keep the basis consistent.
      !     We save the old values first since the rotation reads from all
      !     k vectors but writes to only nsel of them.
   if (isNewtonPO) then
   do i = 1, ksize
   time_old(i) = Q(i)%time
   end do
   do i = 1, nsel
   Q(i)%time = dot_product(time_old, vecs(:, i))
   end do
   end if

      !     =====================================================
      !     Step 7 : Set up for Arnoldi restart
      !     =====================================================
      !
      !     After condensation, the factorization is:
      !        A * Q(:,1:nsel) = Q(:,1:nsel) * H(1:nsel,1:nsel)
      !                        + Q(:,nsel+1) * H(nsel+1, 1:nsel)
      !     where Q(:,nsel+1) is the last generated Krylov vector (residual
      !     direction).  Arnoldi resumes from mstart = nsel + 1, extending
      !     the factorization back to ksize columns.
   mstart = nsel + 1
   call k_copy(Q(mstart), Q(ksize + 1))

   end subroutine schur_condensation

      !-----------------------------------------------------------------------

   subroutine select_eigenvalues(selected, nsel, vals, &
      residuals, delta, nev, n)

      !     Adaptive eigenvalue selection for Krylov-Schur restart.
      !
      !     Selects which eigenvalues to keep during Schur condensation
      !     using three complementary criteria:
      !
      !     1. Near unit circle: |lambda| >= 1 - delta
      !        These are the stability-relevant eigenvalues.
      !
      !     2. Partially converged: Schur residual < sqrt(eigen_tol)
      !        If eigen_tol = 1e-6, this keeps eigenvalues with
      !        residual < 1e-3.  Discarding these wastes ~100 matvecs
      !        of convergence progress per eigenvalue.
      !
      !     3. Minimum guarantee: at least nev + 2 by magnitude
      !        Ensures the target eigenvalues are always retained.
      !
      !     The selection is bounded:  nev+2 <= nsel <= n - nev
      !     to leave room for Arnoldi to make progress each cycle.
      !     Conjugate pairs are kept together (real Schur form).
      !
      !     INPUTS
      !     ------
      !     vals      : complex(n)  — eigenvalues from Schur decomposition
      !     residuals : real(n)     — per-Schur-vector Arnoldi residuals
      !     delta     : real        — magnitude selection radius (schur_del)
      !     nev       : integer     — number of desired eigenvalues (schur_tgt)
      !     n         : integer     — total number of eigenvalues (k_dim)
      !
       !     OUTPUTS
       !     -------
       !     selected  : logical(n)  — which eigenvalues to keep
       !     nsel      : integer     — count of selected eigenvalues

       !     ----- Arguments -----
   integer, intent(in) :: nev, n
   complex(nekStab_dp), dimension(n), intent(in) :: vals
   real, dimension(n), intent(in) :: residuals
   real, intent(in) :: delta
   logical, dimension(n), intent(out) :: selected
   integer, intent(out) :: nsel

      !     ----- Local variables -----
   integer :: i, min_keep, max_keep
   integer, dimension(n) :: idx
   real, dimension(n) :: work_arr
   real :: sqrt_tol
   integer :: n_circle, n_resid

   sqrt_tol = sqrt(eigen_tol)
   selected = .false.

      !     --> Criterion 1: eigenvalues near the unit circle.
   do i = 1, n
   if (abs(vals(i)) >= (1.0d0 - delta))&
      selected(i) = .true.
   end do
   n_circle = count(selected)

      !     --> Criterion 2: partially converged eigenvalues.
      !         Residual < sqrt(tol) means close to convergence;
      !         discarding these wastes the matvecs already invested.
   do i = 1, n
   if (residuals(i) < sqrt_tol)&
      selected(i) = .true.
   end do
   n_resid = count(selected) - n_circle

      !     --> Criterion 3: guarantee at least nev+2 by magnitude.
      !         Use argsort (ascending) — largest magnitudes at end.
   min_keep = nev + 2
   if (count(selected) < min_keep) then
   do i = 1, n
   idx(i) = i
   end do
   work_arr = abs(vals)
   call argsort(n, work_arr, idx)
   do i = n, 1, -1
   if (count(selected) >= min_keep) exit
   selected(idx(i)) = .true.
   end do
   end if

      !     --> Ensure conjugate pairs are not split.
   call ensure_conjugate_pairs(selected, vals, n)

      !     --> Cap: keep at most n - nev to leave room for Arnoldi.
      !         If over budget, drop eigenvalues with largest residual.
   max_keep = n - nev
   if (count(selected) > max_keep) then
   do i = 1, n
   idx(i) = i
   end do
   work_arr = residuals
   call argsort(n, work_arr, idx)
      !        Rebuild selection: keep max_keep with smallest residuals.
   selected = .false.
   do i = 1, max_keep
   selected(idx(i)) = .true.
   end do
   call ensure_conjugate_pairs(selected, vals, n)
   end if

   nsel = count(selected)

   if (nid == 0) then
   write(6, '(A,I4,A,I4,A)')&
      ' Adaptive restart: keeping ', nsel,&
      ' of ', n, ' Ritz values.'
   write(6, '(A,I4,A,I4,A,I4)')&
      '   unit circle: ', n_circle,&
      ', near-converged: ', n_resid,&
      ', floor: ', min_keep
   end if

   return
   end subroutine select_eigenvalues

      !-----------------------------------------------------------------------

   subroutine ensure_conjugate_pairs(selected, vals, n)

      !     Ensure that if one eigenvalue of a complex conjugate pair
      !     is selected, the other is also selected.  In the Schur
      !     decomposition from dgees, conjugate pairs are consecutive:
      !     vals(j) and vals(j+1) have equal real parts and opposite
      !     imaginary parts.  We walk consecutive pairs and enforce
       !     symmetric selection to preserve the real Schur form.

    integer, intent(in) :: n
   complex(nekStab_dp), dimension(n), intent(in) :: vals
   logical, dimension(n), intent(inout) :: selected

   integer :: i
   real :: real_diff, imag_sum, scale, tol
   real, parameter :: REL_TOL = 1.0d-10

   i = 1
   do while (i < n)
      !     --> Check if vals(i) and vals(i+1) form a conjugate pair.
   if (aimag(vals(i)) /= 0.0d0) then
   scale = max(abs(vals(i)), 1.0d0)
   tol = REL_TOL * scale
   real_diff = abs(real(vals(i)) - real(vals(i+1)))
   imag_sum = abs(aimag(vals(i)) + aimag(vals(i+1)))
   if (real_diff < tol .and. imag_sum < tol) then
      !           Conjugate pair: select both if either is selected.
   if (selected(i) .or. selected(i+1)) then
   selected(i) = .true.
   selected(i+1) = .true.
   end if
   i = i + 2
   cycle
   end if
   end if
   i = i + 1
   end do

   return
   end subroutine ensure_conjugate_pairs

   end module nekstab_eigensolvers
