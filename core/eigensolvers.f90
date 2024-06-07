      !-----------------------------------------------------------------------
      
      subroutine inner_product(alpha, px, py, pz, pp, pt, qx, qy, qz, qp, qt)
      
      !     This function provides the user-defined inner product to be used throughout
      !     the computation.
      !
      !     INPUTS
      !     ------
      !
      !     px, py, pz : nek arrays of size lv = lx1*ly1*lz1*lelv.
      !     Arrays containing the velocity fields of the first vector.
      !
      !     pp : nek array of size lp = lx2*ly2*lz2*lelt
      !     Array containing the pressure field of the first vector.
      !
      !     qx, qy, qz : nek arrays of size lv = lx1*ly1*lz1*lelv.
      !     Arrays containing the velocity fields of the second vector.
      !
      !     qp : nek array of size lp = lx2*ly2*lz2*lelt
      !     Array containing the pressure field of the second vector.
      !
      !     RETURN
      !     ------
      !
      !     alpha : real
      !     Value of the inner-product alpha = <p, q>.
      !
         use krylov_subspace
         implicit none
         include "SIZE"
         include "TOTAL"
      
         real, dimension(lv), intent(in) :: px, py, pz
         real, dimension(lv), intent(in) :: qx, qy, qz
         real, dimension(lp), intent(in) :: pp, qp !not used
         real, dimension(lv, ldimt), intent(in) :: pt, qt
      
         real, intent(out) :: alpha
         real :: glsc3
         integer m
      
         nv = nx1*ny1*nz1*nelv
         nt = nx1*ny1*nz1*nelt
      
         alpha = glsc3(px, bm1s, qx, nv) + glsc3(py, bm1s, qy, nv)
         if (if3D) alpha = alpha + glsc3(pz, bm1s, qz, nv)
         if (ifto) alpha = alpha + glsc3(pt(:, 1), bm1s, qt(:, 1), nt)
         if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) alpha = alpha + glsc3(pt(:, m), bm1s, qt(:, m), nt)
         end do
         end if
      
         return
      end subroutine inner_product
      
      !--------------------------------------------------------------------------
      
      subroutine norm(qx, qy, qz, qp, qt, alpha) ! Compute vector norm
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, intent(in), dimension(lv) :: qx, qy, qz
         real, intent(in), dimension(lp) :: qp
         real, intent(in), dimension(lv, ldimt) :: qt
         real, intent(out) :: alpha
      
         call inner_product(alpha, qx, qy, qz, qp, qt, qx, qy, qz, qp, qt)
         alpha = sqrt(alpha)
      
      end subroutine norm
      
      !----------------------------------------------------------------------
      
      subroutine normalize(qx, qy, qz, qp, qt, alpha)
      
      !     This function normalizes the state vector [qx, qy, qz, qp]^T where
      !     qx, qy and qz are the streamwise, cross-stream and spanwise velocity
      !     components while qp is the corresponding pressure field.
      !
      !     INPUTS / OUTPUTS
      !     ----------------
      !
      !     qx, qy, qz : nek arrays of size lv = lx1*ly1*lz1*lelv.
      !     Arrays storing the velocity components.
      !
      !     qp : nek array of size lp = lx2*ly2*lz2*lelt
      !     Array storing the corresponding pressure field.
      !
      !     alpha : real
      !     Norm of the vector.
      !
      !     Last edit : April 2nd 2020 by JC Loiseau.
      
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv), intent(inout) :: qx, qy, qz
         real, dimension(lp), intent(inout) :: qp
         real, dimension(lv, ldimt), intent(inout) :: qt
         real, intent(out) :: alpha
         real :: beta
      
      !     --> Compute the user-defined norm.
         call norm(qx, qy, qz, qp, qt, alpha)
         beta = 1.0d0/alpha
      
      !     --> Normalize the vector.
         call nopcmult(qx, qy, qz, qp, qt, beta)
      
      end subroutine normalize
      
      !-----------------------------------------------------------------------
      
      subroutine krylov_schur
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     -----Krylov basis V for the projection M*V = V*H -----
         type(krylov_vector), allocatable, dimension(:) :: Q
      
      !     ----- Upper Hessenberg matrix -----
         real, allocatable, dimension(:, :) :: H
         real, allocatable, dimension(:, :) :: b_vec
      
      !     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----
         complex(kind=kind(0.0d0)), allocatable, dimension(:) :: vals
         complex(kind=kind(0.0d0)), allocatable, dimension(:, :) :: vecs
      
         real, allocatable, dimension(:) :: residual
      
      !     ----- Miscellaneous -----
         type(krylov_vector) :: wrk, wrk2
      
         integer :: mstart, converged_eigenvalues, m
         real :: alpha
         logical :: converged
         integer :: i, j
         character(len=30) :: filename
      
      !     ----- Allocate arrays -----
         allocate (Q(k_dim + 1))
         allocate (H(k_dim + 1, k_dim), b_vec(1, k_dim), vals(k_dim), vecs(k_dim, k_dim), residual(k_dim))
      
         time = 0.0d0
         H(:, :) = 0.0d0; b_vec = 0.0d0; residual = 0.0d0
         call k_zero(Q(1:k_dim + 1))
      
      !     ----- Loading baseflow from disk (optional) -----
      
         if (ifldbf) then            !skip loading if single run
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
      
         if (istep == 0 .and. (
     $   uparam(1) == 3.11 .or. ! Floquet direct
     $   uparam(1) == 3.21 .or. ! Floquet adjoint
     $   uparam(1) == 3.31    ! Floquet direct-adjoint
     $   )) then
         param(10) = time         ! upo period in field
         if (nid == 0) write (6, *) 'Floquet mode !!!'
         if (nid == 0) write (6, *) ' getting endTime from file: endTime=', param(10)
         end if
         call bcast(param(10), wdsize)
      
      !     ----- First vector (new from noise or restart) -----
      
         if (uparam(2) == 0) then
      
            if (nid == 0) write (6, *) 'Starting first Arnoldi decomposition...'
      
      !     ----- Creates seed vector for the Krylov subspace -----
      
            if (ifseed_nois) then    ! noise as initial seed
      
               if (nid == 0) write (6, *) 'Filling fields with noise...'
               call op_add_noise(wrk2%vx, wrk2%vy, wrk2%vz)
               if (ifto) call add_noise_scal(wrk2%t(:, 1), 9.0e4, 3.0e3, 4.0e5)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call add_noise_scal(wrk2%t(:, m), 9.0e1*m, 3.0e2*m, 4.0e1*m)
               end do
               end if
               call k_normalize(wrk2, alpha)
      !    call outpost2(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%t, nof, 'NOS')
               call matvec(wrk, wrk2)
      !    call outpost2(wrk%vx, wrk%vy, wrk%vz, wrk%pr, wrk%t, nof, 'NOS')
      
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
      ! try this more efficient way
      !    do i = 1, k_dim+1
      !       read(67,"(1E15.7)") H(i, 1:k_dim)
      !    enddo
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
      
            mstart = mstart + 1        !careful here!
            call load_files(Q, mstart, k_dim + 1, 'KRY')
            if (nid == 0) write (6, *) 'Restart fields loaded to memory!'
      
         end if
      
      !     ======================================
      !     =====                            =====
      !     ===== Krylov-Schur decomposition =====
      !     =====                            =====
      !     ======================================
      
         schur_cnt = 0
         converged = .false.
      
         do while (.not. converged)
      
      !     --> Arnoldi factorization.
            call arnoldi_factorization(Q, H, mstart, k_dim, k_dim)
      
      !     if(nid.eq.0) then
      !        open(unit=12345, file="Hessenberg_matrix.dat")
      !        write(12345, *) H(1:k_dim, 1:k_dim)
      !        close(12345)
      !     endif
      
      !     --> Compute the eigenspectrum of the Hessenberg matrix.
            call eig(H(1:k_dim, 1:k_dim), vecs, vals, k_dim)
      
      !     --> Check the residual of the eigenvalues.
      
            residual = abs(H(k_dim + 1, k_dim)*vecs(k_dim, :))
      
            converged_eigenvalues = count(residual < eigen_tol)
      
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
               else                ! Apply Schur condensation before restarting the factorization.
                  schur_cnt = schur_cnt + 1
                  if (nid == 0) write (6, *) 'Starting Schur condensation phase.', schur_cnt
                  call schur_condensation(mstart, H, Q, k_dim)
               end if
      
            end select
      
         end do
      
         if (nid == 0) open (unit=99, file="orthonormality.dat")
         do i = 1, k_dim
            call k_norm(alpha, Q(i))
            if (nid == 0) write (99, '("Norm of the ", I4, "th mode = ", F20.14)') i, alpha
            do j = i + 1, k_dim
               call k_dot(alpha, q(i), q(j))
               if (nid == 0) write (99, '("Orthogonality between mode ", I4, " and mode ", I4, " = ", E15.7)') i, j, alpha
            end do
            if (nid == 0) write (99, *)
         end do
         if (nid == 0) close (unit=99)
      
         if (nid == 0) write (6, *) 'Converged eigenvalues: ', converged_eigenvalues
      
         if (converged_eigenvalues > 0) then
         if (nid == 0) then
            write (6, *) 'Exporting modes...'
      ! write (6, *) 'residual: ', residual
         end if
         call outpost_ks(vals, vecs, Q, residual, converged_eigenvalues)
         end if
      
         if (nid == 0) write (6, *) 'Eigenproblem solver finished.'
      
      !     --> Deallocation.
         deallocate (Q)
         deallocate (H, b_vec, vals, vecs)
      
      end subroutine krylov_schur
      
      !-----------------------------------------------------------------------
      
      subroutine schur_condensation(mstart, H, Q, ksize)
      
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     =================================================
      !     =====                                       =====
      !     ===== Declaration of the required variables =====
      !     =====                                       =====
      !     =================================================
      
         integer, intent(inout) :: mstart
         integer, intent(in) :: ksize
      
      !     ----- Krylov basis V for the projection M*V = V*H -----
      
         type(krylov_vector), dimension(ksize + 1) :: Q
         real, dimension(lv, ksize + 1) :: qx, qy, qz
         real, dimension(lp, ksize + 1) :: qp
         real, dimension(lt, ldimt, ksize + 1) :: qt
      
      !     ----- Upper Hessenberg matrix -----
      
         real, dimension(ksize + 1, ksize) :: H
         real, dimension(ksize) :: b_vec
      
      !     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----
      
         complex(kind=kind(0.0d0)), dimension(ksize) :: vals
      
      !     ----- Miscellaneous -----
         integer :: i, m
         logical, dimension(ksize) :: selected
      
      !     ----- Schur and Hessenberg decomposition -----
         real, dimension(ksize, ksize) :: vecs
      
      !     --> Initialize arrays.
         b_vec = 0.0d0; b_vec(ksize) = H(ksize + 1, ksize)
         vals = (0.0d0, 0.0d0); vecs = 0.0d0
      
      !     --> Perform the Schur decomposition.
         call schur(H(1:ksize, 1:ksize), vecs, vals, ksize)
      
      !     --> Partition the eigenvalues in wanted / unwanted.
         call select_eigenvalues(selected, mstart, vals, schur_del, schur_tgt, ksize)
         if (nid == 0) write (6, *) mstart, 'Ritz eigenpairs have been selected.'
      
      !     --> Re-order the Schur decomposition based on the partition.
         call ordschur(H(1:ksize, 1:ksize), vecs, selected, ksize)
      
      !     --> Zero-out the unwanted blocks of the Schur matrix.
         H(1:mstart, mstart + 1:ksize) = 0.0d0
         H(mstart + 1:ksize + 1, :) = 0.0d0
      
      !     --> Re-order the Krylov basis accordingly.
         do i = 1, k_dim + 1
            qx(:, i) = Q(i)%vx(:)
            qy(:, i) = Q(i)%vy(:)
            if (if3D) qz(:, i) = Q(i)%vz(:)
            if (ifpo) qp(:, i) = Q(i)%pr(:)
            if (ifto) qt(:, i, 1) = Q(i)%t(:, 1)
            if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) qt(:, i, m) = Q(i)%t(:, m)
            end do
            end if
         end do
         qx(:, 1:ksize) = matmul(qx(:, 1:ksize), vecs)
         qy(:, 1:ksize) = matmul(qy(:, 1:ksize), vecs)
         if (if3D) qz(:, 1:ksize) = matmul(qz(:, 1:ksize), vecs)
         if (ifpo) qp(:, 1:ksize) = matmul(qp(:, 1:ksize), vecs)
         if (ifto) qt(:, 1, 1:ksize) = matmul(qt(:, 1, 1:ksize), vecs)
         if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) qt(:, m, 1:ksize) = matmul(qt(:, m, 1:ksize), vecs)
         end do
         end if
      
      !     --> Update the Schur matrix with b.T @ Q corresponding to
      !     the residual beta in the new basis.
         b_vec = matmul(b_vec, vecs)
         H(mstart + 1, :) = b_vec
      
      !     --> Add the last generated Krylov vector as the new starting one.
         mstart = mstart + 1
      
         call nopcopy(qx(:, mstart), qy(:, mstart), qz(:, mstart), qp(:, mstart), qt(:, :, mstart),
     $   qx(:, ksize + 1), qy(:, ksize + 1), qz(:, ksize + 1), qp(:, ksize + 1), qt(:, :, ksize + 1))
         do i = 1, k_dim + 1
            Q(i)%vx(:) = qx(:, i)
            Q(i)%vy(:) = qy(:, i)
            if (if3D) Q(i)%vz(:) = qz(:, i)
            if (ifpo) Q(i)%pr(:) = qp(:, i)
            if (ifto) Q(i)%t(:, 1) = qt(:, 1, i)
            if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) Q(i)%t(:, m) = qt(:, m, i)
            end do
            end if
         end do
      
         return
      end subroutine schur_condensation
      
      !----------------------------------------------------------------------
      
      subroutine outpost_ks(vals, vecs, Q, residual, converged)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      ! Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix
         complex(kind=kind(0.0d0)), dimension(k_dim), intent(in) :: vals
         complex(kind=kind(0.0d0)), dimension(k_dim, k_dim), intent(in) :: vecs
      
      ! Krylov basis V for the projection M*V = V*H
         type(krylov_vector), dimension(k_dim + 1), intent(in) :: Q
         real, dimension(k_dim), intent(in) :: residual
         integer, intent(in) :: converged
      
      ! Krylov vectors
         type(krylov_vector) :: qq, ff
      
      ! Arrays for Krylov basis
         real, dimension(lv, k_dim + 1) :: qx, qy, qz
         real, dimension(lp, k_dim + 1) :: qp
         real, dimension(lt, ldimt, k_dim + 1) :: qt
      
      ! Arrays for the storage/output of a given eigenmode of the NS operator
         complex(kind=kind(0.0d0)), dimension(lv) :: fp_cx, fp_cy, fp_cz
         complex(kind=kind(0.0d0)), dimension(lp) :: fp_cp
         complex(kind=kind(0.0d0)), dimension(lt, ldimt) :: fp_ct
      
      ! Miscellaneous variables
         integer :: i, m
         real :: speriod, trim, spurious_tol!, glmin, glmax
         real :: alpha, alpha_r, alpha_i, beta, old_uparam1
         real :: norma_Re, norma_Im
         complex :: log_transform
         logical ifto_sav, ifpo_sav
      
      ! File handling variables
         character(len=80) :: filename
         character(len=20) :: fich1, fich2, fich3, fmt2, fmt3, fmt4, fmt5, fmt6
         character(len=3) :: nRe, nIm, nRv
         integer :: outp
      
         nv = nx1*ny1*nz1*nelv
         speriod = dt*nsteps ! sampling period
      
      !     evop defined in matrix_vector_product
         nRe = trim(evop)//'Re'
         nIm = trim(evop)//'Im'
         nRv = trim(evop)//'Rv'
      
         fich1 = 'Spectre_H'//trim(evop)//'.dat'
         fich2 = 'Spectre_NS'//trim(evop)//'.dat'
         fich3 = 'Spectre_NS'//trim(evop)//'_conv.dat'
      
         if (nid == 0) then
            open (unit=10, file=fich1, form='formatted')
            open (unit=20, file=fich2, form='formatted')
            open (unit=30, file=fich3, form='formatted')
         end if
      
      ! outpost full spectrum (including spurious modes)
         do i = 1, k_dim
      
            qx(:, i) = Q(i)%vx(:)
            qy(:, i) = Q(i)%vy(:)
            if (if3D) qz(:, i) = Q(i)%vz(:)
            if (ifpo) qp(:, i) = Q(i)%pr(:)
            if (ifto) qt(:, i, 1) = Q(i)%t(:, 1)
            if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) qt(:, i, m) = Q(i)%t(:, m)
            end do
            end if
      
            if (nid == 0) then
      
      !outpost the eigenspectrum of the Hessenberg matrix.
               write (10, "(3E15.7)") real(vals(i)), aimag(vals(i)), residual(i)
      
      !outpost the log-transform spectrum (i.e. eigenspectrum of the linearized Navier-Stokes operator).
               write (20, "(3E15.7)") real(log_transform(vals(i)))/speriod, aimag(log_transform(vals(i)))/speriod, residual(i)
      
            end if ! nid.eq.0
         end do ! i = 1, k_dim
      
         spurious_tol = max(param(21), param(22))
         outp = 0 ! outposted modes counter
         do i = 1, converged ! loop over converged modes
      
            if (outp >= maxmodes) then
      
               if (nid == 0) write (6, *) 'maxmodes reached, skipping converged mode ', i
               cycle ! skip nonconverged modes or if too many modes have been outposted
      
            else ! converged modes
      
      !     ----- Computation of the corresponding eigenmode -----
               fp_cx(:) = matmul(qx(:, 1:k_dim), vecs(:, i))
               fp_cy(:) = matmul(qy(:, 1:k_dim), vecs(:, i))
               if (if3D) fp_cz(:) = matmul(qz(:, 1:k_dim), vecs(:, i))
               if (ifpo) fp_cp(:) = matmul(qp(:, 1:k_dim), vecs(:, i))
               if (ifto) fp_ct(:, 1) = matmul(qt(:, 1, 1:k_dim), vecs(:, i))
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) fp_ct(:, m) = matmul(qt(:, m, 1:k_dim), vecs(:, i))
               end do ! m = 2,ldimt
               end if !ldimt.gt.1
      
      ! normalization to unit-norm (volume integral of FP*conj(FP) = 1.)
               call norm(real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct), alpha_r)
               call norm(aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), aimag(fp_ct), alpha_i)
      
               if (nid == 0) write (6, *) 'Checking eigenvector', i
               if (nid == 0) write (6, *)
               if (nid == 0) write (6, *) '       norm Re/Im:', alpha_r, alpha_i
      
               alpha = alpha_r**2 + alpha_i**2
               beta = 1.0d0/sqrt(alpha)
      
               call norm_grad(real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct), norma_Re)
               call norm_grad(aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), aimag(fp_ct), norma_Im)
               if (nid == 0) write (6, *) '  grad norm Re/Im:', norma_Re, norma_Im
               if (nid == 0) write (6, *)
      
      !    if (norma_Re > 1.1 .or. norma_Im > 1.1) then
      !       if (nid == 0) write (6, *) ' Skipping spurious (non-physical) eigenvector:', i, real(vals(i)), aimag(vals(i))
      !       cycle  ! skip this iteration if all real parts are zero
      !    end if !
      
               if (nid == 0) write (6, *) 'Outposting eigenvector:', i, '/', maxmodes
               if (nid == 0) write (6, *) '  sigma=', real(log_transform(vals(i)))/speriod
               if (nid == 0) write (6, *) '  omega=', aimag(log_transform(vals(i)))/speriod
               if (nid == 0) write (6, *) '      f=', (aimag(log_transform(vals(i)))/speriod)/2.0d0*pi
               outp = outp + 1
      
      !     --> Outpost only the converged part of the log-transformed spectrum.
               if (nid == 0) write (30, "(2E15.7)") real(log_transform(vals(i)))/speriod, aimag(log_transform(vals(i)))/speriod
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
      
            close (10); close (20); close (30)
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
               write (844, fmt4) 'uparam0'//trim(char(i))//'=        ', uparam(i)
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
            write (844, fmt2) 'outp=       ', outp
            close (844)
      
         end if
      
         return
      end subroutine outpost_ks
      
      !-----------------------------------------------------------------------
      
      subroutine select_eigenvalues(selected, converged_eigenvalues, vals, delta, nev, n)
      !     This function selects the eigenvalues to be placed in the upper left corner
      !     during the Schur condensation phase.
      !
      !     INPUTS
      !     ------
      !
      !     vals : n-dimensional complex array.
      !     Array containing the eigenvalues.
      
      !     delta : real
      !     All eigenvalues outside the circle of radius 1-delta will be selected.
      !
      !     nev : integer
      !     Number of desired eigenvalues. At least nev+4 eigenvalues will be selected
      !     to ensure "smooth" convergence of the Krylov-Schur iterations.
      !
      !     n : integer
      !     Total number of eigenvalues.
      !
      !     RETURNS
      !     -------
      !
      !     selected : n-dimensional logical array.
      !     Array indicating which eigenvalue has been selected (.true.).
      !
      !     converged_eigenvalues : integer
      !     Number of selected eigenvalues. converged_eigenvalues >= nev + 4.
      !
      !     Last edit : April 2nd 2020 by JC Loiseau.
      
         implicit none
         include 'SIZE'
         include 'TOTAL'
      !     ----- Input arguments -----
         integer :: nev
         integer :: n
         integer, dimension(n) :: idx
         complex(kind=kind(0.0d0)), dimension(n) :: vals
         real :: delta
      
      !     ----- Output argument -----
         logical, dimension(n) :: selected
         integer :: converged_eigenvalues
      
      !     ----- Miscellaneous -----
         integer :: i
      
      !     --> Sort eigenvalues based on their magnitudes (Note : the array vals itself is not sorted).
         do i = 1, n
            idx(i) = i
         end do
         call quicksort2(n, abs(vals), idx)
      
      !     --> Select eigenvalues closer to the unit circle.
         selected = abs(vals) >= (1.0d0 - delta)
      
      !     --> Select at least the nev+4 largest eigenvalues.
         selected(idx(n - (nev + 3):n)) = .true.
         if (aimag(vals(idx(n - (nev + 3)))) == -aimag(vals(idx(n - (nev + 4))))) then
            selected(idx(n - (nev + 4))) = .true.
         end if
      
         converged_eigenvalues = count(selected)
      
         return
      end subroutine select_eigenvalues
      
      !     ------------------------------------------------------------------------------------
      
      subroutine arnoldi_checkpoint(f_xr, f_yr, f_zr, f_pr, f_tr, H, k)
      
      !     This function implements a fairly simple checkpointing procedure in case one
      !     would need to restart the computation (e.g. in case of cluster shutdown).
      !
      !     INPUTS
      !     ------
      !
      !     f_xr, f_yr, f_zr : nek arrays of size lv = lx1*ly1*lz1*lelv
      !     Velocity components of the latest Krylov vector.
      !
      !     f_pr : nek array of size lp = lx2*ly2*lz2*lelt
      !     Pressure field of the latest Krylov vector.
      !
      !     H : k+1 x k real matrix.
      !     Current upper Hessenberg matrix resulting from the k-step Arnoldi factorization.
      !
      !     k : int
      !     Current iteration of the Arnoldi factorization.
      !
      !     Last edit : April 3rd 2020 by JC Loiseau.
      
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv), intent(in) :: f_xr, f_yr, f_zr
         real, dimension(lp), intent(in) :: f_pr
         real, dimension(lt, ldimt), intent(in) :: f_tr
      
         integer, intent(in) :: k
         real, dimension(k + 1, k), intent(in) :: H
      
         integer :: i, j, converged_eigenvalues
         complex(kind=kind(0.0d0)), dimension(k) :: vals
         complex(kind=kind(0.0d0)), dimension(k, k) :: vecs
         real, dimension(k) :: residual
         character(len=80) filename
         complex :: log_transform
      
         if (nid == 0) write (6, *) 'Outposting Krylov vector to'
      
      !     --> Outpost the latest Krylov vector.
      !     if(uparam(1).gt.3)then
         call whereyouwant("KRY", k + 1) ! skipping one due to the initial condition
      !     else
      !     call whereyouwant("KRY", k) ! for restart of newton solver
      !     endif
      
         time = time*k           !order in ParaView
         call outpost2(f_xr, f_yr, f_zr, f_pr, f_tr, nof, "KRY")
      
      !     --> Compute the eigenvalues and eigenvectors of the current Hessenberg matrix.
         call eig(H(1:k, 1:k), vecs, vals, k)
      
      !     --> Compute the residual (assume H results from classical Arnoldi factorization).
         residual = abs(H(k + 1, k)*vecs(k, :))
         converged_eigenvalues = count(residual < eigen_tol)
      
         if (nid == 0) then
      !     --> Outpost the eigenspectrum and residuals of the current Hessenberg matrix.
            write (filename, '(A,A,i4.4,A)') 'Spectre_H', evop, k, '.dat'
            write (6, *) 'Writing Hessenberg matrix eigenspectrum to', filename
      
            open (67, file=trim(filename), status='unknown', form='formatted')
            write (67, '(3E15.7)') (real(vals(i)), aimag(vals(i)), residual(i), i=1, k)
            close (67)
      
      !     --> Outpost the log-transform spectrum (i.e. eigenspectrum of the linearized Navier-Stokes operator).
            write (filename, '(A,A,i4.4,A)') 'Spectre_NS', evop, k, '.dat'
            write (6, *) 'Writing log-transformed eigenspectrum to', filename
      
            open (67, file=trim(filename), status='unknown', form='formatted')
            write (67, '(3E15.7)') (real(log_transform(vals(i)))/(dt*nsteps),
     $   aimag(log_transform(vals(i)))/(dt*nsteps),
     $   residual(i), i = 1, k)
            close (67)
      
      !     --> Outpost the Hessenberg matrix for restarting purposes (if needed).
            write (filename, '(a, a, i4.4)') 'HES', trim(SESSION), k
            write (6, *) 'Writing Hessenberg matrix to', filename
      
            open (67, file=trim(filename), status='unknown', form='formatted')
            write (67, *) ((H(i, j), j=1, k), i=1, k + 1)
            close (67)
      
      !     --> Write to logfile the current number of converged eigenvalues.
            write (6, *) 'converged eigenvalues:', converged_eigenvalues, 'target:', schur_tgt !keep small caps to ease grep
         end if
      
      !   if (schur_tgt > 0 .and. converged_eigenvalues >= schur_tgt) then
      !         !ifres=.true. is required!
      !         if(nid.eq.0)write(6,*) 'Target reached! exporting and stopping'
      !         call nek_end
      !   endif
      
         return
      end subroutine arnoldi_checkpoint
      
      !     ------------------------------------------------------------------------------------
      function log_transform(x)
         implicit none
         complex, intent(in) :: x
         complex :: log_transform
      
      ! Calculate the natural logarithm of x
         log_transform = log(x)
      ! If x is a real number, return a real value
         if (aimag(x) == 0) log_transform = real(log_transform)
      end function log_transform
