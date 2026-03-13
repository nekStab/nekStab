      !-----------------------------------------------------------------------
      ! krylov_decomposition.f90 — Arnoldi factorization and Hessenberg update
      !
      ! Purpose:
      !   Implements the k-step Arnoldi factorization of the linearized
      !   Navier-Stokes operator with modified Gram-Schmidt orthogonalization
      !   and full re-orthogonalization.
      !
      ! Public interface:
      !   arnoldi_factorization, update_hessenberg_matrix,
      !   arnoldi_checkpoint, log_transform
      !
      ! Dependencies:
      !   krylov_subspace, nekstab_lapack (eig),
      !   nekstab_io (whereyouwant), SIZE, TOTAL
      !-----------------------------------------------------------------------

   module nekstab_krylov_decomposition
   use krylov_subspace
   use nekstab_nek_bridge
   use nekstab_matvec
   use nekstab_lapack
   use nekstab_io
   implicit none
   private
   public :: arnoldi_factorization, &
      update_hessenberg_matrix, &
      arnoldi_checkpoint, log_transform
   contains

      !-----------------------------------------------------------------------
      ! arnoldi_factorization — k-step Arnoldi factorization
      !
      !     This function implements the k-step Arnoldi factorization of the linearized
      !     Navier-Stokes operator. The rank k of the Arnoldi factorization is set as a user
      !     parameter in x_SIZE.f (see parameter k_dim).
      !
      !     INPUT
      !     -----
      !
      !     mstart : integer
      !     Index at which to start the Arnoldi factorization. By default, it should be set to 1.
      !     Note that it changes when the Arnoldi factorization is used as part of the Krylov-Schur
      !     algorithm.
      !
      !     mend : integer
      !     Index at which to stop the Arnoldi factorization. By default, it should be set to kdim.
      !     Note that it changes when the Arnoldi factorization is used as part of the GMRES solver
      !
      !     ksize : integer
      !     Size of the Krylov subspace (same as k_dim).
      !
      !     RETURNS
      !     -------
      !
      !     qx, qy, qz : nek arrays of size (lx1*ly1*lz1*lelt, ksize).
      !     Arrays containing the various Krylov vectors associated to each velocity component.
      !
      !     qp : nek arrays of size (lx2*ly2*lz2*lelt, ksize)
      !     Arrays containing the various Krylov vectors associated to the pressure field.
      !
      !     H : k x k real matrix.
      !     Upper Hessenberg matrix resulting from the Arnoldi factorization of the linearized
      !     Navier-Stokes operator.
      !-----------------------------------------------------------------------
   subroutine arnoldi_factorization(Q, H, mstart, mend, ksize)

   use krylov_subspace

      !     ----- Miscellaneous -----
   real :: alpha
   integer, intent(in) :: mstart, mend, ksize

      !     ----- Timer -----
   real :: eetime0, eetime1, telapsed, tmiss, avg_time

      !     ----- Orthogonal residual f = w - (Q,w)*Q -----
   type(krylov_vector) :: f
      !     ----- Krylov basis V for the projection MQ = QH -----
   type(krylov_vector), dimension(ksize + 1), intent(inout) :: Q

      !     ----- Upper Hessenberg matrix -----
   real, dimension(ksize + 1, ksize), intent(out) :: H

      !     ----- Check k_dim -----
   if (ksize == 0) then
   if (nid == 0) write (6, *) 'Krylov base dimension == 0! Increase it.. STOP'
   call nek_end
   end if

      !     --> Initialize arrays.
   call k_zero(f); alpha = 0.0d0

      !     --> Arnoldi factorization.
   eetime0 = dnekclock() ! Start time for entire process
   do mstep = mstart, mend

   if (mstart < mend) then
   if (nid == 0) write (6, "('        ARNOLDI - Starting iteration ',I3,'/',I3)") mstep, mend
   end if

      !     --> Matrix-vector product f = M * v (e.g. calling the linearized Navier-Stokes solver).
   call matvec(f, Q(mstep))

      !     --> Update Hessenberg matrix and compute the orthogonal residual f.
   call update_hessenberg_matrix(H(1:mstep + 1, 1:mstep), f, Q(1:mstep), mstep)

      !     --> Check for lucky breakdown: if H(k+1,k) is near-zero, we've found an
      !         invariant subspace. The eigenvalues of H(1:k,1:k) are exact.
      !         Continuing would amplify round-off noise into spurious Krylov vectors.
   if (H(mstep + 1, mstep) < 1.0d-12) then
   if (nid == 0) then
   write(6,*) 'ARNOLDI: Early termination at step', mstep
   write(6,*) '         Lucky breakdown - invariant subspace found.'
   end if
   call k_copy(Q(mstep + 1), f)  ! Still copy the (near-zero) residual
   return  ! Exit early - don't pollute basis with noise
   end if

      !     --> Add the residual vector as the new Krylov vector.
   call k_copy(Q(mstep + 1), f)

      !     --> Save checkpoint for restarting/run-time analysis. #  --> not in newton !
   if (ifres .and. (mstart < mend)) call arnoldi_checkpoint(f%vx, f%vy, f%vz, f%pr, f%t, H(1:mstep + 1, 1:mstep), mstep)

   if (nid == 0) then
   if (mstart == mend) then
   write (6, "('        ARNOLDI - Finished iteration ',I3,'/',I3)") mstep, mend
   else
   eetime1 = dnekclock()
   telapsed = (eetime1 - eetime0)/3600.0d0
   avg_time = telapsed/mstep ! Average time per iteration
   tmiss = avg_time*(mend - mstep) ! More accurate ETA based on average
   write (6, "('        ARNOLDI - Finished iteration:',I3,'/',I3,' elapsed:',I3,'h',I2,'m / ETA:',I3,'h',I2,'m')") &
      mstep, mend, int(telapsed), ceiling((telapsed - int(telapsed))*60.0d0), &
      int(tmiss), ceiling((tmiss - int(tmiss))*60.0d0)
   end if
   end if
   end do

   end subroutine arnoldi_factorization

      !-----------------------------------------------------------------------
      ! update_hessenberg_matrix — Orthonormalize and update Hessenberg entries
      !
      !     This function orthonormalizes the latest Krylov vector f w.r.t. all of the
      !     previous ones and updates the entries of the Hessenberg matrix accordingly.
      !
      !     INPUTS
      !     ------
      !
      !     k : int
      !     Current step of the Arnoldi factorization.
      !
      !     f_xr, f_yr, f_zr : nek arrays of size lt = lx1*ly1*lz1*lelt.
      !     Velocity components of the latest Krylov vector.
      !     When returned, it has been orthonormalized w.r.t. to all previous
      !     Krylov vectors.
      !
      !     f_pr : nek array of size lp = lx2*ly2*lz2*lelt.
      !     Pressure component of the latest Krylov vector.
      !
      !     qx, qy, qz : nek arrays of size (lv, k)
      !     Velocity components of the Krylov basis.
      !
      !     qp : nek array of size (lp, k).
      !     Pressure component of the Krylov basis.
      !
      !     H : k x k real matrix.
      !     Upper Hessenberg matrix.
      !
      !     Last edit : April 3rd 2020 by JC Loiseau.
      !-----------------------------------------------------------------------
   subroutine update_hessenberg_matrix(H, f, q, k)

   use krylov_subspace
   use krylov_inner_products

   integer, intent(in) :: k
   real, dimension(k + 1, k), intent(inout) :: H

   type(krylov_vector), dimension(k), intent(in) :: q
   type(krylov_vector), intent(inout) :: f
   type(krylov_vector) :: wrk

   integer :: i
   real :: alpha
   real :: h1(k), h2(k)

   if (use_cgs) then
!        ── CGS2: 2 gop calls total (vs 2k for MGS) ──
!           Pass 1
   call k_project(h1, q, f, k)
   call k_matmul(wrk, q, h1, k)
   call k_sub2(f, wrk)
!           Pass 2 (reorthogonalization)
   call k_project(h2, q, f, k)
   call k_matmul(wrk, q, h2, k)
   call k_sub2(f, wrk)
!           Combine
   do i = 1, k
   H(i, k) = h1(i) + h2(i)
   end do
   else
!        ── MGS with reorthogonalization: 2k gop calls ──
   do i = 1, k
   call k_dot(alpha, f, q(i))
   call k_add2s2(f, q(i), -alpha)
   H(i, k) = alpha
   end do
   do i = 1, k
   call k_dot(alpha, f, q(i))
   call k_add2s2(f, q(i), -alpha)
   H(i, k) = H(i, k) + alpha
   end do
   end if

!     --> Normalise the residual vector.
   call k_normalize(f, alpha)
   H(k + 1, k) = alpha

!     --> Lucky breakdown check.
   if (alpha < 1.0d-12) then
   if (nid == 0) then
   write(6,*) 'ARNOLDI: Lucky breakdown at step', k
   write(6,*) '         Invariant subspace found.'
   end if
   end if

   end subroutine update_hessenberg_matrix

      !-----------------------------------------------------------------------
      ! arnoldi_checkpoint — Save Krylov vector and eigenspectrum checkpoint
      !
      ! Purpose:
      !   Implements checkpointing during the Arnoldi process for restart
      !   capability. Saves the latest Krylov vector, Hessenberg eigenspectrum,
      !   log-transformed spectrum, and the Hessenberg matrix itself.
      !
      ! Arguments:
      !   f_xr, f_yr, f_zr [in] — velocity components of latest Krylov vector
      !   f_pr              [in] — pressure field of latest Krylov vector
      !   f_tr              [in] — temperature/passive scalar fields
      !   H                 [in] — upper Hessenberg matrix (k+1 x k)
      !   k                 [in] — current Arnoldi iteration
      !-----------------------------------------------------------------------
   subroutine arnoldi_checkpoint(f_xr, f_yr, f_zr, f_pr, f_tr, H, k)

   use krylov_subspace

   real, dimension(lv), intent(in) :: f_xr, f_yr, f_zr
   real, dimension(lp), intent(in) :: f_pr
   real, dimension(lt, ldimt), intent(in) :: f_tr

   integer, intent(in) :: k
   real, dimension(k + 1, k), intent(in) :: H

   integer :: i, j, converged_eigenvalues
   complex(nekStab_dp), dimension(k) :: vals
   complex(nekStab_dp), dimension(k, k) :: vecs
   real, dimension(k) :: residual
   character(len=80) filename

   if (nid == 0) write (6, *) 'Outposting Krylov vector to'

      !     --> Outpost the latest Krylov vector.
   call whereyouwant("KRY", k + 1)

   time = time*k !order in ParaView
   call outpost2(f_xr, f_yr, f_zr, f_pr, f_tr, nof, "KRY")

      !     --> Compute the eigenvalues and eigenvectors of the current Hessenberg matrix.
   call eig(H(1:k, 1:k), vecs, vals, k)

      !     --> Compute the residual (assume H results from classical Arnoldi factorization).
   residual = abs(H(k + 1, k)*vecs(k, :))

      ! ------> Enforce minimum value of machine epsilon
   where (residual < epsilon(1.0d0)) residual = epsilon(1.0d0)

   converged_eigenvalues = count(residual < eigen_tol)

   if (nid == 0) then

      !     --> Outpost the eigenspectrum and residuals of the current Hessenberg matrix.
   write (filename, '(A,A,i4.4,A)') 'Spectre_H', evop, k, '.dat'
   write (6, *) 'Writing Hessenberg matrix eigenspectrum to', filename

   open (67, file=trim(filename), status='unknown', form='formatted')
   write (67, '(3E15.7)') (real(vals(i)), aimag(vals(i)), &
      residual(i), i=1, k)
   close (67)

      !     --> Outpost the log-transform spectrum.
   write (filename, '(A,A,i4.4,A)') 'Spectre_NS', evop, k, &
      '.dat'
   write (6, *) 'Writing log-transformed eigenspectrum to', &
      filename

   open (67, file=trim(filename), status='unknown', &
      form='formatted')
   write (67, '(3E15.7)') &
      (real(log_transform(vals(i)))/(dt*nsteps), &
      aimag(log_transform(vals(i)))/(dt*nsteps), &
      residual(i), i = 1, k)
   close (67)

      !     --> Outpost the Hessenberg matrix for restarting purposes.
   write (filename, '(a, a, i4.4)') 'HES', trim(SESSION), k
   write (6, *) 'Writing Hessenberg matrix to', filename

   open (67, file=trim(filename), status='unknown', &
      form='formatted')
   write (67, *) ((H(i, j), j=1, k), i=1, k + 1)
   close (67)

      !     --> Write to logfile the current number of converged eigenvalues.
   write (6, *) 'converged eigenvalues:', &
      converged_eigenvalues, 'target:', schur_tgt

   end if

   end subroutine arnoldi_checkpoint

      !-----------------------------------------------------------------------
      ! log_transform — Complex logarithm for eigenvalue conversion
      !
      ! Purpose:
      !   Computes log(x) for a complex eigenvalue. If the imaginary part
      !   is zero, returns a purely real result (avoids spurious imaginary
      !   component from floating-point noise).
      !-----------------------------------------------------------------------
    function log_transform(x)
    complex(nekStab_dp), intent(in) :: x
    complex(nekStab_dp) :: log_transform
   log_transform = log(x)
   if (aimag(x) == 0) log_transform = real(log_transform)
   end function log_transform

   end module nekstab_krylov_decomposition
