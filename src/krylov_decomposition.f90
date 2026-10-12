!-----------------------------------------------------------------------
! krylov_decomposition.f90 — Arnoldi factorization and Hessenberg update
!
! Purpose:
!   Implements the k-step Arnoldi factorization of the linearized
!   Navier-Stokes operator with modified Gram-Schmidt orthogonalization
!   and full re-orthogonalization.
!
! Public interface:
!   arnoldi_factorization  — k-step Arnoldi factorization
!   update_hessenberg_matrix — orthonormalize and update Hessenberg
!   arnoldi_checkpoint     — save Krylov vector and eigenspectrum checkpoint
!   ritz_pairs             — Ritz values, vectors and residuals of H(1:k,1:k)
!   log_transform          — complex log for eigenvalue (no imag noise)
!
! Dependencies:
!   krylov_subspace, nekstab_lapack (eig), nekstab_ks_select,
!   nekstab_io (whereyouwant), SIZE, TOTAL
!-----------------------------------------------------------------------

module nekstab_krylov_decomposition
   use nekstab_nek_bridge, only: nekStab_dp, nid, ifres, use_cgs, time, dt, &
                                 nsteps, evop, eigen_tol, schur_tgt, SESSION, &
                                 nof, mstep, NEKSTAB_UNIT_HESS, ldimt, isize
   use nekstab_matvec, only: matvec
   use nekstab_lapack, only: eig
   use nekstab_ks_select, only: ks_converged, ks_leading_converged
   use nekstab_io, only: whereyouwant
   implicit none
   private
   public :: arnoldi_factorization, &
             update_hessenberg_matrix, &
             arnoldi_checkpoint, log_transform, ritz_pairs
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
   !     nev : integer, optional
   !     Krylov-Schur only. Stop at the first step where the nev leading Ritz values
   !     have converged. Absent or <= 0: run to mend.
   !
   !     RETURNS
   !     -------
   !
   !     mlast : integer, optional
   !     Last step done: mend, or less after an early stop or a lucky breakdown.
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
   subroutine arnoldi_factorization(Q, H, mstart, mend, ksize, nev, mlast)

      use nekstab_krylov_subspace, only: krylov_vector, k_zero, k_copy, dnekclock

      !     ----- Miscellaneous -----
      real :: alpha
      integer, intent(in) :: mstart, mend, ksize
      integer, intent(in), optional :: nev
      integer, intent(out), optional :: mlast
      integer :: step

      !     ----- Timer -----
      real :: eetime0, eetime1, telapsed, tmiss, avg_time

      !     ----- Orthogonal residual f = w - (Q,w)*Q -----
      type(krylov_vector) :: f
      !     ----- Krylov basis V for the projection MQ = QH -----
      type(krylov_vector), dimension(ksize + 1), intent(inout) :: Q

      !     ----- Upper Hessenberg matrix -----
      !     Krylov-Schur and GMRES resume with columns 1:mstart-1 already filled;
      !     intent(out) would undefine them. Only the new columns are written.
      real, dimension(ksize + 1, ksize), intent(inout) :: H

      !     ----- Check k_dim -----
      if (ksize == 0) then
         if (nid == 0) write (6, *) 'Krylov base dimension == 0! Increase it.. STOP'
         call nek_end
      end if

      !     --> Initialize arrays.
      call k_zero(f); alpha = 0.0d0

      !     --> Arnoldi factorization.
      eetime0 = dnekclock() ! Start time for entire process
      !     Drive the loop with a local index and copy it into mstep.
      !     mstep is what matvec prints. A completed DO would leave it
      !     one past mend, so the residual matvec after a full subspace
      !     logged "from <k_dim+1>/<k_dim>".
      do step = mstart, mend
         mstep = step

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
               write (6, *) 'ARNOLDI: Early termination at step', mstep
               write (6, *) '         Lucky breakdown - invariant subspace found.'
            end if
            call k_copy(Q(mstep + 1), f) ! Still copy the (near-zero) residual
            if (present(mlast)) mlast = mstep
            return ! Exit early - don't pollute basis with noise
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

         !     --> Stop as soon as the nev leading Ritz values have converged.
         !         The test is the one krylov_schur applies after a full pass,
         !         so k_dim is the largest base, not a fixed cost: a long base
         !         costs only the steps that the leading values need.
         !         Rank 0 decides and broadcasts, so all ranks leave the loop
         !         at the same step even if their copies of H differ in the
         !         last bit.
         if (present(nev)) then
            if (nev > 0) then
               if (leading_converged(H, mstep) >= nev) then
                  if (nid == 0) write (6, "('        ARNOLDI - ',I0,' leading Ritz values converged at step ',I0,'/',I0)") &
                     nev, mstep, mend
                  exit
               end if
            end if
         end if
      end do
      if (present(mlast)) mlast = mstep

   end subroutine arnoldi_factorization

   !-----------------------------------------------------------------------
   ! leading_converged — Converged Ritz values at the top of the spectrum of H(1:k,1:k)
   !-----------------------------------------------------------------------
   integer function leading_converged(H, k) result(nlead)
      integer, intent(in) :: k
      real, dimension(:, :), intent(in) :: H
      ! allocatable: called at every step, and k can reach k_dim = 1000 (16 MB for vecs)
      complex(nekStab_dp), allocatable, dimension(:) :: vals
      complex(nekStab_dp), allocatable, dimension(:, :) :: vecs
      real, allocatable, dimension(:) :: residual

      allocate (vals(k), vecs(k, k), residual(k))
      call ritz_pairs(H(1:k + 1, 1:k), k, vals, vecs, residual)
      nlead = ks_leading_converged(vals, residual, k, eigen_tol)
      call bcast(nlead, isize)

   end function leading_converged

   !-----------------------------------------------------------------------
   ! ritz_pairs — Ritz values, vectors and residuals of H(1:k,1:k)
   !
   !     The residual of the Ritz pair (lambda_i, V y_i) is |H(k+1,k) y_i(k)|.
   !     It equals ||M V y_i - lambda_i V y_i|| when M V = V H + f e_k^T holds
   !     for every column. eig sorts by decreasing modulus.
   !-----------------------------------------------------------------------
   subroutine ritz_pairs(H, k, vals, vecs, residual)
      integer, intent(in) :: k
      real, dimension(k + 1, k), intent(in) :: H
      complex(nekStab_dp), dimension(k), intent(out) :: vals
      complex(nekStab_dp), dimension(k, k), intent(out) :: vecs
      real, dimension(k), intent(out) :: residual

      call eig(H(1:k, 1:k), vecs, vals, k)
      residual = abs(H(k + 1, k)*vecs(k, :))

      !     --> Floor the REPORTED residual at float64 relative precision.
      !         The Arnoldi residual r_i = |beta * y_i(k)| is the product of two
      !         independently tiny factors after Krylov-Schur deflation/locking:
      !         beta = H(k+1,k) collapses as the subspace becomes invariant, and
      !         y_i(k) (last component of a converged Ritz vector) -> 0. Their
      !         product can underflow far below machine epsilon (e.g. ~1e-48),
      !         which would write sub-epsilon "accuracy" into the spectrum file -
      !         misleading, since float64 carries only ~16 significant digits.
      !         Floor at eps*max(|lambda|,1) so a fully-locked mode reports
      !         "~1e-16" (converged to machine precision) instead of fake digits.
      !         Classification below is unaffected: the floor eps*max(|lambda|,1)
      !         is many orders of magnitude tighter than eigen_tol*max(|lambda|,1),
      !         so a floored residual still satisfies the convergence criterion
      !         exactly when the raw one did.
      residual = max(residual, epsilon(1.0d0)*max(abs(vals), 1.0d0))

   end subroutine ritz_pairs

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

      use nekstab_krylov_subspace, only: krylov_vector, k_matmul, k_sub2, k_dot, &
                                         k_add2s2, k_normalize
      use nekstab_krylov_inner_products, only: k_project

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
            write (6, *) 'ARNOLDI: Lucky breakdown at step', k
            write (6, *) '         Invariant subspace found.'
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

      use nekstab_krylov_subspace, only: lv, lp, lt

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

      !     --> Compute the eigenvalues, eigenvectors and residuals of the current Hessenberg matrix.
      !         Same residual and test as the run itself, so the checkpoint count
      !         and the stop of the run agree.
      call ritz_pairs(H, k, vals, vecs, residual)
      converged_eigenvalues = count(ks_converged(vals, residual, eigen_tol))

      if (nid == 0) then

         !     --> Outpost the eigenspectrum and residuals of the current Hessenberg matrix.
         write (filename, '(A,A,i4.4,A)') 'Spectre_H', evop, k, '.dat'
         write (6, *) 'Writing Hessenberg matrix eigenspectrum to', filename

         open (NEKSTAB_UNIT_HESS, file=trim(filename), status='unknown', form='formatted')
         write (NEKSTAB_UNIT_HESS, '(3E15.7)') (real(vals(i)), aimag(vals(i)), &
                                                residual(i), i=1, k)
         close (NEKSTAB_UNIT_HESS)

         !     --> Outpost the log-transform spectrum.
         write (filename, '(A,A,i4.4,A)') 'Spectre_NS', evop, k, &
            '.dat'
         write (6, *) 'Writing log-transformed eigenspectrum to', &
            filename

         open (NEKSTAB_UNIT_HESS, file=trim(filename), status='unknown', &
               form='formatted')
         write (NEKSTAB_UNIT_HESS, '(3E15.7)') &
            (real(log_transform(vals(i)))/(dt*nsteps), &
             aimag(log_transform(vals(i)))/(dt*nsteps), &
             residual(i), i=1, k)
         close (NEKSTAB_UNIT_HESS)

         !     --> Outpost the Hessenberg matrix for restarting purposes.
         write (filename, '(a, a, i4.4)') 'HES', trim(SESSION), k
         write (6, *) 'Writing Hessenberg matrix to', filename

         open (NEKSTAB_UNIT_HESS, file=trim(filename), status='unknown', &
               form='formatted')
         write (NEKSTAB_UNIT_HESS, *) ((H(i, j), j=1, k), i=1, k + 1)
         close (NEKSTAB_UNIT_HESS)

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
