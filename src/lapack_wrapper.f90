!-----------------------------------------------------------------------
! lapack_wrapper.f90 -- Wrappers for LAPACK eigenvalue and least-squares routines
!
! Purpose:
!   Provides simplified Fortran interfaces to LAPACK routines for
!   Schur decomposition, eigenvalue problems, reordering, and
!   least-squares solves used throughout the eigensolver.
!
! Public interface:
!   schur, ordschur, eig, sort_eigendecomp, select_eigvals,
!   lstsq, eig_symmetric, eig_hermitian
!
! Dependencies:
!   LAPACK (dgees, dtrsen, dgeev, dgels, dsyev, zheev)
!-----------------------------------------------------------------------

module nekstab_lapack
    use nekstab_nek_bridge
    implicit none
    private
    public :: schur, ordschur, eig, sort_eigendecomp, &
       select_eigvals, lstsq, eig_symmetric, &
       eig_hermitian
contains

!-----------------------------------------------------------------------
! schur -- Schur decomposition of a general matrix
!
!     This function computes the Schur decomposition of a general matrix A.
!     Both the eigenvalues and the corresponding Schur basis are returned.
!     Note that the matrix A is overwritten with its Schur factorization.
!
!     INPUTS
!     ------
!
!     A : n x n real matrix
!     Matrix to be factorized.
!
!     n : integer
!     Number of rows/columns of A.
!
!     RETURNS
!     -------
!
!     A : n x n real matrix
!     Schur decomposition of the input matrix A, in canonical form.
!
!     vecs : n x n real matrix
!     Schur basis associated to A.
!
!     vals : n-dimensional complex array.
!     Unsorted eigenvalues of matrix A.
!-----------------------------------------------------------------------
subroutine schur(A, vecs, vals, n)
   character(len=1) :: jobvs = "V", sort = "S"
   integer, intent(in) :: n
   integer :: lda, sdim, ldvs, lwork, info
   real, dimension(n, n), intent(inout) :: A
   real, dimension(n, n), intent(out) :: vecs
   real, dimension(n) :: wr, wi
   real, dimension(3*n) :: work
   logical, dimension(n) :: bwork
   complex(nekStab_dp), dimension(n), intent(out) :: vals

   ! --> Perform the Schur decomposition.
   lda = max(1, n)
   ldvs = max(1, n)
   lwork = max(1, 3*n)

   call dgees(jobvs, sort, select_eigvals, n, A, lda, sdim, wr, wi, vecs, ldvs, work, lwork, bwork, info)

   ! --> Eigenvalues.
   vals = wr*(1.0d0, 0.0d0) + wi*(0.0d0, 1.0d0)

end subroutine schur

!-----------------------------------------------------------------------
! ordschur -- Reorder Schur decomposition
!
!     Given a matrix T in canonical Schur form and the corresponding Schur basis Q,
!     this function reorder the Schur factorization and returns the reorder Schur
!     matrix and corresponding Schur vectors such that the selected eigenvalues are
!     in the upper-left block of the matrix. Note that, after completion, both T
!     and Q are overwritten by the reordered Schur matrix and vectors.
!
!     INPUTS
!     ------
!
!     T : n x n real matrix
!     Matrix in canonical Schur form to be reordered.
!
!     Q : n x n real matrix
!     Matrix of Schur vectors to be reordered.
!
!     selected : logical n-dimensional array.
!     Logical array indicating which eigenvalues need to be moved to the upper left block.
!
!     n : integer
!     Number of rows/columns of T and Q.
!
!     RETURNS
!     -------
!
!     T : n x n real matrix
!     Reordered Schur matrix.
!
!     Q : n x n real matrix.
!     Reordered Schur vectors.
!
!     Last edit : April 1st 2020 by JC Loiseau.
!-----------------------------------------------------------------------
subroutine ordschur(T, Q, selected, n)
   character(len=1) :: job = "N", compq = "V"
   integer, intent(in) :: n
   integer :: info, ldq, ldt, liwork, lwork, m
    real(nekStab_dp) :: s, sep
   logical, dimension(n), intent(inout) :: selected
   real, dimension(n, n), intent(inout) :: T, Q
   real, dimension(n) :: work, wr, wi
   integer, dimension(1) :: iwork

   ! --> Order the Schur decomposition.
   ldt = max(1, n)
   ldq = n
   lwork = max(1, n)
   liwork = 1

   call dtrsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)

end subroutine ordschur

!-----------------------------------------------------------------------
! eig -- Eigendecomposition of a general matrix
!
!     This function computes the eigendecomposition of a general matrix A.
!     Both the eigenvalues and the right eigenvectors are returned.
!
!     INPUTS
!     ------
!
!     A : n x n real matrix.
!     Matrix to be eigendecomposed.
!
!     n : integer
!     Number of rows/columns of A.
!
!     RETURNS
!     -------
!
!     vecs : n x n complex matrix.
!     Matrix of eigenvectors.
!
!     vals : n-dimensional complex array.
!     Array containing the eigenvalues.
!-----------------------------------------------------------------------
subroutine eig(A, vecs, vals, n)
   character(len=1) :: jobvl = "N", jobvr = "V"
   integer, intent(in) :: n
   integer :: lwork, info, lda, ldvl, ldvr
   real, dimension(n, n), intent(in) :: A
   real, dimension(n, n) :: A_tilde, vr
   real, dimension(1, n) :: vl
   real, dimension(4*n) :: work
   real, dimension(n) :: wr, wi
   complex(nekStab_dp), dimension(n, n), intent(out) :: vecs
   complex(nekStab_dp), dimension(n), intent(out) :: vals
   integer :: i

   ! --> Compute the eigendecomposition of A.
   lda = n
   ldvl = 1
   ldvr = n
   lwork = 4*n
   A_tilde = A

   call dgeev(jobvl, jobvr, n, A_tilde, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

   ! --> Transform from real to complex arithmetic.
   vals = wr*(1.0d0, 0.0d0) + wi*(0.0d0, 1.0d0)
   vecs = vr*(1.0d0, 0.0d0)

   do i = 1, n - 1 ! Process pairs up to n-1 to avoid buffer overflow
      if (wi(i) > 0) then
         vecs(:, i) = vr(:, i)*(1.0d0, 0.0d0) + vr(:, i + 1)*(0.0d0, 1.0d0)
         vecs(:, i + 1) = vr(:, i)*(1.0d0, 0.0d0) - vr(:, i + 1)*(0.0d0, 1.0d0)
      !else if (wi(i) == 0) then ! redundant code since it's handled by the initialization
      !   vecs(:, i) = vr(:, i)*(1.0d0, 0.0d0)
      end if
   end do

   ! --> Sort the eigenvalues and eigenvectors by decreasing magnitudes.
   call sort_eigendecomp(vals, vecs, n)

end subroutine eig

!-----------------------------------------------------------------------
! sort_eigendecomp -- Sort eigenvalues by decreasing magnitude
!
!     This function sorts the eigenvalues in decreasing magnitude using a very
!     naive sorting algorithm.
!
!     INPUTS/OUTPUTS
!     --------------
!
!     vals : n-dimensional complex array.
!     Array containing the eigenvalues to be sorted as input.
!     It is overwritten with the ordered eigenvalues as output.
!
!     vecs : n x n complex matrix.
!     Matrix of corresponding eigenvectors. It is also overwritten with
!     the reordered eigenvectors as output.
!
!     Last edit : April 2nd by JC Loiseau
!-----------------------------------------------------------------------
subroutine sort_eigendecomp(vals, vecs, n)
   integer, intent(in) :: n
   complex(nekStab_dp), dimension(n), intent(inout) :: vals
   complex(nekStab_dp), dimension(n, n), intent(inout) :: vecs
   real, dimension(n) :: norm
   real :: temp_real
   complex(nekStab_dp) :: temp_complex
   complex(nekStab_dp), dimension(n) :: temp_n
   integer :: k, l

   ! ----- Sorting the eigenvalues according to their norm -----
   temp_n = (0.0d0, 0.0d0)
   norm = sqrt(real(vals)**2 + aimag(vals)**2)
   do k = 1, n - 1
      do l = k + 1, n
         if (norm(k) < norm(l)) then
            temp_real = norm(k)
            temp_complex = vals(k)
            temp_n = vecs(:, k)
            norm(k) = norm(l)
            norm(l) = temp_real
            vals(k) = vals(l)
            vals(l) = temp_complex
            vecs(:, k) = vecs(:, l)
            vecs(:, l) = temp_n
         end if
      end do
   end do

end subroutine sort_eigendecomp

!-----------------------------------------------------------------------
! select_eigvals -- Eigenvalue selection function for dgees
!-----------------------------------------------------------------------
!  pure is valid as a LAPACK dgees callback: dgees declares SELECT as
!  EXTERNAL (no explicit interface), and a pure function satisfies a
!  non-pure contract.  Do NOT mark elemental — unusual for callbacks.
!-----------------------------------------------------------------------
pure function select_eigvals(wr, wi)

   ! ----- Miscellaneous declarations -----
   logical :: select_eigvals
   real, intent(in) :: wr, wi
   real, parameter :: EIGVAL_MAG_THRESHOLD = 0.9d0  ! Select eigenvalues near unit circle

   ! --> Select eigenvalues based on its magnitude.

   select_eigvals = (sqrt(wr**2 + wi**2) > EIGVAL_MAG_THRESHOLD)

end function select_eigvals

!-----------------------------------------------------------------------
! lstsq -- Linear least-squares solver
!
!     Wrapper for the LAPACK linear least-squares solver. Given the matrix A
!     and right-hand side vector b, it solves for x that minimizes
!
!     min || Ax - b ||_2
!
!     INPUTS
!     ------
!
!     A : m x n real matrix.
!
!     b : m x 1 real vector.
!
!     m, n : integers.
!
!
!     RETURNS
!     -------
!
!     x : n x 1 real vector.
!
!     Last edit : March 22nd 2021 by JC Loiseau.
!-----------------------------------------------------------------------
subroutine lstsq(A, b, x, m, n)
   character(len=1) :: trans = "N"
   integer, intent(in) :: m, n
   integer :: nrhs, lda, ldb, lwork, info
   real, dimension(m, n), intent(in) :: A
   real, dimension(m, n) :: A_tilde
   real, dimension(m), intent(in) :: b
   real, dimension(m) :: b_tilde
   real, dimension(n), intent(out) :: x
   real, dimension(2*m*n) :: work

   ! --> Solve the least-squares problem min || Ax - b ||_2.
   nrhs = 1
   lda = m
   ldb = m
   lwork = 2*m*n
   A_tilde = A
   b_tilde = b

   call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)

   ! --> Return solution.
   x = b_tilde(1:n)

end subroutine lstsq

!-----------------------------------------------------------------------
! eig_symmetric -- Symmetric eigenvalue problem (DSYEV wrapper)
!     Returns eigenvalues in DESCENDING order (largest first).
!-----------------------------------------------------------------------
subroutine eig_symmetric(A, eigvals, eigvecs, n)
   integer, intent(in) :: n
   real, intent(in) :: A(n, n)
   real, intent(out) :: eigvals(n), eigvecs(n, n)

   real, allocatable :: work(:), Acopy(:,:)
   real :: tmp_val
   real, allocatable :: tmp_vec(:)
   integer :: lwork, info, i, j

   ! Copy input (DSYEV overwrites)
   allocate(Acopy(n,n), tmp_vec(n))
   Acopy = A

   ! Query optimal workspace
   allocate(work(1))
   call dsyev('V', 'U', n, Acopy, n, eigvals, work, -1, info)
   lwork = int(work(1))
   deallocate(work)
   allocate(work(lwork))

   ! Compute eigenvalues/vectors (ascending order from LAPACK)
   call dsyev('V', 'U', n, Acopy, n, eigvals, work, lwork, info)

   if (info /= 0) then
      write(6,*) 'ERROR: DSYEV failed with info =', info
      call nek_end
   end if

   eigvecs = Acopy

   ! Reverse to descending order
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
! eig_hermitian -- Complex Hermitian eigenvalue problem (ZHEEV wrapper)
!     Returns eigenvalues in DESCENDING order (largest first).
!-----------------------------------------------------------------------
subroutine eig_hermitian(A, eigvals, eigvecs, n)
   integer, intent(in) :: n
   complex(nekStab_dp), intent(in) :: A(n, n)
   real, intent(out) :: eigvals(n)
   complex(nekStab_dp), intent(out) :: eigvecs(n, n)

   complex(nekStab_dp), allocatable :: work(:), Acopy(:,:)
   real, allocatable :: rwork(:)
   complex(nekStab_dp), allocatable :: tmp_vec(:)
   real :: tmp_val
   integer :: lwork, info, i, j

   ! Copy input (ZHEEV overwrites)
   allocate(Acopy(n,n), rwork(max(1, 3*n-2)), tmp_vec(n))
   Acopy = A

    ! Query optimal workspace
    allocate(work(1))
    call zheev('V', 'U', n, Acopy, n, eigvals, &
       work, -1, rwork, info)
    lwork = int(real(work(1), nekStab_dp))
    deallocate(work)
    allocate(work(max(1, lwork)))

   ! Compute eigenvalues/vectors (ascending order from LAPACK)
   call zheev('V', 'U', n, Acopy, n, eigvals, &
      work, lwork, rwork, info)

   if (info /= 0) then
      write(6,*) 'ERROR: ZHEEV failed with info =', info
      call nek_end
   end if

   eigvecs = Acopy

   ! Reverse to descending order
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

end module nekstab_lapack
