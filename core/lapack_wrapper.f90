      !-----------------------------------------------------------------------
      
      subroutine schur(A, vecs, vals, n)
      
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
      !     n : float
      !     Number of rows/columns of A.
      !
      !     RETURNS
      !     -------
      !
      !     A : n x n real matrix
      !     Schur decomposition of the input matrix A, in canonical form.
      !
      !     vecs : n x n real matrix
      !     Schur basis asssociated to A.
      !
      !     vals : n-dimensional complex array.
      !     Unsorted eigenvalues of matrix A.
      !
      !
      !     Last edit : April 1st 2020 by JC Loiseau.
      
         implicit none
         character*1 :: jobvs = "V", sort = "S"
         integer :: n, lda, sdim, ldvs, lwork, info
         real, dimension(n, n) :: A, vecs
         real, dimension(n) :: wr, wi
         real, dimension(3*n) :: work
         logical, dimension(n) :: bwork
         complex*16, dimension(n) :: vals
      
         external select_eigvals
      
      !     --> Perform the Schur decomposition.
         lda = max(1, n)
         ldvs = max(1, n)
         lwork = max(1, 3*n)
      
         call dgees(jobvs, sort, select_eigvals, n, A, lda, sdim, wr, wi, vecs, ldvs, work, lwork, bwork, info)
      
      !     --> Eigenvalues.
         vals = wr*(1.0d0, 0.0d0) + wi*(0.0d0, 1.0d0)
      
         return
      end subroutine schur
      
      !-----------------------------------------------------------------------
      
      subroutine ordschur(T, Q, selected, n)
      
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
      
         implicit none
         character*1 :: job = "N", compq = "V"
         integer :: info, ldq, ldt, liwork, lwork, m, n
         double precision :: s, sep
         logical, dimension(n) :: selected
         real, dimension(n, n) :: T, Q
         real, dimension(n) :: work, wr, wi
         integer, dimension(1) :: iwork
      
      !     --> Order the Schur decomposition.
         ldt = max(1, n)
         ldq = n
         lwork = max(1, n)
         liwork = 1
      
         call dtrsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)
      
         return
      end subroutine ordschur
      
      !-----------------------------------------------------------------------
      subroutine eig(A, vecs, vals, n)
      
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
      !
      !     Last edit : April 1st 2020 by JC Loiseau.
      
         implicit none
         character*1 :: jobvl = "N", jobvr = "V"
         integer :: n, lwork, info, lda, ldvl, ldvr
         real, dimension(n, n) :: A, A_tilde, vr
         real, dimension(1, n) :: vl
         real, dimension(4*n) :: work
         real, dimension(n) :: wr, wi
         complex*16, dimension(n, n) :: vecs
         complex*16, dimension(n) :: vals
         integer :: i
         integer, dimension(n) :: idx
      
      !     --> Compute the eigendecomposition of A.
         lda = n
         ldvl = 1
         ldvr = n
         lwork = 4*n
         A_tilde = A
      
         call dgeev(jobvl, jobvr, n, A_tilde, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
      
      !     --> Transform from real to complex arithmetic.
         vals = wr*(1.0d0, 0.0d0) + wi*(0.0d0, 1.0d0)
         vecs = vr*(1.0d0, 0.0d0)
      
         do i = 1, n - 1
         if (wi(i) > 0) then
            vecs(:, i) = vr(:, i)*(1.0d0, 0.0d0) + vr(:, i + 1)*(0.0d0, 1.0d0)
            vecs(:, i + 1) = vr(:, i)*(1.0d0, 0.0d0) - vr(:, i + 1)*(0.0d0, 1.0d0)
         else if (wi(i) == 0) then
            vecs(:, i) = vr(:, i)*(1.0d0, 0.0d0)
         end if
         end do
      
      !     --> Sort the eigenvalues and eigenvectors by decreasing magnitudes.
         call sort_eigendecomp(vals, vecs, n)
      
         return
      end subroutine eig
      
      !-----------------------------------------------------------------------
      
      subroutine sort_eigendecomp(vals, vecs, n)
      
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
      
         implicit none
         integer :: n
         complex*16, dimension(n) :: vals
         complex*16, dimension(n, n) :: vecs
         real, dimension(n) :: norm
         real :: temp_real
         complex*16 :: temp_complex
         complex*16, dimension(n) :: temp_n
         integer :: k, l
      
      !     ----- Sorting the eigenvalues according to their norm -----
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
         return
      end subroutine sort_eigendecomp
      
      !     -------------------------------------------------------------------
      
      function select_eigvals(wr, wi)
         implicit none
      
      !     ----- Miscellaneous declarations     -----
         logical :: select_eigvals
         real :: wr, wi
      
      !     --> Select eigenvalues based on its magnitude.
         select_eigvals = .false.
         if (sqrt(wr**2 + wi**2) > 0.9) select_eigvals = .true.
      
         return
      end function select_eigvals
      
      !     -------------------------------------------------------------------
      
      subroutine lstsq(A, b, x, m, n)
      
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
      
         implicit none
         character*1 :: trans = "N"
         integer :: m, n, nrhs, lda, ldb, lwork, info
         real, dimension(m, n) :: A, A_tilde
         real, dimension(m) :: b, b_tilde
         real, dimension(n) :: x
         real, dimension(2*m*n) :: work
      
      !     --> Solve the least-squares problem min || Ax - b ||_2.
         nrhs = 1
         lda = m
         ldb = m
         lwork = 2*m*n
         A_tilde = A
         b_tilde = b
      
         call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)
      !write (6, *) "Least-Squares solver :", info
      !if (info /= 0) then         ! some error occured
      !write (6, *) "Least-Squares solver UNsuccessful. CHECK!"
      !!http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga225c8efde208eaf246882df48e590eac.html
      !end if
      
      !     --> Return solution.
         x = b_tilde(1:n)
      
         return
      end subroutine lstsq
