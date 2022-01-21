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
      integer  :: n, lda, sdim, ldvs, lwork, info
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
      vals = wr*(1.0D0, 0.0D0) + wi*(0.0D0, 1.0D0)

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
      vals = wr*(1.0D0, 0.0D0) + wi*(0.0D0, 1.0D0)
      vecs = vr*(1.0D0, 0.0D0)

      do i = 1, n-1
         if (wi(i) .gt. 0) then
            vecs(:, i) = vr(:, i)*(1.0D0, 0.0D0) + vr(:, i+1)*(0.0D0, 1.0D0)
            vecs(:, i+1) = vr(:, i)*(1.0D0, 0.0D0) - vr(:, i+1)*(0.0D0, 1.0D0)
         else if (wi(i) .eq. 0) then
            vecs(:, i) = vr(:, i)*(1.0D0, 0.0D0)
         endif
      enddo

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
      integer                    :: n
      complex*16, dimension(n)   :: vals
      complex*16, dimension(n,n) :: vecs
      real, dimension(n)         :: norm
      real                       :: temp_real
      complex*16                 :: temp_complex
      complex*16, dimension(n)   :: temp_n
      integer                    :: k, l

!     ----- Sorting the eigenvalues according to their norm -----
      temp_n   = (0.0D0,0.0D0)
      norm = sqrt( real(vals)**2 + aimag(vals)**2 )
      do k = 1,n-1
         do l = k+1,n
            if (norm(k).LT.norm(l)) then
               temp_real    = norm(k)
               temp_complex = vals(k)
               temp_n       = vecs(:,k)
               norm(k) = norm(l)
               norm(l) = temp_real
               vals(k)  = vals(l)
               vals(l)  = temp_complex
               vecs(:,k) = vecs(:,l)
               vecs(:,l) = temp_n
            endif
         enddo
      enddo
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
      if(sqrt(wr**2 + wi**2) .GT. 0.9) select_eigvals=.true.

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

      write(6,*) "Least-Squares solver :", info

      if(info.ne.0)then         ! some error occured
         write(6,*) "Least-Squares solver UNsuccessful. CHECK!"
!     http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga225c8efde208eaf246882df48e590eac.html
      endif
!     --> Return solution.
      x = b_tilde(1:n)

      return
      end subroutine lstsq

!!!!!!!!!additional LAPACK routines added by Ricardo to remove external dependencies  
!!!!!!!http://www.netlib.org/lapack/#_lapack_version_3_10_0_2
!!!!!!!LAPACK, version 3.10.0 
!!!!!!!Updated: Jun 28, 2021

!     The following files are appended here:
!     dgees.f		dlacn2.f	dlasy2.f	dtrsen.f	dtrtrs.f
!     dgels.f		dlaexc.f	dtrexc.f	dtrsyl.f

*     > \brief <b> DGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DGEES + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgees.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgees.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgees.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
*     VS, LDVS, WORK, LWORK, BWORK, INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER          JOBVS, SORT
*     INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
*     LOGICAL            BWORK( * )
*     DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
*     $                   WR( * )
*     ..
*     .. Function Arguments ..
*     LOGICAL            SELECT
*     EXTERNAL           SELECT
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DGEES computes for an N-by-N real nonsymmetric matrix A, the
*     > eigenvalues, the real Schur form T, and, optionally, the matrix of
*     > Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
*     >
*     > Optionally, it also orders the eigenvalues on the diagonal of the
*     > real Schur form so that selected eigenvalues are at the top left.
*     > The leading columns of Z then form an orthonormal basis for the
*     > invariant subspace corresponding to the selected eigenvalues.
*     >
*     > A matrix is in real Schur form if it is upper quasi-triangular with
*     > 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
*     > form
*     >         [  a  b  ]
*     >         [  c  a  ]
*     >
*     > where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] JOBVS
*     > \verbatim
*     >          JOBVS is CHARACTER*1
*     >          = 'N': Schur vectors are not computed;
*     >          = 'V': Schur vectors are computed.
*     > \endverbatim
*     >
*     > \param[in] SORT
*     > \verbatim
*     >          SORT is CHARACTER*1
*     >          Specifies whether or not to order the eigenvalues on the
*     >          diagonal of the Schur form.
*     >          = 'N': Eigenvalues are not ordered;
*     >          = 'S': Eigenvalues are ordered (see SELECT).
*     > \endverbatim
*     >
*     > \param[in] SELECT
*     > \verbatim
*     >          SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments
*     >          SELECT must be declared EXTERNAL in the calling subroutine.
*     >          If SORT = 'S', SELECT is used to select eigenvalues to sort
*     >          to the top left of the Schur form.
*     >          If SORT = 'N', SELECT is not referenced.
*     >          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
*     >          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
*     >          conjugate pair of eigenvalues is selected, then both complex
*     >          eigenvalues are selected.
*     >          Note that a selected complex eigenvalue may no longer
*     >          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
*     >          ordering may change the value of complex eigenvalues
*     >          (especially if the eigenvalue is ill-conditioned); in this
*     >          case INFO is set to N+2 (see INFO below).
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The order of the matrix A. N >= 0.
*     > \endverbatim
*     >
*     > \param[in,out] A
*     > \verbatim
*     >          A is DOUBLE PRECISION array, dimension (LDA,N)
*     >          On entry, the N-by-N matrix A.
*     >          On exit, A has been overwritten by its real Schur form T.
*     > \endverbatim
*     >
*     > \param[in] LDA
*     > \verbatim
*     >          LDA is INTEGER
*     >          The leading dimension of the array A.  LDA >= max(1,N).
*     > \endverbatim
*     >
*     > \param[out] SDIM
*     > \verbatim
*     >          SDIM is INTEGER
*     >          If SORT = 'N', SDIM = 0.
*     >          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*     >                         for which SELECT is true. (Complex conjugate
*     >                         pairs for which SELECT is true for either
*     >                         eigenvalue count as 2.)
*     > \endverbatim
*     >
*     > \param[out] WR
*     > \verbatim
*     >          WR is DOUBLE PRECISION array, dimension (N)
*     > \endverbatim
*     >
*     > \param[out] WI
*     > \verbatim
*     >          WI is DOUBLE PRECISION array, dimension (N)
*     >          WR and WI contain the real and imaginary parts,
*     >          respectively, of the computed eigenvalues in the same order
*     >          that they appear on the diagonal of the output Schur form T.
*     >          Complex conjugate pairs of eigenvalues will appear
*     >          consecutively with the eigenvalue having the positive
*     >          imaginary part first.
*     > \endverbatim
*     >
*     > \param[out] VS
*     > \verbatim
*     >          VS is DOUBLE PRECISION array, dimension (LDVS,N)
*     >          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
*     >          vectors.
*     >          If JOBVS = 'N', VS is not referenced.
*     > \endverbatim
*     >
*     > \param[in] LDVS
*     > \verbatim
*     >          LDVS is INTEGER
*     >          The leading dimension of the array VS.  LDVS >= 1; if
*     >          JOBVS = 'V', LDVS >= N.
*     > \endverbatim
*     >
*     > \param[out] WORK
*     > \verbatim
*     >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*     >          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
*     > \endverbatim
*     >
*     > \param[in] LWORK
*     > \verbatim
*     >          LWORK is INTEGER
*     >          The dimension of the array WORK.  LWORK >= max(1,3*N).
*     >          For good performance, LWORK must generally be larger.
*     >
*     >          If LWORK = -1, then a workspace query is assumed; the routine
*     >          only calculates the optimal size of the WORK array, returns
*     >          this value as the first entry of the WORK array, and no error
*     >          message related to LWORK is issued by XERBLA.
*     > \endverbatim
*     >
*     > \param[out] BWORK
*     > \verbatim
*     >          BWORK is LOGICAL array, dimension (N)
*     >          Not referenced if SORT = 'N'.
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0: successful exit
*     >          < 0: if INFO = -i, the i-th argument had an illegal value.
*     >          > 0: if INFO = i, and i is
*     >             <= N: the QR algorithm failed to compute all the
*     >                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
*     >                   contain those eigenvalues which have converged; if
*     >                   JOBVS = 'V', VS contains the matrix which reduces A
*     >                   to its partially converged Schur form.
*     >             = N+1: the eigenvalues could not be reordered because some
*     >                   eigenvalues were too close to separate (the problem
*     >                   is very ill-conditioned);
*     >             = N+2: after reordering, roundoff changed values of some
*     >                   complex eigenvalues so that leading eigenvalues in
*     >                   the Schur form no longer satisfy SELECT=.TRUE.  This
*     >                   could also be caused by underflow due to scaling.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleGEeigen
*     
*     =====================================================================
      SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
     $     VS, LDVS, WORK, LWORK, BWORK, INFO )
*     
*     -- LAPACK driver routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      CHARACTER          JOBVS, SORT
      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
     $     WR( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CURSL, LASTSL, LQUERY, LST2SL, SCALEA, WANTST,
     $     WANTVS
      INTEGER            HSWORK, I, I1, I2, IBAL, ICOND, IERR, IEVAL,
     $     IHI, ILO, INXT, IP, ITAU, IWRK, MAXWRK, MINWRK
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY,
     $     DLABAD, DLASCL, DORGHR, DSWAP, DTRSEN, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*     
*     Test the input arguments
*     
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) THEN
         INFO = -11
      END IF
*     
*     Compute workspace
*     (Note: Comments in the code beginning "Workspace:" describe the
*     minimal amount of workspace needed at that point in the code,
*     as well as the preferred amount for good performance.
*     NB refers to the optimal block size for the immediately
*     following subroutine, as returned by ILAENV.
*     HSWORK refers to the workspace preferred by DHSEQR, as
*     calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*     the worst case.)
*     
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 3*N
*     
            CALL DHSEQR( 'S', JOBVS, N, 1, N, A, LDA, WR, WI, VS, LDVS,
     $           WORK, -1, IEVAL )
            HSWORK = WORK( 1 )
*     
            IF( .NOT.WANTVS ) THEN
               MAXWRK = MAX( MAXWRK, N + HSWORK )
            ELSE
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1,
     $              'DORGHR', ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, N + HSWORK )
            END IF
         END IF
         WORK( 1 ) = MAXWRK
*     
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
*     
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEES ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*     
*     Quick return if possible
*     
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*     
*     Get machine constants
*     
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*     
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*     
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $     CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*     
*     Permute the matrix to make it more nearly triangular
*     (Workspace: need N)
*     
      IBAL = 1
      CALL DGEBAL( 'P', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
*     
*     Reduce to upper Hessenberg form
*     (Workspace: need 3*N, prefer 2*N+N*NB)
*     
      ITAU = N + IBAL
      IWRK = N + ITAU
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
     $     LWORK-IWRK+1, IERR )
*     
      IF( WANTVS ) THEN
*     
*     Copy Householder vectors to VS
*     
         CALL DLACPY( 'L', N, N, A, LDA, VS, LDVS )
*     
*     Generate orthogonal matrix in VS
*     (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*     
         CALL DORGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ),
     $        LWORK-IWRK+1, IERR )
      END IF
*     
      SDIM = 0
*     
*     Perform QR iteration, accumulating Schur vectors in VS if desired
*     (Workspace: need N+1, prefer N+HSWORK (see comments) )
*     
      IWRK = ITAU
      CALL DHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, WR, WI, VS, LDVS,
     $     WORK( IWRK ), LWORK-IWRK+1, IEVAL )
      IF( IEVAL.GT.0 )
     $     INFO = IEVAL
*     
*     Sort eigenvalues if desired
*     
      IF( WANTST .AND. INFO.EQ.0 ) THEN
         IF( SCALEA ) THEN
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WR, N, IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WI, N, IERR )
         END IF
         DO 10 I = 1, N
            BWORK( I ) = SELECT( WR( I ), WI( I ) )
 10      CONTINUE
*     
*     Reorder eigenvalues and transform Schur vectors
*     (Workspace: none needed)
*     
         CALL DTRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, WR, WI,
     $        SDIM, S, SEP, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1,
     $        ICOND )
         IF( ICOND.GT.0 )
     $        INFO = N + ICOND
      END IF
*     
      IF( WANTVS ) THEN
*     
*     Undo balancing
*     (Workspace: need N)
*     
         CALL DGEBAK( 'P', 'R', N, ILO, IHI, WORK( IBAL ), N, VS, LDVS,
     $        IERR )
      END IF
*     
      IF( SCALEA ) THEN
*     
*     Undo scaling for the Schur form of A
*     
         CALL DLASCL( 'H', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL DCOPY( N, A, LDA+1, WR, 1 )
         IF( CSCALE.EQ.SMLNUM ) THEN
*     
*     If scaling back towards underflow, adjust WI if an
*     offdiagonal element of a 2-by-2 block in the Schur form
*     underflows.
*     
            IF( IEVAL.GT.0 ) THEN
               I1 = IEVAL + 1
               I2 = IHI - 1
               CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI,
     $              MAX( ILO-1, 1 ), IERR )
            ELSE IF( WANTST ) THEN
               I1 = 1
               I2 = N - 1
            ELSE
               I1 = ILO
               I2 = IHI - 1
            END IF
            INXT = I1 - 1
            DO 20 I = I1, I2
               IF( I.LT.INXT )
     $              GO TO 20
               IF( WI( I ).EQ.ZERO ) THEN
                  INXT = I + 1
               ELSE
                  IF( A( I+1, I ).EQ.ZERO ) THEN
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                  ELSE IF( A( I+1, I ).NE.ZERO .AND. A( I, I+1 ).EQ.
     $                    ZERO ) THEN
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                     IF( I.GT.1 )
     $                    CALL DSWAP( I-1, A( 1, I ), 1, A( 1, I+1 ), 1 )
                     IF( N.GT.I+1 )
     $                    CALL DSWAP( N-I-1, A( I, I+2 ), LDA,
     $                    A( I+1, I+2 ), LDA )
                     IF( WANTVS ) THEN
                        CALL DSWAP( N, VS( 1, I ), 1, VS( 1, I+1 ), 1 )
                     END IF
                     A( I, I+1 ) = A( I+1, I )
                     A( I+1, I ) = ZERO
                  END IF
                  INXT = I + 2
               END IF
 20         CONTINUE
         END IF
*     
*     Undo scaling for the imaginary part of the eigenvalues
*     
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-IEVAL, 1,
     $        WI( IEVAL+1 ), MAX( N-IEVAL, 1 ), IERR )
      END IF
*     
      IF( WANTST .AND. INFO.EQ.0 ) THEN
*     
*     Check if reordering successful
*     
         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 30 I = 1, N
            CURSL = SELECT( WR( I ), WI( I ) )
            IF( WI( I ).EQ.ZERO ) THEN
               IF( CURSL )
     $              SDIM = SDIM + 1
               IP = 0
               IF( CURSL .AND. .NOT.LASTSL )
     $              INFO = N + 2
            ELSE
               IF( IP.EQ.1 ) THEN
*     
*     Last eigenvalue of conjugate pair
*     
                  CURSL = CURSL .OR. LASTSL
                  LASTSL = CURSL
                  IF( CURSL )
     $                 SDIM = SDIM + 2
                  IP = -1
                  IF( CURSL .AND. .NOT.LST2SL )
     $                 INFO = N + 2
               ELSE
*     
*     First eigenvalue of conjugate pair
*     
                  IP = 1
               END IF
            END IF
            LST2SL = LASTSL
            LASTSL = CURSL
 30      CONTINUE
      END IF
*     
      WORK( 1 ) = MAXWRK
      RETURN
*     
*     End of DGEES
*     
      END
*     > \brief <b> DGELS solves overdetermined or underdetermined systems for GE matrices</b>
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DGELS + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgels.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgels.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgels.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*     INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER          TRANS
*     INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DGELS solves overdetermined or underdetermined real linear systems
*     > involving an M-by-N matrix A, or its transpose, using a QR or LQ
*     > factorization of A.  It is assumed that A has full rank.
*     >
*     > The following options are provided:
*     >
*     > 1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     >    an overdetermined system, i.e., solve the least squares problem
*     >                 minimize || B - A*X ||.
*     >
*     > 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     >    an underdetermined system A * X = B.
*     >
*     > 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
*     >    an underdetermined system A**T * X = B.
*     >
*     > 4. If TRANS = 'T' and m < n:  find the least squares solution of
*     >    an overdetermined system, i.e., solve the least squares problem
*     >                 minimize || B - A**T * X ||.
*     >
*     > Several right hand side vectors b and solution vectors x can be
*     > handled in a single call; they are stored as the columns of the
*     > M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*     > matrix X.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] TRANS
*     > \verbatim
*     >          TRANS is CHARACTER*1
*     >          = 'N': the linear system involves A;
*     >          = 'T': the linear system involves A**T.
*     > \endverbatim
*     >
*     > \param[in] M
*     > \verbatim
*     >          M is INTEGER
*     >          The number of rows of the matrix A.  M >= 0.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The number of columns of the matrix A.  N >= 0.
*     > \endverbatim
*     >
*     > \param[in] NRHS
*     > \verbatim
*     >          NRHS is INTEGER
*     >          The number of right hand sides, i.e., the number of
*     >          columns of the matrices B and X. NRHS >=0.
*     > \endverbatim
*     >
*     > \param[in,out] A
*     > \verbatim
*     >          A is DOUBLE PRECISION array, dimension (LDA,N)
*     >          On entry, the M-by-N matrix A.
*     >          On exit,
*     >            if M >= N, A is overwritten by details of its QR
*     >                       factorization as returned by DGEQRF;
*     >            if M <  N, A is overwritten by details of its LQ
*     >                       factorization as returned by DGELQF.
*     > \endverbatim
*     >
*     > \param[in] LDA
*     > \verbatim
*     >          LDA is INTEGER
*     >          The leading dimension of the array A.  LDA >= max(1,M).
*     > \endverbatim
*     >
*     > \param[in,out] B
*     > \verbatim
*     >          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*     >          On entry, the matrix B of right hand side vectors, stored
*     >          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*     >          if TRANS = 'T'.
*     >          On exit, if INFO = 0, B is overwritten by the solution
*     >          vectors, stored columnwise:
*     >          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*     >          squares solution vectors; the residual sum of squares for the
*     >          solution in each column is given by the sum of squares of
*     >          elements N+1 to M in that column;
*     >          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*     >          minimum norm solution vectors;
*     >          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
*     >          minimum norm solution vectors;
*     >          if TRANS = 'T' and m < n, rows 1 to M of B contain the
*     >          least squares solution vectors; the residual sum of squares
*     >          for the solution in each column is given by the sum of
*     >          squares of elements M+1 to N in that column.
*     > \endverbatim
*     >
*     > \param[in] LDB
*     > \verbatim
*     >          LDB is INTEGER
*     >          The leading dimension of the array B. LDB >= MAX(1,M,N).
*     > \endverbatim
*     >
*     > \param[out] WORK
*     > \verbatim
*     >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*     >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*     > \endverbatim
*     >
*     > \param[in] LWORK
*     > \verbatim
*     >          LWORK is INTEGER
*     >          The dimension of the array WORK.
*     >          LWORK >= max( 1, MN + max( MN, NRHS ) ).
*     >          For optimal performance,
*     >          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
*     >          where MN = min(M,N) and NB is the optimum block size.
*     >
*     >          If LWORK = -1, then a workspace query is assumed; the routine
*     >          only calculates the optimal size of the WORK array, returns
*     >          this value as the first entry of the WORK array, and no error
*     >          message related to LWORK is issued by XERBLA.
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0:  successful exit
*     >          < 0:  if INFO = -i, the i-th argument had an illegal value
*     >          > 0:  if INFO =  i, the i-th diagonal element of the
*     >                triangular factor of A is zero, so that A does not have
*     >                full rank; the least squares solution could not be
*     >                computed.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleGEsolve
*     
*     =====================================================================
      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
     $     INFO )
*     
*     -- LAPACK driver routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLABAD, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGELQF, DGEQRF, DLASCL, DLASET, DORMLQ, DORMQR,
     $     DTRTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*     
*     Test the input arguments.
*     
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'T' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY )
     $        THEN
         INFO = -10
      END IF
*     
*     Figure out optimal block size
*     
      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN
*     
         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) )
     $        TPSD = .FALSE.
*     
         IF( M.GE.N ) THEN
            NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LN', M, NRHS, N,
     $              -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N,
     $              -1 ) )
            END IF
         ELSE
            NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M,
     $              -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LN', N, NRHS, M,
     $              -1 ) )
            END IF
         END IF
*     
         WSIZE = MAX( 1, MN+MAX( MN, NRHS )*NB )
         WORK( 1 ) = DBLE( WSIZE )
*     
      END IF
*     
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*     
*     Quick return if possible
*     
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RETURN
      END IF
*     
*     Get machine parameters
*     
      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*     
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*     
      ANRM = DLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*     
*     Scale matrix norm up to SMLNUM
*     
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*     
*     Scale matrix norm down to BIGNUM
*     
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*     
*     Matrix all zero. Return zero solution.
*     
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         GO TO 50
      END IF
*     
      BROW = M
      IF( TPSD )
     $     BROW = N
      BNRM = DLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*     
*     Scale matrix norm up to SMLNUM
*     
         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,
     $        INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*     
*     Scale matrix norm down to BIGNUM
*     
         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,
     $        INFO )
         IBSCL = 2
      END IF
*     
      IF( M.GE.N ) THEN
*     
*     compute QR factorization of A
*     
         CALL DGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,
     $        INFO )
*     
*     workspace at least N, optimally N*NB
*     
         IF( .NOT.TPSD ) THEN
*     
*     Least-Squares Problem min || A * X - B ||
*     
*     B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
*     
            CALL DORMQR( 'Left', 'Transpose', M, NRHS, N, A, LDA,
     $           WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $           INFO )
*     
*     workspace at least NRHS, optimally NRHS*NB
*     
*     B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*     
            CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS,
     $           A, LDA, B, LDB, INFO )
*     
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*     
            SCLLEN = N
*     
         ELSE
*     
*     Underdetermined system of equations A**T * X = B
*     
*     B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
*     
            CALL DTRTRS( 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $           A, LDA, B, LDB, INFO )
*     
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*     
*     B(N+1:M,1:NRHS) = ZERO
*     
            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = ZERO
 10            CONTINUE
 20         CONTINUE
*     
*     B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
*     
            CALL DORMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA,
     $           WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $           INFO )
*     
*     workspace at least NRHS, optimally NRHS*NB
*     
            SCLLEN = M
*     
         END IF
*     
      ELSE
*     
*     Compute LQ factorization of A
*     
         CALL DGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,
     $        INFO )
*     
*     workspace at least M, optimally M*NB.
*     
         IF( .NOT.TPSD ) THEN
*     
*     underdetermined system of equations A * X = B
*     
*     B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*     
            CALL DTRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS,
     $           A, LDA, B, LDB, INFO )
*     
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*     
*     B(M+1:N,1:NRHS) = 0
*     
            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = ZERO
 30            CONTINUE
 40         CONTINUE
*     
*     B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
*     
            CALL DORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA,
     $           WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $           INFO )
*     
*     workspace at least NRHS, optimally NRHS*NB
*     
            SCLLEN = N
*     
         ELSE
*     
*     overdetermined system min || A**T * X - B ||
*     
*     B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
*     
            CALL DORMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA,
     $           WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $           INFO )
*     
*     workspace at least NRHS, optimally NRHS*NB
*     
*     B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
*     
            CALL DTRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS,
     $           A, LDA, B, LDB, INFO )
*     
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*     
            SCLLEN = M
*     
         END IF
*     
      END IF
*     
*     Undo scaling
*     
      IF( IASCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,
     $        INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,
     $        INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $        INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $        INFO )
      END IF
*     
 50   CONTINUE
      WORK( 1 ) = DBLE( WSIZE )
*     
      RETURN
*     
*     End of DGELS
*     
      END
*     > \brief \b DTRSEN
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DTRSEN + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrsen.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrsen.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrsen.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI,
*     M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER          COMPQ, JOB
*     INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
*     DOUBLE PRECISION   S, SEP
*     ..
*     .. Array Arguments ..
*     LOGICAL            SELECT( * )
*     INTEGER            IWORK( * )
*     DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ),
*     $                   WR( * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DTRSEN reorders the real Schur factorization of a real matrix
*     > A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
*     > the leading diagonal blocks of the upper quasi-triangular matrix T,
*     > and the leading columns of Q form an orthonormal basis of the
*     > corresponding right invariant subspace.
*     >
*     > Optionally the routine computes the reciprocal condition numbers of
*     > the cluster of eigenvalues and/or the invariant subspace.
*     >
*     > T must be in Schur canonical form (as returned by DHSEQR), that is,
*     > block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*     > 2-by-2 diagonal block has its diagonal elements equal and its
*     > off-diagonal elements of opposite sign.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] JOB
*     > \verbatim
*     >          JOB is CHARACTER*1
*     >          Specifies whether condition numbers are required for the
*     >          cluster of eigenvalues (S) or the invariant subspace (SEP):
*     >          = 'N': none;
*     >          = 'E': for eigenvalues only (S);
*     >          = 'V': for invariant subspace only (SEP);
*     >          = 'B': for both eigenvalues and invariant subspace (S and
*     >                 SEP).
*     > \endverbatim
*     >
*     > \param[in] COMPQ
*     > \verbatim
*     >          COMPQ is CHARACTER*1
*     >          = 'V': update the matrix Q of Schur vectors;
*     >          = 'N': do not update Q.
*     > \endverbatim
*     >
*     > \param[in] SELECT
*     > \verbatim
*     >          SELECT is LOGICAL array, dimension (N)
*     >          SELECT specifies the eigenvalues in the selected cluster. To
*     >          select a real eigenvalue w(j), SELECT(j) must be set to
*     >          .TRUE.. To select a complex conjugate pair of eigenvalues
*     >          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*     >          either SELECT(j) or SELECT(j+1) or both must be set to
*     >          .TRUE.; a complex conjugate pair of eigenvalues must be
*     >          either both included in the cluster or both excluded.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The order of the matrix T. N >= 0.
*     > \endverbatim
*     >
*     > \param[in,out] T
*     > \verbatim
*     >          T is DOUBLE PRECISION array, dimension (LDT,N)
*     >          On entry, the upper quasi-triangular matrix T, in Schur
*     >          canonical form.
*     >          On exit, T is overwritten by the reordered matrix T, again in
*     >          Schur canonical form, with the selected eigenvalues in the
*     >          leading diagonal blocks.
*     > \endverbatim
*     >
*     > \param[in] LDT
*     > \verbatim
*     >          LDT is INTEGER
*     >          The leading dimension of the array T. LDT >= max(1,N).
*     > \endverbatim
*     >
*     > \param[in,out] Q
*     > \verbatim
*     >          Q is DOUBLE PRECISION array, dimension (LDQ,N)
*     >          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*     >          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*     >          orthogonal transformation matrix which reorders T; the
*     >          leading M columns of Q form an orthonormal basis for the
*     >          specified invariant subspace.
*     >          If COMPQ = 'N', Q is not referenced.
*     > \endverbatim
*     >
*     > \param[in] LDQ
*     > \verbatim
*     >          LDQ is INTEGER
*     >          The leading dimension of the array Q.
*     >          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
*     > \endverbatim
*     >
*     > \param[out] WR
*     > \verbatim
*     >          WR is DOUBLE PRECISION array, dimension (N)
*     > \endverbatim
*     > \param[out] WI
*     > \verbatim
*     >          WI is DOUBLE PRECISION array, dimension (N)
*     >
*     >          The real and imaginary parts, respectively, of the reordered
*     >          eigenvalues of T. The eigenvalues are stored in the same
*     >          order as on the diagonal of T, with WR(i) = T(i,i) and, if
*     >          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
*     >          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
*     >          sufficiently ill-conditioned, then its value may differ
*     >          significantly from its value before reordering.
*     > \endverbatim
*     >
*     > \param[out] M
*     > \verbatim
*     >          M is INTEGER
*     >          The dimension of the specified invariant subspace.
*     >          0 < = M <= N.
*     > \endverbatim
*     >
*     > \param[out] S
*     > \verbatim
*     >          S is DOUBLE PRECISION
*     >          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
*     >          condition number for the selected cluster of eigenvalues.
*     >          S cannot underestimate the true reciprocal condition number
*     >          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
*     >          If JOB = 'N' or 'V', S is not referenced.
*     > \endverbatim
*     >
*     > \param[out] SEP
*     > \verbatim
*     >          SEP is DOUBLE PRECISION
*     >          If JOB = 'V' or 'B', SEP is the estimated reciprocal
*     >          condition number of the specified invariant subspace. If
*     >          M = 0 or N, SEP = norm(T).
*     >          If JOB = 'N' or 'E', SEP is not referenced.
*     > \endverbatim
*     >
*     > \param[out] WORK
*     > \verbatim
*     >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*     >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*     > \endverbatim
*     >
*     > \param[in] LWORK
*     > \verbatim
*     >          LWORK is INTEGER
*     >          The dimension of the array WORK.
*     >          If JOB = 'N', LWORK >= max(1,N);
*     >          if JOB = 'E', LWORK >= max(1,M*(N-M));
*     >          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
*     >
*     >          If LWORK = -1, then a workspace query is assumed; the routine
*     >          only calculates the optimal size of the WORK array, returns
*     >          this value as the first entry of the WORK array, and no error
*     >          message related to LWORK is issued by XERBLA.
*     > \endverbatim
*     >
*     > \param[out] IWORK
*     > \verbatim
*     >          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*     >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*     > \endverbatim
*     >
*     > \param[in] LIWORK
*     > \verbatim
*     >          LIWORK is INTEGER
*     >          The dimension of the array IWORK.
*     >          If JOB = 'N' or 'E', LIWORK >= 1;
*     >          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).
*     >
*     >          If LIWORK = -1, then a workspace query is assumed; the
*     >          routine only calculates the optimal size of the IWORK array,
*     >          returns this value as the first entry of the IWORK array, and
*     >          no error message related to LIWORK is issued by XERBLA.
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0: successful exit
*     >          < 0: if INFO = -i, the i-th argument had an illegal value
*     >          = 1: reordering of T failed because some eigenvalues are too
*     >               close to separate (the problem is very ill-conditioned);
*     >               T may have been partially reordered, and WR and WI
*     >               contain the eigenvalues in the same order as in T; S and
*     >               SEP (if requested) are set to zero.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleOTHERcomputational
*     
*     > \par Further Details:
*     =====================
*     >
*     > \verbatim
*     >
*     >  DTRSEN first collects the selected eigenvalues by computing an
*     >  orthogonal transformation Z to move them to the top left corner of T.
*     >  In other words, the selected eigenvalues are the eigenvalues of T11
*     >  in:
*     >
*     >          Z**T * T * Z = ( T11 T12 ) n1
*     >                         (  0  T22 ) n2
*     >                            n1  n2
*     >
*     >  where N = n1+n2 and Z**T means the transpose of Z. The first n1 columns
*     >  of Z span the specified invariant subspace of T.
*     >
*     >  If T has been obtained from the real Schur factorization of a matrix
*     >  A = Q*T*Q**T, then the reordered real Schur factorization of A is given
*     >  by A = (Q*Z)*(Z**T*T*Z)*(Q*Z)**T, and the first n1 columns of Q*Z span
*     >  the corresponding invariant subspace of A.
*     >
*     >  The reciprocal condition number of the average of the eigenvalues of
*     >  T11 may be returned in S. S lies between 0 (very badly conditioned)
*     >  and 1 (very well conditioned). It is computed as follows. First we
*     >  compute R so that
*     >
*     >                         P = ( I  R ) n1
*     >                             ( 0  0 ) n2
*     >                               n1 n2
*     >
*     >  is the projector on the invariant subspace associated with T11.
*     >  R is the solution of the Sylvester equation:
*     >
*     >                        T11*R - R*T22 = T12.
*     >
*     >  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
*     >  the two-norm of M. Then S is computed as the lower bound
*     >
*     >                      (1 + F-norm(R)**2)**(-1/2)
*     >
*     >  on the reciprocal of 2-norm(P), the true reciprocal condition number.
*     >  S cannot underestimate 1 / 2-norm(P) by more than a factor of
*     >  sqrt(N).
*     >
*     >  An approximate error bound for the computed average of the
*     >  eigenvalues of T11 is
*     >
*     >                         EPS * norm(T) / S
*     >
*     >  where EPS is the machine precision.
*     >
*     >  The reciprocal condition number of the right invariant subspace
*     >  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
*     >  SEP is defined as the separation of T11 and T22:
*     >
*     >                     sep( T11, T22 ) = sigma-min( C )
*     >
*     >  where sigma-min(C) is the smallest singular value of the
*     >  n1*n2-by-n1*n2 matrix
*     >
*     >     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
*     >
*     >  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
*     >  product. We estimate sigma-min(C) by the reciprocal of an estimate of
*     >  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
*     >  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
*     >
*     >  When SEP is small, small changes in T can cause large changes in
*     >  the invariant subspace. An approximate bound on the maximum angular
*     >  error in the computed right invariant subspace is
*     >
*     >                      EPS * norm(T) / SEP
*     > \endverbatim
*     >
*     =====================================================================
      SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI,
     $     M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
*     
*     -- LAPACK computational routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
      DOUBLE PRECISION   S, SEP
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ),
     $     WR( * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, PAIR, SWAP, WANTBH, WANTQ, WANTS,
     $     WANTSP
      INTEGER            IERR, K, KASE, KK, KS, LIWMIN, LWMIN, N1, N2,
     $     NN
      DOUBLE PRECISION   EST, RNORM, SCALE
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLANGE
      EXTERNAL           LSAME, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DLACPY, DTREXC, DTRSYL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*     
*     Decode and test the input parameters
*     
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )
*     
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP )
     $     THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -8
      ELSE
*     
*     Set M to the dimension of the specified invariant subspace,
*     and test LWORK and LIWORK.
*     
         M = 0
         PAIR = .FALSE.
         DO 10 K = 1, N
            IF( PAIR ) THEN
               PAIR = .FALSE.
            ELSE
               IF( K.LT.N ) THEN
                  IF( T( K+1, K ).EQ.ZERO ) THEN
                     IF( SELECT( K ) )
     $                    M = M + 1
                  ELSE
                     PAIR = .TRUE.
                     IF( SELECT( K ) .OR. SELECT( K+1 ) )
     $                    M = M + 2
                  END IF
               ELSE
                  IF( SELECT( N ) )
     $                 M = M + 1
               END IF
            END IF
 10      CONTINUE
*     
         N1 = M
         N2 = N - M
         NN = N1*N2
*     
         IF( WANTSP ) THEN
            LWMIN = MAX( 1, 2*NN )
            LIWMIN = MAX( 1, NN )
         ELSE IF( LSAME( JOB, 'N' ) ) THEN
            LWMIN = MAX( 1, N )
            LIWMIN = 1
         ELSE IF( LSAME( JOB, 'E' ) ) THEN
            LWMIN = MAX( 1, NN )
            LIWMIN = 1
         END IF
*     
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -15
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -17
         END IF
      END IF
*     
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*     
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*     
*     Quick return if possible.
*     
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS )
     $        S = ONE
         IF( WANTSP )
     $        SEP = DLANGE( '1', N, N, T, LDT, WORK )
         GO TO 40
      END IF
*     
*     Collect the selected blocks at the top-left corner of T.
*     
      KS = 0
      PAIR = .FALSE.
      DO 20 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         ELSE
            SWAP = SELECT( K )
            IF( K.LT.N ) THEN
               IF( T( K+1, K ).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT( K+1 )
               END IF
            END IF
            IF( SWAP ) THEN
               KS = KS + 1
*     
*     Swap the K-th block to position KS.
*     
               IERR = 0
               KK = K
               IF( K.NE.KS )
     $              CALL DTREXC( COMPQ, N, T, LDT, Q, LDQ, KK, KS, WORK,
     $              IERR )
               IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
*     
*     Blocks too close to swap: exit.
*     
                  INFO = 1
                  IF( WANTS )
     $                 S = ZERO
                  IF( WANTSP )
     $                 SEP = ZERO
                  GO TO 40
               END IF
               IF( PAIR )
     $              KS = KS + 1
            END IF
         END IF
 20   CONTINUE
*     
      IF( WANTS ) THEN
*     
*     Solve Sylvester equation for R:
*     
*     T11*R - R*T22 = scale*T12
*     
         CALL DLACPY( 'F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 )
         CALL DTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ),
     $        LDT, WORK, N1, SCALE, IERR )
*     
*     Estimate the reciprocal of the condition number of the cluster
*     of eigenvalues.
*     
         RNORM = DLANGE( 'F', N1, N2, WORK, N1, WORK )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )*
     $           SQRT( RNORM ) )
         END IF
      END IF
*     
      IF( WANTSP ) THEN
*     
*     Estimate sep(T11,T22).
*     
         EST = ZERO
         KASE = 0
 30      CONTINUE
         CALL DLACN2( NN, WORK( NN+1 ), WORK, IWORK, EST, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*     
*     Solve  T11*R - R*T22 = scale*X.
*     
               CALL DTRSYL( 'N', 'N', -1, N1, N2, T, LDT,
     $              T( N1+1, N1+1 ), LDT, WORK, N1, SCALE,
     $              IERR )
            ELSE
*     
*     Solve T11**T*R - R*T22**T = scale*X.
*     
               CALL DTRSYL( 'T', 'T', -1, N1, N2, T, LDT,
     $              T( N1+1, N1+1 ), LDT, WORK, N1, SCALE,
     $              IERR )
            END IF
            GO TO 30
         END IF
*     
         SEP = SCALE / EST
      END IF
*     
 40   CONTINUE
*     
*     Store the output eigenvalues in WR and WI.
*     
      DO 50 K = 1, N
         WR( K ) = T( K, K )
         WI( K ) = ZERO
 50   CONTINUE
      DO 60 K = 1, N - 1
         IF( T( K+1, K ).NE.ZERO ) THEN
            WI( K ) = SQRT( ABS( T( K, K+1 ) ) )*
     $           SQRT( ABS( T( K+1, K ) ) )
            WI( K+1 ) = -WI( K )
         END IF
 60   CONTINUE
*     
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*     
      RETURN
*     
*     End of DTRSEN
*     
      END
*     > \brief \b DLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DLACN2 + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacn2.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacn2.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacn2.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
*     
*     .. Scalar Arguments ..
*     INTEGER            KASE, N
*     DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..
*     INTEGER            ISGN( * ), ISAVE( 3 )
*     DOUBLE PRECISION   V( * ), X( * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DLACN2 estimates the 1-norm of a square, real matrix A.
*     > Reverse communication is used for evaluating matrix-vector products.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >         The order of the matrix.  N >= 1.
*     > \endverbatim
*     >
*     > \param[out] V
*     > \verbatim
*     >          V is DOUBLE PRECISION array, dimension (N)
*     >         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*     >         (W is not returned).
*     > \endverbatim
*     >
*     > \param[in,out] X
*     > \verbatim
*     >          X is DOUBLE PRECISION array, dimension (N)
*     >         On an intermediate return, X should be overwritten by
*     >               A * X,   if KASE=1,
*     >               A**T * X,  if KASE=2,
*     >         and DLACN2 must be re-called with all the other parameters
*     >         unchanged.
*     > \endverbatim
*     >
*     > \param[out] ISGN
*     > \verbatim
*     >          ISGN is INTEGER array, dimension (N)
*     > \endverbatim
*     >
*     > \param[in,out] EST
*     > \verbatim
*     >          EST is DOUBLE PRECISION
*     >         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
*     >         unchanged from the previous call to DLACN2.
*     >         On exit, EST is an estimate (a lower bound) for norm(A).
*     > \endverbatim
*     >
*     > \param[in,out] KASE
*     > \verbatim
*     >          KASE is INTEGER
*     >         On the initial call to DLACN2, KASE should be 0.
*     >         On an intermediate return, KASE will be 1 or 2, indicating
*     >         whether X should be overwritten by A * X  or A**T * X.
*     >         On the final return from DLACN2, KASE will again be 0.
*     > \endverbatim
*     >
*     > \param[in,out] ISAVE
*     > \verbatim
*     >          ISAVE is INTEGER array, dimension (3)
*     >         ISAVE is used to save variables between calls to DLACN2
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleOTHERauxiliary
*     
*     > \par Further Details:
*     =====================
*     >
*     > \verbatim
*     >
*     >  Originally named SONEST, dated March 16, 1988.
*     >
*     >  This is a thread safe version of DLACON, which uses the array ISAVE
*     >  in place of a SAVE statement, as follows:
*     >
*     >     DLACON     DLACN2
*     >      JUMP     ISAVE(1)
*     >      J        ISAVE(2)
*     >      ITER     ISAVE(3)
*     > \endverbatim
*     
*     > \par Contributors:
*     ==================
*     >
*     >     Nick Higham, University of Manchester
*     
*     > \par References:
*     ================
*     >
*     >  N.J. Higham, "FORTRAN codes for estimating the one-norm of
*     >  a real or complex matrix, with applications to condition estimation",
*     >  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*     >
*     =====================================================================
      SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
*     
*     -- LAPACK auxiliary routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..
      INTEGER            ISGN( * ), ISAVE( 3 )
      DOUBLE PRECISION   V( * ), X( * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, JLAST
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP, XS
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM
      EXTERNAL           IDAMAX, DASUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT
*     ..
*     .. Executable Statements ..
*     
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
 10      CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      END IF
*     
      GO TO ( 20, 40, 70, 110, 140 )ISAVE( 1 )
*     
*     ................ ENTRY   (ISAVE( 1 ) = 1)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
*     
 20   CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
*     ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
*     
      DO 30 I = 1, N
         IF( X(I).GE.ZERO ) THEN
            X(I) = ONE
         ELSE
            X(I) = -ONE
         END IF
         ISGN( I ) = NINT( X( I ) )
 30   CONTINUE
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN
*     
*     ................ ENTRY   (ISAVE( 1 ) = 2)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*     
 40   CONTINUE
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      ISAVE( 3 ) = 2
*     
*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
*     
 50   CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
 60   CONTINUE
      X( ISAVE( 2 ) ) = ONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN
*     
*     ................ ENTRY   (ISAVE( 1 ) = 3)
*     X HAS BEEN OVERWRITTEN BY A*X.
*     
 70   CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( X(I).GE.ZERO ) THEN
            XS = ONE
         ELSE
            XS = -ONE
         END IF
         IF( NINT( XS ).NE.ISGN( I ) )
     $        GO TO 90
 80   CONTINUE
*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
*     
 90   CONTINUE
*     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )
     $     GO TO 120
*     
      DO 100 I = 1, N
         IF( X(I).GE.ZERO ) THEN
            X(I) = ONE
         ELSE
            X(I) = -ONE
         END IF
         ISGN( I ) = NINT( X( I ) )
 100  CONTINUE
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN
*     
*     ................ ENTRY   (ISAVE( 1 ) = 4)
*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*     
 110  CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND.
     $     ( ISAVE( 3 ).LT.ITMAX ) ) THEN
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      END IF
*     
*     ITERATION COMPLETE.  FINAL STAGE.
*     
 120  CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
 130  CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN
*     
*     ................ ENTRY   (ISAVE( 1 ) = 5)
*     X HAS BEEN OVERWRITTEN BY A*X.
*     
 140  CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
*     
 150  CONTINUE
      KASE = 0
      RETURN
*     
*     End of DLACN2
*     
      END
*     > \brief \b DTREXC
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DTREXC + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrexc.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrexc.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrexc.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
*     INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER          COMPQ
*     INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DTREXC reorders the real Schur factorization of a real matrix
*     > A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
*     > moved to row ILST.
*     >
*     > The real Schur form T is reordered by an orthogonal similarity
*     > transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
*     > is updated by postmultiplying it with Z.
*     >
*     > T must be in Schur canonical form (as returned by DHSEQR), that is,
*     > block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*     > 2-by-2 diagonal block has its diagonal elements equal and its
*     > off-diagonal elements of opposite sign.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] COMPQ
*     > \verbatim
*     >          COMPQ is CHARACTER*1
*     >          = 'V':  update the matrix Q of Schur vectors;
*     >          = 'N':  do not update Q.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The order of the matrix T. N >= 0.
*     >          If N == 0 arguments ILST and IFST may be any value.
*     > \endverbatim
*     >
*     > \param[in,out] T
*     > \verbatim
*     >          T is DOUBLE PRECISION array, dimension (LDT,N)
*     >          On entry, the upper quasi-triangular matrix T, in Schur
*     >          Schur canonical form.
*     >          On exit, the reordered upper quasi-triangular matrix, again
*     >          in Schur canonical form.
*     > \endverbatim
*     >
*     > \param[in] LDT
*     > \verbatim
*     >          LDT is INTEGER
*     >          The leading dimension of the array T. LDT >= max(1,N).
*     > \endverbatim
*     >
*     > \param[in,out] Q
*     > \verbatim
*     >          Q is DOUBLE PRECISION array, dimension (LDQ,N)
*     >          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*     >          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*     >          orthogonal transformation matrix Z which reorders T.
*     >          If COMPQ = 'N', Q is not referenced.
*     > \endverbatim
*     >
*     > \param[in] LDQ
*     > \verbatim
*     >          LDQ is INTEGER
*     >          The leading dimension of the array Q.  LDQ >= 1, and if
*     >          COMPQ = 'V', LDQ >= max(1,N).
*     > \endverbatim
*     >
*     > \param[in,out] IFST
*     > \verbatim
*     >          IFST is INTEGER
*     > \endverbatim
*     >
*     > \param[in,out] ILST
*     > \verbatim
*     >          ILST is INTEGER
*     >
*     >          Specify the reordering of the diagonal blocks of T.
*     >          The block with row index IFST is moved to row ILST, by a
*     >          sequence of transpositions between adjacent blocks.
*     >          On exit, if IFST pointed on entry to the second row of a
*     >          2-by-2 block, it is changed to point to the first row; ILST
*     >          always points to the first row of the block in its final
*     >          position (which may differ from its input value by +1 or -1).
*     >          1 <= IFST <= N; 1 <= ILST <= N.
*     > \endverbatim
*     >
*     > \param[out] WORK
*     > \verbatim
*     >          WORK is DOUBLE PRECISION array, dimension (N)
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0:  successful exit
*     >          < 0:  if INFO = -i, the i-th argument had an illegal value
*     >          = 1:  two adjacent blocks were too close to swap (the problem
*     >                is very ill-conditioned); T may have been partially
*     >                reordered, and ILST points to the first row of the
*     >                current position of the block being moved.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleOTHERcomputational
*     
*     =====================================================================
      SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
     $     INFO )
*     
*     -- LAPACK computational routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAEXC, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*     
*     Decode and test the input arguments.
*     
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF(( IFST.LT.1 .OR. IFST.GT.N ).AND.( N.GT.0 )) THEN
         INFO = -7
      ELSE IF(( ILST.LT.1 .OR. ILST.GT.N ).AND.( N.GT.0 )) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREXC', -INFO )
         RETURN
      END IF
*     
*     Quick return if possible
*     
      IF( N.LE.1 )
     $     RETURN
*     
*     Determine the first row of specified block
*     and find out it is 1 by 1 or 2 by 2.
*     
      IF( IFST.GT.1 ) THEN
         IF( T( IFST, IFST-1 ).NE.ZERO )
     $        IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST.LT.N ) THEN
         IF( T( IFST+1, IFST ).NE.ZERO )
     $        NBF = 2
      END IF
*     
*     Determine the first row of the final block
*     and find out it is 1 by 1 or 2 by 2.
*     
      IF( ILST.GT.1 ) THEN
         IF( T( ILST, ILST-1 ).NE.ZERO )
     $        ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST.LT.N ) THEN
         IF( T( ILST+1, ILST ).NE.ZERO )
     $        NBL = 2
      END IF
*     
      IF( IFST.EQ.ILST )
     $     RETURN
*     
      IF( IFST.LT.ILST ) THEN
*     
*     Update ILST
*     
         IF( NBF.EQ.2 .AND. NBL.EQ.1 )
     $        ILST = ILST - 1
         IF( NBF.EQ.1 .AND. NBL.EQ.2 )
     $        ILST = ILST + 1
*     
         HERE = IFST
*     
 10      CONTINUE
*     
*     Swap block with next one below
*     
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*     
*     Current block either 1 by 1 or 2 by 2
*     
            NBNEXT = 1
            IF( HERE+NBF+1.LE.N ) THEN
               IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO )
     $              NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT,
     $           WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE + NBNEXT
*     
*     Test if 2 by 2 block breaks into two 1 by 1 blocks
*     
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO )
     $              NBF = 3
            END IF
*     
         ELSE
*     
*     Current block consists of two 1 by 1 blocks each of which
*     must be swapped individually
*     
            NBNEXT = 1
            IF( HERE+3.LE.N ) THEN
               IF( T( HERE+3, HERE+2 ).NE.ZERO )
     $              NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT,
     $           WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
*     
*     Swap two 1 by 1 blocks, no problems possible
*     
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT,
     $              WORK, INFO )
               HERE = HERE + 1
            ELSE
*     
*     Recompute NBNEXT in case 2 by 2 split
*     
               IF( T( HERE+2, HERE+1 ).EQ.ZERO )
     $              NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*     
*     2 by 2 Block did not split
*     
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1,
     $                 NBNEXT, WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 2
               ELSE
*     
*     2 by 2 Block did split
*     
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1,
     $                 WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1,
     $                 WORK, INFO )
                  HERE = HERE + 2
               END IF
            END IF
         END IF
         IF( HERE.LT.ILST )
     $        GO TO 10
*     
      ELSE
*     
         HERE = IFST
 20      CONTINUE
*     
*     Swap block with next one above
*     
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*     
*     Current block either 1 by 1 or 2 by 2
*     
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO )
     $              NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT,
     $           NBF, WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE - NBNEXT
*     
*     Test if 2 by 2 block breaks into two 1 by 1 blocks
*     
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO )
     $              NBF = 3
            END IF
*     
         ELSE
*     
*     Current block consists of two 1 by 1 blocks each of which
*     must be swapped individually
*     
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO )
     $              NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT,
     $           1, WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
*     
*     Swap two 1 by 1 blocks, no problems possible
*     
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1,
     $              WORK, INFO )
               HERE = HERE - 1
            ELSE
*     
*     Recompute NBNEXT in case 2 by 2 split
*     
               IF( T( HERE, HERE-1 ).EQ.ZERO )
     $              NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*     
*     2 by 2 Block did not split
*     
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1,
     $                 WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 2
               ELSE
*     
*     2 by 2 Block did split
*     
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1,
     $                 WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1,
     $                 WORK, INFO )
                  HERE = HERE - 2
               END IF
            END IF
         END IF
         IF( HERE.GT.ILST )
     $        GO TO 20
      END IF
      ILST = HERE
*     
      RETURN
*     
*     End of DTREXC
*     
      END
*     > \brief \b DTRSYL
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DTRSYL + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrsyl.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrsyl.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrsyl.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
*     LDC, SCALE, INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER          TRANA, TRANB
*     INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
*     DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DTRSYL solves the real Sylvester matrix equation:
*     >
*     >    op(A)*X + X*op(B) = scale*C or
*     >    op(A)*X - X*op(B) = scale*C,
*     >
*     > where op(A) = A or A**T, and  A and B are both upper quasi-
*     > triangular. A is M-by-M and B is N-by-N; the right hand side C and
*     > the solution X are M-by-N; and scale is an output scale factor, set
*     > <= 1 to avoid overflow in X.
*     >
*     > A and B must be in Schur canonical form (as returned by DHSEQR), that
*     > is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
*     > each 2-by-2 diagonal block has its diagonal elements equal and its
*     > off-diagonal elements of opposite sign.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] TRANA
*     > \verbatim
*     >          TRANA is CHARACTER*1
*     >          Specifies the option op(A):
*     >          = 'N': op(A) = A    (No transpose)
*     >          = 'T': op(A) = A**T (Transpose)
*     >          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
*     > \endverbatim
*     >
*     > \param[in] TRANB
*     > \verbatim
*     >          TRANB is CHARACTER*1
*     >          Specifies the option op(B):
*     >          = 'N': op(B) = B    (No transpose)
*     >          = 'T': op(B) = B**T (Transpose)
*     >          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
*     > \endverbatim
*     >
*     > \param[in] ISGN
*     > \verbatim
*     >          ISGN is INTEGER
*     >          Specifies the sign in the equation:
*     >          = +1: solve op(A)*X + X*op(B) = scale*C
*     >          = -1: solve op(A)*X - X*op(B) = scale*C
*     > \endverbatim
*     >
*     > \param[in] M
*     > \verbatim
*     >          M is INTEGER
*     >          The order of the matrix A, and the number of rows in the
*     >          matrices X and C. M >= 0.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The order of the matrix B, and the number of columns in the
*     >          matrices X and C. N >= 0.
*     > \endverbatim
*     >
*     > \param[in] A
*     > \verbatim
*     >          A is DOUBLE PRECISION array, dimension (LDA,M)
*     >          The upper quasi-triangular matrix A, in Schur canonical form.
*     > \endverbatim
*     >
*     > \param[in] LDA
*     > \verbatim
*     >          LDA is INTEGER
*     >          The leading dimension of the array A. LDA >= max(1,M).
*     > \endverbatim
*     >
*     > \param[in] B
*     > \verbatim
*     >          B is DOUBLE PRECISION array, dimension (LDB,N)
*     >          The upper quasi-triangular matrix B, in Schur canonical form.
*     > \endverbatim
*     >
*     > \param[in] LDB
*     > \verbatim
*     >          LDB is INTEGER
*     >          The leading dimension of the array B. LDB >= max(1,N).
*     > \endverbatim
*     >
*     > \param[in,out] C
*     > \verbatim
*     >          C is DOUBLE PRECISION array, dimension (LDC,N)
*     >          On entry, the M-by-N right hand side matrix C.
*     >          On exit, C is overwritten by the solution matrix X.
*     > \endverbatim
*     >
*     > \param[in] LDC
*     > \verbatim
*     >          LDC is INTEGER
*     >          The leading dimension of the array C. LDC >= max(1,M)
*     > \endverbatim
*     >
*     > \param[out] SCALE
*     > \verbatim
*     >          SCALE is DOUBLE PRECISION
*     >          The scale factor, scale, set <= 1 to avoid overflow in X.
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0: successful exit
*     >          < 0: if INFO = -i, the i-th argument had an illegal value
*     >          = 1: A and B have common or very close eigenvalues; perturbed
*     >               values were used to solve the equation (but the matrices
*     >               A and B are unchanged).
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleSYcomputational
*     
*     =====================================================================
      SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $     LDC, SCALE, INFO )
*     
*     -- LAPACK computational routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,
     $     SMLNUM, SUML, SUMR, XNORM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANGE
      EXTERNAL           LSAME, DDOT, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, DLALN2, DLASY2, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*     
*     Decode and Test input parameters
*     
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
*     
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT.
     $     LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT.
     $        LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN.NE.1 .AND. ISGN.NE.-1 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSYL', -INFO )
         RETURN
      END IF
*     
*     Quick return if possible
*     
      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 )
     $     RETURN
*     
*     Set constants to control overflow
*     
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
*     
      SMIN = MAX( SMLNUM, EPS*DLANGE( 'M', M, M, A, LDA, DUM ),
     $     EPS*DLANGE( 'M', N, N, B, LDB, DUM ) )
*     
      SGN = ISGN
*     
      IF( NOTRNA .AND. NOTRNB ) THEN
*     
*     Solve    A*X + ISGN*X*B = scale*C.
*     
*     The (K,L)th block of X is determined starting from
*     bottom-left corner column by column by
*     
*     A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
*     
*     Where
*     M                         L-1
*     R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
*     I=K+1                       J=1
*     
*     Start column loop (index = L)
*     L1 (L2) : column index of the first (first) row of X(K,L).
*     
         LNEXT = 1
         DO 60 L = 1, N
            IF( L.LT.LNEXT )
     $           GO TO 60
            IF( L.EQ.N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L ).NE.ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
*     
*     Start row loop (index = K)
*     K1 (K2): row index of the first (last) row of X(K,L).
*     
            KNEXT = M
            DO 50 K = M, 1, -1
               IF( K.GT.KNEXT )
     $              GO TO 50
               IF( K.EQ.1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 ).NE.ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
*     
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                 C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*     
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                    SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 10                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*     
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                 LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 20 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 20                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*     
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                 C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*     
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                 C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*     
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                 LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 30 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 30                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLASY2( .FALSE., .FALSE., ISGN, 2, 2,
     $                 A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC,
     $                 2, SCALOC, X, 2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 40 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 40                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*     
 50         CONTINUE
*     
 60      CONTINUE
*     
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
*     
*     Solve    A**T *X + ISGN*X*B = scale*C.
*     
*     The (K,L)th block of X is determined starting from
*     upper-left corner column by column by
*     
*     A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
*     
*     Where
*     K-1        T                    L-1
*     R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
*     I=1                          J=1
*     
*     Start column loop (index = L)
*     L1 (L2): column index of the first (last) row of X(K,L)
*     
         LNEXT = 1
         DO 120 L = 1, N
            IF( L.LT.LNEXT )
     $           GO TO 120
            IF( L.EQ.N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L ).NE.ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
*     
*     Start row loop (index = K)
*     K1 (K2): row index of the first (last) row of X(K,L)
*     
            KNEXT = 1
            DO 110 K = 1, M
               IF( K.LT.KNEXT )
     $              GO TO 110
               IF( K.EQ.M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K ).NE.ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
*     
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*     
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                    SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 70 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 70                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*     
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                 LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 80                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*     
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                 LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 90 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 90                  CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLASY2( .TRUE., .FALSE., ISGN, 2, 2, A( K1, K1 ),
     $                 LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                 2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 100 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 100                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*     
 110        CONTINUE
 120     CONTINUE
*     
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
*     
*     Solve    A**T*X + ISGN*X*B**T = scale*C.
*     
*     The (K,L)th block of X is determined starting from
*     top-right corner column by column by
*     
*     A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
*     
*     Where
*     K-1                            N
*     R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
*     I=1                          J=L+1
*     
*     Start column loop (index = L)
*     L1 (L2): column index of the first (last) row of X(K,L)
*     
         LNEXT = N
         DO 180 L = N, 1, -1
            IF( L.GT.LNEXT )
     $           GO TO 180
            IF( L.EQ.1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 ).NE.ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
*     
*     Start row loop (index = K)
*     K1 (K2): row index of the first (last) row of X(K,L)
*     
            KNEXT = 1
            DO 170 K = 1, M
               IF( K.LT.KNEXT )
     $              GO TO 170
               IF( K.EQ.M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K ).NE.ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
*     
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC,
     $                 B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*     
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                    SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 130 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 130                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*     
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                 LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 140 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 140                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*     
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                 LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 150 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 150                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                 B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLASY2( .TRUE., .TRUE., ISGN, 2, 2, A( K1, K1 ),
     $                 LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                 2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 160 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 160                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*     
 170        CONTINUE
 180     CONTINUE
*     
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
*     
*     Solve    A*X + ISGN*X*B**T = scale*C.
*     
*     The (K,L)th block of X is determined starting from
*     bottom-right corner column by column by
*     
*     A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
*     
*     Where
*     M                          N
*     R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
*     I=K+1                      J=L+1
*     
*     Start column loop (index = L)
*     L1 (L2): column index of the first (last) row of X(K,L)
*     
         LNEXT = N
         DO 240 L = N, 1, -1
            IF( L.GT.LNEXT )
     $           GO TO 240
            IF( L.EQ.1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 ).NE.ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
*     
*     Start row loop (index = K)
*     K1 (K2): row index of the first (last) row of X(K,L)
*     
            KNEXT = M
            DO 230 K = M, 1, -1
               IF( K.GT.KNEXT )
     $              GO TO 230
               IF( K.EQ.1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 ).NE.ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
*     
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                 C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC,
     $                 B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*     
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                    SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 190 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 190                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*     
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                 LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 200 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 200                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*     
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                 C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*     
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                 C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*     
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                 LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                 ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 210 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 210                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*     
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*     
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                 B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                 B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*     
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                 C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                 B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*     
                  CALL DLASY2( .FALSE., .TRUE., ISGN, 2, 2, A( K1, K1 ),
     $                 LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                 2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $                 INFO = 1
*     
                  IF( SCALOC.NE.ONE ) THEN
                     DO 220 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
 220                 CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*     
 230        CONTINUE
 240     CONTINUE
*     
      END IF
*     
      RETURN
*     
*     End of DTRSYL
*     
      END
*     > \brief \b DTRTRS
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DTRTRS + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrtrs.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrtrs.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrtrs.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
*     INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER          DIAG, TRANS, UPLO
*     INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DTRTRS solves a triangular system of the form
*     >
*     >    A * X = B  or  A**T * X = B,
*     >
*     > where A is a triangular matrix of order N, and B is an N-by-NRHS
*     > matrix.  A check is made to verify that A is nonsingular.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] UPLO
*     > \verbatim
*     >          UPLO is CHARACTER*1
*     >          = 'U':  A is upper triangular;
*     >          = 'L':  A is lower triangular.
*     > \endverbatim
*     >
*     > \param[in] TRANS
*     > \verbatim
*     >          TRANS is CHARACTER*1
*     >          Specifies the form of the system of equations:
*     >          = 'N':  A * X = B  (No transpose)
*     >          = 'T':  A**T * X = B  (Transpose)
*     >          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*     > \endverbatim
*     >
*     > \param[in] DIAG
*     > \verbatim
*     >          DIAG is CHARACTER*1
*     >          = 'N':  A is non-unit triangular;
*     >          = 'U':  A is unit triangular.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The order of the matrix A.  N >= 0.
*     > \endverbatim
*     >
*     > \param[in] NRHS
*     > \verbatim
*     >          NRHS is INTEGER
*     >          The number of right hand sides, i.e., the number of columns
*     >          of the matrix B.  NRHS >= 0.
*     > \endverbatim
*     >
*     > \param[in] A
*     > \verbatim
*     >          A is DOUBLE PRECISION array, dimension (LDA,N)
*     >          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
*     >          upper triangular part of the array A contains the upper
*     >          triangular matrix, and the strictly lower triangular part of
*     >          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
*     >          triangular part of the array A contains the lower triangular
*     >          matrix, and the strictly upper triangular part of A is not
*     >          referenced.  If DIAG = 'U', the diagonal elements of A are
*     >          also not referenced and are assumed to be 1.
*     > \endverbatim
*     >
*     > \param[in] LDA
*     > \verbatim
*     >          LDA is INTEGER
*     >          The leading dimension of the array A.  LDA >= max(1,N).
*     > \endverbatim
*     >
*     > \param[in,out] B
*     > \verbatim
*     >          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*     >          On entry, the right hand side matrix B.
*     >          On exit, if INFO = 0, the solution matrix X.
*     > \endverbatim
*     >
*     > \param[in] LDB
*     > \verbatim
*     >          LDB is INTEGER
*     >          The leading dimension of the array B.  LDB >= max(1,N).
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0:  successful exit
*     >          < 0: if INFO = -i, the i-th argument had an illegal value
*     >          > 0: if INFO = i, the i-th diagonal element of A is zero,
*     >               indicating that the matrix is singular and the solutions
*     >               X have not been computed.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleOTHERcomputational
*     
*     =====================================================================
      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $     INFO )
*     
*     -- LAPACK computational routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*     
*     Test the input parameters.
*     
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.
     $        LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRS', -INFO )
         RETURN
      END IF
*     
*     Quick return if possible
*     
      IF( N.EQ.0 )
     $     RETURN
*     
*     Check for singularity.
*     
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $           RETURN
 10      CONTINUE
      END IF
      INFO = 0
*     
*     Solve A * x = b  or  A**T * x = b.
*     
      CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,
     $     LDB )
*     
      RETURN
*     
*     End of DTRTRS
*     
      END
*     > \brief \b DLAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation.
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DLAEXC + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaexc.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaexc.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaexc.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
*     INFO )
*     
*     .. Scalar Arguments ..
*     LOGICAL            WANTQ
*     INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
*     > an upper quasi-triangular matrix T by an orthogonal similarity
*     > transformation.
*     >
*     > T must be in Schur canonical form, that is, block upper triangular
*     > with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
*     > has its diagonal elements equal and its off-diagonal elements of
*     > opposite sign.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] WANTQ
*     > \verbatim
*     >          WANTQ is LOGICAL
*     >          = .TRUE. : accumulate the transformation in the matrix Q;
*     >          = .FALSE.: do not accumulate the transformation.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     >          N is INTEGER
*     >          The order of the matrix T. N >= 0.
*     > \endverbatim
*     >
*     > \param[in,out] T
*     > \verbatim
*     >          T is DOUBLE PRECISION array, dimension (LDT,N)
*     >          On entry, the upper quasi-triangular matrix T, in Schur
*     >          canonical form.
*     >          On exit, the updated matrix T, again in Schur canonical form.
*     > \endverbatim
*     >
*     > \param[in] LDT
*     > \verbatim
*     >          LDT is INTEGER
*     >          The leading dimension of the array T. LDT >= max(1,N).
*     > \endverbatim
*     >
*     > \param[in,out] Q
*     > \verbatim
*     >          Q is DOUBLE PRECISION array, dimension (LDQ,N)
*     >          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
*     >          On exit, if WANTQ is .TRUE., the updated matrix Q.
*     >          If WANTQ is .FALSE., Q is not referenced.
*     > \endverbatim
*     >
*     > \param[in] LDQ
*     > \verbatim
*     >          LDQ is INTEGER
*     >          The leading dimension of the array Q.
*     >          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
*     > \endverbatim
*     >
*     > \param[in] J1
*     > \verbatim
*     >          J1 is INTEGER
*     >          The index of the first row of the first block T11.
*     > \endverbatim
*     >
*     > \param[in] N1
*     > \verbatim
*     >          N1 is INTEGER
*     >          The order of the first block T11. N1 = 0, 1 or 2.
*     > \endverbatim
*     >
*     > \param[in] N2
*     > \verbatim
*     >          N2 is INTEGER
*     >          The order of the second block T22. N2 = 0, 1 or 2.
*     > \endverbatim
*     >
*     > \param[out] WORK
*     > \verbatim
*     >          WORK is DOUBLE PRECISION array, dimension (N)
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          = 0: successful exit
*     >          = 1: the transformed matrix T would be too far from Schur
*     >               form; the blocks are not swapped and T and Q are
*     >               unchanged.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleOTHERauxiliary
*     
*     =====================================================================
      SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
     $     INFO )
*     
*     -- LAPACK auxiliary routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 1.0D+1 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
     $     T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2,
     $     WR1, WR2, XNORM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ),
     $     X( LDX, 2 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLANV2, DLARFG, DLARFX, DLARTG, DLASY2,
     $     DROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*     
      INFO = 0
*     
*     Quick return if possible
*     
      IF( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 )
     $     RETURN
      IF( J1+N1.GT.N )
     $     RETURN
*     
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
*     
      IF( N1.EQ.1 .AND. N2.EQ.1 ) THEN
*     
*     Swap two 1-by-1 blocks.
*     
         T11 = T( J1, J1 )
         T22 = T( J2, J2 )
*     
*     Determine the transformation to perform the interchange.
*     
         CALL DLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
*     
*     Apply transformation to the matrix T.
*     
         IF( J3.LE.N )
     $        CALL DROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS,
     $        SN )
         CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
*     
         T( J1, J1 ) = T22
         T( J2, J2 ) = T11
*     
         IF( WANTQ ) THEN
*     
*     Accumulate transformation in the matrix Q.
*     
            CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
*     
      ELSE
*     
*     Swapping involves at least one 2-by-2 block.
*     
*     Copy the diagonal block of order N1+N2 to the local array D
*     and compute its norm.
*     
         ND = N1 + N2
         CALL DLACPY( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
         DNORM = DLANGE( 'Max', ND, ND, D, LDD, WORK )
*     
*     Compute machine-dependent threshold for test for accepting
*     swap.
*     
         EPS = DLAMCH( 'P' )
         SMLNUM = DLAMCH( 'S' ) / EPS
         THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
*     
*     Solve T11*X - X*T22 = scale*T12 for X.
*     
         CALL DLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD,
     $        D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X,
     $        LDX, XNORM, IERR )
*     
*     Swap the adjacent diagonal blocks.
*     
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
*     
 10      CONTINUE
*     
*     N1 = 1, N2 = 2: generate elementary reflector H so that:
*     
*     ( scale, X11, X12 ) H = ( 0, 0, * )
*     
         U( 1 ) = SCALE
         U( 2 ) = X( 1, 1 )
         U( 3 ) = X( 1, 2 )
         CALL DLARFG( 3, U( 3 ), U, 1, TAU )
         U( 3 ) = ONE
         T11 = T( J1, J1 )
*     
*     Perform swap provisionally on diagonal block in D.
*     
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
*     
*     Test whether to reject swap.
*     
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3,
     $        3 )-T11 ) ).GT.THRESH )GO TO 50
*     
*     Accept swap: apply transformation to the entire matrix T.
*     
         CALL DLARFX( 'L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK )
*     
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J3, J3 ) = T11
*     
         IF( WANTQ ) THEN
*     
*     Accumulate transformation in the matrix Q.
*     
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
*     
 20      CONTINUE
*     
*     N1 = 2, N2 = 1: generate elementary reflector H so that:
*     
*     H (  -X11 ) = ( * )
*     (  -X21 ) = ( 0 )
*     ( scale ) = ( 0 )
*     
         U( 1 ) = -X( 1, 1 )
         U( 2 ) = -X( 2, 1 )
         U( 3 ) = SCALE
         CALL DLARFG( 3, U( 1 ), U( 2 ), 1, TAU )
         U( 1 ) = ONE
         T33 = T( J3, J3 )
*     
*     Perform swap provisionally on diagonal block in D.
*     
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
*     
*     Test whether to reject swap.
*     
         IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1,
     $        1 )-T33 ) ).GT.THRESH )GO TO 50
*     
*     Accept swap: apply transformation to the entire matrix T.
*     
         CALL DLARFX( 'R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK )
*     
         T( J1, J1 ) = T33
         T( J2, J1 ) = ZERO
         T( J3, J1 ) = ZERO
*     
         IF( WANTQ ) THEN
*     
*     Accumulate transformation in the matrix Q.
*     
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
*     
 30      CONTINUE
*     
*     N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
*     that:
*     
*     H(2) H(1) (  -X11  -X12 ) = (  *  * )
*     (  -X21  -X22 )   (  0  * )
*     ( scale    0  )   (  0  0 )
*     (    0  scale )   (  0  0 )
*     
         U1( 1 ) = -X( 1, 1 )
         U1( 2 ) = -X( 2, 1 )
         U1( 3 ) = SCALE
         CALL DLARFG( 3, U1( 1 ), U1( 2 ), 1, TAU1 )
         U1( 1 ) = ONE
*     
         TEMP = -TAU1*( X( 1, 2 )+U1( 2 )*X( 2, 2 ) )
         U2( 1 ) = -TEMP*U1( 2 ) - X( 2, 2 )
         U2( 2 ) = -TEMP*U1( 3 )
         U2( 3 ) = SCALE
         CALL DLARFG( 3, U2( 1 ), U2( 2 ), 1, TAU2 )
         U2( 1 ) = ONE
*     
*     Perform swap provisionally on diagonal block in D.
*     
         CALL DLARFX( 'L', 3, 4, U1, TAU1, D, LDD, WORK )
         CALL DLARFX( 'R', 4, 3, U1, TAU1, D, LDD, WORK )
         CALL DLARFX( 'L', 3, 4, U2, TAU2, D( 2, 1 ), LDD, WORK )
         CALL DLARFX( 'R', 4, 3, U2, TAU2, D( 1, 2 ), LDD, WORK )
*     
*     Test whether to reject swap.
*     
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ),
     $        ABS( D( 4, 2 ) ) ).GT.THRESH )GO TO 50
*     
*     Accept swap: apply transformation to the entire matrix T.
*     
         CALL DLARFX( 'L', 3, N-J1+1, U1, TAU1, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J4, 3, U1, TAU1, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1+1, U2, TAU2, T( J2, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J4, 3, U2, TAU2, T( 1, J2 ), LDT, WORK )
*     
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J4, J1 ) = ZERO
         T( J4, J2 ) = ZERO
*     
         IF( WANTQ ) THEN
*     
*     Accumulate transformation in the matrix Q.
*     
            CALL DLARFX( 'R', N, 3, U1, TAU1, Q( 1, J1 ), LDQ, WORK )
            CALL DLARFX( 'R', N, 3, U2, TAU2, Q( 1, J2 ), LDQ, WORK )
         END IF
*     
 40      CONTINUE
*     
         IF( N2.EQ.2 ) THEN
*     
*     Standardize new 2-by-2 block T11
*     
            CALL DLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ),
     $           T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
            CALL DROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT,
     $           CS, SN )
            CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
            IF( WANTQ )
     $           CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
*     
         IF( N1.EQ.2 ) THEN
*     
*     Standardize new 2-by-2 block T22
*     
            J3 = J1 + N2
            J4 = J3 + 1
            CALL DLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ),
     $           T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
            IF( J3+2.LE.N )
     $           CALL DROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ),
     $           LDT, CS, SN )
            CALL DROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
            IF( WANTQ )
     $           CALL DROT( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN )
         END IF
*     
      END IF
      RETURN
*     
*     Exit with INFO = 1 if swap was rejected.
*     
 50   CONTINUE
      INFO = 1
      RETURN
*     
*     End of DLAEXC
*     
      END
*     > \brief \b DLASY2 solves the Sylvester matrix equation where the matrices are of order 1 or 2.
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download DLASY2 + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasy2.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasy2.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasy2.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
*     LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
*     
*     .. Scalar Arguments ..
*     LOGICAL            LTRANL, LTRANR
*     INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
*     DOUBLE PRECISION   SCALE, XNORM
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
*     $                   X( LDX, * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
*     >
*     >        op(TL)*X + ISGN*X*op(TR) = SCALE*B,
*     >
*     > where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
*     > -1.  op(T) = T or T**T, where T**T denotes the transpose of T.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] LTRANL
*     > \verbatim
*     >          LTRANL is LOGICAL
*     >          On entry, LTRANL specifies the op(TL):
*     >             = .FALSE., op(TL) = TL,
*     >             = .TRUE., op(TL) = TL**T.
*     > \endverbatim
*     >
*     > \param[in] LTRANR
*     > \verbatim
*     >          LTRANR is LOGICAL
*     >          On entry, LTRANR specifies the op(TR):
*     >            = .FALSE., op(TR) = TR,
*     >            = .TRUE., op(TR) = TR**T.
*     > \endverbatim
*     >
*     > \param[in] ISGN
*     > \verbatim
*     >          ISGN is INTEGER
*     >          On entry, ISGN specifies the sign of the equation
*     >          as described before. ISGN may only be 1 or -1.
*     > \endverbatim
*     >
*     > \param[in] N1
*     > \verbatim
*     >          N1 is INTEGER
*     >          On entry, N1 specifies the order of matrix TL.
*     >          N1 may only be 0, 1 or 2.
*     > \endverbatim
*     >
*     > \param[in] N2
*     > \verbatim
*     >          N2 is INTEGER
*     >          On entry, N2 specifies the order of matrix TR.
*     >          N2 may only be 0, 1 or 2.
*     > \endverbatim
*     >
*     > \param[in] TL
*     > \verbatim
*     >          TL is DOUBLE PRECISION array, dimension (LDTL,2)
*     >          On entry, TL contains an N1 by N1 matrix.
*     > \endverbatim
*     >
*     > \param[in] LDTL
*     > \verbatim
*     >          LDTL is INTEGER
*     >          The leading dimension of the matrix TL. LDTL >= max(1,N1).
*     > \endverbatim
*     >
*     > \param[in] TR
*     > \verbatim
*     >          TR is DOUBLE PRECISION array, dimension (LDTR,2)
*     >          On entry, TR contains an N2 by N2 matrix.
*     > \endverbatim
*     >
*     > \param[in] LDTR
*     > \verbatim
*     >          LDTR is INTEGER
*     >          The leading dimension of the matrix TR. LDTR >= max(1,N2).
*     > \endverbatim
*     >
*     > \param[in] B
*     > \verbatim
*     >          B is DOUBLE PRECISION array, dimension (LDB,2)
*     >          On entry, the N1 by N2 matrix B contains the right-hand
*     >          side of the equation.
*     > \endverbatim
*     >
*     > \param[in] LDB
*     > \verbatim
*     >          LDB is INTEGER
*     >          The leading dimension of the matrix B. LDB >= max(1,N1).
*     > \endverbatim
*     >
*     > \param[out] SCALE
*     > \verbatim
*     >          SCALE is DOUBLE PRECISION
*     >          On exit, SCALE contains the scale factor. SCALE is chosen
*     >          less than or equal to 1 to prevent the solution overflowing.
*     > \endverbatim
*     >
*     > \param[out] X
*     > \verbatim
*     >          X is DOUBLE PRECISION array, dimension (LDX,2)
*     >          On exit, X contains the N1 by N2 solution.
*     > \endverbatim
*     >
*     > \param[in] LDX
*     > \verbatim
*     >          LDX is INTEGER
*     >          The leading dimension of the matrix X. LDX >= max(1,N1).
*     > \endverbatim
*     >
*     > \param[out] XNORM
*     > \verbatim
*     >          XNORM is DOUBLE PRECISION
*     >          On exit, XNORM is the infinity-norm of the solution.
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     >          INFO is INTEGER
*     >          On exit, INFO is set to
*     >             0: successful exit.
*     >             1: TL and TR have too close eigenvalues, so TL or
*     >                TR is perturbed to get a nonsingular equation.
*     >          NOTE: In the interests of speed, this routine does not
*     >                check the inputs for errors.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \ingroup doubleSYauxiliary
*     
*     =====================================================================
      SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
     $     LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
*     
*     -- LAPACK auxiliary routine --
*     -- LAPACK is a software package provided by Univ. of Tennessee,    --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     
*     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
     $     X( LDX, * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO, HALF, EIGHT
      PARAMETER          ( TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1,
     $     TEMP, U11, U12, U22, XMAX
*     ..
*     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ),
     $     LOCU22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Data statements ..
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / ,
     $     LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
*     ..
*     .. Executable Statements ..
*     
*     Do not check the input parameters for errors
*     
      INFO = 0
*     
*     Quick return if possible
*     
      IF( N1.EQ.0 .OR. N2.EQ.0 )
     $     RETURN
*     
*     Set constants to control overflow
*     
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SGN = ISGN
*     
      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K
*     
*     1 by 1: TL11*X + SGN*X*TR11 = B11
*     
 10   CONTINUE
      TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
      BET = ABS( TAU1 )
      IF( BET.LE.SMLNUM ) THEN
         TAU1 = SMLNUM
         BET = SMLNUM
         INFO = 1
      END IF
*     
      SCALE = ONE
      GAM = ABS( B( 1, 1 ) )
      IF( SMLNUM*GAM.GT.BET )
     $     SCALE = ONE / GAM
*     
      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN
*     
*     1 by 2:
*     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
*     [TR21 TR22]
*     
 20   CONTINUE
*     
      SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ),
     $     ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ),
     $     SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      IF( LTRANR ) THEN
         TMP( 2 ) = SGN*TR( 2, 1 )
         TMP( 3 ) = SGN*TR( 1, 2 )
      ELSE
         TMP( 2 ) = SGN*TR( 1, 2 )
         TMP( 3 ) = SGN*TR( 2, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40
*     
*     2 by 1:
*     op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
*     [TL21 TL22] [X21]         [X21]         [B21]
*     
 30   CONTINUE
      SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ),
     $     ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ),
     $     SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      IF( LTRANL ) THEN
         TMP( 2 ) = TL( 1, 2 )
         TMP( 3 ) = TL( 2, 1 )
      ELSE
         TMP( 2 ) = TL( 2, 1 )
         TMP( 3 ) = TL( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
 40   CONTINUE
*     
*     Solve 2 by 2 system using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*     
      IPIV = IDAMAX( 4, TMP, 1 )
      U11 = TMP( IPIV )
      IF( ABS( U11 ).LE.SMIN ) THEN
         INFO = 1
         U11 = SMIN
      END IF
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      IF( ABS( U22 ).LE.SMIN ) THEN
         INFO = 1
         U22 = SMIN
      END IF
      IF( BSWAP ) THEN
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      ELSE
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      END IF
      SCALE = ONE
      IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR.
     $     ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      END IF
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      IF( XSWAP ) THEN
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      END IF
      X( 1, 1 ) = X2( 1 )
      IF( N1.EQ.1 ) THEN
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
      ELSE
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
      END IF
      RETURN
*     
*     2 by 2:
*     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
*     [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
*     
*     Solve equivalent 4 by 4 system using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*     
 50   CONTINUE
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
     $     ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
      SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
     $     ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      BTMP( 1 ) = ZERO
      CALL DCOPY( 16, BTMP, 0, T16, 1 )
      T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
      IF( LTRANL ) THEN
         T16( 1, 2 ) = TL( 2, 1 )
         T16( 2, 1 ) = TL( 1, 2 )
         T16( 3, 4 ) = TL( 2, 1 )
         T16( 4, 3 ) = TL( 1, 2 )
      ELSE
         T16( 1, 2 ) = TL( 1, 2 )
         T16( 2, 1 ) = TL( 2, 1 )
         T16( 3, 4 ) = TL( 1, 2 )
         T16( 4, 3 ) = TL( 2, 1 )
      END IF
      IF( LTRANR ) THEN
         T16( 1, 3 ) = SGN*TR( 1, 2 )
         T16( 2, 4 ) = SGN*TR( 1, 2 )
         T16( 3, 1 ) = SGN*TR( 2, 1 )
         T16( 4, 2 ) = SGN*TR( 2, 1 )
      ELSE
         T16( 1, 3 ) = SGN*TR( 2, 1 )
         T16( 2, 4 ) = SGN*TR( 2, 1 )
         T16( 3, 1 ) = SGN*TR( 1, 2 )
         T16( 4, 2 ) = SGN*TR( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )
*     
*     Perform elimination
*     
      DO 100 I = 1, 3
         XMAX = ZERO
         DO 70 IP = I, 4
            DO 60 JP = I, 4
               IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
 60         CONTINUE
 70      CONTINUE
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $        CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T16( I, I ) = SMIN
         END IF
         DO 90 J = I + 1, 4
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
            DO 80 K = I + 1, 4
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
 80         CONTINUE
 90      CONTINUE
 100  CONTINUE
      IF( ABS( T16( 4, 4 ) ).LT.SMIN ) THEN
         INFO = 1
         T16( 4, 4 ) = SMIN
      END IF
      SCALE = ONE
      IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR.
     $     ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR.
     $     ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR.
     $     ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
      SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ),
     $     ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
      BTMP( 1 ) = BTMP( 1 )*SCALE
      BTMP( 2 ) = BTMP( 2 )*SCALE
      BTMP( 3 ) = BTMP( 3 )*SCALE
      BTMP( 4 ) = BTMP( 4 )*SCALE
      END IF
      DO 120 I = 1, 4
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 110 J = K + 1, 4
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
 110     CONTINUE
 120  CONTINUE
      DO 130 I = 1, 3
         IF( JPIV( 4-I ).NE.4-I ) THEN
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         END IF
 130  CONTINUE
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ),
     $     ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
      RETURN
*     
*     End of DLASY2
*     
      END
