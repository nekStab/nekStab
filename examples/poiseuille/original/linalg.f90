!=======================================================================
! Name        : linalg
! Author      : Prabal S. Negi, Mattias Brynjell-Rahkola
! Version     : last modification 2018.02.22
! Copyright   : GPL
! Description : set of linear algebra routines and wrappers for LAPACK
!======================================================================
!---------------------------------------------------------------------- 
!     
!     LAPACK interface for SVD.
!     Upon finishing, SIGMA contains the singular values, VMATX n left
!     singular vectors and VMATXT n transposed right singular vectors.
!
!     Dongarra et al. (1999)
!
      subroutine svd_wrapper(nrow,ncol,dtype)

      implicit none

      include 'SIZE'
      include 'OTD'
      include 'WLAPACK'      
      
      character*1 jobz    !  Specifies options for computing all or part of the matrix U:
                          ! 'A':  all M columns of U and all N rows of V**T are
                          !    returned in the arrays U and VMATXT;
                          ! 'S':  the first min(M,N) columns of U and the first
                          !   min(M,N) rows of V**T are returned in the arrays U and VMATXT;
                          ! 'O':  If M >= N, the first N columns of U are overwritten
                          !   on the array A and all rows of V**T are returned in the array VMATXT;
                          !   otherwise, all columns of U are returned in the array U 
                          !   and the first M rows of V**T are overwritten in the array A;
                          ! 'N':  no columns of U or rows of V**T are computed.

      integer info        ! = 0:  successful exit.
                          ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                          ! > 0:  DBDSDC did not converge, updating process failed.                          

      character*1 job     !  Specifies for which problem the reciprocal condition numbers
                          !  should be computed:
                          ! = 'E':  the eigenvectors of a symmetric/Hermitian matrix;
                          ! = 'L':  the left singular vectors of a general matrix;
                          ! = 'R':  the right singular vectors of a general matrix.
      
      integer m, n, lda, ldu, ldvt,
     $     nrow, ncol
      real u,
     $     epsmch, serrbnd
      character*3 dtype
      character*23 fname

      integer i
      real dlamch       ! function

!     Define LAPACK-variables
      jobz = 'O'                ! return the n first singular vectors 
      m    = nrow               ! no rows in matrix
      n    = ncol               ! no columns in matrix
      lda  = LPERT              ! leading dimension of the matrix
      u    = 0.0                ! left singular vectors (not referenced)
      ldu  = 1                  ! leading dimension of u
      ldvt = LPERT              ! leading dimension of VMATXT

!     VMATX                      ! compute SVD for this matrix.
!     SIGMA                     ! Singular values of VMATX      

!     Compute SVD in double precision with divide-and-conquer
      call dgesdd(jobz,m,n,VMATX,lda,OSIGMA,u,ldu,VMATXT,ldvt,
     $     RWORK,LWORKR,IWORK,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: DBDSDC did not converge, updating process failed.'
     $        , info
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGESDD: successful exit!'
         if (nid.eq.0) write(6,*) '         Optimal LWORKR=',
     $        int(RWORK(1)), LWORKR
      endif

      epsmch = dlamch('E')      ! Machine epsilon in double precision
      if (nid.eq.0) write(6,*) 'Relative machine precision:', epsmch

!     Compute reciprocal condition numbers for singular vectors in double precision
      job = 'l'
      call ddisna(job,m,n,OSIGMA,RCL,info)
!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DDISNA: successful exit!'
      endif

      job = 'r'
      call ddisna(job,m,n,OSIGMA,RCR,info)
!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DDISNA: successful exit!'
      endif
            
!     Output singular values and estimated error bounds
      write(fname,'(A3,A20)') dtype,'_singular_values.txt'
      if (nid.eq.0) then
         open(61,file=fname,action='write',status='unknown')
         write(61,'(TR2,A2,TR18,A6,TR16,A8,TR16,A8,TR16,A8)')
     $        'I;','sigma;','serrbnd;','uerrbnd;','verrbnd'

!     see ch.4.9 Dongarra et al. (1999)
!         serrbnd = epsmch*SIGMA(1)
         do i=1,min(m,n)
            write(61,'(I4,4G24.16)') i, OSIGMA(i), serrbnd,
     $           serrbnd/RCL(i), serrbnd/RCR(i)            
         enddo
         close(61)
      endif
      
      return
      end
!----------------------------------------------------------------------
!      
!     LAPACK interface for the non-symmetric eigenvalue solver.
!     Upon finishing, RITZR and RITZI contain the real and imaginary part,
!     of the computed eigenvalues. Complex conjugate pairs of the
!     eigenvalues appear with the eigenvalue having the positive
!     imaginary part first.
!     The corresponding eigenvectors are stored in EVEC. If the j:th and
!     (j+1):th eigenvalue form a complex conjugate pair, then:
!     v(j) = EVEC(:,j)+i*EVEC(:,j+1), v(j+1) = EVEC(:,j)-i*EVEC(:,j+1)
!
!     Dongarra et al. (1999)
!
      subroutine eig_wrapper(n,kind)

      implicit none

      include 'SIZE'
      include 'OTD'
      include 'WLAPACK'
      
      character*1 jobvl, jobvr, kind
      integer lda, ldvl, ldvr, info,
     $     n,
     $     i, i0

!     Input parameter 'kind' determines whether left and/or right
!     eigenvectors should be computed
!      
!     Define LAPACK-variables
      if (kind.eq.'r') then
         jobvl = 'N'            ! don't compute left eigenvectors
         jobvr = 'V'            ! compute right eigenvectors
      elseif (kind.eq.'l') then
         jobvl = 'V'            ! compute left eigenvectors
         jobvr = 'N'            ! don't compute right eigenvectors
      elseif (kind.eq.'b') then
         jobvl = 'V'            ! compute left eigenvectors
         jobvr = 'V'            ! compute right eigenvectors
      else
         if (nid.eq.0)
     $        write(6,*) 'ERROR: choose left/right/both eigenvectors',
     $        kind
      endif
      lda   = npert              ! leading dimension of LR
      ldvl  = npert              ! leading dimension of EVECL
      ldvr  = npert              ! leading dimension of EVECR
!     Compute the eigenvalues/-vectors in double precision
      call dgeev(jobvl,jobvr,n,LR,lda,EIGR,EIGI,EVL,ldvl,
     $        EVR,ldvr,RWORK,LWORKR,info)
         
!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) then
            write(6,*) 'ERROR: the QR algorithm failed.', info
            write(6,*) '         Converged eigenvalues:'
            i0 = info+1
            do i=i0,n
               write(6,*) EIGR(i), EIGI(i)
            enddo
         endif
         
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGEEV: successful exit!'
         if (nid.eq.0) write(6,*) '        Optimal LWORKR=',
     $        int(RWORK(1)), LWORKR
      endif
      
      return
      end
!----------------------------------------------------------------------      
!                                              
!     LAPACK interface for determining the reciprocal condition number
!     for a square (nflds x nflds) input matrix. The choice of
!     condition number is determined by norm.
!
!     Note: nflds <= LPERT
!      
      subroutine cond_wrapper(rcond,anorm,nflds,norm)

      implicit none

      include 'SIZE'
      include 'OTD'
      include 'WLAPACK'

      character*1 norm    ! Specifies whether the 1-norm condition number
                          ! or the infinity-norm condition number is requested:
                          ! = '1' or 'O': 1-norm
                          ! = 'I'       : Infinity-norm

      integer info        ! = 0:  successful exit.
                          ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                          ! > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                          !       has been completed, but the factor U is exactly singular,
                          !       and division by zero will occur if it is used to solve a
                          !       system of equations. (For DGETRF only)
      
      integer m, n, lda,
     $     nflds
      real anorm, rcond
      real dlange       ! function
      
      if (nflds.ge.LPERT) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: nflds>=LPERT. Increase LPERT.', nflds, LPERT
         call exitt
      endif
      
!     Define LAPACK-variables
      m     = nflds             ! no rows in matrix
      n     = nflds             ! no columns in matrix/the order of the matrix
      lda   = LPERT              ! leading dimension of the matrix

!     Compute the selected norm of the input matrix in double precision
      anorm = dlange(norm,m,n,LR,lda,RWORK)
         
      call copy(ALU,LR,LPERT*LPERT)
!     LU factorize the matrix A in double precision
!     NOTE: Elements 1:min(m,n) of IWORK contain the pivot indices
      call dgetrf(m,n,ALU,lda,IWORK,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) write(6,*) 'WARNING: U(i,i) is exactly zero.',
     $        info
      else
         if (nid.eq.0) write(6,*) 'DGETRF: successful exit!'
      endif
      
!     Compute the resiprocal of the condition number in double precision
      call dgecon(norm,n,ALU,lda,anorm,rcond,RWORK,IWORK,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGECON: successful exit!'
      endif

      return
      end
!----------------------------------------------------------------------

!=======================================================================
*
*     Auxiliary routine: printing eigenvalues.
*
      SUBROUTINE PRINT_EIGENVALUES( DESC, N, WR, WI )
      CHARACTER*(*)    DESC
      INTEGER          N
      REAL             WR( * ), WI( * )
*
      REAL             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO J = 1, N
         IF( WI( J ).EQ.ZERO ) THEN
            WRITE(*,9998,ADVANCE='NO') WR( J )
         ELSE
            WRITE(*,9999,ADVANCE='NO') WR( J ), WI( J )
         END IF
      END DO
      WRITE(*,*)
*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
*
*     Auxiliary routine: printing eigenvectors.
*
      SUBROUTINE PRINT_EIGENVECTORS( DESC, N, WI, V, LDV )
      CHARACTER*(*)    DESC
      INTEGER          N, LDV
      REAL             WI( * ), V( LDV, * )
*
      REAL             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, N
         J = 1
         DO WHILE( J.LE.N )
            IF( WI( J ).EQ.ZERO ) THEN
               WRITE(*,9998,ADVANCE='NO') V( I, J )
               J = J + 1
            ELSE
               WRITE(*,9999,ADVANCE='NO') V( I, J ), V( I, J+1 )
               WRITE(*,9999,ADVANCE='NO') V( I, J ), -V( I, J+1 )
               J = J + 2
            END IF
         END DO
         WRITE(*,*)
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END 
