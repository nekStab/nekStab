c-----------------------------------------------------------------------
      subroutine Schur_decomp(A,VS,VP,n,sdim)
      implicit none

c     ----- Required variables for dgees -----

      character*1           :: JOBVS = 'V', SORT = 'S'
      integer               :: n, LDA, SDIM, LDVS, LWORK, INFO
      real, dimension(n,n)  :: A, VS
      real, dimension(n)    :: WR, WI
      real, dimension(3*n)  :: WORK
      logical, dimension(n) :: BWORK
      complex*16, dimension(n) :: VP

      external select_eigenvalue

c     ----- Schur decomposition -----

      LDA   = max(1,n)
      LDVS  = max(1,n)
      LWORK = max(1,3*n)

      call dgees(JOBVS, SORT, SELECT_EIGENVALUE, N, A, LDA
     $     , SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO)

      VP = WR*(1.0D0,0.0D0) + WI*(0.0D0,1.0D0)

      return
      end
c-----------------------------------------------------------------------
      /* subroutine compute_eigenvec_schur(A,FP,n)

      integer :: n
      complex*16, dimension(n,n) :: FP

c     ----- Required variables for ztrevc -----

      character*1 :: JOB = 'R', HOWMY = 'A'
      real, dimension(n,n) :: A
      real, dimension(n,n) :: T, VR, VL
      integer :: ldt, ldvl, ldvr, mm, m
      real, dimension(3*n) :: WORK
      integer :: INFO, i
      real :: norme

c     ----- Compute the eigenvectors of T -----

      T    = A
      ldt  = max(1,n)
      ldvl = max(1,n)
      ldvr = max(1,n)
      m    = n
      mm   = m

      call dtrevc( JOB, HOWMY, select_eigenvalue, N, T, LDT, VL,
     $     LDVL, VR, LDVR, MM, M, WORK, INFO )

      i = 1
      do while ( i.LT.n )

         if ( abs(T(i+1,i)).LT.1e-8 ) then
            FP(:,i) = VR(:,i)*(1.0D0,0.0D0)
            norme   = sqrt( sum( VR(:,i)**2. ) )
            FP(:,i) = FP(:,i)/norme
            i = i + 1
         else
            FP(:,i)   = VR(:,i)*(1.0D0,0.0D0) + VR(:,i+1)*(0.0D0,1.0D0)
            FP(:,i+1) = VR(:,i)*(1.0D0,0.0D0) - VR(:,i+1)*(0.0D0,1.0D0)
            norme     = sqrt( sum( VR(:,i)**2 + VR(:,i+1)**2) )
            FP(:,i)   = FP(:,i)/norme
            FP(:,i+1) = FP(:,i+1)/norme
            i = i + 2
         endif

      enddo

      return
      end */
c-----------------------------------------------------------------------
      subroutine eigen_decomp(A,n,VP,VR)
      implicit none
c     ----- Required variables for zgeev (eigenpairs computations) -----
      integer    :: INFO
      integer    :: n, LWORK
      real       :: A(n,n)
      real       :: RWORK(2*n)
      complex*16 :: VL(1,n)
      complex*16 :: VR(n,n)
      complex*16 :: WORK(4*n)
      complex*16 :: VP(n)
c     ----- Computation of eig(A) -----
      LWORK = 4*n
      call ZGEEV('N','V',n,A*(1.0,0.0),n,VP,VL,1,VR,n,
     $     WORK,LWORK,RWORK,INFO)
c     ----- Sort the eigedecomposition -----
      call sort_eigendecomp(VP,VR,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine matrix_matrix(A,B,nra,nc,ncb)
      implicit none
      integer :: nra, nc, ncb
c     ----- Required variables for dgeem -----
      character :: TRANSA = 'N', TRANSB = 'N'
      integer   :: m, n, k
      integer   :: lda, ldb, ldc
      real, dimension(nra,nc)  :: A
      real, dimension(nc,ncb)  :: B
      real, dimension(nra,ncb) :: C
      real :: alpha, beta
c     ----- Matrix-matrix multiplication -----
      m = nra;      n = ncb;      k = nc
      lda = max(1,m); ldb = max(1,k); ldc = max(1,m)
      alpha = 1.0D0; beta  = 0.0D0
      call dgemm(TRANSA,TRANSB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
      A(:,1:ncb) = C
      return
      end subroutine
c-----------------------------------------------------------------------
