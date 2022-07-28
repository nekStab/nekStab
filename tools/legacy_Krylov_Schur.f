      subroutine Krylov_Schur

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'

c     =================================================
c     =====                                       =====
c     ===== Declaration of the required variables =====
c     =====                                       =====
c     =================================================

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

c     ----- Krylov basis V for the projection M*V = V*H -----

      real, dimension(lt,k_dim+1)        :: V_x, V_y, V_z
      real, dimension(lt2,k_dim+1)       :: Pressure

c     ----- Upper Hessenberg matrix -----

      real, dimension(k_dim,k_dim)       :: H_mat
      real, dimension(1,k_dim)           :: b_vec

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: VP
      complex*16, dimension(k_dim,k_dim) :: FP

c     ----- Miscellaneous -----

      integer :: mstart

      real, dimension(k_dim)             :: residu
      real                               :: alpha, beta
      logical                            :: is_converged
      integer                            :: have_converged

c     =========================================
c     =====                               =====
c     ===== Initialization of most arrays =====
c     =====                               =====
c     =========================================

      n    = nx1*ny1*nz1*nelt
      n2   = nx2*ny2*nz2*nelt

      time = 0.0

      V_x = 0.0D+00
      V_y = 0.0D+00

      H_mat  = 0.0D+00
      b_vec  = 0.0D+00
      residu = 0.0D+00

      if ( if3D ) then
         V_z = 0.0D+00
      endif

      if ( ifpo ) then 
         Pressure = 0.0D+00
      endif

c     ==========================
c     =====                =====
c     ===== PRE-PROCESSING =====
c     =====                =====
c     ==========================

c     ----- Load base flow -----

      if(ifldbf)then
       if(nid.eq.0)write(*,*)'Loading base flow: ',baseflow_handle
        call load_fld(baseflow_handle)
      endif
 
      ifto = .true. !to output scalar field t
      call comp_vort3(t(1,1,1,1,2),wo1,wo2,vx,vy,vz)
      call lambda2(t(1,1,1,1,3))
      call outpost(vx,vy,vz,pr,t,'BFL')

      call opcopy(ubase,vbase,wbase,vx,vy,vz)

c     ----- Creates seed vector for the Krylov subspace -----

      if(uparam(1).gt.0.0)then

      !noise will be generated in the vxp,vyp,vzp in useric
      call outpost(vxp,vyp,vzp,prp,t,'IV_')

      call dsavg(vxp)
      call dsavg(vyp)
      if ( if3D ) call dsavg(vzp)
      if ( ifpo ) call dsavg(prp)
      call bcdirvc(vxp,vyp,vzp,v1mask,v2mask,v3mask)

      alpha = glsc3(vxp,BM1,vxp,n) + glsc3(vyp,BM1,vyp,n)
      if( if3D ) alpha = alpha + glsc3(vzp,BM1,vzp,n)
      alpha = 1.d0/sqrt(alpha)
      call opcmult(vxp,vyp,vzp,alpha)
      if ( ifpo ) call cmult(prp,alpha,n2)

      call outpost(vxp,vyp,vzp,prp,t,'IV_') !initial condition VEC t=0

      DO ISTEP = 1, NSTEPS
         call opcopy(vx,vy,vz,ubase,vbase,wbase)
         call nek_advance
         !call userchk !compute buffer zone forces
      ENDDO

      else
       call opcopy(vxp,vyp,vzp,ubase,vbase,wbase)
      endif

c     ----- Normalized to unit-norm -----
      
      alpha = glsc3(vxp,BM1,vxp,n)
     $     +  glsc3(vyp,BM1,vyp,n)

      if( if3D ) alpha = alpha + glsc3(vzp,BM1,vzp,n)
      
      alpha = 1.d0/sqrt(alpha)
      call opcmult(vxp,vyp,vzp,alpha)
      if ( ifpo ) call cmult(prp,alpha,n2)

c     ======================================
c     =====                            =====
c     ===== Krylov-Schur decomposition =====
c     =====                            =====
c     ======================================

c     ----- Storing the starting vector in the Krylov basis -----

      call opcopy(v_x(:,1),v_y(:,1),v_z(:,1),vxp,vyp,vzp)
      if ( ifpo ) call copy(Pressure(:,1),prp(:,1),n2)

      mstart = 1
      is_converged = .false.
      
      do while ( .not. is_converged )

c     ----- m-step Arnoldi factorization -----

         call Arnoldi_factorization( V_x , V_y , V_z , Pressure
     $        , H_mat , beta , mstart )

c     ----- Check the convergence of the Ritz eigenpairs -----

         call check_convergence( residu , have_converged , beta
     $        , tolerance , H_mat , k_dim )

         if ( nid.EQ.0 ) then
            write(*,*) 'Number of converged eigenvalues : ' 
     $           , have_converged
         endif

         if ( ifschur .EQV. .false. ) have_converged = k_dim        
 
         if ( have_converged .GE. wanted_converged ) then

            is_converged = .true.

         else

c     ----- Performs Schur factorization of the Krylov basis -----

            call Schur_factorization( mstart , H_mat 
     $           , V_x , V_y , V_z , Pressure , residu , beta)
         endif

      enddo

c     ----- Final computation of the eigenpairs of the Hessenberg matrix -----

      call eigendecomposition( H_mat , k_dim , VP , FP )
      
c     ----- Output all the spectrums and converged eigenmodes -----

      call output_KrylovSchur( VP , FP , 
     $     V_x , V_y , V_z , Pressure , residu )

      call nek_end
      return
      end

c--------------------------------------------------------------------------------





      subroutine Arnoldi_factorization(Qx,Qy,Qz,Qp,H,beta,mstart)
      
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'

      real*4 papi_mflops
      integer*8 papi_flops

c     =================================================
c     =====                                       =====
c     ===== Declaration of the required variables =====
c     =====                                       =====
c     =================================================

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt
cc      integer, parameter                 :: k_dim_tot = k_dim+1
      
c     ----- Miscellaneous -----

      real                               :: alpha, beta
      integer                            :: mstep, mstart

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: VP
      complex*16, dimension(k_dim,k_dim) :: FP

c     ----- Krylov basis V for the projection MV = VH -----

      real, dimension(lt,k_dim+1)        :: Qx, Qy, Qz
      real, dimension(lt2,k_dim+1)       :: Qp

c     ----- Orthogonal residual f = w - (V,w)*V -----

      real, dimension(lt)                :: F_x, F_y, F_z
      real, dimension(lt2)               :: F_p

c     ----- Upper Hessenberg matrix -----

      real, dimension(k_dim,k_dim)       :: H
      real, dimension(k_dim)             :: h_vec

c     ----- Working arrays -----

      real, dimension(lt)                :: work1, work2, work3
      real, dimension(lt2)               :: workp

c     =========================================
c     =====                               =====
c     ===== Initialization of most arrays =====
c     =====                               =====
c     =========================================

      n  = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

      F_x = 0.0D+00
      F_y = 0.0D+00

      h_vec = 0.0D+00
      
      if ( if3D ) then
         F_z = 0.0D+00
      endif

      if ( ifpo ) then
         W_p = 0.0D+00
         F_p = 0.0D+00
      endif

c     ----- Arnoldi factorization -----
         
      do mstep = mstart, k_dim
         
         if (nid.EQ.0) then
            write(6,*) 'Begin iteration :', mstep , '/' , k_dim
         endif

         call opcopy(vxp,vyp,vzp,
     $        Qx(:,mstep),Qy(:,mstep),Qz(:,mstep))
         if ( ifpo ) call copy(prp,Qp(:,mstep),n2)
         
c     ----- Matrix-vector product w = M*v -----

         call matrix_vector_product(f_x, f_y, f_z, f_p,
     $        qx(:,mstep), qy(:,mstep), qz(:,mstep))

c     ----- Compute the orthogonal residual f with respect to all previous Krylov vectors -----

         h_vec = 0.0D+00

         do i = 1, mstep

            call opcopy(work1,work2,work3
     $           ,Qx(:,i),Qy(:,i),Qz(:,i))
            if ( ifpo ) call copy(workp,Qp(:,i),n2)
            
            h_vec(i) = glsc3(work1,BM1,f_x,n)
     $           +     glsc3(work2,BM1,f_y,n)
            if( if3D ) h_vec(i) = h_vec(i) + glsc3(work3,BM1,f_z,n)

            call opcmult(work1,work2,work3,h_vec(i))
            if ( ifpo ) call cmult(workp,h_vec(i),n2)

            call opsub2(f_x,f_y,f_z,work1,work2,work3)
            if ( ifpo ) call sub2(f_p,workp,n2)

         enddo

c     ----- Fills in the upper Hessenberg matrix -----

         H(1:mstep,mstep) = h_vec(1:mstep)
         
c     ----- Compute the norm of the orthogonal residual f -----
         
         beta = glsc3(f_x,BM1,f_x,n)
     $        + glsc3(f_y,BM1,f_y,n)
         
         if( if3D ) beta = beta + glsc3(f_z,BM1,f_z,n)
         
         beta = sqrt(beta)
                           
         if ( mstep.LT.k_dim ) then

c     ----- Normalizes the orthogonal residual and uses it as a new Krylov vector -----
            
            H(mstep+1,mstep) = beta
                        
         endif

         beta = 1.d0/beta

         call opcmult(f_x,f_y,f_z,beta)
         if ( ifpo ) call cmult(f_p,beta,n2)
         
         call opcopy(Qx(:,mstep+1),Qy(:,mstep+1),Qz(:,mstep+1)
     $        ,f_x,f_y,f_z)
         if ( ifpo ) call copy(Qp(:,mstep+1),f_p,n2)

         call outpost(f_x,f_y,f_z,f_p,t,'AR_')

         if(nid.EQ.0) write(6,*) "End of iteration:", mstep

      enddo

      return
      end subroutine





c----------------------------------------------------------------------





      subroutine matrix_matrix(A,B,nra,nc,ncb)

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

      m = nra
      n = ncb
      k = nc

      lda = max(1,m)
      ldb = max(1,k)
      ldc = max(1,m)

      alpha = 1.0D0
      beta  = 0.0D0
    
      call dgemm(TRANSA,TRANSB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)

      A(:,1:ncb) = C

      return
      end subroutine





c----------------------------------------------------------------------





      subroutine Schur_decomposition(A,VS,VP,n,sdim)

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

      VP = WR*(1.0,0.0) + WI*(0.0,1.0)

      return
      end subroutine

      subroutine compute_eigenvec_schur(A,FP,n)

      complex*16, dimension(n,n) :: FP

c     ----- Required variables for ztrevc -----

      character*1 :: JOB = 'R', HOWMY = 'A'
      integer :: n
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
            FP(:,i) = VR(:,i)*(1.0,0.0)

            norme   = sqrt( sum( VR(:,i)**2. ) )

            FP(:,i) = FP(:,i)/norme

            i = i + 1
         else
            FP(:,i)   = VR(:,i)*(1.0,0.0) + VR(:,i+1)*(0.0,1.0)
            FP(:,i+1) = VR(:,i)*(1.0,0.0) - VR(:,i+1)*(0.0,1.0)

            norme     = sqrt( sum( VR(:,i)**2. + VR(:,i+1)**2.) )
            FP(:,i)   = FP(:,i)/norme
            FP(:,i+1) = FP(:,i+1)/norme

            i = i + 2
         endif

      enddo
      
      return
      end

      function select_eigenvalue(wr,wi)

      logical :: select_eigenvalue
      real :: wr, wi
      real :: delta = 0.1D0
      
      if ( sqrt(wr**2. + wi**2.) .GT. 1.d0-delta ) then
         select_eigenvalue = .true.
      else
         select_eigenvalue = .false.
      endif

      end





c---------------------------------------------------------------





      subroutine sort_eigendecomp(VP,FP,n)
      
      integer                    :: n
      complex*16, dimension(n)   :: VP
      complex*16, dimension(n,n) :: FP
      real, dimension(n)         :: norm
      real                       :: temp_real
      complex*16                 :: temp_complex
      complex*16, dimension(n)   :: temp_n
      integer                    :: i, j, k, l
      
c     ----- Sorting the eigenvalues according to their norm -----
      
      temp_n   = (0.0,0.0)
      norm     = 0.0

      norm = sqrt(real(VP)**2. + aimag(VP)**2.)
      
      do k = 1,n-1
         do l = k+1,n
            
            if (norm(k).LT.norm(l)) then
               
               temp_real    = norm(k)
               temp_complex = VP(k)
               temp_n       = FP(:,k)
               
               norm(k) = norm(l)
               norm(l) = temp_real
               
               VP(k)  = VP(l)
               VP(l)  = temp_complex
               
               FP(:,k) = FP(:,l)
               FP(:,l) = temp_n

            endif
            
         enddo
      enddo
      
      return
      end





c-----------------------------------------------------------------------





      subroutine eigendecomposition(A,n,VP,VR)

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





c----------------------------------------------------------------------





      subroutine matvec(FP_Cx,Qx,FP,n)!,k_dim)
      include 'SIZE'
      include 'TOTAL'
      
      complex*16, dimension(n)       :: FP_Cx
      complex*16, dimension(k_dim)   :: FP
      real      , dimension(n,k_dim) :: Qx
      
      FP_Cx = (0.,0.)
      
      do i = 1,k_dim
         FP_Cx = FP_Cx + Qx(:,i)*FP(i)
      enddo
      
      return
      end





c----------------------------------------------------------------------





      subroutine output_KrylovSchur(VP,FP,Qx,Qy,Qz,Qp,residu)
      
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

c     ----- Krylov basis V for the projection M*V = V*H -----

      real, dimension(lt,k_dim+1)        :: Qx, Qy, Qz
      real, dimension(lt2,k_dim+1)       :: Qp

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: VP
      complex*16, dimension(k_dim,k_dim) :: FP

c     ----- Arrays for the storage/output of a given eigenmode of the NS operator -----

      complex*16, dimension(lt)          :: FP_Cx, FP_Cy, FP_Cz
      complex*16, dimension(lt2)         :: FP_Cp

c     ----- Miscellaneous -----

      integer :: ifich1 = 10, ifich2 = 20, ifich3 = 30
      integer :: mstart

      real                               :: sampling_period
      real, dimension(k_dim)             :: residu
      real                               :: alpha, beta

      sampling_period = dt*nsteps

      n = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nez2*nelt

c     ----- Output all the spectrums and converged eigenmodes -----
      
      if(nid.EQ.0) then
         
         open( unit = ifich1, 
     $        file  = 'Spectre_H.dat',
     $        form  = 'formatted')
         
         open( unit = ifich2,
     $        file  = 'Spectre_NS.dat',
     $        form  = 'formatted')
         
         open( unit = ifich3, 
     $        file  = 'Spectre_NS_conv.dat',
     $        form  = 'formatted')
                  
      endif
      
      do i = 1,k_dim
 
         time = real(i)       
  
         if(nid.EQ.0) then
            
            write(ifich1,*) real(VP(i)),
     $           aimag(VP(i)),
     $           residu(i)
            
            write(ifich2,*) log(abs(VP(i)))/sampling_period,
     $           ATAN2(aimag(VP(i)),real(VP(i)))/sampling_period,
     $           residu(i)
            
         endif
         
         if( residu(i).LT.tolerance ) then
            
            if(nid.EQ.0) then
               
               write(ifich3,*) log(abs(VP(i)))/sampling_period,
     $              ATAN2(aimag(VP(i)),real(VP(i)))/sampling_period
               
            endif
            
c     ----- Computation of the corresponding eigenmode -----
            
            FP_Cx = (0.0,0.0)
            FP_Cy = (0.0,0.0)
            if ( if3D ) FP_Cz = (0.0,0.0)
            if ( ifpo ) FP_Cp = (0.0,0.0)
            
            call matvec(FP_Cx,Qx(:,1:k_dim),FP(:,i),lt)!,k_dim)
            call matvec(FP_Cy,Qy(:,1:k_dim),FP(:,i),lt)!,k_dim)
            if ( if3D ) then
               call matvec(FP_Cz,Qz(:,1:k_dim),FP(:,i),lt)!,k_dim)
            endif
            if ( ifpo ) then
               call matvec(FP_Cp,Qp(:,1:k_dim),FP(:,i),lt2)!,k_dim)
            endif
            
c     ----- Normalization to be unit-norm -----
c     Note: volume integral of FP*conj(FP) = 1.
            
            alpha = glsc3(real(FP_Cx),bm1,real(FP_Cx),n)
     $           +  glsc3(real(FP_Cy),bm1,real(FP_Cy),n)
     $           +  glsc3(aimag(FP_Cx),bm1,aimag(FP_Cx),n)
     $           +  glsc3(aimag(FP_Cy),bm1,aimag(FP_Cy),n)
            
            if ( if3D ) alpha = alpha 
     $           +  glsc3(real(FP_Cz),bm1,real(FP_Cz),n)
     $           +  glsc3(aimag(FP_Cz),bm1,aimag(FP_Cz),n)
            
            alpha = 1./sqrt(alpha)
            
c     ----- Output the imaginary part -----
            
            call opcopy(vx,vy,vz
     $           ,aimag(FP_Cx)
     $           ,aimag(FP_Cy)
     $           ,aimag(FP_Cz))
            if ( ifpo ) call copy(pr,aimag(FP_Cp),n2)
            
            call opcmult(vx,vy,vz,alpha)
            if ( ifpo ) call cmult(pr,alpha,n2)
            call outpost(vx,vy,vz,pr,t,'Im_')
            
c     ----- Output the real part -----
            
            call opcopy(vx,vy,vz
     $           ,real(FP_Cx)
     $           ,real(FP_Cy)
     $           ,real(FP_Cz))
            if ( ifpo ) call copy(pr,real(FP_Cp),n2)
            
            call opcmult(vx,vy,vz,alpha)
            if ( ifpo ) call cmult(pr,alpha,n2)
            call outpost(vx,vy,vz,pr,t,'Re_')
            
         endif
         
      enddo

      if ( nid.EQ.0 ) then
         
         close(ifich1)
         close(ifich2)
         close(ifich3)

      endif

      return
      end





c----------------------------------------------------------------------





      subroutine check_convergence( residual , have_converged
     $     , beta , tolerance , H , n )

      !----- Inputs -----!

      integer              :: n
      real, dimension(n,n) :: H
      real                 :: tolerance, beta

      !----- Outputs -----!

      real, dimension(n) :: residual
      integer            :: have_converged

      !----- Miscellaneous -----!

      complex*16, dimension(n)   :: VP
      complex*16, dimension(n,n) :: FP

      integer :: i

      !----- Compute the eigedecomposition of H -----!

      VP = ( 0.0D+00 , 0.0D+00 )
      FP = ( 0.0D+00 , 0.0D+00 )

      call eigendecomposition(H,n,VP,FP)
      
      !----- Compute the number of converged Ritz pairs -----!

      have_converged = 0
      residual      = 1.0D+00

      do i = 1,n
         residual(i) = abs( beta * FP(n,i) )
         if ( residual(i) .LT. tolerance ) then
            have_converged = have_converged + 1
         endif
      enddo

      return
      end




c----------------------------------------------------------------------





      subroutine Schur_factorization( mstart , H_mat 
     $     , V_x , V_y , V_z , Pressure , residu , beta)

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
cc      include 'LINEAR_STABILITY'

c     =================================================
c     =====                                       =====
c     ===== Declaration of the required variables =====
c     =====                                       =====
c     =================================================

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

c     ----- Krylov basis V for the projection M*V = V*H -----

      real, dimension(lt,k_dim+1)        :: V_x, V_y, V_z
      real, dimension(lt2,k_dim+1)       :: Pressure

c     ----- Upper Hessenberg matrix -----

      real, dimension(k_dim,k_dim)       :: H_mat
      real, dimension(1,k_dim)           :: b_vec

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: VP
      complex*16, dimension(k_dim,k_dim) :: FP

c     ----- Miscellaneous -----

      integer :: mstart
      real, dimension(k_dim)             :: residu

c     ----- Schur and Hessenberg decomposition -----

      real, dimension(k_dim,k_dim)       :: Q1

c     =========================================
c     =====                               =====
c     ===== Initialization of most arrays =====
c     =====                               =====
c     =========================================

      n  = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

      Q1 = 0.0D+00

      b_vec          = 0.0D+00
      b_vec(1,k_dim) = beta

c     ----- Compute the Schur decomposition -----
            
      VP = ( 0.0D+00 , 0.0D+00 )
      Q1 = 0.0D+00
      
      call Schur_decomposition(H_mat,Q1,VP,k_dim,mstart)
      
      H_mat(1:mstart,mstart+1:k_dim)       = 0.0D+00
      H_mat(mstart+1:k_dim,1:mstart)       = 0.0D+00
      H_mat(mstart+1:k_dim,mstart+1:k_dim) = 0.0D+00
      
      if ( mstart.EQ.0 ) then
         
         mstart = 1
         H_mat  = 0.0D+00
         
         call opcopy( V_x(:,mstart), V_y(:,mstart), V_z(:,mstart),
     $        V_x(:,k_dim+1), V_y(:,k_dim+1), V_z(:,k_dim+1) )
         
         if ( ifpo ) call copy( Pressure(:,mstart) ,
     $        Pressure(:,k_dim+1), n2)
         
      else
         
         if ( nid.EQ.0 ) then
            write(*,*) mstart,
     $           'Ritz eigenpairs have been selected'
            write(*,*)

            do i = 1, k_dim
               write(*,*) 'Residual of eigenvalue', 
     $              i, ' : ', residu(i)
            enddo
            write(*,*)
         endif
         
c     ----- Re-order the Krylov basis -----
         
         call matrix_matrix(
     $        V_x(:,1:k_dim),
     $        Q1,
     $        lt,
     $        k_dim,
     $        k_dim)
         
         call matrix_matrix(
     $        V_y(:,1:k_dim),
     $        Q1,
     $        lt,
     $        k_dim,
     $        k_dim)
         
         if ( if3D ) call matrix_matrix(
     $        V_z(:,1:k_dim),
     $        Q1,
     $        lt,
     $        k_dim,
     $        k_dim)
         
         if ( ifpo ) call matrix_matrix(
     $        Pressure(:,1:k_dim),
     $        Q1,
     $        lt2,
     $        k_dim,
     $        k_dim)
         
c     ----- Update the matrix with b(prime)*Q -----
         
         call matrix_matrix(b_vec,Q1,1,k_dim,k_dim)
         
         H_mat(mstart+1,1:mstart) = b_vec(1,1:mstart)
         mstart = mstart + 1
         
         call opcopy(V_x(:,mstart), V_y(:,mstart), V_z(:,mstart),
     $        V_x(:,k_dim+1), V_y(:,k_dim+1), V_z(:,k_dim+1))
         
      endif
      
      return
      end




c----------------------------------------------------------------------




      subroutine matrix_vector_product(fx, fy, fz, fp, qx, qy, qz)

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'ADJOINT'

c     =================================================
c     =====                                       =====
c     ===== Declaration of the required variables =====
c     =====                                       =====
c     =================================================

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

!     ----- Working arrays -----

      real, dimension(lt)                :: qx, qy, qz
      real, dimension(lt)                :: fx, fy, fz
      real, dimension(lt2)               :: fp

!     ====================================================
!     =====                                          =====
!     ===== Computation of the matrix-vector product =====
!     =====                                          =====
!     ====================================================

      n = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

!     ----- Initial condition -----

      call opcopy(vxp, vyp, vzp, qx, qy, qz)
      call opcopy(vx, vy, vz, ubase, vbase, wbase)

!     ----- Time-stepper matrix-vector product -----

      do istep = 1, nsteps
         call nek_advance
c         call userchk
      enddo

!     ----- Save the result into the output arrays -----
      
      call opcopy(fx, fy, fz, vxp, vyp, vzp)
      if (ifpo) call copy(fp, prp, n2)
      
      return
      end
