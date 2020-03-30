!!!!!!!!!
!!! norm           -> return |u|
!!! orthonormalize -> return u/|u|
!!! switch_to_lnse_steady -> move somewhere else

!!! krylov_schur
!!! Arnoldi_fact
!!! check_conv
!!! Schur_decomp
!!! Schur_fact
!!! matrix_vector_product (fx, fy, fz, fp, qx, qy, qz)
!!! outpost_ks(VP,FP,Qx,Qy,Qz,Qp,residu)

!!! select_eigenvalue
!!! matrix_matrix
!!! eigen_decomp
!!! sort_eigendecomp
!!! matvec

c----------------------------------------------------------------------
      subroutine norm(qx, qy, qz, alpha) ! Compute vector norm
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      real, intent(in), dimension(lt)  :: qx, qy, qz
      real, intent(out)                              :: alpha
      real                                           :: glsc3, alpha_tmp
      integer n
      n      = nx1*ny1*nz1*nelt

      alpha_tmp = 0.0d0 ; n = nx1*ny1*nz1*nelt
      alpha_tmp = glsc3(qx,bm1,qx,n) + glsc3(qy,bm1,qy,n)
      if( if3d ) alpha_tmp = alpha_tmp + glsc3(qz,bm1,qz,n)

      alpha = sqrt(alpha_tmp)

      end subroutine
c----------------------------------------------------------------------
      subroutine orthonormalize2(qx, qy, qz, qp, beta) ! Normalize to the unit-norm
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt)  :: qx, qy, qz
      real, dimension(lt2)  :: qp
      real                               :: alpha,beta_tmp
      real, intent(out), optional        :: beta
      integer n2
      n2     = nx2*ny2*nz2*nelt

      call norm(qx, qy, qz, alpha)
      beta_tmp = 1.0D0/alpha

      call opcmult(qx, qy, qz, beta_tmp)
      if(ifpo) call cmult(qp,beta_tmp,n2)
      beta = beta_tmp

      end subroutine
c----------------------------------------------------------------------
      subroutine orthonormalize(qx, qy, qz, qp, beta) ! Normalize to the unit-norm
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt)  :: qx, qy, qz
      real, dimension(lt2)  :: qp
      real                               :: alpha,beta,glsc3
      integer n,n2
      n      = nx1*ny1*nz1*nelt
      n2     = nx2*ny2*nz2*nelt

      alpha = 0.0d0
      alpha = glsc3(qx,bm1,qx,n) + glsc3(qy,bm1,qy,n)
      if( if3d ) alpha = alpha + glsc3(qz,bm1,qz,n)
      alpha = sqrt(alpha)

      beta = 1.0D0/alpha
      call opcmult(qx, qy, qz, beta)
      if(ifpo) call cmult(qp,beta,n2)

      end
c-----------------------------------------------------------------------
      subroutine switch_to_lnse_steady(steps_in)
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'OPCTR'
      include 'CTIMER'

      integer, intent(in), optional :: steps_in
      integer steps
      real sampling_period

      if(.not.present(steps_in))then
        steps = steps_in
      else
        if(nid.eq.0)write(6,*)' nsteps not specified: computing with uparam(04)'
        sampling_period = real(1.0d0/uparam(04))
        steps = int(sampling_period/dt)
      endif

      if(nid.eq.0)then
       write(6,*)' Switching to linearized solver...'
       write(6,*)' Base flow stability: nsteps=',steps
       write(6,*)
      endif

      time = 0
      istep = 0
      nsteps = steps
      param(10) = 0     !endTime
      param(11) = steps !numSteps to overwrite final time
      param(14) = 999. !write interval runtime
      param(12) = -1*abs(param(12)) !freeze dt
      param(31) = 1 !numberOfPerturbations
      param(63) = 1 ! Enforce 64-bit output
      npert = param(31)

      call bcast(param,200*wdsize) !broadcast all parameters to processors

      ifpert = .true.
      call bcast(ifpert, lsize) !start linearized

      ifbase = .false.
      call bcast(ifbase, lsize) !stop dns

      call chkParam !sanity check

      return
      end
c----------------------------------------------------------------------
      subroutine krylov_schur
      implicit none
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

      real, dimension(lt,k_dim+1)        :: V_x, V_y, V_z
      real, dimension(lt2,k_dim+1)       :: Pressure

c     ----- Upper Hessenberg matrix -----

      real, dimension(k_dim,k_dim)       :: H_mat
      real, dimension(1,k_dim)           :: b_vec

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: VP
      complex*16, dimension(k_dim,k_dim) :: FP

c     ----- Miscellaneous -----

      real wo1(lt),wo2(lt),wo3(lt)
      common /ugrad/ wo1,wo2,wo3

      integer :: mstart
      real, dimension(k_dim)             :: residu
      real                               :: alpha, beta, dumm, glsc3
      logical                            :: is_converged
      integer                            :: have_converged, n, n2, i, j
      integer                            :: ifich3 = 30
      character(len=30)                  :: filename
      n      = nx1*ny1*nz1*nelt
      n2     = nx2*ny2*nz2*nelt
      time   = 0.0D0
      H_mat  = 0.0D0
      b_vec  = 0.0D0
      residu = 0.0D0

      call oprzero(wo1,wo2,wo3)
      call oprzero(V_x,V_y,V_z)
      if(ifpo)Pressure = 0.0D0

      !if(ifldbf)then
         if(.not. ifbfcv)then !skip loading if single run
          if(nid.eq.0)write(*,*)'Loading base flow: ',bf_handle
          call load_fld(bf_handle)
         endif
         ifto = .true. ;ifpo = .true.
         call comp_vort3(t(1,1,1,1,2),wo1,wo2,vx,vy,vz)
         call lambda2(t(1,1,1,1,3))
         call outpost(vx,vy,vz,pr,t,'BF_') !outpost for sanity check
         call opcopy(ubase,vbase,wbase,vx,vy,vz)
       !endif!ifldbf

c     ----- Creates seed vector for the Krylov subspace -----

        if(ifnois)then

         call add_noise(vxp,vyp,vzp)
         call orthonormalize(vxp,vyp,vzp,prp,alpha)
         call matrix_vector_product(vxp,vyp,vzp,prp, vxp,vyp,vzp)

	      else !ifnois.eq..false.

          call opcopy(vxp,vyp,vzp,ubase,vbase,wbase)

        endif

c     ----- Normalized to unit-norm -----

      call orthonormalize(vxp,vyp,vzp,prp,alpha)

c     ======================================
c     =====                            =====
c     ===== Krylov-Schur decomposition =====
c     =====                            =====
c     ======================================

c     ----- Storing the starting vector in the Krylov basis -----


      if(nid.eq.0)write(6,*)'Restarting from:',uparam(2)


      if (uparam(2).gt.0) then

        mstart = uparam(2)

         if(nid.eq.0)then
          write(6,'(a,a,i4.4)')' Loading Hessenberg matrix: HES',trim(SESSION),mstart
          write(filename,'(a,a,i4.4)')'HES',trim(SESSION),mstart
          open(67,file=trim(filename),
     $            status='unknown',form='formatted')
          read(67,*) beta
          write(6,*)'Loading beta',beta

          if (k_dim.lt.mstart) then !subsampling
            do i = 1,k_dim
              do j = 1,mstart
                  if (j.le.k_dim) then
                     read(67,"(1E15.7)") H_mat(i,j)
                  else
                     read(67,"(1E15.7)") dumm
                  endif
               enddo
            enddo
          else
            do i = 1,mstart+1
              do j = 1,mstart
                 read(67,*) H_mat(i,j)
              enddo
            enddo
          endif

          close(67)
         endif!nid.eq.0

         !if(nid.eq.0)then !-> debug only
         ! write(filename,'(a,a,i4.4)')'HESloaded',trim(SESSION),mstart
         ! write(6,*) ''
         ! open(67,file=trim(filename),status='unknown',form='formatted')
         ! write(67,*) beta
         ! do i = 1,mstart+1
         !   do j = 1,mstart
         !      write(67,*) H_mat(i,j)
         !   enddo
         ! enddo
         ! close(67)
         !endif

         mstart=mstart+1  !careful here!
         call load_files(V_x, V_y, V_z, Pressure, mstart, k_dim+1, 'KRY')
         !k_dim+1 is the dim of V_x

      else !uparam(2).eq.0


        if(nid.eq.0)write(6,*)'Starting new Arnoldi decomposition...'

         mstart = 1; istep = 1; time = 0.0D0
         call opcopy(v_x(:,1),v_y(:,1),v_z(:,1),vxp,vyp,vzp)
         if(ifpo) call copy(Pressure(:,1),prp(:,1),n2)
         call whereyouwant('KRY',1)
         if(nid.eq.0)write(6,*)'Sanity check: ensure line 1094 of prepost changed to common!'
         ifto=.false.
         call outpost(v_x(:,1),v_y(:,1),v_z(:,1), Pressure(:,1),t,'KRY')

      endif

      is_converged = .false.
      do while ( .not. is_converged )

c     ----- m-step Arnoldi factorization -----


         call Arnoldi_fact(V_x,V_y,V_z,Pressure,H_mat,beta,mstart)


c     ----- Check the convergence of the Ritz eigenpairs -----

         call check_conv(residu,have_converged,beta,eigen_tol,H_mat,k_dim)

         if(nid.eq.0)write(6,*)'Number of converged eigenvalues: ' ,have_converged

         if( schur_tgt .le. 1 ) then

          have_converged = k_dim; is_converged = .true.
          if(nid.eq.0)write(6,*)'Schur step OFF. Outposting the converged eigenvalues: ',have_converged

         else
                if ( have_converged .GE. schur_tgt ) then
                    is_converged = .true.
                else
                    if(nid.eq.0)write(6,*)'Starting Schur factorization...'
                    call Schur_fact(mstart,H_mat,V_x,V_y,V_z,Pressure,residu,beta)
                endif

         endif

      enddo

c     ----- Final computation of the eigenpairs of the Hessenberg matrix -----

      if(nid.eq.0)write(6,*)'Final eigendecomposition...'
      call eigen_decomp( H_mat , k_dim , VP , FP )

c     ----- Output all the spectrums and converged eigenmodes -----

      if(nid.eq.0)write(6,*)'Exporting modes...'
      call outpost_ks( VP , FP , V_x , V_y , V_z , Pressure , residu )

      if(nid.eq.0)write(6,*)'Stopping code...'
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine Arnoldi_fact(Qx,Qy,Qz,Qp,H,beta,mstart)
      implicit none
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

c     ----- Miscellaneous -----

      real                               :: glsc3,beta,resd
      integer                            :: mstep, mstart, n, n2, i, j, cvrgd

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: VP
      complex*16, dimension(k_dim,k_dim) :: FP

c     ----- Krylov basis V for the projection MV = VH -----

      real, dimension(lt,k_dim+1)        :: Qx, Qy, Qz
      real, dimension(lt2,k_dim+1)       :: Qp

c     ----- Orthogonal residual f = w - (V,w)*V -----

      real, dimension(lt)                :: F_xr, F_yr, F_zr
      real, dimension(lt2)               :: F_pr

c     ----- Upper Hessenberg matrix -----

      real, dimension(k_dim,k_dim)       :: H
      real, dimension(k_dim)             :: h_vec

c     ----- Working arrays -----

      real, dimension(lt)                :: work1, work2, work3
      real, dimension(lt2)               :: workp

      character(len=30) :: filename
      integer :: have_cnv

c     =========================================
c     =====                               =====
c     ===== Initialization of most arrays =====
c     =====                               =====
c     =========================================

      n  = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt
      h_vec(:) = 0.0D0
      call oprzero(F_xr,F_yr,F_zr)
      if(ifpo)F_pr = 0.0D0

c     ----- Arnoldi factorization -----

      do mstep = mstart, k_dim

         if (nid.EQ.0) write(6,*) 'Begin iteration :', mstep , '/' , k_dim

         call opcopy(vxp,vyp,vzp, qx(:,mstep),qy(:,mstep),qz(:,mstep))
         if(ifpo) call copy(prp,qp(:,mstep),n2)

c     ----- Matrix-vector product w = M*v -----

         call matrix_vector_product(f_xr, f_yr, f_zr, f_pr, qx(:,mstep), qy(:,mstep), qz(:,mstep))

c     ----- Compute the orthogonal residual f with respect to all previous Krylov vectors -----

         h_vec(:) = 0.0D0
         do i = 1, mstep

            call opcopy(work1,work2,work3, qx(:,i),qy(:,i),qz(:,i))
            if(ifpo) call copy(workp,      qp(:,i),n2)

            h_vec(i) = glsc3(work1,BM1,f_xr,n) + glsc3(work2,BM1,f_yr,n)
            if( if3D ) h_vec(i) = h_vec(i) + glsc3(work3,BM1,f_zr,n)

            call opcmult(work1,work2,work3, h_vec(i))
            if(ifpo) call cmult(workp,h_vec(i),n2)

            call opsub2(f_xr,f_yr,f_zr, work1,work2,work3)
            if(ifpo) call sub2(f_pr,workp,n2)

         enddo

c     ----- Fills in the upper Hessenberg matrix -----

         H(1:mstep,mstep) = h_vec(1:mstep)

c     ----- Compute the norm of the orthogonal residual f -----

         beta = glsc3(f_xr,BM1,f_xr,n)+glsc3(f_yr,BM1,f_yr,n)
         if(if3D)beta = beta + glsc3(f_zr,BM1,f_zr,n)
         beta = sqrt(beta)

c     ----- Normalizes the orthogonal residual and uses it as a new Krylov vector

         if ( mstep.lt.k_dim ) H(mstep+1,mstep) = beta
         beta = 1.0D0/beta
         call opcmult(f_xr,f_yr,f_zr,beta)
         if(ifpo) call cmult(f_pr,beta,n2)

         call opcopy(Qx(:,mstep+1),Qy(:,mstep+1),Qz(:,mstep+1),f_xr,f_yr,f_zr)
         if(ifpo) call copy(Qp(:,mstep+1),f_pr,n2)

         if(nid.eq.0)write(6,*)
         if(nid.eq.0)write(6,*) ' Outposting Krylov subspace to'

         time=real(time*mstep)
         call whereyouwant('KRY',mstep+1); ifto = .false.
         call outpost(f_xr,f_yr,f_zr,f_pr,t,'KRY')

         VP(1:mstep)=0.0D0; FP(1:mstep,1:mstep)=0.0D0
         call eigen_decomp( H(1:mstep,1:mstep) , mstep , VP(1:mstep) , FP(1:mstep,1:mstep) )

         have_cnv = 0; resd = 0.0D0
         if(nid.eq.0)then ! outpost Hessenberg matrix spectre
          !write(6,'(a,a,i4.4)') ' Ouposting Hessenberg matrix spectre to: HSP',trim(SESSION),mstep
          !write(filename,'(a,a,i4.4)')'HSP',trim(SESSION),mstep

            write(6,'(a,a,i4.4)') ' Ouposting Hessenberg matrix spectre to: H',trim(SESSION),mstep
            write(filename,'(a,i4.4)')'H',mstep

          write(6,*) ''
          open(67,file=trim(filename),status='unknown',form='formatted')
          do i = 1,mstep
          resd = abs(beta*FP(mstep,i))
          if( (resd .LT. eigen_tol)) have_cnv = have_cnv + 1
          !write(67,"(3E15.7)")log(abs(VP(i)))/time,ATAN2(aimag(VP(i)),real(VP(i)))/time,resd
          write(67,"(3E15.7)")real(VP(i)),aimag(VP(i)),resd
          enddo
          close(67)
         endif

         if(nid.eq.0)then
          if(mstep.eq.mstart)then
               write(filename,'(a,i4.4,a,i4.4,a)')'Convegence',mstart,'_',k_dim,'.dat'
          open(99,file=trim(filename),status='unknown',form='formatted')
          write(99,*)'Tolerance:',eigen_tol
          endif
          write(99,'(2i8)')mstep,have_cnv
         endif

         if(nid.eq.0)then ! outpost Hessenberg matrix
          write(6,'(a,a,i4.4)') ' Ouposting Hessenberg matrix to: HES',trim(SESSION),mstep
          write(filename,'(a,a,i4.4)')'HES',trim(SESSION),mstep
          write(6,*) ''
          open(67,file=trim(filename),status='unknown',form='formatted')
          write(67,*) beta
          do i = 1,mstep+1
            do j = 1,mstep
               write(67,*) H(i,j)
            enddo
          enddo
          close(67)
         endif

        if(nid.EQ.0) write(6,*) "End of iteration:", mstep

      enddo !mstep = mstart, k_dim

      return
      end
c-----------------------------------------------------------------------
      subroutine check_conv(res,have_cnv,beta,tol,H,n) !check convergence
      implicit none

!-----Inputs -----!

      integer              :: n
      real, dimension(n,n) :: H
      real                 :: tol, beta

!-----Outputs -----!

      real, dimension(n) :: res
      integer            :: have_cnv

!-----Miscellaneous -----!

      complex*16, dimension(n)   :: VP
      complex*16, dimension(n,n) :: FP
      integer :: i

!-----Compute the eigedecomposition of H -----!

      VP = ( 0.0D0 , 0.0D0 ); FP = ( 0.0D0 , 0.0D0 )
      call eigen_decomp(H,n,VP,FP)

!-----Compute the number of converged Ritz pairs -----!

      have_cnv = 0; res = 1.0D0
      do i = 1,n
         res(i) = abs( beta * FP(n,i) )
         if ( res(i) .LT. tol ) then
            have_cnv = have_cnv + 1
         endif
      enddo

      return
      end
c----------------------------------------------------------------------
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
      end subroutine

      subroutine compute_eigenvec_schur(A,FP,n)

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
      end
c-----------------------------------------------------------------------
      subroutine Schur_fact(mstart,H_mat,V_x,V_y,V_z,Pressure,residu,beta)

      implicit none
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

c     ----- Miscellaneous -----

      integer :: mstart,n,n2,i,j
      real beta
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

      Q1 = 0.0D0

      b_vec          = 0.0D0
      b_vec(1,k_dim) = beta

c     ----- Compute the Schur decomposition -----

      VP = ( 0.0D0 , 0.0D0 )
      Q1 = 0.0D0

      call Schur_decomp(H_mat,Q1,VP,k_dim,mstart)

      H_mat(1:mstart,mstart+1:k_dim)       = 0.0D0
      H_mat(mstart+1:k_dim,1:mstart)       = 0.0D0
      H_mat(mstart+1:k_dim,mstart+1:k_dim) = 0.0D0

      if ( mstart.EQ.0 ) then

         mstart = 1
         H_mat  = 0.0D0

         call opcopy( V_x(:,mstart), V_y(:,mstart), V_z(:,mstart),
     $        V_x(:,k_dim+1), V_y(:,k_dim+1), V_z(:,k_dim+1) )

         if(ifpo) call copy( Pressure(:,mstart) ,
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

         call matrix_matrix(V_x(:,1:k_dim),     Q1,lt ,k_dim,k_dim)
         call matrix_matrix(V_y(:,1:k_dim),     Q1,lt ,k_dim,k_dim)
          if(if3D)
     $   call matrix_matrix(V_z(:,1:k_dim),     Q1,lt ,k_dim,k_dim)
          if(ifpo)
     $   call matrix_matrix(Pressure(:,1:k_dim),Q1,lt2,k_dim,k_dim)

c     ----- Update the matrix with b(prime)*Q -----

         call matrix_matrix(b_vec,Q1,1,k_dim,k_dim)

         H_mat(mstart+1,1:mstart) = b_vec(1,1:mstart)

         mstart = mstart + 1

         call opcopy(V_x(:,mstart),  V_y(:,mstart),  V_z(:,mstart),
     $               V_x(:,k_dim+1), V_y(:,k_dim+1), V_z(:,k_dim+1))

!           FIX RESTART WITH SCHUR STEP!!!
!         if (nid.eq.0) then
!            open(unit=ifich1,file='Hessenberg.dat',form='formatted')
!            write(ifich1,*) k_dim
!            !write(ifich1,*) mstart
!            write(ifich1,*) beta
!            do i = 1,k_dim
!               do j = 1,k_dim
!                  write(ifich1,*) H_mat(i,j)
!               enddo
!            enddo
!            close(ifich1)
!            open( unit = ifich2,
!     $           file  = 'Hessenberg_info.dat',
!     $           form  = 'formatted')
!            write(ifich2,*) mstart
!            write(ifich2,*) beta
!            close(ifich2)
!         endif

         do j = 1,mstart
            call  whereyouwant('KRY',j)
            call  outpost(V_x(:,j),V_y(:,j),V_z(:,j),Pressure(:,j),t,'KRY')
         enddo

      endif

      return
      end
c----------------------------------------------------------------------
      subroutine matrix_vector_product (fx, fy, fz, fp, qx, qy, qz)
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'ADJOINT'
      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt
      real, dimension(lt)                :: fx, fy, fz, qx, qy, qz
      real, dimension(lt2)               :: fp
      integer n, n2
      n = nx1*ny1*nz1*nelt; n2 = nx2*ny2*nz2*nelt

!     ----- Initial condition -----
      call opcopy(vxp, vyp, vzp, qx, qy, qz)

!     ----- Time-stepper matrix-vector product -----
      time = 0.0D0
      do istep = 1, nsteps

!     ----- Prepare the base flow to vx,vy,vz

        if(ifldbf)then

          if(nid.eq.0)write(6,*)'Copying base flow!'
          call opcopy(vx,vy,vz,ubase,vbase,wbase)

        else

          if(nid.eq.0)then
            if(ifbase.eqv..true.)then
              write(6,*)'Running DNS alongside stability!'
            else
             write(6,*)'Base flow prescribed in useric!'
            endif
          endif

        endif

!      ----- Integrate in time vxp,vyp,vzp on top of vx,vy,vz

         call nek_advance

      enddo

      call opcopy(fx, fy, fz, vxp, vyp, vzp)
      if(ifpo) call copy(fp, prp, n2)

      return
      end
c----------------------------------------------------------------------
      subroutine outpost_ks(VP,FP,Qx,Qy,Qz,Qp,residu) !outposting vectors
      implicit none
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

      real wo1(lt),wo2(lt),wo3(lt),vort(lt,3)
      common /ugrad/ wo1,wo2,wo3,vort

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
      integer :: n, n2, i, j

      real                               :: sampling_period
      real, dimension(k_dim)             :: residu
      real                               :: glsc3,alpha!, beta

      sampling_period = dt*nsteps

      n = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

c     ----- Output all the spectrums and converged eigenmodes -----

      if(nid.EQ.0) then
       open(unit=ifich1,file='Spectre_H.dat'      ,form='formatted')
       open(unit=ifich2,file='Spectre_NS.dat'     ,form='formatted')
       open(unit=ifich3,file='Spectre_NS_conv.dat',form='formatted')
      endif

      do i = 1,k_dim

         time = real(i) !here for the outposted files have different times

         if(nid.EQ.0) then

            !Spectre_H.dat
           write(ifich1,"(3E15.7)") real(VP(i)),
     $                              aimag(VP(i)),
     $                              residu(i)

            !Spectre_NS.dat
           write(ifich2,"(3E15.7)") log(abs(VP(i)))/sampling_period,
     $             ATAN2(aimag(VP(i)),real(VP(i)))/sampling_period,
     $             residu(i)

         endif

         if(residu(i).LT.eigen_tol)then !just the converged oned

            !Spectre_NS_conv.dat
          if(nid.EQ.0)write(ifich3,"(3E15.7)")
     $     log(abs(VP(i)))/sampling_period,
     $     ATAN2(aimag(VP(i)),real(VP(i)))/sampling_period

c     ----- Computation of the corresponding eigenmode -----

            FP_Cx = (0.0D0,0.0D0)
            FP_Cy = (0.0D0,0.0D0)
            if(if3D) FP_Cz = (0.0D0,0.0D0)
            if(ifpo) FP_Cp = (0.0D0,0.0D0)

            call matvec(FP_Cx,Qx(:,1:k_dim),FP(:,i),lt)
            call matvec(FP_Cy,Qy(:,1:k_dim),FP(:,i),lt)
            if(if3D)
     $      call matvec(FP_Cz,Qz(:,1:k_dim),FP(:,i),lt)
            if(ifpo)
     $      call matvec(FP_Cp,Qp(:,1:k_dim),FP(:,i),lt2)

c     ----- Normalization to be unit-norm -----
c     Note: volume integral of FP*conj(FP) = 1.

            alpha = glsc3(real(FP_Cx),bm1,real(FP_Cx),n)
     $           +  glsc3(real(FP_Cy),bm1,real(FP_Cy),n)
     $           +  glsc3(aimag(FP_Cx),bm1,aimag(FP_Cx),n)
     $           +  glsc3(aimag(FP_Cy),bm1,aimag(FP_Cy),n)

            if(if3D) alpha = alpha
     $           +  glsc3(real(FP_Cz),bm1,real(FP_Cz),n)
     $           +  glsc3(aimag(FP_Cz),bm1,aimag(FP_Cz),n)

            alpha = 1.0D0/sqrt(alpha)

c     ----- Output the imaginary part -----

            call opcopy(vx,vy,vz
     $           ,aimag(FP_Cx)
     $           ,aimag(FP_Cy)
     $           ,aimag(FP_Cz))
            if(ifpo) call copy(pr,aimag(FP_Cp),n2)

            call opcmult(vx,vy,vz,alpha)
            if(ifpo) call cmult(pr,alpha,n2)
            ifto = .false.; ifpo = .true.
            call outpost(vx,vy,vz,pr,t,'Im_')

c     ----- Output the real part -----

            call opcopy(vx,vy,vz
     $           ,real(FP_Cx)
     $           ,real(FP_Cy)
     $           ,real(FP_Cz))
            if(ifpo) call copy(pr,real(FP_Cp),n2)

            call opcmult(vx,vy,vz,alpha)
            if(ifpo) call cmult(pr,alpha,n2)
            ifto = .false.; ifpo = .true.
            call outpost(vx,vy,vz,pr,t,'Re_')

c     ----- Output vorticity from real part -----

            call oprzero(wo1,wo2,wo3)
            call comp_vort3(vort,wo1,wo2,vx,vy,vz)
            ifto = .false.; ifpo = .false.
            call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t, 'Rev')
            ifto = .true.; ifpo = .true.

       endif

      enddo

      if(nid.EQ.0)then
         close(ifich1);         close(ifich2);         close(ifich3)
      endif

      return
      end
c---------------------------------------------------------------
      function select_eigenvalue(wr,wi)
      implicit none
      include 'SIZE'
      logical :: select_eigenvalue
      real :: wr, wi
      select_eigenvalue = .false.
      if(sqrt(wr**2+wi**2).GT.(1.0D0-schur_del))select_eigenvalue=.true.
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
      subroutine sort_eigendecomp(VP,FP,n)
      implicit none
      integer                    :: n
      complex*16, dimension(n)   :: VP
      complex*16, dimension(n,n) :: FP
      real, dimension(n)         :: norm
      real                       :: temp_real
      complex*16                 :: temp_complex
      complex*16, dimension(n)   :: temp_n
      integer                    :: k, l
c     ----- Sorting the eigenvalues according to their norm -----
      temp_n   = (0.0D0,0.0D0)
      norm = sqrt( real(VP)**2 + aimag(VP)**2 )
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
      subroutine matvec(FP_Cx,Qx,FP,n)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer i,n
      complex*16, dimension(n)       :: FP_Cx
      complex*16, dimension(k_dim)   :: FP
      real      , dimension(n,k_dim) :: Qx
      FP_Cx = (0.0D0,0.0D0)
      do i = 1,k_dim
         FP_Cx = FP_Cx + Qx(:,i)*FP(i)
      enddo
      return
      end
