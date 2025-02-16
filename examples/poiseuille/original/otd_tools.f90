!
!     Tools used in the computation of the OTD modes
!     Based on the method proposed in
!     Babaee,H. and Sapsis,T.P., "A minimization principle for the
!      descripertion of modes associated with finite-time instabilities",
!      Proc. R. Soc. A 474: 20150779,
!      https://dx.doi.org/10.1098/rspa.2015.0779  
!      
!     Simon Kern
!     Email:     skern@mech.kth.se 
!
!     List of subroutines:
!       Lu_op_conv(uxp,uyp,uzp,ipert)
!                  IN:  uxp,uyp,uzp         perturbation velocity fields   
!                  IN:  ipert               index
!                  OUT: conv[xyz] (OTD)     convective part of Lu
!       Lu_op_gradp(prpert,ipert)
!                  IN:  prpert              perturbation pressure field
!                  IN:  ipert               index
!                  OUT: gradp[xyz] (OTD)    pressure gradient part of Lu
!       Lu_op_diff(uxp,uyp,uzp,ipert)    
!                  IN:  uxp,uyp,uzp         perturbation velocity fields   
!                  IN:  ipert               index
!                  OUT: diff[xyz] (OTD)     diffusive part of Lu
!       laplacian(lapu,up)
!                  IN:  up                  perturbation velocity field
!                  OUT: lapu                laplacian of up
!       op_ip_vp(uxp,uyp,uzp,jpert)
!                  IN:  uxp,uyp,uzp         input velocity fields   
!                  IN:  jpert               index
!       op_ip(ipert,jpert,iflag)
!                  IN:  ipert               index of first pert. field   
!                  IN:  jpert               index of second pert. field
!                  IN:  iflag               flag to switch u_i / Lu_i
!       op_norm(uxp,uyp,uzp)
!                  IN:  uxp,uyp,uzp         input velocity fields   
!       CGS (classical Gram-Schmidt orthonormalisation)  
!                  IN:  ---
!                  OUT: ---
!       MGS (modified Gram-Schmidt orthonormalisation)
!                  IN:  ---
!                  OUT: ---
!       computeNO(N,O,info,flag) (compute measures for normality and orthogonality)
!                  IN:  info (6 chars)      info
!                  IN:  flag (logical)      print N,L
!                  OUT: N,O
!       sorteigs (sort eigenvalues in decreasing order)
!                  IN:  ---
!                  OUT: ---
!     
!======================================================================

!======================================================================
!> @brief Construct the convective terms for Lu
!
!     Lu_conv = (u.grad) Ub + (Ub.grad) u
!
!     Using the convop routine takes care of the dealiasing.
!
      subroutine Lu_op_conv(uxp,uyp,uzp,ipert)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'SOLN'            ! v[xyz]
      include 'OTD'             ! conv[xyz]

      ! argument list
      real uxp (lx1*ly1*lz1*lelv)   ! perturbation velocity components
     $ ,   uyp (lx1*ly1*lz1*lelv)
     $ ,   uzp (lx1*ly1*lz1*lelv)
      integer ipert                 ! index of the considered perturbation

      ! local variables
      integer ntot,i
      real TA1 (LX1,LY1,LZ1,LELV)
     $ ,   TA2 (LX1,LY1,LZ1,LELV)
     $ ,   TA3 (LX1,LY1,LZ1,LELV)
     $ ,   TB1 (LX1,LY1,LZ1,LELV)
     $ ,   TB2 (LX1,LY1,LZ1,LELV)
     $ ,   TB3 (LX1,LY1,LZ1,LELV)
!-----------------------------------------------------------------------

      ntot = lx1*ly1*lz1*lelv
!
      if (if3d) then
        call opcopy  (tb1,tb2,tb3,vx,vy,vz)         ! Save velocity
        call opcopy  (vx,vy,vz,uxp,uyp,uzp)         ! U <-- u
! convop(conv,fld): builds the convective term for the scalar field fld
!       conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
        call convop  (ta1,tb1)                      ! (u.grad) Ub
        call convop  (ta2,tb2)                                
        call convop  (ta3,tb3)
        ! Copy fields into the correct variables
        call opcopy  (convx(1,ipert),convy(1,ipert),convz(1,ipert)
     $ ,              ta1,ta2,ta3)
        call opcopy  (vx,vy,vz,tb1,tb2,tb3)         ! Restore velocity
c
!       conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
        call convop  (tb1,uxp)                      ! (Ub.grad) u
        call convop  (tb2,uyp)
        call convop  (tb3,uzp)
        ! Add fields to the convective term
        call opadd2  (convx(1,ipert),convy(1,ipert),convz(1,ipert)
     $ ,              tb1,tb2,tb3)

!   DIAGNOSTICS      
!        call opcopy  (convx1(1,ipert),convy1(1,ipert),convz1(1,ipert)
!     $ ,              ta1,ta2,ta3)
!        call opcopy  (convx2(1,ipert),convy2(1,ipert),convz2(1,ipert)
!     $ ,              tb1,tb2,tb3)
      else ! 2D
        call opcopy  (tb1,tb2,tb3,vx,vy,vz)         ! Save velocity
        call opcopy  (vx,vy,vz,uxp,uyp,uzp)         ! U <-- u
! convop(conv,fld): builds the convective term for the scalar field fld
!       conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
        call convop  (ta1,tb1)                      ! (u.grad) Ub
        call convop  (ta2,tb2)                                
        ! Copy fields into the correct variables
        call opcopy  (convx(1,ipert),convy(1,ipert),ta3,ta1,ta2,ta3)
        call opcopy  (vx,vy,vz,tb1,tb2,tb3)         ! Restore velocity
c
!       conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
        call convop  (tb1,uxp)                      ! (Ub.grad) u
        call convop  (tb2,uyp)
        ! Add fields to the convective term
        call opadd2  (convx(1,ipert),convy(1,ipert),tb3,tb1,tb2,tb3)
!   DIAGNOSTICS      
!        call opcopy  (convx1(1,ipert),convy1(1,ipert),ta3,ta1,ta2,ta3)
!        call opcopy  (convx2(1,ipert),convy2(1,ipert),tb1,tb1,tb2,tb3)
      endif ! if3d

      return
      end subroutine Lu_op_conv

!======================================================================
!> @brief Construct the pressure gradient term for Lu
!
!     Lu_gradp = grad p
!
!     Note: The pressure gradient is computed directly on the v-mesh!
!
      subroutine Lu_op_gradp(prpert,ipert)
      implicit none

      include 'SIZE'
      include 'OTD'             ! gradp[xyz]

      ! argument list
      real prpert (lx2*ly2*lz2*lelv,1) ! perturbation pressure field
      integer ipert                    ! number of the considered pert.

      ! local variables
      integer ntot
      real ta1 (lx1,ly1,lz1,lelv)  
     $ ,   ta2 (lx1,ly1,lz1,lelv)  
     $ ,   wrk (lx1,ly1,lz1,lelv)  
!-----------------------------------------------------------------------

      ntot = lx2*ly2*lz2*lelv
!     Map the perturbation pressure to the velocity mesh
      call mappr(wrk,prpert,ta1,ta2)
!     compute the gradient on the velocity mesh directly
      call gradm1(gradpx(1,ipert),gradpy(1,ipert),gradpz(1,ipert),wrk)

      return
      end subroutine Lu_op_gradp

!======================================================================
!> @brief Construct the diffusive term for Lu
!
!     Lu_op_diff = 1/Re*grad^2 u
!
      subroutine Lu_op_diff(uxp,uyp,uzp,ipert)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'SOLN'            ! vdiff
      include 'OTD'             ! diff[xyz]

      ! argument list
      real uxp (lx1*ly1*lz1*lelv,1) ! perturbation velocity components
     $ ,   uyp (lx1*ly1*lz1*lelv,1)
     $ ,   uzp (lx1*ly1*lz1*lelv,1)
      integer ipert                 ! number of the considered pert.
      ! local variables
      integer ntot
!-----------------------------------------------------------------------

      ntot = lx1*ly1*lz1*lelv
      ! compute laplacian
      call laplacian(diffx(1,ipert),uxp)
      call laplacian(diffy(1,ipert),uyp)
      if (if3d) call laplacian(diffz(1,ipert),uzp) 
      ! multiply by 1/Re                > remove for operator diagnostics
      call col2(diffx(1,ipert),vdiff,ntot) 
      call col2(diffy(1,ipert),vdiff,ntot) 
      if (if3d) call col2(diffz(1,ipert),vdiff,ntot) 

      return
      end subroutine Lu_op_diff

!======================================================================
!> @brief Construct the diffusion term (laplacian of u) for direction i
!
      subroutine laplacian(lapu,up)
      implicit none
!
      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'DXYZ'            ! dxm1,d[xy]tm1
      include 'GEOM'            ! r[xy]m1,s[xy]m1,t[xy]m1,jacmi
!
      ! argument list
      real up (lx1*ly1*lz1*lelv,1)       ! perturbation velocity component
!
      ! output
      real lapu (lx1*ly1*lz1,lelv)
!
      ! local variables
      real ux  (lx1*ly1*lz1,lelv)
     $ ,   uy  (lx1*ly1*lz1,lelv)
     $ ,   uz  (lx1*ly1*lz1,lelv)
     $ ,   ur  (lx1*ly1*lz1)
     $ ,   us  (lx1*ly1*lz1)
     $ ,   ut  (lx1*ly1*lz1)
!
      common /ctmp1/ ur,us,ut
!
      integer e,i,lxyz,nel,N
!-----------------------------------------------------------------------

      lxyz = lx1*ly1*lz1
      nel = nx1-1
      call gradm1(ux,uy,uz,up)
      do e=1,lelt
        if (if3d) then
          call local_grad3(ur,us,ut,ux,nel,e,dxm1,dxtm1)
          do i=1,lxyz
            lapu(i,e) = jacmi(i,e)*(  ur(i)*rxm1(i,1,1,e)
     $                              + us(i)*sxm1(i,1,1,e)
     $                              + ut(i)*txm1(i,1,1,e) )
          enddo
          call local_grad3(ur,us,ut,uy,nel,e,dxm1,dxtm1)
          do i=1,lxyz
            lapu(i,e) = lapu(i,e) + jacmi(i,e)*(  ur(i)*rym1(i,1,1,e)
     $                                          + us(i)*sym1(i,1,1,e)
     $                                          + ut(i)*tym1(i,1,1,e) )
          enddo
          call local_grad3(ur,us,ut,uz,nel,e,dxm1,dxtm1)
          do i=1,lxyz   
            lapu(i,e) = lapu(i,e) + jacmi(i,e)*(  ur(i)*rzm1(i,1,1,e)
     $                                          + us(i)*szm1(i,1,1,e)
     $                                          + ut(i)*tzm1(i,1,1,e) )
          enddo
        else ! 2D
          call local_grad2(ur,us,ux,nel,e,dxm1,dytm1)
          do i=1,lxyz
            lapu(i,e) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                            + us(i)*sxm1(i,1,1,e) )
          enddo
          call local_grad2(ur,us,uy,nel,e,dxm1,dytm1)
          do i=1,lxyz
            lapu(i,e) = lapu(i,e)
     $                  + jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                              + us(i)*sym1(i,1,1,e) )
          enddo
        endif ! if3d
      enddo
!      
      return
      end subroutine laplacian
!
!======================================================================
!> @brief Compute the global inner product with the pert. velocities 
!
!     < vc_i, v[xyz]p_j >
!      
      real function op_ip_vp(vcx,vcy,vcz,jpert)
      implicit none

      include 'SIZE'
      include 'SOLN'            ! v[xyz]p, jp
      include 'TSTEP'           ! ifield
      include 'FRAMELP'         ! mntr_abort
      include 'MASS'            ! bm1

      ! argument list
      real vcx(lx1*ly1*lz1*lelv)
      real vcy(lx1*ly1*lz1*lelv)
      real vcz(lx1*ly1*lz1*lelv)
      integer jpert
      ! functions and local variables
      real op_glsc2_wt

      ifield = 1
      op_ip_vp = 0.5*op_glsc2_wt(VCX,         VCY,         VCZ,
     $                           VXP(1,jpert),VYP(1,jpert),VZP(1,jpert),
     $                           bm1)

      return
      end function op_ip_vp
!
!======================================================================
!> @brief Compute the global inner product for the perturbations 
!
!     < v[xyz]p_i , v[xyz]p_j >  or  < LU[xyz]_i , v[xyz]p_j >
!      
      real function op_ip(ipert,jpert,iflag)
      implicit none

      include 'SIZE'
      include 'SOLN'            ! v[xyz]p, jp
      include 'FRAMELP'         ! mntr_abort
      include 'OTD'             ! LU[xyz], otd_id

      ! argument list
      integer ipert,jpert       ! [IN]  perturbation indices
      integer iflag             ! [IN]  iflag = 1: < U_i, U_j >
      !                                 iflag = 2: < LU_i, U_j >
      ! functions and local variables
      real op_ip_vp

      if (iflag.eq.1) then
        op_ip = op_ip_vp(VXP(1,ipert),VYP(1,ipert),VZP(1,ipert),jpert)
      elseif (iflag.eq.2) then
        op_ip = op_ip_vp(LUx(1,ipert),LUy(1,ipert),LUz(1,ipert),jpert)
      else
        call mntr_abort(otd_id,'Error in op_ip!')
      endif

      return
      end function op_ip
!
!======================================================================
!> @brief Normalize vector field 
!
!     u_i = v_i/||v_i||
!      
      subroutine op_norm(uxp,uyp,uzp)
      implicit none

      include 'SIZE'
      include 'TSTEP'           ! ifield
      include 'MASS'            ! bm1
      include 'FRAMELP'         ! mntr_abort
      include 'INPUT'           ! if3d
      include 'OTD'

      ! argument list
      real uxp (lx1*ly1*lz1*lelv,1) ! perturbation velocity components
     $ ,   uyp (lx1*ly1*lz1*lelv,1)
     $ ,   uzp (lx1*ly1*lz1*lelv,1)

      ! functions and local variables
      real invnorm,n2,ntot
      real op_glsc2_wt

      ifield = 1
      ntot = lx1*ly1*lz1*lelv
      n2 = 0.5*op_glsc2_wt(uxp,uyp,uzp,uxp,uyp,uzp,bm1)
      if (n2.le.0.0) then
        call mntr_abort(otd_id,'Error in op_norm!')
      endif
      invnorm = 1./sqrt(n2)
      call cmult(uxp,invnorm,ntot)
      call cmult(uyp,invnorm,ntot)
      if (if3d) call cmult(uzp,invnorm,ntot)

      return
      end subroutine op_norm
!
!======================================================================
!> @brief Perform Classical Gram-Schmidt (CGS) orthonormalization on the 
!         perturbation velocity field 
!
!     u_k = v_k - sum_{j = 1 -> k-1} proj_{u_j} (u_k)
!
!        with proj_{u_j} (u_k) = < u_k , u_j >/||u_j|| * u_j
!
      subroutine CGS()
      implicit none
!
      include 'SIZE'
      include 'INPUT'   ! if3d
      include 'TSTEP'   ! istep
      include 'SOLN'    ! V[XYZ]P
      include 'FRAMELP' ! lp_inf
      include 'OTD'     ! otd_id
!
      ! local variables
      integer i,j,ntot
      real invnorm, proj

      ! functions
      real op_ip
!-----------------------------------------------------------------------
      ntot = lx1*ly1*lz1*lelv
      ! orthonormalize
      do i=1,npert
        do j=1,i-1
          proj = -op_ip(i,j,1)/op_ip(j,j,1)
          call add2s2(vxp(1,i),vxp(1,j),proj,ntot)
          call add2s2(vyp(1,i),vyp(1,j),proj,ntot)
          if (if3d) then
            call add2s2(vzp(1,i),vzp(1,j),proj,ntot)
          endif
        enddo
        invnorm = 1/sqrt(op_ip(i,i,1))
        call cmult(vxp(1,i),invnorm,ntot)
        call cmult(vyp(1,i),invnorm,ntot)
        if (if3d) then
          call cmult(vzp(1,i),invnorm,ntot)
        endif
      enddo
! 
!     stamp logs
      call mntr_logi(otd_id,lp_inf,
     $               'V[XYZ]P orthonormalized (CGS)',istep)

      return 
      end subroutine CGS
!
!======================================================================
!> @brief Perform Modified Gram-Schmidt (MGS) orthonormalization on the 
!         perturbation velocity field for improved numerical stability 
!
!     do i=1,npert
!       u_i = v_i/||v_i||
!       do j=i+1,npert
!         u_j = v_j - proj_{u_i} (v_j)
!       enddo
!     enddo
!
!        with proj_{u_i} (v_j) = < v_j , u_i >/||u_i|| * v_j
!                              = < v_j , u_i > * v_j   since ||u_i|| = 1
!
      subroutine MGS()
      implicit none
!
      include 'SIZE'
      include 'INPUT'   ! if3d
      include 'TSTEP'   ! istep
      include 'SOLN'    ! V[XYZ]P
      include 'FRAMELP' ! lp_inf
      include 'OTD'     ! otd_id
!
      ! local variables
      integer i,j,ntot
      real invnorm, proj

      ! functions
      real op_ip
!-----------------------------------------------------------------------
      ntot = lx1*ly1*lz1*lelv
      ! orthonormalize
      do i=1,npert
        invnorm = 1/sqrt(op_ip(i,i,1))
        call cmult(vxp(1,i),invnorm,ntot)
        call cmult(vyp(1,i),invnorm,ntot)
        if (if3d) then
          call cmult(vzp(1,i),invnorm,ntot)
        endif
        do j=i+1,npert
          proj = -op_ip(i,j,1)
          call add2s2(vxp(1,j),vxp(1,i),proj,ntot)
          call add2s2(vyp(1,j),vyp(1,i),proj,ntot)
          if (if3d) then
            call add2s2(vzp(1,j),vzp(1,i),proj,ntot)
          endif
        enddo
      enddo
! 
!     stamp logs
      call mntr_logi(otd_id,lp_inf,
     $               'V[XYZ]P orthonormalized (MGS)',istep)

      return 
      end subroutine MGS
!
!======================================================================
!> @brief Compute measures of orthogonality and normality of the
!         perturbations 
!
      subroutine compute_NO(N,O,info,flag)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'

      ! argument list
      real        N,O
      logical     flag
      character*6 info

      ! local variables
      integer     i,j,n9
      real        ip(npert,npert)

      ! function
      real        op_ip

      do i=1,npert
        do j=1,npert
          ip(i,j) = op_ip(i,j,1)
        enddo
      enddo

      if (nid.eq.0) then
        ! compute measure for basis vector normality
        N = 0
        do i=1,npert
          N = N + ip(i,i)**2
        enddo
        N = sqrt(N/npert)

        ! compute measure for basis vector orthogonality
        if (npert.gt.1) then
         O = 0
          do i=1,npert
            do j=i+1,npert
              O = O + ip(i,j)**2
            enddo
          enddo
          O = sqrt(2*O)/(npert*(npert-1))
        endif

        ! output
        if (otd_debug) then            ! debug output pre MGS
          n9 = min(npert,9)
          write(6,*) ' [OTD] ipout   istp ipert
     $                <u_i,u_j> ON '//trim(info)
          do i=1,npert
            write(6,100) istep,i,(ip(i,j),j=1,n9)
          enddo
        endif
        if (flag) write(6,101) istep,N-1.0,O
      endif
 100  format('  [OTD] ipout',I7,1X,I5,2X,9(1X,F9.6))
 101  format('  [OTD] NOout',I7,1X,'N-1.0',1X,E15.7,1X,'O',1X,E15.7)
    
      return
      end subroutine compute_NO

!======================================================================
!> @brief 1. Sort the eigenvalues l_i such that their real parts are 
!            ranked in decreasing order
!
!            Re(l_1) .ge. Re(l_i) .ge. Re(l_r), i = 1,...,r
!
!         2. Apply the same sorting to the columns of the right 
!            eigenvector matrix and separate real and imaginary parts
!
!            EVR => EVRR + i*EVRI
! 
      subroutine sorteigs(str)
      implicit none
!
      include 'SIZE'
      include 'OTD'     ! EIG[IR], EVR, EVR[RI]
      include 'TSTEP'

      ! argument list
      character*3 str
!
      ! local variables
      integer i,j,id
      logical mk(lpert)
      real    wrk1(lpert)
      real    wrk2(lpert)
!-----------------------------------------------------------------------
      ! zero out indices, output and mask
      call izero(idx,lpert)
      call rzero(EVRR,lpert*lpert)
      call rzero(EVRI,lpert*lpert)
      do i=1,npert
        mk(i) = .true.
      enddo
      call copy(wrk1,EIGR,lpert)
      call copy(wrk2,EIGI,lpert)

      ! we need to exclude the trailing zeros for the sorting to work
      if (lpert.gt.npert) then
        do i=npert+1,lpert
          mk(i) = .false.
        enddo
      endif
      if (otd_debug) then
        if (nid.eq.0) then
          call print_eigenvalues(' [OTD debug] e-vals: '//trim(str)
     $ ,                         npert,EIGR,EIGI)
        call print_eigenvectors(' [OTD debug] right e-vecs: '//trim(str)
     $ ,                         npert,EIGI,EVR,npert)
        endif
      endif
      ! sorting
      j = 1
      do while (j.le.npert)
        EIGR(j)    = maxval(wrk1,mask=mk)       ! find largest real eigenvalue in remaining list
        id         = maxloc(wrk1,1,mk)          ! find its index
        EIGI(j)    = wrk2(id)                   ! extract corresponding imaginary part
        mk(id)     = .false.                    ! update mask
        idx(j)     = id
        if (EIGI(j).eq.0.0) then
          do i=1,npert
            EVRR(i,j) = EVR(i,id)
            ! EVRI(i,j) = 0.0
          enddo
          j=j+1
        else  ! complex conjugate eigenvectors!
          EIGI(j+1)    = -EIGI(j)
          EIGR(j+1)    =  EIGR(j)
          mk(id+1)     = .false.
          idx(j+1)     = id+1
          do i=1,npert
            EVRR(i,j)   =  EVR(i,id)
            EVRR(i,j+1) =  EVR(i,id)
            EVRI(i,j)   =  EVR(i,id+1)
            EVRI(i,j+1) = -EVR(i,id+1)
          enddo
          j=j+2
        endif
      enddo

! debugging output
      if (otd_debug) then
        if (nid.eq.0) then
        call print_eigenvalues(' [OTD debug] sorted e-vals: '//trim(str)
     $ ,                         npert,EIGR,EIGI)
          write(6,*)
          write(6,*) ' [OTD debug] sorted right e-vecs: '//trim(str)
          do i=1,npert
            do j=1,npert
              write(6,200,ADVANCE='NO') EVRR(i,j), EVRI(i,j)
            enddo
            write(6,*)
          enddo
          write(6,*)
        endif
      endif
  200 format( 9(:,3X,F6.2,' + i*',F6.2) )

      return
      end subroutine sorteigs

