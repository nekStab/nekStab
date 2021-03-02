c-----------------------------------------------------------------------
      subroutine bf_sensitivity

         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, parameter :: lt  = lx1*ly1*lz1*lelt

         real, dimension(lt) :: vx_dRe,vy_dRe,vz_dRe
         real, dimension(lt) :: vx_dIm,vy_dIm,vz_dIm
         real, dimension(lt) :: vx_aRe,vy_aRe,vz_aRe
         real, dimension(lt) :: vx_aIm,vy_aIm,vz_aIm

         real, dimension(lt) :: dudx_dRe, dudy_dRe, dudz_dRe
         real, dimension(lt) :: dvdx_dRe, dvdy_dRe, dvdz_dRe
         real, dimension(lt) :: dwdx_dRe, dwdy_dRe, dwdz_dRe

         real, dimension(lt) :: dudx_dIm, dudy_dIm, dudz_dIm
         real, dimension(lt) :: dvdx_dIm, dvdy_dIm, dvdz_dIm
         real, dimension(lt) :: dwdx_dIm, dwdy_dIm, dwdz_dIm

         real, dimension(lt) :: dudx_aRe, dudy_aRe, dudz_aRe
         real, dimension(lt) :: dvdx_aRe, dvdy_aRe, dvdz_aRe
         real, dimension(lt) :: dwdx_aRe, dwdy_aRe, dwdz_aRe
         
         real, dimension(lt) :: dudx_aIm, dudy_aIm, dudz_aIm
         real, dimension(lt) :: dvdx_aIm, dvdy_aIm, dvdz_aIm
         real, dimension(lt) :: dwdx_aIm, dwdy_aIm, dwdz_aIm

         real, dimension(lt) :: vx_tr,vy_tr,vz_tr
         real, dimension(lt) :: vx_ti,vy_ti,vz_ti
         !real, dimension(lt) :: vx_tr,vy_tr,vz_tr
         !real, dimension(lt) :: vx_tr,vy_tr,vz_tr

         ifto=.false.;ifpo=.false.

         call load_fld('dRe_1cyl0.f00001')
         call opcopy(vx_dRe,vy_dRe,vz_dRe,vx,vy,vz)
         call gradm1(dudx_dRe, dudy_dRe, dudz_dRe, vx_dRe ,nelv)
         call gradm1(dvdx_dRe, dvdy_dRe, dvdz_dRe, vy_dRe ,nelv)
         call gradm1(dwdx_dRe, dwdy_dRe, dwdz_dRe, vz_dRe ,nelv)
         call dsavg(dudx_dRe);call dsavg(dudy_dRe); call dsavg(dudz_dRe)
         call dsavg(dvdx_dRe);call dsavg(dvdy_dRe); call dsavg(dvdz_dRe)
         call dsavg(dwdx_dRe);call dsavg(dwdy_dRe); call dsavg(dwdz_dRe)

         call load_fld('dIm_1cyl0.f00001')
         call opcopy(vx_dIm,vy_dIm,vz_dIm,vx,vy,vz)
         call gradm1(dudx_dIm, dudy_dIm, dudz_dIm, vx_dIm ,nelv)
         call gradm1(dvdx_dIm, dvdy_dIm, dvdz_dIm, vy_dIm ,nelv)
         call gradm1(dwdx_dIm, dwdy_dIm, dwdz_dIm, vz_dIm ,nelv)
         call dsavg(dudx_dIm);call dsavg(dudy_dIm); call dsavg(dudz_dIm)
         call dsavg(dvdx_dIm);call dsavg(dvdy_dIm); call dsavg(dvdz_dIm)
         call dsavg(dwdx_dIm);call dsavg(dwdy_dIm); call dsavg(dwdz_dIm)

         call load_fld('aRe_1cyl0.f00001')
         call opcopy(vx_aRe,vy_aRe,vz_aRe,vx,vy,vz)
         call gradm1(dudx_aRe, dudy_aRe, dudz_aRe, vx_aRe ,nelv)
         call gradm1(dvdx_aRe, dvdy_aRe, dvdz_aRe, vy_aRe ,nelv)
         call gradm1(dwdx_aRe, dwdy_aRe, dwdz_aRe, vz_aRe ,nelv)
         call dsavg(dudx_aRe);call dsavg(dudy_aRe); call dsavg(dudz_aRe)
         call dsavg(dvdx_aRe);call dsavg(dvdy_aRe); call dsavg(dvdz_aRe)
         call dsavg(dwdx_aRe);call dsavg(dwdy_aRe); call dsavg(dwdz_aRe)

         call load_fld('aIm_1cyl0.f00001') 
         call opcopy(vx_aIm,vy_aIm,vz_aIm,vx,vy,vz)
         call gradm1(dudx_aIm, dudy_aIm, dudz_aIm, vx_aIm ,nelv)
         call gradm1(dvdx_aIm, dvdy_aIm, dvdz_aIm, vy_aIm ,nelv)
         call gradm1(dwdx_aIm, dwdy_aIm, dwdz_aIm, vz_aIm ,nelv)

         call dsavg(dudx_aIm);call dsavg(dudy_aIm); call dsavg(dudz_aIm)
         call dsavg(dvdx_aIm);call dsavg(dvdy_aIm); call dsavg(dvdz_aIm)
         call dsavg(dwdx_aIm);call dsavg(dwdy_aIm); call dsavg(dwdz_aIm)

         call oprzero  (vx_tr,vy_tr,vz_tr)
         call opaddcol3(vx_tr,vy_tr,vz_tr,-vx_aRe,-vx_aRe,-vx_aRe,dudx_dRe,dvdx_dRe,dwdx_dRe)
         call opaddcol3(vx_tr,vy_tr,vz_tr,-vy_aRe,-vy_aRe,-vy_aRe,dudx_dRe,dvdx_dRe,dwdx_dRe)
         call opaddcol3(vx_tr,vy_tr,vz_tr,-vz_aRe,-vz_aRe,-vz_aRe,dudx_dRe,dvdx_dRe,dwdx_dRe)
         call opaddcol3(vx_tr,vy_tr,vz_tr,-vx_aIm,-vx_aIm,-vx_aIm,dudx_dIm,dvdx_dIm,dwdx_dIm)
         call opaddcol3(vx_tr,vy_tr,vz_tr,-vy_aIm,-vy_aIm,-vy_aIm,dudx_dIm,dvdx_dIm,dwdx_dIm)
         call opaddcol3(vx_tr,vy_tr,vz_tr,-vz_aIm,-vz_aIm,-vz_aIm,dudx_dIm,dvdx_dIm,dwdx_dIm)

      !vx_tr = -vx_aRe*dudx_dRe  -vy_aRe*dudy_dRe  -vz_aRe*dudz_dRe  -vx_aIm*dudx_dIm  -vy_aIm*dudy_dIm   -vz_aIm*dudz_dIm
      !vy_tr = -vx_aRe*dvdx_dRe  -vy_aRe*dvdy_dRe  -vz_aRe*dvdz_dRe  -vx_aIm*dvdx_dIm  -vy_aIm*dvdy_dIm   -vz_aIm*dvdz_dIm
      !vz_tr = -vx_aRe*dwdx_dRe  -vy_aRe*dwdy_dRe  -vz_aRe*dwdz_dRe  -vx_aIm*dwdx_dIm  -vy_aIm*dwdy_dIm   -vz_aIm*dwdz_dIm

      ! vx_pr =  vx_aRe*dudx_dIm+vy_aRe*dudy_dIm+vz_aRe*dudz_dIm-vx_aIm*dudx_dRe-vy_aIm*dudy_dRe-vz_aIm*dudz_dRe
      ! vy_pr =  vx_aRe*dvdx_dIm+vy_aRe*dvdy_dIm+vz_aRe*dvdz_dIm-vx_aIm*dvdx_dRe-vy_aIm*dvdy_dRe-vz_aIm*dvdz_dRe
      ! vz_pr =  vx_aRe*dwdx_dIm+vy_aRe*dwdy_dIm+vz_aRe*dwdz_dIm-vx_aIm*dwdx_dRe-vy_aIm*dwdy_dRe-vz_aIm*dwdz_dRe
! 
      ! vx_ti =  vx_dRe*dudx_dIm+v_r*dvdx_aRe+w_r*dwdx_aRe+u_i*dudx_aIm+v_i*dvdx_aIm+w_i*dwdx_aIm
      ! vy_ti =  vx_dRe*dudy_dIm+v_r*dvdy_aRe+w_r*dwdy_aRe+u_i*dudy_aIm+v_i*dvdy_aIm+w_i*dwdy_aIm
      ! vz_ti =  vx_dRe*dudz_dIm+v_r*dvdz_aRe+w_r*dwdz_aRe+u_i*dudz_aIm+v_i*dvdz_aIm+w_i*dwdz_aIm
! 
      ! vx_pi =  vx_dRe*dudx_aIm+v_r*dvdx_aIm+w_r*dwdx_aIm-u_i*dudx_aRe-v_i*dvdx_aRe-w_i*dwdx_aRe
      ! vy_pi =  vx_dRe*dudy_aIm+v_r*dvdy_aIm+w_r*dwdy_aIm-u_i*dudy_aRe-v_i*dvdy_aRe-w_i*dwdy_aRe
      ! vz_pi =  vx_dRe*dudz_aIm+v_r*dvdz_aIm+w_r*dwdz_aIm-u_i*dudz_aRe-v_i*dvdz_aRe-w_i*dwdz_aRe

         call filter_s0(vx_tr,0.5,1,'vortx') 
         call filter_s0(vy_tr,0.5,1,'vortx') 
         call filter_s0(vz_tr,0.5,1,'vortx') 
         call outpost(vx_tr, vy_tr, vz_tr, pr,t, 'tr_')

         call filter_s0(vx_ti,0.5,1,'vortx') 
         call filter_s0(vy_ti,0.5,1,'vortx') 
         call filter_s0(vz_ti,0.5,1,'vortx') 
         call outpost(vx_ti, vy_ti, vz_ti, pr,t, 'pr_')

         call oprzero(vx,vy,vz)
         call opaddcol3(vx,vy,vz, vx_tr,vy_tr,vz_tr, vx_tr,vy_tr,vz_tr)
         call opaddcol3(vx,vy,vz, vx_ti,vy_ti,vz_ti, vx_ti,vy_ti,vz_ti)
         call vsqrt(vx); call vsqrt(vy); call vsqrt(vz)
         call outpost(vx,vy,vz,pr,t,'st_')

         ! call filter_s0(vx_tr,0.5,1,'vortx') 
         ! call filter_s0(vy_tr,0.5,1,'vortx') 
         ! call filter_s0(vz_tr,0.5,1,'vortx') 
         ! call outpost(vx_ti, vy_ti, vz_ti, pr,t, 'ti_')

         ! call filter_s0(vx_tr,0.5,1,'vortx')
         ! call filter_s0(vy_tr,0.5,1,'vortx')
         ! call filter_s0(vz_tr,0.5,1,'vortx')
         ! call outpost(vx_pi, vy_pi, vz_pi, pr,t, 'pi_')

      return
      end
c-----------------------------------------------------------------------
      subroutine vortex_core(l2, vortex) 
      include 'SIZE'
      include 'TOTAL'
      real l2(lx1,ly1,lz1,1)
      character*(*) vortex

      if     (vortex.EQ."lambda2") then
         call lambda2(l2)
      elseif (vortex.EQ."q") then
         call compute_q(l2)
      elseif (vortex.EQ."delta") then
         call compute_delta(l2)
      elseif (vortex.EQ."swirling") then
         call compute_swirling(l2)
      elseif (vortex.EQ."omega") then
         call compute_omega_jc(l2)
      elseif (vortex.EQ."symmetric") then
         call compute_symmetricVec(l2)
      elseif (vortex.EQ."assymetric") then
         call compute_assymetricVec(l2)
      else
            if(nid.eq.0) write(6,*)"ABORT:unknown vortex determinant", vortex,
     &". Please use one of: 'lambda2', 'q', 'delta', or 'swirling'"
            call exitt
      endif

      return
      end
c---------------------------------------------------------
      subroutine compute_omega_jc(omega)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: n = lx1 * ly1 * lz1 * lelv
      integer, parameter :: nxyz = lx1 * ly1 * lz1
      real, parameter :: eps = 1.0D-5
 
      ! --> Element-wise velocity pseudo-gradient.
      real gije(lx1*ly1*lz1, ldim, ldim)

      ! --> Point-wise symmetric and anti-symmetric parts.
      real ss(ldim, ldim), oo(ldim, ldim)
      real omega(lx1, ly1, lz1, lelv)
      real norm_a, norm_b

      ! --> Miscellaneous.
      integer ie, l, i, j

      ! --> Loop through the elements.
      do ie = 1, nelv

        ! --> Compute the velocity pseudo-gradient.
        call comp_gije(gije, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
     
        ! --> Point-wise computations.
        do l = 1, nxyz
        
        do j = 1, ldim
           do i = 1, ldim
              ! --> Compute the symmetric and antisymmetric component.
              ss(i, j) = 0.50d0 * (gije(l, i, j) + gije(l, j, i))
              oo(i, j) = 0.50d0 * (gije(l, i, j) - gije(l, j, i))
           enddo
        enddo

        ! --> Compute the Frobenius norm
        norm_a = (norm2(ss) ** 2) ** 2 
        norm_b = (norm2(oo) ** 2) ** 2

        ! --> Compute omega.
        omega(l, 1, 1, ie) =  norm_b / (norm_a + norm_b + eps)
        
       enddo   
      enddo
      call filter_s0(omega,0.5,1,'vortx') !filtering is necessary here!
      return
      end
c---------------------------------------------------------
      subroutine compute_omega(l2)
c
c     Generate \Omega criterion vortex of 
c
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lxyz=lx1*ly1*lz1
      real l2(lx1,ly1,lz1,1)
      real mygi(lxyz,ldim,ldim)
      integer n,ie,l,nxyz
      real a,b
      common /mygrad/ mygi
      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv
      do ie=1,nelv ! Compute velocity gradient tensor
         call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
         do l=1,nxyz
          call compute_symmetric    (A,l)
          call compute_antisymmetric(B,l)
          l2(l,1,1,ie) = (B**2)/(B**2+A**2+0.0020d0*glmax_qc)
         enddo
      enddo
      call filter_s0(l2,0.5,1,'vortx')
      return
      end
c---------------------------------------------------------
      subroutine compute_symmetricVec(l2)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lxyz=lx1*ly1*lz1
      real l2(lx1,ly1,lz1,1)
      real mygi(lxyz,ldim,ldim)
      integer n,ie,l,nxyz
      common /mygrad/ mygi
      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv
      do ie=1,nelv ! Compute velocity gradient tensor
         call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
         do l=1,nxyz
          call compute_symmetric    (l2(l,1,1,ie),l)
         enddo
      enddo
      call filter_s0(l2,0.5,1,'vortx')
      return
      end
c---------------------------------------------------------
      subroutine compute_assymetricVec(l2)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lxyz=lx1*ly1*lz1
      real l2(lx1,ly1,lz1,1)
      real mygi(lxyz,ldim,ldim)
      integer n,ie,l,nxyz
      common /mygrad/ mygi
      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv
      do ie=1,nelv ! Compute velocity gradient tensor
         call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
         do l=1,nxyz
          call compute_antisymmetric(l2(l,1,1,ie),l)
         enddo
      enddo
      call filter_s0(l2,0.5,1,'vortx')
      return
      end
c---------------------------------------------------------
       subroutine compute_q(l2)
c
c     Generate Q criterion vortex of Hunt, Wray & Moin, CTR-S88 1988
c     positive second invariant of velocity gradient tensor
c
      include 'SIZE'
      include 'TOTAL'
      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lxyz=lx1*ly1*lz1)
      real l2(lx1,ly1,lz1,1)
      real mygi(lxyz,ldim,ldim)
      real Q1
      common /mygrad/ mygi

      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv

      do ie=1,nelv
         call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
         do l=1,nxyz
            call compute_secondInv(Q1,l)         
            l2(l,1,1,ie) = Q1
         enddo
      enddo
      call filter_s0(l2,0.5,1,'vortx') 

      return
      end
c---------------------------------------------------------
      subroutine compute_delta(l2)
c
c     Generate  Discriminant (DELTA) criterion vortex of Chong, Perry & Cantwell, Phys. Fluids 1990
c     complex eigenvalues of velocity gradient tensor
c
      include 'SIZE'
      include 'TOTAL'
      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lxyz=lx1*ly1*lz1)
      real l2(lx1,ly1,lz1,1)
      real mygi(lxyz,ldim,ldim)
      real P1,Q1,R1,Q,R,Delta
      common /mygrad/ mygi

      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv

      do ie=1,nelv
         ! Compute velocity gradient tensor
         call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
         do l=1,nxyz
            call compute_firstInv(P1,l)
            P1=-P1 !negative sign
            call compute_secondInv(Q1,l) 
            call compute_thirdInv(R1,l)
            R1=-R1 !negative sign
            Q=Q1-(P1**2)/3
            R=R1+(P1**3)*2/27-P1*Q1/3   
            l2(l,1,1,ie) = (R/2)**2+(Q/3)**3
         enddo
      enddo
      call filter_s0(l2,0.5,1,'vortx')
      return
      end
c---------------------------------------------------------
       subroutine compute_swirling(l2)
c
c     Generate Swirling Strength criterion vortex of 
c     Zhou, Adrian, Balachandar and Kendall, JFM 1999
c
c     imaginary part of complex eigenvalues of velocity gradient tensor
c                   presented as lambda_ci^2 
c
c   for the 3x3 case
c                |d11, d12, d13|  |du/dx,du/dy,du/dz|
c   D = [d_ij] = |d21, d22, d23|= |dv/dx,dv/dy,dv/dz|
c                |d31, d32, d33|  |dw/dx,dw/dy,dw/dz|
c
c                       |lambda_r,  0,       0      |
c   [d_ij]= [vr,vcr,vci]|  0,   lambda_cr, lambda_ci|[vr,vcr,vci]^-1
c                       |  0,  -lambda_ci, lambda_cr| 
c
c    λ^3+P*λ^2+Qλ+R=0
c    
c   where P = -I_D   = -tr(D)
c         Q = II_D  =  0.5[P*P - tr(DD)]
c         R = -III_D = 1/3.[-P+3QP-tr(DDD)]
c   for the 2x2 case
c
c   λ^2+Qλ+R=0
c
c   where Q = II_D  =  0.5[P*P - tr(DD)]
c         R = -III_D = 1/3.[-P+3QP-tr(DDD)]    
c
c
      include 'SIZE'
      include 'TOTAL'
      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (lxyz=lx1*ly1*lz1)
      real l2(lx1,ly1,lz1,1)
      real gije(lxyz,ldim,ldim),mygi(lxyz,ldim,ldim)
      real P1,Q1,R1,Q,R,Delta,a1,a2,lambdaCi
      common /mygrad/ mygi

      nxyz = lx1*ly1*lz1
      n    = nxyz*nelv
      if (if3d) then     ! 3D CASE
         do ie=1,nelv
         ! Compute velocity gradient tensor
        call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
           do l=1,nxyz
              call compute_firstInv(P1,l)
              P1=-P1 !negative sign
              call compute_secondInv(Q1,l) 
              call compute_thirdInv(R1,l)
              R1=-R1 !negative sign
              !Q=Q1-(P1**2)/3
              R=R1+(P1**3)*2/27-P1*Q1/3
              !Delta = (R/2)**2+(Q/3)**3          
              call  cubicLambdaCi(P1, Q1, R1, lambdaCi)
                  
            l2(l,1,1,ie) =lambdaCi
           enddo
        enddo
      elseif (ifaxis) then   ! AXISYMMETRIC CASE
            if(nid.eq.0) write(6,*) 
     &      'ABORT:no compute_swirling axialsymmetric support for now'
            call exitt
       else              ! 2D CASE
         do ie=1,nelv
         ! Compute velocity gradient tensor
        call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
           do l=1,nxyz
              call compute_firstInv(P1,l)
              P1=-P1 !negative sign
              !call compute_secondInv(Q1,l) 
              !Q1=0.
              call compute_thirdInv(R1,l)
              R1=-R1 !negative sign
              !Q=Q1-(P1**2)/3
              !R=R1+(P1**3)*2/27-P1*Q1/3
              !Delta = (R/2)**2+(Q/3)**3          
              call  quadLambdaCi(P1, R1, lambdaCi)

                  
            l2(l,1,1,ie) =lambdaCi
           enddo
        enddo

      endif
     
      call filter_s0(l2,0.5,1,'vortx') 
        
      do ie=1,nelv
       do l=1,nxyz
        l2(l,1,1,ie) = l2(l,1,1,ie)*l2(l,1,1,ie)  
       enddo
      enddo

      return
      end
c---------------------------------------------------------
       subroutine compute_antisymmetric(B,l) !B=0.5(\Nabla u - (\Nabla u)^T)
c      |11(dudx),12(dudy),13(dudz)|
c      |21(dvdx),22(dvdy),23(dvdz)|
c      |31(dwdx),32(dwdy),33(dwdz)|
       implicit none
       include 'SIZE'
       include 'TOTAL'
       integer, parameter :: lxyz=lx1*ly1*lz1
       real mygi(lxyz,ldim,ldim),B
       integer l
       common /mygrad/ mygi
       B=0.0d0
       if(if3d)then ! 3D CASE
        B=((mygi(l,1,2)+mygi(l,2,1))**2+(mygi(l,1,3)+mygi(l,3,1))**2+(mygi(l,2,3)+mygi(l,3,2))**2)/4
       else         ! 2D CASE
        B=((mygi(l,1,2)+mygi(l,2,1))**2)/4
       endif 
       return
       end
c---------------------------------------------------------
       subroutine compute_symmetric(A,l) !A=0.5(\Nabla u + (\Nabla u)^T)
       implicit none
       include 'SIZE'
       include 'TOTAL'
       integer, parameter :: lxyz=lx1*ly1*lz1
       real mygi(lxyz,ldim,ldim),A,B
       integer l
       common /mygrad/ mygi
       A=0.0d0;B=0.0d0
       if(if3d)then ! 3D CASE
        B=((mygi(l,1,2)+mygi(l,2,1))**2+(mygi(l,1,3)+mygi(l,3,1))**2+(mygi(l,2,3)+mygi(l,3,2))**2)/4
        A=B+(mygi(l,1,1)**2+mygi(l,2,2)**2+mygi(l,3,3)**2)/2
       else         ! 2D CASE
        B=((mygi(l,1,2)+mygi(l,2,1))**2)/4
        A=B+(mygi(l,1,1)**2+mygi(l,2,2)**2)/2
       endif 
       return
       end
c---------------------------------------------------------
       subroutine compute_firstInv(a, l)
c
c   for the 3x3 case
c                |d11, d12, d13|
c   D = [d_ij] = |d21, d22, d23|
c                |d31, d32, d33|
c   compute_firstInv returns tr(d_ij)=(d11+d22+d33)
c
       include 'SIZE'
       include 'TOTAL'
       parameter(lxyz=lx1*ly1*lz1)
       real mygi(lxyz,ldim,ldim),a
       integer l
       common /mygrad/ mygi
       if(if3d)then ! 3D CASE
        a=(mygi(l,1,1)+mygi(l,2,2)+mygi(l,3,3))
       elseif(ifaxis)then ! AXISYMMETRIC CASE
        if(nid.eq.0)write(6,*)'ABORT: compute_firstInv axialsymmetric support for now'
        call exitt
       else ! 2D CASE
        a=(mygi(l,1,1)+mygi(l,2,2))
       endif 
       return 
       end
c-------------------------------------------------------------------
       subroutine compute_secondInv(a, l)
c
c   for the 3x3 case
c                |d11, d12, d13|
c   D = [d_ij] = |d21, d22, d23|
c                |d31, d32, d33|
c   compute_SecondInv returns 
c    0.5*(tr(d_ij)^2+tr(d_ij*d_ij))=-(d22*d33-d23*d32)-(d11*d22-d12*d21)-(d33*d11-d13*d31)
c
       include 'SIZE'
       include 'TOTAL'
       parameter (lxyz=lx1*ly1*lz1)
       real mygi(lxyz,ldim,ldim)
       real a
       integer l
       common /mygrad/ mygi
       if (if3d) then     ! 3D CASE
        a=(mygi(l,2,2)*mygi(l,3,3)-mygi(l,2,3)*mygi(l,3,2))
     &   +(mygi(l,1,1)*mygi(l,2,2)-mygi(l,1,2)*mygi(l,2,1))
     &   +(mygi(l,3,3)*mygi(l,1,1)-mygi(l,1,3)*mygi(l,3,1))
       elseif (ifaxis) then   ! AXISYMMETRIC CASE
        if(nid.eq.0) write(6,*)'ABORT: compute_secondInv axialsymmetric support for now'
        call exitt
       else              ! 2D CASE
        a=( mygi(l,1,1)+mygi(l,2,2) )*( mygi(l,1,1)+mygi(l,2,2) )
     &  -2*mygi(l,1,2)*mygi(l,2,1)-mygi(l,1,1)*mygi(l,1,1)-mygi(l,2,2)*mygi(l,2,2)
       endif 

       return 
       end
c-------------------------------------------------------------------
       subroutine compute_thirdInv(a, l)
c
c   for the 3x3 case
c                |d11, d12, d13|
c   D = [d_ij] = |d21, d22, d23|
c                |d31, d32, d33|
c
c    compute_thirdInv returns det(D)
c    =-d11*(d23*d32-d22*d33)-d12*(d21*d33-d31*d23)-d13*(d31*d22-d21*d32)
c
       include 'SIZE'
       include 'TOTAL'
       parameter (lxyz=lx1*ly1*lz1)
       real mygi(lxyz,ldim,ldim)
       real a
       integer l
       common /mygrad/ mygi
       if (if3d) then     ! 3D CASE
       a=-mygi(l,1,1)*(mygi(l,2,3)*mygi(l,3,2)
     &         -mygi(l,2,2)*mygi(l,3,3))
     &         -mygi(l,1,2)*(mygi(l,2,1)*mygi(l,3,3)
     &         -mygi(l,3,1)*mygi(l,2,3))
     &         -mygi(l,1,3)*(mygi(l,3,1)*mygi(l,2,2)
     &         -mygi(l,2,1)*mygi(l,3,2))
       elseif (ifaxis) then   ! AXISYMMETRIC CASE
            if(nid.eq.0) write(6,*) 
     &        'ABORT: compute_thirdInv axialsymmetric support for now'
            call exitt
       else              ! 2D CASE
            a=mygi(l,1,1)*mygi(l,2,2)-mygi(l,1,2)*mygi(l,2,1)

       endif 
       return 
       end
c-------------------------------------------------------------------
      subroutine cubicLambdaCi(b, c, d, lci)
      include 'SIZE'
      include 'TOTAL'
      parameter (lxyz=lx1*ly1*lz1)
      real mygi(lxyz,ldim,ldim)
      real b,c,d
      real f,g,h,r,s,t2,u
      real lci      
      complex ci
      complex x1

      ci =sqrt(cmplx(-1.))

      a = 1
      f = c/a - 1/3.*(b/a)**2.
      g = ((2.*(b**3.)/(a**3.))-(9.*b*c/(a**2.))+(27.*d/a))/27.
      h = (g/2.)**2.+(f/3.)**3. 
        
      if (h .LE. 0.) then
          lci = 0.
      else
          r = -(g/2.) + sqrt(h)
          if (r .le. 0.) then
              s = sign(abs(r)**(1.0/3.0), r)
          else
              s = (r)**(1./3.)
          endif
          t2 = -(g/2.) - sqrt(h)
          if (t2 .le. 0.) then 
              u = sign(abs(t2)**(1.0/3.0), t2)
         else
              u = ((t2)**(1/3.))
         endif
         x1 = -(s+u)/2.+ (b/3./a)+ci*(s-u)*sqrt(3.)/2.
         lci = aimag(x1)
      endif
    
      end 
c-------------------------------------------------------------------
      subroutine quadLambdaCi(b, c, lci)
      include 'SIZE'
      include 'TOTAL'
      parameter (lxyz=lx1*ly1*lz1)
      real mygi(lxyz,ldim,ldim)
      real b,c
      real  f
      real lci      
      complex ci
      complex x1

      ci =sqrt(cmplx(-1.))

      a = 1
      d = b**2.-4.*a*c
        
      if (d .GE. 0.) then
          lci = 0.
      else
          f = sqrt(abs(d))!sign(abs(d)**(1.0/2.0), d)
          x1 = -b/2./a+f*ci/2./a     
         lci = aimag(x1)
      endif
    
      end 
c-------------------------------------------------------------------
      subroutine LambdaCi2D(lci)
      include 'SIZE'
      include 'TOTAL'
      parameter (lxyz=lx1*ly1*lz1)
      real mygi(lxyz,ldim,ldim)
      real d
      real lci      
      complex ci
      complex x1
      common /mygrad/ mygi
      ci =sqrt(cmplx(-1.))


      
      d = (mygi(l,1,2)-mygi(l,2,1))**2-4*mygi(l,1,2)*mygi(l,2,1)

      if (d .GE. 0.) then
         lci = 0.
      else
         lci = sqrt(abs(d))/2.
      endif

      end 
c-------------------------------------------------------------------
c-----------------------------------------------------------------------
        subroutine avg_all_time
        include 'SIZE'  
        include 'TOTAL' 
        include 'AVG'

        integer, parameter :: lt=lx1*ly1*lz1*lelt
        real do1(lt),do2(lt),do3(lt)

        real,save :: x0(3)
        data x0 /0.0, 0.0, 0.0/ 

        integer,save :: icalld
        data    icalld /0/

        logical ifverbose

        if (ax1.ne.lx1 .or. ay1.ne.ly1 .or. az1.ne.lz1) then
         if(nid.eq.0) write(6,*)'ABORT: wrong size of ax1,ay1,az1 in avg_all(), check SIZE!'
         call exitt
        endif
        if (ax2.ne.lx2 .or. ay2.ne.ay2 .or. az2.ne.lz2) then
         if(nid.eq.0) write(6,*)'ABORT: wrong size of ax2,ay2,az2 in avg_all(), check SIZE!'
         call exitt
        endif

        ntot  = lx1*ly1*lz1*nelv
        ntott = lx1*ly1*lz1*nelt
        nto2  = lx2*ly2*lz2*nelv

        ! initialization
        if (icalld.eq.0) then
         icalld = icalld + 1
         atime  = 0.; timel  = time
         call oprzero(uavg,vavg,wavg)
         call rzero(pavg,nto2)
         call oprzero(urms,vrms,wrms)
         call rzero(prms,nto2)
         call oprzero(vwms,wums,uvms)
        endif

        dtime = time  - timel
        atime = atime + dtime

        ! dump freq
        iastep = param(68)
        if  (iastep.eq.0) iastep=param(15)   ! same as iostep
        if  (iastep.eq.0) iastep=500

        ifverbose=.false.
        if (istep.le.10) ifverbose=.true.
        if  (mod(istep,iastep).eq.0) ifverbose=.true.

        if (atime.ne.0..and.dtime.ne.0.) then
        if(nio.eq.0) write(6,*) 'Compute statistics ...'
        beta  = dtime/atime
        alpha = 1.-beta

        ! compute averages E(X) !avg(k) = alpha*avg(k) + beta*f(k)
        call avg1    (uavg,vx,alpha,beta,ntot ,'um  ',ifverbose)
        call avg1    (vavg,vy,alpha,beta,ntot ,'vm  ',ifverbose)
        call avg1    (wavg,vz,alpha,beta,ntot ,'wm  ',ifverbose)
        call avg1    (pavg,pr,alpha,beta,nto2 ,'prm ',ifverbose)
        
        ! compute averages E(X^2) 
        call avg2    (urms,vx,alpha,beta,ntot ,'ums ',ifverbose)
        call avg2    (vrms,vy,alpha,beta,ntot ,'vms ',ifverbose)
        call avg2    (wrms,vz,alpha,beta,ntot ,'wms ',ifverbose)
        !call avg2    (prms,pr,alpha,beta,nto2 ,'prms',ifverbose)
        
        ! compute averages E(X*Y) (for now just for the velocities)
        !call avg3    (uvms,vx,vy,alpha,beta,ntot,'uvm ',ifverbose)
        !call avg3    (vwms,vy,vz,alpha,beta,ntot,'vwm ',ifverbose)
        !call avg3    (wums,vz,vx,alpha,beta,ntot,'wum ',ifverbose)
        
        !call torque_calc(1.0,x0,.false.,.false.) ! compute wall shear
        !dragx_avg = alpha*dragx_avg + beta*dragx(iobj_wall)

        endif

        if ( (mod(istep,iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then

        time_temp = time
        time      = atime   ! Output the duration of this avg
        dtmp      = param(63)
        param(63) = 1       ! Enforce 64-bit output


        !mean fluctuation fields
        ifto=.false.;ifpo=.false.
        call opsub3 (do1,do2,do3, uavg,vavg,wavg, ubase,vbase,wbase)
        call outpost(do1,do2,do3,pr,t,'avt')
            

        call outpost2(uavg,vavg,wavg,pavg,tavg,ldimt,'avg')
        call outpost2(urms,vrms,wrms,prms,trms,ldimt,'rms')
        call outpost (uvms,vwms,wums,prms,trms,      'rm2')

        !urms2 = urms-uavg*uavg
        !vrms2 = vrms-vavg*vavg
        !wrms2 = wrms-wavg*wavg
        !uvms2 = uvms-uavg*vavg
        !vwms2 = vwms-vavg*wavg
        !wums2 = wums-uavg*wavg
        !call outpost(urms2,vrms2,wrms2,pr,t,'rm1')
        !call outpost(uvms2,vwms2,wums2,pr,t,'rm2')

        param(63) = dtmp
        atime = 0.
        time  = time_temp  ! Restore clock

        endif
        timel = time

        return
        end
c----------------------------------------------------------------------
