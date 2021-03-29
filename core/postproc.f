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
     &        ". Please use one of: 'lambda2', 'q', 'delta', or 'swirling'"
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

!     --> Element-wise velocity pseudo-gradient.
      real gije(lx1*ly1*lz1, ldim, ldim)

!     --> Point-wise symmetric and anti-symmetric parts.
      real ss(ldim, ldim), oo(ldim, ldim)
      real omega(lx1, ly1, lz1, lelv)
      real norm_a, norm_b

!     --> Miscellaneous.
      integer ie, l, i, j

!     --> Loop through the elements.
      do ie = 1, nelv

!     --> Compute the velocity pseudo-gradient.
         call comp_gije(gije, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)

!     --> Point-wise computations.
         do l = 1, nxyz

            do j = 1, ldim
               do i = 1, ldim
!     --> Compute the symmetric and antisymmetric component.
                  ss(i, j) = 0.50d0 * (gije(l, i, j) + gije(l, j, i))
                  oo(i, j) = 0.50d0 * (gije(l, i, j) - gije(l, j, i))
               enddo
            enddo

!     --> Compute the Frobenius norm
            norm_a = (norm2(ss) ** 2) ** 2
            norm_b = (norm2(oo) ** 2) ** 2

!     --> Compute omega.
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
      do ie=1,nelv              ! Compute velocity gradient tensor
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
      do ie=1,nelv              ! Compute velocity gradient tensor
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
      do ie=1,nelv              ! Compute velocity gradient tensor
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
!     Compute velocity gradient tensor
         call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
         do l=1,nxyz
            call compute_firstInv(P1,l)
            P1=-P1              !negative sign
            call compute_secondInv(Q1,l)
            call compute_thirdInv(R1,l)
            R1=-R1              !negative sign
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
c     presented as lambda_ci^2
c     
c     for the 3x3 case
c     |d11, d12, d13|  |du/dx,du/dy,du/dz|
c     D = [d_ij] = |d21, d22, d23|= |dv/dx,dv/dy,dv/dz|
c     |d31, d32, d33|  |dw/dx,dw/dy,dw/dz|
c     
c     |lambda_r,  0,       0      |
c     [d_ij]= [vr,vcr,vci]|  0,   lambda_cr, lambda_ci|[vr,vcr,vci]^-1
c     |  0,  -lambda_ci, lambda_cr|
c     
c     λ^3+P*λ^2+Qλ+R=0
c     
c     where P = -I_D   = -tr(D)
c     Q = II_D  =  0.5[P*P - tr(DD)]
c     R = -III_D = 1/3.[-P+3QP-tr(DDD)]
c     for the 2x2 case
c     
c     λ^2+Qλ+R=0
c     
c     where Q = II_D  =  0.5[P*P - tr(DD)]
c     R = -III_D = 1/3.[-P+3QP-tr(DDD)]
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
      if (if3d) then            ! 3D CASE
         do ie=1,nelv
!     Compute velocity gradient tensor
            call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
            do l=1,nxyz
               call compute_firstInv(P1,l)
               P1=-P1           !negative sign
               call compute_secondInv(Q1,l)
               call compute_thirdInv(R1,l)
               R1=-R1           !negative sign
!     Q=Q1-(P1**2)/3
               R=R1+(P1**3)*2/27-P1*Q1/3
!     Delta = (R/2)**2+(Q/3)**3
               call  cubicLambdaCi(P1, Q1, R1, lambdaCi)

               l2(l,1,1,ie) =lambdaCi
            enddo
         enddo
      elseif (ifaxis) then      ! AXISYMMETRIC CASE
         if(nid.eq.0) write(6,*)
     &        'ABORT:no compute_swirling axialsymmetric support for now'
         call exitt
      else                      ! 2D CASE
         do ie=1,nelv
!     Compute velocity gradient tensor
            call comp_gije(mygi,vx(1,1,1,ie),vy(1,1,1,ie),vz(1,1,1,ie),ie)
            do l=1,nxyz
               call compute_firstInv(P1,l)
               P1=-P1           !negative sign
!     call compute_secondInv(Q1,l)
!     Q1=0.
               call compute_thirdInv(R1,l)
               R1=-R1           !negative sign
!     Q=Q1-(P1**2)/3
!     R=R1+(P1**3)*2/27-P1*Q1/3
!     Delta = (R/2)**2+(Q/3)**3
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
c     |11(dudx),12(dudy),13(dudz)|
c     |21(dvdx),22(dvdy),23(dvdz)|
c     |31(dwdx),32(dwdy),33(dwdz)|
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lxyz=lx1*ly1*lz1
      real mygi(lxyz,ldim,ldim),B
      integer l
      common /mygrad/ mygi
      B=0.0d0
      if(if3d)then              ! 3D CASE
         B=((mygi(l,1,2)+mygi(l,2,1))**2+(mygi(l,1,3)+mygi(l,3,1))**2+(mygi(l,2,3)+mygi(l,3,2))**2)/4
      else                      ! 2D CASE
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
      if(if3d)then              ! 3D CASE
         B=((mygi(l,1,2)+mygi(l,2,1))**2+(mygi(l,1,3)+mygi(l,3,1))**2+(mygi(l,2,3)+mygi(l,3,2))**2)/4
         A=B+(mygi(l,1,1)**2+mygi(l,2,2)**2+mygi(l,3,3)**2)/2
      else                      ! 2D CASE
         B=((mygi(l,1,2)+mygi(l,2,1))**2)/4
         A=B+(mygi(l,1,1)**2+mygi(l,2,2)**2)/2
      endif
      return
      end
c---------------------------------------------------------
      subroutine compute_firstInv(a, l)
c     
c     for the 3x3 case
c     |d11, d12, d13|
c     D = [d_ij] = |d21, d22, d23|
c     |d31, d32, d33|
c     compute_firstInv returns tr(d_ij)=(d11+d22+d33)
c     
      include 'SIZE'
      include 'TOTAL'
      parameter(lxyz=lx1*ly1*lz1)
      real mygi(lxyz,ldim,ldim),a
      integer l
      common /mygrad/ mygi
      if(if3d)then              ! 3D CASE
         a=(mygi(l,1,1)+mygi(l,2,2)+mygi(l,3,3))
      elseif(ifaxis)then        ! AXISYMMETRIC CASE
         if(nid.eq.0)write(6,*)'ABORT: compute_firstInv axialsymmetric support for now'
         call exitt
      else                      ! 2D CASE
         a=(mygi(l,1,1)+mygi(l,2,2))
      endif
      return
      end
c-------------------------------------------------------------------
      subroutine compute_secondInv(a, l)
c     
c     for the 3x3 case
c     |d11, d12, d13|
c     D = [d_ij] = |d21, d22, d23|
c     |d31, d32, d33|
c     compute_SecondInv returns
c     0.5*(tr(d_ij)^2+tr(d_ij*d_ij))=-(d22*d33-d23*d32)-(d11*d22-d12*d21)-(d33*d11-d13*d31)
c     
      include 'SIZE'
      include 'TOTAL'
      parameter (lxyz=lx1*ly1*lz1)
      real mygi(lxyz,ldim,ldim)
      real a
      integer l
      common /mygrad/ mygi
      if (if3d) then            ! 3D CASE
         a=(mygi(l,2,2)*mygi(l,3,3)-mygi(l,2,3)*mygi(l,3,2))
     &        +(mygi(l,1,1)*mygi(l,2,2)-mygi(l,1,2)*mygi(l,2,1))
     &        +(mygi(l,3,3)*mygi(l,1,1)-mygi(l,1,3)*mygi(l,3,1))
      elseif (ifaxis) then      ! AXISYMMETRIC CASE
         if(nid.eq.0) write(6,*)'ABORT: compute_secondInv axialsymmetric support for now'
         call exitt
      else                      ! 2D CASE
         a=( mygi(l,1,1)+mygi(l,2,2) )*( mygi(l,1,1)+mygi(l,2,2) )
     &        -2*mygi(l,1,2)*mygi(l,2,1)-mygi(l,1,1)*mygi(l,1,1)-mygi(l,2,2)*mygi(l,2,2)
      endif

      return
      end
c-------------------------------------------------------------------
      subroutine compute_thirdInv(a, l)
c     
c     for the 3x3 case
c     |d11, d12, d13|
c     D = [d_ij] = |d21, d22, d23|
c     |d31, d32, d33|
c     
c     compute_thirdInv returns det(D)
c     =-d11*(d23*d32-d22*d33)-d12*(d21*d33-d31*d23)-d13*(d31*d22-d21*d32)
c     
      include 'SIZE'
      include 'TOTAL'
      parameter (lxyz=lx1*ly1*lz1)
      real mygi(lxyz,ldim,ldim)
      real a
      integer l
      common /mygrad/ mygi
      if (if3d) then            ! 3D CASE
         a=-mygi(l,1,1)*(mygi(l,2,3)*mygi(l,3,2)
     &        -mygi(l,2,2)*mygi(l,3,3))
     &        -mygi(l,1,2)*(mygi(l,2,1)*mygi(l,3,3)
     &        -mygi(l,3,1)*mygi(l,2,3))
     &        -mygi(l,1,3)*(mygi(l,3,1)*mygi(l,2,2)
     &        -mygi(l,2,1)*mygi(l,3,2))
      elseif (ifaxis) then      ! AXISYMMETRIC CASE
         if(nid.eq.0) write(6,*)
     &        'ABORT: compute_thirdInv axialsymmetric support for now'
         call exitt
      else                      ! 2D CASE
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
         f = sqrt(abs(d))       !sign(abs(d)**(1.0/2.0), d)
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
      subroutine nekStab_avg(ifstatis)
      include 'SIZE'
      include 'TOTAL'
      include 'AVG'

      logical, intent(in) :: ifstatis
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

!     initialization
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

!     dump freq
      iastep = param(68)
      if  (iastep.eq.0) iastep=param(15) ! same as iostep
      if  (iastep.eq.0) iastep=500

      ifverbose=.false.
      if (istep.le.10) ifverbose=.true.
      if  (mod(istep,iastep).eq.0) ifverbose=.true.

      if (atime.ne.0..and.dtime.ne.0.) then
         if(nio.eq.0) write(6,*) 'Computing statistics ...'
         beta  = dtime/atime
         alpha = 1.-beta

!     compute averages E(X) !avg(k) = alpha*avg(k) + beta*f(k)
         call avg1    (uavg,vx,alpha,beta,ntot ,'um  ',ifverbose)
         call avg1    (vavg,vy,alpha,beta,ntot ,'vm  ',ifverbose)
         call avg1    (wavg,vz,alpha,beta,ntot ,'wm  ',ifverbose)
         call avg1    (pavg,pr,alpha,beta,nto2 ,'prm ',ifverbose)
         call avg1  (tavg,t(1,1,1,1,1),alpha,beta,ntot ,'tm ',ifverbose)

         if(ifstatis)then       !compute fluctuations
!     compute averages E(X^2)
            call avg2    (urms,vx,alpha,beta,ntot ,'ums ',ifverbose)
            call avg2    (vrms,vy,alpha,beta,ntot ,'vms ',ifverbose)
            call avg2    (wrms,vz,alpha,beta,ntot ,'wms ',ifverbose)
            call avg2  (prms,pr,alpha,beta,nto2 ,'prms',ifverbose)
            call avg2  (trms,t(1,1,1,1,1),alpha,beta,ntot ,'trms',ifverbose)

!     compute averages E(X*Y)
!     call avg3    (uvms,vx,vy,alpha,beta,ntot,'uvm ',ifverbose)
!     call avg3    (vwms,vy,vz,alpha,beta,ntot,'vwm ',ifverbose)
!     call avg3    (wums,vz,vx,alpha,beta,ntot,'wum ',ifverbose)
!     if(ifheat)then
!     call avg3  (uvms,vx,t(1,1,1,1,1),alpha,beta,ntot,'utm ',ifverbose)
!     call avg3  (vwms,vy,t(1,1,1,1,1),alpha,beta,ntot,'vtm ',ifverbose)
!     call avg3  (wums,vz,t(1,1,1,1,1),alpha,beta,ntot,'wtm ',ifverbose)
!     endif

         endif

!     call torque_calc(1.0,x0,.false.,.false.) ! compute wall shear
!     dragx_avg = alpha*dragx_avg + beta*dragx(iobj_wall)

      endif

      if ( (mod(istep,iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then

         time_temp = time
         time      = atime      ! Output the duration of this avg
         dtmp      = param(63)
         param(63) = 1          ! Enforce 64-bit output


!     mean fluctuation fields
         ifto=.false.;ifpo=.false.
         call opsub3 (do1,do2,do3, uavg,vavg,wavg, ubase,vbase,wbase)
         call outpost(do1,do2,do3,pr,t,'avt')


         call outpost2(uavg,vavg,wavg,pavg,tavg,ldimt,'avg')
         call outpost2(urms,vrms,wrms,prms,trms,ldimt,'rms')
         call outpost (uvms,vwms,wums,prms,trms,      'rm2')

!     urms2 = urms-uavg*uavg
!     vrms2 = vrms-vavg*vavg
!     wrms2 = wrms-wavg*wavg
!     uvms2 = uvms-uavg*vavg
!     vwms2 = vwms-vavg*wavg
!     wums2 = wums-uavg*wavg
!     call outpost(urms2,vrms2,wrms2,pr,t,'rm1')
!     call outpost(uvms2,vwms2,wums2,pr,t,'rm2')

         param(63) = dtmp
         atime = 0.
         time  = time_temp      ! Restore clock

      endif
      timel = time

      return
      end
c----------------------------------------------------------------------
