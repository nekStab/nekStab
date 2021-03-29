c-----------------------------------------------------------------------
      subroutine SFD
!     SFD with higher-order temporal scheme and variable time-step
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt=lx1*ly1*lz1*lelt

      real, dimension(lt), save   :: qxo,qyo,qzo,
     $     vxo,vyo,vzo,
     $     uta,vta,wta,
     $     utb,vtb,wtb,
     $     utc,vtc,wtc

      real, dimension(lt)         :: do1,do2,do3,
     $     utmp,vtmp,wtmp

      real adt,bdt,cdt
      real residu,h1,l2,semi,linf,rate,residu0,cutoff,gain,frq,sig,umax
      real glmin,glmax,glsum,glsc3,tol
      integer n,i
      save n,residu0,tol

      if(uparam(4).gt.0. OR. uparam(5).gt.0) then

!     cutoff = uparam(04)
!     gain = -uparam(05)
         
         frq = uparam(04)*8*atan(1.0D0) ! transform St to omega
         sig = uparam(05)

!     Akervik 2006
!     cutoff = 0.5*frq
!     gain = -2*sig

!     Casacuberta 2018 JCP 375:481-497
         cutoff = 0.5*(sqrt(frq**2+sig**2)-sig)
         gain = -0.5*(sqrt(frq**2+sig**2)+sig)

         if(istep.eq.0)then

            n = nx1*ny1*nz1*nelt
            call oprzero(utmp,vtmp,wtmp)
            call oprzero(uta,vta,wta)
            call oprzero(utb,vtb,wtb)
            call oprzero(utc,vtc,wtc)
            residu0=0.0D0
            if(nid.eq.0)open(unit=10,file='residu.dat')
            call opcopy(qxo,qyo,qzo,vx,vy,vz)
            call opcopy(vxo,vyo,vzo,vx,vy,vz)
            tol = max(param(21), param(22))

         else
            
            call setab3(adt,bdt,cdt)
            call opcopy(utc,vtc,wtc, utb,vtb,wtb)
            call opcopy(utb,vtb,wtb, uta,vta,wta)
            call opsub3(uta,vta,wta, vxo,vyo,vzo, qxo,qyo,qzo) !uta=vxo-qxo

            call opcopy (do1,do2,do3, uta,vta,wta)
            call opcmult(do1,do2,do3, adt)
            call opcopy (utmp,vtmp,wtmp, do1,do2,do3)

            call opcopy (do1,do2,do3, utb,vtb,wtb)
            call opcmult(do1,do2,do3, bdt)
            call opadd2 (utmp,vtmp,wtmp, do1,do2,do3)

            call opcopy (do1,do2,do3, utc,vtc,wtc)
            call opcmult(do1,do2,do3, cdt)
            call opadd2 (utmp,vtmp,wtmp, do1,do2,do3)

!     call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo) ! Euler for reference 

            call opcmult(utmp,vtmp,wtmp, cutoff*dt)
            call opadd2(qxo,qyo,qzo, utmp,vtmp,wtmp)

         endif

         call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo) !fc=vo-qo
         call opcmult(do1,do2,do3, gain) !f=fc*-chi
         call opadd2 (fcx,fcy,fcz, do1,do2,do3) !FORCE HERE DO NOT COPY, ADD!
      else
         if(nid.eq.0)open(unit=10,file='residu.dat')
         if(nid.eq.0)write(6,*)' SFD in continuation mode'
      endif

      if(istep.ge.1)then

         call opsub2(vxo,vyo,vzo,vx,vy,vz)
         call normvc(h1,semi,l2,linf,vxo,vyo,vzo)
         residu = l2; rate = (residu-residu0); residu0 = residu
         call opcopy(vxo,vyo,vzo,vx,vy,vz)

         if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,rate
            write(6,"(A,2E15.7)")' SFD AB3 residu =',residu,rate
            write(6,*)
         endif

         if( istep.gt.100 .and. residu .lt. tol )then !save to disk and change flag
            if(nid.eq.0)write(6,*)' Converged base flow to:',tol

!     call compute_cfl(umax,vx,vy,vz,1.0)
!     if (nid.eq.0) write(6,*) 'CFL=0.5 dt max=',0.50d0/umax
!     if (nid.eq.0) write(6,*) 'CFL=5 dt max=',5.0d0/umax

            ifbfcv = .true.
            call bcast(ifbfcv  , lsize)
            param(63) = 1       ! Enforce 64-bit output
            call bcast(param,200*wdsize)
            call outpost(vx,vy,vz,pr,t,'BF_')
            param(63) = 0       ! Enforce 32-bit output
            call bcast(param,200*wdsize)
         endif

         if(istep.eq.nsteps)close(10)
      endif
      return
      end subroutine SFD
c-----------------------------------------------------------------------c
      subroutine BoostConv
!     boostconv core subroutine
      implicit none
      include 'SIZE'
      include 'TOTAL'
      
      real, dimension(lx1*ly1*lz1*lelt) ::  dvx,dvy,dvz
      real                              ::  residu,h1,semi,linf,rate
      real, save                        ::  residu0,tol
      logical, save :: init
      data             init /.false./

      if(mod(istep,bst_skp).eq.0)then

         if(.not.init)then
            residu = 0.0d0; rate = 0.0d0; residu0 = 0.0d0
            open(unit=10,file='residu.dat')
            tol = max(param(21), param(22))
            init=.true.
         endif

         call opsub3(dvx,dvy,dvz,vx,vy,vz,vxlag(1,1,1,1,1),vylag(1,1,1,1,1),vzlag(1,1,1,1,1)) !dv=v-vold
         call normvc(h1,semi,residu,linf,dvx,dvy,dvz); rate = (residu-residu0); residu0 = residu
         call boostconv_core(dvx,dvy,dvz)
         call opadd3(vx,vy,vz,vxlag(1,1,1,1,1),vylag(1,1,1,1,1),vzlag(1,1,1,1,1),dvx,dvy,dvz) !v=vold+dv

         if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,rate
            write(6,"(' BoostConv residu=',1pE11.4,' delta= ',1pE11.4)")residu,rate
            write(6,*)' '
         endif

         if(residu.lt.tol)then
            if(nid.eq.0)write(6,*)' Converged base flow to:',tol
            ifbfcv = .true.
            call bcast(ifbfcv  , lsize)
            param(63) = 1       ! Enforce 64-bit output
            call bcast(param,200*wdsize)
            call outpost(vx,vy,vz,pr,t,'BF_')
            param(63) = 0       ! Enforce 32-bit output
            call bcast(param,200*wdsize)
         endif

      endif

      return
      end subroutine BoostConv
c-----------------------------------------------------------------------
      subroutine boostconv_core (rbx, rby, rbz)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter   :: lt  = lx1*ly1*lz1*lelt
      integer, save        :: rot, n
      integer              :: j

      logical, save :: init
      data             init /.FALSE./

      real, dimension(bst_snp, bst_snp) :: dd
      real, dimension(bst_snp) :: cc, ccb
      real, dimension(lt) :: rbx,rby,rbz,dumx,dumy,dumz

      real, dimension(lt,bst_snp) :: q_x, q_y, q_z
      real, dimension(lt,bst_snp) :: x_x, x_y, x_z, y_x, y_y, y_z

      real :: glsc3
      n = nx1*ny1*nz1*nelt

      if (.not.init) then

         call oprzero(x_x(:,:), x_y(:,:), x_z(:,:))
         call oprzero(y_x(:,:), y_y(:,:), y_z(:,:))
         call oprzero(q_x(:,:), q_y(:,:), q_z(:,:))
         call opcopy(y_x(:,1),y_y(:,1),y_z(:,1),rbx,rby,rbz)
         call opcopy(x_x(:,1),x_y(:,1),x_z(:,1),rbx,rby,rbz)
         dd(:,:) = 1.0d0; rot = 1; init = .true.

      else

         call opsub2(y_x(:,rot),y_y(:,rot),y_z(:,rot),rbx,rby,rbz)
         call opsub2(x_x(:,rot),x_y(:,rot),x_z(:,rot),y_x(:,rot),y_y(:,rot),y_z(:,rot))
         call qr_dec(dd,q_x,q_y,q_z,y_x,y_y,y_z)

         do j = 1, bst_snp
            cc(j) = glsc3(rbx,bm1,q_x(:,j),n) + glsc3(rby,bm1,q_y(:,j),n)
            if(if3d) cc(j) = cc(j) + glsc3(rbz,bm1,q_z(:,j),n)
         enddo

         call linear_system(ccb,cc,dd,bst_snp); rot = mod(rot,bst_snp)+1
         call opcopy(y_x(:,rot),y_y(:,rot),y_z(:,rot),rbx,rby,rbz)

         do j = 1, bst_snp
            call opcopy(dumx,dumy,dumz,x_x(:,j),x_y(:,j),x_z(:,j))
            call opcmult(dumx,dumy,dumz,ccb(j))
            call opadd2(rbx,rby,rbz,dumx,dumy,dumz)
         enddo
         call opcopy(x_x(:,rot),x_y(:,rot),x_z(:,rot),rbx,rby,rbz)

      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine qr_dec (rr,q_x,q_y,q_z,x_x,x_y,x_z)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer,parameter :: lt=lx1*ly1*lz1*lelt
      integer i, j, n
      real, dimension(lt,bst_snp) :: x_x,x_y,x_z,q_x,q_y,q_z !res subspaces
      real, dimension(lt)         :: dum_x1,dum_y1,dum_z1,dum_x,dum_y,dum_z
      real, dimension(bst_snp,bst_snp) ::  rr
      real norma,glsc3

      n = nx1*ny1*nz1*nelt
      rr = 0.0d0; norma = 0.0d0
      call oprzero(Q_x(:,:), Q_y(:,:), Q_z(:,:))

      call opcopy(dum_x,dum_y,dum_z,x_x(:,1),x_y(:,1),x_z(:,1))

      norma = glsc3(dum_x,bm1,dum_x,n) + glsc3(dum_y,bm1,dum_y,n)
      if(if3d) norma = norma + glsc3(dum_z,bm1,dum_z,n)
      norma = sqrt(norma)

      call opcmult(dum_x,dum_y,dum_z,1./norma)
      call opcopy(q_x(:,1),q_y(:,1),q_z(:,1),dum_x,dum_y,dum_z)
      rr(1,1) = norma

      do j = 2,bst_snp
         call opcopy(dum_x,dum_y,dum_z,x_x(:,j),x_y(:,j),x_z(:,j))
         do i = 1,j-1

            rr(i,j)=glsc3(dum_x,bm1,q_x(:,i),n)+glsc3(dum_y,bm1,q_y(:,i),n)
            if(if3d)rr(i,j)=rr(i,j)+ glsc3(dum_z,bm1,q_z(:,i),n)

            call opcopy(dum_x1,dum_y1,dum_z1,q_x(:,i),q_y(:,i),q_z(:,i))
            call opcmult(dum_x1,dum_y1,dum_z1,rr(i,j))
            call opsub2(dum_x,dum_y,dum_z,dum_x1,dum_y1,dum_z1)

         end do

         norma = glsc3(dum_x,bm1,dum_x,n)+glsc3(dum_y,bm1,dum_y,n)
         if(if3d) norma = norma + glsc3(dum_z,bm1,dum_z,n)

         if (norma.lt.1e-60) then
            norma = 1.0d0
            q_x(:,j) = 0.0d0; q_y(:,j) = 0.0d0
            if(if3d)q_z(:,j) = 0.0d0
         else
            call opcmult(dum_x,dum_y,dum_z,1.0d0/sqrt(norma))
            call opcopy(q_x(:,j),q_y(:,j),q_z(:,j),dum_x,dum_y,dum_z)
         endif

         rr(j,j) = sqrt(norma)

      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine linear_system (outp,inp,m,size_m)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer size_m,j,k
      real m(size_m,size_m)
      real inp(size_m)
      real outp(size_m)
      outp = 0.0d0
      do j = size_m,1,-1
         outp(j) = inp(j)
         do k = j+1,size_m
            outp(j) = outp(j) - m(j,k)*outp(k)
         enddo
         outp(j) = outp(j)/m(j,j)
      enddo
      return
      end
c----------------------------------------------------------------------
