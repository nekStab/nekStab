c-----------------------------------------------------------------------c
      subroutine tdf            !Time-delayed Feedback
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      
      real, dimension(lv) :: do1,do2,do3
      real                :: h1,semi,l2,linf,rate,tol
      real, save          :: residu0, gain, porbit
      integer, save       :: i,norbit

      real vort(lv,3),wo1(lv),wo2(lv) ! to outpost vorticity BFV
      common /ugrad/ vort,wo1,wo2
      
      logical, save :: init
      data             init /.FALSE./

      if (.not.init) then
         
         if(nid.eq.0)write(6,*)'   Target frequency specified in uparam(5)=',uparam(5)
         porbit = 1.0d0/uparam(5) ! user input from .usr -> also used on the inflow BC
         
         call compute_cfl(ctarg, vx, vy, vz, 1.0d0) ! ctarg contains the sum ( ux_i / dx_i )
         dt = param(26) / ctarg ! dt given target CFL
         norbit = ceiling(porbit / dt) ! computing a safe value of norbit
         if(nid.eq.0)write(6,*)' Computing norbit=',norbit    
         
         dt = porbit / norbit   ! reducing dt to match forced orbit to machine accuracy   
         param(12) = dt

         if(nid.eq.0)write(6,*)' Adjusting timeStep dt=',dt
         call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
         if(nid.eq.0)write(6,*)' veryfing current CFL and target=',ctarg,param(26)
         param(12) = -abs(param(12))       
         
         gain = -0.04432D0*8*atan(1.0D0)/porbit ! Theoretical optimal feedback parameter see reference
         
         if(nid .eq. 0) write(6,*) 'Allocating TDF orbit with nsteps:',norbit,norbit*dt
         allocate(uor(lv, norbit), vor(lv, norbit))
         if(if3d)then
            allocate(wor(lv,norbit))
         else
            allocate(wor(1,1))
         endif
         call oprzero(uor(:,:), vor(:,:), wor(:,:))
         
         if(ifheat)then
            allocate(tor(lv,norbit))
            tor(:,:)=0.0d0
         endif
         
         rate = 0.0d0; residu0 = 0.0d0
         open(unit=10,file='residu.dat')
         
         init = .true.
         
      else

         if(istep.le.norbit)then !t<T->save solutions
            
            if(nid.eq.0)write(6,*)' Storing initial solution in memory:',istep,'/',norbit
            call opcopy(uor(:,istep),vor(:,istep),wor(:,istep),vx,vy,vz)
            if(ifheat)call copy(tor(1,istep),t(1,1,1,1,1),nx1*ny1*nz1*nelv)

         else                   !t>T->compute forcing !f(t)= - \Lambda * 2*pi*St * ( u(t) - u(t-T) )

            call opsub3 (do1,do2,do3, vx,vy,vz, uor(:,1),vor(:,1),wor(:,1)) !ub=v-vold
            if(istep.gt.norbit+1) call normvc(h1,semi,l2,linf,do1,do2,do3); rate=(l2-residu0); residu0=l2
            
            call opcmult(do1,do2,do3, gain) !f=fc*-chi
            call opadd2 (fcx,fcy,fcz, do1,do2,do3) !FORCE HERE DO NOT COPY, ADD!

            do i=1,norbit-1     !discard the i=1 solution
               
               uor(:,i)=uor(:,i+1)
               vor(:,i)=vor(:,i+1)
               if (if3d)   wor(:,i)=wor(:,i+1)
               if (ifheat) tor(:,i)=tor(:,i+1)
               
            enddo               !store the last one
            call opcopy(uor(1,norbit),vor(1,norbit),wor(1,norbit),vx,vy,vz) !store solution
            if (ifheat)call copy(tor(1,norbit),t(1,1,1,1,1),nx1*ny1*nz1*nelv)
            
            if(nid.eq.0)then
               write(10,"(3E15.7)")time,l2,rate
               write(6,"(' TDF residu=',1pE11.4,' rate of change= ',1pE11.4)")l2,rate
               write(6,*)' '
            endif
            
            tol =  max(param(21), param(22))
            if(l2 .gt. 0.0d0 .and. l2 .lt. tol)then
               if(nid.eq.0)write(6,*)' Converged base flow to:',tol
               ifbfcv = .true.
               call bcast(ifbfcv  , lsize)
               param(63) = 1    ! Enforce 64-bit output
               call bcast(param,200*wdsize)
               call outpost(vx,vy,vz,pr,t,'BF_')
               param(63) = 0    ! Enforce 32-bit output
               call bcast(param,200*wdsize)
               
               if(ifvor)then    ! outpost vorticity
                  call comp_vort3(vort,wo1,wo2,vx,vy,vz);ifto=.false.;ifpo=.false.
                  call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t,'BFV')
               endif
            endif
            
         endif                  ! else
      endif                     ! not init

      return
      end subroutine TDF
c-----------------------------------------------------------------------
      subroutine SFD
!     SFD with higher-order temporal scheme and variable time-step
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real, dimension(lv), save :: qxo,qyo,qzo,vxo,vyo,vzo,uta,vta,wta,utb,vtb,wtb,utc,vtc,wtc
      real, dimension(lv) :: do1,do2,do3,utmp,vtmp,wtmp
      real adt,bdt,cdt,desired_tolerance
      real residu,h1,l2,semi,linf,residu0,cutoff,gain,frq,sig,umax,rate
      real glmin,glmax,glsum,glsc3,tol
      integer i
      real vort(lv,3),wo1(lv),wo2(lv) ! to outpost vorticity BFV
      common /ugrad/ vort,wo1,wo2

      save residu0,desired_tolerance
      logical, save :: tic,toc
      data             tic /.true./
      data             toc /.true./

      if(uparam(5).gt.0) then
         
         frq = abs(uparam(04))*8*atan(1.0d0) ! transform St to omega
         sig = abs(uparam(05))

         if(uparam(4).gt.0)then ! Akervik 2006
            
            cutoff = 0.5*frq
            gain  = -2*sig
            if(nid.eq.0)write(6,*)' Akervik  cutoff,gain:',cutoff,gain

         else                   ! Casacuberta 2018 JCP 375:481-497
            
            cutoff = 0.5*(sqrt(frq**2+sig**2)-sig)
            gain  = -0.5*(sqrt(frq**2+sig**2)+sig)
            if(nid.eq.0)write(6,*)' Casacub. cutoff,gain:',cutoff,gain
            
         endif

         tol = max(param(21), param(22))
         if(istep.eq.0)then

            n = nx1*ny1*nz1*nelv
            
            call oprzero(utmp,vtmp,wtmp)
            call oprzero(uta,vta,wta)
            call oprzero(utb,vtb,wtb)
            call oprzero(utc,vtc,wtc)
            residu0=0.0D0
            if(nid.eq.0)open(unit=10,file='residu.dat')
            call opcopy(qxo,qyo,qzo,vx,vy,vz)
            call opcopy(vxo,vyo,vzo,vx,vy,vz)
            desired_tolerance = max(param(21),param(22))

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
         if(nid.eq.0)then
            open(unit=10,file='residu.dat')
            write(6,*)' SFD in continuation mode'
         endif
         if(istep.eq.0)desired_tolerance = max(param(21),param(22))
      endif

      if(istep.ge.1)then

         call opsub2(vxo,vyo,vzo,vx,vy,vz)
         call normvc(h1,semi,l2,linf,vxo,vyo,vzo)
         residu = l2; rate = (residu-residu0)/dt; residu0 = residu
         call opcopy(vxo,vyo,vzo,vx,vy,vz)

         tol = 0.0d0
         if(nid.eq.0)then
            write(10,"(4E15.7)")time,residu,rate,param(21)
            write(6,"(A,3E15.7)")' SFD AB3 residu =',residu,rate
            write(6,*)
            tol = residu
         endif
         if ( ifdyntol .and. mod(istep,20)==0 ) call set_solv_tole(abs(tol)/20)           

         if( istep.gt.100 .and. residu .lt. desired_tolerance )then !save to disk and change flag
            if(nid.eq.0)write(6,*)' Converged base flow to:',desired_tolerance
            ifbfcv = .true.
            call bcast(ifbfcv  , lsize)
            param(63) = 1       ! Enforce 64-bit output
            call bcast(param,200*wdsize)
            call outpost(vx,vy,vz,pr,t,'BF_')
            param(63) = 0       ! Enforce 32-bit output
            call bcast(param,200*wdsize)
            
            if(ifvor)then       ! outpost vorticity
               call comp_vort3(vort,wo1,wo2,vx,vy,vz);ifto=.false.;ifpo=.false.
               call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t,'BFV')
            endif
            
         endif

         if(istep.eq.nsteps)close(10)
      endif
      return
      end subroutine SFD
c-----------------------------------------------------------------------c
      subroutine spec_tole_sfd(i)
!     Subroutine to progressively tight tolerances to maximise computational time
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer :: i
      real, save :: desired_tolerance
      logical, save :: init
      data             init /.false./

      if(.not.init)then         ! save user specified tolerance
         desired_tolerance = max(param(21),param(22))
         init = .true.
      endif
!     if restarting we need to move i previous iteration to match the tolerances of the previous case!
      if (uparam(2).gt.0)i = i + 1
      if     (i == 1 .and. desired_tolerance.le.1e-6) then
         call set_solv_tole(1e-6)
      elseif (i == 2 .and. desired_tolerance.le.1e-7) then
         call set_solv_tole(1e-7)
      elseif (i == 3 .and. desired_tolerance.le.1e-8) then
         call set_solv_tole(1e-8)
      elseif (i == 4 .and. desired_tolerance.le.1e-9) then
         call set_solv_tole(1e-9)
      elseif (i == 5 .and. desired_tolerance.le.1e-10) then
         call set_solv_tole(1e-10)
      elseif (i == 6 .and. desired_tolerance.le.1e-11) then
         call set_solv_tole(1e-11)
      elseif (i == 7 .and. desired_tolerance.le.1e-12) then
         call set_solv_tole(1e-12)
      elseif (i == 8 .and. desired_tolerance.le.1e-13) then
         call set_solv_tole(1e-13)
      else
         call set_solv_tole(desired_tolerance)
      endif
      return
      end subroutine spec_tole_SFD
c-----------------------------------------------------------------------c
      subroutine BoostConv
!     boostconv core subroutine
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      
      real, dimension(lv) ::  dvx,dvy,dvz
      real                              ::  residu,h1,semi,linf,rate,tol
      real, save                        ::  residu0
      logical, save :: init
      data             init /.false./

      if(mod(istep,bst_skp).eq.0)then

         if(.not.init)then
            residu = 0.0d0; rate = 0.0d0; residu0 = 0.0d0
            open(unit=10,file='residu.dat')
            init=.true.
         endif

         call opsub3(dvx,dvy,dvz,vx,vy,vz,vxlag(1,1,1,1,1),vylag(1,1,1,1,1),vzlag(1,1,1,1,1)) !dv=v-vold
         call normvc(h1,semi,residu,linf,dvx,dvy,dvz); rate = (residu-residu0)/dt; residu0 = residu
         call boostconv_core(dvx,dvy,dvz)
         call opadd3(vx,vy,vz,vxlag(1,1,1,1,1),vylag(1,1,1,1,1),vzlag(1,1,1,1,1),dvx,dvy,dvz) !v=vold+dv

         if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,rate
            write(6,"(' BoostConv residu=',1pE11.4,' delta= ',1pE11.4)")residu,rate
            write(6,*)' '
         endif

         tol = max(param(21), param(22))
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
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, save        :: rot
      integer              :: j

      logical, save :: init
      data             init /.FALSE./

      real, allocatable, save, dimension(:)   :: cc, ccb
      real, allocatable, save, dimension(:, :) :: dd
      real, allocatable, save, dimension(:, :) :: q_x, q_y, q_z
      real, allocatable, save, dimension(:, :) :: x_x, x_y, x_z, y_x, y_y, y_z

      real, dimension(lv) :: rbx,rby,rbz,dumx,dumy,dumz

      real :: glsc3
      n = nx1*ny1*nz1*nelv

      if (.not.init) then

         allocate(cc(bst_snp),ccb(bst_snp),dd(bst_snp, bst_snp))
         allocate(q_x(lv,bst_snp),q_y(lv,bst_snp),q_z(lv,bst_snp))
         allocate(x_x(lv,bst_snp),x_y(lv,bst_snp),x_z(lv,bst_snp))
         allocate(y_x(lv,bst_snp),y_y(lv,bst_snp),y_z(lv,bst_snp))
         
         if(nid.eq.0)write(6,*)'Allocating BoostConv variables with:',bst_snp
         if(nid.eq.0)write(6,*)'                     skipping every:',bst_skp

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
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer i, j
      real, dimension(lv,bst_snp) :: x_x,x_y,x_z,q_x,q_y,q_z !res subspaces
      real, dimension(lv)         :: dum_x1,dum_y1,dum_z1,dum_x,dum_y,dum_z
      real, dimension(bst_snp,bst_snp) ::  rr
      real norma,glsc3

      n = nx1*ny1*nz1*nelv
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
