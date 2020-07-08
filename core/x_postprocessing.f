c-----------------------------------------------------------------------
        subroutine avg_all_time
        include 'SIZE'  
        include 'TOTAL' 
        include 'AVG'

        integer,parameter :: lt=lx1*ly1*lz1*lelt
        real do1(lt),do2(lt),do3(lt)

        real,save :: x0(3)
        data x0 /0.0,0.0,0.0/

        integer,save :: icalld
        data    icalld /0/

        logical ifverbose

        if (ax1.ne.lx1 .or. ay1.ne.ly1 .or. az1.ne.lz1) then
         if(nid.eq.0) write(6,*)'ABORT: wrong size of ax1,ay1,
     $      az1 in avg_all(),check SIZE!'
         call exitt
        endif
        if (ax2.ne.lx2 .or. ay2.ne.ay2 .or. az2.ne.lz2) then
         if(nid.eq.0) write(6,*)'ABORT: wrong size of ax2,ay2,
     $      az2 in avg_all(),check SIZE!'
         call exitt
        endif

        ntot  = lx1*ly1*lz1*nelv
        ntott = lx1*ly1*lz1*nelt
        nto2  = lx2*ly2*lz2*nelv

        ! initialization
        if (icalld.eq.0) then
         call rzero(prms,nto2)
         call oprzero(vwms,wums,uvms)
        endif

        dtime = time-timel
        atime = atime+dtime

        ! dump freq
        iastep = 10
        if  (iastep.eq.0) iastep=param(15)   ! same as iostep
        if  (iastep.eq.0) iastep=500

        ifverbose=.false.
        if (istep.le.10) ifverbose=.true.
        if  (mod(istep,iastep).eq.0) ifverbose=.true.

        if (atime.ne.0..and.dtime.ne.0.) then
        if(nio.eq.0) write(6,*) 'Compute statistics ...'
        beta  = dtime/atime
        alpha = 1.-beta

        ! compute averages E(X) !avg(k) = alpha*avg(k)+beta*f(k)
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
        !dragx_avg = alpha*dragx_avg+beta*dragx(iobj_wall)

        endif

        if ( (mod(istep,
     $     iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then

        time_temp = time
        time      = atime   ! Output the duration of this avg
        dtmp      = param(63)
        param(63) = 1       ! Enforce 64-bit output


        !mean fluctuation fields
        ifto=.false.;ifpo=.false.
        call opsub3 (do1,do2,do3,uavg,vavg,wavg,ubase,vbase,wbase)
        call outpost(do1,do2,do3,pr,t,'avt')
            

        call outpost2(uavg,vavg,wavg,pavg,tavg,ldimt,'avg')
        call outpost2(urms,vrms,wrms,prms,trms,ldimt,'rms')
        call outpost (uvms,vwms,wums,prms,trms,'rm2')

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
        subroutine user_stat_trnsv(lvel,dudx,dvdx,dwdx,vort)
        implicit none
        include 'SIZE'
        include 'SOLN'
        include 'INPUT'               ! if3d
        real lvel(LX1,LY1,LZ1,LELT,3) ! velocity array
        real dudx(LX1,LY1,LZ1,LELT,3) ! velocity derivatives; U
        real dvdx(LX1,LY1,LZ1,LELT,3) ! V
        real dwdx(LX1,LY1,LZ1,LELT,3) ! W
        real vort(LX1,LY1,LZ1,LELT,3) ! vorticity
        integer itmp                  ! dummy variable
        ! Velocity transformation; simple copy
        itmp = NX1*NY1*NZ1*NELV
        call copy(lvel(1,1,1,1,1),VX,itmp)
        call copy(lvel(1,1,1,1,2),VY,itmp)
        call copy(lvel(1,1,1,1,3),VZ,itmp)
        ! Derivative transformation  ! No transformation
        call gradm1(dudx(1,1,1,1,1),dudx(1,1,1,1,2),dudx(1,1,1,1,3),
     $     lvel(1,1,1,1,1))
        call gradm1(dvdx(1,1,1,1,1),dvdx(1,1,1,1,2),dvdx(1,1,1,1,3),
     $     lvel(1,1,1,1,2))
        call gradm1(dwdx(1,1,1,1,1),dwdx(1,1,1,1,2),dwdx(1,1,1,1,3),
     $     lvel(1,1,1,1,3))
        ! Vorticity
        if (IF3D) then
           call sub3(vort(1,1,1,1,1),dwdx(1,1,1,1,2),dvdx(1,1,1,1,3),
     $        itmp)
           call sub3(vort(1,1,1,1,2),dudx(1,1,1,1,3),dwdx(1,1,1,1,1),
     $        itmp)
        endif
        call sub3(vort(1,1,1,1,3),dvdx(1,1,1,1,1),dudx(1,1,1,1,2),itmp)
        return
        end
c----------------------------------------------------------------------

      subroutine stat_avg
      implicit none 
      include 'SIZE'
      include 'TOTAL' 

      if (ISTEP.eq.0) then ! Initialize statistics

         call rzero(STAT,lx1*ly1*lz1*lelt*STAT_LVAR)
         call rzero(STAT_TEMP,lx1*ly1*lz1*lelt)
         call rzero(STAT_UU,lx1*ly1*lz1*lelt)
         call rzero(STAT_VV,lx1*ly1*lz1*lelt)
         call rzero(STAT_WW,lx1*ly1*lz1*lelt)
         call rzero(STAT_PP,lx1*ly1*lz1*lelt)
         call rzero(STAT_UUU,lx1*ly1*lz1*lelt)
         call rzero(STAT_VVV,lx1*ly1*lz1*lelt)
         call rzero(STAT_WWW,lx1*ly1*lz1*lelt)
         call rzero(STAT_PPP,lx1*ly1*lz1*lelt)

         stat_atime = 0.
         stat_tstart = time

      else ! averaging

        if (mod(istep,10).eq.0) call stat_compute

        if (ifoutfld) then

        ifpo=.false.
        ifto=.true.
        
        ! fields to outpost: <u>t,<v>t,<w>t,<p>t
        call outpost(stat(1,1),stat(1,2),stat(1,3),pr,stat(1,4),'s01')
        
        ! fields to outpost: <uu>t,<vv>t,<ww>t,<pp>t
        call outpost(stat(1,5),stat(1,6),stat(1,7),pr,stat(1,8),'s02')
        
        ! fields to outpost: <uv>t,<vw>t,<uw>t,<pu>t
        call outpost(stat(1,9),stat(1,10),stat(1,11),pr,stat(1,12),'s03')

        ! fields to outpost: <pv>t,<pw>t,<pdudx>t,<pdudy>t
        call outpost(stat(1,13),stat(1,14),stat(1,15),pr,stat(1,16),'s04')
        
        ! fields to outpost: <pdudz>t,<pdvdx>t,<pdvdy>t,<pdvdz>t
        call outpost(stat(1,17),stat(1,18),stat(1,19),pr,stat(1,20),'s05')

        ! fields to outpost: <pdwdx>t,<pdwdy>t,<pdwdz>t,<uuu>t
        call outpost(stat(1,21),stat(1,22),stat(1,23),pr,stat(1,24),'s06')

        ! fields to outpost:  <vvv>t,<www>t,<uuv>t,<uuw>t
        call outpost(stat(1,25),stat(1,26),stat(1,27),pr,stat(1,28),'s07')
        
        ! fields to outpost: <vvu>t,<vvw>t,<wwu>t,<wwv>t
        call outpost(stat(1,29),stat(1,30),stat(1,31),pr,stat(1,32),'s08')
        
        ! fields to outpost:  <ppp>t,<pppp>t,<uvw>t,<uuuu>t
        call outpost(stat(1,33),stat(1,34),stat(1,35),pr,stat(1,36),'s09')
        
        ! fields to outpost: <vvvv>t,<wwww>t,<e11>t,<e22>t
        call outpost(stat(1,37),stat(1,38),stat(1,39),pr,stat(1,40),'s10')
        
        ! fields to outpost: <e33>t,<e12>t,<e13>t,<e23>t
        call outpost(stat(1,41),stat(1,42),stat(1,43),pr,stat(1,44),'s11')

      ifpo=.true.
      ifto=.true.
      stat_atime = 0.
      stat_tstart = time
        
      endif

      endif      
      
      return
      end subroutine stat_avg
c----------------------------------------------------------------------
      subroutine stat_compute
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real ltim_init
      integer npos              ! position in STAT_RUAVG
      real alpha,beta,dtime    ! time averaging parameters
      integer lnvar             ! count number of variables
      integer i                 ! loop index
      integer itmp              ! dummy variable
      real rtmp                 ! dummy variable
      integer ntot
      
      real slvel(LX1,LY1,LZ1,LELT,3),slp(LX1,LY1,LZ1,LELT)
      common /SCRMG/ slvel,slp
      real tmpvel(LX1,LY1,LZ1,LELT,3),tmppr(LX1,LY1,LZ1,LELT)
      common /SCRUZ/ tmpvel,tmppr

      real dudx(LX1,LY1,LZ1,LELT,3) ! du/dx,du/dy and du/dz
      real dvdx(LX1,LY1,LZ1,LELT,3) ! dv/dx,dv/dy and dv/dz
      real dwdx(LX1,LY1,LZ1,LELT,3) ! dw/dx,dw/dy and dw/dz
      common /SCRNS/ dudx,dvdx
      common /SCRSF/ dwdx
      real dnekclock
     
!     Calculate time span of current statistical sample
      dtime=time-stat_atime-stat_tstart
      if(nid.eq.0)write(6,*)'dtime=',dtime

!     Update total time over which the current stat file is averaged
      stat_atime=time-stat_tstart
      if(nid.eq.0)write(6,*)'stat_atime=',stat_atime

!     Time average compuated as: Accumulated=alpha*Accumulated+beta*New
      beta=dtime/stat_atime
      alpha=1.0d0-beta
      
!     Map pressure to velocity mesh
      call mappr(tmppr,PR,tmpvel(1,1,1,1,2),tmpvel(1,1,1,1,3))

!     Compute derivative tensor
      call user_stat_trnsv(tmpvel,dudx,dvdx,dwdx,slvel)
      
!     Number of points per processor in each of the 3D fields
!      ntot=lx1*ly1*lz1*lelt
      ntot  = lx1*ly1*lz1*nelv

!     Computation of statistics
      
!     <u>t
      call add2sxy(STAT(1,1),alpha,vx,beta,ntot)

!     <v>t
      call add2sxy(STAT(1,2),alpha,vy,beta,ntot)

!     <w>t
      call add2sxy(STAT(1,3),alpha,vz,beta,ntot)

!     <p>t
      call add2sxy(STAT(1,4),alpha,tmppr,beta,ntot)

c------------------------------------------------------------
      
!     <uu>t
      call col3(STAT_UU,vx,vx,ntot)
      call add2sxy(STAT(1,5),alpha,STAT_UU,beta,ntot)

!     <vv>t
         call col3(STAT_VV,vy,vy,ntot)
         call add2sxy(STAT(1,6),alpha,STAT_VV,beta,ntot)

!     <ww>t
         call col3(STAT_WW,vz,vz,ntot)
         call add2sxy(STAT(1,7),alpha,STAT_WW,beta,ntot)

!     <pp>t
         call col3(STAT_PP,tmppr,tmppr,ntot)
         call add2sxy(STAT(1,8),alpha,STAT_PP,beta,ntot)

c------------------------------------------------------------ 

!     <uv>t
         call col3(STAT_TEMP,vx,vy,ntot)
         call add2sxy(STAT(1,9),alpha,STAT_TEMP,beta,ntot)

!     <vw>t
         call col3(STAT_TEMP,vy,vz,ntot)
         call add2sxy(STAT(1,10),alpha,STAT_TEMP,beta,ntot)

!     <uw>t
         call col3(STAT_TEMP,vx,vz,ntot)
         call add2sxy(STAT(1,11),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <pu>t
         call col3(STAT_TEMP,tmppr,vx,ntot)
         call add2sxy(STAT(1,12),alpha,STAT_TEMP,beta,ntot)

!     <pv>t
         call col3(STAT_TEMP,tmppr,vy,ntot)
         call add2sxy(STAT(1,13),alpha,STAT_TEMP,beta,ntot)

!     <pw>t
         call col3(STAT_TEMP,tmppr,vz,ntot)
         call add2sxy(STAT(1,14),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------

!     <pdudx>t

         call col3(STAT_TEMP,tmppr,dudx(1,1,1,1,1),ntot)
         call add2sxy(STAT(1,15),alpha,STAT_TEMP,beta,ntot)

!     <pdudy>t
         call col3(STAT_TEMP,tmppr,dudx(1,1,1,1,2),ntot)
         call add2sxy(STAT(1,16),alpha,STAT_TEMP,beta,ntot)

!     <pdudz>t
         call col3(STAT_TEMP,tmppr,dudx(1,1,1,1,3),ntot)
         call add2sxy(STAT(1,17),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <pdvdx>t
         call col3(STAT_TEMP,tmppr,dvdx(1,1,1,1,1),ntot)
         call add2sxy(STAT(1,18),alpha,STAT_TEMP,beta,ntot)

!     <pdvdy>t
         call col3(STAT_TEMP,tmppr,dvdx(1,1,1,1,2),ntot)
         call add2sxy(STAT(1,19),alpha,STAT_TEMP,beta,ntot)

!     <pdvdz>t
         call col3(STAT_TEMP,tmppr,dvdx(1,1,1,1,3),ntot)
         call add2sxy(STAT(1,20),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <pdwdx>t
         call col3(STAT_TEMP,tmppr,dwdx(1,1,1,1,1),ntot)
         call add2sxy(STAT(1,21),alpha,STAT_TEMP,beta,ntot)

!     <pdwdy>t
         call col3(STAT_TEMP,tmppr,dwdx(1,1,1,1,2),ntot)
         call add2sxy(STAT(1,22),alpha,STAT_TEMP,beta,ntot)

!     <pdwdz>t
         call col3(STAT_TEMP,tmppr,dwdx(1,1,1,1,3),ntot)
         call add2sxy(STAT(1,23),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 
! !     UU, VV, WW
!       if (STATS3D.eq.0) then
!          itmp = LX1*LY1*LZ1*LELT*LDIM
!          call col3(tmpvel(1,1,1,1,1),slvel(1,1,1,1,1),slvel(1,1,1,1,1),itmp)
!       endif

!     <uuu>t
         call col3(STAT_UUU,STAT_UU,vx,ntot)
         call add2sxy(STAT(1,24),alpha,STAT_UUU,beta,ntot)

!     <vvv>t
         call col3(STAT_VVV,STAT_VV,vy,ntot)
         call add2sxy(STAT(1,25),alpha,STAT_VVV,beta,ntot)

!     <www>t
         call col3(STAT_WWW,STAT_WW,vz,ntot)
         call add2sxy(STAT(1,26),alpha,STAT_WWW,beta,ntot)

c------------------------------------------------------------ 

!     <uuv>t
         call col3(STAT_TEMP,STAT_UU,vy,ntot)
         call add2sxy(STAT(1,27),alpha,STAT_TEMP,beta,ntot)

!     <uuw>t
         call col3(STAT_TEMP,STAT_UU,vz,ntot)
         call add2sxy(STAT(1,28),alpha,STAT_TEMP,beta,ntot)

!     <vvu>t
         call col3(STAT_TEMP,STAT_VV,vx,ntot)
         call add2sxy(STAT(1,29),alpha,STAT_TEMP,beta,ntot)

!     <vvw>t
         call col3(STAT_TEMP,STAT_VV,vz,ntot)
         call add2sxy(STAT(1,30),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <wwu>t
         call col3(STAT_TEMP,STAT_WW,vx,ntot)
         call add2sxy(STAT(1,31),alpha,STAT_TEMP,beta,ntot)

!     <wwv>t
         call col3(STAT_TEMP,STAT_WW,vy,ntot)
         call add2sxy(STAT(1,32),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <ppp>t
         call col3(STAT_PPP,STAT_PP,tmppr,ntot)
         call add2sxy(STAT(1,33),alpha,STAT_PPP,beta,ntot)

c------------------------------------------------------------ 

!     <pppp>t
         call col3(STAT_TEMP,STAT_PPP,tmppr,ntot)
         call add2sxy(STAT(1,34),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <uvw>t
!     uv
         call col3(STAT_TEMP,vx,vy,ntot)
!     uvw
         call col2(STAT_TEMP,vz,ntot)
!     <uvw>t
         call add2sxy(STAT(1,35),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------ 

!     <uuuu>t
         call col3(STAT_TEMP,STAT_UUU,vx,ntot)
         call add2sxy(STAT(1,36),alpha,STAT_TEMP,beta,ntot)

!     <vvvv>t
         call col3(STAT_TEMP,STAT_VVV,vy,ntot)
         call add2sxy(STAT(1,37),alpha,STAT_TEMP,beta,ntot)

         !     <wwww>t
         call col3(STAT_TEMP,STAT_WWW,vz,ntot)
         call add2sxy(STAT(1,38),alpha,STAT_TEMP,beta,ntot)

c------------------------------------------------------------
      
!     <e11>t : (du/dx)^2 + (du/dy)^2 + (du/dz)^2
         call col3(tmpvel(1,1,1,1,1),dudx(1,1,1,1,1),dudx(1,1,1,1,1),
     $      3*ntot)
!     (du/dx)^2+(du/dy)^2
         call add3(STAT_TEMP,tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),ntot)
!     e11: (du/dx)^2+(du/dy)^2+(du/dz)^2
         call add2(STAT_TEMP,tmpvel(1,1,1,1,3),ntot)
!     <e11>t
         call add2sxy(STAT(1,39),alpha,STAT_TEMP,beta,ntot)

!     <e22>t: (dv/dx)^2 + (dv/dy)^2 + (dv/dz)^2

         call col3(tmpvel(1,1,1,1,1),dvdx(1,1,1,1,1),dvdx(1,1,1,1,1),
     $      3*ntot)
!     (dv/dx)^2+(dv/dy)^2
         call add3(STAT_TEMP,tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),ntot)
!     e22: (dv/dx)^2+(dv/dy)^2+(dv/dz)^2
         call add2(STAT_TEMP,tmpvel(1,1,1,1,3),ntot)
!     <e22>t
         call add2sxy(STAT(1,40),alpha,STAT_TEMP,beta,ntot)
      
c!     <e33>t: (dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2
         call col3(tmpvel(1,1,1,1,1),dwdx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $      3*ntot)
!      (dw/dx)^2+(dw/dy)^2
         call add3(STAT_TEMP,tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),ntot)
!     e33: (dw/dx)^2+(dw/dy)^2+(dw/dz)^2
         call add2(STAT_TEMP,tmpvel(1,1,1,1,3),ntot)
!     <e33>t
         call add2sxy(STAT(1,41),alpha,STAT_TEMP,beta,ntot)
      
c------------------------------------------------------------ 

!     <e12>t: (du/dx)*(dv/dx) + (du/dy)*(dv/dy) + (du/dz)*(dv/dz)
         call col3(tmpvel(1,1,1,1,1),dudx(1,1,1,1,1),dvdx(1,1,1,1,1),
     $      3*ntot)
!     du/dx*dv/dx+du/dy*dv/dy
         call add3(STAT_TEMP,tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),ntot)
!     e12: du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz
         call add2(STAT_TEMP,tmpvel(1,1,1,1,3),ntot)
!     <e12>t
         call add2sxy(STAT(1,42),alpha,STAT_TEMP,beta,ntot)

!     <e13>t: (du/dx)*(dw/dx) + (du/dy)*(dw/dy) + (du/dz)*(dw/dz)
         call col3(tmpvel(1,1,1,1,1),dudx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $      3*ntot)
!     du/dx*dw/dx+du/dy*dw/dy
         call add3(STAT_TEMP,tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),ntot)
!     e13: du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz
         call add2(STAT_TEMP,tmpvel(1,1,1,1,3),ntot)
!     <e13>t
         call add2sxy(STAT(1,43),alpha,STAT_TEMP,beta,ntot)
      
!     <e23>t: (dv/dx)*(dw/dx) + (dv/dy)*(dw/dy) + (dv/dz)*(dw/dz)
         call col3(tmpvel(1,1,1,1,1),dvdx(1,1,1,1,1),dwdx(1,1,1,1,1),
     $      3*ntot)
!     dv/dx*dw/dx+dv/dy*dw/dy
         call add3(STAT_TEMP,tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),ntot)
!     e23: dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz
         call add2(STAT_TEMP,tmpvel(1,1,1,1,3),ntot)
!     <e23>t
         call add2sxy(STAT(1,44),alpha,STAT_TEMP,beta,ntot)
     
      return
      end subroutine stat_compute

c----------------------------------------------------------------------
        subroutine stats_3d_full
        implicit none
        include 'SIZE'  
        include 'TOTAL' 
        integer,parameter                  :: lt  = lx1*ly1*lz1*lelt
        integer,parameter                  :: lt2 = lx2*ly2*lz2*lelt
      
        character*132 inputname1,hdr,pippo,pippa
        integer nfiles,ssf,sfn,ntot,ifld,j
        real nu,rho,mtimes,ftime,ltime,fdt,fttime,alpha

        real,dimension(lt) :: dum1,dum2,U,V,W,P,uu,vv,ww,uv,uw,vw,pp,
     $     ppp,pppp,uuu,vvv,www,uuv,uuw,uvv,vvw,uww,vww,uvw,dUdx,dUdy,
     $     dUdz,dVdx,dVdy,dVdz,dWdx,dWdy,dWdz,Pxx,Pyy,Pzz,Pxy,Pyz,Pxz,
     $     Dxx,Dyy,Dzz,Dxy,Dyz,Dxz,duudx,duudy,duudz,dvvdx,dvvdy,dvvdz,
     $     dwwdx,dwwdy,dwwdz,duvdx,duvdy,duvdz,duwdx,duwdy,duwdz,dvwdx,
     $     dvwdy,dvwdz,Cxx,Cyy,Czz,Cxy,Cyz,Cxz,d2uudx2,d2uudy2,d2uudz2,
     $     d2vvdx2,d2vvdy2,d2vvdz2,d2wwdx2,d2wwdy2,d2wwdz2,d2uvdx2,
     $     d2uvdy2,d2uvdz2,d2uwdx2,d2uwdy2,d2uwdz2,d2vwdx2,d2vwdy2,
     $     d2vwdz2,VDxx,VDyy,VDzz,VDxy,VDyz,VDxz,duuudx,duvvdx,duvvdy,
     $     duwwdx,duwwdz,duuvdx,duuvdy,duuwdx,duuwdz,duvwdx,duvwdy,
     $     duvwdz,dvvvdy,dvwwdy,dvwwdz,dvvwdy,dvvwdz,dwwwdz,Txx,Tyy,Tzz,
     $     Txy,Tyz,Txz,PTxx,PTyy,PTzz,PTxy,PTyz,PTxz,dpudx,dpudy,dpudz,
     $     dpvdx,dpvdy,dpvdz,dpwdx,dpwdy,dpwdz,dPdx,dPdy,dPdz,PSxx,PSyy,
     $     PSzz,PSxy,PSyz,PSxz,pdudx,pdudy,pdudz,pdvdx,pdvdy,pdvdz,
     $     pdwdx,pdwdy,pdwdz,Pixx,Piyy,Pizz,Pixy,Piyz,Pixz,Pk,Dk,Tk,VDk,
     $     Pik,Ck,Resk

        real stats_temp(lt,stat_lvar)
        real,save :: time_zero

        nfiles=stat_comp         ! Number of stat tiles
        nu=param(2)               ! Kinematic viscosity
        rho=1.0                   ! Fluid density
        mtimes=0.0                ! Starting time CHANGE
        
        ntot=lx1*ly1*lz1*lelt
        alpha=1.0
        
        call rzero(stat,lx1*ly1*lz1*lelt*stat_lvar)
        call rzero(stats_temp,lx1*ly1*lz1*lelt*stat_lvar)
  
        do sfn=1,11
           write(pippa,'(i2.2)') sfn
           ltime=mtimes           ! Time of last field
           ftime=0.0              ! Time of current field
           fttime=0.0             ! Total accumulated time
           do ssf = 1,nfiles
           
              write(pippo,'(i5.5)') ssf
              inputname1 = 's'//trim(pippa)//trim(SESSION)//'0.f'//trim(
     $           pippo)
              if(nid.eq.0)write(6,*)'Reading: ',inputname1
              call read_hdr(inputname1,
     $           ftime) ! We read header to get the times
  
              fdt=ftime-ltime     ! Time span of this field
              fttime=fttime+fdt   ! Total averaging time
           
              if(nid.eq.0)write(6,*) '**FIELD,Ts,Tf,Ta',ssf,ltime,ftime,
     $           fdt
           
              ltime=ftime         ! Update last field time
           
              !call load_field(inputname1)
              call load_fld(inputname1)
  
              call add2sxy(STATS_TEMP(1,4*(sfn-1)+1),alpha,vx,fdt,ntot)
              call add2sxy(STATS_TEMP(1,4*(sfn-1)+2),alpha,vy,fdt,ntot)
              call add2sxy(STATS_TEMP(1,4*(sfn-1)+3),alpha,vz,fdt,ntot)
              call add2sxy(STATS_TEMP(1,4*(sfn-1)+4),alpha,t,fdt,ntot)
  
           enddo
        enddo
  
  !     Reorder fields. No change for the first 26 variables
        do ifld=1,26
           do j=1,ntot
              STAT(j,ifld) = STATS_TEMP(j,ifld)
           enddo
        enddo

        !     Various shifts in position
        do ifld=27,32
                do j=1,ntot
                STAT(j,ifld+1) = STATS_TEMP(j,ifld)
                enddo
        enddo
        
        do j=1,ntot
                STAT(j,27) = STATS_TEMP(j,33)
        enddo

        do j=1,ntot
                STAT(j,38) = STATS_TEMP(j,34)
        enddo
        
        do ifld=35,38
                do j=1,ntot
                STAT(j,ifld-1) = STATS_TEMP(j,ifld)
                enddo
        enddo
        
        !     No change for the last variables
        do ifld=39,44
                do j=1,ntot
                STAT(j,ifld) = STATS_TEMP(j,ifld)
                enddo
        enddo
        
        !     Divide by the total averaging time
        do ifld=1,STAT_LVAR
                do j=1,ntot
                STAT(j,ifld) = STAT(j,ifld)/fttime
                enddo
        enddo


        ! Mean velocities. Tensors of Rank 1.
        do j=1,ntot
                U(j)=STAT(j,1)
                V(j)=STAT(j,2)
                W(j)=STAT(j,3)
        enddo

        ! Reynolds-stress tensor. Tensor of Rank 2.
        do j=1,ntot
                uu(j)=STAT(j,5)-U(j)*U(j)
                vv(j)=STAT(j,6)-V(j)*V(j)
                ww(j)=STAT(j,7)-W(j)*W(j)
                uv(j)=STAT(j,9)-U(j)*V(j)
                uw(j)=STAT(j,11)-U(j)*W(j)
                vw(j)=STAT(j,10)-V(j)*W(j)
        enddo

        ! Mean,RMS,skewness and flatness of pressure
        do j=1,ntot
                P(j)=STAT(j,4)
                pp(j)=STAT(j,8)-P(j)*P(j)
                ppp(j)=STAT(j,27)-3*P(j)*pp(j)-P(j)*P(j)*P(j)
                pppp(j)=STAT(j,38)-4*P(j)*ppp(j)-6*P(j)*P(j)*pp(j)-P(j)*
     $             P(j)*P(j)*P(j)
        enddo

        !Skewness tensor. Tensor of Rank 3. Tensor form:
        ![ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ]
        ![ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw vww ]
        ![ uuw uvw uww ] [ uvw vvw vww ] [ uww vww www ]
        do j=1,ntot
         uuu(j)=STAT(j,24)-3*U(j)*uu(j)-U(j)*U(j)*U(j)
         vvv(j)=STAT(j,25)-3*V(j)*vv(j)-V(j)*V(j)*V(j)
         www(j)=STAT(j,26)-3*W(j)*ww(j)-W(j)*W(j)*W(j)
         uuv(j)=STAT(j,28)-2*U(j)*uv(j)-V(j)*uu(j)-U(j)*U(j)*V(j)
         uuw(j)=STAT(j,29)-2*U(j)*uw(j)-W(j)*uu(j)-U(j)*U(j)*W(j)
         uvv(j)=STAT(j,30)-2*V(j)*uv(j)-U(j)*vv(j)-V(j)*V(j)*U(j)
         vvw(j)=STAT(j,31)-2*V(j)*vw(j)-W(j)*vv(j)-V(j)*V(j)*W(j)
         uww(j)=STAT(j,32)-2*W(j)*uw(j)-U(j)*ww(j)-W(j)*W(j)*U(j)
         vww(j)=STAT(j,33)-2*W(j)*vw(j)-V(j)*ww(j)-W(j)*W(j)*V(j)
         uvw(j)=STAT(j,34)-U(j)*vw(j)-V(j)*uw(j)-W(j)*uv(j)-U(j)*V(j)*
     $      W(j)
        enddo

        !Velocity gradient tensor. Tensor of Rank 2.
        call gradm1(dUdx,dUdy,dUdz,U)
        call gradm1(dVdx,dVdy,dVdz,V)
        call gradm1(dWdx,dWdy,dWdz,W)

        !Production tensor. Tensor of Rank 2.
        do j=1,ntot
         Pxx(j)=-2*(uu(j)*dUdx(j)+uv(j)*dUdy(j)+uw(j)*dUdz(j))
         Pyy(j)=-2*(uv(j)*dVdx(j)+vv(j)*dVdy(j)+vw(j)*dVdz(j))
         Pzz(j)=-2*(uw(j)*dWdx(j)+vw(j)*dWdy(j)+ww(j)*dWdz(j))
         Pxy(j)=-(uu(j)*dVdx(j)+uv(j)*dVdy(j)+uv(j)*dUdx(j)+vv(j)*
     $      dUdy(j)+uw(j)*dVdz(j)+vw(j)*dUdz(j))
         Pxz(j)=-(uu(j)*dWdx(j)+uv(j)*dWdy(j)+uw(j)*dUdx(j)+vw(j)*
     $      dUdy(j)+uw(j)*dWdz(j)+ww(j)*dUdz(j))
         Pyz(j)=-(uv(j)*dWdx(j)+vv(j)*dWdy(j)+uw(j)*dVdx(j)+vw(j)*
     $      dVdy(j)+vw(j)*dWdz(j)+ww(j)*dVdz(j))
        enddo

        !     Dissipation tensor. Tensor of Rank 2.
        do j=1,ntot
         Dxx(j)=-2*nu*(STAT(j,39)-dUdx(j)*dUdx(j)-dUdy(j)*dUdy(j)-
     $      dUdz(j)*dUdz(j))
         Dyy(j)=-2*nu*(STAT(j,40)-dVdx(j)*dVdx(j)-dVdy(j)*dVdy(j)-
     $      dVdz(j)*dVdz(j))
         Dzz(j)=-2*nu*(STAT(j,41)-dWdx(j)*dWdx(j)-dWdy(j)*dWdy(j)-
     $      dWdz(j)*dWdz(j))
         Dxy(j)=-2*nu*(STAT(j,42)-dUdx(j)*dVdx(j)-dUdy(j)*dVdy(j)-
     $      dUdz(j)*dVdz(j))
         Dxz(j)=-2*nu*(STAT(j,43)-dUdx(j)*dWdx(j)-dUdy(j)*dWdy(j)-
     $      dUdz(j)*dWdz(j))
         Dyz(j)=-2*nu*(STAT(j,44)-dVdx(j)*dWdx(j)-dVdy(j)*dWdy(j)-
     $      dVdz(j)*dWdz(j))
        enddo

        ! Derivatives of the Reynolds-stress tensor components
        call gradm1(duudx,duudy,duudz,uu)
        call gradm1(dvvdx,dvvdy,dvvdz,vv)
        call gradm1(dwwdx,dwwdy,dwwdz,ww)
        call gradm1(duvdx,duvdy,duvdz,uv)
        call gradm1(duwdx,duwdy,duwdz,uw)
        call gradm1(dvwdx,dvwdy,dvwdz,vw)

        ! Mean convection tensor. Tensor of Rank 2.
        do j=1,ntot
         Cxx(j)=U(j)*duudx(j)+V(j)*duudy(j)+W(j)*duudz(j)
         Cyy(j)=U(j)*dvvdx(j)+V(j)*dvvdy(j)+W(j)*dvvdz(j)
         Czz(j)=U(j)*dwwdx(j)+V(j)*dwwdy(j)+W(j)*dwwdz(j)
         Cxy(j)=U(j)*duvdx(j)+V(j)*duvdy(j)+W(j)*duvdz(j)
         Cxz(j)=U(j)*duwdx(j)+V(j)*duwdy(j)+W(j)*duwdz(j)
         Cyz(j)=U(j)*dvwdx(j)+V(j)*dvwdy(j)+W(j)*dvwdz(j)
        enddo
        
        ! Second derivatives of the Reynolds-stress tensor components
        call gradm1(d2uudx2,dum1,dum2,duudx)
        call gradm1(dum1,d2uudy2,dum2,duudy)
        call gradm1(dum1,dum2,d2uudz2,duudz)
        call gradm1(d2vvdx2,dum1,dum2,dvvdx)
        call gradm1(dum1,d2vvdy2,dum2,dvvdy)
        call gradm1(dum1,dum2,d2vvdz2,dvvdz)
        call gradm1(d2wwdx2,dum1,dum2,dwwdx)
        call gradm1(dum1,d2wwdy2,dum2,dwwdy)
        call gradm1(dum1,dum2,d2wwdz2,dwwdz)
        call gradm1(d2uvdx2,dum1,dum2,duvdx)
        call gradm1(dum1,d2uvdy2,dum2,duvdy)
        call gradm1(dum1,dum2,d2uvdz2,duvdz)
        call gradm1(d2uwdx2,dum1,dum2,duwdx)
        call gradm1(dum1,d2uwdy2,dum2,duwdy)
        call gradm1(dum1,dum2,d2uwdz2,duwdz)
        call gradm1(d2vwdx2,dum1,dum2,dvwdx)
        call gradm1(dum1,d2vwdy2,dum2,dvwdy)
        call gradm1(dum1,dum2,d2vwdz2,dvwdz)
        
        ! Viscous diffusion tensor. Tensor of Rank 2.
        do j=1,ntot
         VDxx(j)=nu*(d2uudx2(j)+d2uudy2(j)+d2uudz2(j))
         VDyy(j)=nu*(d2vvdx2(j)+d2vvdy2(j)+d2vvdz2(j))
         VDzz(j)=nu*(d2wwdx2(j)+d2wwdy2(j)+d2wwdz2(j))
         VDxy(j)=nu*(d2uvdx2(j)+d2uvdy2(j)+d2uvdz2(j))
         VDxz(j)=nu*(d2uwdx2(j)+d2uwdy2(j)+d2uwdz2(j))
         VDyz(j)=nu*(d2vwdx2(j)+d2vwdy2(j)+d2vwdz2(j))
        enddo
        
        ! Derivatives of the triple-product terms
        call gradm1(duuudx,dum1,dum2,uuu)
        call gradm1(duvvdx,duvvdy,dum2,uvv)
        call gradm1(duwwdx,dum1,duwwdz,uww)
        call gradm1(duuvdx,duuvdy,dum2,uuv)
        call gradm1(duuwdx,dum1,duuwdz,uuw)
        call gradm1(duvwdx,duvwdy,duvwdz,uvw)
        call gradm1(dum1,dvvvdy,dum2,vvv)
        call gradm1(dum1,dvwwdy,dvwwdz,vww)
        call gradm1(dum1,dvvwdy,dvvwdz,vvw)
        call gradm1(dum1,dum2,dwwwdz,www)

        ! Turbulent transport tensor. Tensor of Rank 2.
        do j=1,ntot
         Txx(j)=-(duuudx(j)+duuvdy(j)+duuwdz(j))
         Tyy(j)=-(duvvdx(j)+dvvvdy(j)+dvvwdz(j))
         Tzz(j)=-(duwwdx(j)+dvwwdy(j)+dwwwdz(j))
         Txy(j)=-(duuvdx(j)+duvvdy(j)+duvwdz(j))
         Txz(j)=-(duuwdx(j)+duvwdy(j)+duwwdz(j))
         Tyz(j)=-(duvwdx(j)+dvvwdy(j)+dvwwdz(j))
        enddo

        ! Derivatives of the pressure-velocity products
        call gradm1(dpudx,dpudy,dpudz,STAT(1,12))
        call gradm1(dpvdx,dpvdy,dpvdz,STAT(1,13))
        call gradm1(dpwdx,dpwdy,dpwdz,STAT(1,14))

        ! Derivatives of the mean pressure field
        call gradm1(dPdx,dPdy,dPdz,P)
        
        ! Pressure transport tensor. Tensor of Rank 2.
        do j=1,ntot
         dpudx(j)=dpudx(j)-P(j)*dUdx(j)-U(j)*dPdx(j)
         dpvdx(j)=dpvdx(j)-P(j)*dVdx(j)-V(j)*dPdx(j)
         dpwdx(j)=dpwdx(j)-P(j)*dWdx(j)-W(j)*dPdx(j)
         dpudy(j)=dpudy(j)-P(j)*dUdy(j)-U(j)*dPdy(j)
         dpvdy(j)=dpvdy(j)-P(j)*dVdy(j)-V(j)*dPdy(j)
         dpwdy(j)=dpwdy(j)-P(j)*dWdy(j)-W(j)*dPdy(j)
         dpudz(j)=dpudz(j)-P(j)*dUdz(j)-U(j)*dPdz(j)
         dpvdz(j)=dpvdz(j)-P(j)*dVdz(j)-V(j)*dPdz(j)
         dpwdz(j)=dpwdz(j)-P(j)*dWdz(j)-W(j)*dPdz(j)
        enddo

        do j=1,ntot
         PTxx(j)=-2.d0/rho*dpudx(j)
         PTyy(j)=-2.d0/rho*dpvdy(j)
         PTzz(j)=-2.d0/rho*dpwdz(j)
         PTxy(j)=-1.d0/rho*(dpudy(j)+dpvdx(j))
         PTxz(j)=-1.d0/rho*(dpudz(j)+dpwdx(j))
         PTyz(j)=-1.d0/rho*(dpvdz(j)+dpwdy(j))
        enddo

        ! Pressure strain tensor. Tensor of Rank 2.
        do j=1,ntot
         pdudx(j)=STAT(j,15)-P(j)*dUdx(j)
         pdudy(j)=STAT(j,16)-P(j)*dUdy(j)
         pdudz(j)=STAT(j,17)-P(j)*dUdz(j)
         pdvdx(j)=STAT(j,18)-P(j)*dVdx(j)
         pdvdy(j)=STAT(j,19)-P(j)*dVdy(j)
         pdvdz(j)=STAT(j,20)-P(j)*dVdz(j)
         pdwdx(j)=STAT(j,21)-P(j)*dWdx(j)
         pdwdy(j)=STAT(j,22)-P(j)*dWdy(j)
         pdwdz(j)=STAT(j,23)-P(j)*dWdz(j)
        enddo

        do j=1,ntot
         PSxx(j)=-2.d0/rho*pdudx(j)
         PSyy(j)=-2.d0/rho*pdvdy(j)
         PSzz(j)=-2.d0/rho*pdwdz(j)
         PSxy(j)=-1.d0/rho*(pdudy(j)+pdvdx(j))
         PSxz(j)=-1.d0/rho*(pdudz(j)+pdwdx(j))
         PSyz(j)=-1.d0/rho*(pdvdz(j)+pdwdy(j))
        enddo

        ! Velocity-pressure-gradient tensor. Tensor of Rank 2.
        do j=1,ntot
         Pixx(j)=PTxx(j)-PSxx(j)
         Piyy(j)=PTyy(j)-PSyy(j)
         Pizz(j)=PTzz(j)-PSzz(j)
         Pixy(j)=PTxy(j)-PSxy(j)
         Pixz(j)=PTxz(j)-PSxz(j)
         Piyz(j)=PTyz(j)-PSyz(j)
        enddo

        ! Calculation of TKE budget
        do j=1,ntot
         Pk(j)=0.5d0*(Pxx(j)+Pyy(j)+Pzz(j))
         Dk(j)=0.5d0*(Dxx(j)+Dyy(j)+Dzz(j))
         Tk(j)=0.5d0*(Txx(j)+Tyy(j)+Tzz(j))
         VDk(j)=0.5d0*(VDxx(j)+VDyy(j)+VDzz(j))
         Pik(j)=0.5d0*(Pixx(j)+Piyy(j)+Pizz(j))
         Ck(j)=0.5d0*(Cxx(j)+Cyy(j)+Czz(j))
         Resk(j)=Pk(j)+Dk(j)+Tk(j)+VDk(j)+Pik(j)-Ck(j)
        enddo

        ifpo=.false.
        ifto=.false.


        call whereyouwant('a01',1)
        call outpost(U,V,W,pr,uu,'a01')
        call outpost(vv,ww,uv,pr,uw,'a02')
        call outpost(vw,P,pp,pr,ppp,'a03')
        call outpost(pppp,uuu,vvv,pr,www,'a04')
        call outpost(uuv,uuw,uvv,pr,vvw,'a05')
        call outpost(uww,vww,uvw,pr,Pxx,'a06')
        call outpost(Pyy,Pzz,Pxy,pr,Pxz,'a07')
        call outpost(Pyz,Dxx,Dyy,pr,Dzz,'a08')
        call outpost(Dxy,Dxz,Dyz,pr,Txx,'a09')
        call outpost(Tyy,Tzz,Txy,pr,Txz,'a10')
        call outpost(Tyz,VDxx,VDyy,pr,VDzz,'a11')
        call outpost(VDxy,VDxz,VDyz,pr,Pixx,'a12')
        call outpost(Piyy,Pizz,Pixy,pr,Pixz,'a13')
        call outpost(Piyz,Cxx,Cyy,pr,Czz,'a14')
        call outpost(Cxy,Cxz,Cyz,pr,Pk,'a15')
        call outpost(Dk,Tk,VDk,pr,Pik,'a16')
        call outpost(Ck,Resk,PTxx,pr,PTyy,'a17')
        call outpost(PTzz,PTxy,PTxz,pr,PTyz,'a18')
        call outpost(PSxx,PSyy,PSzz,pr,PSxy,'a19')
        call outpost(PSxz,PTyz,dUdx,pr,dUdy,'a20')
        call outpost(dUdz,dVdx,dVdy,pr,dVdz,'a21')
        call outpost(dWdx,dWdy,dWdz,pr,Tk,'a22')
        
        ifpo=.true.
        ifto=.true.
        
        return
        end

        ! subroutine load_field(field)
        ! implicit none
        ! include 'SIZE'  
        ! include 'TOTAL' 

        ! character*132 field    
        ! call load_fld(field)

        ! return
        ! end

        subroutine read_hdr(field,mtimee)
        implicit none
        include 'SIZE'  
        include 'TOTAL'
          
        character*132 field,hdr,fmt1
        character*10 tmpchar
        integer twdsize,mnelx,mnely,mnelz,nelo,isteps,fid0,nfileoo,
     $     stat_gnum
        real mtimee
          
        open(unit=33,file=field,form='unformatted')
        read(33) hdr
        close(33)
          
        fmt1 = '(1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,1x,i9,
     $     1x,i6,1x,i6,1x,10a)'
          
        read(hdr,fmt1)twdsize,mnelx,mnely,mnelz,nelo,stat_gnum,mtimee,
     $     isteps,fid0,nfileoo,tmpchar
          
        return
        end