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
