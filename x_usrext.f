c-----------------------------------------------------------------------
      subroutine mycomment
      include 'SIZE'
      include 'TOTAL'
      real*8,save :: eetime0,eetime1,eetime2
      data           eetime0,eetime1,eetime2 /0.0d0, 0.0d0, 0.0d0/
      real, save :: deltatime
      real telapsed,tpernondt,tmiss

      if (courno.gt.2.0d0) then
        if (nio.eq.0)then
          write(6,*)' '
          write(6,*)'    cfl > 2. stopping'
          write(6,*)' '
        endif
        lastep = 1
      endif

      if (nio.ne.0) return

      if (eetime0.eq.0.0 .and. istep.eq.1)then
        eetime0=dnekclock()
        deltatime=time
      endif

      eetime1=eetime2
      eetime2=dnekclock()

      if (istep.gt.0 .and. lastep.eq.0 .and. iftran) then

        ttime_stp = eetime2-eetime1 ! time per timestep
        ttime     = eetime2-eetime0 ! sum of all timesteps

        if(istep.eq.1)then
          ttime_stp = 0.0d0; ttime = 0.0d0
        endif

        if (mod(istep,5).eq.0) then

          telapsed = ttime/3600.0d0
          tpernondt = (ttime/(time-deltatime))
          tmiss = (param(10)-time)*tpernondt/3600.0d0

          write(6,*)' '
          write(6,103)1.d0/param(2)
          write(6,102)ttime/istep,int(((ttime/istep)-ttime_stp)*1000) !to ms
          write(6,104)int(telapsed),int((telapsed-int(telapsed))*60.)
          write(6,108)int(tmiss),int((tmiss-int(tmiss))*60.)
          if(tpernondt.gt.60.)then
            write(6,105)tpernondt
          else
            write(6,106)tpernondt/60.0d0
          endif
          write(6,107)time-deltatime, int((time-deltatime)/param(14))+1 !,deltatime,param(10)
          write(6,*)
          write(6,*)

          !write(6,*)'ramp, pert = ',1./(1.+exp(4.-0.25*time)),(1.+0.05*cos(time*8.*atan(1.)*uparam(8)))

        endif

      endif

  102 format('      Mean time per timestep: ',F8.4,'  dev:',I8,'ms')
  103 format('      Re=',F8.2,F8.2)
  104 format('      Elapsed time: ',I8,' h ',I2,' min')
  105 format('      Time per nondimensional time: ',F8.2,' sec')
  106 format('      Time per nondimensional time: ',F8.2,' min ')
  107 format('      Local time: ',F8.4,'  File:',I8) !,' StartFrom=',F8.4,'endTime=',F8.4)
  108 format('      Remaining time: ',I8,' h ',I2,' min')

      return
      end
c-----------------------------------------------------------------------
      subroutine estimate_strouhal

      include 'SIZE'
      include 'TOTAL'

      real tlast,vlast,tcurr,vcurr,t0,t1
      save tlast,vlast,tcurr,vcurr,t0,t1
      data tlast,vlast,tcurr,vcurr,t0,t1 / 6*0 /

      integer e,eg,eg0,e0

      eg0 = 622          ! Identify element/processor in wake
      mid = gllnid(eg0)
      e0  = gllel (eg0)

      st  = 0

      if (nid.eq.mid) then

         tlast = tcurr
         vlast = vcurr

         tcurr = time
         vcurr = vy (1,ny1,1,e0)

         xcurr = xm1(1,ny1,1,e0)
         ycurr = ym1(1,ny1,1,e0)

         write(6,2) istep,time,vcurr,xcurr,ycurr
    2    format(i9,1p4e13.5,' vcurr')

         if (vlast.gt.0.and.vcurr.le.0) then ! zero crossing w/ negative slope
            t0  = t1
            t1  = tlast + (tcurr-tlast)*(vlast-0)/(vlast-vcurr)
            per = t1-t0
            if (per.gt.0) st = 1./per
         endif
      endif

      st = glmax(st,1)

      n  = nx1*ny1*nz1*nelv
      ux = glamax(vx,n)
      uy = glamax(vy,n)

      if (nid.eq.0.and.st.gt.0) write(6,1) istep,time,st,ux,uy
    1 format(i5,1p4e12.4,' Strouhal')

      return
      end
c-----------------------------------------------------------------------
