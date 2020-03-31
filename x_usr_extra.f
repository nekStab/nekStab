c-----------------------------------------------------------------------
      subroutine nekStab
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lt=lx1*ly1*lz1*lelt

      real wo1(lt),wo2(lt),wo3(lt),vort(lt,3)
      common /ugrad/ wo1,wo2,wo3,vort

      real xmn,xmx,ymn,ymx,zmn,zmx,period
      real uek,vek,wek,eek
      real glsc3,glmin,glmax,glsum,re,viscos,v_now
      logical ifto_sav
      integer i,n,ntot,steps
      save eek,n,re,ntot

      call oprzero(fcx,fcy,fcz) ! never comment this!
      if(istep.eq.0)call print_parameters

      if(uparam(1).ne.3)then

         if(istep.eq.0)then

            if(nid.eq.0)open(unit=11,file='stats.dat')
            ifbfcv=.false.

            n = nx1*ny1*nz1*nelv
            re = 1.0d0/param(2) !inital value
            xmn = glmin(xm1,n); xmx = glmax(xm1,n)
            ymn = glmin(ym1,n); ymx = glmax(ym1,n)
            zmn = glmin(zm1,n); zmx = glmax(zm1,n)
            eek = 0.50d0/volvm1

            if(nid.eq.0)then
               write(6,*)' x min max =',xmn,xmx
               write(6,*)' y min max =',ymn,ymx
               write(6,*)' z min max =',zmn,zmx
               write(6,*)' eek =',eek
            endif

         !viscous sponge
         !do i=1,lx1*ly1*lz1*lelt!*ldimt1
         !vdiff(i,1,1,1,:)=(1./(0.6+0.4*tanh(0.8*xmx-xm1(i,1,1,1))))*param(2)
         !enddo
         !call outpost(vdiff(1,1,1,1,1),vdiff(1,1,1,1,2),vdiff(1,1,1,1,3),pr,t,'VDF')

         endif !istep.eq.0

         if((istep.eq.0).OR.ifoutfld)then

            call oprzero(wo1,wo2,wo3)
            call oprzero(vort(:,1),vort(:,2),vort(:,3))
            call comp_vort3(vort,wo1,wo2,vx,vy,vz)
            ifto = .false.; ifpo = .false.
            call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t, 'vor')
            ifto = .true.; ifpo = .true.
            call lambda2(t(1,1,1,1,1))

         endif
         if(istep.eq.0)then
         ifpo=.true.;call outpost(vx,vy,vz,pr,t,'   ')
         endif

         if (mod(istep,10).eq.0) then
            uek = glsc3(vx,bm1,vx,n)*eek
            vek = glsc3(vy,bm1,vy,n)*eek
            if(if3d)then
               wek = glsc3(vz,bm1,vz,n)*eek
               write(11,"(5E15.7)")time,uek,vek,wek,uek+vek+wek
            else
               write(11,"(4E15.7)")time,uek,vek,uek+vek
            endif
         endif

        ifto_sav=ifto
        ifpo = .false.; ifto = .false.
         call hpts
        ifpo = ifto_sav; ifto = ifto_sav

         if(uparam(01).eq.1)then !compose forcings to fcx,fcy,fcz
            if(uparam(03).eq.1)call sfd
            if(uparam(03).eq.1.1)call sfd_ab3
         endif

         call mycomment

         if( ifbfcv )then

          call switch_to_lnse_steady
          call krylov_schur
          call exitt

         endif !ifbfcv

      else !uparam(01)==3

         !call switch_to_lnse_steady
         call krylov_schur
         call exitt

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sponge !here we compose fcx,fcy,fcz
      implicit none
      include 'SIZE'
      include 'TOTAL'
      !integer, parameter :: lt=lx1*ly1*lz1*lelt

      !if(uparam(1).gt.1)then    !pre composed forced on
      !   e = gllel(eg)
      !   ffx = fcx(ix,iy,iz,e)
      !   ffy = fcy(ix,iy,iz,e)
      !   ffz = fcz(ix,iy,iz,e)
      !else
      !   ffx=0.0d0; ffy=0.0d0; ffz=0.0d0
      !endif


      spng_str = 0.0
      if(spng_str.gt.0.0)then
        if(istep.eq.0)then

                  spng_wl(1)=1.0 ! Sponge left section width; dimension X
                  spng_wl(2)=0.0 ! Sponge left section width; dimension Y
          if(IF3D)spng_wl(3)=0.0 ! Sponge left section width; dimension Z

                 spng_wr(1)=40.0 ! Sponge right section width; dimension X
                 spng_wr(2)= 0.0
          if(IF3D)spng_wr(3)=0.0

                 spng_dl(1)=0.333*spng_wl(1) ! Sponge left drop/rise section width; dimension X
                 spng_dl(2)=0.333*spng_wl(2)
         if(IF3D)spng_dl(3)=0.333*spng_wl(3)

                spng_dr(1)=0.333*spng_wr(1) ! Sponge right drop/rise section width; dimension X
                spng_dr(2)=0.333*spng_wr(2)
        if(IF3D)spng_dr(3)=0.333*spng_wr(3)

        call spng_init

       endif

!       DO 100 IEL=1,NELV
!          ielg = lglel(iel)
!          DO 100 K=1,lz1
!          DO 100 J=1,ly1
!          DO 100 I=1,lx1
             !CALL USERF   (I,J,K,IELG)
             !fcx(I,J,K,IEL) = FFX
             !fcy(I,J,K,IEL) = FFY
             !fcz(I,J,K,IEL) = FFZ
!  100  CONTINUE

       !call spng_forcing(ffx,ffy,ffz,ix,iy,iz,eg) !add forcing to ff
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mycomment
      include 'SIZE'
      include 'TOTAL'
      real*8,save :: eetime0,eetime1,eetime2
      data           eetime0,eetime1,eetime2 /0.0d0, 0.0d0, 0.0d0/
      real, save :: deltatime
      real telapsed,tpernondt,tmiss

      if (courno.gt.4.0d0) then
        if (nio.eq.0)then
          write(6,*)
          write(6,*)'    cfl > 4. stopping'
          write(6,*)
         endif
        call exitt
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
      subroutine print_parameters
      include 'SIZE'
      include 'TOTAL'
      if(nid.eq.0)then
         write(6,*)'P10=',param(10),''
         write(6,*)'P11=',param(11),''
         write(6,*)'P12=',param(12),''
         write(6,*)'P13=',param(13),''
         write(6,*)'P21=',param(21),'pressure sol tol'
         write(6,*)'P22=',param(22),'velocity sol tol'
         write(6,*)'P24=',param(24),'rel tol Helmholtz sol'
         write(6,*)'P25=',param(25),'abs tol Helmholtz sol'
         write(6,*)'P27=',param(27),'temporal integration order'
         write(6,*)'P31=',param(31),'numberOfPerturbations'
         write(6,*)'P41=',param(41),'1 for multiplicative SEMG'
         write(6,*)'P42=',param(42),'lin solv for the pres equation 0:GMRES,1:CG'
         write(6,*)'P43=',param(43),'0:additive multilevel scheme 1:orig 2lvl sch'
         write(6,*)'P44=',param(44),'0=E-based addit Schwarz PnPn-2;1=A-based'
         write(6,*)'P93=',param(93),'n of prev sol to use for res proj'
         write(6,*)'P94=',param(94),'n steps star res proj for vel and pas.scal'
         write(6,*)'P95=',param(95),'projection for pressure solve on/off'
         write(6,*)'uparam1=',uparam(1)
         write(6,*)'uparam01=',uparam(01)
         write(6,*)'uparam02=',uparam(02)
         write(6,*)'uparam03=',uparam(03)
         write(6,*)'uparam04=',uparam(04)
         write(6,*)'uparam05=',uparam(05)
         write(6,*)'uparam06=',uparam(06)
         write(6,*)'uparam07=',uparam(07)
         write(6,*)'uparam08=',uparam(08)
         write(6,*)'uparam09=',uparam(09)
         write(6,*)'uparam10=',uparam(10)
      endif
      end
c-----------------------------------------------------------------------
