c-----------------------------------------------------------------------
      subroutine userchk
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

      if(uparam(1).ne.3)then

         if(istep.eq.0)then

            if(nid.eq.0)open(unit=11,file='stats.dat')
            ifbfcv=.false.

            n = nx1*ny1*nz1*nelv
            re = 1.0d0/param(2) !inital value
            xmn = glmin(xm1,n); xmx = glmax(xm1,n)
            ymn = glmin(ym1,n); ymx = glmax(ym1,n)
            zmn = glmin(zm1,n); zmx = glmax(zm1,n)
            eek = 0.5/volvm1

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
      subroutine userf (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,eg
      if(uparam(1).gt.1)then    !pre composed forced on
         e = gllel(eg)
         ffx = fcx(ix,iy,iz,e)
         ffy = fcy(ix,iy,iz,e)
         ffz = fcz(ix,iy,iz,e)
      else
         ffx=0.0d0; ffy=0.0d0; ffz=0.0d0
      endif

      
      spng_str = 1.0
      if(spng_str.gt.0.0)then
        if(istep.eq.0)then 

                  spng_wl(1)=1.0 ! Sponge left section width; dimension X
                  spng_wl(2)=0.0  ! Sponge left section width; dimension Y
          if(IF3D)spng_wl(3)=0.0  ! Sponge left section width; dimension Z

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
       call spng_forcing(ffx,ffy,ffz,ix,iy,iz,eg) !add forcing to ff
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,eg

      if (JP.eq.0) then
         e  = gllel(eg)
         ux=1.0d0
         uy=0.0d0
         uz=0.0d0
         temp=0.0d0
      else                      ! perturbation
         ux = 0.0d0
         uy = 0.0d0
         uz = 0.0d0
         temp = 0.0d0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,eg

      if (JP.eq.0) then         ! base flow
         e  = gllel(eg)
         ux = 1.0d0
         uy = 0.0d0
         uz = 0.0d0
         temp=0.0d0
      else                      ! perturbation
         ux = 0.0d0
         uy = 0.0d0
         uz = 0.0d0
         temp = 0.0d0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      udiff = 0.0d0
      utrans = 0.0d0
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      qvol = 0.0d0
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

      if(nid.eq.0)then
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
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      return
      end

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
