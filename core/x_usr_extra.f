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
      logical ifto_sav, ifpo_sav
      integer i,n,ntot,steps
      save eek,n,re,ntot

      call oprzero(fcx,fcy,fcz) ! never comment this!
      if(istep.eq.0) then
         call print_parameters
         if(uparam(10).gt.0)then
            if(nid.eq.0)write(6,*)' Initializing sponge...'
            spng_str = 1.0d0
            call spng_init
         endif
      endif

      if(uparam(1).le.2)then

         if(istep.eq.0)then

            !if(nid.eq.0)open(unit=11,file='stats.dat')
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

         endif !istep.eq.0

         if((istep.eq.0).OR.ifoutfld)then

            call oprzero(wo1,wo2,wo3)
            call oprzero(vort(:,1),vort(:,2),vort(:,3))
            call comp_vort3(vort,wo1,wo2,vx,vy,vz)
            ifto_sav = ifto; ifpo_sav = ifpo; ifto = .false.; ifpo = .false.
            call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t, 'vor')
            ifto = ifto_sav ; ifpo = ifpo_sav

         endif
         if(istep.eq.0)call outpost(vx,vy,vz,pr,t,'   ')


        ! if (mod(istep,10).eq.0) then
        !    uek = glsc3(vx,bm1,vx,n)*eek
        !    vek = glsc3(vy,bm1,vy,n)*eek
        !    if(if3d)then
        !       wek = glsc3(vz,bm1,vz,n)*eek
        !       write(11,"(5E15.7)")time,uek,vek,wek,uek+vek+wek
        !    else
        !       write(11,"(4E15.7)")time,uek,vek,uek+vek
        !    endif
        ! endif
        !if(istep.eq.lastep)close(11)

        call hpts

         if(uparam(01).eq.1)then !compose forcings to fcx,fcy,fcz
            !if(uparam(03).eq.1)call sfd
            if(uparam(03).eq.1)call sfd_ab3
         endif

         call mycomment

         if( ifbfcv )then ! after converging base fow...

         !call switch_to_lnse_steady !in utilities.f
         !call krylov_schur ! in eigensolvers.f
         if(nid.eq.0)write(6,*)'Stopping code...'
         call nek_end

         endif !ifbfcv

      elseif(uparam(01).ge.3)then !3:direct,4:adj,5:dir-adj,6:adj-dir

         param(12) = -abs(param(12)) !freeze dt
         param(31) = 1 ; npert = param(31)
         call bcast(param,200*wdsize) !broadcast all parameters to processors
         call krylov_schur ! in eigensolvers.f
         if(nid.eq.0)write(6,*)'Stopping code...'
         call nek_end

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
      subroutine estimate_strouhal !original routine from NekExamples

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
      subroutine set_rjet(ub) !round jet profile for axissymetric jet
      include 'SIZE'
      include 'TOTAL'
      real ub(1),theta_0
      theta_0=0.0250d0
      do i=1,nx1*ny1*nz1*nelv
         x = xm1(i,1,1,1)
         y = ym1(i,1,1,1)
         ub(i)=0.50d0*(1.0d0-tanh((1.0d0/(4.0d0*theta_0))*(y-(1.0d0/(4.0d0*y)))))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_sb (v_jet)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      real v_jet(1),A_1,x_m,H,f_x
      integer n

      n = nx1*ny1*nz1*nelv
      A_1 = uparam(06)
      x_m = uparam(07)

      do i=1,n
      x = xm1(i,1,1,1)
      y = ym1(i,1,1,1)
      if(y.eq.0.)then
         H = exp( -((x-x_m)**2)/(3.10d0**2))
         f_x = 15.18750d0*H**5 -35.43750d0*H**4 +20.250d0*H**3
         v_jet(i)=A_1*f_x
      else
         v_jet(i)=0.0d0
      endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine outpost_blayer_pert
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter   :: lt=lx1*ly1*lz1*lelv
      real, dimension(lt)  :: do1,do2,do3
      real ampx, ampy, glamax
      integer n
      logical ifto_sav, ifpo_sav
      
      n = nx1*ny1*nz1*nelv

      if((istep.eq.0).OR.ifoutfld)then

         call opsub3 ( do1,do2,do3, vx,vy,vz, ubase,vbase,wbase)

         ifto_sav = ifto; ifpo_sav = ifpo; ifto = .false.; ifpo = .false.
         call outpost( do1,do2,do3,pr,t,'per')
         ifto = ifto_sav ; ifpo = ifpo_sav

         ampx = glamax(do1,n)      
         ampy = glamax(do2,n)

         if(nid.eq.0)then
            if(istep.eq.0)then
            open(unit=111,file='ts_amp.dat',status='unknown',form='formatted')
            write(112,'(A)')'#  t  A  up  vp  up2  vp2'
            endif
            write(111,"(6E15.7)")time,uparam(06),ampx,ampy,ampx**2,ampy**2
            if(istep.eq.nsteps)close(111)
         endif
      endif      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_blasius(ub,vb) ! Compute the Blasius profile
         include 'SIZE'
         include 'TOTAL'

         real ub(1),vb(1),x_m
         real,save :: xmx,xmn
         n = nx1*ny1*nz1*nelv

         xmx = glmax(xm1,n)
         xmn = glmin(xm1,n)
         x_m = uparam(07)

         visc = param(2)/param(1) !density / dynamic viscosity
         delta99_0 = 5.0d0/1.72080d0 !2.9
         delta_star= 1.0d0
         u_0   = 1.0d0
         x_0 = (delta_star/1.7208d0)**2 / visc * u_0  ! Reference x

         x_inflow = (605.0d0/740.0d0)**2 * x_0 !original blasius
         x_inflow = x_0

         if(nid.eq.0)then
         write(6,*)'visc=',visc
         write(6,*)'delta99_0=',delta99_0
         write(6,*)'delta_star=',delta_star
         write(6,*)'u_0=',u_0
         write(6,*)'x_0=',x_0 
         write(6,*)'x_inflow=',x_inflow
         write(6,*)'x_m=',x_m
         endif

         do i=1,n
            x = xm1(i,1,1,1)
            y = ym1(i,1,1,1)

            x_t = x_inflow + x
            rex = u_0 * x_t / visc

            if(x.eq.xmn)write(6,*)'if x,rex=',real(x,4),real(rex,4),real(sqrt(rex),4)
            if(x.eq.x_m)write(6,*)'sb x,rex=',real(x,4),real(rex,4),real(sqrt(rex),4)
            if(x.eq.xmx)write(6,*)'of x,rex=',real(x,4),real(rex,4),real(sqrt(rex),4)
                     
            eta = y*sqrt(rex)/x_t
            call blasius(ub(i),vb(i),eta,rex)

         enddo

         return
         end
c-----------------------------------------------------------------------
      subroutine blasius(u,v,eta,rex)

         integer icalld
         save    icalld
         data    icalld /0/

         parameter (lb=100)
         real blasius_soln(0:4,0:lb)
         save blasius_soln
c
c     Algorithm found in Lighthills monograph on mathematical fluid mech.
c     (c. of M. Choudhari)
c
         real  w(4)

         twok2 =  1.6551903602308323382003140460740
         rk2   =  0.5*twok2
         rk    =  sqrt(rk2)

         if (icalld.eq.0) then
            icalld = 1

            call set_ics (blasius_soln(0,0))

            dt=.05
            rr=1.0725
            do i=1,lb
               blasius_soln(0,i) = blasius_soln(0,i-1) + dt
               dt = dt*rr
            enddo

            do i=1,lb
               eta0 = blasius_soln(0,i-1)
               eta1 = blasius_soln(0,i  )
               t0 = 0.5*eta0/rk
               t1 = 0.5*eta1/rk
               dt = .0005    !  Note, this is good to about 12 digits
               m=3
               call copy(blasius_soln(1,i),blasius_soln(1,i-1),m)
               call rk4_integrate(blasius_soln(1,i),3,t1,t0,dt)
            enddo
         endif

         if (eta.gt.blasius_soln(0,lb)) then

            call copy(w,blasius_soln(1,lb),2)

         else

            i = interval_find(eta,blasius_soln,5,lb)

            eta0 = blasius_soln(0,i)
            t0   = 0.5*eta0/rk
            t1   = 0.5*eta/rk
            dt   = .0005    !  Note, this is good to about 12 digits
            m    = 3
            call copy(w,blasius_soln(1,i),m)
            call rk4_integrate(w,3,t1,t0,dt)

         endif

         g  = w(1)
         gp = w(2)

         f  = g  / rk
         fp = gp / twok2

         u  = fp
         v  = 0.5*(eta*fp-f)/sqrt(rex)

c     write(6,1) eta,u,v,f,fp,rex
c  1  format(1p6e14.6,' eta')

         return
         end
c-----------------------------------------------------------------------
      subroutine rk4_integrate(w,n,tfinal,tstart,dti) !Program to integrate dW/dt = F(W,t)
!     Input:   w() is initial condition at t=tstart
!     Output:  w() is solution at t = tfinal
!     n = length of vector
      real  w(n)
      if (tfinal.gt.tstart .and. dti.gt.0.) then

         tdelta = tfinal-tstart
         dt     = dti
         nsteps = tdelta/dt
         nsteps = max(nsteps,1)
         dt     = tdelta/nsteps

         t = tstart
         do k=1,nsteps !  TIME STEPPING

            call rk4 (w,t,dt,n) ! Single RK4 step (nmax=4)
            t = t+dt

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_ics (w)!Initial conditions for modified Blasius equation g''' + g g'' = 0
      real  w(0:3)
      w(0) = 0.0d0    ! eta = 0
      w(1) = 0.0d0    ! g
      w(2) = 0.0d0    ! g'
      w(3) = 1.0d0    ! g"
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_f(f,w,t) !Compute RHS of ODE:
      real  f(4),w(4)
      f(1) = w(2)
      f(2) = w(3)
      f(3) = -w(1)*w(3)
      return
      end
c-----------------------------------------------------------------------
      subroutine add3s2nd(x,y,z,c,n)
      real  x(1),y(1),z(1),c
      do i=1,n
         x(i) = y(i) + c*z(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rk4(w,t,dt,n)
      real  w(1),t,dt
      real  wh(4),f1(4),f2(4),f3(4),f4(4)

      dt2 = dt/2.0d0
      dt3 = dt/3.0d0
      dt6 = dt/6.0d0

      t2 = t+dt2
      tt = t+dt

      call compute_f (f1,w ,t )
      call add3s2nd  (wh,w,f1,dt2,n)

      call compute_f (f2,wh,t2)
      call add3s2nd  (wh,w,f2,dt2,n)

      call compute_f (f3,wh,t2)
      call add3s2nd  (wh,w,f3,dt ,n)

      call compute_f (f4,wh,tt)

      do i=1,n
         w(i) = w(i) + dt6*(f1(i)+f4(i)) + dt3*(f2(i)+f3(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      function interval_find(x,xa,m,n) !Find interval. p. 88-89, numerical recipes
      real xa(m,0:n)

      if (x.ge.xa(1,n)) then
         interval_find = n
      elseif (x.le.xa(1,0)) then
         interval_find = 0
      else

         klo=0
         khi=n
   1    if ((khi-klo).gt.1) then
         k=(khi+klo)/2
         if (xa(1,k).gt.x) then
            khi=k
         else
            klo=k
         endif
         goto 1
         endif

         h=xa(1,khi)-xa(1,klo)
         if (h.eq.0) then
            write(6,*) xa(1,klo),xa(1,khi),klo,khi,'ERROR: Zero jump in interval_find.'
            return
         endif
         interval_find = klo
      endif

      return
      end
c-----------------------------------------------------------------------
