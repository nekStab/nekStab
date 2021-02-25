c-----------------------------------------------------------------------
      subroutine sfd !selective frequency damping with Euler
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt=lx1*ly1*lz1*lelt
      real, save, dimension(lt)   :: qxo,qyo,qzo
      real, save, dimension(lt)   :: vxo,vyo,vzo
      real, dimension(lt)         :: do1,do2,do3

      real residu,h1,l2,semi,linf,rate,residu0,cutoff,gain,frq,sig,umax
      save residu0

      if(istep.eq.0) then

        call opcopy (qxo,qyo,qzo, vx,vy,vz)
        call opcopy (vxo,vyo,vzo, vx,vy,vz)

        residu0=0.0D0
        if(nid.eq.0)open(unit=10,file='residu.dat')

      else

        !cutoff = uparam(04)
        !gain = -uparam(05)
   
        frq = uparam(04)*8*atan(1.0D0) ! transform St to omega
        sig = uparam(05)

        !Akervik 2006
        !cutoff = 0.5*frq
        !gain = -2*sig

        !Casacuberta 2018 JCP 375:481-497
        cutoff = 0.5*(sqrt(frq**2+sig**2)-sig)
        gain = -0.5*(sqrt(frq**2+sig**2)+sig)

        call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo)
        call opcmult(do1,do2,do3, cutoff*dt)
        call opadd2 (qxo,qyo,qzo, do1,do2,do3)

        call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo)
        call opcmult(do1,do2,do3, gain)
        call opadd2 (fcx,fcy,fcz, do1,do2,do3)

        call opsub2 (vxo,vyo,vzo,vx,vy,vz)
        call normvc (h1,semi,l2,linf,vxo,vyo,vzo)
        residu = l2!/dt; 
        rate = (residu-residu0); residu0 = residu

        call opcopy (vxo,vyo,vzo, vx,vy,vz)

          if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,residu/dt
            write(6,"(A,2E15.7)")' Rs, rate =',residu,rate
            write(6,*)
          endif

          if( istep.gt.100 .and. residu .lt. max(param(21),param(22)) )then !save to disk and change flag
            if(nid.eq.0)write(6,*)' Converged base flow to:',max(param(21),param(22))
            ifbfcv = .true.

            call compute_cfl(umax,vx,vy,vz,1.0)
            if (nid.eq.0) write(6,*) 'CFL=0.5 dt max=',0.50d0/umax
            if (nid.eq.0) write(6,*) 'CFL=5 dt max=',5.0d0/umax

            param(63) = 1 ! Enforce 64-bit output
            call bcast(param,200*wdsize)
            call outpost(vx,vy,vz,pr,t,'BF_')
            param(63) = 0 ! Enforce 32-bit output
            call bcast(param,200*wdsize)
          endif

          if(istep.eq.nsteps)close(10)
        endif
      return
      end
c-----------------------------------------------------------------------
      subroutine sfd_ab3 !SFD with higher-order temporal scheme and variable time-step
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt=lx1*ly1*lz1*lelt

      real, dimension(lt), save   :: qxo,qyo,qzo,
     $                               vxo,vyo,vzo,
     $                               uta,vta,wta,
     $                               utb,vtb,wtb,
     $                               utc,vtc,wtc

      real, dimension(lt)         :: do1,do2,do3,
     $                               utmp,vtmp,wtmp

      real adt,bdt,cdt
      real residu,h1,l2,semi,linf,rate,residu0,cutoff,gain,frq,sig,umax
      real glmin,glmax,glsum,glsc3
      integer n,i
      save n,residu0

      if(uparam(4).gt.0. OR. uparam(5).gt.0) then

        !cutoff = uparam(04)
        !gain = -uparam(05)
   
        frq = uparam(04)*8*atan(1.0D0) ! transform St to omega
        sig = uparam(05)

        !Akervik 2006
        !cutoff = 0.5*frq
        !gain = -2*sig

        !Casacuberta 2018 JCP 375:481-497
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

         else

            call setab3(adt,bdt,cdt)
            call opcopy(utc,vtc,wtc, utb,vtb,wtb)
            call opcopy(utb,vtb,wtb, uta,vta,wta)
            call opsub3(uta,vta,wta, vxo,vyo,vzo, qxo,qyo,qzo)!uta=vxo-qxo

            call opcopy (do1,do2,do3, uta,vta,wta)
            call opcmult(do1,do2,do3, adt)
            call opcopy (utmp,vtmp,wtmp, do1,do2,do3)

            call opcopy (do1,do2,do3, utb,vtb,wtb)
            call opcmult(do1,do2,do3, bdt)
            call opadd2 (utmp,vtmp,wtmp, do1,do2,do3)

            call opcopy (do1,do2,do3, utc,vtc,wtc)
            call opcmult(do1,do2,do3, cdt)
            call opadd2 (utmp,vtmp,wtmp, do1,do2,do3)

            call opcmult(utmp,vtmp,wtmp, cutoff*dt)
            call opadd2(qxo,qyo,qzo, utmp,vtmp,wtmp)

          endif

        call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo)!fc=vo-qo
        call opcmult(do1,do2,do3, gain)!f=fc*-chi
        call opadd2 (fcx,fcy,fcz, do1,do2,do3)!FORCE HERE DO NOT COPY, ADD!
        else
            if(nid.eq.0)open(unit=10,file='residu.dat')
            if(nid.eq.0)write(6,*)' SFD in continuation mode'
        endif

        if(istep.ge.1)then

          call opsub2(vxo,vyo,vzo,vx,vy,vz)
          call normvc(h1,semi,l2,linf,vxo,vyo,vzo)
          residu = l2!/dt; 
          rate = (residu-residu0); residu0 = residu
          call opcopy(vxo,vyo,vzo,vx,vy,vz)

          if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,rate
            write(6,"(A,2E15.7)")' Rs, rate =',residu,rate
            write(6,*)
          endif

          if( istep.gt.100 .and. residu .lt. max(param(21),param(22)) )then !save to disk and change flag
            if(nid.eq.0)write(6,*)' Converged base flow to:',max(param(21),param(22))
            ifbfcv = .true.

            call compute_cfl(umax,vx,vy,vz,1.0)
            if (nid.eq.0) write(6,*) 'CFL=0.5 dt max=',0.50d0/umax
            if (nid.eq.0) write(6,*) 'CFL=5 dt max=',5.0d0/umax

            param(63) = 1 ! Enforce 64-bit output
            call bcast(param,200*wdsize)
            call outpost(vx,vy,vz,pr,t,'BF_')
            param(63) = 0 ! Enforce 32-bit output
            call bcast(param,200*wdsize)
          endif

        if(istep.eq.nsteps)close(10)
        endif
      return
      end
c-----------------------------------------------------------------------
