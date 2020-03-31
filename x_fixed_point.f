c-----------------------------------------------------------------------
      subroutine sfd !original selective frequency damping
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt=lx1*ly1*lz1*lelt

      !real, save, dimension(lt)   :: qxo,qyo,qzo
      !real, save, dimension(lt)   :: vxo,vyo,vzo

      real, dimension(lt)   :: qxo,qyo,qzo
      common /teste1/ qxo,qyo,qzo

      real, dimension(lt)   :: vxo,vyo,vzo
      common /teste2/ vxo,vyo,vzo

      real, dimension(lt)         :: do1,do2,do3

      real residu,h1,l2,semi,linf,rate,residu0,cutoff,gain
      save residu0

      if(istep.eq.0) then

        call opcopy (qxo,qyo,qzo, vx,vy,vz)
        call opcopy (vxo,vyo,vzo, vx,vy,vz)

        residu0=0.0D0
        if(nid.eq.0)open(unit=10,file='residu.dat')

      else

        cutoff = uparam(04)*4*atan(1.0D0)*dt
        gain = -uparam(05)

        !cutoff = uparam(04)*4*atan(1.0D0)*dt
        cutoff = 0.2*dt
        gain = -0.125 !uparam(05)

        call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo)
        call opcmult(do1,do2,do3, cutoff)
        call opadd2 (qxo,qyo,qzo, do1,do2,do3)

        call opsub3 (do1,do2,do3, vxo,vyo,vzo, qxo,qyo,qzo)
        call opcmult(do1,do2,do3, gain)
        call opadd2 (fcx,fcy,fcz, do1,do2,do3)

        call opsub2 (vxo,vyo,vzo,vx,vy,vz)
        call normvc (h1,semi,l2,linf,vxo,vyo,vzo)
        residu = l2/dt; rate = (residu-residu0); residu0 = residu

        call opcopy(vxo,vyo,vzo, vx,vy,vz)

          if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,rate
            write(6,"(A,3E15.7)")' Rs, rate =',residu,rate,l2
            write(6,*)
          endif

          if( residu .lt. 1.0e-09 )then !save to disk and change flag
            if(nid.eq.0)write(6,*)' Converged base flow to 1.0E-09...'
            ifbfcv = .true.
            call outpost(vx,vy,vz,pr,t,'BF_')
          endif

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sfd_ab3
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt=lx1*ly1*lz1*lelt

      !real, dimension(lt), save   :: qxo,qyo,qzo
      !real, dimension(lt), save   :: vxo,vyo,vzo
      real, dimension(lt)         :: do1,do2,do3

      real, dimension(lt)   :: qxo,qyo,qzo
      common /teste1/ qxo,qyo,qzo

      real, dimension(lt)   :: vxo,vyo,vzo
      common /teste2/ vxo,vyo,vzo

      real utmp(lt),vtmp(lt),wtmp(lt),
     $     uta(lt), vta(lt), wta(lt),
     $     utb(lt), vtb(lt), wtb(lt),
     $     utc(lt), vtc(lt), wtc(lt)

      common /tempo/ utmp,vtmp,wtmp,
     $               uta,vta,wta,
     $               utb,vtb,wtb,
     $               utc,vtc,wtc

      real adt,bdt,cdt
      real residu,h1,l2,semi,linf,residu0,rate
      real glmin,glmax,glsum,glsc3,cutoff,gain
      integer n,i
      save n,residu0

      !cutoff = uparam(04)*4*atan(1.0D0)*dt
      cutoff = 0.2*dt
      gain = -0.125 !uparam(05)

          if(istep.eq.0)then
            n = nx1*ny1*nz1*nelt
            call oprzero(fcx,fcy,fcz)
            call oprzero(utmp,vtmp,wtmp)
            call oprzero(uta,vta,wta)
            call oprzero(utb,vtb,wtb)
            call oprzero(utc,vtc,wtc)
            residu0=0.0D0
            if(nid.eq.0)open(unit=10,file='residu.dat')
            call opcopy(qxo,qyo,qzo,vx,vy,vz)
            call opcopy(vxo,vyo,vzo,vx,vy,vz)
          endif

          if(istep.eq.1)then !start with Euler

            call opsub3 (uta,vta,wta, vxo,vyo,vzo, qxo,qyo,qzo)
            call opcopy (utmp,vtmp,wtmp, uta,vta,wta)
            call opcmult(uta,vta,wta, cutoff)
            call opadd2 (qxo,qyo,qzo, utmp,vtmp,wtmp)

          elseif(istep.eq.2 .AND. param(27).gt.1)then !then step with AB2 fixed step

            call opcopy(utb,vtb,wtb, uta,vta,wta)
            call opsub3(uta,vta,wta, vxo,vyo,vzo, qxo,qyo,qzo)
            utmp = cutoff*(1.5*uta-0.5*utb)
            vtmp = cutoff*(1.5*vta-0.5*vtb)
            wtmp = cutoff*(1.5*wta-0.5*wtb)
            call opadd2(qxo,qyo,qzo,utmp,vtmp,wtmp)

          elseif(istep.ge.3 .AND. param(27).gt.2 )then

            call setab3(adt,bdt,cdt)
            call opcopy(utc,vtc,wtc,utb,vtb,wtb)
            call opcopy(utb,vtb,wtb,uta,vta,wta)
            call opsub3 (uta,vta,wta,vxo,vyo,vzo,qxo,qyo,qzo)!uta=vxo-qxo

            !call opcopy (do1,do2,do3, uta,vta,wta)
            !call opcmult(do1,do2,do3, adt)
            !call opcopy (utmp,vtmp,wtmp, do1,do2,do3)
            !call opcopy (do1,do2,do3, utb,vtb,wtb)
            !call opcmult(do1,do2,do3, bdt)
            !call opadd2 (utmp,vtmp,wtmp, do1,do2,do3)
            !call opcopy (do1,do2,do3, utc,vtc,wtc)
            !!call opcmult(do1,do2,do3, cdt)
            !call opadd2 (utmp,vtmp,wtmp, do1,do2,do3)

            !call opcmult(utmp,vtmp,wtmp, cutoff)

            utmp = cutoff*(adt*uta+bdt*utb+cdt*utc)
            vtmp = cutoff*(adt*vta+bdt*vtb+cdt*vtc)
            wtmp = cutoff*(adt*wta+bdt*wtb+cdt*wtc)

            call opadd2(qxo,qyo,qzo,utmp,vtmp,wtmp)

          endif

        call opsub3 (do1,do2,do3,vxo,vyo,vzo,qxo,qyo,qzo)!fc=vo-qo
        call opcmult(do1,do2,do3,gain)!f=-chi*uo
        call opadd2 (fcx,fcy,fcz,do1,do2,do3)!f=uo !FORCE HERE DONT COPY! ADD!

        if(istep.ge.1)then

          call opsub2(vxo,vyo,vzo,vx,vy,vz)
          call normvc(h1,semi,l2,linf,vxo,vyo,vzo)
          residu = l2/dt; rate = (residu-residu0); residu0 = residu
          call opcopy(vxo,vyo,vzo,vx,vy,vz)

          if(nid.eq.0)then
            write(10,"(3E15.7)")time,residu,rate
            write(6,"(A,3E15.7)")' Rs, rate =',residu,rate,l2
            write(6,*)
          endif

          if( residu .lt. 1.0e-09 )then !save to disk and change flag
            if(nid.eq.0)write(6,*)' Converged base flow to 1.0E-09...'
            ifbfcv = .true.
            call outpost(vx,vy,vz,pr,t,'BF_')
          endif

        elseif(istep.eq.nsteps)then
          close(10)
        endif

      return
      end
c-----------------------------------------------------------------------
