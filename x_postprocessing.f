!!!!! sponge region and noise generation routines
!!!!! adapted from KTH Toolbox https://github.com/KTH-Nek5000/KTH_Toolbox
! spng_init ->
! spng_set ->
! mth_stepf ->
! spng_forcing ->
!
! add_noise
! mth_rand
c-----------------------------------------------------------------------
      subroutine switch_to_lnse_steady(steps_in)
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'OPCTR'
      include 'CTIMER'

      integer, intent(in), optional :: steps_in
      integer steps
      real sampling_period

      if(.not.present(steps_in))then
        steps = steps_in
      else
        if(nid.eq.0)write(6,*)' nsteps not specified: computing with uparam(04)'
        sampling_period = real(1.0d0/uparam(04))
        steps = int(sampling_period/dt)
      endif

      if(nid.eq.0)then
       write(6,*)' Switching to linearized solver...'
       write(6,*)' Base flow stability: nsteps=',steps
       write(6,*)
      endif

      time = 0
      istep = 0
      nsteps = steps
      param(10) = 0     !endTime
      param(11) = steps !numSteps to overwrite final time
      param(14) = 999. !write interval runtime
      param(12) = -1*abs(param(12)) !freeze dt
      param(31) = 1 !numberOfPerturbations
      param(63) = 1 ! Enforce 64-bit output
      npert = param(31)

      call bcast(param,200*wdsize) !broadcast all parameters to processors

      ifpert = .true.
      call bcast(ifpert, lsize) !start linearized

      ifbase = .false.
      call bcast(ifbase, lsize) !stop dns

      call chkParam !sanity check

      return
      end
c-----------------------------------------------------------------------
      subroutine spng_init
      implicit none
      include 'SIZE'
      integer nn
      real vtmp(lx1*ly1*lz1*lelt,ldim)
      common /CTMP1/ vtmp

      nn = lx1*ly1*lz1*lelt*ldim
      call rzero(vtmp,nn)
      call spng_set(vtmp(1,1),vtmp(1,2),vtmp(1,ndim))

      return
      end
c-----------------------------------------------------------------------
      subroutine spng_set(lvx,lvy,lvz) !set sponge function and refernece fields
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      real lvx   (LX1*LY1*LZ1*LELV),lvy(LX1*LY1*LZ1*LELV),lvz(LX1*LY1*LZ1*LELV)
      real lcoord(LX1*LY1*LZ1*LELV)
      common /SCRUZ/ lcoord

      integer ierr, nhour, nmin, itmp, ntot, il, jl
      real rtmp, ltim, bmin(LDIM), bmax(LDIM)
      real xxmax, xxmax_c, xxmin, xxmin_c, arg
      real glmin, glmax, mth_stepf
      logical ltmp

      ntot = NX1*NY1*NZ1*NELV
      bmin(1) = glmin(XM1,ntot)
      bmax(1) = glmax(XM1,ntot)
      bmin(2) = glmin(YM1,ntot)
      bmax(2) = glmax(YM1,ntot)
      if(IF3D) then
        bmin(NDIM) = glmin(ZM1,ntot)
        bmax(NDIM) = glmax(ZM1,ntot)
      endif

      !spng_wl(1)=WIDTHLX
      !spng_wl(2)=WIDTHLY
      !if(IF3D)spng_wl(3)=WIDTHLZ

      !spng_wr(1)=WIDTHRX
      !spng_wr(2)=WIDTHRY
      !if(IF3D)spng_wr(3)=WIDTHRZ

      !spng_dl(1)=DROPLX
      !spng_dl(2)=DROPLY
      !if(IF3D)spng_dl(3)=DROPLZ

      !spng_dr(1)=DROPRX
      !spng_dr(2)=DROPRY
      !if(IF3D)spng_dr(3)=DROPRZ

      call rzero(spng_fun,ntot)

         ! save reference field
         call copy(spng_vr(1,1),lvx,ntot)
         call copy(spng_vr(1,2),lvy,ntot)
         if (IF3D) call copy(spng_vr(1,NDIM),lvz, ntot)

         ! for every dimension
         do il=1,NDIM

          if (spng_wl(il).gt.0.0.or.spng_wr(il).gt.0.0) then

             if (spng_wl(il).lt.spng_dl(il).or.spng_wr(il).lt.spng_dr(il)) then
                write(6,*)'Wrong sponge parameters!'
             endif

             xxmax   = bmax(il) - spng_wr(il)! sponge beginning (rise at xmax; right)
             xxmin   = bmin(il) + spng_wl(il)! end (drop at xmin; left)
             xxmax_c = xxmax    + spng_dr(il)! beginnign of constant part (right)
             xxmin_c = xxmin    - spng_dl(il)! beginnign of constant part (left)

             ! get SPNG_FUN
             if (xxmax.le.xxmin) then
                write(6,*)'Sponge too wide'
             else
                ! this should be done by pointers, but for now I avoid it
                if (il.eq.1) then
                   call copy(lcoord,XM1, ntot)
                elseif (il.eq.2) then
                   call copy(lcoord,YM1, ntot)
                elseif (il.eq.3) then
                   call copy(lcoord,ZM1, ntot)
                endif

                do jl=1,ntot
                   rtmp = lcoord(jl)
                   if(rtmp.le.xxmin_c) then ! constant; xmin
                      rtmp=spng_str
                   elseif(rtmp.lt.xxmin) then ! fall; xmin
                      arg = (xxmin-rtmp)/spng_wl(il)
                      rtmp = mth_stepf(arg)
                   elseif (rtmp.le.xxmax) then ! zero
                      rtmp = 0.0
                   elseif (rtmp.lt.xxmax_c) then ! rise
                      arg = (rtmp-xxmax)/spng_wr(il)
                      rtmp = mth_stepf(arg)
                   else    ! constant
                      rtmp = spng_str
                   endif
                   spng_fun(jl)=max(spng_fun(jl),rtmp)
                enddo

             endif  ! xxmax.le.xxmin

          endif  ! spng_w(il).gt.0.0
       enddo

      ltmp = ifto
      ifto = .TRUE.
      call outpost2(spng_vr,spng_vr(1,2),spng_vr(1,NDIM),spng_fun,spng_fun,1,'SPG')
      ifto = ltmp

      return
      end
c-----------------------------------------------------------------------
      real function mth_stepf(x) !compute sponge function
      implicit none
      real x, xdmin, xdmax
      parameter (xdmin = 0.001, xdmax = 0.999)
      if (x.le.xdmin) then
         mth_stepf = 0.0d0
      else if (x.le.xdmax) then
         mth_stepf = 1.0d0/( 1.0d0 + exp(1.0d0/(x - 1.0d0) + 1.0d0/x) )
      else
         mth_stepf = 1.0d0
      end if
      end
c-----------------------------------------------------------------------
      subroutine spng_forcing(ffx,ffy,ffz,ix,iy,iz,ieg) !compute sponge
      implicit none
      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg,iel,ip

      iel=GLLEL(ieg)
      if (SPNG_STR.gt.0.0) then
         ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(iel-1)))
         if (JP.eq.0) then
            ! dns
            ffx = ffx + SPNG_FUN(ip)*(SPNG_VR(ip,1) - VX(ix,iy,iz,iel))
            ffy = ffy + SPNG_FUN(ip)*(SPNG_VR(ip,2) - VY(ix,iy,iz,iel))
            if (IF3D) ffz = ffz + SPNG_FUN(ip)*(SPNG_VR(ip,NDIM) - VZ(ix,iy,iz,iel))
         else
            ! perturbation
            ffx = ffx - SPNG_FUN(ip)*VXP(ip,JP)
            ffy = ffy - SPNG_FUN(ip)*VYP(ip,JP)
            if(IF3D) ffz = ffz - SPNG_FUN(ip)*VZP(ip,JP)
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine add_noise(qx, qy, qz) !input random number to fields
      implicit none
      include 'SIZE'            ! NX1, NY1, NZ1, NELV, NID
      include 'TSTEP'           ! TIME, DT
      include 'PARALLEL'        ! LGLEL
      include 'INPUT'           ! IF3D
      include 'SOLN'            ! VX, VY, VZ, VMULT
      include 'GEOM'            ! XM1, YM1, ZM1

      real, dimension(lx1,ly1,lz1,lelv) :: qx, qy, qz
      integer iel,ieg,il,jl,kl,nl
      real xl(LDIM),mth_rand,fcoeff(3) !< coefficients for random distribution

      do iel=1,NELV
       do kl=1,NZ1
        do jl=1,NY1
         do il=1,NX1

         ieg = LGLEL(iel)
         xl(1) = XM1(il,jl,kl,iel)
         xl(2) = YM1(il,jl,kl,iel)
         if (if3D) xl(NDIM) = ZM1(il,jl,kl,iel)

         fcoeff(1)=  3.0e4;fcoeff(2)= -1.5e3;fcoeff(3)=  0.5e5
         qx(il,jl,kl,iel)=qx(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fcoeff)

         fcoeff(1)=  2.3e4;fcoeff(2)=  2.3e3;fcoeff(3)= -2.0e5
         qy(il,jl,kl,iel)=qy(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fcoeff)

         if (IF3D) then
           fcoeff(1)= 2.e4;fcoeff(2)= 1.e3;fcoeff(3)= 1.e5
           qz(il,jl,kl,iel)=qz(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fcoeff)
         endif

          enddo
         enddo
        enddo
       enddo

      ! face averaging
      call opdssum(qx, qy, qz)
      call opcolv (qx, qy, qz, VMULT)

      call dsavg(qx)
      call dsavg(qy)
      if(if3D) call dsavg(qz)

      call bcdirvc(qx, qy, qz,v1mask,v2mask,v3mask)

      return
      end
c-----------------------------------------------------------------------
      real function mth_rand(ix,iy,iz,ieg,xl,fcoeff) !generate random number
      implicit none
      include 'SIZE'
      include 'INPUT'       ! IF3D
      integer ix,iy,iz,ieg
      real xl(LDIM), fcoeff(3)
      mth_rand = fcoeff(1)*(ieg+xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy+fcoeff(3)*ix
      if (IF3D) mth_rand = fcoeff(1)*(ieg +xl(NDIM)*sin(mth_rand))+fcoeff(2)*iz*ix+fcoeff(3)*iz
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = cos(mth_rand)
      return
      end
c-----------------------------------------------------------------------
      subroutine opadd3 (a1,a2,a3,b1,b2,b3,c1,c2,c3)
      implicit none
      include 'SIZE'
      integer ntot1
      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),c1(1),c2(1),c3(1)
      ntot1=nx1*ny1*nz1*nelv
      call add3(a1,b1,c1,ntot1)
      call add3(a2,b2,c2,ntot1)
      if (ndim.eq.3) call add3(a3,b3,c3,ntot1)
      return
      end
c-----------------------------------------------------------------------
