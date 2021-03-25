!!!!! sponge region and noise generation routines
!!!!! adapted from KTH Toolbox https://github.com/KTH-Nek5000/KTH_Toolbox
! spng_init ->
! spng_set ->
! mth_stepf ->
! spng_forcing ->
!
! add_noise
! mth_rand

      subroutine quicksort2(n, arr, idx)

      integer n, m, nstack
      real arr(n)
      integer idx(n)
      parameter (m=7, nstack=50)
      integer i, ir, j, jstack, k, l, istack(nstack)
      real a, temp
      integer b, temp_int

      jstack = 0
      l = 1
      ir = n

 1    if (ir-l .lt. m) then
         do 12 j = l+1, ir
            a = arr(j)
            b = idx(j)
            do 11 i = j-1, l, -1
               if (arr(i) .le. a) goto 2
               arr(i+1) = arr(i)
               idx(i+1) = idx(i)
 11         continue
            i = l-1
 2          arr(i+1) = a
            idx(i+1) = b
 12      continue

         if (jstack .eq. 0) return

         ir = istack(jstack)
         l = istack(jstack - 1)
         jstack = jstack - 2
      else
         k = (l + ir) / 2
         temp = arr(k)
         arr(k) = arr(l+1)
         arr(l+1) = temp
         temp_int = idx(k)
         idx(k) = idx(l+1)
         idx(l+1) = temp_int

         if (arr(l) .gt. arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp

            temp_int = idx(l)
            idx(l) = idx(ir)
            idx(ir) = temp_int
         endif

         if (arr(l+1) .gt. arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp

            temp_int = idx(l+1)
            idx(l+1) = idx(ir)
            idx(ir) = temp_int
         endif

         if (arr(l) .gt. arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
            temp_int = idx(l)
            idx(l) = idx(l+1)
            idx(l+1) = temp_int
         endif

         i = l + 1
         j = ir
         a = arr(l+1)
         b = idx(l+1)

 3       continue

         i = i+1

         if (arr(i) .lt. a) goto 3

 4       continue

         j = j -1

         if (arr(j) .gt. a) goto 4
         if (j .lt. i) goto 5

         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp

         temp_int = idx(i)
         idx(i) = idx(j)
         idx(j) = temp_int
         goto 3

 5       arr(l+1) = arr(j)
         arr(j) = a
         idx(l+1) = idx(j)
         idx(j) = b
         jstack = jstack + 2

         if (jstack .gt. nstack) pause "..." !'NSTACK too small in quicksort2'

         if (ir-i+1 .ge. j-1) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j-1
         else
            istack(jstack) = j-1
            istack(jstack-1) = l
            l = i
         endif
      endif
      goto 1
      end
c-----------------------------------------------------------------------
      subroutine force_INCOMPLINNS
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'OPCTR'
      include 'CTIMER'
      time = 0.0d0
      istep = 0
      param(12) = -1*abs(param(12)) !freeze dt
      param(31) = 1 !numberOfPerturbations
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
      subroutine force_INCOMPLINADJNS
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'OPCTR'
      include 'CTIMER'
      time = 0.0d0
      istep = 0
      param(12) = -1*abs(param(12)) !freeze dt
      param(31) = 1 !numberOfPerturbations
      npert = param(31)
      call bcast(param,200*wdsize) !broadcast all parameters to processors
      ifpert = .true.
      call bcast(ifpert, lsize) !start linearized
      ifadj = .true.
      call bcast(ifadj, lsize) !start linearized
      ifbase = .false.
      call bcast(ifbase, lsize) !stop dns
      call chkParam !sanity check
      return
      end
c-----------------------------------------------------------------------
      subroutine switch_to_lnse_steady!(steps_in)
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'OPCTR'
      include 'CTIMER'

      !integer, intent(in), optional :: steps_in
      integer steps
      real sampling_period

      !if(.not.present(steps_in))then
      !  steps = steps_in
      !else
      !  if(nid.eq.0)write(6,*)' nsteps not specified: computing with uparam(04)'
      !  sampling_period = real(1.0d0/uparam(04))
      !  steps = int(sampling_period/dt)
      !endif

      steps = int(uparam(08))

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

      param(21) = 1e-8 ! pressure tolerance
      param(22) = 1e-10 ! velocity tolerance
      restol(0) = param(22)
      restol(1) = param(22)
      do i=1,ldimt
        restol(1+i) = param(22)
      enddo

      param(31) = 1 !numberOfPerturbations
      param(63) = 0 ! Enforce 32-bit output
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
      subroutine nekStab_forcing (ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none
      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg,iel,ip

      iel=gllel(ieg)
      if (jp.eq.0) then
         ffx = ffx + fcx(ix,iy,iz,iel)
         ffy = ffy + fcy(ix,iy,iz,iel)
         if (if3d) ffz = ffz + fcz(ix,iy,iz,iel)
      endif

      if (spng_str.gt.0) then
      ip=ix+nx1*(iy-1+ny1*(iz-1+nz1*(iel-1)))
         if (jp.eq.0) then
            ! dns
            ffx = ffx + spng_fun(ip)*(spng_vr(ip,1) - vx(ix,iy,iz,iel))*spng_str
            ffy = ffy + spng_fun(ip)*(spng_vr(ip,2) - vy(ix,iy,iz,iel))*spng_str
            if (if3d) ffz = ffz + spng_fun(ip)*(spng_vr(ip,ndim) - vz(ix,iy,iz,iel))*spng_str
         else
            ! perturbation
            ffx = ffx - spng_fun(ip)*vxp(ip,jp)
            ffy = ffy - spng_fun(ip)*vyp(ip,jp)
            if(if3d) ffz = ffz - spng_fun(ip)*vzp(ip,jp)
         endif
      endif
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine nekStab_forcing_temp(temp,ix,iy,iz,ieg)
      implicit none
      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      real temp
      integer ix,iy,iz,ieg,iel,ip

      iel=gllel(ieg) !SUM HERE NEKSTAB CUSTOM TEMP FORCING
      if (jp.eq.0) temp  = temp + fct(ix,iy,iz,iel)
     
      if (spng_str.gt.0) then !!!HERE SPONGE STRENGHT ALWAYS UNITY!
      ip=ix+nx1*(iy-1+ny1*(iz-1+nz1*(iel-1)))
         if (jp.eq.0) then ! dns ! t(1,1,1,1,ifield-1)
            temp = temp + spng_fun(ip)*(spng_vt(ip) - t(ix,iy,iz,iel,1))
         else ! perturbation   ! tp(lpx1*lpy1*lpz1*lpelt,ldimt,lpert)
            temp = temp - spng_fun(ip)*tp(ip,1,jp)
         endif
      endif
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine spng_init
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real glmin,glmax,deltax,xmn,xmx
      real left_sponge,right_sponge
      integer n

      n = nx1*ny1*nz1*nelv
      xmn = glmin(xm1,n); xmx = glmax(xm1,n)
      deltax = xmx - xmn

      left_sponge = real(uparam(08)) !0.1*deltax
      right_sponge = real(uparam(09)) !0.15*deltax

                  spng_wl(1)=(2./3.)*left_sponge ! Sponge left section width; dimension X
                  spng_wl(2)=0.0
          if(IF3D)spng_wl(3)=0.0

                  spng_wr(1)=(2./3.)*right_sponge! Sponge right section width; dimension X
                  spng_wr(2)=0.0
          if(IF3D)spng_wr(3)=0.0

                 spng_dl(1)=(1./3.)*left_sponge ! Sponge left drop/rise section width; dimension X
                 spng_dl(2)=0.0
         if(IF3D)spng_dl(3)=0.0

                spng_dr(1)=(1./3.)*right_sponge ! Sponge right drop/rise section width; dimension X
                spng_dr(2)=0.0
        if(IF3D)spng_dr(3)=0.0

         if(nid.eq.0)then
            write(6,*)' Left sponge width: ',left_sponge
            write(6,*)' Right sponge width: ',right_sponge
            write(6,*)' Sponge left section width: ',spng_wl
            write(6,*)' Sponge right section width: ',spng_wr
            write(6,*)' Sponge left drop/rise section width: ',spng_dl
            write(6,*)' Sponge right drop/rise section width: ',spng_dr
         endif

      ! save reference field -> sponge value reference
      call opcopy(spng_vr(1,1),spng_vr(1,2),spng_vr(1,NDIM),vx,vy,vz) !only DNS
      if(ifheat)call copy(spng_vt,t(1,1,1,1,1),n) !only DNS - temperature
      call spng_set ! -> compute spng_fun

      return
      end
c-----------------------------------------------------------------------
      subroutine spng_set !set sponge function and refernece fields
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real lcoord(LX1*LY1*LZ1*LELV)
      common /SCRUZ/ lcoord

      integer ierr, nhour, nmin, itmp, ntot, il, jl
      real rtmp, ltim, bmin(LDIM), bmax(LDIM)
      real xxmax, xxmax_c, xxmin, xxmin_c, arg
      real glmin, glmax, mth_stepf
      logical ltmp, ltmp2

      ntot = NX1*NY1*NZ1*NELV
      bmin(1) = glmin(XM1,ntot)
      bmax(1) = glmax(XM1,ntot)
      bmin(2) = glmin(YM1,ntot)
      bmax(2) = glmax(YM1,ntot)
      if(IF3D) then
        bmin(NDIM) = glmin(ZM1,ntot)
        bmax(NDIM) = glmax(ZM1,ntot)
      endif

         call rzero(spng_fun,ntot)
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
                      rtmp=1.0d0 !spng_str
                   elseif(rtmp.lt.xxmin) then ! fall; xmin
                      arg = (xxmin-rtmp)/spng_wl(il)
                      rtmp = mth_stepf(arg)
                   elseif (rtmp.le.xxmax) then ! zero
                      rtmp = 0.0
                   elseif (rtmp.lt.xxmax_c) then ! rise
                      arg = (rtmp-xxmax)/spng_wr(il)
                      rtmp = mth_stepf(arg)
                   else    ! constant
                      rtmp = 1.0d0 !spng_str
                   endif
                   spng_fun(jl)=max(spng_fun(jl),rtmp)
                enddo

             endif  ! xxmax.le.xxmin
          endif  ! spng_w(il).gt.0.0
       enddo

      ltmp = ifto; ltmp2 = ifpo
      ifto = .true.; ifpo= .false.
      call outpost2(spng_vr(1,1),spng_vr(1,2),spng_vr(1,NDIM),spng_fun,spng_fun,1,'SPG')
      ifto = ltmp; ifpo = ltmp2

      return
      end
c-----------------------------------------------------------------------
      real function mth_stepf(x) !compute sponge function
      implicit none
      real x, xdmin, xdmax
      parameter (xdmin = 0.0010d0, xdmax = 0.9990d0)
      if (x.le.xdmin) then
         mth_stepf = 0.0d0
      else if (x.le.xdmax) then
         mth_stepf = 1.0d0/( 1.0d0 + exp(1.0d0/(x - 1.0d0) + 1.0d0/x) )
      else
         mth_stepf = 1.0d0
      end if
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
      subroutine add_symmetric_seed(qx, qy, qz) !generate symmetric seed to fields
      implicit none
      include "SIZE"
      include "TOTAL"
      !include 'TSTEP'           ! TIME, DT
      !include 'PARALLEL'        ! LGLEL
      !include 'INPUT'           ! IF3D
      !include 'SOLN'            ! VX, VY, VZ, VMULT
      !include 'GEOM'            ! XM1, YM1, ZM1

      real, dimension(lx1,ly1,lz1,lelv) :: qx, qy, qz
      integer iel,ieg,il,jl,kl,nl,ntot
      real xmn,xmx,xlx,ymn,ymx,yly,zmn,zmx,zlz,alpha,x,y,z
      real glmin,glmax,glsc3,amp

      ntot = NX1*NY1*NZ1*NELV
      ! --> Get the dimensions of the computational box
      xmn = glmin(xm1, ntot)
      xmx = glmax(xm1, ntot)
      xlx = xmx - xmn
 
      ymn = glmin(ym1, ntot)
      ymx = glmax(ym1, ntot)
      yly = ymx - ymn
 
      zmn = glmin(zm1, ntot)
      zmx = glmax(zm1, ntot)
      zlz = zmx - zmn
 
      alpha = 2*pi/zlz
 
      ! --> Create the initial velocity perturbation.
      
      do iel=1,NELV
       do kl=1,NZ1
        do jl=1,NY1
         do il=1,NX1

         ieg = LGLEL(iel)
         x = XM1(il,jl,kl,iel)
         y = YM1(il,jl,kl,iel)
         if (if3D) z = ZM1(il,jl,kl,iel)
 
      ! -> Construct the perturbation. ! Note: Spanwise invariant.
       qx(il,jl,kl,iel) = cos(alpha*z)*sin(2.*pi*y)
       qz(il,jl,kl,iel) = -(2.*pi)/(alpha)*cos(alpha*z)*cos(2.*pi*y)
  
      enddo
      enddo
      enddo
      enddo
 
      amp = glsc3(qx, bm1, qx, ntot)+ glsc3(qy, bm1, qy, ntot)
      if(if3d) amp = amp + glsc3(qz, bm1, qz, ntot)
      amp = 1e-6/(0.50d0*amp)
      call opcmult(qx, qy, qz, amp)

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
      subroutine opaddcol3 (a1,a2,a3,b1,b2,b3,c1,c2,c3)
      implicit none
      include 'SIZE'
      integer ntot1
      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),c1(1),c2(1),c3(1)
      ntot1=nx1*ny1*nz1*nelv
      call addcol3(a1,b1,c1,ntot1)
      call addcol3(a2,b2,c2,ntot1)
      if (ndim.eq.3) call addcol3(a3,b3,c3,ntot1)
      return
      end
c-----------------------------------------------------------------------