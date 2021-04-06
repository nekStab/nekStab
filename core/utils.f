
      subroutine quicksort2(n, arr, idx)
      implicit none
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
      end subroutine quicksort2
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
      ffx = ffx + fcx(ix,iy,iz,iel)
      ffy = ffy + fcy(ix,iy,iz,iel)
      if (if3d) ffz = ffz + fcz(ix,iy,iz,iel)

      if (spng_str.gt.0) then
         ip=ix+nx1*(iy-1+ny1*(iz-1+nz1*(iel-1)))
         if (jp.eq.0) then
!     dns
            ffx = ffx + spng_fun(ip)*(spng_vr(ip,1) - vx(ix,iy,iz,iel))*spng_str
            ffy = ffy + spng_fun(ip)*(spng_vr(ip,2) - vy(ix,iy,iz,iel))*spng_str
            if (if3d) ffz = ffz + spng_fun(ip)*(spng_vr(ip,ndim) - vz(ix,iy,iz,iel))*spng_str
         else
!     perturbation
            ffx = ffx - spng_fun(ip)*vxp(ip,jp)
            ffy = ffy - spng_fun(ip)*vyp(ip,jp)
            if(if3d) ffz = ffz - spng_fun(ip)*vzp(ip,jp)
         endif
      endif
      return
      end subroutine nekStab_forcing
c-----------------------------------------------------------------------
      subroutine nekStab_forcing_temp(temp,ix,iy,iz,ieg)
      implicit none
      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      real temp
      integer ix,iy,iz,ieg,iel,ip

      iel=gllel(ieg)            !SUM HERE NEKSTAB CUSTOM TEMP FORCING
      if (jp.eq.0) temp  = temp + fct(ix,iy,iz,iel)

      if (spng_str.gt.0) then   !!!HERE SPONGE STRENGHT ALWAYS UNITY!
         ip=ix+nx1*(iy-1+ny1*(iz-1+nz1*(iel-1)))
         if (jp.eq.0) then      ! dns ! t(1,1,1,1,ifield-1)
            temp = temp + spng_fun(ip)*(spng_vt(ip) - t(ix,iy,iz,iel,1))
         else                   ! perturbation   ! tp(lpx1*lpy1*lpz1*lpelt,ldimt,lpert)
            temp = temp - spng_fun(ip)*tp(ip,1,jp)
         endif
      endif
      return
      end subroutine nekStab_forcing_temp
c-----------------------------------------------------------------------
      subroutine spng_init
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real glmin,glmax,deltax
      real left_sponge,right_sponge
      integer n

      n = nx1*ny1*nz1*nelv
      acc_spg = ABS(acc_spg)

      spng_wl(1)=(1.0d0-acc_spg)*xLspg ! Sponge left section width; dimension X
      spng_wl(2)=(1.0d0-acc_spg)*yLspg
      if(IF3D)spng_wl(3)=(1.0d0-acc_spg)*zLspg

      spng_wr(1)=(1.0d0-acc_spg)*xRspg ! Sponge right section width; dimension X
      spng_wr(2)=(1.0d0-acc_spg)*yRspg
      if(IF3D)spng_wr(3)=(1.0d0-acc_spg)*zRspg

      spng_dl(1)=(acc_spg)*xLspg ! Sponge left drop/rise section width; dimension X
      spng_dl(2)=(acc_spg)*yLspg
      if(IF3D)spng_dl(3)=(acc_spg)*zLspg

      spng_dr(1)=(acc_spg)*xRspg !Sponge right drop/rise section width; dimension X
      spng_dr(2)=(acc_spg)*yRspg
      if(IF3D)spng_dr(3)=(acc_spg)*zRspg

      if(nid.eq.0)then
         write(6,*)' Left sponge width: ',left_sponge
         write(6,*)' Right sponge width: ',right_sponge
         write(6,*)' Sponge left section width: ',spng_wl
         write(6,*)' Sponge right section width: ',spng_wr
         write(6,*)' Sponge left drop/rise section width: ',spng_dl
         write(6,*)' Sponge right drop/rise section width: ',spng_dr
      endif

!     save reference field -> sponge value reference
      call opcopy(spng_vr(1,1),spng_vr(1,2),spng_vr(1,NDIM),vx,vy,vz) !only DNS
      if(ifheat)call copy(spng_vt,t(1,1,1,1,1),n) !only DNS - temperature
      call spng_set             ! -> compute spng_fun

      return
      end subroutine spng_init
c-----------------------------------------------------------------------
      subroutine spng_set
!     set sponge function and refernece fields
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
      bmin(1) = xmn             !glmin(XM1,ntot)
      bmax(1) = xmx             !glmax(XM1,ntot)
      bmin(2) = ymn             !glmin(YM1,ntot)
      bmax(2) = ymx             !glmax(YM1,ntot)
      if(IF3D) then
         bmin(NDIM) = zmn       !glmin(ZM1,ntot)
         bmax(NDIM) = zmx       !glmax(ZM1,ntot)
      endif

      call rzero(spng_fun,ntot)
!     for every dimension
      do il=1,NDIM

         if (spng_wl(il).gt.0.0.or.spng_wr(il).gt.0.0) then

            if (spng_wl(il).lt.spng_dl(il).or.spng_wr(il).lt.spng_dr(il)) then
            write(6,*)'Wrong sponge parameters!'
         endif

         xxmax   = bmax(il) - spng_wr(il) ! sponge beginning (rise at xmax; right)
         xxmin   = bmin(il) + spng_wl(il) ! end (drop at xmin; left)
         xxmax_c = xxmax    + spng_dr(il) ! beginnign of constant part (right)
         xxmin_c = xxmin    - spng_dl(il) ! beginnign of constant part (left)

!     get SPNG_FUN
         if (xxmax.le.xxmin) then
            write(6,*)'Sponge too wide'
         else
!     this should be done by pointers, but for now I avoid it
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
                  rtmp=1.0d0    !spng_str
               elseif(rtmp.lt.xxmin) then ! fall; xmin
                  arg = (xxmin-rtmp)/spng_wl(il)
                  rtmp = mth_stepf(arg)
               elseif (rtmp.le.xxmax) then ! zero
                  rtmp = 0.0
               elseif (rtmp.lt.xxmax_c) then ! rise
                  arg = (rtmp-xxmax)/spng_wr(il)
                  rtmp = mth_stepf(arg)
               else             ! constant
                  rtmp = 1.0d0  !spng_str
               endif
               spng_fun(jl)=max(spng_fun(jl),rtmp)
            enddo

         endif                  ! xxmax.le.xxmin
      endif                     ! spng_w(il).gt.0.0
      enddo

      ltmp = ifto; ltmp2 = ifpo
      ifto = .true.; ifpo= .false.
      call outpost2(spng_vr(1,1),spng_vr(1,2),spng_vr(1,NDIM),spng_fun,spng_fun,1,'SPG')
      ifto = ltmp; ifpo = ltmp2

      return
      end subroutine spng_set
c-----------------------------------------------------------------------
      real function mth_stepf(x)
!     compute sponge function
      implicit none
      real x, xdmin, xdmax
      parameter (xdmin = 0.0010d0, xdmax = 0.9990d0)
      if(x.le.xdmin) then
         mth_stepf = 0.0d0
      elseif (x.le.xdmax) then
         mth_stepf = 1.0d0/( 1.0d0 + exp(1.0d0/(x - 1.0d0) + 1.0d0/x) )
      else
         mth_stepf = 1.0d0
      endif
      end function mth_stepf
c-----------------------------------------------------------------------
      subroutine add_noise(qx, qy, qz, qp)
!     input random number to fields
      implicit none
      include 'SIZE'            ! NX1, NY1, NZ1, NELV, NID
      include 'TSTEP'           ! TIME, DT
      include 'PARALLEL'        ! LGLEL
      include 'INPUT'           ! IF3D
      include 'SOLN'            ! VX, VY, VZ, VMULT
      include 'GEOM'            ! XM1, YM1, ZM1

      real, dimension(lx1,ly1,lz1,lelv) :: qx, qy, qz, qp
      integer iel,ieg,il,jl,kl,nl,n
      real xl(LDIM),mth_rand,fcoeff(3) !< coefficients for random distribution
      n = nx1*ny1*nz1*nelv

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

                  if (if3d) then
                     fcoeff(1)= 2.e4;fcoeff(2)= 1.e3;fcoeff(3)= 1.e5
                     qz(il,jl,kl,iel)=qz(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fcoeff)
                  endif

                  if (ifheat) then
                     fcoeff(1)= 9.e4;fcoeff(2)= 3.e3;fcoeff(3)= 4.e5
                     qp(il,jl,kl,iel)=qp(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fcoeff)
                  endif

               enddo
            enddo
         enddo
      enddo

!     face averaging
      call opdssum(qx, qy, qz)
      call opcolv (qx, qy, qz, VMULT)

      call dsavg(qx)
      call dsavg(qy)
      if(if3D) call dsavg(qz)
      call bcdirVC(qx, qy, qz,v1mask,v2mask,v3mask)

      if(ifheat)then
         call dssum(qp,lx1,ly1,lz1)
         call col2(qp, VMULT, n)
         call dsavg(qp)
         call bcdirSC(qp)
      endif

      return
      end subroutine add_noise
c-----------------------------------------------------------------------
      subroutine add_symmetric_seed(qx, qy, qz, qp) !generate symmetric seed to fields
      implicit none
      include "SIZE"
      include "TOTAL"
      real, dimension(lx1,ly1,lz1,lelv) :: qx, qy, qz, qp
      integer iel,ieg,il,jl,kl,nl,ntot
      real xlx,yly,zlz,alpha,x,y,z
      real glmin,glmax,glsc3,amp

      ntot = NX1*NY1*NZ1*NELV
      xlx = xmx - xmn
      yly = ymx - ymn
      zlz = zmx - zmn

      alpha = 2*pi/zlz

!     --> Create the initial velocity perturbation.

      do iel=1,NELV
         do kl=1,NZ1
            do jl=1,NY1
               do il=1,NX1

                  ieg = LGLEL(iel)
                  x = XM1(il,jl,kl,iel)
                  y = YM1(il,jl,kl,iel)
                  if (if3D) z = ZM1(il,jl,kl,iel)

!     -> Construct the perturbation. ! Note: Spanwise invariant.
                  qx(il,jl,kl,iel) = cos(alpha*z)*sin(2.*pi*y)
                  qz(il,jl,kl,iel) = -(2.*pi)/(alpha)*cos(alpha*z)*cos(2.*pi*y)
                  qp(il,jl,kl,iel) = cos(alpha*z)*cos(2.*pi*y)

               enddo
            enddo
         enddo
      enddo

      amp = glsc3(qx, bm1, qx, ntot)+ glsc3(qy, bm1, qy, ntot)
      if(if3d) amp = amp + glsc3(qz, bm1, qz, ntot)
      amp = 1e-6/(0.50d0*amp)
      call opcmult(qx, qy, qz, amp)
      call cmult(qp, amp, ntot)

      return
      end subroutine add_symmetric_seed
c-----------------------------------------------------------------------
      real function mth_rand(ix,iy,iz,ieg,xl,fcoeff) !generate random number
      implicit none
      include 'SIZE'
      include 'INPUT'           ! IF3D
      integer ix,iy,iz,ieg
      real xl(LDIM), fcoeff(3)
      mth_rand = fcoeff(1)*(ieg+xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy+fcoeff(3)*ix
      if (IF3D) mth_rand = fcoeff(1)*(ieg +xl(NDIM)*sin(mth_rand))+fcoeff(2)*iz*ix+fcoeff(3)*iz
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = cos(mth_rand)
      return
      end function mth_rand
c-----------------------------------------------------------------------
      subroutine nopcopy(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)
      implicit none
      include 'SIZE'
      include 'INPUT'
      integer n,n2
      real a1(1),a2(1),a3(1),a4(1),a5(1),b1(1),b2(1),b3(1),b4(1),b5(1)
      n=nx1*ny1*nz1*nelv
      n2=nx2*ny2*nz2*nelv
      call copy(a1,b1,n)
      call copy(a2,b2,n)
      if (ndim.eq.3)call copy(a3,b3,n)
      if (ifpo)call copy(a4,b4,n2)
      if (ifheat)call copy(a5,b5,n)
      return
      end subroutine nopcopy
c-----------------------------------------------------------------------
      subroutine nopsub2(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)
      implicit none
      include 'SIZE'
      include 'INPUT'
      integer n,n2
      real a1(1),a2(1),a3(1),a4(1),a5(1),b1(1),b2(1),b3(1),b4(1),b5(1)
      n=nx1*ny1*nz1*nelv
      n2=nx2*ny2*nz2*nelv
      call sub2(a1,b1,n)
      call sub2(a2,b2,n)
      if (ndim.eq.3)call sub2(a3,b3,n)
      if (ifpo)call sub2(a4,b4,n2)
      if (ifheat)call sub2(a5,b5,n)

      return
      end subroutine nopsub2
c-----------------------------------------------------------------------
      subroutine nopcmult(a1,a2,a3,a4,a5,c)
      implicit none
      include 'SIZE'
      include 'INPUT'
      integer n,n2
      real a1(1),a2(1),a3(1),a4(1),a5(1),c
      n=nx1*ny1*nz1*nelv
      n2=nx2*ny2*nz2*nelv
      call cmult(a1,c,n)
      call cmult(a2,c,n)
      if (ndim.eq.3)call cmult(a3,c,n)
      if (ifpo)call cmult(a4,c,n2)
      if (ifheat)call cmult(a5,c,n)
      return
      end subroutine nopcmult
c-----------------------------------------------------------------------
      subroutine nopchsign(a1,a2,a3,a4,a5)
      implicit none
      include 'SIZE'
      include 'INPUT'
      integer n,n2
      real a1(1),a2(1),a3(1),a4(1),a5(1)
      n=nx1*ny1*nz1*nelv
      n2=nx2*ny2*nz2*nelv
      call chsign(a1,n)
      call chsign(a2,n)
      if (ndim.eq.3)call chsign(a3,n)
      if (ifpo)call chsign(a4,n2)
      if (ifheat)call chsign(a5,n)
      return
      end subroutine nopchsign
c-----------------------------------------------------------------------
      subroutine noprzero(a1,a2,a3,a4,a5)
      implicit none
      include 'SIZE'
      include 'INPUT'
      integer n,n2
      real a1(1),a2(1),a3(1),a4(1),a5(1)
      n=nx1*ny1*nz1*nelv
      n2=nx2*ny2*nz2*nelv
      call rzero(a1,n)
      call rzero(a2,n)
      if (ndim.eq.3)call rzero(a3,n)
      if (ifpo)call rzero(a4,n2)
      if (ifheat)call rzero(a5,n)
      return
      end subroutine noprzero
c-----------------------------------------------------------------------
      subroutine opadd2 (a1,a2,a3,b1,b2,b3)
      implicit none
      include 'SIZE'
      integer ntot1
      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)
      ntot1=nx1*ny1*nz1*nelv
      call add2(a1,b1,ntot1)
      call add2(a2,b2,ntot1)
      if (ndim.eq.3) call add2(a3,b3,ntot1)
      return
      end subroutine opadd2
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
      end subroutine opadd3
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
      end subroutine opaddcol3
c-----------------------------------------------------------------------
      subroutine add4(a,b,c,d,n)
      implicit none
      include 'SIZE'
      integer i,n
      real a(1),b(1),c(1),d(1)
      do i=1,n
         a(i)=b(i)+c(i)+d(i)
      enddo
      return
      end subroutine add4
c-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)
      implicit none
      include 'SIZE'
      integer i,n
      real a(1),b(1),c(1)
      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo
      return
      end subroutine addcol3
c-----------------------------------------------------------------------
      subroutine set_rjet(ub)   !round jet profile for axissymetric jet
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
      end subroutine set_rjet
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
      end subroutine compute_sb
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
      n = nx1*ny1*nz1*nelv

      x_m = uparam(07)

      visc = param(2)/param(1)  !density / dynamic viscosity
      delta99_0 = 5.0d0/1.72080d0 !2.9
      delta_star= 1.0d0
      u_0   = 1.0d0
      x_0 = (delta_star/1.7208d0)**2 / visc * u_0 ! Reference x

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
            dt = .0005          !  Note, this is good to about 12 digits
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
         dt   = .0005           !  Note, this is good to about 12 digits
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
c     1  format(1p6e14.6,' eta')

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
         do k=1,nsteps          !  TIME STEPPING

            call rk4 (w,t,dt,n) ! Single RK4 step (nmax=4)
            t = t+dt

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_ics (w)    !Initial conditions for modified Blasius equation g''' + g g'' = 0
      real  w(0:3)
      w(0) = 0.0d0              ! eta = 0
      w(1) = 0.0d0              ! g
      w(2) = 0.0d0              ! g'
      w(3) = 1.0d0              ! g"
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
 1       if ((khi-klo).gt.1) then
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
