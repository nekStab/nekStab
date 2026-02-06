      !-----------------------------------------------------------------------
      ! fst.f90 -- Freestream turbulence inflow generation
      !
      ! Purpose:
      !   Generates synthetic turbulent inflow boundary conditions using
      !   a superposition of discrete modes with a von Karman energy
      !   spectrum. Original implementation by M. A. Bucci.
      !
      ! Public interface:
      !   fst              -- main driver (init + update each step)
      !   initWavenumbers  -- read wavenumber data from files
      !   initModes        -- read velocity mode shapes from files
      !   defineBC         -- build pointer to inlet boundary points
      !   interpolateModes -- spline-interpolate modes onto inlet mesh
      !   computeBC        -- generate turbulent velocity at inlet
      !   computeTurbu     -- synthesize turbulent fluctuations
      !   spline           -- cubic spline second-derivative setup
      !   splint           -- cubic spline evaluation
      !
      ! Dependencies:
      !   SIZE, TOTAL
      !-----------------------------------------------------------------------

      module nekstab_fst
         private
         public :: fst, initWavenumbers, initModes,
     $             defineBC, interpolateModes, computeBC,
     $             computeTurbu, spline, splint

         include 'SIZE'

         ! FST parameters (override in SIZE if using FST)
         integer, parameter :: fst_numk_default = 20
         integer, parameter :: fst_nmodes_default = 10
         integer, parameter :: fst_maxpoints = 9999

         ! FST field arrays
         integer :: fst_numk, fst_nmodes, npointModes, npointBC
         real :: fst_okini, fst_okfin, fst_length, fst_tu
         real :: fst_uin(lx1,ly1,lz1,lelv)
         real :: fst_vin(lx1,ly1,lz1,lelv)
         real :: fst_win(lx1,ly1,lz1,lelv)
         real :: frec(2,fst_numk_default*fst_nmodes_default)
         real :: pointBC(fst_maxpoints,7)
         real :: umodes(fst_maxpoints,7,fst_numk_default*fst_nmodes_default)
         real :: umodesBC(fst_maxpoints,6,fst_numk_default*fst_nmodes_default)

      contains

      !-----------------------------------------------------------------------
      ! fst -- Main driver for freestream turbulence generation
      !
      ! Purpose:
      !   On first call initializes wavenumbers, modes, inlet points,
      !   and interpolation. Every call updates the turbulent BC.
      !-----------------------------------------------------------------------
      subroutine fst
      ! Freestream turbulence (original implementation by M A Bucci)
         include 'SIZE'
         include 'TOTAL'
      
         if (istep == 0) then
            call oprzero(fst_uin, fst_vin, fst_win)
            call initWavenumbers
            call initModes
            call defineBC
            call interpolateModes
         end if
      
         call computeBC
      
      end subroutine fst

      !-----------------------------------------------------------------------
      ! initWavenumbers -- Read wavenumber/frequency data from FST_data/
      !-----------------------------------------------------------------------
      subroutine initWavenumbers
         include 'SIZE'
         include 'TOTAL'
         integer :: itervp, i
         character(len=60) :: filename
         itervp = 0
         do i = 1, fst_numk*fst_nmodes
            itervp = itervp + 1
            write (filename, '(A,I3.3,A)') 'FST_data/wavenumber', itervp, '.dat'
            open (299, file=trim(filename), form='formatted', status='old')
            read (299, *) frec(1, i) !omega
            read (299, *)
            read (299, *) frec(2, i) !beta
            close (299)
         end do
      end subroutine initWavenumbers

      !-----------------------------------------------------------------------
      ! initModes -- Read velocity mode shapes from FST_data/
      !-----------------------------------------------------------------------
      subroutine initModes
         include 'SIZE'
         include 'TOTAL'
         integer :: i, j, k, itervp
         character(len=60) :: filename
         itervp = 0
         do k = 1, fst_numk
            do j = 1, fst_nmodes
               itervp = itervp + 1
               write (filename, '(A,I3.3,A)') 'FST_data/velocity', itervp, '.dat'
               open (200, file=trim(filename), form='formatted', status='unknown')
               read (200, *) npointModes
               do i = 1, npointModes
                  read (200, *) umodes(i, 1, itervp), umodes(i, 2, itervp), umodes(i, 3, itervp),
     $   umodes(i, 4, itervp), umodes(i, 5, itervp), umodes(i, 6, itervp), umodes(i, 7, itervp)
               end do
               close (200)
            end do
         end do
      end subroutine initModes

      !-----------------------------------------------------------------------
      ! defineBC -- Identify inlet boundary points for FST injection
      !-----------------------------------------------------------------------
      subroutine defineBC
         include 'SIZE'
         include 'TOTAL'
         include 'NEKUSE'
         integer i, j, k, e, f, i0, i1, j0, j1, k0, k1
      
         npointBC = 0
         do e = 1, nelv
            do f = 1, 2*ndim
      
               if (cbc(f, e, 1) == 'v  ') then
                  call facind(i0, i1, j0, j1, k0, k1, nx1, ny1, nz1, f)
                  do i = i0, i1
                     do j = j0, j1
                        do k = k0, k1
                           npointBC = npointBC + 1
                           pointBC(npointBC, 1) = i
                           pointBC(npointBC, 2) = j
                           pointBC(npointBC, 3) = k
                           pointBC(npointBC, 4) = e
                           pointBC(npointBC, 5) = xm1(i, j, k, e)
                           pointBC(npointBC, 6) = ym1(i, j, k, e)
                           pointBC(npointBC, 7) = zm1(i, j, k, e)
                        end do
                     end do
                  end do
               end if
      
            end do
         end do
      end subroutine defineBC

      !-----------------------------------------------------------------------
      ! interpolateModes -- Spline-interpolate mode shapes onto inlet mesh
      !-----------------------------------------------------------------------
      subroutine interpolateModes
         include 'SIZE'
         include 'TOTAL'
         real y1(npointModes), y2(npointModes), y3(npointModes),
     $   y4(npointModes), y5(npointModes), y6(npointModes),
     $   uint1, uint2, uint3, uint4, uint5, uint6, yint
         integer itervp, kvalue, j, i
      
         itervp = 0
         do kvalue = 1, fst_numk
            do j = 1, fst_nmodes
               itervp = itervp + 1
      
      !SPLINE SECOND DERIVATE
               call spline(umodes(1, 1, 1), umodes(1, 2, itervp), npointModes, 1e30, 1e30, y1(1))
               call spline(umodes(1, 1, 1), umodes(1, 3, itervp), npointModes, 1e30, 1e30, y2(1))
               call spline(umodes(1, 1, 1), umodes(1, 4, itervp), npointModes, 1e30, 1e30, y3(1))
               call spline(umodes(1, 1, 1), umodes(1, 5, itervp), npointModes, 1e30, 1e30, y4(1))
               call spline(umodes(1, 1, 1), umodes(1, 6, itervp), npointModes, 1e30, 1e30, y5(1))
               call spline(umodes(1, 1, 1), umodes(1, 7, itervp), npointModes, 1e30, 1e30, y6(1))
      
      !     INTERPOLATING VALUES
               do i = 1, npointBC
                  yint = pointBC(i, 6)
                  call splint(umodes(1, 1, 1), umodes(1, 2, itervp), y1(1), npointModes, yint, uint1)
                  call splint(umodes(1, 1, 1), umodes(1, 3, itervp), y2(1), npointModes, yint, uint2)
                  call splint(umodes(1, 1, 1), umodes(1, 4, itervp), y3(1), npointModes, yint, uint3)
                  call splint(umodes(1, 1, 1), umodes(1, 5, itervp), y4(1), npointModes, yint, uint4)
                  call splint(umodes(1, 1, 1), umodes(1, 6, itervp), y5(1), npointModes, yint, uint5)
                  call splint(umodes(1, 1, 1), umodes(1, 7, itervp), y6(1), npointModes, yint, uint6)
                  umodesBC(i, 1, itervp) = uint1
                  umodesBC(i, 2, itervp) = uint2
                  umodesBC(i, 3, itervp) = uint3
                  umodesBC(i, 4, itervp) = uint4
                  umodesBC(i, 5, itervp) = uint5
                  umodesBC(i, 6, itervp) = uint6
               end do
            end do
         end do
      end subroutine interpolateModes

      !-----------------------------------------------------------------------
      ! computeBC -- Generate turbulent velocity field at inlet points
      !-----------------------------------------------------------------------
      subroutine computeBC
         include 'SIZE'
         include 'TOTAL'
         real u_turbu(npointBC, 3)
         integer i, j, k, e, ieg, index
         e = gllel(ieg)
      
         call computeTurbu(u_turbu) ! Generate the turbulent profile
      
         do index = 1, npointBC
            i = pointBC(index, 1)
            j = pointBC(index, 2)
            k = pointBC(index, 3)
            e = pointBC(index, 4)
      
            fst_uin(i, j, k, e) = u_turbu(index, 1)
            fst_vin(i, j, k, e) = u_turbu(index, 2)
            fst_win(i, j, k, e) = u_turbu(index, 3)
         end do
      
      end subroutine computeBC

      !-----------------------------------------------------------------------
      ! computeTurbu -- Synthesize turbulent fluctuations from mode basis
      !
      ! Arguments:
      !   u_turbu [out] -- turbulent velocity at inlet points (npointBC,3)
      !-----------------------------------------------------------------------
      subroutine computeTurbu(u_turbu)
         include 'SIZE'
         include 'TOTAL'
         real u_turbu(npointBC, 3), turbu_aux(npointBC, 3)
         real kk, dkk, dkke, kk1, kk2, enspect1, enspect2, integral
         real ampli, enspect, auxcos(npointBC), auxsin(npointBC), aa, bb
         integer i, j, k, itervp, kvalue
      
         aa = 1.6060d0 ! von Karman Spectrum constant 1
         bb = 1.3500d0 ! von Karman Spectrum constant 2
         dkk = (fst_okfin - fst_okini)/(fst_numk - 1)
         kk1 = fst_okini - dkk/2
         kk2 = fst_okfin + dkk/2
         dkke = (kk2 - kk1)/fst_numk
         kk = kk1
         integral = 0.0d0
      
         do i = 1, fst_numk
            enspect1 = 2.0d0/3.0d0*fst_length*((aa*(kk*fst_length)**4)/(bb + (kk*fst_length)**2)**(17.0d0/6.0d0))
            kk = kk + dkke
            enspect2 = 2.0d0/3.0d0*fst_length*((aa*(kk*fst_length)**4)/(bb + (kk*fst_length)**2)**(17.0d0/6.0d0))
            integral = integral + (enspect1 + enspect2)*dkke/2
         end do
         integral = 1.0d0/integral
         itervp = 0.0d0
      
         call rzero(u_turbu(1, 1), npointBC)
         call rzero(u_turbu(1, 2), npointBC)
         call rzero(u_turbu(1, 3), npointBC)
      
         kk = fst_okini
         do kvalue = 1, fst_numk
      
            do i = 1, fst_nmodes
               enspect = integral*fst_tu**2*fst_length*((aa*(kk*fst_length)**4)/(bb + (kk*fst_length)**2)**(17.0d0/6.0d0))
               ampli = sqrt(enspect*dkk/(fst_nmodes*2)*2)
               itervp = itervp + 1
      
               do j = 1, npointBC
                  auxcos(j) = +cos(+frec(1, itervp)*time + frec(2, itervp)*pointBC(j, 7))
     $   +cos(-frec(1, itervp)*time + frec(2, itervp)*pointBC(j, 7))
                  auxsin(j) = -sin(+frec(1, itervp)*time + frec(2, itervp)*pointBC(j, 7))
     $   -sin(-frec(1, itervp)*time + frec(2, itervp)*pointBC(j, 7))
               end do
      
               call rzero(turbu_aux(1, 1), npointBC)
               call addcol3(turbu_aux(1, 1), umodesBC(1, 1, itervp), auxcos(1), npointBC) !a = a+b*c
               call addcol3(turbu_aux(1, 1), umodesBC(1, 2, itervp), auxsin(1), npointBC)
               call rzero(turbu_aux(1, 2), npointBC)
               call addcol3(turbu_aux(1, 2), umodesBC(1, 3, itervp), auxcos(1), npointBC)
               call addcol3(turbu_aux(1, 2), umodesBC(1, 4, itervp), auxsin(1), npointBC)
               call rzero(turbu_aux(1, 3), npointBC)
               call addcol3(turbu_aux(1, 3), umodesBC(1, 5, itervp), auxcos(1), npointBC)
               call addcol3(turbu_aux(1, 3), umodesBC(1, 6, itervp), auxsin(1), npointBC)
      
               call add2s2(u_turbu(1, 1), turbu_aux(1, 1), ampli, npointBC)
               call add2s2(u_turbu(1, 2), turbu_aux(1, 2), ampli, npointBC)
               call add2s2(u_turbu(1, 3), turbu_aux(1, 3), ampli, npointBC)
      
            end do
            kk = kk + dkk
         end do
      
      end subroutine computeTurbu
      ! !----------------------------------------------------------------------
      !       SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      !       implicit none
      !       INTEGER i,k,n,NMAX
      !       REAL yp1,ypn,x(n),y(n),y2(n)
      !       PARAMETER (NMAX=9999)
      !       REAL p,qn,sig,un,u(NMAX)
      
      !       if (yp1.gt..99e30) then ! NATURAL
      !          y2(1)=0.0d0
      !          u(1)=0.0d0
      !       else ! FIRST DERIVATIVE
      !          y2(1)=-0.5
      !          u (1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      !       endif
      !       do  i=2,n-1
      !          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      !          p=sig*y2(i-1)+2.
      !          y2(i)=(sig-1.)/p
      !          u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      !       enddo
      !       if (ypn.gt..99e30) then
      !         qn=0.0d0
      !         un=0.0d0
      !       else
      !         qn=0.50d0
      !         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      !       endif
      !       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      !       do  k=n-1,1,-1
      !       y2(k)=y2(k)*y2(k+1)+u(k)
      !       enddo
      !       return
      !       END SUBROUTINE SPLINE
      ! !----------------------------------------------------------------------
      !       SUBROUTINE splint(xa,ya,y2a,n,x,y)
      !       implicit none
      !       INTEGER n,k,khi,klo
      !       REAL a,b,h,x,y,xa(n),y2a(n),ya(n)
      !       klo=1; khi=n
      
      ! 1     if (khi-klo.gt.1) then
      !          k=(khi+klo)/2
      !          if(xa(k).gt.x)then
      !              khi=k
      !           else
      !              klo=k
      !           endif
      !       goto 1
      !       endif
      !       h=xa(khi)-xa(klo)
      !       if (h.eq.0.) then
      !          !if (nid.eq.0) print *, 'bad xa input in splint'
      !          call exitt
      !       endif
      !       a=(xa(khi)-x)/h
      !       b=(x-xa(klo))/h
      !       y=a*ya(klo)+b*ya(khi)+
      !     $                        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      
      !       return
      !       END SUBROUTINE splint
      !-----------------------------------------------------------------------
      ! spline -- Cubic spline second-derivative setup
      !
      ! Arguments:
      !   x, y   [in]  -- data arrays (n points)
      !   n      [in]  -- number of data points
      !   yp1    [in]  -- first derivative at x(1) (>0.99e30 = natural)
      !   ypn    [in]  -- first derivative at x(n) (>0.99e30 = natural)
      !   y2     [out] -- second derivatives at data points
      !-----------------------------------------------------------------------
      subroutine spline(x, y, n, yp1, ypn, y2)
         implicit none
         integer, intent(in) :: n
         real, intent(in) :: yp1, ypn, x(n), y(n)
         real, intent(out) :: y2(n)
         real, allocatable, dimension(:) :: u
         real :: dx_forward, dx_backward, dx_total, dy_forward_norm, dy_backward_norm
         real :: p, sig, qn, un
         integer :: i
      
         allocate (u(n))
      
         if (yp1 > 0.99e30) then ! Natural spline
            y2(1) = 0.0
            u(1) = 0.0
         else ! First derivative
            y2(1) = -0.5
            u(1) = (3.0/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
         end if
      
         do i = 2, n - 1
            dx_forward = x(i + 1) - x(i)
            dx_backward = x(i) - x(i - 1)
            dx_total = x(i + 1) - x(i - 1)
      
            sig = dx_backward/dx_total
      
            p = sig*y2(i - 1) + 2.0
      
            y2(i) = (sig - 1.0)/p
      
            dy_forward_norm = (y(i + 1) - y(i))/dx_forward
            dy_backward_norm = (y(i) - y(i - 1))/dx_backward
      
            u(i) = (6.0*(dy_forward_norm - dy_backward_norm)/dx_total - sig*u(i - 1))/p
         end do
      
         if (ypn > 0.99e30) then
            qn = 0.0
            un = 0.0
         else
            qn = 0.50
            un = (3.0/(x(n) - x(n - 1)))*(ypn - (y(n) - y(n - 1))/(x(n) - x(n - 1)))
         end if
      
         y2(n) = (un - qn*u(n - 1))/(qn*y2(n - 1) + 1.0)
      
         do i = n - 1, 1, -1
            y2(i) = y2(i)*y2(i + 1) + u(i)
         end do
      
         deallocate (u)
      end subroutine spline

      !-----------------------------------------------------------------------
      ! splint -- Cubic spline evaluation at a single point
      !
      ! Arguments:
      !   xa, ya [in]  -- data arrays (n points)
      !   y2a    [in]  -- second derivatives from spline()
      !   n      [in]  -- number of data points
      !   x      [in]  -- evaluation point
      !   y      [out] -- interpolated value
      !-----------------------------------------------------------------------
      subroutine splint(xa, ya, y2a, n, x, y)
         implicit none
         integer, intent(in) :: n
         real, intent(in) :: x, xa(n), ya(n), y2a(n)
         real, intent(out) :: y
         integer :: k, khi, klo
         real :: a, b, h
      
         klo = 1
         khi = n
      
         do
            if (khi - klo > 1) then
               k = (khi + klo)/2
               if (xa(k) > x) then
                  khi = k
               else
                  klo = k
               end if
            else
               exit
            end if
         end do
      
         h = xa(khi) - xa(klo)
      
         if (h == 0.0) then
            call exitt
         end if
      
         a = (xa(khi) - x)/h
         b = (x - xa(klo))/h
      
         y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*h**2/6.0

      end subroutine splint

      end module nekstab_fst
