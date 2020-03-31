c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      call nekStab

      return
      end
c-----------------------------------------------------------------------
      subroutine userf (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e

      if(uparam(10).gt.0)call spng_forcing
      
      e = gllel(ieg)
      ffx = fcx(ix,iy,iz,e)
      ffy = fcy(ix,iy,iz,e)
      ffz = fcz(ix,iy,iz,e)

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e

      if (JP.eq.0) then         ! velocity
         e  = gllel(ieg)
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
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e

      if (JP.eq.0) then         ! velocity
         e  = gllel(ieg)
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
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      udiff = 0.0d0
      utrans = 0.0d0
      return
      end
c-----------------------------------------------------------------------
      subroutine userq (ix,iy,iz,ieg)
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
