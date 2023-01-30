c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      call nekStab

      return
      end
c-----------------------------------------------------------------------
      subroutine nekstab_usrchk
      include 'SIZE'
      include 'TOTAL'

      if(istep.eq.0)then !change defaults
        k_dim = int(uparam(7)) ; call bcast(k_dim,isize)
        schur_tgt = 2 ; call bcast(schur_tgt,isize)
        maxmodes = 10 ; call bcast(maxmodes,isize)
        ifres = .false. ; call bcast(ifres,lsize)
        ifvor = .false. ; call bcast(ifvor,lsize)
        ifvox = .false. ; call bcast(ifvox,lsize)
        ifdyntol = .true. ; call bcast(ifdyntol,lsize)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userf (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ffx = 0.0d0
      ffy = 0.0d0
      ffz = 0.0d0

      call nekStab_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)

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
         ux=ubb(ix,iy,iz,e)
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
      real ramp,pert


      ! following Leshaft JFM 2019 
      !ramp = 1./(1.+exp(4.-0.25*time))
      pert = (1.+0.05*cos(time*8.*atan(1.)*uparam(05)))!*ramp

      if (JP.eq.0) then         ! velocity
         e  = gllel(ieg)
         ux = ubb(ix,iy,iz,e)*pert
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
      subroutine usrdat
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'
      call set_rjet(ubb)
      return
      end

c automatically added by makenek
      subroutine uservp(ix,iy,iz,eg)

      return
      end

c automatically added by makenek
      subroutine userq(ix,iy,iz,eg)

      return
      end

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrdat3 

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
