c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      !if(istep.eq.0) call gfldr("1cyl0.fXXXXX")
      !if(istep.eq.0)t=0.0d0
      call nekStab
      !call hpts

      !call stat_avg
      !call stats_3d_full

      return
      end
c-----------------------------------------------------------------------
      subroutine nekStab_usrchk
      implicit none
      include 'SIZE'
      include 'TOTAL'

      if(istep.eq.0)then !change defaults

        xLspg = uparam(8); call bcast(xLspg , wdsize)
        xRspg = uparam(9); call bcast(xRspg , wdsize)

        schur_tgt = 4 ; call bcast(schur_tgt,isize)
        maxmodes = 2 ; call bcast(maxmodes,isize)
        ifres = .false. ; call bcast(ifres,lsize)
        !ifvor = .true. ; call bcast(ifvor,lsize)
        !ifvox = .true. ; call bcast(ifvox,lsize)

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
         ux= ubb(ix,iy,iz,e)
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
      integer iel,ifc

                 !ID, FIELD, BC
      call setbc(4,1,'v  ')
      call setbc(2,1,'O  ')
      call setbc(3,1,'W  ')

      !this mesh was generatd with genbox
      !with only BCs for velicity
      !here we chenge outflow to dirichlet if adjoint
      !and initilize BCs for scalar field

      if(uparam(1)==3.2)then !if adjoint, change BCs
      do iel=1,nelt
      do ifc = 1, 2*ndim
        if(cbc(ifc,iel,1).eq.'O  ')cbc(ifc,iel,1)='v  '
      enddo
      enddo
      endif

      call compute_inflow(ubb) !compute inflow profile to variables
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_inflow(ubb1) !compute parabolic profile for jet in crossflow
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real ubb1(1)
      integer i
      do i=1,nx1*ny1*nz1*nelv
       ubb1(i)=4.0d0*ym1(i,1,1,1)*(1.0d0-ym1(i,1,1,1))
      enddo
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