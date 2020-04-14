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

      include 'SIZE'
      include 'TOTAL'
      common /USRBC/ bcflag(lx1,ly1,lz1,lelt)
      integer iel,ifc,id_face
      integer isideset1
      ntot = nx1*ny1*nz1*nelv	  
      call rzero(bcflag,ntot)

      do iel=1,nelv
      do ifc=1,2*ndim
        id_face = bc(5,ifc,iel,1)
        if (id_face.eq.1) then                   ! surface 1 for inlet
           cbc(ifc,iel,1) = 'v  '
        elseif (id_face.eq.2) then               ! surface 2 for outlet
           cbc(ifc,iel,1) = 'O  '
        elseif (id_face.eq.3) then               ! surface 3 for upper
           cbc(ifc,iel,1) = 'SYM'
        elseif (id_face.eq.4) then               ! surface 4 for lower
           cbc(ifc,iel,1) = 'SYM'
        elseif (id_face.eq.5) then               ! surface 5 for wall
           cbc(ifc,iel,1) = 'W  '

        endif
      enddo
      enddo


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
