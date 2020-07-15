c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
    
      !if(istep.eq.0) call gfldr("1cyl0.fXXXXX")

      call nekStab

      !call stat_avg
      !call stats_3d_full

      return
      end
c-----------------------------------------------------------------------
      subroutine userf (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      real pra,ray
      pra  = 1.0d0
      ray  = 1720.d0

      ffx = 0.0d0
      ffy = 0.0d0
      if(ifto) ffy = temp*(ray*pra) !coupling the scalar field
      ffz = 0.0d0

      call nekStab_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,idum
      save    idum
      data    idum /99/

      if (JP.eq.0) then         ! velocity
         e  = gllel(ieg)
         ux=0.0d0
         uy=0.0d0
         uz=0.0d0
      
         !ran = 2.e4*(ieg+x*sin(y)) + 1.e3*ix*iy + 1.e5*ix 
         !ran = 1.e3*sin(ran)
         !ran = 1.e3*sin(ran)
         !ran = cos(ran)
         !amp = .0010d0
         temp = 1.0d0-y !+ ran*amp*(1.0d0-y)*y*x*(9.0d0-x) 
         
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
         ux = 0.0d0
         uy = 0.0d0
         uz = 0.0d0
         temp=1.0d0-y
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
      !common /USRBC/ bcflag(lx1,ly1,lz1,lelt)
      !integer iel,ifc,id_face
      !ntot = nx1*ny1*nz1*nelv
      
      !cbc(:,:,2)=cbc(:,:,1)
      !nface = 2*ndim
      !do iel=1,nelt
      ! do iface = 1, nface
      !  if(cbc(iface,iel,1).eq.'W  ')cbc(iface,iel,2)='t  '
      !  if(cbc(iface,iel,1).eq.'v  ')cbc(iface,iel,2)='t  '
      !enddo
      !enddo

      !do iel=1,nelv
      ! do ifc=1,2*ndim
      !  id_face = bc(5,ifc,iel,1)
      !  write(6,*) id_face
      !  if (id_face.eq.1) then                   ! surface 1 for inlet¬
      !     cbc(ifc,iel,1) = 'v  '
      !     cbc(ifc,iel,2) = 't  '
      !  elseif (id_face.eq.2) then               ! surface 2 for outlet¬
      !     cbc(ifc,iel,1) = 'W  '
      !     cbc(ifc,iel,2) = 'I  '
      !  elseif (id_face.eq.3) then               ! surface 3 for upper¬
      !     cbc(ifc,iel,1) = 'O  '
      !     cbc(ifc,iel,2) = 'O  '
      !  elseif (id_face.eq.4) then               ! surface 4 for lower¬
      !     cbc(ifc,iel,1) = 'SYM'
      !     cbc(ifc,iel,2) = 'SYM'
      !  elseif (id_face.eq.5) then               ! surface 5 for wall¬
      !     cbc(ifc,iel,1) = 'SYM'
      !     cbc(ifc,iel,2) = 'SYM'
      !  endif
      ! enddo
      !enddo
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
