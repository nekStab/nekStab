c-----------------------------------------------------------------------
      subroutine nekStab_usrchk
      include 'SIZE'
      include 'TOTAL'
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      udiff =0.0d0
      utrans=0.0d0
      end
c-----------------------------------------------------------------------
      subroutine userf(ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      ffx = 0.0d0
      ffy = 0.0d0
      ffz = 0.0d0
      call makebf_str_pulse(ffx,ffy,ffz,ix,iy,iz,ieg)
      end
c-----------------------------------------------------------------------
      subroutine userq(ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      qvol = 0.0d0
      source = 0.0d0
      end
c-----------------------------------------------------------------------
      subroutine userchk
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'HELIXD'

      if (ISTEP.eq.0) then
         ! omega = (womersley**2)*cpfld(1,1)          ! pulsation frequency
         ! pulse_T = 2.0d0*pi/omega                   ! pulsation period
         ! an_phi = atan2(pitch_s,curv_radius)
         ! sweep = length*cos(an_phi)/curv_radius     ! sweep angle in radians

         call nekstab_init
         ! k_dim = 200
         ! schur_tgt = 2
         ! maxmodes = 2 
         ! ifdyntol = .true.

         ! isNewtonFP  = .true.
         ! call newton_krylov
         ! isNewtonFP  = .false.
   
         ! call linear_stability_analysis('steady','direct')
      endif

      call hpts
      end
c----------------------------------------------------------------------------------
      subroutine userbc(i,j,k,iside,eg)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'HELIXD'
      integer e,eg, i, j, k, iside
      real xo, yo, zo, rr, r2, uaxial
      real angle_t, xstream, ystream, zstream, stream_mag 
      e = gllel(eg)

      ! xo = xax(i,j,k,e)
      ! yo = yax(i,j,k,e)
      ! zo = zax(i,j,k,e)

      ! rr = sqrt(xo**2+yo**2)
      ! r2 = (rr-curv_radius)**2 + zo**2
      ! uaxial = 2.0d0*(1 - r2/(rad**2))
       
      ! angle_t = atan2(xo,yo)
      ! xstream = cos(angle_t)*cos(an_phi)
      ! ystream = -sin(angle_t)*cos(an_phi)
      ! zstream = sin(an_phi)
      ! stream_mag = sqrt(xstream**2+ystream**2+zstream**2)
      ! xstream = xstream/stream_mag
      ! ystream = ystream/stream_mag
      ! zstream = zstream/stream_mag

!      Express velocity components in Cartesian coordinates
      ux = 0.0d0 ! uaxial*xstream
      uy = 0.0d0 !uaxial*ystream
      uz = 0.0d0 !uaxial*zstream
      temp=0.0d0
      end
c-----------------------------------------------------------------------
      subroutine useric(i,j,k,eg)
      implicit none 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer i, j, k, eg
      ux = 0.0d0
      uy = 0.0d0
      uz = 0.0d0
      temp = 0.0d0
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      implicit none
      include 'SIZE'
      include 'TOTAL'
      call setbc(1,1,'W  ')
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'HELIXD'

      if (istep.eq.0) then
         omega = (womersley**2)*cpfld(1,1)      ! pulsation frequency
         pulse_T = 2.0D0*pi/omega               ! pulsation period
         an_phi = atan2(pitch_s,curv_radius)
         sweep = length*cos(an_phi)/curv_radius ! sweep angle in radians
      endif

      call helix_pipe
      if (nid.eq.0) write(6,*) 'USERDAT2: phi = ', an_phi*180./3.14
      if (nid.eq.0) write(6,*) 'USERDAT2: curv_radius = ', curv_radius
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      end
c-----------------------------------------------------------------------

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
