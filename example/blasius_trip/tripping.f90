
!! Tripping function supporting conformal and AMR version of nek5000
!!  The tripping is based on a similar implementation in the
!!   SIMSON code (Chevalier et al. 2007, KTH Mechanics), and is described
!!   in detail in the paper Schlatter & Örlü, JFM 2012, DOI 10.1017/jfm.2012.324.
!! @author Adam Peplinski
!! @date May 03, 2018
!=======================================================================

      subroutine tripl_init()
         implicit none
         include 'SIZE'
         include 'INPUT'
         include 'GEOM'
         include 'TRIP'
         integer il, jl, kl
         real :: rtmp, rtmpv(ldim)
         logical :: ltmp
         integer :: itmp

         if (nid == 0) write (6, *) 'Initialisation started'
         ! Check if the module was already initialized
         if (tripl_ifinit) then
            if (nid == 0) write (6, *) 'module ['//trim(tripl_name)//'] already initialized.'
            return
         endif

         tripl_nline = 1  ! Number of tripping lines

         ! Set parameters for tripping lines
         ! Line 1
         tripl_tiamp(1) = 0.0d0 ! Time independent amplitude
         tripl_tdamp(1) = 10.0d0 ! Time dependent amplitude
         tripl_spos(1, 1) = 8.300d0 ! Starting point X
         tripl_spos(2, 1) = 0.05d0 ! Starting point Y
         tripl_spos(3, 1) = 0.000d0 ! Starting point Z
         tripl_epos(1, 1) = 8.300d0 ! Ending point X
         tripl_epos(2, 1) = 0.05d0 ! Ending point Y
         tripl_epos(3, 1) = 4.500d0 ! Ending point Z
         tripl_smth(1, 1) = 1.360d0 ! Smoothing length X
         tripl_smth(2, 1) = 0.340d0 ! Smoothing length Y
         tripl_smth(3, 1) = 0.300d0 ! Smoothing length Z
         tripl_lext(1) = .true. ! Line extension
         tripl_rota(1) = 0.0d0 ! Rotation angle in radians
         tripl_nmode(1) = 16 ! Number of Fourier modes
         tripl_tdt(1) = 0.140 ! Time step for tripping

         ! ! Line 2
         ! tripl_tiamp(2) = 0.00000000E+00 ! Time independent amplitude
         ! tripl_tdamp(2) = 10.0d0 ! Time dependent amplitude
         ! tripl_spos(1, 2) = 8.80000000E+00 ! Starting point X
         ! tripl_spos(2, 2) = 0.99430000E+00 ! Starting point Y
         ! tripl_spos(3, 2) = 0.00000000E+00 ! Starting point Z
         ! tripl_epos(1, 2) = 8.80000000E+00 ! Ending point X
         ! tripl_epos(2, 2) = 0.99430000E+00 ! Ending point Y
         ! tripl_epos(3, 2) = 4.50000000E+00 ! Ending point Z
         ! tripl_smth(1, 2) = 0.28000000E+00 ! Smoothing length X
         ! tripl_smth(2, 2) = 0.07000000E+00 ! Smoothing length Y
         ! tripl_smth(3, 2) = 0.20000000E+00 ! Smoothing length Z
         ! tripl_lext(2) = .false. ! Line extension
         ! tripl_rota(2) = 0.12300000E+00 ! Rotation angle in radians
         ! tripl_nmode(2) = 76 ! Number of Fourier modes
         ! tripl_tdt(2) = 0.14000000E+00 ! Time step for tripping

         ! Check simulation dimension
         if (.not.IF3D) then
            if (nid == 0) write (6, *) '2D simulation is not supported.'
            call nek_end
         endif

         ! Calculate line versors, inverse line lengths, and scaled smoothing lengths
         do il = 1, tripl_nline
            call mntr_logi(tripl_id, lp_inf, 'Line info; line nr: ', il)
            tripl_ilngt(il) = 0.0

            do jl = 1, LDIM
                  ! The last (third) versor is parallel to the line
                  tripl_vrs(jl, ldim, il) = tripl_epos(jl, il) - tripl_spos(jl, il)
                  tripl_ilngt(il) = tripl_ilngt(il) + tripl_vrs(jl, ldim, il)**2
            end do
            tripl_ilngt(il) = sqrt(tripl_ilngt(il))
            if (nid == 0) write (6, *) 'Line length: ', tripl_ilngt(il)
            if (tripl_ilngt(il) > 0.0) then
                  tripl_ilngt(il) = 1.0 / tripl_ilngt(il)
                  do jl = 1, LDIM
                     tripl_vrs(jl, ldim, il) = tripl_vrs(jl, ldim, il) * tripl_ilngt(il)
                  end do
            else
                  if (nid == 0) write (6, *) 'Line with zero length is not supported.'
            endif

            ! The rest of versors given by cross product starting with the second versor
            call rzero(rtmpv, ldim)
            ! The first versor guess depends on the last versor coordinates
            if (tripl_vrs(1, ldim, il) < 0.95) then
                  rtmpv(1) = 1.0
            else
                  rtmpv(2) = 1.0
            endif
            ! The second versor
            call cross(tripl_vrs(1, 2, il), tripl_vrs(1, ldim, il), rtmpv)
            ! The first versor
            call cross(tripl_vrs(1, 1, il), tripl_vrs(1, 2, il), tripl_vrs(1, ldim, il))
            ! Correct versor length
            do jl = 1, ldim
                  rtmp = 0.0
                  do kl = 1, ldim
                     rtmp = rtmp + tripl_vrs(kl, jl, il) * tripl_vrs(kl, jl, il)
                  end do
                  if (rtmp > 0.0) then
                     rtmp = 1.0 / sqrt(rtmp)
                     do kl = 1, ldim
                        tripl_vrs(kl, jl, il) = tripl_vrs(kl, jl, il) * rtmp
                     end do
                  else
                     if (nid == 0) write (6, *) 'Line versor with zero length.'
                  end if
            end do
            ! Rotate the first and second versors along the third versor
            do jl = 1, 2
                  call math_rot3da(rtmpv, tripl_vrs(1, jl, il), tripl_vrs(1, ldim, il), tripl_rota(il))
                  do kl = 1, ldim
                     tripl_vrs(kl, jl, il) = rtmpv(kl)
                  end do
            end do

            ! Stump the log
            do jl = 1, ldim
                  if (nid == 0) write (6, *) 'Line versor: ', jl
                  if (nid == 0) write (6, *) 'Coordinates:', tripl_vrs(1, jl, il), ldim
            end do

            ! Rescale smoothing lengths
            do jl = 1, LDIM
                  tripl_smth(jl, il) = tripl_smth(jl, il) * tripl_ilngt(il)
            end do
            ! Get inverse smoothing length
            do jl = 1, LDIM
                  if (tripl_smth(jl, il) > 0.0) then
                     tripl_ismth(jl, il) = 1.0 / tripl_smth(jl, il)
                  else
                     tripl_ismth(jl, il) = 1.0
                  endif
            end do
         end do
         if (nid == 0) write (6, *) 'Local base calculated'

         ! Get 1D projection and array mapping
         call tripl_1dprj
         if (nid == 0) write (6, *) '1D projection finalized'

         ! Initialize random generator seed and number of time intervals
         do il = 1, tripl_nline
            tripl_seed(il) = -32 * il
            tripl_ntdt(il) = 1 - tripl_nset_max
            tripl_ntdt_old(il) = tripl_ntdt(il)
         end do

         ! Generate random phases (time independent and time dependent)
         il = tripl_nmode_max * tripl_nset_max * tripl_nline_max
         call rzero(tripl_rphs, il)
         call tripl_rphs_get
         if (nid == 0) write (6, *) 'Random phases calculated'

         ! Get forcing
         call tripl_frcs_get(.true.)
         if (nid == 0) write (6, *) 'Forcing calculated'

         ! Everything is initialized
         tripl_ifinit = .true.

         if (nid == 0) write (6, *) 'Initialization finalized'

         return
      end subroutine tripl_init
!=======================================================================
      logical function tripl_is_initialised()
      implicit none
      include 'SIZE'
      include 'TRIP'
      tripl_is_initialised = tripl_ifinit
      return
      end function
!=======================================================================

      subroutine tripl_update() ! Update tripping
      implicit none
      include 'SIZE'
      include 'TRIP'
      ! update random phases (time independent and time dependent)
      call tripl_rphs_get
      ! update forcing
      call tripl_frcs_get(.false.)

      return
      end subroutine      
!=======================================================================
      subroutine tripl_forcing(ffx,ffy,ffz,ix,iy,iz,ieg) ! Compute tripping forcing
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'TRIP'

      ! argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

      ! local variables
      integer ipos,iel,il
      real ffn
!-----------------------------------------------------------------------
      iel=GLLEL(ieg)

      do il= 1, tripl_nline
         ffn = tripl_fsmth(ix,iy,iz,iel,il)
         if (ffn.gt.0.0) then
            ipos = tripl_map(ix,iy,iz,iel,il)
            ffn = tripl_ftrp(ipos,il)*ffn

            ! I assume forcing direction is given by the second versor
            ffx = ffx + ffn*tripl_vrs(1,2,il)
            ffy = ffy + ffn*tripl_vrs(2,2,il)
            ffz = ffz + ffn*tripl_vrs(ldim,2,il)
         endif
      enddo

      return
      end subroutine
!=======================================================================

      subroutine tripl_reset() ! Reset tripping
      implicit none
      include 'SIZE'
      include 'TRIP'

      ! get 1D projection and array mapping
      call tripl_1dprj
      ! update forcing
      call tripl_frcs_get(.true.)

      return
      end subroutine
!=======================================================================
!> Get 1D projection, array mapping and forcing smoothing
!! @ingroup tripl
!!  This routine supports straight lines given by their starting
!!    and ending points. Additional flagg allows to introuduce forcing
!!    periodicity or contain it between starting and ending points + smooting
!!    lenght in z
!! @remark This routine uses global scratch space \a CTMP0 and \a CTMP1
      subroutine tripl_1dprj()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TRIP'

      ! global memory access
      real lcoord(LX1*LY1*LZ1*LELT)
      common /CTMP0/ lcoord
      integer lmap(LX1*LY1*LZ1*LELT), lmap_el(4,LX1*LY1*LZ1*LELT)
      common /CTMP1/ lmap, lmap_el

      ! local variables
      integer nv, itmp, itmp2
      integer il, jl, kl, ll, ml, nl
      real rota, rtmp, epsl
      parameter (epsl = 1.0e-10)
      real rtmpv(ldim), rtmpc(ldim)
      ! functions
      real dot
!-----------------------------------------------------------------------
      nv = nx1*ny1*nz1*nelv

      ! for each line
      do il=1,tripl_nline
         ! reset mapping array
         call ifill(tripl_map(1,1,1,1,il),-1,nv)
         ! initialise number of points per line
         tripl_npoint(il) = 0
         ! initialize smoothing factor
         call rzero(tripl_fsmth(1,1,1,1,il),nv)
         ! initialize projected point position
         call rzero(tripl_prj(1,il),nv)

         ! Projection onto the line
         ! count points on the line
         itmp = 0
         do jl = 1, nelv
            do kl = 1, lz1
               do ll = 1, ly1
                  do ml = 1, lx1
                     ! get point position relative to the line start
                     rtmpv(1) = xm1(ml,ll,kl,jl)-tripl_spos(1,il)
                     rtmpv(2) = ym1(ml,ll,kl,jl)-tripl_spos(2,il)
                     rtmpv(ldim) = zm1(ml,ll,kl,jl)-tripl_spos(ldim,il)
                     ! get point coordinates in the local line system
                     do nl = 1, ldim
                        rtmpc(nl) = dot(rtmpv,tripl_vrs(1,nl,il),ldim)
                        rtmpc(nl) = rtmpc(nl)*tripl_ilngt(il)
                     end do
                     ! distance from the line
                     ! 2D
                     rtmp = (rtmpc(1)*tripl_ismth(1,il))**2+(rtmpc(2)*tripl_ismth(2,il))**2
                     ! do we extend a line beyond its ends
                     if (.not.tripl_lext(il)) then
                        if (rtmpc(ldim).lt.0.0) then
                           rtmp = rtmp +(rtmpc(ldim)*tripl_ismth(ldim,il))**2
                        elseif (rtmpc(ldim).gt.1.0) then
                           rtmp = rtmp +((rtmpc(ldim)-1.0)*tripl_ismth(ldim,il))**2
                        end if
                     end if

                     ! get smoothing profile
                     ! Gauss; cannot be used with lines not extended beyond their ending points
                     !tripl_fsmth(itmp,jtmp,ktmp,eltmp,il) = exp(-4.0*rtmp)
                     ! limited support
                     if (rtmp.lt.1.0) then
                        tripl_fsmth(ml,ll,kl,jl,il) =exp(-rtmp)*(1-rtmp)**2
                        ! add the point to the list
                        itmp = itmp + 1
                        ! save data
                        ! coordinate along the line in line length unit
                        lcoord(itmp) = rtmpc(ldim)
                        ! point position
                        lmap_el(1,itmp) = jl
                        lmap_el(2,itmp) = kl
                        lmap_el(3,itmp) = ll
                        lmap_el(4,itmp) = ml
                     else
                        tripl_fsmth(ml,ll,kl,jl,il) = 0.0d0
                     endif
                  end do
               end do
            end do
         end do

         if (itmp.ge.1) then
            ! point sorting acording to the last coordinate
            call sort(lcoord,lmap,itmp)

            ! identify unique points
            tripl_npoint(il) = 1
            tripl_prj(tripl_npoint(il),il) = lcoord(1)
            ! generate mapping
            itmp2 = lmap(1)
            tripl_map(lmap_el(4,itmp2),lmap_el(3,itmp2),lmap_el(2,itmp2),lmap_el(1,itmp2),il) = tripl_npoint(il)
            do jl = 2, itmp
               ! compare positions along the line
               if((lcoord(jl)-tripl_prj(tripl_npoint(il),il)).gt.max(epsl,abs(epsl*lcoord(jl)))) then
                  tripl_npoint(il) = tripl_npoint(il) + 1
                  tripl_prj(tripl_npoint(il),il) = lcoord(jl)
               endif
               ! generate mapping
               itmp2 = lmap(jl)
               tripl_map(lmap_el(4,itmp2),lmap_el(3,itmp2),lmap_el(2,itmp2),lmap_el(1,itmp2),il) = tripl_npoint(il)
            end do
         end if
      enddo

      return
      end subroutine      
!=======================================================================
      subroutine tripl_rphs_get ! Generate set of random phases
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'TRIP'
      
      ! local variables
      integer il, jl, kl
      integer itmp
      real tripl_ran2

!-----------------------------------------------------------------------
      ! time independent part
      do il = 1, tripl_nline
         if (tripl_tiamp(il).gt.0.0.and..not.tripl_ifinit) then
            do jl=1, tripl_nmode(il)
               tripl_rphs(jl,1,il) = 2.0*pi*tripl_ran2(il)
            end do
         end if
      end do

      ! time dependent part
      do il = 1, tripl_nline
         itmp = int(time/tripl_tdt(il))
         call bcast(itmp,ISIZE) ! just for safety
         do kl= tripl_ntdt(il)+1, itmp
            do jl= tripl_nset_max,3,-1
               call copy(tripl_rphs(1,jl,il),tripl_rphs(1,jl-1,il),tripl_nmode(il))
            enddo
            do jl=1, tripl_nmode(il)
               tripl_rphs(jl,2,il) = 2.0*pi*tripl_ran2(il)
            enddo
         enddo
         ! update time interval
         tripl_ntdt_old(il) = tripl_ntdt(il)
         tripl_ntdt(il) = itmp
      enddo

      return
      end subroutine
!=======================================================================

      real function tripl_ran2(il) ! Random number generator
!! @param[in]   il      line number
!! @return      ran
!!   Requires 32-bit integer arithmetic. Taken from Numerical
!!   Recipes, William Press et al. Gives correlation free random
!!   numbers but does not have a very large dynamic range, i.e only
!!   generates 714025 different numbers. Set seed negative for
!!   initialization
      implicit none
      include 'SIZE'
      include 'TRIP'
      integer il
      integer iff(tripl_nline_max), iy(tripl_nline_max)
      integer ir(97,tripl_nline_max)
      integer m,ia,ic,j
      real rm
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      save iff,ir,iy
      data iff /tripl_nline_max*0/
!-----------------------------------------------------------------------
      ! initialise
      if (tripl_seed(il).lt.0.or.iff(il).eq.0) then
         iff(il)=1
         tripl_seed(il)=mod(ic-tripl_seed(il),m)
         do j=1,97
            tripl_seed(il)=mod(ia*tripl_seed(il)+ic,m)
            ir(j,il)=tripl_seed(il)
         end do
         tripl_seed(il)=mod(ia*tripl_seed(il)+ic,m)
         iy(il)=tripl_seed(il)
      end if
      
      ! generate random number
      j=1+(97*iy(il))/m
      iy(il)=ir(j,il)
      tripl_ran2=iy(il)*rm
      tripl_seed(il)=mod(ia*tripl_seed(il)+ic,m)
      ir(j,il)=tripl_seed(il)

      end function
!=======================================================================
!> Generate forcing along 1D line
!! @ingroup tripl
!! @param[in] ifreset    reset flag
      subroutine tripl_frcs_get(ifreset)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'TRIP'

      ! argument list
      logical ifreset

#ifdef TRIPL_PR_RST
      ! variables necessary to reset pressure projection for P_n-P_n-2
      integer nprv(2)
      common /orthbi/ nprv

      ! variables necessary to reset velocity projection for P_n-P_n-2
      include 'VPROJ'
#endif      
      ! local variables
      integer il, jl, kl, ll
      integer istart
      real theta0, theta
      logical ifntdt_dif

#ifdef DEBUG
      character*3 str1, str2
      integer iunit, ierr
      ! call number
      integer icalldl
      save icalldl
      data icalldl /0/
#endif
!-----------------------------------------------------------------------
      ! reset all
      if (ifreset) then
         do il= 1, tripl_nline
            ! do we include time independent part?
            if (tripl_tiamp(il).gt.0.0) then
               istart = 1
            else
               istart = 2
            endif
            ! get forcing
            do jl = istart, tripl_nset_max
               call rzero(tripl_frcs(1,jl,il),tripl_npoint(il))
               do kl= 1, tripl_npoint(il)
                  theta0 = 2*pi*tripl_prj(kl,il)
                  do ll= 1, tripl_nmode(il)
                     theta = theta0*ll
                     tripl_frcs(kl,jl,il) = tripl_frcs(kl,jl,il) + sin(theta+tripl_rphs(ll,jl,il))
                  enddo
               enddo
            enddo
            ! rescale time independent part
            if (tripl_tiamp(il).gt.0.0) call cmult(tripl_frcs(1,1,il),tripl_tiamp(il),tripl_npoint(il))
         enddo

      else
         ! reset only time dependent part if needed
         ifntdt_dif = .FALSE.
         do il= 1, tripl_nline
            if (tripl_ntdt(il).ne.tripl_ntdt_old(il)) then
               ifntdt_dif = .TRUE.
               do jl= tripl_nset_max,3,-1
                  call copy(tripl_frcs(1,jl,il),tripl_frcs(1,jl-1,il),tripl_npoint(il))
               enddo
               call rzero(tripl_frcs(1,2,il),tripl_npoint(il))
               do jl= 1, tripl_npoint(il)
                  theta0 = 2*pi*tripl_prj(jl,il)
                  do kl= 1, tripl_nmode(il)
                     theta = theta0*kl
                     tripl_frcs(jl,2,il) = tripl_frcs(jl,2,il) + sin(theta+tripl_rphs(kl,2,il))
                  enddo
               enddo
            endif
         enddo
         if (ifntdt_dif) then
#ifdef TRIPL_PR_RST
            ! reset projection space
            ! pressure
            if (int(PARAM(95)).gt.0) then
               PARAM(95) = ISTEP
               nprv(1) = 0      ! veloctiy field only
            endif
            ! velocity
            if (int(PARAM(94)).gt.0) then
               PARAM(94) = ISTEP!+2
               ivproj(2,1) = 0
               ivproj(2,2) = 0
               if (IF3D) ivproj(2,3) = 0
            endif
#endif
         endif
      endif
      
      ! get tripping for current time step
      do il= 1, tripl_nline
         if (tripl_tiamp(il).gt.0.0) then
            call copy(tripl_ftrp(1,il),tripl_frcs(1,1,il),
     $           tripl_npoint(il))
         else
            call rzero(tripl_ftrp(1,il),tripl_npoint(il))
         end if
      end do
      ! interpolation in time
      do il = 1, tripl_nline
         theta0= time/tripl_tdt(il)-real(tripl_ntdt(il))
         if (theta0.gt.0.0) then
            theta0=theta0*theta0*(3.0-2.0*theta0)
            !theta0=theta0*theta0*theta0*(10.0+(6.0*theta0-15.0)*theta0)
            do jl= 1, tripl_npoint(il)
               tripl_ftrp(jl,il) = tripl_ftrp(jl,il) +
     $              tripl_tdamp(il)*((1.0-theta0)*tripl_frcs(jl,3,il) +
     $              theta0*tripl_frcs(jl,2,il))
            enddo
         else
            theta0=theta0+1.0
            theta0=theta0*theta0*(3.0-2.0*theta0)
            !theta0=theta0*theta0*theta0*(10.0+(6.0*theta0-15.0)*theta0)
            do jl= 1, tripl_npoint(il)
               tripl_ftrp(jl,il) = tripl_ftrp(jl,il) +
     $              tripl_tdamp(il)*((1.0-theta0)*tripl_frcs(jl,4,il) +
     $              theta0*tripl_frcs(jl,3,il))
            enddo
         endif
      enddo

#ifdef DEBUG
      ! for testing, to output refinement
      icalldl = icalldl+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalldl
      open(unit=iunit,file='trp_fcr.txt'//str1//'i'//str2)
      do il=1,tripl_npoint(1)
         write(iunit,*) il,tripl_prj(il,1),tripl_ftrp(il,1),tripl_frcs(il,1:4,1)
      enddo
      close(iunit)
#endif
      
      return
      end subroutine
!=======================================================================