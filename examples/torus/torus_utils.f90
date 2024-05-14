      subroutine interpolate_spectral(th_sec)

         implicit none
         include 'SIZE'
         include 'TOTAL'
         include 'HELIXD'
         real, intent(in) :: th_sec
         real vel_s(lx1,ly1,lz1,lelv)
         real vel_rho(lx1,ly1,lz1,lelv)
         real vel_th(lx1,ly1,lz1,lelv)
         real pres(lx1,ly1,lz1,lelv)
         real tmp(lx1*ly1*lz1,lelt,2)

         call compute_usrt(vel_s,vel_rho,vel_th)
         call interp_cross(vel_s,   th_sec, 'interp_data_us_', 12)
         call interp_cross(vel_rho, th_sec, 'interp_data_ur_', 13)
         call interp_cross(vel_th,  th_sec, 'interp_data_ut_', 14)
         call mappr(pres,PR,tmp(1,1,1),tmp(1,1,2))
         call interp_cross(pres,    th_sec, 'interp_data_pr_', 15)

      end subroutine interpolate_spectral

      subroutine save_ubar
         implicit none
         include 'SIZE'
         include 'TOTAL'
         include 'HELIXD'
         real current_ubar, get_ubar
         logical, save :: initialized
         data initialized/.false./

         if (.not. initialized) then
            if (nid == 0) write (6, *) 'Initializing save_ubar routine...'
            if (nid == 0) open (9833, file='ubar.dat', action='write', status='replace')
            initialized = .true.
         end if

         current_ubar = get_ubar(vx,vy,vz)

         if (nid == 0) then
            write(9833, "(2E15.7)") time, current_ubar
            write(6,*) 't =', time, 'ubar = ', current_ubar, 'Re = ', current_ubar/param(2)
         endif
      end subroutine save_ubar

      subroutine helix_pipe()
         implicit none

         include 'SIZE'            ! NX1, NY1, NZ1, NELV, NID
         include 'TSTEP'           ! TIME, DT
         include 'PARALLEL'        ! LGLEL
         include 'INPUT'           ! IF3D
         include 'SOLN'            ! VX, VY, VZ, VMULT
         include 'GEOM'            ! XM1, YM1, ZM1
         include 'HELIXD'

         real rtmp, xmin, xmax, scale, angle
         integer iel, il, jl, kl, nn
         real glmax, glmin

!     Geometry modification for helical pipe

         rad  = diam/2.
         radm = -rad
         pi = 4.*atan(1.)

         if (nid.eq.0) write(6,*) 'Start rescaling'
         call rescale_x(xm1,radm,rad)
         call rescale_x(ym1,radm,rad)
         call rescale_x(zm1,0.,1.) ! rescale the length of the pipe
         if (nid.eq.0) write(6,*) 'Done rescaling'

         if (nid.eq.0) write(6,*) 'USERDAT2: phi = ', an_phi*180./pi
         if (nid.eq.0) write(6,*) 'USERDAT2: curv_radius = ', curv_radius

!     rotate mesh to set the centre of the pipe along x-axis
         do iel = 1, nelt
            do kl = 1,lz1
               do jl = 1,ly1
                  do il = 1,lx1
                     rtmp = xm1(il,jl,kl,iel)
                     xm1(il,jl,kl,iel) = zm1(il,jl,kl,iel)
                     zm1(il,jl,kl,iel) = - rtmp
                  enddo
               enddo
            enddo
         enddo

!     rescale the mesh
         nn = lx1*ly1*lz1*nelt
         xmin = glmin(xm1,nn)
         xmax = glmax(xm1,nn)
         scale = sweep/(xmax-xmin)
         do il =1,nn
            xm1(il,1,1,1) = scale*xm1(il,1,1,1)
         enddo
         ! call outpost(vx,vy,vz,pr,t,'str')
         call copy(zax,zm1,nn)

!     new centre of the plane
         do iel = 1, nelt
            do kl = 1,lz1
               do jl = 1,ly1
                  do il = 1,lx1
                     angle = xm1(il,jl,kl,iel)
                     ox(il,jl,kl,iel) = curv_radius*sin(angle)
                     oy(il,jl,kl,iel) = curv_radius*cos(angle)
                     xm1(il,jl,kl,iel) = ox(il,jl,kl,iel) +
     $                    ym1(il,jl,kl,iel)*sin(angle)
                     ym1(il,jl,kl,iel) = oy(il,jl,kl,iel) +
     $                    ym1(il,jl,kl,iel)*cos(angle)
                     xax(il,jl,kl,iel) = xm1(il,jl,kl,iel)
                     yax(il,jl,kl,iel) = ym1(il,jl,kl,iel)
                     xm1(il,jl,kl,iel) = xm1(il,jl,kl,iel)
     $                             -zax(il,jl,kl,iel)*sin(an_phi)*cos(angle)
                     ym1(il,jl,kl,iel) = ym1(il,jl,kl,iel)
     $                             +zax(il,jl,kl,iel)*sin(an_phi)*sin(angle)
                     zm1(il,jl,kl,iel) = angle*pitch_s
     $                             +zax(il,jl,kl,iel)*cos(an_phi)
                  enddo
               enddo
            enddo
         enddo

         param(59) = 1.   !  All elements deformed

         return
      end subroutine

!=======================================================================

      subroutine check_periodic(norm_periodic)
         implicit none

         include 'SIZE'            ! LX1, LY1, LZ1, LELT, NFIELD
         include 'TSTEP'           ! ISTEP, LASTEP
         include 'SOLN'            ! V[XYZ], T
         include 'MASS'            ! BM1
         include 'HELIXD'          ! VX_TPER,VY_TPER,VZ_TPER,T_TPER

         real cnht_glsc2_wt, norm_periodic

         cnht_sv = 1.

         if (istep.eq.0) then
            call cnht_oprzero(VX_TPER,VY_TPER,VZ_TPER,T_TPER)
         endif

         if (istep.ge.0) then
            call cnht_opsub2(VX_TPER,VY_TPER,VZ_TPER,T_TPER,VX,VY,VZ,T)

            ! norm of the difference
            norm_periodic = cnht_glsc2_wt(VX_TPER,VY_TPER,VZ_TPER,T_TPER,
     $                          VX_TPER,VY_TPER,VZ_TPER,T_TPER,BM1)

            if (nid.eq.0) write(6,*) 'NORM_PER = ', istep, norm_periodic

            norm_periodic = sqrt(norm_periodic)

            ! updating v[xyz]old with v[xyz] at present iteration
            call cnht_opcopy(VX_TPER,VY_TPER,VZ_TPER,T_TPER,VX,VY,VZ,T)

         endif

         if (nid.eq.0) write(6,*) 'normper = ', ISTEP, TIME, norm_periodic

         return
      end

!=========================================================================

      subroutine makebf_str_pulse(ffx_p,ffy_p,ffz_p,ix,iy,iz,ieg)
         implicit none

         include 'SIZE'
         include 'TSTEP'
         include 'PARALLEL'
         include 'INPUT'
         include 'HELIXD'

         integer ix,iy,iz,iel,nbes,ierr,e,ieg
         real dpds, dp1ds, dpdt, dpds_steady
         real ffx_p, ffy_p, ffz_p
         real p, pp, alpha, rr, angle_t, delta, tau, fshape

         complex imag
         parameter (imag=(0.0,1.0))
         complex sqrtmi, wwi, b_pulse, j0_b, j1_b, jn_b(2)
         parameter (sqrtmi=(1.-imag)/sqrt(2.))
         complex binv_pulse

         e = gllel(ieg)
         pi = 4.0d0*atan(1.0d0)
         angle_t  = atan2(xax(ix,iy,iz,e),yax(ix,iy,iz,e))  ! measured clockwise from 12 noon
         pp=xax(ix,iy,iz,e)*xax(ix,iy,iz,e)+yax(ix,iy,iz,e)*yax(ix,iy,iz,e)
         p   = sqrt(pp)-curv_radius
         alpha = atan2(p,zax(ix,iy,iz,e))
         rr = p*p + zax(ix,iy,iz,e)*zax(ix,iy,iz,e)
         delta = curv_radius/(curv_radius**2+pitch_s**2)
         tau = pitch_s/(curv_radius**2+pitch_s**2)
         wwi   = sqrtmi*sqrt(1.)*womersley
         call cbesj(wwi, 0., 1, 2, jn_b, nbes, ierr)
         j0_b   = jn_b(1) ; j1_b   = jn_b(2)
         binv_pulse   = 1.0d0/(2.0d0*j1_b/(j0_b*wwi)-1.0d0)
         fshape = 1.0d0/sqrt((1+delta*sqrt(rr)*sin(alpha))**2)
         dp1ds = 2.0d0*real( (dp1dsr + dp1dsi*imag)*cexp(imag*omega*time) )
         dpds = ( dp0ds + dp1ds )*fshape/curv_radius

         ffx_p =  dpds*cos(an_phi)*cos(angle_t)
         ffy_p = -dpds*cos(an_phi)*sin(angle_t)
         ffz_p =  dpds*sin(an_phi)

         return
      end

!=========================================================================

      function get_ubar(u,v,w)
         implicit none

         include 'SIZE'
         include 'TOTAL'
         include 'HELIXD'

         real u(lx1,ly1,lz1,lelv),v(lx1,ly1,lz1,lelv),w(lx1,ly1,lz1,lelv)
         real num, den, s_L, angle_t, delta, tau, glsum
         real x, y, z, p, pp, alpha, rr, ri, us, usr, ubar, get_ubar
         integer i, j, k, e, eg, n

         n=nx1*ny1*nz1*nelv

         num = 0.0d0
         den = 0.0d0

         !if (nid.eq.0) write(6,*) 'ps = ', pitch_s
         delta = curv_radius/(curv_radius**2+pitch_s**2)
         !if (nid.eq.0) write(6,*) 'delta = ',delta
         tau = pitch_s/(curv_radius**2+pitch_s**2)
         !if (nid.eq.0) write(6,*) 'tau = ', tau

         do e = 1,nelv
            eg = lglel(e)
            do k = 1,lz1
               do j = 1,ly1
                  do i = 1,lx1
                     x  = xax(i,j,k,e)
                     y  = yax(i,j,k,e)
                     z  = zax(i,j,k,e)
                     pi = atan(1.0d0)*4.0d0
                     pp   = x*x + y*y
                     p   = sqrt(pp)-curv_radius
                     alpha = atan2(p,z)
                     rr = p*p + z*z
                     angle_t = atan2(x,y) ! clockwise from y axis
                     ri = 1./sqrt((1+delta*sqrt(rr)*sin(alpha))**2)

                     us = cos(an_phi)*(cos(angle_t)*u(i,j,k,e)
     $                    -sin(angle_t)*v(i,j,k,e))+sin(an_phi)*w(i,j,k,e)

                     usr = us*ri      ! Streamwise u/r
                     num = num + usr*bm1(i,j,k,e)
                     den = den +  ri*bm1(i,j,k,e)
                  enddo
               enddo
            enddo
         enddo
         num=glsum(num,1)
         den=glsum(den,1)
         s_L = sweep*sqrt(curv_radius**2+pitch_s**2)
         ubar = num/den  ! "1/r"-weighted volumetric average of streamwise velocity
         ! if (nid.eq.0) write(6,*) 'Num = ', num
         ! if (nid.eq.0) write(6,*) 'Den = ', den/s_L, den
         get_ubar = ubar
         return
      end

!=========================================================================

      subroutine compute_usrt(vel_s, vel_rho, vel_th)

         implicit none

         include 'SIZE'
         include 'TOTAL'             ! v[xyz]
         include 'HELIXD'            ! curv_radius, an_phi

         ! input parameters
         ! internal variables
         integer n, i, j, k, e
         real x,y,z
         real angle_t, tt, rr, r2, alpha
         real vel_p(lx1,ly1,lz1,lelv), vel_v(lx1,ly1,lz1,lelv)
         ! output variables
         real vel_s(lx1,ly1,lz1,lelv), vel_rho(lx1,ly1,lz1,lelv),
     $        vel_th(lx1,ly1,lz1,lelv)

         ! Compute velocity components in helical coordinates
         n = nx1*ny1*nz1*nelv

         do e=1,nelv
            do k=1,lz1
               do j=1,ly1
                  do i=1,lx1
                     x  = xax(i,j,k,e)
                     y  = yax(i,j,k,e)
                     z  = zax(i,j,k,e)

                     angle_t = atan2(x,y) ! clockwise from y axis
                     rr = sqrt(x**2+y**2)
                     r2 = (rr-curv_radius)**2 + z**2

                     vel_s(i,j,k,e) = cos(an_phi)*(cos(angle_t)*vx(i,j,k,e)
     $                    -sin(angle_t)*vy(i,j,k,e)) + sin(an_phi)*vz(i,j,k,e)
                     vel_p(i,j,k,e) = sin(angle_t)*vx(i,j,k,e)
     $                               +cos(angle_t)*vy(i,j,k,e)
                     vel_v(i,j,k,e) = sin(an_phi)*(-cos(angle_t)*vx(i,j,k,e)
     $                    -sin(angle_t)*vy(i,j,k,e)) + cos(an_phi)*vz(i,j,k,e)

                     alpha = atan2(z,sqrt(x*x + y*y)-curv_radius)

                     vel_rho(i,j,k,e) = + cos(alpha)*vel_p(i,j,k,e)
     $                                  + sin(alpha)*vel_v(i,j,k,e)
                     vel_th(i,j,k,e)  = + sin(alpha)*vel_p(i,j,k,e)
     $                                  - cos(alpha)*vel_v(i,j,k,e)
                  enddo
               enddo
            enddo
         enddo

         return
      end

!=========================================================================

      subroutine compute_area(vel_s,ubulk,xarea)

         implicit none

         include 'SIZE'
         include 'TOTAL'             ! area
         include 'HELIXD'            ! nslices, nElf

         ! input parameters
         real vel_s(lx1,ly1,lz1,lelv)
         ! internal variables
         integer n, i, j, e, eg, islicez
         real vel_p(lx1,ly1,lz1,lelv), vel_v(lx1,ly1,lz1,lelv)
         real ffone(lx1,ly1,lz1,lelv)
         real work(nslices)
         ! output variables
         real ubulk(nslices), xarea(nslices)

         ! initialisation
         n = nx1*ny1*nz1*nelv
         call rzero(ubulk,nslices)
         call rzero(xarea,nslices)
         call rone(ffone,n)

         do e=1,nelv
            eg = lglel(e)
            islicez = (eg-1)/nElf+1
            do j=1,ly1
               do i=1,lx1
                  ubulk(islicez) = ubulk(islicez) + vel_s(i,j,8,e)*area(i,j,5,e)
                  xarea(islicez) = xarea(islicez) + ffone(i,j,1,e)*area(i,j,5,e)
               enddo
            enddo
         enddo
         call gop(ubulk,work,'+  ',nslices)
         call gop(xarea,work,'+  ',nslices)

         return
      end

      subroutine compute_norm(norm2_steady)
         implicit none

         include 'SIZE'            ! LX1, LY1, LZ1, LELT, NFIELD
         include 'TSTEP'           ! ISTEP, LASTEP
         include 'SOLN'            ! V[XYZ], T
         include 'MASS'            ! BM1
         include 'ADJOINT'         ! IFADJ
         include 'HELIXD'

         real cnht_glsc2_wt, norm2_steady

         cnht_sv = 1.

         if (istep.eq.1) then
            call cnht_oprzero(VXOLD,VYOLD,VZOLD,TOLD)
         endif

         if (istep.ge.1) then
            if(.not.IFADJ) then    !non-linear
               call cnht_opsub2(VXOLD,VYOLD,VZOLD,TOLD,VX,VY,VZ,T)
            else                   !linear
               call cnht_opsub2(VXOLD,VYOLD,VZOLD,TOLD,VXP,VYP,VZP,TP)
            endif

            !     norm of the difference
            norm2_steady = cnht_glsc2_wt(VXOLD,VYOLD,VZOLD,TOLD,
     $                            VXOLD,VYOLD,VZOLD,TOLD,BM1)

            if (nid.eq.0) write(6,*) 'NORM2 = ', istep, dt, norm2_steady

            norm2_steady=sqrt(norm2_steady)

            !     normalization by DT (estimate of the norm of time derivative)
            norm2_steady = norm2_steady/DT

            !     updating v[xyz]old with v[xyz] at present iteration
            if(.not.IFADJ) then    !non-linear
               call cnht_opcopy(VXOLD,VYOLD,VZOLD,TOLD,VX,VY,VZ,T)
            else                   !linear
               call cnht_opcopy(VXOLD,VYOLD,VZOLD,TOLD,VXP,VYP,VZP,TP)
            endif
         endif

         if (nid.eq.0) write(6,*) 'norm = ', ISTEP, TIME, norm2_steady

         return
      end
      subroutine interp_cross(vel_interp,th_sec,interp_prefix,file_id)
         implicit none

         include 'SIZE'
         include 'TSTEP'
         include 'PARALLEL'
         include 'GEOM'
         include 'HELIXD'

         real vel_interp(lx1,ly1,lz1,lelv)
         real th_sec
         integer i, j, k, e, il, file_id

         integer nidd,npp,nekcomm,nekgroup,nekreal
         common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

         integer nt,nxf,nyf,nzf,nfail
         integer npt_max,npoints,n,npts,npts1
         real tol,bb_t,toldist

         integer M1, M2, M1_d
         parameter(M1=30,M2=70)

         real pts      (ldim,M1,M2),
     $        fieldout (M1,M2),
     $        dist     (lhis),
     $        rst      (lhis*ldim),
     $        rrv      (2*M1),
     $        rrm      (M1,M2),
     $        thetam   (M1,M2)
         integer rcode(lhis),elid(lhis),proc(lhis)
         real ex(M1,M2)
         integer inth_hpts

         integer icalld
         save    icalld
         data    icalld  /0/

         character*132 interp_fname,  strstep
         character*14 interp_prefix

         write(*,*) 'INTERPOLATION ROUTINE: '
         write(*,*) 'time step and time: ',istep,time
         write(*,*) 'number of elements on this rank: ',nelv
         write(*,*) 'total number of elements: ',nelgv
         write(*,*) 'Max number of inteprolation points per rank : ',lhis
         write(*,*) 'Dimension: ',ldim
         write(*,*) 'Processor: ',nid,' of ',np

         ! Preparing file name
         if (nid.eq.0) then
            write(strstep,'(I0.5)') istep
            print *, trim(strstep)
            interp_fname = trim(interp_prefix) // trim(strstep)
            print *, trim(interp_fname)
         endif

         write(*,*) 'Done prep fname: ',nid,' of ',np

         ! Build the polar mesh for interpolation
         nt=lx1*ly1*lz1*nelv
         pi = 4.*atan(1.)

         if (icalld.eq.0) then
            M1_d = 2*M1
            do il=1,M1_d
               rrv(il) = rad*sin(pi*(M1_d+1-2*il)/(2*M1_d-2))! GLC nodes (Weideman $ Reddy)
            end do

            ! just on the master node
            if (nid.eq.0) then
               npts = M1*M2
               do i=1,M1
                  do j=1,M2
                     rrm(i,j) = rrv(i)
                     thetam(i,j) = 2*pi*(j-1)/M2
                     pts(1,i,j) = rrm(i,j)*sin(thetam(i,j))*sin(th_sec) ! x-coordinate
     $                          + curv_radius*sin(th_sec)
                     pts(2,i,j) = rrm(i,j)*sin(thetam(i,j))*cos(th_sec) ! y-coordinate
     $                          + curv_radius*cos(th_sec)
                     pts(3,i,j) = rrm(i,j)*cos(thetam(i,j))

                     pts(1,i,j) = pts(1,i,j)
     $                          - pts(3,i,j)*sin(an_phi)*cos(th_sec)
                     pts(2,i,j) = pts(2,i,j)
     $                          - pts(3,i,j)*sin(an_phi)*sin(th_sec)
                     pts(3,i,j) = pts(3,i,j)*cos(an_phi)+th_sec*pitch_s
                  end do
               end do
            else
               npts = 0 ! zero on all other ranks
            end if

            ! compute total number of interpolation points
            npts1 = npts
            call igop(npts1,npoints,'+  ',1)
            write(*,*) 'Number of interpolation pts on this rank: ',npts
            write(*,*) 'total number of interpolation points: ',npoints

            if (lhis.lt.npts) then
               write(*,*) 'increase lhis to at least ',npts
               call exitt()
            end if

            ! set up interpolation handle
            toldist = 5e-6
            tol     = 1e-14
            npt_max = npoints
            nxf     = 2*lx1        ! fine mesh for bb-test
            nyf     = 2*ly1
            nzf     = 2*lz1
            bb_t    = 0.01         ! relative size to expand bounding boxes by

            call fgslib_findpts_setup(inth_hpts,nekcomm,np,ldim,
     $           xm1,ym1,zm1,lx1,ly1,lz1,
     $           nelt,nxf,nyf,nzf,bb_t,nt,nt,
     $           npt_max,tol)
            call fgslib_findpts(inth_hpts,rcode,1,
     $           proc,1,
     $           elid,1,
     $           rst,ldim,
     $           dist,1,
     $           pts(1,1,1),ldim,
     $           pts(2,1,1),ldim,
     $           pts(3,1,1),ldim,npts)

            ! error checks using return code rcode
            nfail = 0
            do i=1,npts
               if(rcode(i).eq.1) then
                  if(sqrt(dist(i)).gt.toldist) then
                     nfail = nfail + 1
                     IF (NFAIL.LE.5) WRITE(6,'(a,1p4e15.7)')
     $                    'ERROR: point outside the mesh xy[z]d^2:',
     $                    (pts(k,i,1),k=1,ldim),dist(i)
                  endif
               else if(rcode(i).eq.2) then
                  nfail = nfail + 1
                  if (nfail.le.5) write(6,'(a,1p3e15.7)')
     $                 ' ERROR: point not within mesh xy[z]: !',
     $                 (pts(k,i,1),k=1,ldim)
               endif
            enddo
            if (nfail.gt.0) then
               call exitt()
            end if

            icalld = 1
         end if


         ! doing the actual interpolation
         call fgslib_findpts_eval(inth_hpts,fieldout,1,
     $        rcode,1,
     $        proc,1,
     $        elid,1,
     $        rst,ldim,npts,
     $        vel_interp)

         ! printing on screen
         if (nid.eq.0) then
!         do j=1,M2
!            do i=1,M1
!               write(*,*) i,j,pts(1,i,j),pts(2,i,j),pts(3,i,j),
!     $              fieldout(i,j),nid
!            end do
!         end do

            ! file output
            open(file=trim(interp_fname),form='unformatted',unit=file_id)
            write(file_id) ldim,istep,time
            write(file_id) M1,M2
            do k=1,ldim
               write(file_id) ((pts(k,i,j),i=1,M1),j=1,M2)
            end do
            write(file_id) ((fieldout(i,j),i=1,M1),j=1,M2)
            close(file_id)
         end if

         return
      end
!> @file cnht_tools.f
!! @ingroup cnht
!! @brief Set of utilities related to conjugated heat transfer to build
!!    single scalar product for velocity nad temperature.
!! @author Clio Saglietti, Adam Peplinski
!! @date Mar 4, 2019
!=======================================================================
!> @brief Calcualte forcing ralted to conjugated heat transfer
!! @ingroup cnht
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    ieg             global element number
      subroutine cnht_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
         implicit none

         include 'SIZE'            !
         include 'INPUT'           ! IF3D, IFHEAT, CPFLD
         include 'PARALLEL'        ! GLLEL
         include 'TSTEP'           ! IFIELD
         include 'SOLN'            ! JP, T, TP
         include 'ADJOINT'         ! IFADJ, G_ADJ, DTD[XYZ]
         include 'HELIXD'           ! CHGR[XYZ]

!     argument list
         real ffx, ffy, ffz
         integer ix,iy,iz,ieg

!     local variables
         integer iel, ip
         real rtmp
!-----------------------------------------------------------------------
         if (IFHEAT) then
            iel=GLLEL(ieg)
            if (JP.eq.0) then
               rtmp = T(ix,iy,iz,iel,IFIELD)/CPFLD(1,2)
               ffx = ffx + cnht_gx*rtmp
               ffy = ffy + cnht_gy*rtmp
               if (IF3D) ffz = ffz + cnht_gz*rtmp
            else
               ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(iel-1)))
               if (.not.IFADJ) then
                  rtmp = TP(ip,IFIELD,JP)/CPFLD(1,2)
                  ffx = ffx + G_ADJ(1)*rtmp
                  ffy = ffy + G_ADJ(2)*rtmp
                  if (IF3D) ffz = ffz + G_ADJ(3)*rtmp
               else
                  ffx = ffx - DTDX(ip)*TP(ip,IFIELD,JP)
                  ffy = ffy - DTDY(ip)*TP(ip,IFIELD,JP)
                  if (IF3D) ffz = ffz - DTDZ(ip)*TP(ip,IFIELD,JP)
               end if
            end if
         endif

         return
      end subroutine
!=======================================================================
!> @brief Set cpfld coefficient for given type of simulation
!! @ingroup cnht
      subroutine cnht_cpfld_set()
         implicit none

         include 'SIZE'            !
         include 'INPUT'           ! CPFLD, PARAM
         include 'ADJOINT'         ! IFADJ
         include 'HELIXD'          ! cnht_Ra, cnht_Ra
!-----------------------------------------------------------------------
         if (IFHEAT) then
            if (IFADJ) then
               CPFLD(1,1)=cnht_Ra/sqrt(cnht_Ra)
               CPFLD(1,2)=1.0

               CPFLD(2,1)=1.0/sqrt(cnht_Ra)
               CPFLD(2,2)=1.0
            else
               CPFLD(1,1)=1.0/sqrt(cnht_Ra)
               CPFLD(1,2)=1.0/cnht_Ra

               CPFLD(2,1)=1.0/sqrt(cnht_Ra)
               CPFLD(2,2)=1.0
            endif
         else
            if (PARAM(2).lt.0.0) then
               CPFLD(1,1) = -1.0/PARAM(2)
            else
               CPFLD(1,1) = PARAM(2)
            endif

            if (PARAM(1).lt.0.0) then
               CPFLD(1,2) = -1.0/PARAM(1)
            else
               CPFLD(1,2) = PARAM(1)
            endif
         endif

         return
      end subroutine
!=======================================================================
!> @brief Zero velocity and temperature vectors
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
      subroutine cnht_oprzero (a1,a2,a3,a4)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1)

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call rzero(a1,ntotv)
            call rzero(a2,ntotv)
            if(IF3D) call rzero(a3,ntotv)
         endif
         if (IFHEAT) call rzero(a4,ntott)
         return
      end subroutine
!=======================================================================
!> @brief Copy vectors A=B (velocity and temperature)
!! @ingroup cnht
!! @param[out] a1, a2, a3    vlocity field 1
!! @param[out] a4            temperature field 1
!! @param[in]  b1, b2, b3    vlocity field 2
!! @param[in]  b4            temperature field 2
      subroutine cnht_opcopy (a1,a2,a3,a4,b1,b2,b3,b4)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call copy(a1,b1,ntotv)
            call copy(a2,b2,ntotv)
            if(IF3D) call copy(a3,b3,ntotv)
         endif
         if (IFHEAT) call copy(a4,b4,ntott)

         return
      end subroutine
!=======================================================================
!> @brief Add velocity and temperature vectors A = A+B
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
      subroutine cnht_opadd2 (a1,a2,a3,a4,b1,b2,b3,b4)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call add2(a1,b1,ntotv)
            call add2(a2,b2,ntotv)
            if(IF3D) call add2(a3,b3,ntotv)
         endif
         if (IFHEAT) call add2(a4,b4,ntott)

         return
      end subroutine
!=======================================================================
!> @brief Subtract vectors A = A-B (velocity and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
      subroutine cnht_opsub2 (a1,a2,a3,a4,b1,b2,b3,b4)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call sub2(a1,b1,ntotv)
            call sub2(a2,b2,ntotv)
            if(IF3D) call sub2(a3,b3,ntotv)
         endif
         if (IFHEAT) call sub2(a4,b4,ntott)

         return
      end subroutine
!=======================================================================
!> @brief Subtract vectors A = B-C (velocity and temperature)
!! @ingroup cnht
!! @param[out] a1, a2, a3    vlocity field 1
!! @param[out] a4            temperature field 1
!! @param[in]  b1, b2, b3    vlocity field 2
!! @param[in]  b4            temperature field 2
!! @param[in]  c1, c2, c3    vlocity field 3
!! @param[in]  c4            temperature field 3
      subroutine cnht_opsub3 (a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
         real c1(1),c2(1),c3(1),c4(1)

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call sub3(a1,b1,c1,ntotv)
            call sub3(a2,b2,c2,ntotv)
            if(IF3D) call sub3(a3,b3,c3,ntotv)
         endif
         if (IFHEAT) call sub3(a4,b4,c4,ntott)
         return
      end subroutine
!=======================================================================
!> @brief Multiply vector by constant A = c*A (single coeff. for velocity
!!    and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity fields
!! @param[inout] a4            temperature field
!! @param[in]    const         coefficient
      subroutine cnht_opcmult (a1,a2,a3,a4,const)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1)
         real const

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call cmult(a1,const,ntotv)
            call cmult(a2,const,ntotv)
            if(IF3D) call cmult(a3,const,ntotv)
         endif
         if (IFHEAT) call cmult(a4,const,ntott)
         return
      end subroutine
!=======================================================================
!> @brief Multiply vector by constant A = c*A with separate const. for
!!    velocity and temperature
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity fields
!! @param[inout] a4            temperature field
!! @param[in]    const1        velocity coefficient
!! @param[in]    const2        temperature coefficient
      subroutine cnht_opcmult2c (a1,a2,a3,a4,const1, const2)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1)
         real const1, const2

!     local variables
         integer ntotv, ntott
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            call cmult(a1,const1,ntotv)
            call cmult(a2,const1,ntotv)
            if(IF3D) call cmult(a3,const1,ntotv)
         endif
         if (IFHEAT) call cmult(a4,const2,ntott)
         return
      end subroutine
!=======================================================================
!> @brief  Vector summation with scaling A = A+c*B (velocity and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
!! @param[in]    coeff         scaling coefficient
      subroutine cnht_opadd2cm (a1,a2,a3,a4,b1,b2,b3,b4,coeff)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
         real coeff

!     local variables
         integer ntotv, ntott
         integer il
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            if (IF3D) then
               do il=1,ntotv
                  a1(il) = a1(il) + b1(il)*coeff
                  a2(il) = a2(il) + b2(il)*coeff
                  a3(il) = a3(il) + b3(il)*coeff
               enddo
            else
               do il=1,ntotv
                  a1(il) = a1(il) + b1(il)*coeff
                  a2(il) = a2(il) + b2(il)*coeff
               enddo
            endif
         endif
         if (IFHEAT) then
            do il=1,ntott
               a4(il) = a4(il) + b4(il)*coeff
            enddo
         endif
         return
      end subroutine
!=======================================================================
!> @brief  Vector subtraction with scaling A = A-c*B (velocity and temperature)
!! @ingroup cnht
!! @param[inout] a1, a2, a3    vlocity field 1
!! @param[inout] a4            temperature field 1
!! @param[in]    b1, b2, b3    vlocity field 2
!! @param[in]    b4            temperature field 2
!! @param[in]    coeff         scaling coefficient
      subroutine cnht_opsub2cm (a1,a2,a3,a4,b1,b2,b3,b4,coeff)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D

!     argument list
         real a1(1),a2(1),a3(1),a4(1),b1(1),b2(1),b3(1),b4(1)
         real coeff

!     local variables
         integer ntotv, ntott
         integer il
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         if (IFFLOW) then
            if (IF3D) then
               do il=1,ntotv
                  a1(il) = a1(il) - b1(il)*coeff
                  a2(il) = a2(il) - b2(il)*coeff
                  a3(il) = a3(il) - b3(il)*coeff
               enddo
            else
               do il=1,ntotv
                  a1(il) = a1(il) - b1(il)*coeff
                  a2(il) = a2(il) - b2(il)*coeff
               enddo
            endif
         endif
         if (IFHEAT) then
            do il=1,ntott
               a4(il) = a4(il) - b4(il)*coeff
            enddo
         endif
         return
      end subroutine
!=======================================================================
!> @brief Weigth velocity and temperature fields
!! @ingroup cnht
!! @param[inout] lvx, lvy, lvz    vlocity fields
!! @param[inout] lt               temperature field
!! @param[in]    coeff            scaling coefficient
      subroutine cnht_weight_fun (lvx,lvy,lvz,lt,coeff)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'MASS'            ! VOLVM1, VOLTM1
         include 'HELIXD'           ! cnht_sc,cnht_sv,cnht_st

         ! argument list
         real lvx(1),lvy(1),lvz(1),lt(1)
         real coeff

         ! local variables
         real f1, f2
!-----------------------------------------------------------------------
         f1=cnht_sv/VOLVM1/coeff
         f2=cnht_st*cnht_sc/VOLTM1/coeff

         !rescale
         call cnht_opcmult2c (lvx,lvy,lvz,lt,f1,f2)

         return
      end subroutine
!=======================================================================
!> @brief Global inner product of velocity and temperature fields
!! @ingroup cnht
!! @param[in] b1, b2, b3    vlocity field 1
!! @param[in] b4            temperature field 1
!! @param[in] x1, x2, x3    vlocity field 2
!! @param[in] x4            temperature field 2
!! @param[in] wt            mass matrix
!! @return cnht_glsc2_wt
      real function cnht_glsc2_wt (b1,b2,b3,b4,x1,x2,x3,x4,wt)
         implicit none

         include 'SIZE'            ! N[XYZ]1, NEL[VT]
         include 'INPUT'           ! IFFLOW, IFHEAT, IF3D
         include 'MASS'            ! VOLVM1, VOLTM1
         include 'HELIXD'           ! cnht_sc,cnht_sv,cnht_st

!     argument list
         real b1(1),b2(1),b3(1),b4(1),x1(1),x2(1),x3(1),x4(1),wt(1)

!     local variables
         integer ntotv, ntott
         real sum, f1, f2
         integer il
!     functions
         real glsum
!-----------------------------------------------------------------------
         ntotv = NX1*NY1*NZ1*NELV
         ntott = NX1*NY1*NZ1*NELT

         ! scaling factor velocity vs temperature
         ! veorsion for newton
         !  f1 = coeff_v
         !  f2 = coeff_T
         ! version for oic
         f1=cnht_sv/VOLVM1
         f2=cnht_st*cnht_sc/VOLTM1

         sum = 0.
         if (IFFLOW) then          !if vel
            if (IFHEAT) then       !if temp $ vel
               if (IF3D) then
                  do il=1,ntotv
                     sum = sum + wt(il)*(f1*(b1(il)*x1(il)+b2(il)*x2(il)
     $                    +b3(il)*x3(il))+f2*b4(il)*x4(il))
                  end do
               else
                  do il=1,ntotv
                     sum =sum + wt(il)*(f1*(b1(il)*x1(il)+b2(il)*x2(il))
     $                    +f2*b4(il)*x4(il))
                  end do
               end if

               ! for conjugate heat transfer
               if (ntott.gt.ntotv) then
                  do il=ntotv+1,ntott
                     sum = sum + wt(il)*f2*b4(il)*x4(il)
                  end do
               end if
            else                   !just vel
               if (IF3D) then
                  do il=1,ntotv
                     sum = sum + wt(il)*f1*(b1(il)*x1(il)+
     $                    b2(il)*x2(il)+b3(il)*x3(il))
                  end do
               else
                  do il=1,ntotv
                     sum = sum + wt(il)*f1*(b1(il)*x1(il)+b2(il)*x2(il))
                  end do
               end if
            end if
         else                      !just temp
            if (IFHEAT) then
               do il=1,ntott
                  sum = sum + wt(il)*(f2*b4(il)*x4(il))
               end do
            end if
         end if

         cnht_glsc2_wt = glsum(sum,1)

         return
      end function
!=======================================================================
      SUBROUTINE CACAI(Z, FNU, KODE, MR, N, Y, NZ, RL, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CACAI
C***REFER TO  CAIRY
C
CC CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO CACON CAN RESULT IF CACON
C     IS CALLED FROM CAIRY.
C
C***ROUTINES CALLED  CASYI,CBKNU,CMLRI,CSERI,CS1S2,R1MACH
C***END PROLOGUE  CACAI
         COMPLEX CSGN, CSPN, C1, C2, Y, Z, ZN, CY
         REAL ALIM, ARG, ASCLE, AZ, CPN, DFNU, ELIM, FMR, FNU, PI, RL,
     *    SGN, SPN, TOL, YY, R1MACH
         INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
         DIMENSION Y(N), CY(2)
         DATA PI / 3.14159265358979324E0 /
         NZ = 0
         ZN = -Z
         AZ = CABS(Z)
         NN = N
         DFNU = FNU + FLOAT(N-1)
         IF (AZ.LE.2.0E0) GO TO 10
         IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
         CALL CSERI(ZN, FNU, KODE, NN, Y, NW, TOL, ELIM, ALIM)
         GO TO 40
   20    CONTINUE
         IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
         CALL CASYI(ZN, FNU, KODE, NN, Y, NW, RL, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 70
         GO TO 40
   30    CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
         CALL CMLRI(ZN, FNU, KODE, NN, Y, NW, TOL)
         IF(NW.LT.0) GO TO 70
   40    CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
         CALL CBKNU(ZN, FNU, KODE, 1, CY, NW, TOL, ELIM, ALIM)
         IF (NW.NE.0) GO TO 70
         FMR = FLOAT(MR)
         SGN = -SIGN(PI,FMR)
         CSGN = CMPLX(0.0E0,SGN)
         IF (KODE.EQ.1) GO TO 50
         YY = -AIMAG(ZN)
         CPN = COS(YY)
         SPN = SIN(YY)
         CSGN = CSGN*CMPLX(CPN,SPN)
   50    CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(FNU)
         ARG = (FNU-FLOAT(INU))*SGN
         CPN = COS(ARG)
         SPN = SIN(ARG)
         CSPN = CMPLX(CPN,SPN)
         IF (MOD(INU,2).EQ.1) CSPN = -CSPN
         C1 = CY(1)
         C2 = Y(1)
         IF (KODE.EQ.1) GO TO 60
         IUF = 0
         ASCLE = 1.0E+3*R1MACH(1)/TOL
         CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
         NZ = NZ + NW
   60    CONTINUE
         Y(1) = CSPN*C1 + CSGN*C2
         RETURN
   70    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      SUBROUTINE CACON(Z, FNU, KODE, MR, N, Y, NZ, RL, FNUL, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CACON
C***REFER TO  CBESK,CBESH
C
CC CACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***ROUTINES CALLED  CBINU,CBKNU,CS1S2,R1MACH
C***END PROLOGUE  CACON
         COMPLEX CK, CONE, CS, CSCL, CSCR, CSGN, CSPN, CSS, CSR, C1, C2,
     *    RZ, SC1, SC2, ST, S1, S2, Y, Z, ZN, CY
         REAL ALIM, ARG, ASCLE, AS2, BSCLE, BRY, CPN, C1I, C1M, C1R, ELIM,
     *    FMR, FNU, FNUL, PI, RL, SGN, SPN, TOL, YY, R1MACH
         INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
         DIMENSION Y(N), CY(2), CSS(3), CSR(3), BRY(3)
         DATA PI / 3.14159265358979324E0 /
         DATA CONE / (1.0E0,0.0E0) /
         NZ = 0
         ZN = -Z
         NN = N
         CALL CBINU(ZN, FNU, KODE, NN, Y, NW, RL, FNUL, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 80
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
         NN = MIN0(2,N)
         CALL CBKNU(ZN, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
         IF (NW.NE.0) GO TO 80
         S1 = CY(1)
         FMR = FLOAT(MR)
         SGN = -SIGN(PI,FMR)
         CSGN = CMPLX(0.0E0,SGN)
         IF (KODE.EQ.1) GO TO 10
         YY = -AIMAG(ZN)
         CPN = COS(YY)
         SPN = SIN(YY)
         CSGN = CSGN*CMPLX(CPN,SPN)
   10    CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(FNU)
         ARG = (FNU-FLOAT(INU))*SGN
         CPN = COS(ARG)
         SPN = SIN(ARG)
         CSPN = CMPLX(CPN,SPN)
         IF (MOD(INU,2).EQ.1) CSPN = -CSPN
         IUF = 0
         C1 = S1
         C2 = Y(1)
         ASCLE = 1.0E+3*R1MACH(1)/TOL
         IF (KODE.EQ.1) GO TO 20
         CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
         NZ = NZ + NW
         SC1 = C1
   20    CONTINUE
         Y(1) = CSPN*C1 + CSGN*C2
         IF (N.EQ.1) RETURN
         CSPN = -CSPN
         S2 = CY(2)
         C1 = S2
         C2 = Y(2)
         IF (KODE.EQ.1) GO TO 30
         CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
         NZ = NZ + NW
         SC2 = C1
   30    CONTINUE
         Y(2) = CSPN*C1 + CSGN*C2
         IF (N.EQ.2) RETURN
         CSPN = -CSPN
         RZ = CMPLX(2.0E0,0.0E0)/ZN
         CK = CMPLX(FNU+1.0E0,0.0E0)*RZ
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         CSCR = CMPLX(TOL,0.0E0)
         CSS(1) = CSCL
         CSS(2) = CONE
         CSS(3) = CSCR
         CSR(1) = CSCR
         CSR(2) = CONE
         CSR(3) = CSCL
         BRY(1) = ASCLE
         BRY(2) = 1.0E0/ASCLE
         BRY(3) = R1MACH(2)
         AS2 = CABS(S2)
         KFLAG = 2
         IF (AS2.GT.BRY(1)) GO TO 40
         KFLAG = 1
         GO TO 50
   40    CONTINUE
         IF (AS2.LT.BRY(2)) GO TO 50
         KFLAG = 3
   50    CONTINUE
         BSCLE = BRY(KFLAG)
         S1 = S1*CSS(KFLAG)
         S2 = S2*CSS(KFLAG)
         CS = CSR(KFLAG)
         DO 70 I=3,N
            ST = S2
            S2 = CK*S2 + S1
            S1 = ST
            C1 = S2*CS
            ST = C1
            C2 = Y(I)
            IF (KODE.EQ.1) GO TO 60
            IF (IUF.LT.0) GO TO 60
            CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
            NZ = NZ + NW
            SC1 = SC2
            SC2 = C1
            IF (IUF.NE.3) GO TO 60
            IUF = -4
            S1 = SC1*CSS(KFLAG)
            S2 = SC2*CSS(KFLAG)
            ST = SC2
   60       CONTINUE
            Y(I) = CSPN*C1 + CSGN*C2
            CK = CK + RZ
            CSPN = -CSPN
            IF (KFLAG.GE.3) GO TO 70
            C1R = REAL(C1)
            C1I = AIMAG(C1)
            C1R = ABS(C1R)
            C1I = ABS(C1I)
            C1M = AMAX1(C1R,C1I)
            IF (C1M.LE.BSCLE) GO TO 70
            KFLAG = KFLAG + 1
            BSCLE = BRY(KFLAG)
            S1 = S1*CS
            S2 = ST
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            CS = CSR(KFLAG)
   70    CONTINUE
         RETURN
   80    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      SUBROUTINE CAIRY(Z, ID, KODE, AI, NZ, IERR)

c*********************************************************************72
c
cc CAIRY computes the complex Airy function AI(Z) or its derivative.
c
C***BEGIN PROLOGUE  CAIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C         ON KODE=1, CAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
C         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
C         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
C         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
C         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z)
C
C         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
C         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
C         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
C         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             AI=AI(Z)                ON ID=0 OR
C                             AI=DAI(Z)/DZ            ON ID=1
C                        = 2  RETURNS
C                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
C                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)
C
C         OUTPUT
C           AI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           NZ     - UNDERFLOW INDICATOR
C                    NZ= 0   , NORMAL RETURN
C                    NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
C                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
C                            TOO LARGE WITH KODE=1.
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C
C***LONG DESCRIPTION
C
C         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
C         FUNCTIONS BY
C
C            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
C                           C=1.0/(PI*SQRT(3.0))
C                           ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CACAI,CBKNU,I1MACH,R1MACH
C***END PROLOGUE  CAIRY
         COMPLEX AI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
         REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CK, COEF, C1, C2, DIG,
     *    DK, D1, D2, ELIM, FID, FNU, RL, R1M5, SFAC, TOL, TTH, ZI, ZR,
     *    Z3I, Z3R, R1MACH, BB, ALAZ
         INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
         DIMENSION CY(1)
         DATA TTH, C1, C2, COEF /6.66666666666666667E-01,
     *    3.55028053887817240E-01,2.58819403792806799E-01,
     *    1.83776298473930683E-01/
         DATA  CONE / (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CAIRY
         IERR = 0
         NZ=0
         IF (ID.LT.0 .OR. ID.GT.1) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (IERR.NE.0) RETURN
         AZ = CABS(Z)
         TOL = AMAX1(R1MACH(4),1.0E-18)
         FID = FLOAT(ID)
         IF (AZ.GT.1.0E0) GO TO 60
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
         S1 = CONE
         S2 = CONE
         IF (AZ.LT.TOL) GO TO 160
         AA = AZ*AZ
         IF (AA.LT.TOL/AZ) GO TO 40
         TRM1 = CONE
         TRM2 = CONE
         ATRM = 1.0E0
         Z3 = Z*Z*Z
         AZ3 = AZ*AA
         AK = 2.0E0 + FID
         BK = 3.0E0 - FID - FID
         CK = 4.0E0 - FID
         DK = 3.0E0 + FID + FID
         D1 = AK*DK
         D2 = BK*CK
         AD = AMIN1(D1,D2)
         AK = 24.0E0 + 9.0E0*FID
         BK = 30.0E0 - 9.0E0*FID
         Z3R = REAL(Z3)
         Z3I = AIMAG(Z3)
         DO 30 K=1,25
            TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
            S1 = S1 + TRM1
            TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
            S2 = S2 + TRM2
            ATRM = ATRM*AZ3/AD
            D1 = D1 + AK
            D2 = D2 + BK
            AD = AMIN1(D1,D2)
            IF (ATRM.LT.TOL*AD) GO TO 40
            AK = AK + 18.0E0
            BK = BK + 18.0E0
   30    CONTINUE
   40    CONTINUE
         IF (ID.EQ.1) GO TO 50
         AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
         IF (KODE.EQ.1) RETURN
         ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
         AI = AI*CEXP(ZTA)
         RETURN
   50    CONTINUE
         AI = -S2*CMPLX(C2,0.0E0)
         IF (AZ.GT.TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
         IF (KODE.EQ.1) RETURN
         ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
         AI = AI*CEXP(ZTA)
         RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   60    CONTINUE
         FNU = (1.0E0+FID)/3.0E0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         R1M5 = R1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         K1 = I1MACH(11) - 1
         AA = R1M5*FLOAT(K1)
         DIG = AMIN1(AA,18.0E0)
         AA = AA*2.303E0
         ALIM = ELIM + AMAX1(-AA,-41.45E0)
         RL = 1.2E0*DIG + 3.0E0
         ALAZ=ALOG(AZ)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA=0.5E0/TOL
         BB=FLOAT(I1MACH(9))*0.5E0
         AA=AMIN1(AA,BB)
         AA=AA**TTH
         IF (AZ.GT.AA) GO TO 260
         AA=SQRT(AA)
         IF (AZ.GT.AA) IERR=3
         CSQ=CSQRT(Z)
         ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
         IFLAG = 0
         SFAC = 1.0E0
         ZI = AIMAG(Z)
         ZR = REAL(Z)
         AK = AIMAG(ZTA)
         IF (ZR.GE.0.0E0) GO TO 70
         BK = REAL(ZTA)
         CK = -ABS(BK)
         ZTA = CMPLX(CK,AK)
   70    CONTINUE
         IF (ZI.NE.0.0E0) GO TO 80
         IF (ZR.GT.0.0E0) GO TO 80
         ZTA = CMPLX(0.0E0,AK)
   80    CONTINUE
         AA = REAL(ZTA)
         IF (AA.GE.0.0E0 .AND. ZR.GT.0.0E0) GO TO 100
         IF (KODE.EQ.2) GO TO 90
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         IF (AA.GT.(-ALIM)) GO TO 90
         AA = -AA + 0.25E0*ALAZ
         IFLAG = 1
         SFAC = TOL
         IF (AA.GT.ELIM) GO TO 240
   90    CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
         MR = 1
         IF (ZI.LT.0.0E0) MR = -1
         CALL CACAI(ZTA, FNU, KODE, MR, 1, CY, NN, RL, TOL, ELIM, ALIM)
         IF (NN.LT.0) GO TO 250
         NZ = NZ + NN
         GO TO 120
  100    CONTINUE
         IF (KODE.EQ.2) GO TO 110
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
         IF (AA.LT.ALIM) GO TO 110
         AA = -AA - 0.25E0*ALAZ
         IFLAG = 2
         SFAC = 1.0E0/TOL
         IF (AA.LT.(-ELIM)) GO TO 180
  110    CONTINUE
         CALL CBKNU(ZTA, FNU, KODE, 1, CY, NZ, TOL, ELIM, ALIM)
  120    CONTINUE
         S1 = CY(1)*CMPLX(COEF,0.0E0)
         IF (IFLAG.NE.0) GO TO 140
         IF (ID.EQ.1) GO TO 130
         AI = CSQ*S1
         RETURN
  130    AI = -Z*S1
         RETURN
  140    CONTINUE
         S1 = S1*CMPLX(SFAC,0.0E0)
         IF (ID.EQ.1) GO TO 150
         S1 = S1*CSQ
         AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
         RETURN
  150    CONTINUE
         S1 = -S1*Z
         AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
         RETURN
  160    CONTINUE
         AA = 1.0E+3*R1MACH(1)
         S1 = CMPLX(0.0E0,0.0E0)
         IF (ID.EQ.1) GO TO 170
         IF (AZ.GT.AA) S1 = CMPLX(C2,0.0E0)*Z
         AI = CMPLX(C1,0.0E0) - S1
         RETURN
  170    CONTINUE
         AI = -CMPLX(C2,0.0E0)
         AA = SQRT(AA)
         IF (AZ.GT.AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
         AI = AI + S1*CMPLX(C1,0.0E0)
         RETURN
  180    CONTINUE
         NZ = 1
         AI = CMPLX(0.0E0,0.0E0)
         RETURN
  240    CONTINUE
         NZ = 0
         IERR=2
         RETURN
  250    CONTINUE
         IF(NN.EQ.(-1)) GO TO 240
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         IERR=4
         NZ=0
         RETURN
      END
      SUBROUTINE CASYI(Z, FNU, KODE, N, Y, NZ, RL, TOL, ELIM, ALIM)

c*********************************************************************72
c
cc CASYI computs the I Bessel function.
c
C***BEGIN PROLOGUE  CASYI
C***REFER TO  CBESI,CBESK
C
C     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***ROUTINES CALLED  R1MACH
C***END PROLOGUE  CASYI
         COMPLEX AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1, RZ, S2,
     *    Y, Z
         REAL AA, ACZ, AEZ, AK, ALIM, ARG, ARM, ATOL, AZ, BB, BK, DFNU,
     *    DNU2, ELIM, FDN, FNU, PI, RL, RTPI, RTR1, S, SGN, SQK, TOL, X,
     *    YY, R1MACH
         INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
         DIMENSION Y(N)
         DATA PI, RTPI  /3.14159265358979324E0 , 0.159154943091895336E0 /
         DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
         NZ = 0
         AZ = CABS(Z)
         X = REAL(Z)
         ARM = 1.0E+3*R1MACH(1)
         RTR1 = SQRT(ARM)
         IL = MIN0(2,N)
         DFNU = FNU + FLOAT(N-IL)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         AK1 = CMPLX(RTPI,0.0E0)/Z
         AK1 = CSQRT(AK1)
         CZ = Z
         IF (KODE.EQ.2) CZ = Z - CMPLX(X,0.0E0)
         ACZ = REAL(CZ)
         IF (ABS(ACZ).GT.ELIM) GO TO 80
         DNU2 = DFNU + DFNU
         KODED = 1
         IF ((ABS(ACZ).GT.ALIM) .AND. (N.GT.2)) GO TO 10
         KODED = 0
         AK1 = AK1*CEXP(CZ)
   10    CONTINUE
         FDN = 0.0E0
         IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
         EZ = Z*CMPLX(8.0E0,0.0E0)
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
         AEZ = 8.0E0*AZ
         S = TOL/AEZ
         JL = INT(RL+RL) + 2
         YY = AIMAG(Z)
         P1 = CZERO
         IF (YY.EQ.0.0E0) GO TO 20
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
         INU = INT(FNU)
         ARG = (FNU-FLOAT(INU))*PI
         INU = INU + N - IL
         AK = -SIN(ARG)
         BK = COS(ARG)
         IF (YY.LT.0.0E0) BK = -BK
         P1 = CMPLX(AK,BK)
         IF (MOD(INU,2).EQ.1) P1 = -P1
   20    CONTINUE
         DO 50 K=1,IL
            SQK = FDN - 1.0E0
            ATOL = S*ABS(SQK)
            SGN = 1.0E0
            CS1 = CONE
            CS2 = CONE
            CK = CONE
            AK = 0.0E0
            AA = 1.0E0
            BB = AEZ
            DK = EZ
            DO 30 J=1,JL
               CK = CK*CMPLX(SQK,0.0E0)/DK
               CS2 = CS2 + CK
               SGN = -SGN
               CS1 = CS1 + CK*CMPLX(SGN,0.0E0)
               DK = DK + EZ
               AA = AA*ABS(SQK)/BB
               BB = BB + AEZ
               AK = AK + 8.0E0
               SQK = SQK - AK
               IF (AA.LE.ATOL) GO TO 40
   30       CONTINUE
            GO TO 90
   40       CONTINUE
            S2 = CS1
            IF (X+X.LT.ELIM) S2 = S2 + P1*CS2*CEXP(-Z-Z)
            FDN = FDN + 8.0E0*DFNU + 4.0E0
            P1 = -P1
            M = N - IL + K
            Y(M) = S2*AK1
   50    CONTINUE
         IF (N.LE.2) RETURN
         NN = N
         K = NN - 2
         AK = FLOAT(K)
         RZ = (CONE+CONE)/Z
         IB = 3
         DO 60 I=IB,NN
            Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
            AK = AK - 1.0E0
            K = K - 1
   60    CONTINUE
         IF (KODED.EQ.0) RETURN
         CK = CEXP(CZ)
         DO 70 I=1,NN
            Y(I) = Y(I)*CK
   70    CONTINUE
         RETURN
   80    CONTINUE
         NZ = -1
         RETURN
   90    CONTINUE
         NZ=-2
         RETURN
      END
      SUBROUTINE CBESH(Z, FNU, KODE, M, N, CY, NZ, IERR)

c*********************************************************************72
c
cc CBESH computes a sequence of complex Hankel functions.
c
C***BEGIN PROLOGUE  CBESH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
C         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
C         Z.NE.CMPLX(0.0E0,0.0E0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI.
C         ON KODE=2, CBESH COMPUTES THE SCALED HANKEL FUNCTIONS
C
C         CY(I)=H(M,FNU+J-1,Z)*EXP(-MM*Z*I)       MM=3-2M,      I**2=-1.
C
C         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER
C         AND LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN
C         THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z),      J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
C                                  J=1,...,N  ,  I**2=-1
C           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(J)=H(M,FNU+J-1,Z)  OR
C                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
C                    DEPENDING ON KODE, I**2=-1.
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
C                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR
C                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY
C                              HALF PLANES, NZ STATES ONLY THE NUMBER
C                              OF UNDERFLOWS.
C           IERR    -ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 TOO
C                            LARGE OR CABS(Z) TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE RELATION
C
C         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
C             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
C
C         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
C         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED
C         TO THE LEFT HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
C         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
C         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
C         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
C         WHOLE Z PLANE FOR Z TO INFINITY.
C
C         FOR NEGATIVE ORDERS,THE FORMULAE
C
C               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
C               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
C                         I**2=-1
C
C         CAN BE USED.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
C***END PROLOGUE  CBESH
C
         COMPLEX CY, Z, ZN, ZT, CSGN
         REAL AA, ALIM, ALN, ARG, AZ, CPN, DIG, ELIM, FMM, FN, FNU, FNUL,
     *    HPI, RHPI, RL, R1M5, SGN, SPN, TOL, UFL, XN, XX, YN, YY, R1MACH,
     *    BB, ASCLE, RTOL, ATOL
         INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     *    MM, MR, N, NN, NUF, NW, NZ, I1MACH
         DIMENSION CY(N)
C
         DATA HPI /1.57079632679489662E0/
C
C***FIRST EXECUTABLE STATEMENT  CBESH
         NZ=0
         XX = REAL(Z)
         YY = AIMAG(Z)
         IERR = 0
         IF (XX.EQ.0.0E0 .AND. YY.EQ.0.0E0) IERR=1
         IF (FNU.LT.0.0E0) IERR=1
         IF (M.LT.1 .OR. M.GT.2) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
         TOL = AMAX1(R1MACH(4),1.0E-18)
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         R1M5 = R1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         K1 = I1MACH(11) - 1
         AA = R1M5*FLOAT(K1)
         DIG = AMIN1(AA,18.0E0)
         AA = AA*2.303E0
         ALIM = ELIM + AMAX1(-AA,-41.45E0)
         FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
         RL = 1.2E0*DIG + 3.0E0
         FN = FNU + FLOAT(NN-1)
         MM = 3 - M - M
         FMM = FLOAT(MM)
         ZN = Z*CMPLX(0.0E0,-FMM)
         XN = REAL(ZN)
         YN = AIMAG(ZN)
         AZ = CABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA = 0.5E0/TOL
         BB=FLOAT(I1MACH(9))*0.5E0
         AA=AMIN1(AA,BB)
         IF(AZ.GT.AA) GO TO 240
         IF(FN.GT.AA) GO TO 240
         AA=SQRT(AA)
         IF(AZ.GT.AA) IERR=3
         IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
         UFL = R1MACH(1)*1.0E+3
         IF (AZ.LT.UFL) GO TO 220
         IF (FNU.GT.FNUL) GO TO 90
         IF (FN.LE.1.0E0) GO TO 70
         IF (FN.GT.2.0E0) GO TO 60
         IF (AZ.GT.TOL) GO TO 70
         ARG = 0.5E0*AZ
         ALN = -FN*ALOG(ARG)
         IF (ALN.GT.ELIM) GO TO 220
         GO TO 70
   60    CONTINUE
         CALL CUOIK(ZN, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
         IF (NUF.LT.0) GO TO 220
         NZ = NZ + NUF
         NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
         IF (NN.EQ.0) GO TO 130
   70    CONTINUE
         IF ((XN.LT.0.0E0) .OR. (XN.EQ.0.0E0 .AND. YN.LT.0.0E0 .AND.
     *    M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
         CALL CBKNU(ZN, FNU, KODE, NN, CY, NZ, TOL, ELIM, ALIM)
         GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80    CONTINUE
         MR = -MM
         CALL CACON(ZN, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 230
         NZ=NW
         GO TO 110
   90    CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
         MR = 0
         IF ((XN.GE.0.0E0) .AND. (XN.NE.0.0E0 .OR. YN.GE.0.0E0 .OR.
     *    M.NE.2)) GO TO 100
         MR = -MM
         IF (XN.EQ.0.0E0 .AND. YN.LT.0.0E0) ZN = -ZN
  100    CONTINUE
         CALL CBUNK(ZN, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 230
         NZ = NZ + NW
  110    CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
         SGN = SIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(FNU)
         INUH = INU/2
         IR = INU - 2*INUH
         ARG = (FNU-FLOAT(INU-IR))*SGN
         RHPI = 1.0E0/SGN
         CPN = RHPI*COS(ARG)
         SPN = RHPI*SIN(ARG)
C     ZN = CMPLX(-SPN,CPN)
         CSGN = CMPLX(-SPN,CPN)
C     IF (MOD(INUH,2).EQ.1) ZN = -ZN
         IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
         ZT = CMPLX(0.0E0,-FMM)
         RTOL = 1.0E0/TOL
         ASCLE = UFL*RTOL
         DO 120 I=1,NN
C       CY(I) = CY(I)*ZN
C       ZN = ZN*ZT
            ZN=CY(I)
            AA=REAL(ZN)
            BB=AIMAG(ZN)
            ATOL=1.0E0
            IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 125
            ZN = ZN*CMPLX(RTOL,0.0E0)
            ATOL = TOL
  125       CONTINUE
            ZN = ZN*CSGN
            CY(I) = ZN*CMPLX(ATOL,0.0E0)
            CSGN = CSGN*ZT
  120    CONTINUE
         RETURN
  130    CONTINUE
         IF (XN.LT.0.0E0) GO TO 220
         RETURN
  220    CONTINUE
         IERR=2
         NZ=0
         RETURN
  230    CONTINUE
         IF(NW.EQ.(-1)) GO TO 220
         NZ=0
         IERR=5
         RETURN
  240    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE CBESI(Z, FNU, KODE, N, CY, NZ, IERR)

c*********************************************************************72
c
cc CBESI computes a sequence of complex Bessel I functions.
c
C***BEGIN PROLOGUE  CBESI
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESI RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)
C
C         WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF.1)
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL I FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=I(FNU+J-1,Z), J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(J)=I(FNU+J-1,Z)  OR
C                    CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
C                    DEPENDING ON KODE, X=REAL(Z)
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0),
C                              J = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
C                            LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
C         SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
C         THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
C         NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
C         UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
C         FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
C         SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
C
C         THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
C         CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
C
C         I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z).GT.0.0
C                       M = +I OR -I,  I**2=-1
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
C         NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBINU,I1MACH,R1MACH
C***END PROLOGUE  CBESI
         COMPLEX CONE, CSGN, CY, Z, ZN
         REAL AA, ALIM, ARG, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, S1, S2,
     *    TOL, XX, YY, R1MACH, AZ, FN, BB, ASCLE, RTOL, ATOL
         INTEGER I, IERR, INU, K, KODE, K1, K2, N, NN, NZ, I1MACH
         DIMENSION CY(N)
         DATA PI /3.14159265358979324E0/
         DATA CONE / (1.0E0,0.0E0) /
C
C***FIRST EXECUTABLE STATEMENT  CBESI
         IERR = 0
         NZ=0
         IF (FNU.LT.0.0E0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         XX = REAL(Z)
         YY = AIMAG(Z)
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
         TOL = AMAX1(R1MACH(4),1.0E-18)
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         R1M5 = R1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         K1 = I1MACH(11) - 1
         AA = R1M5*FLOAT(K1)
         DIG = AMIN1(AA,18.0E0)
         AA = AA*2.303E0
         ALIM = ELIM + AMAX1(-AA,-41.45E0)
         RL = 1.2E0*DIG + 3.0E0
         FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
         AZ = CABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA = 0.5E0/TOL
         BB=FLOAT(I1MACH(9))*0.5E0
         AA=AMIN1(AA,BB)
         IF(AZ.GT.AA) GO TO 140
         FN=FNU+FLOAT(N-1)
         IF(FN.GT.AA) GO TO 140
         AA=SQRT(AA)
         IF(AZ.GT.AA) IERR=3
         IF(FN.GT.AA) IERR=3
         ZN = Z
         CSGN = CONE
         IF (XX.GE.0.0E0) GO TO 40
         ZN = -Z
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(FNU)
         ARG = (FNU-FLOAT(INU))*PI
         IF (YY.LT.0.0E0) ARG = -ARG
         S1 = COS(ARG)
         S2 = SIN(ARG)
         CSGN = CMPLX(S1,S2)
         IF (MOD(INU,2).EQ.1) CSGN = -CSGN
   40    CONTINUE
         CALL CBINU(ZN, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
         IF (NZ.LT.0) GO TO 120
         IF (XX.GE.0.0E0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
C-----------------------------------------------------------------------
         NN = N - NZ
         IF (NN.EQ.0) RETURN
         RTOL = 1.0E0/TOL
         ASCLE = R1MACH(1)*RTOL*1.0E+3
         DO 50 I=1,NN
C       CY(I) = CY(I)*CSGN
            ZN=CY(I)
            AA=REAL(ZN)
            BB=AIMAG(ZN)
            ATOL=1.0E0
            IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
            ZN = ZN*CMPLX(RTOL,0.0E0)
            ATOL = TOL
   55       CONTINUE
            ZN = ZN*CSGN
            CY(I) = ZN*CMPLX(ATOL,0.0E0)
            CSGN = -CSGN
   50    CONTINUE
         RETURN
  120    CONTINUE
         IF(NZ.EQ.(-2)) GO TO 130
         NZ = 0
         IERR=2
         RETURN
  130    CONTINUE
         NZ=0
         IERR=5
         RETURN
  140    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE CBESJ(Z, FNU, KODE, N, CY, NZ, IERR)

c*********************************************************************72
c
cc CBESJ computes a sequence of complex Bessel J functions.
c
C***BEGIN PROLOGUE  CBESJ
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESJ RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=J(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=J(FNU+I-1,Z)  OR
C                    CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE, Y=AIMAG(Z).
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0
C
C         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0
C
C         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
C         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBINU,I1MACH,R1MACH
C***END PROLOGUE  CBESJ
C
         COMPLEX CI, CSGN, CY, Z, ZN
         REAL AA, ALIM, ARG, DIG, ELIM, FNU, FNUL, HPI, RL, R1, R1M5, R2,
     *    TOL, YY, R1MACH, AZ, FN, BB, ASCLE, RTOL, ATOL
         INTEGER I, IERR, INU, INUH, IR, KODE, K1, K2, N, NL, NZ, I1MACH, K
         DIMENSION CY(N)
         DATA HPI /1.57079632679489662E0/
C
C***FIRST EXECUTABLE STATEMENT  CBESJ
         IERR = 0
         NZ=0
         IF (FNU.LT.0.0E0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
         TOL = AMAX1(R1MACH(4),1.0E-18)
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         R1M5 = R1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         K1 = I1MACH(11) - 1
         AA = R1M5*FLOAT(K1)
         DIG = AMIN1(AA,18.0E0)
         AA = AA*2.303E0
         ALIM = ELIM + AMAX1(-AA,-41.45E0)
         RL = 1.2E0*DIG + 3.0E0
         FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
         CI = CMPLX(0.0E0,1.0E0)
         YY = AIMAG(Z)
         AZ = CABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA = 0.5E0/TOL
         BB=FLOAT(I1MACH(9))*0.5E0
         AA=AMIN1(AA,BB)
         FN=FNU+FLOAT(N-1)
         IF(AZ.GT.AA) GO TO 140
         IF(FN.GT.AA) GO TO 140
         AA=SQRT(AA)
         IF(AZ.GT.AA) IERR=3
         IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(FNU)
         INUH = INU/2
         IR = INU - 2*INUH
         ARG = (FNU-FLOAT(INU-IR))*HPI
         R1 = COS(ARG)
         R2 = SIN(ARG)
         CSGN = CMPLX(R1,R2)
         IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
         ZN = -Z*CI
         IF (YY.GE.0.0E0) GO TO 40
         ZN = -ZN
         CSGN = CONJG(CSGN)
         CI = CONJG(CI)
   40    CONTINUE
         CALL CBINU(ZN, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
         IF (NZ.LT.0) GO TO 120
         NL = N - NZ
         IF (NL.EQ.0) RETURN
         RTOL = 1.0E0/TOL
         ASCLE = R1MACH(1)*RTOL*1.0E+3
         DO 50 I=1,NL
C       CY(I)=CY(I)*CSGN
            ZN=CY(I)
            AA=REAL(ZN)
            BB=AIMAG(ZN)
            ATOL=1.0E0
            IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
            ZN = ZN*CMPLX(RTOL,0.0E0)
            ATOL = TOL
   55       CONTINUE
            ZN = ZN*CSGN
            CY(I) = ZN*CMPLX(ATOL,0.0E0)
            CSGN = CSGN*CI
   50    CONTINUE
         RETURN
  120    CONTINUE
         IF(NZ.EQ.(-2)) GO TO 130
         NZ = 0
         IERR = 2
         RETURN
  130    CONTINUE
         NZ=0
         IERR=5
         RETURN
  140    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE CBESK(Z, FNU, KODE, N, CY, NZ, IERR)

c*********************************************************************72
c
cc CBESK computes a sequence of complex Bessel K functions.
c
C***BEGIN PROLOGUE  CBESK
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
C             BESSEL FUNCTION OF THE THIRD KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
C         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESK
C         RETURNS THE SCALED K FUNCTIONS,
C
C         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
C
C         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0E0
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=K(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
C                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C                    DEPENDING ON KODE
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
C                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
C                              IN THE SEQUENCE.
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
C         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
C         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
C         HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
C         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
C
C         FOR NEGATIVE ORDERS, THE FORMULA
C
C                       K(-FNU,Z) = K(FNU,Z)
C
C         CAN BE USED.
C
C         CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
C         AVAILABLE.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
C***END PROLOGUE  CBESK
C
         COMPLEX CY, Z
         REAL AA, ALIM, ALN, ARG, AZ, DIG, ELIM, FN, FNU, FNUL, RL, R1M5,
     *    TOL, UFL, XX, YY, R1MACH, BB
         INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, I1MACH
         DIMENSION CY(N)
C***FIRST EXECUTABLE STATEMENT  CBESK
         IERR = 0
         NZ=0
         XX = REAL(Z)
         YY = AIMAG(Z)
         IF (YY.EQ.0.0E0 .AND. XX.EQ.0.0E0) IERR=1
         IF (FNU.LT.0.0E0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
         TOL = AMAX1(R1MACH(4),1.0E-18)
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         R1M5 = R1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         K1 = I1MACH(11) - 1
         AA = R1M5*FLOAT(K1)
         DIG = AMIN1(AA,18.0E0)
         AA = AA*2.303E0
         ALIM = ELIM + AMAX1(-AA,-41.45E0)
         FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
         RL = 1.2E0*DIG + 3.0E0
         AZ = CABS(Z)
         FN = FNU + FLOAT(NN-1)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA = 0.5E0/TOL
         BB=FLOAT(I1MACH(9))*0.5E0
         AA=AMIN1(AA,BB)
         IF(AZ.GT.AA) GO TO 210
         IF(FN.GT.AA) GO TO 210
         AA=SQRT(AA)
         IF(AZ.GT.AA) IERR=3
         IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
C     UFL = EXP(-ELIM)
         UFL = R1MACH(1)*1.0E+3
         IF (AZ.LT.UFL) GO TO 180
         IF (FNU.GT.FNUL) GO TO 80
         IF (FN.LE.1.0E0) GO TO 60
         IF (FN.GT.2.0E0) GO TO 50
         IF (AZ.GT.TOL) GO TO 60
         ARG = 0.5E0*AZ
         ALN = -FN*ALOG(ARG)
         IF (ALN.GT.ELIM) GO TO 180
         GO TO 60
   50    CONTINUE
         CALL CUOIK(Z, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
         IF (NUF.LT.0) GO TO 180
         NZ = NZ + NUF
         NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
         IF (NN.EQ.0) GO TO 100
   60    CONTINUE
         IF (XX.LT.0.0E0) GO TO 70
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
C-----------------------------------------------------------------------
         CALL CBKNU(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 200
         NZ=NW
         RETURN
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
C-----------------------------------------------------------------------
   70    CONTINUE
         IF (NZ.NE.0) GO TO 180
         MR = 1
         IF (YY.LT.0.0E0) MR = -1
         CALL CACON(Z, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 200
         NZ=NW
         RETURN
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
   80    CONTINUE
         MR = 0
         IF (XX.GE.0.0E0) GO TO 90
         MR = 1
         IF (YY.LT.0.0E0) MR = -1
   90    CONTINUE
         CALL CBUNK(Z, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 200
         NZ = NZ + NW
         RETURN
  100    CONTINUE
         IF (XX.LT.0.0E0) GO TO 180
         RETURN
  180    CONTINUE
         NZ = 0
         IERR=2
         RETURN
  200    CONTINUE
         IF(NW.EQ.(-1)) GO TO 180
         NZ=0
         IERR=5
         RETURN
  210    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE CBESY(Z, FNU, KODE, N, CY, NZ, CWRK, IERR)

c*********************************************************************72
c
cc CBESY computes a sequence of complex Bessel Y functions.
c
C***BEGIN PROLOGUE  CBESY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101  (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF SECOND KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESY RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0E0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRK   - A COMPLEX WORK VECTOR OF DIMENSION AT LEAST N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON KODE=2)
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT IN TERMS OF THE I(FNU,Z) AND
C         K(FNU,Z) BESSEL FUNCTIONS IN THE RIGHT HALF PLANE BY
C
C             Y(FNU,Z) = I*CC*I(FNU,ARG) - (2/PI)*CONJG(CC)*K(FNU,ARG)
C
C             Y(FNU,Z) = CONJG(Y(FNU,CONJG(Z)))
C
C         FOR AIMAG(Z).GE.0 AND AIMAG(Z).LT.0 RESPECTIVELY, WHERE
C         CC=EXP(I*PI*FNU/2), ARG=Z*EXP(-I*PI/2) AND I**2=-1.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C             Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBESI,CBESK,I1MACH,R1MACH
C***END PROLOGUE  CBESY
C
         COMPLEX CWRK, CY, CI, CIP, CSGN, CSPN, EX, Z, ZU, ZV, ZZ, ZN
         REAL ARG, ELIM, EY, FNU, R1, R2, TAY, XX, YY, R1MACH, ASCLE, RTOL,
     *    ATOL, TOL, AA, BB, FFNU, HPI, RHPI, R1M5
         INTEGER I, IERR, IFNU, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH,
     *   I4
         DIMENSION CY(N), CWRK(N), CIP(4)
         DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     *    (1.0E0,0.0E0) , (0.0E0,1.0E0) , (-1.0E0,0.0E0) , (0.0E0,-1.0E0) /
         DATA HPI / 1.57079632679489662E0 /
C***FIRST EXECUTABLE STATEMENT  CBESY
         XX = REAL(Z)
         YY = AIMAG(Z)
         IERR = 0
         NZ=0
         IF (XX.EQ.0.0E0 .AND. YY.EQ.0.0E0) IERR=1
         IF (FNU.LT.0.0E0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         CI = CMPLX(0.0E0,1.0E0)
         ZZ=Z
         IF (YY.LT.0.0E0) ZZ=CONJG(Z)
         ZN = -CI*ZZ
         CALL CBESI(ZN, FNU, KODE, N, CY, NZ1, IERR)
         IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
         CALL CBESK(ZN, FNU, KODE, N, CWRK, NZ2, IERR)
         IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
         NZ = MIN0(NZ1,NZ2)
         IFNU = INT(FNU)
         FFNU = FNU - FLOAT(IFNU)
         ARG = HPI*FFNU
         CSGN = CMPLX(COS(ARG),SIN(ARG))
         I4 = MOD(IFNU,4) + 1
         CSGN = CSGN*CIP(I4)
         RHPI = 1.0E0/HPI
         CSPN = CONJG(CSGN)*CMPLX(RHPI,0.0E0)
         CSGN = CSGN*CI
         IF (KODE.EQ.2) GO TO 60
         DO 50 I=1,N
            CY(I) = CSGN*CY(I)-CSPN*CWRK(I)
            CSGN =  CI*CSGN
            CSPN = -CI*CSPN
   50    CONTINUE
         IF (YY.LT.0.0E0) THEN
            DO 55 I=1,N
               CY(I)=CONJG(CY(I))
   55       CONTINUE
         ENDIF
         RETURN
   60    CONTINUE
         R1 = COS(XX)
         R2 = SIN(XX)
         EX = CMPLX(R1,R2)
         TOL = AMAX1(R1MACH(4),1.0E-18)
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         K = MIN0(IABS(K1),IABS(K2))
         R1M5 = R1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         EY = 0.0E0
         TAY = ABS(YY+YY)
         IF (TAY.LT.ELIM) EY = EXP(-TAY)
         CSPN = EX*CMPLX(EY,0.0E0)*CSPN
         NZ = 0
         RTOL = 1.0E0/TOL
         ASCLE = R1MACH(1)*RTOL*1.0E+3
         DO 80 I=1,N
C----------------------------------------------------------------------
C       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN
C       SCALED MODE IF CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO
C       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.
C----------------------------------------------------------------------
            ZV = CWRK(I)
            AA=REAL(ZV)
            BB=AIMAG(ZV)
            ATOL=1.0E0
            IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 75
            ZV = ZV*CMPLX(RTOL,0.0E0)
            ATOL = TOL
   75       CONTINUE
            ZV = ZV*CSPN
            ZV = ZV*CMPLX(ATOL,0.0E0)
            ZU = CY(I)
            AA=REAL(ZU)
            BB=AIMAG(ZU)
            ATOL=1.0E0
            IF (AMAX1(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 85
            ZU = ZU*CMPLX(RTOL,0.0E0)
            ATOL = TOL
   85       CONTINUE
            ZU = ZU*CSGN
            ZU = ZU*CMPLX(ATOL,0.0E0)
            CY(I) = ZU - ZV
            IF (YY.LT.0.0E0) CY(I)=CONJG(CY(I))
            IF (CY(I).EQ.CMPLX(0.0E0,0.0E0) .AND. EY.EQ.0.0E0) NZ = NZ + 1
            CSGN =  CI*CSGN
            CSPN = -CI*CSPN
   80    CONTINUE
         RETURN
   90    CONTINUE
         NZ = 0
         RETURN
      END
      SUBROUTINE CBINU(Z, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CBINU
C***REFER TO  CBESH,CBESI,CBESJ,CBESK,CAIRY,CBIRY
C
cc CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE.
C
C***ROUTINES CALLED  CASYI,CBUNI,CMLRI,CSERI,CUOIK,CWRSK
C***END PROLOGUE  CBINU
         COMPLEX CW, CY, CZERO, Z
         REAL ALIM, AZ, DFNU, ELIM, FNU, FNUL, RL, TOL
         INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
         DIMENSION CY(N), CW(2)
         DATA CZERO / (0.0E0,0.0E0) /
C
         NZ = 0
         AZ = CABS(Z)
         NN = N
         DFNU = FNU + FLOAT(N-1)
         IF (AZ.LE.2.0E0) GO TO 10
         IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
         CALL CSERI(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
         INW = IABS(NW)
         NZ = NZ + INW
         NN = NN - INW
         IF (NN.EQ.0) RETURN
         IF (NW.GE.0) GO TO 120
         DFNU = FNU + FLOAT(NN-1)
   20    CONTINUE
         IF (AZ.LT.RL) GO TO 40
         IF (DFNU.LE.1.0E0) GO TO 30
         IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30    CONTINUE
         CALL CASYI(Z, FNU, KODE, NN, CY, NW, RL, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 130
         GO TO 120
   40    CONTINUE
         IF (DFNU.LE.1.0E0) GO TO 70
   50    CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
         CALL CUOIK(Z, FNU, KODE, 1, NN, CY, NW, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 130
         NZ = NZ + NW
         NN = NN - NW
         IF (NN.EQ.0) RETURN
         DFNU = FNU+FLOAT(NN-1)
         IF (DFNU.GT.FNUL) GO TO 110
         IF (AZ.GT.FNUL) GO TO 110
   60    CONTINUE
         IF (AZ.GT.RL) GO TO 80
   70    CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
         CALL CMLRI(Z, FNU, KODE, NN, CY, NW, TOL)
         IF(NW.LT.0) GO TO 130
         GO TO 120
   80    CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
         CALL CUOIK(Z, FNU, KODE, 2, 2, CW, NW, TOL, ELIM, ALIM)
         IF (NW.GE.0) GO TO 100
         NZ = NN
         DO 90 I=1,NN
            CY(I) = CZERO
   90    CONTINUE
         RETURN
  100    CONTINUE
         IF (NW.GT.0) GO TO 130
         CALL CWRSK(Z, FNU, KODE, NN, CY, NW, CW, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 130
         GO TO 120
  110    CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
         NUI = INT(FNUL-DFNU) + 1
         NUI = MAX0(NUI,0)
         CALL CBUNI(Z, FNU, KODE, NN, CY, NW, NUI, NLAST, FNUL, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 130
         NZ = NZ + NW
         IF (NLAST.EQ.0) GO TO 120
         NN = NLAST
         GO TO 60
  120    CONTINUE
         RETURN
  130    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      SUBROUTINE CBIRY(Z, ID, KODE, BI, IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CBIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C         ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
C         ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
C         DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
C         BOTH THE LEFT AND RIGHT HALF PLANES WHERE
C         ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
C         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             BI=BI(Z)                 ON ID=0 OR
C                             BI=DBI(Z)/DZ             ON ID=1
C                        = 2  RETURNS
C                             BI=CEXP(-AXZTA)*BI(Z)     ON ID=0 OR
C                             BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
C                             AND AXZTA=ABS(XZTA)
C
C         OUTPUT
C           BI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
C                            TOO LARGE WITH KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         BI AND DBI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE I BESSEL
C         FUNCTIONS BY
C
C                BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
C               DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
C                               C=1.0/SQRT(3.0)
C                               ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  CBINU,I1MACH,R1MACH
C***END PROLOGUE  CBIRY
         COMPLEX BI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
         REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BK, CK, COEF, C1, C2,
     *    DIG, DK, D1, D2, ELIM, FID, FMR, FNU, FNUL, PI, RL, R1M5, SFAC,
     *    TOL, TTH, ZI, ZR, Z3I, Z3R, R1MACH
         INTEGER ID, IERR, K, KODE, K1, K2, NZ, I1MACH
         DIMENSION CY(2)
         DATA TTH, C1, C2, COEF, PI /6.66666666666666667E-01,
     *    6.14926627446000736E-01,4.48288357353826359E-01,
     *    5.77350269189625765E-01,3.14159265358979324E+00/
         DATA  CONE / (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CBIRY
         IERR = 0
         NZ=0
         IF (ID.LT.0 .OR. ID.GT.1) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (IERR.NE.0) RETURN
         AZ = CABS(Z)
         TOL = AMAX1(R1MACH(4),1.0E-18)
         FID = FLOAT(ID)
         IF (AZ.GT.1.0E0) GO TO 60
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
         S1 = CONE
         S2 = CONE
         IF (AZ.LT.TOL) GO TO 110
         AA = AZ*AZ
         IF (AA.LT.TOL/AZ) GO TO 40
         TRM1 = CONE
         TRM2 = CONE
         ATRM = 1.0E0
         Z3 = Z*Z*Z
         AZ3 = AZ*AA
         AK = 2.0E0 + FID
         BK = 3.0E0 - FID - FID
         CK = 4.0E0 - FID
         DK = 3.0E0 + FID + FID
         D1 = AK*DK
         D2 = BK*CK
         AD = AMIN1(D1,D2)
         AK = 24.0E0 + 9.0E0*FID
         BK = 30.0E0 - 9.0E0*FID
         Z3R = REAL(Z3)
         Z3I = AIMAG(Z3)
         DO 30 K=1,25
            TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
            S1 = S1 + TRM1
            TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
            S2 = S2 + TRM2
            ATRM = ATRM*AZ3/AD
            D1 = D1 + AK
            D2 = D2 + BK
            AD = AMIN1(D1,D2)
            IF (ATRM.LT.TOL*AD) GO TO 40
            AK = AK + 18.0E0
            BK = BK + 18.0E0
   30    CONTINUE
   40    CONTINUE
         IF (ID.EQ.1) GO TO 50
         BI = S1*CMPLX(C1,0.0E0) + Z*S2*CMPLX(C2,0.0E0)
         IF (KODE.EQ.1) RETURN
         ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
         AA = REAL(ZTA)
         AA = -ABS(AA)
         BI = BI*CMPLX(EXP(AA),0.0E0)
         RETURN
   50    CONTINUE
         BI = S2*CMPLX(C2,0.0E0)
         IF (AZ.GT.TOL) BI = BI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
         IF (KODE.EQ.1) RETURN
         ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
         AA = REAL(ZTA)
         AA = -ABS(AA)
         BI = BI*CMPLX(EXP(AA),0.0E0)
         RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   60    CONTINUE
         FNU = (1.0E0+FID)/3.0E0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
         K1 = I1MACH(12)
         K2 = I1MACH(13)
         R1M5 = R1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
         K1 = I1MACH(11) - 1
         AA = R1M5*FLOAT(K1)
         DIG = AMIN1(AA,18.0E0)
         AA = AA*2.303E0
         ALIM = ELIM + AMAX1(-AA,-41.45E0)
         RL = 1.2E0*DIG + 3.0E0
         FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA=0.5E0/TOL
         BB=FLOAT(I1MACH(9))*0.5E0
         AA=AMIN1(AA,BB)
         AA=AA**TTH
         IF (AZ.GT.AA) GO TO 190
         AA=SQRT(AA)
         IF (AZ.GT.AA) IERR=3
         CSQ=CSQRT(Z)
         ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
         SFAC = 1.0E0
         ZI = AIMAG(Z)
         ZR = REAL(Z)
         AK = AIMAG(ZTA)
         IF (ZR.GE.0.0E0) GO TO 70
         BK = REAL(ZTA)
         CK = -ABS(BK)
         ZTA = CMPLX(CK,AK)
   70    CONTINUE
         IF (ZI.EQ.0.0E0 .AND. ZR.LE.0.0E0) ZTA = CMPLX(0.0E0,AK)
         AA = REAL(ZTA)
         IF (KODE.EQ.2) GO TO 80
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         BB = ABS(AA)
         IF (BB.LT.ALIM) GO TO 80
         BB = BB + 0.25E0*ALOG(AZ)
         SFAC = TOL
         IF (BB.GT.ELIM) GO TO 170
   80    CONTINUE
         FMR = 0.0E0
         IF (AA.GE.0.0E0 .AND. ZR.GT.0.0E0) GO TO 90
         FMR = PI
         IF (ZI.LT.0.0E0) FMR = -PI
         ZTA = -ZTA
   90    CONTINUE
C-----------------------------------------------------------------------
C     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
C     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBINU
C-----------------------------------------------------------------------
         CALL CBINU(ZTA, FNU, KODE, 1, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
         IF (NZ.LT.0) GO TO 180
         AA = FMR*FNU
         Z3 = CMPLX(SFAC,0.0E0)
         S1 = CY(1)*CMPLX(COS(AA),SIN(AA))*Z3
         FNU = (2.0E0-FID)/3.0E0
         CALL CBINU(ZTA, FNU, KODE, 2, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
         CY(1) = CY(1)*Z3
         CY(2) = CY(2)*Z3
C-----------------------------------------------------------------------
C     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
C-----------------------------------------------------------------------
         S2 = CY(1)*CMPLX(FNU+FNU,0.0E0)/ZTA + CY(2)
         AA = FMR*(FNU-1.0E0)
         S1 = (S1+S2*CMPLX(COS(AA),SIN(AA)))*CMPLX(COEF,0.0E0)
         IF (ID.EQ.1) GO TO 100
         S1 = CSQ*S1
         BI = S1*CMPLX(1.0E0/SFAC,0.0E0)
         RETURN
  100    CONTINUE
         S1 = Z*S1
         BI = S1*CMPLX(1.0E0/SFAC,0.0E0)
         RETURN
  110    CONTINUE
         AA = C1*(1.0E0-FID) + FID*C2
         BI = CMPLX(AA,0.0E0)
         RETURN
  170    CONTINUE
         NZ=0
         IERR=2
         RETURN
  180    CONTINUE
         IF(NZ.EQ.(-1)) GO TO 170
         NZ=0
         IERR=5
         RETURN
  190    CONTINUE
         IERR=4
         NZ=0
         RETURN
      END
      SUBROUTINE CBKNU(Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CBKNU
C***REFER TO  CBESI,CBESK,CAIRY,CBESH
C
C     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  CKSCL,CSHCH,GAMLN,I1MACH,R1MACH,CUCHK
C***END PROLOGUE  CBKNU
C
         COMPLEX CCH, CK, COEF, CONE, CRSC, CS, CSCL, CSH, CSR, CSS, CTWO,
     *    CZ, CZERO, F, FMU, P, PT, P1, P2, Q, RZ, SMU, ST, S1, S2, Y, Z,
     *    ZD, CELM, CY
         REAL AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, CC, DNU,
     *    DNU2, ELIM, ETEST, FC, FHS, FK, FKS, FNU, FPI, G1, G2, HPI, PI,
     *    P2I, P2M, P2R, RK, RTHPI, R1, S, SPI, TM, TOL, TTH, T1, T2, XX,
     *    YY, GAMLN, R1MACH, HELIM, ELM, XD, YD, ALAS, AS
         INTEGER I, IDUM, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N,
     *    NZ, I1MACH, NW, J, IC, INUB
         DIMENSION BRY(3), CC(8), CSS(3), CSR(3), Y(N), CY(2)
C
         DATA KMAX / 30 /
         DATA R1 / 2.0E0 /
         DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
C
         DATA PI, RTHPI, SPI ,HPI, FPI, TTH /
     1        3.14159265358979324E0,       1.25331413731550025E0,
     2        1.90985931710274403E0,       1.57079632679489662E0,
     3        1.89769999331517738E0,       6.66666666666666666E-01/
C
         DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1        5.77215664901532861E-01,    -4.20026350340952355E-02,
     2       -4.21977345555443367E-02,     7.21894324666309954E-03,
     3       -2.15241674114950973E-04,    -2.01348547807882387E-05,
     4        1.13302723198169588E-06,     6.11609510448141582E-09/
C
         XX = REAL(Z)
         YY = AIMAG(Z)
         CAZ = CABS(Z)
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         CRSC = CMPLX(TOL,0.0E0)
         CSS(1) = CSCL
         CSS(2) = CONE
         CSS(3) = CRSC
         CSR(1) = CRSC
         CSR(2) = CONE
         CSR(3) = CSCL
         BRY(1) = 1.0E+3*R1MACH(1)/TOL
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = R1MACH(2)
         NZ = 0
         IFLAG = 0
         KODED = KODE
         RZ = CTWO/Z
         INU = INT(FNU+0.5E0)
         DNU = FNU - FLOAT(INU)
         IF (ABS(DNU).EQ.0.5E0) GO TO 110
         DNU2 = 0.0E0
         IF (ABS(DNU).GT.TOL) DNU2 = DNU*DNU
         IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR CABS(Z).LE.R1
C-----------------------------------------------------------------------
         FC = 1.0E0
         SMU = CLOG(RZ)
         FMU = SMU*CMPLX(DNU,0.0E0)
         CALL CSHCH(FMU, CSH, CCH)
         IF (DNU.EQ.0.0E0) GO TO 10
         FC = DNU*PI
         FC = FC/SIN(FC)
         SMU = CSH*CMPLX(1.0E0/DNU,0.0E0)
   10    CONTINUE
         A2 = 1.0E0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
         T2 = EXP(-GAMLN(A2,IDUM))
         T1 = 1.0E0/(T2*FC)
         IF (ABS(DNU).GT.0.1E0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
         AK = 1.0E0
         S = CC(1)
         DO 20 K=2,8
            AK = AK*DNU2
            TM = CC(K)*AK
            S = S + TM
            IF (ABS(TM).LT.TOL) GO TO 30
   20    CONTINUE
   30    G1 = -S
         GO TO 50
   40    CONTINUE
         G1 = (T1-T2)/(DNU+DNU)
   50    CONTINUE
         G2 = 0.5E0*(T1+T2)*FC
         G1 = G1*FC
         F = CMPLX(G1,0.0E0)*CCH + SMU*CMPLX(G2,0.0E0)
         PT = CEXP(FMU)
         P = CMPLX(0.5E0/T2,0.0E0)*PT
         Q = CMPLX(0.5E0/T1,0.0E0)/PT
         S1 = F
         S2 = P
         AK = 1.0E0
         A1 = 1.0E0
         CK = CONE
         BK = 1.0E0 - DNU2
         IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
         IF (CAZ.LT.TOL) GO TO 70
         CZ = Z*Z*CMPLX(0.25E0,0.0E0)
         T1 = 0.25E0*CAZ*CAZ
   60    CONTINUE
         F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
         P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
         Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
         RK = 1.0E0/AK
         CK = CK*CZ*CMPLX(RK,0.0)
         S1 = S1 + CK*F
         A1 = A1*T1*RK
         BK = BK + AK + AK + 1.0E0
         AK = AK + 1.0E0
         IF (A1.GT.TOL) GO TO 60
   70    CONTINUE
         Y(1) = S1
         IF (KODED.EQ.1) RETURN
         Y(1) = S1*CEXP(Z)
         RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80    CONTINUE
         IF (CAZ.LT.TOL) GO TO 100
         CZ = Z*Z*CMPLX(0.25E0,0.0E0)
         T1 = 0.25E0*CAZ*CAZ
   90    CONTINUE
         F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
         P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
         Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
         RK = 1.0E0/AK
         CK = CK*CZ*CMPLX(RK,0.0E0)
         S1 = S1 + CK*F
         S2 = S2 + CK*(P-F*CMPLX(AK,0.0E0))
         A1 = A1*T1*RK
         BK = BK + AK + AK + 1.0E0
         AK = AK + 1.0E0
         IF (A1.GT.TOL) GO TO 90
  100    CONTINUE
         KFLAG = 2
         BK = REAL(SMU)
         A1 = FNU + 1.0E0
         AK = A1*ABS(BK)
         IF (AK.GT.ALIM) KFLAG = 3
         P2 = S2*CSS(KFLAG)
         S2 = P2*RZ
         S1 = S1*CSS(KFLAG)
         IF (KODED.EQ.1) GO TO 210
         F = CEXP(Z)
         S1 = S1*F
         S2 = S2*F
         GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110    CONTINUE
         COEF = CMPLX(RTHPI,0.0E0)/CSQRT(Z)
         KFLAG = 2
         IF (KODED.EQ.2) GO TO 120
         IF (XX.GT.ALIM) GO TO 290
C     BLANK LINE
         A1 = EXP(-XX)*REAL(CSS(KFLAG))
         PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
         COEF = COEF*PT
  120    CONTINUE
         IF (ABS(DNU).EQ.0.5E0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR CABS(Z).GT.R1
C-----------------------------------------------------------------------
         AK = COS(PI*DNU)
         AK = ABS(AK)
         IF (AK.EQ.0.0E0) GO TO 300
         FHS = ABS(0.25E0-DNU2)
         IF (FHS.EQ.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(11))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
         T1 = FLOAT(I1MACH(11)-1)*R1MACH(5)*3.321928094E0
         T1 = AMAX1(T1,12.0E0)
         T1 = AMIN1(T1,60.0E0)
         T2 = TTH*T1 - 6.0E0
         IF (XX.NE.0.0E0) GO TO 130
         T1 = HPI
         GO TO 140
  130    CONTINUE
         T1 = ATAN(YY/XX)
         T1 = ABS(T1)
  140    CONTINUE
         IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
C-----------------------------------------------------------------------
         ETEST = AK/(PI*CAZ*TOL)
         FK = 1.0E0
         IF (ETEST.LT.1.0E0) GO TO 180
         FKS = 2.0E0
         RK = CAZ + CAZ + 2.0E0
         A1 = 0.0E0
         A2 = 1.0E0
         DO 150 I=1,KMAX
            AK = FHS/FKS
            BK = RK/(FK+1.0E0)
            TM = A2
            A2 = BK*A2 - AK*A1
            A1 = TM
            RK = RK + 2.0E0
            FKS = FKS + FK + FK + 2.0E0
            FHS = FHS + FK + FK
            FK = FK + 1.0E0
            TM = ABS(A2)*FK
            IF (ETEST.LT.TM) GO TO 160
  150    CONTINUE
         GO TO 310
  160    CONTINUE
         FK = FK + SPI*T1*SQRT(T2/CAZ)
         FHS = ABS(0.25E0-DNU2)
         GO TO 180
  170    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
C-----------------------------------------------------------------------
         A2 = SQRT(CAZ)
         AK = FPI*AK/(TOL*SQRT(A2))
         AA = 3.0E0*T1/(1.0E0+CAZ)
         BB = 14.7E0*T1/(28.0E0+CAZ)
         AK = (ALOG(AK)+CAZ*COS(AA)/(1.0E0+0.008E0*CAZ))/COS(BB)
         FK = 0.12125E0*AK*AK/CAZ + 1.5E0
  180    CONTINUE
         K = INT(FK)
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
         FK = FLOAT(K)
         FKS = FK*FK
         P1 = CZERO
         P2 = CMPLX(TOL,0.0E0)
         CS = P2
         DO 190 I=1,K
            A1 = FKS - FK
            A2 = (FKS+FK)/(A1+FHS)
            RK = 2.0E0/(FK+1.0E0)
            T1 = (FK+XX)*RK
            T2 = YY*RK
            PT = P2
            P2 = (P2*CMPLX(T1,T2)-P1)*CMPLX(A2,0.0E0)
            P1 = PT
            CS = CS + P2
            FKS = A1 - FK + 1.0E0
            FK = FK - 1.0E0
  190    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
         TM = CABS(CS)
         PT = CMPLX(1.0E0/TM,0.0E0)
         S1 = PT*P2
         CS = CONJG(CS)*PT
         S1 = COEF*S1*CS
         IF (INU.GT.0 .OR. N.GT.1) GO TO 200
         ZD = Z
         IF(IFLAG.EQ.1) GO TO 270
         GO TO 240
  200    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
C-----------------------------------------------------------------------
         TM = CABS(P2)
         PT = CMPLX(1.0E0/TM,0.0E0)
         P1 = PT*P1
         P2 = CONJG(P2)*PT
         PT = P1*P2
         S2 = S1*(CONE+(CMPLX(DNU+0.5E0,0.0E0)-PT)/Z)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210    CONTINUE
         CK = CMPLX(DNU+1.0E0,0.0E0)*RZ
         IF (N.EQ.1) INU = INU - 1
         IF (INU.GT.0) GO TO 220
         IF (N.EQ.1) S1=S2
         ZD = Z
         IF(IFLAG.EQ.1) GO TO 270
         GO TO 240
  220    CONTINUE
         INUB = 1
         IF (IFLAG.EQ.1) GO TO 261
  225    CONTINUE
         P1 = CSR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 230 I=INUB,INU
            ST = S2
            S2 = CK*S2 + S1
            S1 = ST
            CK = CK + RZ
            IF (KFLAG.GE.3) GO TO 230
            P2 = S2*P1
            P2R = REAL(P2)
            P2I = AIMAG(P2)
            P2R = ABS(P2R)
            P2I = ABS(P2I)
            P2M = AMAX1(P2R,P2I)
            IF (P2M.LE.ASCLE) GO TO 230
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1 = S1*P1
            S2 = P2
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            P1 = CSR(KFLAG)
  230    CONTINUE
         IF (N.EQ.1) S1 = S2
  240    CONTINUE
         Y(1) = S1*CSR(KFLAG)
         IF (N.EQ.1) RETURN
         Y(2) = S2*CSR(KFLAG)
         IF (N.EQ.2) RETURN
         KK = 2
  250    CONTINUE
         KK = KK + 1
         IF (KK.GT.N) RETURN
         P1 = CSR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 260 I=KK,N
            P2 = S2
            S2 = CK*S2 + S1
            S1 = P2
            CK = CK + RZ
            P2 = S2*P1
            Y(I) = P2
            IF (KFLAG.GE.3) GO TO 260
            P2R = REAL(P2)
            P2I = AIMAG(P2)
            P2R = ABS(P2R)
            P2I = ABS(P2I)
            P2M = AMAX1(P2R,P2I)
            IF (P2M.LE.ASCLE) GO TO 260
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1 = S1*P1
            S2 = P2
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            P1 = CSR(KFLAG)
  260    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261    CONTINUE
         HELIM = 0.5E0*ELIM
         ELM = EXP(-ELIM)
         CELM = CMPLX(ELM,0.0)
         ASCLE = BRY(1)
         ZD = Z
         XD = XX
         YD = YY
         IC = -1
         J = 2
         DO 262 I=1,INU
            ST = S2
            S2 = CK*S2+S1
            S1 = ST
            CK = CK+RZ
            AS = CABS(S2)
            ALAS = ALOG(AS)
            P2R = -XD+ALAS
            IF(P2R.LT.(-ELIM)) GO TO 263
            P2 = -ZD+CLOG(S2)
            P2R = REAL(P2)
            P2I = AIMAG(P2)
            P2M = EXP(P2R)/TOL
            P1 = CMPLX(P2M,0.0E0)*CMPLX(COS(P2I),SIN(P2I))
            CALL CUCHK(P1,NW,ASCLE,TOL)
            IF(NW.NE.0) GO TO 263
            J=3-J
            CY(J) = P1
            IF(IC.EQ.(I-1)) GO TO 264
            IC = I
            GO TO 262
  263       CONTINUE
            IF(ALAS.LT.HELIM) GO TO 262
            XD = XD-ELIM
            S1 = S1*CELM
            S2 = S2*CELM
            ZD = CMPLX(XD,YD)
  262    CONTINUE
         IF(N.EQ.1) S1 = S2
         GO TO 270
  264    CONTINUE
         KFLAG = 1
         INUB = I+1
         S2 = CY(J)
         J = 3 - J
         S1 = CY(J)
         IF(INUB.LE.INU) GO TO 225
         IF(N.EQ.1) S1 = S2
         GO TO 240
  270    CONTINUE
         Y(1) = S1
         IF (N.EQ.1) GO TO 280
         Y(2) = S2
  280    CONTINUE
         ASCLE = BRY(1)
         CALL CKSCL(ZD, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
         INU = N - NZ
         IF (INU.LE.0) RETURN
         KK = NZ + 1
         S1 = Y(KK)
         Y(KK) = S1*CSR(1)
         IF (INU.EQ.1) RETURN
         KK = NZ + 2
         S2 = Y(KK)
         Y(KK) = S2*CSR(1)
         IF (INU.EQ.2) RETURN
         T2 = FNU + FLOAT(KK-1)
         CK = CMPLX(T2,0.0E0)*RZ
         KFLAG = 1
         GO TO 250
  290    CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY EXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
         KODED = 2
         IFLAG = 1
         KFLAG = 2
         GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300    CONTINUE
         S1 = COEF
         S2 = COEF
         GO TO 210
  310    CONTINUE
         NZ=-2
         RETURN
      END
      SUBROUTINE CBUNI(Z, FNU, KODE, N, Y, NZ, NUI, NLAST, FNUL, TOL,
     * ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CBUNI
C***REFER TO  CBESI,CBESK
C
C     CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***ROUTINES CALLED  CUNI1,CUNI2,R1MACH
C***END PROLOGUE  CBUNI
         COMPLEX CSCL, CSCR, CY, RZ, ST, S1, S2, Y, Z
         REAL ALIM, AX, AY, DFNU, ELIM, FNU, FNUI, FNUL, GNU, TOL, XX, YY,
     *    ASCLE, BRY, STR, STI, STM, R1MACH
         INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
         DIMENSION Y(N), CY(2), BRY(3)
         NZ = 0
         XX = REAL(Z)
         YY = AIMAG(Z)
         AX = ABS(XX)*1.7321E0
         AY = ABS(YY)
         IFORM = 1
         IF (AY.GT.AX) IFORM = 2
         IF (NUI.EQ.0) GO TO 60
         FNUI = FLOAT(NUI)
         DFNU = FNU + FLOAT(N-1)
         GNU = DFNU + FNUI
         IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
         CALL CUNI1(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
         GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
         CALL CUNI2(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   20    CONTINUE
         IF (NW.LT.0) GO TO 50
         IF (NW.NE.0) GO TO 90
         AY = CABS(CY(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
         BRY(1) = 1.0E+3*R1MACH(1)/TOL
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = BRY(2)
         IFLAG = 2
         ASCLE = BRY(2)
         AX = 1.0E0
         CSCL = CMPLX(AX,0.0E0)
         IF (AY.GT.BRY(1)) GO TO 21
         IFLAG = 1
         ASCLE = BRY(1)
         AX = 1.0E0/TOL
         CSCL = CMPLX(AX,0.0E0)
         GO TO 25
   21    CONTINUE
         IF (AY.LT.BRY(2)) GO TO 25
         IFLAG = 3
         ASCLE = BRY(3)
         AX = TOL
         CSCL = CMPLX(AX,0.0E0)
   25    CONTINUE
         AY = 1.0E0/AX
         CSCR = CMPLX(AY,0.0E0)
         S1 = CY(2)*CSCL
         S2 = CY(1)*CSCL
         RZ = CMPLX(2.0E0,0.0E0)/Z
         DO 30 I=1,NUI
            ST = S2
            S2 = CMPLX(DFNU+FNUI,0.0E0)*RZ*S2 + S1
            S1 = ST
            FNUI = FNUI - 1.0E0
            IF (IFLAG.GE.3) GO TO 30
            ST = S2*CSCR
            STR = REAL(ST)
            STI = AIMAG(ST)
            STR = ABS(STR)
            STI = ABS(STI)
            STM = AMAX1(STR,STI)
            IF (STM.LE.ASCLE) GO TO 30
            IFLAG = IFLAG+1
            ASCLE = BRY(IFLAG)
            S1 = S1*CSCR
            S2 = ST
            AX = AX*TOL
            AY = 1.0E0/AX
            CSCL = CMPLX(AX,0.0E0)
            CSCR = CMPLX(AY,0.0E0)
            S1 = S1*CSCL
            S2 = S2*CSCL
   30    CONTINUE
         Y(N) = S2*CSCR
         IF (N.EQ.1) RETURN
         NL = N - 1
         FNUI = FLOAT(NL)
         K = NL
         DO 40 I=1,NL
            ST = S2
            S2 = CMPLX(FNU+FNUI,0.0E0)*RZ*S2 + S1
            S1 = ST
            ST = S2*CSCR
            Y(K) = ST
            FNUI = FNUI - 1.0E0
            K = K - 1
            IF (IFLAG.GE.3) GO TO 40
            STR = REAL(ST)
            STI = AIMAG(ST)
            STR = ABS(STR)
            STI = ABS(STI)
            STM = AMAX1(STR,STI)
            IF (STM.LE.ASCLE) GO TO 40
            IFLAG = IFLAG+1
            ASCLE = BRY(IFLAG)
            S1 = S1*CSCR
            S2 = ST
            AX = AX*TOL
            AY = 1.0E0/AX
            CSCL = CMPLX(AX,0.0E0)
            CSCR = CMPLX(AY,0.0E0)
            S1 = S1*CSCL
            S2 = S2*CSCL
   40    CONTINUE
         RETURN
   50    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
   60    CONTINUE
         IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
         CALL CUNI1(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
         GO TO 80
   70    CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
         CALL CUNI2(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   80    CONTINUE
         IF (NW.LT.0) GO TO 50
         NZ = NW
         RETURN
   90    CONTINUE
         NLAST = N
         RETURN
      END
      SUBROUTINE CBUNK(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CBUNK
C***REFER TO  CBESK,CBESH
C
C     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
C
C***ROUTINES CALLED  CUNK1,CUNK2
C***END PROLOGUE  CBUNK
         COMPLEX Y, Z
         REAL ALIM, AX, AY, ELIM, FNU, TOL, XX, YY
         INTEGER KODE, MR, N, NZ
         DIMENSION Y(N)
         NZ = 0
         XX = REAL(Z)
         YY = AIMAG(Z)
         AX = ABS(XX)*1.7321E0
         AY = ABS(YY)
         IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
         CALL CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
         GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
         CALL CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
   20    CONTINUE
         RETURN
      END
      SUBROUTINE CKSCL(ZR, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CKSCL
C***REFER TO  CBKNU,CUNK1,CUNK2
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***ROUTINES CALLED  CUCHK
C***END PROLOGUE  CKSCL
         COMPLEX CK, CS, CY, CZERO, RZ, S1, S2, Y, ZR, ZD, CELM
         REAL AA, ASCLE, ACS, AS, CSI, CSR, ELIM, FN, FNU, TOL, XX, ZRI,
     *    ELM, ALAS, HELIM
         INTEGER I, IC, K, KK, N, NN, NW, NZ
         DIMENSION Y(N), CY(2)
         DATA CZERO / (0.0E0,0.0E0) /
C
         NZ = 0
         IC = 0
         XX = REAL(ZR)
         NN = MIN0(2,N)
         DO 10 I=1,NN
            S1 = Y(I)
            CY(I) = S1
            AS = CABS(S1)
            ACS = -XX + ALOG(AS)
            NZ = NZ + 1
            Y(I) = CZERO
            IF (ACS.LT.(-ELIM)) GO TO 10
            CS = -ZR + CLOG(S1)
            CSR = REAL(CS)
            CSI = AIMAG(CS)
            AA = EXP(CSR)/TOL
            CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
            CALL CUCHK(CS, NW, ASCLE, TOL)
            IF (NW.NE.0) GO TO 10
            Y(I) = CS
            NZ = NZ - 1
            IC = I
   10    CONTINUE
         IF (N.EQ.1) RETURN
         IF (IC.GT.1) GO TO 20
         Y(1) = CZERO
         NZ = 2
   20    CONTINUE
         IF (N.EQ.2) RETURN
         IF (NZ.EQ.0) RETURN
         FN = FNU + 1.0E0
         CK = CMPLX(FN,0.0E0)*RZ
         S1 = CY(1)
         S2 = CY(2)
         HELIM = 0.5E0*ELIM
         ELM = EXP(-ELIM)
         CELM = CMPLX(ELM,0.0E0)
         ZRI =AIMAG(ZR)
         ZD = ZR
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
         DO 30 I=3,N
            KK = I
            CS = S2
            S2 = CK*S2 + S1
            S1 = CS
            CK = CK + RZ
            AS = CABS(S2)
            ALAS = ALOG(AS)
            ACS = -XX + ALAS
            NZ = NZ + 1
            Y(I) = CZERO
            IF (ACS.LT.(-ELIM)) GO TO 25
            CS = -ZD + CLOG(S2)
            CSR = REAL(CS)
            CSI = AIMAG(CS)
            AA = EXP(CSR)/TOL
            CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
            CALL CUCHK(CS, NW, ASCLE, TOL)
            IF (NW.NE.0) GO TO 25
            Y(I) = CS
            NZ = NZ - 1
            IF (IC.EQ.(KK-1)) GO TO 40
            IC = KK
            GO TO 30
   25       CONTINUE
            IF(ALAS.LT.HELIM) GO TO 30
            XX = XX-ELIM
            S1 = S1*CELM
            S2 = S2*CELM
            ZD = CMPLX(XX,ZRI)
   30    CONTINUE
         NZ = N
         IF(IC.EQ.N) NZ=N-1
         GO TO 45
   40    CONTINUE
         NZ = KK - 2
   45    CONTINUE
         DO 50 K=1,NZ
            Y(K) = CZERO
   50    CONTINUE
         RETURN
      END
      SUBROUTINE CMLRI(Z, FNU, KODE, N, Y, NZ, TOL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CMLRI
C***REFER TO  CBESI,CBESK
C
C     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***ROUTINES CALLED  GAMLN,R1MACH
C***END PROLOGUE  CMLRI
         COMPLEX CK, CNORM, CONE, CTWO, CZERO, PT, P1, P2, RZ, SUM, Y, Z
         REAL ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF, FNU, RHO,
     *    RHO2, SCLE, TFNF, TOL, TST, X, GAMLN, R1MACH
         INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N
         DIMENSION Y(N)
         DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
         SCLE = 1.0E+3*R1MACH(1)/TOL
         NZ=0
         AZ = CABS(Z)
         X = REAL(Z)
         IAZ = INT(AZ)
         IFNU = INT(FNU)
         INU = IFNU + N - 1
         AT = FLOAT(IAZ) + 1.0E0
         CK = CMPLX(AT,0.0E0)/Z
         RZ = CTWO/Z
         P1 = CZERO
         P2 = CONE
         ACK = (AT+1.0E0)/AZ
         RHO = ACK + SQRT(ACK*ACK-1.0E0)
         RHO2 = RHO*RHO
         TST = (RHO2+RHO2)/((RHO2-1.0E0)*(RHO-1.0E0))
         TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
         AK = AT
         DO 10 I=1,80
            PT = P2
            P2 = P1 - CK*P2
            P1 = PT
            CK = CK + RZ
            AP = CABS(P2)
            IF (AP.GT.TST*AK*AK) GO TO 20
            AK = AK + 1.0E0
   10    CONTINUE
         GO TO 110
   20    CONTINUE
         I = I + 1
         K = 0
         IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
         P1 = CZERO
         P2 = CONE
         AT = FLOAT(INU) + 1.0E0
         CK = CMPLX(AT,0.0E0)/Z
         ACK = AT/AZ
         TST = SQRT(ACK/TOL)
         ITIME = 1
         DO 30 K=1,80
            PT = P2
            P2 = P1 - CK*P2
            P1 = PT
            CK = CK + RZ
            AP = CABS(P2)
            IF (AP.LT.TST) GO TO 30
            IF (ITIME.EQ.2) GO TO 40
            ACK = CABS(CK)
            FLAM = ACK + SQRT(ACK*ACK-1.0E0)
            FKAP = AP/CABS(P1)
            RHO = AMIN1(FLAM,FKAP)
            TST = TST*SQRT(RHO/(RHO*RHO-1.0E0))
            ITIME = 2
   30    CONTINUE
         GO TO 110
   40    CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
         K = K + 1
         KK = MAX0(I+IAZ,K+INU)
         FKK = FLOAT(KK)
         P1 = CZERO
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
         P2 = CMPLX(SCLE,0.0E0)
         FNF = FNU - FLOAT(IFNU)
         TFNF = FNF + FNF
         BK = GAMLN(FKK+TFNF+1.0E0,IDUM) - GAMLN(FKK+1.0E0,IDUM)
     *        -GAMLN(TFNF+1.0E0,IDUM)
         BK = EXP(BK)
         SUM = CZERO
         KM = KK - INU
         DO 50 I=1,KM
            PT = P2
            P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
            P1 = PT
            AK = 1.0E0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
            BK = ACK
            FKK = FKK - 1.0E0
   50    CONTINUE
         Y(N) = P2
         IF (N.EQ.1) GO TO 70
         DO 60 I=2,N
            PT = P2
            P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
            P1 = PT
            AK = 1.0E0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
            BK = ACK
            FKK = FKK - 1.0E0
            M = N - I + 1
            Y(M) = P2
   60    CONTINUE
   70    CONTINUE
         IF (IFNU.LE.0) GO TO 90
         DO 80 I=1,IFNU
            PT = P2
            P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
            P1 = PT
            AK = 1.0E0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
            BK = ACK
            FKK = FKK - 1.0E0
   80    CONTINUE
   90    CONTINUE
         PT = Z
         IF (KODE.EQ.2) PT = PT - CMPLX(X,0.0E0)
         P1 = -CMPLX(FNF,0.0E0)*CLOG(RZ) + PT
         AP = GAMLN(1.0E0+FNF,IDUM)
         PT = P1 - CMPLX(AP,0.0E0)
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
         P2 = P2 + SUM
         AP = CABS(P2)
         P1 = CMPLX(1.0E0/AP,0.0E0)
         CK = CEXP(PT)*P1
         PT = CONJG(P2)*P1
         CNORM = CK*PT
         DO 100 I=1,N
            Y(I) = Y(I)*CNORM
  100    CONTINUE
         RETURN
  110    CONTINUE
         NZ=-2
         RETURN
      END
      SUBROUTINE CRATI(Z, FNU, N, CY, TOL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CRATI
C***REFER TO  CBESI,CBESK,CBESH
C
C     CRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CRATI
         COMPLEX CDFNU, CONE, CY, CZERO, PT, P1, P2, RZ, T1, Z
         REAL AK, AMAGZ, AP1, AP2, ARG, AZ, DFNU, FDNU, FLAM, FNU, FNUP,
     *    RAP1, RHO, TEST, TEST1, TOL
         INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
         DIMENSION CY(N)
         DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
         AZ = CABS(Z)
         INU = INT(FNU)
         IDNU = INU + N - 1
         FDNU = FLOAT(IDNU)
         MAGZ = INT(AZ)
         AMAGZ = FLOAT(MAGZ+1)
         FNUP = AMAX1(AMAGZ,FDNU)
         ID = IDNU - MAGZ - 1
         ITIME = 1
         K = 1
         RZ = (CONE+CONE)/Z
         T1 = CMPLX(FNUP,0.0E0)*RZ
         P2 = -T1
         P1 = CONE
         T1 = T1 + RZ
         IF (ID.GT.0) ID = 0
         AP2 = CABS(P2)
         AP1 = CABS(P1)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
         ARG = (AP2+AP2)/(AP1*TOL)
         TEST1 = SQRT(ARG)
         TEST = TEST1
         RAP1 = 1.0E0/AP1
         P1 = P1*CMPLX(RAP1,0.0E0)
         P2 = P2*CMPLX(RAP1,0.0E0)
         AP2 = AP2*RAP1
   10    CONTINUE
         K = K + 1
         AP1 = AP2
         PT = P2
         P2 = P1 - T1*P2
         P1 = PT
         T1 = T1 + RZ
         AP2 = CABS(P2)
         IF (AP1.LE.TEST) GO TO 10
         IF (ITIME.EQ.2) GO TO 20
         AK = CABS(T1)*0.5E0
         FLAM = AK + SQRT(AK*AK-1.0E0)
         RHO = AMIN1(AP2/AP1,FLAM)
         TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0E0))
         ITIME = 2
         GO TO 10
   20    CONTINUE
         KK = K + 1 - ID
         AK = FLOAT(KK)
         DFNU = FNU + FLOAT(N-1)
         CDFNU = CMPLX(DFNU,0.0E0)
         T1 = CMPLX(AK,0.0E0)
         P1 = CMPLX(1.0E0/AP2,0.0E0)
         P2 = CZERO
         DO 30 I=1,KK
            PT = P1
            P1 = RZ*(CDFNU+T1)*P1 + P2
            P2 = PT
            T1 = T1 - CONE
   30    CONTINUE
         IF (REAL(P1).NE.0.0E0 .OR. AIMAG(P1).NE.0.0E0) GO TO 40
         P1 = CMPLX(TOL,TOL)
   40    CONTINUE
         CY(N) = P2/P1
         IF (N.EQ.1) RETURN
         K = N - 1
         AK = FLOAT(K)
         T1 = CMPLX(AK,0.0E0)
         CDFNU = CMPLX(FNU,0.0E0)*RZ
         DO 60 I=2,N
            PT = CDFNU + T1*RZ + CY(K+1)
            IF (REAL(PT).NE.0.0E0 .OR. AIMAG(PT).NE.0.0E0) GO TO 50
            PT = CMPLX(TOL,TOL)
   50       CONTINUE
            CY(K) = CONE/PT
            T1 = T1 - CONE
            K = K - 1
   60    CONTINUE
         RETURN
      END
      SUBROUTINE CS1S2(ZR, S1, S2, NZ, ASCLE, ALIM, IUF)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CS1S2
C***REFER TO  CBESK,CAIRY
C
C     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CS1S2
         COMPLEX CZERO, C1, S1, S1D, S2, ZR
         REAL AA, ALIM, ALN, ASCLE, AS1, AS2, XX
         INTEGER IUF, NZ
         DATA CZERO / (0.0E0,0.0E0) /
         NZ = 0
         AS1 = CABS(S1)
         AS2 = CABS(S2)
         AA = REAL(S1)
         ALN = AIMAG(S1)
         IF (AA.EQ.0.0E0 .AND. ALN.EQ.0.0E0) GO TO 10
         IF (AS1.EQ.0.0E0) GO TO 10
         XX = REAL(ZR)
         ALN = -XX - XX + ALOG(AS1)
         S1D = S1
         S1 = CZERO
         AS1 = 0.0E0
         IF (ALN.LT.(-ALIM)) GO TO 10
         C1 = CLOG(S1D) - ZR - ZR
         S1 = CEXP(C1)
         AS1 = CABS(S1)
         IUF = IUF + 1
   10    CONTINUE
         AA = AMAX1(AS1,AS2)
         IF (AA.GT.ASCLE) RETURN
         S1 = CZERO
         S2 = CZERO
         NZ = 1
         IUF = 0
         RETURN
      END
      SUBROUTINE CSERI(Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CSERI
C***REFER TO  CBESI,CBESK
C
C     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***ROUTINES CALLED  CUCHK,GAMLN,R1MACH
C***END PROLOGUE  CSERI
         COMPLEX AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ, S1, S2, W,
     *    Y, Z
         REAL AA, ACZ, AK, ALIM, ARM, ASCLE, ATOL, AZ, DFNU, ELIM, FNU,
     *    FNUP, RAK1, RS, RTR1, S, SS, TOL, X, GAMLN, R1MACH
         INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NW, NZ
         DIMENSION Y(N), W(2)
         DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
         NZ = 0
         AZ = CABS(Z)
         IF (AZ.EQ.0.0E0) GO TO 150
         X = REAL(Z)
         ARM = 1.0E+3*R1MACH(1)
         RTR1 = SQRT(ARM)
         CRSC = CMPLX(1.0E0,0.0E0)
         IFLAG = 0
         IF (AZ.LT.ARM) GO TO 140
         HZ = Z*CMPLX(0.5E0,0.0E0)
         CZ = CZERO
         IF (AZ.GT.RTR1) CZ = HZ*HZ
         ACZ = CABS(CZ)
         NN = N
         CK = CLOG(HZ)
   10    CONTINUE
         DFNU = FNU + FLOAT(NN-1)
         FNUP = DFNU + 1.0E0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
         AK1 = CK*CMPLX(DFNU,0.0E0)
         AK = GAMLN(FNUP,IDUM)
         AK1 = AK1 - CMPLX(AK,0.0E0)
         IF (KODE.EQ.2) AK1 = AK1 - CMPLX(X,0.0E0)
         RAK1 = REAL(AK1)
         IF (RAK1.GT.(-ELIM)) GO TO 30
   20    CONTINUE
         NZ = NZ + 1
         Y(NN) = CZERO
         IF (ACZ.GT.DFNU) GO TO 170
         NN = NN - 1
         IF (NN.EQ.0) RETURN
         GO TO 10
   30    CONTINUE
         IF (RAK1.GT.(-ALIM)) GO TO 40
         IFLAG = 1
         SS = 1.0E0/TOL
         CRSC = CMPLX(TOL,0.0E0)
         ASCLE = ARM*SS
   40    CONTINUE
         AK = AIMAG(AK1)
         AA = EXP(RAK1)
         IF (IFLAG.EQ.1) AA = AA*SS
         COEF = CMPLX(AA,0.0E0)*CMPLX(COS(AK),SIN(AK))
         ATOL = TOL*ACZ/FNUP
         IL = MIN0(2,NN)
         DO 80 I=1,IL
            DFNU = FNU + FLOAT(NN-I)
            FNUP = DFNU + 1.0E0
            S1 = CONE
            IF (ACZ.LT.TOL*FNUP) GO TO 60
            AK1 = CONE
            AK = FNUP + 2.0E0
            S = FNUP
            AA = 2.0E0
   50       CONTINUE
            RS = 1.0E0/S
            AK1 = AK1*CZ*CMPLX(RS,0.0E0)
            S1 = S1 + AK1
            S = S + AK
            AK = AK + 2.0E0
            AA = AA*ACZ*RS
            IF (AA.GT.ATOL) GO TO 50
   60       CONTINUE
            M = NN - I + 1
            S2 = S1*COEF
            W(I) = S2
            IF (IFLAG.EQ.0) GO TO 70
            CALL CUCHK(S2, NW, ASCLE, TOL)
            IF (NW.NE.0) GO TO 20
   70       CONTINUE
            Y(M) = S2*CRSC
            IF (I.NE.IL) COEF = COEF*CMPLX(DFNU,0.0E0)/HZ
   80    CONTINUE
         IF (NN.LE.2) RETURN
         K = NN - 2
         AK = FLOAT(K)
         RZ = (CONE+CONE)/Z
         IF (IFLAG.EQ.1) GO TO 110
         IB = 3
   90    CONTINUE
         DO 100 I=IB,NN
            Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
            AK = AK - 1.0E0
            K = K - 1
  100    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  110    CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3
C-----------------------------------------------------------------------
         S1 = W(1)
         S2 = W(2)
         DO 120 L=3,NN
            CK = S2
            S2 = S1 + CMPLX(AK+FNU,0.0E0)*RZ*S2
            S1 = CK
            CK = S2*CRSC
            Y(K) = CK
            AK = AK - 1.0E0
            K = K - 1
            IF (CABS(CK).GT.ASCLE) GO TO 130
  120    CONTINUE
         RETURN
  130    CONTINUE
         IB = L + 1
         IF (IB.GT.NN) RETURN
         GO TO 90
  140    CONTINUE
         NZ = N
         IF (FNU.EQ.0.0E0) NZ = NZ - 1
  150    CONTINUE
         Y(1) = CZERO
         IF (FNU.EQ.0.0E0) Y(1) = CONE
         IF (N.EQ.1) RETURN
         DO 160 I=2,N
            Y(I) = CZERO
  160    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
C-----------------------------------------------------------------------
  170    CONTINUE
         NZ = -NZ
         RETURN
      END
      SUBROUTINE CSHCH(Z, CSH, CCH)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CSHCH
C***REFER TO  CBESK,CBESH
C
C     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CSHCH
         COMPLEX CCH, CSH, Z
         REAL CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y, COSH, SINH
         X = REAL(Z)
         Y = AIMAG(Z)
         SH = SINH(X)
         CH = COSH(X)
         SN = SIN(Y)
         CN = COS(Y)
         CSHR = SH*CN
         CSHI = CH*SN
         CSH = CMPLX(CSHR,CSHI)
         CCHR = CH*CN
         CCHI = SH*SN
         CCH = CMPLX(CCHR,CCHI)
         RETURN
      END
      SUBROUTINE CUCHK(Y, NZ, ASCLE, TOL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUCHK
C***REFER TO CSERI,CUOIK,CUNK1,CUNK2,CUNI1,CUNI2,CKSCL
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*R1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDER FLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUCHK
C
         COMPLEX Y
         REAL ASCLE, SS, ST, TOL, YR, YI
         INTEGER NZ
         NZ = 0
         YR = REAL(Y)
         YI = AIMAG(Y)
         YR = ABS(YR)
         YI = ABS(YI)
         ST = AMIN1(YR,YI)
         IF (ST.GT.ASCLE) RETURN
         SS = AMAX1(YR,YI)
         ST=ST/TOL
         IF (SS.LT.ST) NZ = 1
         RETURN
      END
      SUBROUTINE CUNHJ(Z, FNU, IPMTR, TOL, PHI, ARG, ZETA1, ZETA2,
     * ASUM, BSUM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUNHJ
C***REFER TO  CBESI,CBESK
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         CUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUNHJ
         COMPLEX ARG, ASUM, BSUM, CFNU, CONE, CR, CZERO, DR, P, PHI,
     *    PRZTH, PTFN, RFN13, RTZTA, RZTH, SUMA, SUMB, TFN, T2, UP, W, W2,
     *    Z, ZA, ZB, ZC, ZETA, ZETA1, ZETA2, ZTH
         REAL ALFA, ANG, AP, AR, ATOL, AW2, AZTH, BETA, BR, BTOL, C, EX1,
     *    EX2, FNU, FN13, FN23, GAMA, HPI, PI, PP, RFNU, RFNU2, THPI, TOL,
     *    WI, WR, ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR, ASUMR, ASUMI, BSUMR,
     *    BSUMI, TEST, TSTR, TSTI, AC
         INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     *    LRP1, L1, L2, M
         DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     *    AP(30), P(30), UP(14), CR(14), DR(14)
         DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1        AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2        1.00000000000000000E+00,     1.04166666666666667E-01,
     3        8.35503472222222222E-02,     1.28226574556327160E-01,
     4        2.91849026464140464E-01,     8.81627267443757652E-01,
     5        3.32140828186276754E+00,     1.49957629868625547E+01,
     6        7.89230130115865181E+01,     4.74451538868264323E+02,
     7        3.20749009089066193E+03,     2.40865496408740049E+04,
     8        1.98923119169509794E+05,     1.79190200777534383E+06/
         DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1        BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2        1.00000000000000000E+00,    -1.45833333333333333E-01,
     3       -9.87413194444444444E-02,    -1.43312053915895062E-01,
     4       -3.17227202678413548E-01,    -9.42429147957120249E-01,
     5       -3.51120304082635426E+00,    -1.57272636203680451E+01,
     6       -8.22814390971859444E+01,    -4.92355370523670524E+02,
     7       -3.31621856854797251E+03,    -2.48276742452085896E+04,
     8       -2.04526587315129788E+05,    -1.83844491706820990E+06/
         DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1        C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2        C(19), C(20), C(21), C(22), C(23), C(24)/
     3        1.00000000000000000E+00,    -2.08333333333333333E-01,
     4        1.25000000000000000E-01,     3.34201388888888889E-01,
     5       -4.01041666666666667E-01,     7.03125000000000000E-02,
     6       -1.02581259645061728E+00,     1.84646267361111111E+00,
     7       -8.91210937500000000E-01,     7.32421875000000000E-02,
     8        4.66958442342624743E+00,    -1.12070026162229938E+01,
     9        8.78912353515625000E+00,    -2.36408691406250000E+00,
     A        1.12152099609375000E-01,    -2.82120725582002449E+01,
     B        8.46362176746007346E+01,    -9.18182415432400174E+01,
     C        4.25349987453884549E+01,    -7.36879435947963170E+00,
     D        2.27108001708984375E-01,     2.12570130039217123E+02,
     E       -7.65252468141181642E+02,     1.05999045252799988E+03/
         DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1        C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2        C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3       -6.99579627376132541E+02,     2.18190511744211590E+02,
     4       -2.64914304869515555E+01,     5.72501420974731445E-01,
     5       -1.91945766231840700E+03,     8.06172218173730938E+03,
     6       -1.35865500064341374E+04,     1.16553933368645332E+04,
     7       -5.30564697861340311E+03,     1.20090291321635246E+03,
     8       -1.08090919788394656E+02,     1.72772750258445740E+00,
     9        2.02042913309661486E+04,    -9.69805983886375135E+04,
     A        1.92547001232531532E+05,    -2.03400177280415534E+05,
     B        1.22200464983017460E+05,    -4.11926549688975513E+04,
     C        7.10951430248936372E+03,    -4.93915304773088012E+02,
     D        6.07404200127348304E+00,    -2.42919187900551333E+05,
     E        1.31176361466297720E+06,    -2.99801591853810675E+06/
         DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1        C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2        C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3        3.76327129765640400E+06,    -2.81356322658653411E+06,
     4        1.26836527332162478E+06,    -3.31645172484563578E+05,
     5        4.52187689813627263E+04,    -2.49983048181120962E+03,
     6        2.43805296995560639E+01,     3.28446985307203782E+06,
     7       -1.97068191184322269E+07,     5.09526024926646422E+07,
     8       -7.41051482115326577E+07,     6.63445122747290267E+07,
     9       -3.75671766607633513E+07,     1.32887671664218183E+07,
     A       -2.78561812808645469E+06,     3.08186404612662398E+05,
     B       -1.38860897537170405E+04,     1.10017140269246738E+02,
     C       -4.93292536645099620E+07,     3.25573074185765749E+08,
     D       -9.39462359681578403E+08,     1.55359689957058006E+09,
     E       -1.62108055210833708E+09,     1.10684281682301447E+09/
         DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1        C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2        C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3       -4.95889784275030309E+08,     1.42062907797533095E+08,
     4       -2.44740627257387285E+07,     2.24376817792244943E+06,
     5       -8.40054336030240853E+04,     5.51335896122020586E+02,
     6        8.14789096118312115E+08,    -5.86648149205184723E+09,
     7        1.86882075092958249E+10,    -3.46320433881587779E+10,
     8        4.12801855797539740E+10,    -3.30265997498007231E+10,
     9        1.79542137311556001E+10,    -6.56329379261928433E+09,
     A        1.55927986487925751E+09,    -2.25105661889415278E+08,
     B        1.73951075539781645E+07,    -5.49842327572288687E+05,
     C        3.03809051092238427E+03,    -1.46792612476956167E+10,
     D        1.14498237732025810E+11,    -3.99096175224466498E+11,
     E        8.19218669548577329E+11,    -1.09837515608122331E+12/
         DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1        C(105)/
     2        1.00815810686538209E+12,    -6.45364869245376503E+11,
     3        2.87900649906150589E+11,    -8.78670721780232657E+10,
     4        1.76347306068349694E+10,    -2.16716498322379509E+09,
     5        1.43157876718888981E+08,    -3.87183344257261262E+06,
     6        1.82577554742931747E+04/
         DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1        ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2        ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3        ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4       -4.44444444444444444E-03,    -9.22077922077922078E-04,
     5       -8.84892884892884893E-05,     1.65927687832449737E-04,
     6        2.46691372741792910E-04,     2.65995589346254780E-04,
     7        2.61824297061500945E-04,     2.48730437344655609E-04,
     8        2.32721040083232098E-04,     2.16362485712365082E-04,
     9        2.00738858762752355E-04,     1.86267636637545172E-04,
     A        1.73060775917876493E-04,     1.61091705929015752E-04,
     B        1.50274774160908134E-04,     1.40503497391269794E-04,
     C        1.31668816545922806E-04,     1.23667445598253261E-04,
     D        1.16405271474737902E-04,     1.09798298372713369E-04,
     E        1.03772410422992823E-04,     9.82626078369363448E-05/
         DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1        ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2        ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3        ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4        9.32120517249503256E-05,     8.85710852478711718E-05,
     5        8.42963105715700223E-05,     8.03497548407791151E-05,
     6        7.66981345359207388E-05,     7.33122157481777809E-05,
     7        7.01662625163141333E-05,     6.72375633790160292E-05,
     8        6.93735541354588974E-04,     2.32241745182921654E-04,
     9       -1.41986273556691197E-05,    -1.16444931672048640E-04,
     A       -1.50803558053048762E-04,    -1.55121924918096223E-04,
     B       -1.46809756646465549E-04,    -1.33815503867491367E-04,
     C       -1.19744975684254051E-04,    -1.06184319207974020E-04,
     D       -9.37699549891194492E-05,    -8.26923045588193274E-05,
     E       -7.29374348155221211E-05,    -6.44042357721016283E-05/
         DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1        ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2        ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3        ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4       -5.69611566009369048E-05,    -5.04731044303561628E-05,
     5       -4.48134868008882786E-05,    -3.98688727717598864E-05,
     6       -3.55400532972042498E-05,    -3.17414256609022480E-05,
     7       -2.83996793904174811E-05,    -2.54522720634870566E-05,
     8       -2.28459297164724555E-05,    -2.05352753106480604E-05,
     9       -1.84816217627666085E-05,    -1.66519330021393806E-05,
     A       -1.50179412980119482E-05,    -1.35554031379040526E-05,
     B       -1.22434746473858131E-05,    -1.10641884811308169E-05,
     C       -3.54211971457743841E-04,    -1.56161263945159416E-04,
     D        3.04465503594936410E-05,     1.30198655773242693E-04,
     E        1.67471106699712269E-04,     1.70222587683592569E-04/
         DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1        ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2        ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3        ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4        1.56501427608594704E-04,     1.36339170977445120E-04,
     5        1.14886692029825128E-04,     9.45869093034688111E-05,
     6        7.64498419250898258E-05,     6.07570334965197354E-05,
     7        4.74394299290508799E-05,     3.62757512005344297E-05,
     8        2.69939714979224901E-05,     1.93210938247939253E-05,
     9        1.30056674793963203E-05,     7.82620866744496661E-06,
     A        3.59257485819351583E-06,     1.44040049814251817E-07,
     B       -2.65396769697939116E-06,    -4.91346867098485910E-06,
     C       -6.72739296091248287E-06,    -8.17269379678657923E-06,
     D       -9.31304715093561232E-06,    -1.02011418798016441E-05,
     E       -1.08805962510592880E-05,    -1.13875481509603555E-05/
         DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1        ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2        ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3        ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4       -1.17519675674556414E-05,    -1.19987364870944141E-05,
     5        3.78194199201772914E-04,     2.02471952761816167E-04,
     6       -6.37938506318862408E-05,    -2.38598230603005903E-04,
     7       -3.10916256027361568E-04,    -3.13680115247576316E-04,
     8       -2.78950273791323387E-04,    -2.28564082619141374E-04,
     9       -1.75245280340846749E-04,    -1.25544063060690348E-04,
     A       -8.22982872820208365E-05,    -4.62860730588116458E-05,
     B       -1.72334302366962267E-05,     5.60690482304602267E-06,
     C        2.31395443148286800E-05,     3.62642745856793957E-05,
     D        4.58006124490188752E-05,     5.24595294959114050E-05,
     E        5.68396208545815266E-05,     5.94349820393104052E-05/
         DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1        ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2        ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3        ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4        6.06478527578421742E-05,     6.08023907788436497E-05,
     5        6.01577894539460388E-05,     5.89199657344698500E-05,
     6        5.72515823777593053E-05,     5.52804375585852577E-05,
     7        5.31063773802880170E-05,     5.08069302012325706E-05,
     8        4.84418647620094842E-05,     4.60568581607475370E-05,
     9       -6.91141397288294174E-04,    -4.29976633058871912E-04,
     A        1.83067735980039018E-04,     6.60088147542014144E-04,
     B        8.75964969951185931E-04,     8.77335235958235514E-04,
     C        7.49369585378990637E-04,     5.63832329756980918E-04,
     D        3.68059319971443156E-04,     1.88464535514455599E-04/
         DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1        ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2        ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3        ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4        3.70663057664904149E-05,    -8.28520220232137023E-05,
     5       -1.72751952869172998E-04,    -2.36314873605872983E-04,
     6       -2.77966150694906658E-04,    -3.02079514155456919E-04,
     7       -3.12594712643820127E-04,    -3.12872558758067163E-04,
     8       -3.05678038466324377E-04,    -2.93226470614557331E-04,
     9       -2.77255655582934777E-04,    -2.59103928467031709E-04,
     A       -2.39784014396480342E-04,    -2.20048260045422848E-04,
     B       -2.00443911094971498E-04,    -1.81358692210970687E-04,
     C       -1.63057674478657464E-04,    -1.45712672175205844E-04,
     D       -1.29425421983924587E-04,    -1.14245691942445952E-04/
         DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1        ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2        ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3        ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4        1.92821964248775885E-03,     1.35592576302022234E-03,
     5       -7.17858090421302995E-04,    -2.58084802575270346E-03,
     6       -3.49271130826168475E-03,    -3.46986299340960628E-03,
     7       -2.82285233351310182E-03,    -1.88103076404891354E-03,
     8       -8.89531718383947600E-04,     3.87912102631035228E-06,
     9        7.28688540119691412E-04,     1.26566373053457758E-03,
     A        1.62518158372674427E-03,     1.83203153216373172E-03,
     B        1.91588388990527909E-03,     1.90588846755546138E-03,
     C        1.82798982421825727E-03,     1.70389506421121530E-03,
     D        1.55097127171097686E-03,     1.38261421852276159E-03/
         DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1        ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2        1.20881424230064774E-03,     1.03676532638344962E-03,
     3        8.71437918068619115E-04,     7.16080155297701002E-04,
     4        5.72637002558129372E-04,     4.42089819465802277E-04,
     5        3.24724948503090564E-04,     2.20342042730246599E-04,
     6        1.28412898401353882E-04,     4.82005924552095464E-05/
         DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1        BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2        BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3        BETA(19), BETA(20), BETA(21), BETA(22)/
     4        1.79988721413553309E-02,     5.59964911064388073E-03,
     5        2.88501402231132779E-03,     1.80096606761053941E-03,
     6        1.24753110589199202E-03,     9.22878876572938311E-04,
     7        7.14430421727287357E-04,     5.71787281789704872E-04,
     8        4.69431007606481533E-04,     3.93232835462916638E-04,
     9        3.34818889318297664E-04,     2.88952148495751517E-04,
     A        2.52211615549573284E-04,     2.22280580798883327E-04,
     B        1.97541838033062524E-04,     1.76836855019718004E-04,
     C        1.59316899661821081E-04,     1.44347930197333986E-04,
     D        1.31448068119965379E-04,     1.20245444949302884E-04,
     E        1.10449144504599392E-04,     1.01828770740567258E-04/
         DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1        BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2        BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3        BETA(41), BETA(42), BETA(43), BETA(44)/
     4        9.41998224204237509E-05,     8.74130545753834437E-05,
     5        8.13466262162801467E-05,     7.59002269646219339E-05,
     6        7.09906300634153481E-05,     6.65482874842468183E-05,
     7        6.25146958969275078E-05,     5.88403394426251749E-05,
     8       -1.49282953213429172E-03,    -8.78204709546389328E-04,
     9       -5.02916549572034614E-04,    -2.94822138512746025E-04,
     A       -1.75463996970782828E-04,    -1.04008550460816434E-04,
     B       -5.96141953046457895E-05,    -3.12038929076098340E-05,
     C       -1.26089735980230047E-05,    -2.42892608575730389E-07,
     D        8.05996165414273571E-06,     1.36507009262147391E-05,
     E        1.73964125472926261E-05,     1.98672978842133780E-05/
         DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1        BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2        BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3        BETA(63), BETA(64), BETA(65), BETA(66)/
     4        2.14463263790822639E-05,     2.23954659232456514E-05,
     5        2.28967783814712629E-05,     2.30785389811177817E-05,
     6        2.30321976080909144E-05,     2.28236073720348722E-05,
     7        2.25005881105292418E-05,     2.20981015361991429E-05,
     8        2.16418427448103905E-05,     2.11507649256220843E-05,
     9        2.06388749782170737E-05,     2.01165241997081666E-05,
     A        1.95913450141179244E-05,     1.90689367910436740E-05,
     B        1.85533719641636667E-05,     1.80475722259674218E-05,
     C        5.52213076721292790E-04,     4.47932581552384646E-04,
     D        2.79520653992020589E-04,     1.52468156198446602E-04,
     E        6.93271105657043598E-05,     1.76258683069991397E-05/
         DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1        BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2        BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3        BETA(85), BETA(86), BETA(87), BETA(88)/
     4       -1.35744996343269136E-05,    -3.17972413350427135E-05,
     5       -4.18861861696693365E-05,    -4.69004889379141029E-05,
     6       -4.87665447413787352E-05,    -4.87010031186735069E-05,
     7       -4.74755620890086638E-05,    -4.55813058138628452E-05,
     8       -4.33309644511266036E-05,    -4.09230193157750364E-05,
     9       -3.84822638603221274E-05,    -3.60857167535410501E-05,
     A       -3.37793306123367417E-05,    -3.15888560772109621E-05,
     B       -2.95269561750807315E-05,    -2.75978914828335759E-05,
     C       -2.58006174666883713E-05,    -2.41308356761280200E-05,
     D       -2.25823509518346033E-05,    -2.11479656768912971E-05,
     E       -1.98200638885294927E-05,    -1.85909870801065077E-05/
         DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1        BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2        BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3        BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4       -1.74532699844210224E-05,    -1.63997823854497997E-05,
     5       -4.74617796559959808E-04,    -4.77864567147321487E-04,
     6       -3.20390228067037603E-04,    -1.61105016119962282E-04,
     7       -4.25778101285435204E-05,     3.44571294294967503E-05,
     8        7.97092684075674924E-05,     1.03138236708272200E-04,
     9        1.12466775262204158E-04,     1.13103642108481389E-04,
     A        1.08651634848774268E-04,     1.01437951597661973E-04,
     B        9.29298396593363896E-05,     8.40293133016089978E-05,
     C        7.52727991349134062E-05,     6.69632521975730872E-05,
     D        5.92564547323194704E-05,     5.22169308826975567E-05,
     E        4.58539485165360646E-05,     4.01445513891486808E-05/
         DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1        BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2        BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3        BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4        3.50481730031328081E-05,     3.05157995034346659E-05,
     5        2.64956119950516039E-05,     2.29363633690998152E-05,
     6        1.97893056664021636E-05,     1.70091984636412623E-05,
     7        1.45547428261524004E-05,     1.23886640995878413E-05,
     8        1.04775876076583236E-05,     8.79179954978479373E-06,
     9        7.36465810572578444E-04,     8.72790805146193976E-04,
     A        6.22614862573135066E-04,     2.85998154194304147E-04,
     B        3.84737672879366102E-06,    -1.87906003636971558E-04,
     C       -2.97603646594554535E-04,    -3.45998126832656348E-04,
     D       -3.53382470916037712E-04,    -3.35715635775048757E-04/
         DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1        BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2        BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3        BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4       -3.04321124789039809E-04,    -2.66722723047612821E-04,
     5       -2.27654214122819527E-04,    -1.89922611854562356E-04,
     6       -1.55058918599093870E-04,    -1.23778240761873630E-04,
     7       -9.62926147717644187E-05,    -7.25178327714425337E-05,
     8       -5.22070028895633801E-05,    -3.50347750511900522E-05,
     9       -2.06489761035551757E-05,    -8.70106096849767054E-06,
     A        1.13698686675100290E-06,     9.16426474122778849E-06,
     B        1.56477785428872620E-05,     2.08223629482466847E-05,
     C        2.48923381004595156E-05,     2.80340509574146325E-05,
     D        3.03987774629861915E-05,     3.21156731406700616E-05/
         DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1        BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2        BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3        BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4       -1.80182191963885708E-03,    -2.43402962938042533E-03,
     5       -1.83422663549856802E-03,    -7.62204596354009765E-04,
     6        2.39079475256927218E-04,     9.49266117176881141E-04,
     7        1.34467449701540359E-03,     1.48457495259449178E-03,
     8        1.44732339830617591E-03,     1.30268261285657186E-03,
     9        1.10351597375642682E-03,     8.86047440419791759E-04,
     A        6.73073208165665473E-04,     4.77603872856582378E-04,
     B        3.05991926358789362E-04,     1.60315694594721630E-04,
     C        4.00749555270613286E-05,    -5.66607461635251611E-05,
     D       -1.32506186772982638E-04,    -1.90296187989614057E-04/
         DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1        BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2        BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3        BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4       -2.32811450376937408E-04,    -2.62628811464668841E-04,
     5       -2.82050469867598672E-04,    -2.93081563192861167E-04,
     6       -2.97435962176316616E-04,    -2.96557334239348078E-04,
     7       -2.91647363312090861E-04,    -2.83696203837734166E-04,
     8       -2.73512317095673346E-04,    -2.61750155806768580E-04,
     9        6.38585891212050914E-03,     9.62374215806377941E-03,
     A        7.61878061207001043E-03,     2.83219055545628054E-03,
     B       -2.09841352012720090E-03,    -5.73826764216626498E-03,
     C       -7.70804244495414620E-03,    -8.21011692264844401E-03,
     D       -7.65824520346905413E-03,    -6.47209729391045177E-03/
         DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1        BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2        BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3        BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4       -4.99132412004966473E-03,    -3.45612289713133280E-03,
     5       -2.01785580014170775E-03,    -7.59430686781961401E-04,
     6        2.84173631523859138E-04,     1.10891667586337403E-03,
     7        1.72901493872728771E-03,     2.16812590802684701E-03,
     8        2.45357710494539735E-03,     2.61281821058334862E-03,
     9        2.67141039656276912E-03,     2.65203073395980430E-03,
     A        2.57411652877287315E-03,     2.45389126236094427E-03,
     B        2.30460058071795494E-03,     2.13684837686712662E-03,
     C        1.95896528478870911E-03,     1.77737008679454412E-03,
     D        1.59690280765839059E-03,     1.42111975664438546E-03/
         DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1        GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2        GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3        GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4        6.29960524947436582E-01,     2.51984209978974633E-01,
     5        1.54790300415655846E-01,     1.10713062416159013E-01,
     6        8.57309395527394825E-02,     6.97161316958684292E-02,
     7        5.86085671893713576E-02,     5.04698873536310685E-02,
     8        4.42600580689154809E-02,     3.93720661543509966E-02,
     9        3.54283195924455368E-02,     3.21818857502098231E-02,
     A        2.94646240791157679E-02,     2.71581677112934479E-02,
     B        2.51768272973861779E-02,     2.34570755306078891E-02,
     C        2.19508390134907203E-02,     2.06210828235646240E-02,
     D        1.94388240897880846E-02,     1.83810633800683158E-02,
     E        1.74293213231963172E-02,     1.65685837786612353E-02/
         DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1        GAMA(29), GAMA(30)/
     2        1.57865285987918445E-02,     1.50729501494095594E-02,
     3        1.44193250839954639E-02,     1.38184805735341786E-02,
     4        1.32643378994276568E-02,     1.27517121970498651E-02,
     5        1.22761545318762767E-02,     1.18338262398482403E-02/
         DATA EX1, EX2, HPI, PI, THPI /
     1        3.33333333333333333E-01,     6.66666666666666667E-01,
     2        1.57079632679489662E+00,     3.14159265358979324E+00,
     3        4.71238898038468986E+00/
         DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
         RFNU = 1.0E0/FNU
C     ZB = Z*CMPLX(RFNU,0.0E0)
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
         TSTR = REAL(Z)
         TSTI = AIMAG(Z)
         TEST = R1MACH(1)*1.0E+3
         AC = FNU*TEST
         IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
         AC = 2.0E0*ABS(ALOG(TEST))+FNU
         ZETA1 = CMPLX(AC,0.0E0)
         ZETA2 = CMPLX(FNU,0.0E0)
         PHI=CONE
         ARG=CONE
         RETURN
   15    CONTINUE
         ZB = Z*CMPLX(RFNU,0.0E0)
         RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
         FN13 = FNU**EX1
         FN23 = FN13*FN13
         RFN13 = CMPLX(1.0E0/FN13,0.0E0)
         W2 = CONE - ZB*ZB
         AW2 = CABS(W2)
         IF (AW2.GT.0.25E0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(W2).LE.0.25E0
C-----------------------------------------------------------------------
         K = 1
         P(1) = CONE
         SUMA = CMPLX(GAMA(1),0.0E0)
         AP(1) = 1.0E0
         IF (AW2.LT.TOL) GO TO 20
         DO 10 K=2,30
            P(K) = P(K-1)*W2
            SUMA = SUMA + P(K)*CMPLX(GAMA(K),0.0E0)
            AP(K) = AP(K-1)*AW2
            IF (AP(K).LT.TOL) GO TO 20
   10    CONTINUE
         K = 30
   20    CONTINUE
         KMAX = K
         ZETA = W2*SUMA
         ARG = ZETA*CMPLX(FN23,0.0E0)
         ZA = CSQRT(SUMA)
         ZETA2 = CSQRT(W2)*CMPLX(FNU,0.0E0)
         ZETA1 = ZETA2*(CONE+ZETA*ZA*CMPLX(EX2,0.0E0))
         ZA = ZA + ZA
         PHI = CSQRT(ZA)*RFN13
         IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
         SUMB = CZERO
         DO 30 K=1,KMAX
            SUMB = SUMB + P(K)*CMPLX(BETA(K),0.0E0)
   30    CONTINUE
         ASUM = CZERO
         BSUM = SUMB
         L1 = 0
         L2 = 30
         BTOL = TOL*CABS(BSUM)
         ATOL = TOL
         PP = 1.0E0
         IAS = 0
         IBS = 0
         IF (RFNU2.LT.TOL) GO TO 110
         DO 100 IS=2,7
            ATOL = ATOL/RFNU2
            PP = PP*RFNU2
            IF (IAS.EQ.1) GO TO 60
            SUMA = CZERO
            DO 40 K=1,KMAX
               M = L1 + K
               SUMA = SUMA + P(K)*CMPLX(ALFA(M),0.0E0)
               IF (AP(K).LT.ATOL) GO TO 50
   40       CONTINUE
   50       CONTINUE
            ASUM = ASUM + SUMA*CMPLX(PP,0.0E0)
            IF (PP.LT.TOL) IAS = 1
   60       CONTINUE
            IF (IBS.EQ.1) GO TO 90
            SUMB = CZERO
            DO 70 K=1,KMAX
               M = L2 + K
               SUMB = SUMB + P(K)*CMPLX(BETA(M),0.0E0)
               IF (AP(K).LT.ATOL) GO TO 80
   70       CONTINUE
   80       CONTINUE
            BSUM = BSUM + SUMB*CMPLX(PP,0.0E0)
            IF (PP.LT.BTOL) IBS = 1
   90       CONTINUE
            IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
            L1 = L1 + 30
            L2 = L2 + 30
  100    CONTINUE
  110    CONTINUE
         ASUM = ASUM + CONE
         PP = RFNU*REAL(RFN13)
         BSUM = BSUM*CMPLX(PP,0.0E0)
  120    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     CABS(W2).GT.0.25E0
C-----------------------------------------------------------------------
  130    CONTINUE
         W = CSQRT(W2)
         WR = REAL(W)
         WI = AIMAG(W)
         IF (WR.LT.0.0E0) WR = 0.0E0
         IF (WI.LT.0.0E0) WI = 0.0E0
         W = CMPLX(WR,WI)
         ZA = (CONE+W)/ZB
         ZC = CLOG(ZA)
         ZCR = REAL(ZC)
         ZCI = AIMAG(ZC)
         IF (ZCI.LT.0.0E0) ZCI = 0.0E0
         IF (ZCI.GT.HPI) ZCI = HPI
         IF (ZCR.LT.0.0E0) ZCR = 0.0E0
         ZC = CMPLX(ZCR,ZCI)
         ZTH = (ZC-W)*CMPLX(1.5E0,0.0E0)
         CFNU = CMPLX(FNU,0.0E0)
         ZETA1 = ZC*CFNU
         ZETA2 = W*CFNU
         AZTH = CABS(ZTH)
         ZTHR = REAL(ZTH)
         ZTHI = AIMAG(ZTH)
         ANG = THPI
         IF (ZTHR.GE.0.0E0 .AND. ZTHI.LT.0.0E0) GO TO 140
         ANG = HPI
         IF (ZTHR.EQ.0.0E0) GO TO 140
         ANG = ATAN(ZTHI/ZTHR)
         IF (ZTHR.LT.0.0E0) ANG = ANG + PI
  140    CONTINUE
         PP = AZTH**EX2
         ANG = ANG*EX2
         ZETAR = PP*COS(ANG)
         ZETAI = PP*SIN(ANG)
         IF (ZETAI.LT.0.0E0) ZETAI = 0.0E0
         ZETA = CMPLX(ZETAR,ZETAI)
         ARG = ZETA*CMPLX(FN23,0.0E0)
         RTZTA = ZTH/ZETA
         ZA = RTZTA/W
         PHI = CSQRT(ZA+ZA)*RFN13
         IF (IPMTR.EQ.1) GO TO 120
         TFN = CMPLX(RFNU,0.0E0)/W
         RZTH = CMPLX(RFNU,0.0E0)/ZTH
         ZC = RZTH*CMPLX(AR(2),0.0E0)
         T2 = CONE/W2
         UP(2) = (T2*CMPLX(C(2),0.0E0)+CMPLX(C(3),0.0E0))*TFN
         BSUM = UP(2) + ZC
         ASUM = CZERO
         IF (RFNU.LT.TOL) GO TO 220
         PRZTH = RZTH
         PTFN = TFN
         UP(1) = CONE
         PP = 1.0E0
         BSUMR = REAL(BSUM)
         BSUMI = AIMAG(BSUM)
         BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
         KS = 0
         KP1 = 2
         L = 3
         IAS = 0
         IBS = 0
         DO 210 LR=2,12,2
            LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
            DO 160 K=LR,LRP1
               KS = KS + 1
               KP1 = KP1 + 1
               L = L + 1
               ZA = CMPLX(C(L),0.0E0)
               DO 150 J=2,KP1
                  L = L + 1
                  ZA = ZA*T2 + CMPLX(C(L),0.0E0)
  150          CONTINUE
               PTFN = PTFN*TFN
               UP(KP1) = PTFN*ZA
               CR(KS) = PRZTH*CMPLX(BR(KS+1),0.0E0)
               PRZTH = PRZTH*RZTH
               DR(KS) = PRZTH*CMPLX(AR(KS+2),0.0E0)
  160       CONTINUE
            PP = PP*RFNU2
            IF (IAS.EQ.1) GO TO 180
            SUMA = UP(LRP1)
            JU = LRP1
            DO 170 JR=1,LR
               JU = JU - 1
               SUMA = SUMA + CR(JR)*UP(JU)
  170       CONTINUE
            ASUM = ASUM + SUMA
            ASUMR = REAL(ASUM)
            ASUMI = AIMAG(ASUM)
            TEST = ABS(ASUMR) + ABS(ASUMI)
            IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180       CONTINUE
            IF (IBS.EQ.1) GO TO 200
            SUMB = UP(LR+2) + UP(LRP1)*ZC
            JU = LRP1
            DO 190 JR=1,LR
               JU = JU - 1
               SUMB = SUMB + DR(JR)*UP(JU)
  190       CONTINUE
            BSUM = BSUM + SUMB
            BSUMR = REAL(BSUM)
            BSUMI = AIMAG(BSUM)
            TEST = ABS(BSUMR) + ABS(BSUMI)
            IF (PP.LT.BTOL .AND. TEST.LT.TOL) IBS = 1
  200       CONTINUE
            IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210    CONTINUE
  220    CONTINUE
         ASUM = ASUM + CONE
         BSUM = -BSUM*RFN13/RTZTA
         GO TO 120
      END
      SUBROUTINE CUNI1(Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUNI1
C***REFER TO  CBESI,CBESK
C
C     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  CUCHK,CUNIK,CUOIK,R1MACH
C***END PROLOGUE  CUNI1
         COMPLEX CFN, CONE, CRSC, CSCL, CSR, CSS, CWRK, CZERO, C1, C2,
     *    PHI, RZ, SUM, S1, S2, Y, Z, ZETA1, ZETA2, CY
         REAL ALIM, APHI, ASCLE, BRY, C2I, C2M, C2R, ELIM, FN, FNU, FNUL,
     *    RS1, TOL, YY, R1MACH
         INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
         DIMENSION BRY(3), Y(N), CWRK(16), CSS(3), CSR(3), CY(2)
         DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
         NZ = 0
         ND = N
         NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         CRSC = CMPLX(TOL,0.0E0)
         CSS(1) = CSCL
         CSS(2) = CONE
         CSS(3) = CRSC
         CSR(1) = CRSC
         CSR(2) = CONE
         CSR(3) = CSCL
         BRY(1) = 1.0E+3*R1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
         FN = AMAX1(FNU,1.0E0)
         INIT = 0
         CALL CUNIK(Z, FN, 1, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
         IF (KODE.EQ.1) GO TO 10
         CFN = CMPLX(FN,0.0E0)
         S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
         GO TO 20
   10    CONTINUE
         S1 = -ZETA1 + ZETA2
   20    CONTINUE
         RS1 = REAL(S1)
         IF (ABS(RS1).GT.ELIM) GO TO 130
   30    CONTINUE
         NN = MIN0(2,ND)
         DO 80 I=1,NN
            FN = FNU + FLOAT(ND-I)
            INIT = 0
            CALL CUNIK(Z, FN, 1, 0, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
            IF (KODE.EQ.1) GO TO 40
            CFN = CMPLX(FN,0.0E0)
            YY = AIMAG(Z)
            S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + CMPLX(0.0E0,YY)
            GO TO 50
   40       CONTINUE
            S1 = -ZETA1 + ZETA2
   50       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = REAL(S1)
            IF (ABS(RS1).GT.ELIM) GO TO 110
            IF (I.EQ.1) IFLAG = 2
            IF (ABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = CABS(PHI)
            RS1 = RS1 + ALOG(APHI)
            IF (ABS(RS1).GT.ELIM) GO TO 110
            IF (I.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0E0) GO TO 60
            IF (I.EQ.1) IFLAG = 3
   60       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF CABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
            S2 = PHI*SUM
            C2R = REAL(S1)
            C2I = AIMAG(S1)
            C2M = EXP(C2R)*REAL(CSS(IFLAG))
            S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (IFLAG.NE.1) GO TO 70
            CALL CUCHK(S2, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 110
   70       CONTINUE
            M = ND - I + 1
            CY(I) = S2
            Y(M) = S2*CSR(IFLAG)
   80    CONTINUE
         IF (ND.LE.2) GO TO 100
         RZ = CMPLX(2.0E0,0.0E0)/Z
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = R1MACH(2)
         S1 = CY(1)
         S2 = CY(2)
         C1 = CSR(IFLAG)
         ASCLE = BRY(IFLAG)
         K = ND - 2
         FN = FLOAT(K)
         DO 90 I=3,ND
            C2 = S2
            S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
            S1 = C2
            C2 = S2*C1
            Y(K) = C2
            K = K - 1
            FN = FN - 1.0E0
            IF (IFLAG.GE.3) GO TO 90
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = AMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 90
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1 = S1*C1
            S2 = C2
            S1 = S1*CSS(IFLAG)
            S2 = S2*CSS(IFLAG)
            C1 = CSR(IFLAG)
   90    CONTINUE
  100    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110    CONTINUE
         IF (RS1.GT.0.0E0) GO TO 120
         Y(ND) = CZERO
         NZ = NZ + 1
         ND = ND - 1
         IF (ND.EQ.0) GO TO 100
         CALL CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
         IF (NUF.LT.0) GO TO 120
         ND = ND - NUF
         NZ = NZ + NUF
         IF (ND.EQ.0) GO TO 100
         FN = FNU + FLOAT(ND-1)
         IF (FN.GE.FNUL) GO TO 30
         NLAST = ND
         RETURN
  120    CONTINUE
         NZ = -1
         RETURN
  130    CONTINUE
         IF (RS1.GT.0.0E0) GO TO 120
         NZ = N
         DO 140 I=1,N
            Y(I) = CZERO
  140    CONTINUE
         RETURN
      END
      SUBROUTINE CUNI2(Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUNI2
C***REFER TO  CBESI,CBESK
C
C     CUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  CAIRY,CUCHK,CUNHJ,CUOIK,R1MACH
C***END PROLOGUE  CUNI2
         COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CID, CIP, CONE, CRSC, CSCL,
     *    CSR, CSS, CY, CZERO, C1, C2, DAI, PHI, RZ, S1, S2, Y, Z, ZB,
     *    ZETA1, ZETA2, ZN, ZAR
         REAL AARG, AIC, ALIM, ANG, APHI, ASCLE, AY, BRY, CAR, C2I, C2M,
     *    C2R, ELIM, FN, FNU, FNUL, HPI, RS1, SAR, TOL, YY, R1MACH
         INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     *    NN, NUF, NW, NZ, IDUM
         DIMENSION BRY(3), Y(N), CIP(4), CSS(3), CSR(3), CY(2)
         DATA CZERO,CONE,CI/(0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0)/
         DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     1    (1.0E0,0.0E0), (0.0E0,1.0E0), (-1.0E0,0.0E0), (0.0E0,-1.0E0)/
         DATA HPI, AIC  /
     1         1.57079632679489662E+00,     1.265512123484645396E+00/
C
         NZ = 0
         ND = N
         NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         CRSC = CMPLX(TOL,0.0E0)
         CSS(1) = CSCL
         CSS(2) = CONE
         CSS(3) = CRSC
         CSR(1) = CRSC
         CSR(2) = CONE
         CSR(3) = CSCL
         BRY(1) = 1.0E+3*R1MACH(1)/TOL
         YY = AIMAG(Z)
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
         ZN = -Z*CI
         ZB = Z
         CID = -CI
         INU = INT(FNU)
         ANG = HPI*(FNU-FLOAT(INU))
         CAR = COS(ANG)
         SAR = SIN(ANG)
         C2 = CMPLX(CAR,SAR)
         ZAR = C2
         IN = INU + N - 1
         IN = MOD(IN,4)
         C2 = C2*CIP(IN+1)
         IF (YY.GT.0.0E0) GO TO 10
         ZN = CONJG(-ZN)
         ZB = CONJG(ZB)
         CID = -CID
         C2 = CONJG(C2)
   10    CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
         FN = AMAX1(FNU,1.0E0)
         CALL CUNHJ(ZN, FN, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
         IF (KODE.EQ.1) GO TO 20
         CFN = CMPLX(FNU,0.0E0)
         S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2))
         GO TO 30
   20    CONTINUE
         S1 = -ZETA1 + ZETA2
   30    CONTINUE
         RS1 = REAL(S1)
         IF (ABS(RS1).GT.ELIM) GO TO 150
   40    CONTINUE
         NN = MIN0(2,ND)
         DO 90 I=1,NN
            FN = FNU + FLOAT(ND-I)
            CALL CUNHJ(ZN, FN, 0, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
            IF (KODE.EQ.1) GO TO 50
            CFN = CMPLX(FN,0.0E0)
            AY = ABS(YY)
            S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2)) + CMPLX(0.0E0,AY)
            GO TO 60
   50       CONTINUE
            S1 = -ZETA1 + ZETA2
   60       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = REAL(S1)
            IF (ABS(RS1).GT.ELIM) GO TO 120
            IF (I.EQ.1) IFLAG = 2
            IF (ABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
            APHI = CABS(PHI)
            AARG = CABS(ARG)
            RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
            IF (ABS(RS1).GT.ELIM) GO TO 120
            IF (I.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0E0) GO TO 70
            IF (I.EQ.1) IFLAG = 3
   70       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
            CALL CAIRY(ARG, 0, 2, AI, NAI, IDUM)
            CALL CAIRY(ARG, 1, 2, DAI, NDAI, IDUM)
            S2 = PHI*(AI*ASUM+DAI*BSUM)
            C2R = REAL(S1)
            C2I = AIMAG(S1)
            C2M = EXP(C2R)*REAL(CSS(IFLAG))
            S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (IFLAG.NE.1) GO TO 80
            CALL CUCHK(S2, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 120
   80       CONTINUE
            IF (YY.LE.0.0E0) S2 = CONJG(S2)
            J = ND - I + 1
            S2 = S2*C2
            CY(I) = S2
            Y(J) = S2*CSR(IFLAG)
            C2 = C2*CID
   90    CONTINUE
         IF (ND.LE.2) GO TO 110
         RZ = CMPLX(2.0E0,0.0E0)/Z
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = R1MACH(2)
         S1 = CY(1)
         S2 = CY(2)
         C1 = CSR(IFLAG)
         ASCLE = BRY(IFLAG)
         K = ND - 2
         FN = FLOAT(K)
         DO 100 I=3,ND
            C2 = S2
            S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
            S1 = C2
            C2 = S2*C1
            Y(K) = C2
            K = K - 1
            FN = FN - 1.0E0
            IF (IFLAG.GE.3) GO TO 100
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = AMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 100
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1 = S1*C1
            S2 = C2
            S1 = S1*CSS(IFLAG)
            S2 = S2*CSS(IFLAG)
            C1 = CSR(IFLAG)
  100    CONTINUE
  110    CONTINUE
         RETURN
  120    CONTINUE
         IF (RS1.GT.0.0E0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
         Y(ND) = CZERO
         NZ = NZ + 1
         ND = ND - 1
         IF (ND.EQ.0) GO TO 110
         CALL CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
         IF (NUF.LT.0) GO TO 140
         ND = ND - NUF
         NZ = NZ + NUF
         IF (ND.EQ.0) GO TO 110
         FN = FNU + FLOAT(ND-1)
         IF (FN.LT.FNUL) GO TO 130
C      FN = AIMAG(CID)
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1 = CIP(K)
C      IF (FN.LT.0.0E0) S1 = CONJG(S1)
C      C2 = C2*S1
         IN = INU + ND - 1
         IN = MOD(IN,4) + 1
         C2 = ZAR*CIP(IN)
         IF (YY.LE.0.0E0)C2=CONJG(C2)
         GO TO 40
  130    CONTINUE
         NLAST = ND
         RETURN
  140    CONTINUE
         NZ = -1
         RETURN
  150    CONTINUE
         IF (RS1.GT.0.0E0) GO TO 140
         NZ = N
         DO 160 I=1,N
            Y(I) = CZERO
  160    CONTINUE
         RETURN
      END
      SUBROUTINE CUNIK(ZR, FNU, IKFLG, IPMTR, TOL, INIT, PHI, ZETA1,
     * ZETA2, SUM, CWRK)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUNIK
C***REFER TO  CBESI,CBESK
C
C        CUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUNIK
         COMPLEX CFN, CON, CONE, CRFN, CWRK, CZERO, PHI, S, SR, SUM, T,
     *    T2, ZETA1, ZETA2, ZN, ZR
         REAL AC, C, FNU, RFN, TEST, TOL, TSTR, TSTI
         INTEGER I, IKFLG, INIT, IPMTR, J, K, L
         DIMENSION C(120), CWRK(16), CON(2)
         DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
         DATA CON(1), CON(2)  /
     1   (3.98942280401432678E-01,0.0E0),(1.25331413731550025E+00,0.0E0)/
         DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1        C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2        C(19), C(20), C(21), C(22), C(23), C(24)/
     3        1.00000000000000000E+00,    -2.08333333333333333E-01,
     4        1.25000000000000000E-01,     3.34201388888888889E-01,
     5       -4.01041666666666667E-01,     7.03125000000000000E-02,
     6       -1.02581259645061728E+00,     1.84646267361111111E+00,
     7       -8.91210937500000000E-01,     7.32421875000000000E-02,
     8        4.66958442342624743E+00,    -1.12070026162229938E+01,
     9        8.78912353515625000E+00,    -2.36408691406250000E+00,
     A        1.12152099609375000E-01,    -2.82120725582002449E+01,
     B        8.46362176746007346E+01,    -9.18182415432400174E+01,
     C        4.25349987453884549E+01,    -7.36879435947963170E+00,
     D        2.27108001708984375E-01,     2.12570130039217123E+02,
     E       -7.65252468141181642E+02,     1.05999045252799988E+03/
         DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1        C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2        C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3       -6.99579627376132541E+02,     2.18190511744211590E+02,
     4       -2.64914304869515555E+01,     5.72501420974731445E-01,
     5       -1.91945766231840700E+03,     8.06172218173730938E+03,
     6       -1.35865500064341374E+04,     1.16553933368645332E+04,
     7       -5.30564697861340311E+03,     1.20090291321635246E+03,
     8       -1.08090919788394656E+02,     1.72772750258445740E+00,
     9        2.02042913309661486E+04,    -9.69805983886375135E+04,
     A        1.92547001232531532E+05,    -2.03400177280415534E+05,
     B        1.22200464983017460E+05,    -4.11926549688975513E+04,
     C        7.10951430248936372E+03,    -4.93915304773088012E+02,
     D        6.07404200127348304E+00,    -2.42919187900551333E+05,
     E        1.31176361466297720E+06,    -2.99801591853810675E+06/
         DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1        C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2        C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3        3.76327129765640400E+06,    -2.81356322658653411E+06,
     4        1.26836527332162478E+06,    -3.31645172484563578E+05,
     5        4.52187689813627263E+04,    -2.49983048181120962E+03,
     6        2.43805296995560639E+01,     3.28446985307203782E+06,
     7       -1.97068191184322269E+07,     5.09526024926646422E+07,
     8       -7.41051482115326577E+07,     6.63445122747290267E+07,
     9       -3.75671766607633513E+07,     1.32887671664218183E+07,
     A       -2.78561812808645469E+06,     3.08186404612662398E+05,
     B       -1.38860897537170405E+04,     1.10017140269246738E+02,
     C       -4.93292536645099620E+07,     3.25573074185765749E+08,
     D       -9.39462359681578403E+08,     1.55359689957058006E+09,
     E       -1.62108055210833708E+09,     1.10684281682301447E+09/
         DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1        C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2        C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3       -4.95889784275030309E+08,     1.42062907797533095E+08,
     4       -2.44740627257387285E+07,     2.24376817792244943E+06,
     5       -8.40054336030240853E+04,     5.51335896122020586E+02,
     6        8.14789096118312115E+08,    -5.86648149205184723E+09,
     7        1.86882075092958249E+10,    -3.46320433881587779E+10,
     8        4.12801855797539740E+10,    -3.30265997498007231E+10,
     9        1.79542137311556001E+10,    -6.56329379261928433E+09,
     A        1.55927986487925751E+09,    -2.25105661889415278E+08,
     B        1.73951075539781645E+07,    -5.49842327572288687E+05,
     C        3.03809051092238427E+03,    -1.46792612476956167E+10,
     D        1.14498237732025810E+11,    -3.99096175224466498E+11,
     E        8.19218669548577329E+11,    -1.09837515608122331E+12/
         DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1        C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2        C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3        1.00815810686538209E+12,    -6.45364869245376503E+11,
     4        2.87900649906150589E+11,    -8.78670721780232657E+10,
     5        1.76347306068349694E+10,    -2.16716498322379509E+09,
     6        1.43157876718888981E+08,    -3.87183344257261262E+06,
     7        1.82577554742931747E+04,     2.86464035717679043E+11,
     8       -2.40629790002850396E+12,     9.10934118523989896E+12,
     9       -2.05168994109344374E+13,     3.05651255199353206E+13,
     A       -3.16670885847851584E+13,     2.33483640445818409E+13,
     B       -1.23204913055982872E+13,     4.61272578084913197E+12,
     C       -1.19655288019618160E+12,     2.05914503232410016E+11,
     D       -2.18229277575292237E+10,     1.24700929351271032E+09/
         DATA C(119), C(120)/
     1       -2.91883881222208134E+07,     1.18838426256783253E+05/
C
         IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
         RFN = 1.0E0/FNU
         CRFN = CMPLX(RFN,0.0E0)
C     T = ZR*CRFN
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
         TSTR = REAL(ZR)
         TSTI = AIMAG(ZR)
         TEST = R1MACH(1)*1.0E+3
         AC = FNU*TEST
         IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
         AC = 2.0E0*ABS(ALOG(TEST))+FNU
         ZETA1 = CMPLX(AC,0.0E0)
         ZETA2 = CMPLX(FNU,0.0E0)
         PHI=CONE
         RETURN
   15    CONTINUE
         T=ZR*CRFN
         S = CONE + T*T
         SR = CSQRT(S)
         CFN = CMPLX(FNU,0.0E0)
         ZN = (CONE+SR)/T
         ZETA1 = CFN*CLOG(ZN)
         ZETA2 = CFN*SR
         T = CONE/SR
         SR = T*CRFN
         CWRK(16) = CSQRT(SR)
         PHI = CWRK(16)*CON(IKFLG)
         IF (IPMTR.NE.0) RETURN
         T2 = CONE/S
         CWRK(1) = CONE
         CRFN = CONE
         AC = 1.0E0
         L = 1
         DO 20 K=2,15
            S = CZERO
            DO 10 J=1,K
               L = L + 1
               S = S*T2 + CMPLX(C(L),0.0E0)
   10       CONTINUE
            CRFN = CRFN*SR
            CWRK(K) = CRFN*S
            AC = AC*RFN
            TSTR = REAL(CWRK(K))
            TSTI = AIMAG(CWRK(K))
            TEST = ABS(TSTR) + ABS(TSTI)
            IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20    CONTINUE
         K = 15
   30    CONTINUE
         INIT = K
   40    CONTINUE
         IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
         S = CZERO
         DO 50 I=1,INIT
            S = S + CWRK(I)
   50    CONTINUE
         SUM = S
         PHI = CWRK(16)*CON(1)
         RETURN
   60    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
         S = CZERO
         T = CONE
         DO 70 I=1,INIT
            S = S + T*CWRK(I)
            T = -T
   70    CONTINUE
         SUM = S
         PHI = CWRK(16)*CON(2)
         RETURN
      END
      SUBROUTINE CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUNK1
C***REFER TO  CBESK
C
C     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  CS1S2,CUCHK,CUNIK,R1MACH
C***END PROLOGUE  CUNK1
         COMPLEX CFN, CK, CONE, CRSC, CS, CSCL, CSGN, CSPN, CSR, CSS,
     *    CWRK, CY, CZERO, C1, C2, PHI,  RZ, SUM,  S1, S2, Y, Z,
     *    ZETA1,  ZETA2,  ZR, PHID, ZETA1D, ZETA2D, SUMD
         REAL ALIM, ANG, APHI, ASC, ASCLE, BRY, CPN, C2I, C2M, C2R, ELIM,
     *    FMR, FN, FNF, FNU, PI, RS1, SGN, SPN, TOL, X, R1MACH
         INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     *    KK, KODE, MR, N, NW, NZ, J, IPARD, INITD, IC
         DIMENSION BRY(3), INIT(2), Y(N), SUM(2), PHI(2), ZETA1(2),
     *    ZETA2(2), CY(2), CWRK(16,3), CSS(3), CSR(3)
         DATA CZERO, CONE / (0.0E0,0.0E0) , (1.0E0,0.0E0) /
         DATA PI / 3.14159265358979324E0 /
C
         KDFLG = 1
         NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         CRSC = CMPLX(TOL,0.0E0)
         CSS(1) = CSCL
         CSS(2) = CONE
         CSS(3) = CRSC
         CSR(1) = CRSC
         CSR(2) = CONE
         CSR(3) = CSCL
         BRY(1) = 1.0E+3*R1MACH(1)/TOL
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = R1MACH(2)
         X = REAL(Z)
         ZR = Z
         IF (X.LT.0.0E0) ZR = -Z
         J=2
         DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
            J = 3 - J
            FN = FNU + FLOAT(I-1)
            INIT(J) = 0
            CALL CUNIK(ZR, FN, 2, 0, TOL, INIT(J), PHI(J), ZETA1(J),
     *       ZETA2(J), SUM(J), CWRK(1,J))
            IF (KODE.EQ.1) GO TO 20
            CFN = CMPLX(FN,0.0E0)
            S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
            GO TO 30
   20       CONTINUE
            S1 = ZETA1(J) - ZETA2(J)
   30       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = REAL(S1)
            IF (ABS(RS1).GT.ELIM) GO TO 60
            IF (KDFLG.EQ.1) KFLAG = 2
            IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = CABS(PHI(J))
            RS1 = RS1 + ALOG(APHI)
            IF (ABS(RS1).GT.ELIM) GO TO 60
            IF (KDFLG.EQ.1) KFLAG = 1
            IF (RS1.LT.0.0E0) GO TO 40
            IF (KDFLG.EQ.1) KFLAG = 3
   40       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
            S2 = PHI(J)*SUM(J)
            C2R = REAL(S1)
            C2I = AIMAG(S1)
            C2M = EXP(C2R)*REAL(CSS(KFLAG))
            S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (KFLAG.NE.1) GO TO 50
            CALL CUCHK(S2, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 60
   50       CONTINUE
            CY(KDFLG) = S2
            Y(I) = S2*CSR(KFLAG)
            IF (KDFLG.EQ.2) GO TO 75
            KDFLG = 2
            GO TO 70
   60       CONTINUE
            IF (RS1.GT.0.0E0) GO TO 290
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
            IF (X.LT.0.0E0) GO TO 290
            KDFLG = 1
            Y(I) = CZERO
            NZ=NZ+1
            IF (I.EQ.1) GO TO 70
            IF (Y(I-1).EQ.CZERO) GO TO 70
            Y(I-1) = CZERO
            NZ=NZ+1
   70    CONTINUE
         I=N
   75    CONTINUE
         RZ = CMPLX(2.0E0,0.0E0)/ZR
         CK = CMPLX(FN,0.0E0)*RZ
         IB = I+1
         IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
C     ON UNDERFLOW
C-----------------------------------------------------------------------
         FN = FNU+FLOAT(N-1)
         IPARD = 1
         IF (MR.NE.0) IPARD = 0
         INITD = 0
         CALL CUNIK(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD,
     *   CWRK(1,3))
         IF (KODE.EQ.1) GO TO 80
         CFN=CMPLX(FN,0.0E0)
         S1=ZETA1D-CFN*(CFN/(ZR+ZETA2D))
         GO TO 90
   80    CONTINUE
         S1=ZETA1D-ZETA2D
   90    CONTINUE
         RS1=REAL(S1)
         IF (ABS(RS1).GT.ELIM) GO TO 95
         IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
         APHI=CABS(PHID)
         RS1=RS1+ALOG(APHI)
         IF (ABS(RS1).LT.ELIM) GO TO 100
   95    CONTINUE
         IF (RS1.GT.0.0E0) GO TO 290
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
         IF (X.LT.0.0E0) GO TO 290
         NZ=N
         DO 96 I=1,N
            Y(I) = CZERO
   96    CONTINUE
         RETURN
  100    CONTINUE
C-----------------------------------------------------------------------
C     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
         S1 = CY(1)
         S2 = CY(2)
         C1 = CSR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 120 I=IB,N
            C2 = S2
            S2 = CK*S2 + S1
            S1 = C2
            CK = CK + RZ
            C2 = S2*C1
            Y(I) = C2
            IF (KFLAG.GE.3) GO TO 120
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = AMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 120
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1 = S1*C1
            S2 = C2
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            C1 = CSR(KFLAG)
  120    CONTINUE
  160    CONTINUE
         IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C-----------------------------------------------------------------------
         NZ = 0
         FMR = FLOAT(MR)
         SGN = -SIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
C-----------------------------------------------------------------------
         CSGN = CMPLX(0.0E0,SGN)
         INU = INT(FNU)
         FNF = FNU - FLOAT(INU)
         IFN = INU + N - 1
         ANG = FNF*SGN
         CPN = COS(ANG)
         SPN = SIN(ANG)
         CSPN = CMPLX(CPN,SPN)
         IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
         ASC = BRY(1)
         KK = N
         IUF = 0
         KDFLG = 1
         IB = IB-1
         IC = IB-1
         DO 260 K=1,N
            FN = FNU + FLOAT(KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
            M=3
            IF (N.GT.2) GO TO 175
  170       CONTINUE
            INITD = INIT(J)
            PHID = PHI(J)
            ZETA1D = ZETA1(J)
            ZETA2D = ZETA2(J)
            SUMD = SUM(J)
            M = J
            J = 3 - J
            GO TO 180
  175       CONTINUE
            IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
            IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 170
            INITD = 0
  180       CONTINUE
            CALL CUNIK(ZR, FN, 1, 0, TOL, INITD, PHID, ZETA1D,
     *       ZETA2D, SUMD, CWRK(1,M))
            IF (KODE.EQ.1) GO TO 190
            CFN = CMPLX(FN,0.0E0)
            S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
            GO TO 200
  190       CONTINUE
            S1 = -ZETA1D + ZETA2D
  200       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = REAL(S1)
            IF (ABS(RS1).GT.ELIM) GO TO 250
            IF (KDFLG.EQ.1) IFLAG = 2
            IF (ABS(RS1).LT.ALIM) GO TO 210
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = CABS(PHID)
            RS1 = RS1 + ALOG(APHI)
            IF (ABS(RS1).GT.ELIM) GO TO 250
            IF (KDFLG.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0E0) GO TO 210
            IF (KDFLG.EQ.1) IFLAG = 3
  210       CONTINUE
            S2 = CSGN*PHID*SUMD
            C2R = REAL(S1)
            C2I = AIMAG(S1)
            C2M = EXP(C2R)*REAL(CSS(IFLAG))
            S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (IFLAG.NE.1) GO TO 220
            CALL CUCHK(S2, NW, BRY(1), TOL)
            IF (NW.NE.0) S2 = CMPLX(0.0E0,0.0E0)
  220       CONTINUE
            CY(KDFLG) = S2
            C2 = S2
            S2 = S2*CSR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
            S1 = Y(KK)
            IF (KODE.EQ.1) GO TO 240
            CALL CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  240       CONTINUE
            Y(KK) = S1*CSPN + S2
            KK = KK - 1
            CSPN = -CSPN
            IF (C2.NE.CZERO) GO TO 245
            KDFLG = 1
            GO TO 260
  245       CONTINUE
            IF (KDFLG.EQ.2) GO TO 265
            KDFLG = 2
            GO TO 260
  250       CONTINUE
            IF (RS1.GT.0.0E0) GO TO 290
            S2 = CZERO
            GO TO 220
  260    CONTINUE
         K = N
  265    CONTINUE
         IL = N - K
         IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
         S1 = CY(1)
         S2 = CY(2)
         CS = CSR(IFLAG)
         ASCLE = BRY(IFLAG)
         FN = FLOAT(INU+IL)
         DO 280 I=1,IL
            C2 = S2
            S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
            S1 = C2
            FN = FN - 1.0E0
            C2 = S2*CS
            CK = C2
            C1 = Y(KK)
            IF (KODE.EQ.1) GO TO 270
            CALL CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  270       CONTINUE
            Y(KK) = C1*CSPN + C2
            KK = KK - 1
            CSPN = -CSPN
            IF (IFLAG.GE.3) GO TO 280
            C2R = REAL(CK)
            C2I = AIMAG(CK)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = AMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 280
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1 = S1*CS
            S2 = CK
            S1 = S1*CSS(IFLAG)
            S2 = S2*CSS(IFLAG)
            CS = CSR(IFLAG)
  280    CONTINUE
         RETURN
  290    CONTINUE
         NZ = -1
         RETURN
      END
      SUBROUTINE CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUNK2
C***REFER TO  CBESK
C
C     CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  CAIRY,CS1S2,CUCHK,CUNHJ,R1MACH
C***END PROLOGUE  CUNK2
         COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CIP,
     *    CK, CONE, CRSC, CR1, CR2, CS, CSCL, CSGN, CSPN, CSR, CSS, CY,
     *    CZERO, C1, C2, DAI, PHI,  RZ, S1, S2, Y, Z, ZB, ZETA1,
     *    ZETA2, ZN, ZR, PHID, ARGD, ZETA1D, ZETA2D, ASUMD, BSUMD
         REAL AARG, AIC, ALIM, ANG, APHI, ASC, ASCLE, BRY, CAR, CPN, C2I,
     *    C2M, C2R, ELIM, FMR, FN, FNF, FNU, HPI, PI, RS1, SAR, SGN, SPN,
     *    TOL, X, YY, R1MACH
         INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     *    KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
         DIMENSION BRY(3), Y(N), ASUM(2), BSUM(2), PHI(2), ARG(2),
     *    ZETA1(2), ZETA2(2), CY(2), CIP(4), CSS(3), CSR(3)
         DATA CZERO, CONE, CI, CR1, CR2 /
     1            (0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0),
     1   (1.0E0,1.73205080756887729E0),(-0.5E0,-8.66025403784438647E-01)/
         DATA HPI, PI, AIC /
     1        1.57079632679489662E+00,     3.14159265358979324E+00,
     1        1.26551212348464539E+00/
         DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     1    (1.0E0,0.0E0), (0.0E0,-1.0E0), (-1.0E0,0.0E0), (0.0E0,1.0E0)/
C
         KDFLG = 1
         NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         CRSC = CMPLX(TOL,0.0E0)
         CSS(1) = CSCL
         CSS(2) = CONE
         CSS(3) = CRSC
         CSR(1) = CRSC
         CSR(2) = CONE
         CSR(3) = CSCL
         BRY(1) = 1.0E+3*R1MACH(1)/TOL
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = R1MACH(2)
         X = REAL(Z)
         ZR = Z
         IF (X.LT.0.0E0) ZR = -Z
         YY = AIMAG(ZR)
         ZN = -ZR*CI
         ZB = ZR
         INU = INT(FNU)
         FNF = FNU - FLOAT(INU)
         ANG = -HPI*FNF
         CAR = COS(ANG)
         SAR = SIN(ANG)
         CPN = -HPI*CAR
         SPN = -HPI*SAR
         C2 = CMPLX(-SPN,CPN)
         KK = MOD(INU,4) + 1
         CS = CR1*C2*CIP(KK)
         IF (YY.GT.0.0E0) GO TO 10
         ZN = CONJG(-ZN)
         ZB = CONJG(ZB)
   10    CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
         J = 2
         DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
            J = 3 - J
            FN = FNU + FLOAT(I-1)
            CALL CUNHJ(ZN, FN, 0, TOL, PHI(J), ARG(J), ZETA1(J), ZETA2(J),
     *       ASUM(J), BSUM(J))
            IF (KODE.EQ.1) GO TO 20
            CFN = CMPLX(FN,0.0E0)
            S1 = ZETA1(J) - CFN*(CFN/(ZB+ZETA2(J)))
            GO TO 30
   20       CONTINUE
            S1 = ZETA1(J) - ZETA2(J)
   30       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = REAL(S1)
            IF (ABS(RS1).GT.ELIM) GO TO 60
            IF (KDFLG.EQ.1) KFLAG = 2
            IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = CABS(PHI(J))
            AARG = CABS(ARG(J))
            RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
            IF (ABS(RS1).GT.ELIM) GO TO 60
            IF (KDFLG.EQ.1) KFLAG = 1
            IF (RS1.LT.0.0E0) GO TO 40
            IF (KDFLG.EQ.1) KFLAG = 3
   40       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
            C2 = ARG(J)*CR2
            CALL CAIRY(C2, 0, 2, AI, NAI, IDUM)
            CALL CAIRY(C2, 1, 2, DAI, NDAI, IDUM)
            S2 = CS*PHI(J)*(AI*ASUM(J)+CR2*DAI*BSUM(J))
            C2R = REAL(S1)
            C2I = AIMAG(S1)
            C2M = EXP(C2R)*REAL(CSS(KFLAG))
            S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (KFLAG.NE.1) GO TO 50
            CALL CUCHK(S2, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 60
   50       CONTINUE
            IF (YY.LE.0.0E0) S2 = CONJG(S2)
            CY(KDFLG) = S2
            Y(I) = S2*CSR(KFLAG)
            CS = -CI*CS
            IF (KDFLG.EQ.2) GO TO 75
            KDFLG = 2
            GO TO 70
   60       CONTINUE
            IF (RS1.GT.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
            IF (X.LT.0.0E0) GO TO 300
            KDFLG = 1
            Y(I) = CZERO
            CS = -CI*CS
            NZ=NZ+1
            IF (I.EQ.1) GO TO 70
            IF (Y(I-1).EQ.CZERO) GO TO 70
            Y(I-1) = CZERO
            NZ=NZ+1
   70    CONTINUE
         I=N
   75    CONTINUE
         RZ = CMPLX(2.0E0,0.0E0)/ZR
         CK = CMPLX(FN,0.0E0)*RZ
         IB = I + 1
         IF (N.LT.IB) GO TO 170
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
C     ON UNDERFLOW
C-----------------------------------------------------------------------
         FN = FNU+FLOAT(N-1)
         IPARD = 1
         IF (MR.NE.0) IPARD = 0
         CALL CUNHJ(ZN,FN,IPARD,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD,BSUMD)
         IF (KODE.EQ.1) GO TO 80
         CFN=CMPLX(FN,0.0E0)
         S1=ZETA1D-CFN*(CFN/(ZB+ZETA2D))
         GO TO 90
   80    CONTINUE
         S1=ZETA1D-ZETA2D
   90    CONTINUE
         RS1=REAL(S1)
         IF (ABS(RS1).GT.ELIM) GO TO 95
         IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
         APHI=CABS(PHID)
         AARG = CABS(ARGD)
         RS1=RS1+ALOG(APHI)-0.25E0*ALOG(AARG)-AIC
         IF (ABS(RS1).LT.ELIM) GO TO 100
   95    CONTINUE
         IF (RS1.GT.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
         IF (X.LT.0.0E0) GO TO 300
         NZ=N
         DO 96 I=1,N
            Y(I) = CZERO
   96    CONTINUE
         RETURN
  100    CONTINUE
C-----------------------------------------------------------------------
C     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
         S1 = CY(1)
         S2 = CY(2)
         C1 = CSR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 120 I=IB,N
            C2 = S2
            S2 = CK*S2 + S1
            S1 = C2
            CK = CK + RZ
            C2 = S2*C1
            Y(I) = C2
            IF (KFLAG.GE.3) GO TO 120
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = AMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 120
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1 = S1*C1
            S2 = C2
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            C1 = CSR(KFLAG)
  120    CONTINUE
  170    CONTINUE
         IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C-----------------------------------------------------------------------
         NZ = 0
         FMR = FLOAT(MR)
         SGN = -SIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
         CSGN = CMPLX(0.0E0,SGN)
         IF (YY.LE.0.0E0) CSGN = CONJG(CSGN)
         IFN = INU + N - 1
         ANG = FNF*SGN
         CPN = COS(ANG)
         SPN = SIN(ANG)
         CSPN = CMPLX(CPN,SPN)
         IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
         CS = CMPLX(CAR,-SAR)*CSGN
         IN = MOD(IFN,4) + 1
         C2 = CIP(IN)
         CS = CS*CONJG(C2)
         ASC = BRY(1)
         KK = N
         KDFLG = 1
         IB = IB-1
         IC = IB-1
         IUF = 0
         DO 270 K=1,N
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
            FN = FNU+FLOAT(KK-1)
            IF (N.GT.2) GO TO 180
  175       CONTINUE
            PHID = PHI(J)
            ARGD = ARG(J)
            ZETA1D = ZETA1(J)
            ZETA2D = ZETA2(J)
            ASUMD = ASUM(J)
            BSUMD = BSUM(J)
            J = 3 - J
            GO TO 190
  180       CONTINUE
            IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 190
            IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 175
            CALL CUNHJ(ZN, FN, 0, TOL, PHID, ARGD, ZETA1D, ZETA2D,
     *       ASUMD, BSUMD)
  190       CONTINUE
            IF (KODE.EQ.1) GO TO 200
            CFN = CMPLX(FN,0.0E0)
            S1 = -ZETA1D + CFN*(CFN/(ZB+ZETA2D))
            GO TO 210
  200       CONTINUE
            S1 = -ZETA1D + ZETA2D
  210       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = REAL(S1)
            IF (ABS(RS1).GT.ELIM) GO TO 260
            IF (KDFLG.EQ.1) IFLAG = 2
            IF (ABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = CABS(PHID)
            AARG = CABS(ARGD)
            RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
            IF (ABS(RS1).GT.ELIM) GO TO 260
            IF (KDFLG.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0E0) GO TO 220
            IF (KDFLG.EQ.1) IFLAG = 3
  220       CONTINUE
            CALL CAIRY(ARGD, 0, 2, AI, NAI, IDUM)
            CALL CAIRY(ARGD, 1, 2, DAI, NDAI, IDUM)
            S2 = CS*PHID*(AI*ASUMD+DAI*BSUMD)
            C2R = REAL(S1)
            C2I = AIMAG(S1)
            C2M = EXP(C2R)*REAL(CSS(IFLAG))
            S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (IFLAG.NE.1) GO TO 230
            CALL CUCHK(S2, NW, BRY(1), TOL)
            IF (NW.NE.0) S2 = CMPLX(0.0E0,0.0E0)
  230       CONTINUE
            IF (YY.LE.0.0E0) S2 = CONJG(S2)
            CY(KDFLG) = S2
            C2 = S2
            S2 = S2*CSR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
            S1 = Y(KK)
            IF (KODE.EQ.1) GO TO 250
            CALL CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  250       CONTINUE
            Y(KK) = S1*CSPN + S2
            KK = KK - 1
            CSPN = -CSPN
            CS = -CS*CI
            IF (C2.NE.CZERO) GO TO 255
            KDFLG = 1
            GO TO 270
  255       CONTINUE
            IF (KDFLG.EQ.2) GO TO 275
            KDFLG = 2
            GO TO 270
  260       CONTINUE
            IF (RS1.GT.0.0E0) GO TO 300
            S2 = CZERO
            GO TO 230
  270    CONTINUE
         K = N
  275    CONTINUE
         IL = N-K
         IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
         S1 = CY(1)
         S2 = CY(2)
         CS = CSR(IFLAG)
         ASCLE = BRY(IFLAG)
         FN = FLOAT(INU+IL)
         DO 290 I=1,IL
            C2 = S2
            S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
            S1 = C2
            FN = FN - 1.0E0
            C2 = S2*CS
            CK = C2
            C1 = Y(KK)
            IF (KODE.EQ.1) GO TO 280
            CALL CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  280       CONTINUE
            Y(KK) = C1*CSPN + C2
            KK = KK - 1
            CSPN = -CSPN
            IF (IFLAG.GE.3) GO TO 290
            C2R = REAL(CK)
            C2I = AIMAG(CK)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = AMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 290
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1 = S1*CS
            S2 = CK
            S1 = S1*CSS(IFLAG)
            S2 = S2*CSS(IFLAG)
            CS = CSR(IFLAG)
  290    CONTINUE
         RETURN
  300    CONTINUE
         NZ = -1
         RETURN
      END
      SUBROUTINE CUOIK(Z, FNU, KODE, IKFLG, N, Y, NUF, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CUOIK
C***REFER TO  CBESI,CBESK,CBESH
C
C     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***ROUTINES CALLED  CUCHK,CUNHJ,CUNIK,R1MACH
C***END PROLOGUE  CUOIK
         COMPLEX ARG, ASUM, BSUM, CWRK, CZ, CZERO, PHI, SUM, Y, Z, ZB,
     *    ZETA1, ZETA2, ZN, ZR
         REAL AARG, AIC, ALIM, APHI, ASCLE, AX, AY, ELIM, FNN, FNU, GNN,
     *    GNU, RCZ, TOL, X, YY
         INTEGER I, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
         DIMENSION Y(N), CWRK(16)
         DATA CZERO / (0.0E0,0.0E0) /
         DATA AIC / 1.265512123484645396E+00 /
         NUF = 0
         NN = N
         X = REAL(Z)
         ZR = Z
         IF (X.LT.0.0E0) ZR = -Z
         ZB = ZR
         YY = AIMAG(ZR)
         AX = ABS(X)*1.7321E0
         AY = ABS(YY)
         IFORM = 1
         IF (AY.GT.AX) IFORM = 2
         GNU = AMAX1(FNU,1.0E0)
         IF (IKFLG.EQ.1) GO TO 10
         FNN = FLOAT(NN)
         GNN = FNU + FNN - 1.0E0
         GNU = AMAX1(GNN,FNN)
   10    CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
         IF (IFORM.EQ.2) GO TO 20
         INIT = 0
         CALL CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM,
     *    CWRK)
         CZ = -ZETA1 + ZETA2
         GO TO 40
   20    CONTINUE
         ZN = -ZR*CMPLX(0.0E0,1.0E0)
         IF (YY.GT.0.0E0) GO TO 30
         ZN = CONJG(-ZN)
   30    CONTINUE
         CALL CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
         CZ = -ZETA1 + ZETA2
         AARG = CABS(ARG)
   40    CONTINUE
         IF (KODE.EQ.2) CZ = CZ - ZB
         IF (IKFLG.EQ.2) CZ = -CZ
         APHI = CABS(PHI)
         RCZ = REAL(CZ)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         IF (RCZ.GT.ELIM) GO TO 170
         IF (RCZ.LT.ALIM) GO TO 50
         RCZ = RCZ + ALOG(APHI)
         IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
         IF (RCZ.GT.ELIM) GO TO 170
         GO TO 100
   50    CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
         IF (RCZ.LT.(-ELIM)) GO TO 60
         IF (RCZ.GT.(-ALIM)) GO TO 100
         RCZ = RCZ + ALOG(APHI)
         IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
         IF (RCZ.GT.(-ELIM)) GO TO 80
   60    CONTINUE
         DO 70 I=1,NN
            Y(I) = CZERO
   70    CONTINUE
         NUF = NN
         RETURN
   80    CONTINUE
         ASCLE = 1.0E+3*R1MACH(1)/TOL
         CZ = CZ + CLOG(PHI)
         IF (IFORM.EQ.1) GO TO 90
         CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
   90    CONTINUE
         AX = EXP(RCZ)/TOL
         AY = AIMAG(CZ)
         CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
         CALL CUCHK(CZ, NW, ASCLE, TOL)
         IF (NW.EQ.1) GO TO 60
  100    CONTINUE
         IF (IKFLG.EQ.2) RETURN
         IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  110    CONTINUE
         GNU = FNU + FLOAT(NN-1)
         IF (IFORM.EQ.2) GO TO 120
         INIT = 0
         CALL CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM,
     *    CWRK)
         CZ = -ZETA1 + ZETA2
         GO TO 130
  120    CONTINUE
         CALL CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
         CZ = -ZETA1 + ZETA2
         AARG = CABS(ARG)
  130    CONTINUE
         IF (KODE.EQ.2) CZ = CZ - ZB
         APHI = CABS(PHI)
         RCZ = REAL(CZ)
         IF (RCZ.LT.(-ELIM)) GO TO 140
         IF (RCZ.GT.(-ALIM)) RETURN
         RCZ = RCZ + ALOG(APHI)
         IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
         IF (RCZ.GT.(-ELIM)) GO TO 150
  140    CONTINUE
         Y(NN) = CZERO
         NN = NN - 1
         NUF = NUF + 1
         IF (NN.EQ.0) RETURN
         GO TO 110
  150    CONTINUE
         ASCLE = 1.0E+3*R1MACH(1)/TOL
         CZ = CZ + CLOG(PHI)
         IF (IFORM.EQ.1) GO TO 160
         CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
  160    CONTINUE
         AX = EXP(RCZ)/TOL
         AY = AIMAG(CZ)
         CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
         CALL CUCHK(CZ, NW, ASCLE, TOL)
         IF (NW.EQ.1) GO TO 140
         RETURN
  170    CONTINUE
         NUF = -1
         RETURN
      END
      SUBROUTINE CWRSK(ZR, FNU, KODE, N, Y, NZ, CW, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  CWRSK
C***REFER TO  CBESI,CBESK
C
C     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN
C
C***ROUTINES CALLED  CBKNU,CRATI,R1MACH
C***END PROLOGUE  CWRSK
         COMPLEX CINU, CSCL, CT, CW, C1, C2, RCT, ST, Y, ZR
         REAL ACT, ACW, ALIM, ASCLE, ELIM, FNU, S1, S2, TOL, YY
         INTEGER I, KODE, N, NW, NZ
         DIMENSION Y(N), CW(2)
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
         NZ = 0
         CALL CBKNU(ZR, FNU, KODE, 2, CW, NW, TOL, ELIM, ALIM)
         IF (NW.NE.0) GO TO 50
         CALL CRATI(ZR, FNU, N, Y, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
         CINU = CMPLX(1.0E0,0.0E0)
         IF (KODE.EQ.1) GO TO 10
         YY = AIMAG(ZR)
         S1 = COS(YY)
         S2 = SIN(YY)
         CINU = CMPLX(S1,S2)
   10    CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
         ACW = CABS(CW(2))
         ASCLE = 1.0E+3*R1MACH(1)/TOL
         CSCL = CMPLX(1.0E0,0.0E0)
         IF (ACW.GT.ASCLE) GO TO 20
         CSCL = CMPLX(1.0E0/TOL,0.0E0)
         GO TO 30
   20    CONTINUE
         ASCLE = 1.0E0/ASCLE
         IF (ACW.LT.ASCLE) GO TO 30
         CSCL = CMPLX(TOL,0.0E0)
   30    CONTINUE
         C1 = CW(1)*CSCL
         C2 = CW(2)*CSCL
         ST = Y(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0E0/CABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
C-----------------------------------------------------------------------
         CT = ZR*(C2+ST*C1)
         ACT = CABS(CT)
         RCT = CMPLX(1.0E0/ACT,0.0E0)
         CT = CONJG(CT)*RCT
         CINU = CINU*RCT*CT
         Y(1) = CINU*CSCL
         IF (N.EQ.1) RETURN
         DO 40 I=2,N
            CINU = ST*CINU
            ST = Y(I)
            Y(I) = CINU*CSCL
   40    CONTINUE
         RETURN
   50    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c    D1MACH ( 3) = B**(-T), the smallest relative spacing.
c    D1MACH ( 4) = B**(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
         implicit none

         double precision d1mach
         integer i

         if ( i < 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'D1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
            write ( *, '(a,i12)' ) '  I = ', i
            d1mach = 0.0D+00
            stop
         else if ( i == 1 ) then
            d1mach = 4.450147717014403D-308
         else if ( i == 2 ) then
            d1mach = 8.988465674311579D+307
         else if ( i == 3 ) then
            d1mach = 1.110223024625157D-016
         else if ( i == 4 ) then
            d1mach = 2.220446049250313D-016
         else if ( i == 5 ) then
            d1mach = 0.301029995663981D+000
         else if ( 5 < i ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'D1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
            write ( *, '(a,i12)' ) '  I = ', i
            d1mach = 0.0D+00
            stop
         end if

         return
      end
      DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C               **** A DOUBLE PRECISION ROUTINE ****
C         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT      Z IS D0UBLE PRECISION
C           Z      - ARGUMENT, Z.GT.0.0D0
C
C         OUTPUT      DGAMLN IS DOUBLE PRECISION
C           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
C           IERR    - ERROR FLAG
C                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
C
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,D1MACH
C***END PROLOGUE  DGAMLN
         DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
     *    T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, D1MACH
         INTEGER I, IERR, I1M, K, MZ, NZ, I1MACH
         DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
         DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1        GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2        GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3        GLN(21), GLN(22)/
     4        0.00000000000000000D+00,     0.00000000000000000D+00,
     5        6.93147180559945309D-01,     1.79175946922805500D+00,
     6        3.17805383034794562D+00,     4.78749174278204599D+00,
     7        6.57925121201010100D+00,     8.52516136106541430D+00,
     8        1.06046029027452502D+01,     1.28018274800814696D+01,
     9        1.51044125730755153D+01,     1.75023078458738858D+01,
     A        1.99872144956618861D+01,     2.25521638531234229D+01,
     B        2.51912211827386815D+01,     2.78992713838408916D+01,
     C        3.06718601060806728D+01,     3.35050734501368889D+01,
     D        3.63954452080330536D+01,     3.93398841871994940D+01,
     E        4.23356164607534850D+01,     4.53801388984769080D+01/
         DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1        GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2        GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3        GLN(41), GLN(42), GLN(43), GLN(44)/
     4        4.84711813518352239D+01,     5.16066755677643736D+01,
     5        5.47847293981123192D+01,     5.80036052229805199D+01,
     6        6.12617017610020020D+01,     6.45575386270063311D+01,
     7        6.78897431371815350D+01,     7.12570389671680090D+01,
     8        7.46582363488301644D+01,     7.80922235533153106D+01,
     9        8.15579594561150372D+01,     8.50544670175815174D+01,
     A        8.85808275421976788D+01,     9.21361756036870925D+01,
     B        9.57196945421432025D+01,     9.93306124547874269D+01,
     C        1.02968198614513813D+02,     1.06631760260643459D+02,
     D        1.10320639714757395D+02,     1.14034211781461703D+02,
     E        1.17771881399745072D+02,     1.21533081515438634D+02/
         DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1        GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2        GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3        GLN(63), GLN(64), GLN(65), GLN(66)/
     4        1.25317271149356895D+02,     1.29123933639127215D+02,
     5        1.32952575035616310D+02,     1.36802722637326368D+02,
     6        1.40673923648234259D+02,     1.44565743946344886D+02,
     7        1.48477766951773032D+02,     1.52409592584497358D+02,
     8        1.56360836303078785D+02,     1.60331128216630907D+02,
     9        1.64320112263195181D+02,     1.68327445448427652D+02,
     A        1.72352797139162802D+02,     1.76395848406997352D+02,
     B        1.80456291417543771D+02,     1.84533828861449491D+02,
     C        1.88628173423671591D+02,     1.92739047287844902D+02,
     D        1.96866181672889994D+02,     2.01009316399281527D+02,
     E        2.05168199482641199D+02,     2.09342586752536836D+02/
         DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1        GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2        GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3        GLN(85), GLN(86), GLN(87), GLN(88)/
     4        2.13532241494563261D+02,     2.17736934113954227D+02,
     5        2.21956441819130334D+02,     2.26190548323727593D+02,
     6        2.30439043565776952D+02,     2.34701723442818268D+02,
     7        2.38978389561834323D+02,     2.43268849002982714D+02,
     8        2.47572914096186884D+02,     2.51890402209723194D+02,
     9        2.56221135550009525D+02,     2.60564940971863209D+02,
     A        2.64921649798552801D+02,     2.69291097651019823D+02,
     B        2.73673124285693704D+02,     2.78067573440366143D+02,
     C        2.82474292687630396D+02,     2.86893133295426994D+02,
     D        2.91323950094270308D+02,     2.95766601350760624D+02,
     E        3.00220948647014132D+02,     3.04686856765668715D+02/
         DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1        GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2        3.09164193580146922D+02,     3.13652829949879062D+02,
     3        3.18152639620209327D+02,     3.22663499126726177D+02,
     4        3.27185287703775217D+02,     3.31717887196928473D+02,
     5        3.36261181979198477D+02,     3.40815058870799018D+02,
     6        3.45379407062266854D+02,     3.49954118040770237D+02,
     7        3.54539085519440809D+02,     3.59134205369575399D+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
         DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1        CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2        CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3        8.33333333333333333D-02,    -2.77777777777777778D-03,
     4        7.93650793650793651D-04,    -5.95238095238095238D-04,
     5        8.41750841750841751D-04,    -1.91752691752691753D-03,
     6        6.41025641025641026D-03,    -2.95506535947712418D-02,
     7        1.79644372368830573D-01,    -1.39243221690590112D+00,
     8        1.34028640441683920D+01,    -1.56848284626002017D+02,
     9        2.19310333333333333D+03,    -3.61087712537249894D+04,
     A        6.91472268851313067D+05,    -1.52382215394074162D+07,
     B        3.82900751391414141D+08,    -1.08822660357843911D+10,
     C        3.47320283765002252D+11,    -1.23696021422692745D+13,
     D        4.88788064793079335D+14,    -2.13203339609193739D+16/
C
C             LN(2*PI)
         DATA CON                    /     1.83787706640934548D+00/
C
C***FIRST EXECUTABLE STATEMENT  DGAMLN
         IERR=0
         IF (Z.LE.0.0D0) GO TO 70
         IF (Z.GT.101.0D0) GO TO 10
         NZ = INT(Z)
         FZ = Z - FLOAT(NZ)
         IF (FZ.GT.0.0D0) GO TO 10
         IF (NZ.GT.100) GO TO 10
         DGAMLN = GLN(NZ)
         RETURN
   10    CONTINUE
         WDTOL = D1MACH(4)
         WDTOL = DMAX1(WDTOL,0.5D-18)
         I1M = I1MACH(14)
         RLN = D1MACH(5)*FLOAT(I1M)
         FLN = DMIN1(RLN,20.0D0)
         FLN = DMAX1(FLN,3.0D0)
         FLN = FLN - 3.0D0
         ZM = 1.8000D0 + 0.3875D0*FLN
         MZ = INT(SNGL(ZM)) + 1
         ZMIN = FLOAT(MZ)
         ZDMY = Z
         ZINC = 0.0D0
         IF (Z.GE.ZMIN) GO TO 20
         ZINC = ZMIN - FLOAT(NZ)
         ZDMY = Z + ZINC
   20    CONTINUE
         ZP = 1.0D0/ZDMY
         T1 = CF(1)*ZP
         S = T1
         IF (ZP.LT.WDTOL) GO TO 40
         ZSQ = ZP*ZP
         TST = T1*WDTOL
         DO 30 K=2,22
            ZP = ZP*ZSQ
            TRM = CF(K)*ZP
            IF (DABS(TRM).LT.TST) GO TO 40
            S = S + TRM
   30    CONTINUE
   40    CONTINUE
         IF (ZINC.NE.0.0D0) GO TO 50
         TLG = DLOG(Z)
         DGAMLN = Z*(TLG-1.0D0) + 0.5D0*(CON-TLG) + S
         RETURN
   50    CONTINUE
         ZP = 1.0D0
         NZ = INT(SNGL(ZINC))
         DO 60 I=1,NZ
            ZP = ZP*(Z+FLOAT(I-1))
   60    CONTINUE
         TLG = DLOG(ZDMY)
         DGAMLN = ZDMY*(TLG-1.0D0) - DLOG(ZP) + 0.5D0*(CON-TLG) + S
         RETURN
C
C
   70    CONTINUE
         IERR=1
         RETURN
      END
      FUNCTION GAMLN(Z,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  GAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT
C           Z      - REAL ARGUMENT, Z.GT.0.0E0
C
C         OUTPUT
C           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                    IERR=1, Z.LE.0.0E0,    NO COMPUTATION
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,R1MACH
C***END PROLOGUE  GAMLN
C
         INTEGER I, I1M, K, MZ, NZ, IERR, I1MACH
         REAL CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST, T1, WDTOL, Z,
     *    ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ
         REAL R1MACH
         DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
         DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1        GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2        GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3        GLN(21), GLN(22)/
     4        0.00000000000000000E+00,     0.00000000000000000E+00,
     5        6.93147180559945309E-01,     1.79175946922805500E+00,
     6        3.17805383034794562E+00,     4.78749174278204599E+00,
     7        6.57925121201010100E+00,     8.52516136106541430E+00,
     8        1.06046029027452502E+01,     1.28018274800814696E+01,
     9        1.51044125730755153E+01,     1.75023078458738858E+01,
     A        1.99872144956618861E+01,     2.25521638531234229E+01,
     B        2.51912211827386815E+01,     2.78992713838408916E+01,
     C        3.06718601060806728E+01,     3.35050734501368889E+01,
     D        3.63954452080330536E+01,     3.93398841871994940E+01,
     E        4.23356164607534850E+01,     4.53801388984769080E+01/
         DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1        GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2        GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3        GLN(41), GLN(42), GLN(43), GLN(44)/
     4        4.84711813518352239E+01,     5.16066755677643736E+01,
     5        5.47847293981123192E+01,     5.80036052229805199E+01,
     6        6.12617017610020020E+01,     6.45575386270063311E+01,
     7        6.78897431371815350E+01,     7.12570389671680090E+01,
     8        7.46582363488301644E+01,     7.80922235533153106E+01,
     9        8.15579594561150372E+01,     8.50544670175815174E+01,
     A        8.85808275421976788E+01,     9.21361756036870925E+01,
     B        9.57196945421432025E+01,     9.93306124547874269E+01,
     C        1.02968198614513813E+02,     1.06631760260643459E+02,
     D        1.10320639714757395E+02,     1.14034211781461703E+02,
     E        1.17771881399745072E+02,     1.21533081515438634E+02/
         DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1        GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2        GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3        GLN(63), GLN(64), GLN(65), GLN(66)/
     4        1.25317271149356895E+02,     1.29123933639127215E+02,
     5        1.32952575035616310E+02,     1.36802722637326368E+02,
     6        1.40673923648234259E+02,     1.44565743946344886E+02,
     7        1.48477766951773032E+02,     1.52409592584497358E+02,
     8        1.56360836303078785E+02,     1.60331128216630907E+02,
     9        1.64320112263195181E+02,     1.68327445448427652E+02,
     A        1.72352797139162802E+02,     1.76395848406997352E+02,
     B        1.80456291417543771E+02,     1.84533828861449491E+02,
     C        1.88628173423671591E+02,     1.92739047287844902E+02,
     D        1.96866181672889994E+02,     2.01009316399281527E+02,
     E        2.05168199482641199E+02,     2.09342586752536836E+02/
         DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1        GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2        GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3        GLN(85), GLN(86), GLN(87), GLN(88)/
     4        2.13532241494563261E+02,     2.17736934113954227E+02,
     5        2.21956441819130334E+02,     2.26190548323727593E+02,
     6        2.30439043565776952E+02,     2.34701723442818268E+02,
     7        2.38978389561834323E+02,     2.43268849002982714E+02,
     8        2.47572914096186884E+02,     2.51890402209723194E+02,
     9        2.56221135550009525E+02,     2.60564940971863209E+02,
     A        2.64921649798552801E+02,     2.69291097651019823E+02,
     B        2.73673124285693704E+02,     2.78067573440366143E+02,
     C        2.82474292687630396E+02,     2.86893133295426994E+02,
     D        2.91323950094270308E+02,     2.95766601350760624E+02,
     E        3.00220948647014132E+02,     3.04686856765668715E+02/
         DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1        GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2        3.09164193580146922E+02,     3.13652829949879062E+02,
     3        3.18152639620209327E+02,     3.22663499126726177E+02,
     4        3.27185287703775217E+02,     3.31717887196928473E+02,
     5        3.36261181979198477E+02,     3.40815058870799018E+02,
     6        3.45379407062266854E+02,     3.49954118040770237E+02,
     7        3.54539085519440809E+02,     3.59134205369575399E+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
         DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1        CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2        CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3        8.33333333333333333E-02,    -2.77777777777777778E-03,
     4        7.93650793650793651E-04,    -5.95238095238095238E-04,
     5        8.41750841750841751E-04,    -1.91752691752691753E-03,
     6        6.41025641025641026E-03,    -2.95506535947712418E-02,
     7        1.79644372368830573E-01,    -1.39243221690590112E+00,
     8        1.34028640441683920E+01,    -1.56848284626002017E+02,
     9        2.19310333333333333E+03,    -3.61087712537249894E+04,
     A        6.91472268851313067E+05,    -1.52382215394074162E+07,
     B        3.82900751391414141E+08,    -1.08822660357843911E+10,
     C        3.47320283765002252E+11,    -1.23696021422692745E+13,
     D        4.88788064793079335E+14,    -2.13203339609193739E+16/
C
C             LN(2*PI)
         DATA CON                    /     1.83787706640934548E+00/
C
C***FIRST EXECUTABLE STATEMENT  GAMLN
         IERR=0
         IF (Z.LE.0.0E0) GO TO 70
         IF (Z.GT.101.0E0) GO TO 10
         NZ = INT(Z)
         FZ = Z - FLOAT(NZ)
         IF (FZ.GT.0.0E0) GO TO 10
         IF (NZ.GT.100) GO TO 10
         GAMLN = GLN(NZ)
         RETURN
   10    CONTINUE
         WDTOL = R1MACH(4)
         WDTOL = AMAX1(WDTOL,0.5E-18)
         I1M = I1MACH(11)
         RLN = R1MACH(5)*FLOAT(I1M)
         FLN = AMIN1(RLN,20.0E0)
         FLN = AMAX1(FLN,3.0E0)
         FLN = FLN - 3.0E0
         ZM = 1.8000E0 + 0.3875E0*FLN
         MZ = INT(ZM) + 1
         ZMIN = FLOAT(MZ)
         ZDMY = Z
         ZINC = 0.0E0
         IF (Z.GE.ZMIN) GO TO 20
         ZINC = ZMIN - FLOAT(NZ)
         ZDMY = Z + ZINC
   20    CONTINUE
         ZP = 1.0E0/ZDMY
         T1 = CF(1)*ZP
         S = T1
         IF (ZP.LT.WDTOL) GO TO 40
         ZSQ = ZP*ZP
         TST = T1*WDTOL
         DO 30 K=2,22
            ZP = ZP*ZSQ
            TRM = CF(K)*ZP
            IF (ABS(TRM).LT.TST) GO TO 40
            S = S + TRM
   30    CONTINUE
   40    CONTINUE
         IF (ZINC.NE.0.0E0) GO TO 50
         TLG = ALOG(Z)
         GAMLN = Z*(TLG-1.0E0) + 0.5E0*(CON-TLG) + S
         RETURN
   50    CONTINUE
         ZP = 1.0E0
         NZ = INT(ZINC)
         DO 60 I=1,NZ
            ZP = ZP*(Z+FLOAT(I-1))
   60    CONTINUE
         TLG = ALOG(ZDMY)
         GAMLN = ZDMY*(TLG-1.0E0) - ALOG(ZP) + 0.5E0*(CON-TLG) + S
         RETURN
C
C
   70    CONTINUE
         IERR=1
         RETURN
      END
      function i1mach ( i )

c*********************************************************************72
c
cc I1MACH returns integer machine dependent constants.
c
c  Discussion:
c
c    Input/output unit numbers.
c
c      I1MACH(1) = the standard input unit.
c      I1MACH(2) = the standard output unit.
c      I1MACH(3) = the standard punch unit.
c      I1MACH(4) = the standard error message unit.
c
c    Words.
c
c      I1MACH(5) = the number of bits per integer storage unit.
c      I1MACH(6) = the number of characters per integer storage unit.
c
c    Integers.
c
c    Assume integers are represented in the S digit base A form:
c
c      Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
c
c    where 0 <= X(1:S-1) < A.
c
c      I1MACH(7) = A, the base.
c      I1MACH(8) = S, the number of base A digits.
c      I1MACH(9) = A**S-1, the largest integer.
c
c    Floating point numbers
c
c    Assume floating point numbers are represented in the T digit
c    base B form:
c
c      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
c
c    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
c
c      I1MACH(10) = B, the base.
c
c    Single precision
c
c      I1MACH(11) = T, the number of base B digits.
c      I1MACH(12) = EMIN, the smallest exponent E.
c      I1MACH(13) = EMAX, the largest exponent E.
c
c    Double precision
c
c      I1MACH(14) = T, the number of base B digits.
c      I1MACH(15) = EMIN, the smallest exponent E.
c      I1MACH(16) = EMAX, the largest exponent E.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 16.
c
c    Output, integer I1MACH, the value of the chosen parameter.
c
         implicit none

         integer i
         integer i1mach

         if ( i < 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'I1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
            write ( *, '(a,i12)' ) '  I = ', i
            i1mach = 0
            stop
         else if ( i == 1 ) then
            i1mach = 5
         else if ( i == 2 ) then
            i1mach = 6
         else if ( i == 3 ) then
            i1mach = 7
         else if ( i == 4 ) then
            i1mach = 6
         else if ( i == 5 ) then
            i1mach = 32
         else if ( i == 6 ) then
            i1mach = 4
         else if ( i == 7 ) then
            i1mach = 2
         else if ( i == 8 ) then
            i1mach = 31
         else if ( i == 9 ) then
            i1mach = 2147483647
         else if ( i == 10 ) then
            i1mach = 2
         else if ( i == 11 ) then
            i1mach = 24
         else if ( i == 12 ) then
            i1mach = -125
         else if ( i == 13 ) then
            i1mach = 128
         else if ( i == 14 ) then
            i1mach = 53
         else if ( i == 15 ) then
            i1mach = -1021
         else if ( i == 16 ) then
            i1mach = 1024
         else if ( 16 < i ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'I1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
            write ( *, '(a,i12)' ) '  I = ', i
            i1mach = 0
            stop
         end if

         return
      end
      function r1mach ( i )

c*********************************************************************72
c
cc R1MACH returns single precision real machine constants.
c
c  Discussion:
c
c    Assume that single precision real numbers are stored with a mantissa
c    of T digits in base B, with an exponent whose value must lie
c    between EMIN and EMAX.  Then for values of I between 1 and 5,
c    R1MACH will return the following values:
c
c      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
c      R1MACH(3) = B**(-T), the smallest relative spacing.
c      R1MACH(4) = B**(1-T), the largest relative spacing.
c      R1MACH(5) = log10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 5.
c
c    Output, real R1MACH, the value of the chosen parameter.
c
         implicit none

         integer i
         real r1mach

         if ( i < 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
            write ( *, '(a,i12)' ) '  I = ', i
            r1mach = 0.0E+00
            stop
         else if ( i == 1 ) then
            r1mach = 1.1754944E-38
         else if ( i == 2 ) then
            r1mach = 3.4028235E+38
         else if ( i == 3 ) then
            r1mach = 5.9604645E-08
         else if ( i == 4 ) then
            r1mach = 1.1920929E-07
         else if ( i == 5 ) then
            r1mach = 0.3010300E+00
         else if ( 5 < i ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
            write ( *, '(a,i12)' ) '  I = ', i
            r1mach = 0.0E+00
            stop
         end if

         return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
         implicit none

         character * ( 8 ) ampm
         integer d
         character * ( 8 ) date
         integer h
         integer m
         integer mm
         character * ( 9 ) month(12)
         integer n
         integer s
         character * ( 10 ) time
         integer y

         save month

         data month /
     $     'January  ', 'February ', 'March    ', 'April    ',
     $     'May      ', 'June     ', 'July     ', 'August   ',
     $     'September', 'October  ', 'November ', 'December ' /

         call date_and_time ( date, time )

         read ( date, '(i4,i2,i2)' ) y, m, d
         read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

         if ( h .lt. 12 ) then
            ampm = 'AM'
         else if ( h .eq. 12 ) then
            if ( n .eq. 0 .and. s .eq. 0 ) then
               ampm = 'Noon'
            else
               ampm = 'PM'
            end if
         else
            h = h - 12
            if ( h .lt. 12 ) then
               ampm = 'PM'
            else if ( h .eq. 12 ) then
               if ( n .eq. 0 .and. s .eq. 0 ) then
                  ampm = 'Midnight'
               else
                  ampm = 'AM'
               end if
            end if
         end if

         write ( *,
     $     '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     $     d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

         return
      end
      SUBROUTINE XERROR(MESS,NMESS,L1,L2)

c*********************************************************************72
c
C     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
C     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
C     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
C     ROUTINE.
C
         INTEGER NMESS, L1, L2, NN, NR, K, I, KMIN
         CHARACTER*(*) MESS
         NN=NMESS/70
         NR=NMESS-70*NN
         IF(NR.NE.0) NN=NN+1
         K=1
         PRINT 900
  900    FORMAT(/)
         DO 10 I=1,NN
            KMIN=MIN0(K+69,NMESS)
            PRINT *, MESS(K:KMIN)
            K=K+70
   10    CONTINUE
         PRINT 900
         RETURN
      END
      DOUBLE PRECISION FUNCTION ZABS(ZR, ZI)

c*********************************************************************72
c
cc ZABS carries out double precision complex absolute values.
c
C***BEGIN PROLOGUE  ZABS
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
C     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZABS
         DOUBLE PRECISION ZR, ZI, U, V, Q, S
         U = DABS(ZR)
         V = DABS(ZI)
         S = U + V
C-----------------------------------------------------------------------
C     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
C     TRUE FLOATING ZERO
C-----------------------------------------------------------------------
         S = S*1.0D+0
         IF (S.EQ.0.0D+0) GO TO 20
         IF (U.GT.V) GO TO 10
         Q = U/V
         ZABS = V*DSQRT(1.D+0+Q*Q)
         RETURN
   10    Q = V/U
         ZABS = U*DSQRT(1.D+0+Q*Q)
         RETURN
   20    ZABS = 0.0D+0
         RETURN
      END
      SUBROUTINE ZACAI(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     * ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZACAI
C***REFER TO  ZAIRY
C
C     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
C     IS CALLED FROM ZAIRY.
C
C***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,ZABS
C***END PROLOGUE  ZACAI
C     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
         EXTERNAL ZABS
         DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
     *    CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
     *    RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, D1MACH, ZABS
         INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
         DIMENSION YR(N), YI(N), CYR(2), CYI(2)
         DATA PI / 3.14159265358979324D0 /
         NZ = 0
         ZNR = -ZR
         ZNI = -ZI
         AZ = ZABS(ZR,ZI)
         NN = N
         DFNU = FNU + DBLE(FLOAT(N-1))
         IF (AZ.LE.2.0D0) GO TO 10
         IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
         CALL ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
         GO TO 40
   20    CONTINUE
         IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
         CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 80
         GO TO 40
   30    CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
         CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
         IF(NW.LT.0) GO TO 80
   40    CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
         CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
         IF (NW.NE.0) GO TO 80
         FMR = DBLE(FLOAT(MR))
         SGN = -DSIGN(PI,FMR)
         CSGNR = 0.0D0
         CSGNI = SGN
         IF (KODE.EQ.1) GO TO 50
         YY = -ZNI
         CSGNR = -CSGNI*DSIN(YY)
         CSGNI = CSGNI*DCOS(YY)
   50    CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(SNGL(FNU))
         ARG = (FNU-DBLE(FLOAT(INU)))*SGN
         CSPNR = DCOS(ARG)
         CSPNI = DSIN(ARG)
         IF (MOD(INU,2).EQ.0) GO TO 60
         CSPNR = -CSPNR
         CSPNI = -CSPNI
   60    CONTINUE
         C1R = CYR(1)
         C1I = CYI(1)
         C2R = YR(1)
         C2I = YI(1)
         IF (KODE.EQ.1) GO TO 70
         IUF = 0
         ASCLE = 1.0D+3*D1MACH(1)/TOL
         CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
         NZ = NZ + NW
   70    CONTINUE
         YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
         YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
         RETURN
   80    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      SUBROUTINE ZACON(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL,
     * TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZACON
C***REFER TO  ZBESK,ZBESH
C
C     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***ROUTINES CALLED  ZBINU,ZBKNU,ZS1S2,D1MACH,ZABS,ZMLT
C***END PROLOGUE  ZACON
C     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
C    *S1,S2,Y,Z,ZN
         EXTERNAL ZABS
         DOUBLE PRECISION ALIM, ARG, ASCLE, AS2, AZN, BRY, BSCLE, CKI,
     *    CKR, CONER, CPN, CSCL, CSCR, CSGNI, CSGNR, CSPNI, CSPNR,
     *    CSR, CSRR, CSSR, CYI, CYR, C1I, C1M, C1R, C2I, C2R, ELIM, FMR,
     *    FN, FNU, FNUL, PI, PTI, PTR, RAZN, RL, RZI, RZR, SC1I, SC1R,
     *    SC2I, SC2R, SGN, SPN, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR,
     *    YY, ZEROR, ZI, ZNI, ZNR, ZR, D1MACH, ZABS
         INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
         DIMENSION YR(N), YI(N), CYR(2), CYI(2), CSSR(3), CSRR(3), BRY(3)
         DATA PI / 3.14159265358979324D0 /
         DATA ZEROR,CONER / 0.0D0,1.0D0 /
         NZ = 0
         ZNR = -ZR
         ZNI = -ZI
         NN = N
         CALL ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL,
     *    ELIM, ALIM)
         IF (NW.LT.0) GO TO 90
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
         NN = MIN0(2,N)
         CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
         IF (NW.NE.0) GO TO 90
         S1R = CYR(1)
         S1I = CYI(1)
         FMR = DBLE(FLOAT(MR))
         SGN = -DSIGN(PI,FMR)
         CSGNR = ZEROR
         CSGNI = SGN
         IF (KODE.EQ.1) GO TO 10
         YY = -ZNI
         CPN = DCOS(YY)
         SPN = DSIN(YY)
         CALL ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI)
   10    CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(SNGL(FNU))
         ARG = (FNU-DBLE(FLOAT(INU)))*SGN
         CPN = DCOS(ARG)
         SPN = DSIN(ARG)
         CSPNR = CPN
         CSPNI = SPN
         IF (MOD(INU,2).EQ.0) GO TO 20
         CSPNR = -CSPNR
         CSPNI = -CSPNI
   20    CONTINUE
         IUF = 0
         C1R = S1R
         C1I = S1I
         C2R = YR(1)
         C2I = YI(1)
         ASCLE = 1.0D+3*D1MACH(1)/TOL
         IF (KODE.EQ.1) GO TO 30
         CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
         NZ = NZ + NW
         SC1R = C1R
         SC1I = C1I
   30    CONTINUE
         CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
         CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
         YR(1) = STR + PTR
         YI(1) = STI + PTI
         IF (N.EQ.1) RETURN
         CSPNR = -CSPNR
         CSPNI = -CSPNI
         S2R = CYR(2)
         S2I = CYI(2)
         C1R = S2R
         C1I = S2I
         C2R = YR(2)
         C2I = YI(2)
         IF (KODE.EQ.1) GO TO 40
         CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
         NZ = NZ + NW
         SC2R = C1R
         SC2I = C1I
   40    CONTINUE
         CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
         CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
         YR(2) = STR + PTR
         YI(2) = STI + PTI
         IF (N.EQ.2) RETURN
         CSPNR = -CSPNR
         CSPNI = -CSPNI
         AZN = ZABS(ZNR,ZNI)
         RAZN = 1.0D0/AZN
         STR = ZNR*RAZN
         STI = -ZNI*RAZN
         RZR = (STR+STR)*RAZN
         RZI = (STI+STI)*RAZN
         FN = FNU + 1.0D0
         CKR = FN*RZR
         CKI = FN*RZI
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
         CSCL = 1.0D0/TOL
         CSCR = TOL
         CSSR(1) = CSCL
         CSSR(2) = CONER
         CSSR(3) = CSCR
         CSRR(1) = CSCR
         CSRR(2) = CONER
         CSRR(3) = CSCL
         BRY(1) = ASCLE
         BRY(2) = 1.0D0/ASCLE
         BRY(3) = D1MACH(2)
         AS2 = ZABS(S2R,S2I)
         KFLAG = 2
         IF (AS2.GT.BRY(1)) GO TO 50
         KFLAG = 1
         GO TO 60
   50    CONTINUE
         IF (AS2.LT.BRY(2)) GO TO 60
         KFLAG = 3
   60    CONTINUE
         BSCLE = BRY(KFLAG)
         S1R = S1R*CSSR(KFLAG)
         S1I = S1I*CSSR(KFLAG)
         S2R = S2R*CSSR(KFLAG)
         S2I = S2I*CSSR(KFLAG)
         CSR = CSRR(KFLAG)
         DO 80 I=3,N
            STR = S2R
            STI = S2I
            S2R = CKR*STR - CKI*STI + S1R
            S2I = CKR*STI + CKI*STR + S1I
            S1R = STR
            S1I = STI
            C1R = S2R*CSR
            C1I = S2I*CSR
            STR = C1R
            STI = C1I
            C2R = YR(I)
            C2I = YI(I)
            IF (KODE.EQ.1) GO TO 70
            IF (IUF.LT.0) GO TO 70
            CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
            NZ = NZ + NW
            SC1R = SC2R
            SC1I = SC2I
            SC2R = C1R
            SC2I = C1I
            IF (IUF.NE.3) GO TO 70
            IUF = -4
            S1R = SC1R*CSSR(KFLAG)
            S1I = SC1I*CSSR(KFLAG)
            S2R = SC2R*CSSR(KFLAG)
            S2I = SC2I*CSSR(KFLAG)
            STR = SC2R
            STI = SC2I
   70       CONTINUE
            PTR = CSPNR*C1R - CSPNI*C1I
            PTI = CSPNR*C1I + CSPNI*C1R
            YR(I) = PTR + CSGNR*C2R - CSGNI*C2I
            YI(I) = PTI + CSGNR*C2I + CSGNI*C2R
            CKR = CKR + RZR
            CKI = CKI + RZI
            CSPNR = -CSPNR
            CSPNI = -CSPNI
            IF (KFLAG.GE.3) GO TO 80
            PTR = DABS(C1R)
            PTI = DABS(C1I)
            C1M = DMAX1(PTR,PTI)
            IF (C1M.LE.BSCLE) GO TO 80
            KFLAG = KFLAG + 1
            BSCLE = BRY(KFLAG)
            S1R = S1R*CSR
            S1I = S1I*CSR
            S2R = STR
            S2I = STI
            S1R = S1R*CSSR(KFLAG)
            S1I = S1I*CSSR(KFLAG)
            S2R = S2R*CSSR(KFLAG)
            S2I = S2I*CSSR(KFLAG)
            CSR = CSRR(KFLAG)
   80    CONTINUE
         RETURN
   90    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      SUBROUTINE ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)

c*********************************************************************72
c
cc ZAIRY computes a sequence of complex Airy Ai functions.
c
C***BEGIN PROLOGUE  ZAIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
C         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
C         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
C         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
C         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
C
C         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
C         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
C         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
C         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             AI=AI(Z)                ON ID=0 OR
C                             AI=DAI(Z)/DZ            ON ID=1
C                        = 2  RETURNS
C                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
C                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)
C
C         OUTPUT     AIR,AII ARE DOUBLE PRECISION
C           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           NZ     - UNDERFLOW INDICATOR
C                    NZ= 0   , NORMAL RETURN
C                    NZ= 1   , AI=CMPLX(0.0D0,0.0D0) DUE TO UNDERFLOW IN
C                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
C         FUNCTIONS BY
C
C            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
C                           C=1.0/(PI*SQRT(3.0))
C                            ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZACAI,ZBKNU,ZEXP,ZSQRT,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZAIRY
C     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
         EXTERNAL ZABS
         DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
     *    CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
     *    DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
     *    S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
     *    ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS, ALAZ, BB
         INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
         DIMENSION CYR(1), CYI(1)
         DATA TTH, C1, C2, COEF /6.66666666666666667D-01,
     *    3.55028053887817240D-01,2.58819403792806799D-01,
     *    1.83776298473930683D-01/
         DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZAIRY
         IERR = 0
         NZ=0
         IF (ID.LT.0 .OR. ID.GT.1) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (IERR.NE.0) RETURN
         AZ = ZABS(ZR,ZI)
         TOL = DMAX1(D1MACH(4),1.0D-18)
         FID = DBLE(FLOAT(ID))
         IF (AZ.GT.1.0D0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
         S1R = CONER
         S1I = CONEI
         S2R = CONER
         S2I = CONEI
         IF (AZ.LT.TOL) GO TO 170
         AA = AZ*AZ
         IF (AA.LT.TOL/AZ) GO TO 40
         TRM1R = CONER
         TRM1I = CONEI
         TRM2R = CONER
         TRM2I = CONEI
         ATRM = 1.0D0
         STR = ZR*ZR - ZI*ZI
         STI = ZR*ZI + ZI*ZR
         Z3R = STR*ZR - STI*ZI
         Z3I = STR*ZI + STI*ZR
         AZ3 = AZ*AA
         AK = 2.0D0 + FID
         BK = 3.0D0 - FID - FID
         CK = 4.0D0 - FID
         DK = 3.0D0 + FID + FID
         D1 = AK*DK
         D2 = BK*CK
         AD = DMIN1(D1,D2)
         AK = 24.0D0 + 9.0D0*FID
         BK = 30.0D0 - 9.0D0*FID
         DO 30 K=1,25
            STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
            TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
            TRM1R = STR
            S1R = S1R + TRM1R
            S1I = S1I + TRM1I
            STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
            TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
            TRM2R = STR
            S2R = S2R + TRM2R
            S2I = S2I + TRM2I
            ATRM = ATRM*AZ3/AD
            D1 = D1 + AK
            D2 = D2 + BK
            AD = DMIN1(D1,D2)
            IF (ATRM.LT.TOL*AD) GO TO 40
            AK = AK + 18.0D0
            BK = BK + 18.0D0
   30    CONTINUE
   40    CONTINUE
         IF (ID.EQ.1) GO TO 50
         AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
         AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
         IF (KODE.EQ.1) RETURN
         CALL ZSQRT(ZR, ZI, STR, STI)
         ZTAR = TTH*(ZR*STR-ZI*STI)
         ZTAI = TTH*(ZR*STI+ZI*STR)
         CALL ZEXP(ZTAR, ZTAI, STR, STI)
         PTR = AIR*STR - AII*STI
         AII = AIR*STI + AII*STR
         AIR = PTR
         RETURN
   50    CONTINUE
         AIR = -S2R*C2
         AII = -S2I*C2
         IF (AZ.LE.TOL) GO TO 60
         STR = ZR*S1R - ZI*S1I
         STI = ZR*S1I + ZI*S1R
         CC = C1/(1.0D0+FID)
         AIR = AIR + CC*(STR*ZR-STI*ZI)
         AII = AII + CC*(STR*ZI+STI*ZR)
   60    CONTINUE
         IF (KODE.EQ.1) RETURN
         CALL ZSQRT(ZR, ZI, STR, STI)
         ZTAR = TTH*(ZR*STR-ZI*STI)
         ZTAI = TTH*(ZR*STI+ZI*STR)
         CALL ZEXP(ZTAR, ZTAI, STR, STI)
         PTR = STR*AIR - STI*AII
         AII = STR*AII + STI*AIR
         AIR = PTR
         RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70    CONTINUE
         FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         R1M5 = D1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
         K1 = I1MACH(14) - 1
         AA = R1M5*DBLE(FLOAT(K1))
         DIG = DMIN1(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + DMAX1(-AA,-41.45D0)
         RL = 1.2D0*DIG + 3.0D0
         ALAZ = DLOG(AZ)
C--------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
         AA=0.5D0/TOL
         BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
         AA=DMIN1(AA,BB)
         AA=AA**TTH
         IF (AZ.GT.AA) GO TO 260
         AA=DSQRT(AA)
         IF (AZ.GT.AA) IERR=3
         CALL ZSQRT(ZR, ZI, CSQR, CSQI)
         ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
         ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
         IFLAG = 0
         SFAC = 1.0D0
         AK = ZTAI
         IF (ZR.GE.0.0D0) GO TO 80
         BK = ZTAR
         CK = -DABS(BK)
         ZTAR = CK
         ZTAI = AK
   80    CONTINUE
         IF (ZI.NE.0.0D0) GO TO 90
         IF (ZR.GT.0.0D0) GO TO 90
         ZTAR = 0.0D0
         ZTAI = AK
   90    CONTINUE
         AA = ZTAR
         IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
         IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         IF (AA.GT.(-ALIM)) GO TO 100
         AA = -AA + 0.25D0*ALAZ
         IFLAG = 1
         SFAC = TOL
         IF (AA.GT.ELIM) GO TO 270
  100    CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
         MR = 1
         IF (ZI.LT.0.0D0) MR = -1
         CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
     *    ELIM, ALIM)
         IF (NN.LT.0) GO TO 280
         NZ = NZ + NN
         GO TO 130
  110    CONTINUE
         IF (KODE.EQ.2) GO TO 120
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
         IF (AA.LT.ALIM) GO TO 120
         AA = -AA - 0.25D0*ALAZ
         IFLAG = 2
         SFAC = 1.0D0/TOL
         IF (AA.LT.(-ELIM)) GO TO 210
  120    CONTINUE
         CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
     *    ALIM)
  130    CONTINUE
         S1R = CYR(1)*COEF
         S1I = CYI(1)*COEF
         IF (IFLAG.NE.0) GO TO 150
         IF (ID.EQ.1) GO TO 140
         AIR = CSQR*S1R - CSQI*S1I
         AII = CSQR*S1I + CSQI*S1R
         RETURN
  140    CONTINUE
         AIR = -(ZR*S1R-ZI*S1I)
         AII = -(ZR*S1I+ZI*S1R)
         RETURN
  150    CONTINUE
         S1R = S1R*SFAC
         S1I = S1I*SFAC
         IF (ID.EQ.1) GO TO 160
         STR = S1R*CSQR - S1I*CSQI
         S1I = S1R*CSQI + S1I*CSQR
         S1R = STR
         AIR = S1R/SFAC
         AII = S1I/SFAC
         RETURN
  160    CONTINUE
         STR = -(S1R*ZR-S1I*ZI)
         S1I = -(S1R*ZI+S1I*ZR)
         S1R = STR
         AIR = S1R/SFAC
         AII = S1I/SFAC
         RETURN
  170    CONTINUE
         AA = 1.0D+3*D1MACH(1)
         S1R = ZEROR
         S1I = ZEROI
         IF (ID.EQ.1) GO TO 190
         IF (AZ.LE.AA) GO TO 180
         S1R = C2*ZR
         S1I = C2*ZI
  180    CONTINUE
         AIR = C1 - S1R
         AII = -S1I
         RETURN
  190    CONTINUE
         AIR = -C2
         AII = 0.0D0
         AA = DSQRT(AA)
         IF (AZ.LE.AA) GO TO 200
         S1R = 0.5D0*(ZR*ZR-ZI*ZI)
         S1I = ZR*ZI
  200    CONTINUE
         AIR = AIR + C1*S1R
         AII = AII + C1*S1I
         RETURN
  210    CONTINUE
         NZ = 1
         AIR = ZEROR
         AII = ZEROI
         RETURN
  270    CONTINUE
         NZ = 0
         IERR=2
         RETURN
  280    CONTINUE
         IF(NN.EQ.(-1)) GO TO 270
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         IERR=4
         NZ=0
         RETURN
      END
      SUBROUTINE ZASYI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZASYI
C***REFER TO  ZBESI,ZBESK
C
C     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***ROUTINES CALLED  D1MACH,ZABS,ZDIV,ZEXP,ZMLT,ZSQRT
C***END PROLOGUE  ZASYI
C     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
         EXTERNAL ZABS
         DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL,
     *    AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI,
     *    CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I,
     *    P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I,
     *    S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR, D1MACH, ZABS
         INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
         DIMENSION YR(N), YI(N)
         DATA PI, RTPI  /3.14159265358979324D0 , 0.159154943091895336D0 /
         DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
         NZ = 0
         AZ = ZABS(ZR,ZI)
         ARM = 1.0D+3*D1MACH(1)
         RTR1 = DSQRT(ARM)
         IL = MIN0(2,N)
         DFNU = FNU + DBLE(FLOAT(N-IL))
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         RAZ = 1.0D0/AZ
         STR = ZR*RAZ
         STI = -ZI*RAZ
         AK1R = RTPI*STR*RAZ
         AK1I = RTPI*STI*RAZ
         CALL ZSQRT(AK1R, AK1I, AK1R, AK1I)
         CZR = ZR
         CZI = ZI
         IF (KODE.NE.2) GO TO 10
         CZR = ZEROR
         CZI = ZI
   10    CONTINUE
         IF (DABS(CZR).GT.ELIM) GO TO 100
         DNU2 = DFNU + DFNU
         KODED = 1
         IF ((DABS(CZR).GT.ALIM) .AND. (N.GT.2)) GO TO 20
         KODED = 0
         CALL ZEXP(CZR, CZI, STR, STI)
         CALL ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I)
   20    CONTINUE
         FDN = 0.0D0
         IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
         EZR = ZR*8.0D0
         EZI = ZI*8.0D0
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
         AEZ = 8.0D0*AZ
         S = TOL/AEZ
         JL = INT(SNGL(RL+RL)) + 2
         P1R = ZEROR
         P1I = ZEROI
         IF (ZI.EQ.0.0D0) GO TO 30
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
         INU = INT(SNGL(FNU))
         ARG = (FNU-DBLE(FLOAT(INU)))*PI
         INU = INU + N - IL
         AK = -DSIN(ARG)
         BK = DCOS(ARG)
         IF (ZI.LT.0.0D0) BK = -BK
         P1R = AK
         P1I = BK
         IF (MOD(INU,2).EQ.0) GO TO 30
         P1R = -P1R
         P1I = -P1I
   30    CONTINUE
         DO 70 K=1,IL
            SQK = FDN - 1.0D0
            ATOL = S*DABS(SQK)
            SGN = 1.0D0
            CS1R = CONER
            CS1I = CONEI
            CS2R = CONER
            CS2I = CONEI
            CKR = CONER
            CKI = CONEI
            AK = 0.0D0
            AA = 1.0D0
            BB = AEZ
            DKR = EZR
            DKI = EZI
            DO 40 J=1,JL
               CALL ZDIV(CKR, CKI, DKR, DKI, STR, STI)
               CKR = STR*SQK
               CKI = STI*SQK
               CS2R = CS2R + CKR
               CS2I = CS2I + CKI
               SGN = -SGN
               CS1R = CS1R + CKR*SGN
               CS1I = CS1I + CKI*SGN
               DKR = DKR + EZR
               DKI = DKI + EZI
               AA = AA*DABS(SQK)/BB
               BB = BB + AEZ
               AK = AK + 8.0D0
               SQK = SQK - AK
               IF (AA.LE.ATOL) GO TO 50
   40       CONTINUE
            GO TO 110
   50       CONTINUE
            S2R = CS1R
            S2I = CS1I
            IF (ZR+ZR.GE.ELIM) GO TO 60
            TZR = ZR + ZR
            TZI = ZI + ZI
            CALL ZEXP(-TZR, -TZI, STR, STI)
            CALL ZMLT(STR, STI, P1R, P1I, STR, STI)
            CALL ZMLT(STR, STI, CS2R, CS2I, STR, STI)
            S2R = S2R + STR
            S2I = S2I + STI
   60       CONTINUE
            FDN = FDN + 8.0D0*DFNU + 4.0D0
            P1R = -P1R
            P1I = -P1I
            M = N - IL + K
            YR(M) = S2R*AK1R - S2I*AK1I
            YI(M) = S2R*AK1I + S2I*AK1R
   70    CONTINUE
         IF (N.LE.2) RETURN
         NN = N
         K = NN - 2
         AK = DBLE(FLOAT(K))
         STR = ZR*RAZ
         STI = -ZI*RAZ
         RZR = (STR+STR)*RAZ
         RZI = (STI+STI)*RAZ
         IB = 3
         DO 80 I=IB,NN
            YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
            YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
            AK = AK - 1.0D0
            K = K - 1
   80    CONTINUE
         IF (KODED.EQ.0) RETURN
         CALL ZEXP(CZR, CZI, CKR, CKI)
         DO 90 I=1,NN
            STR = YR(I)*CKR - YI(I)*CKI
            YI(I) = YR(I)*CKI + YI(I)*CKR
            YR(I) = STR
   90    CONTINUE
         RETURN
  100    CONTINUE
         NZ = -1
         RETURN
  110    CONTINUE
         NZ=-2
         RETURN
      END
      SUBROUTINE ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)

c*********************************************************************72
c
cc ZBESH computes a sequence of complex Hankel functions.
c
C***BEGIN PROLOGUE  ZBESH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
C         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
C         Z.NE.CMPLX(0.0,0.0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI.
C         ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS
C
C         CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1.
C
C         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND
C         LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE
C         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PT.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z),   J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
C                                  J=1,...,N  ,  I**2=-1
C           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
C           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(J)=H(M,FNU+J-1,Z)  OR
C                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
C                    DEPENDING ON KODE, I**2=-1.
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
C                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0)
C                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR
C                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY
C                              HALF PLANES, NZ STATES ONLY THE NUMBER
C                              OF UNDERFLOWS.
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU TOO
C                            LARGE OR CABS(Z) TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE RELATION
C
C         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
C             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
C
C         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
C         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED
C         TO THE LEFT HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
C         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
C         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
C         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
C         WHOLE Z PLANE FOR Z TO INFINITY.
C
C         FOR NEGATIVE ORDERS,THE FORMULAE
C
C               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
C               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
C                         I**2=-1
C
C         CAN BE USED.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0D-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESH
C
C     COMPLEX CY,Z,ZN,ZT,CSGN
         EXTERNAL ZABS
         DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM,
     *    FMM, FN, FNU, FNUL, HPI, RHPI, RL, R1M5, SGN, STR, TOL, UFL, ZI,
     *    ZNI, ZNR, ZR, ZTI, D1MACH, ZABS, BB, ASCLE, RTOL, ATOL, STI,
     *    CSGNR, CSGNI
         INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     *    MM, MR, N, NN, NUF, NW, NZ, I1MACH
         DIMENSION CYR(N), CYI(N)
C
         DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESH
         IERR = 0
         NZ=0
         IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
         IF (FNU.LT.0.0D0) IERR=1
         IF (M.LT.1 .OR. M.GT.2) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
         TOL = DMAX1(D1MACH(4),1.0D-18)
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         R1M5 = D1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
         K1 = I1MACH(14) - 1
         AA = R1M5*DBLE(FLOAT(K1))
         DIG = DMIN1(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + DMAX1(-AA,-41.45D0)
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
         RL = 1.2D0*DIG + 3.0D0
         FN = FNU + DBLE(FLOAT(NN-1))
         MM = 3 - M - M
         FMM = DBLE(FLOAT(MM))
         ZNR = FMM*ZI
         ZNI = -FMM*ZR
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
         AZ = ZABS(ZR,ZI)
         AA = 0.5D0/TOL
         BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
         AA = DMIN1(AA,BB)
         IF (AZ.GT.AA) GO TO 260
         IF (FN.GT.AA) GO TO 260
         AA = DSQRT(AA)
         IF (AZ.GT.AA) IERR=3
         IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
         UFL = D1MACH(1)*1.0D+3
         IF (AZ.LT.UFL) GO TO 230
         IF (FNU.GT.FNUL) GO TO 90
         IF (FN.LE.1.0D0) GO TO 70
         IF (FN.GT.2.0D0) GO TO 60
         IF (AZ.GT.TOL) GO TO 70
         ARG = 0.5D0*AZ
         ALN = -FN*DLOG(ARG)
         IF (ALN.GT.ELIM) GO TO 230
         GO TO 70
   60    CONTINUE
         CALL ZUOIK(ZNR, ZNI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     *    ALIM)
         IF (NUF.LT.0) GO TO 230
         NZ = NZ + NUF
         NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
         IF (NN.EQ.0) GO TO 140
   70    CONTINUE
         IF ((ZNR.LT.0.0D0) .OR. (ZNR.EQ.0.0D0 .AND. ZNI.LT.0.0D0 .AND.
     *    M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
         CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NZ, TOL, ELIM, ALIM)
         GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80    CONTINUE
         MR = -MM
         CALL ZACON(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     *    TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 240
         NZ=NW
         GO TO 110
   90    CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
         MR = 0
         IF ((ZNR.GE.0.0D0) .AND. (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0 .OR.
     *    M.NE.2)) GO TO 100
         MR = -MM
         IF (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0) GO TO 100
         ZNR = -ZNR
         ZNI = -ZNI
  100    CONTINUE
         CALL ZBUNK(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 240
         NZ = NZ + NW
  110    CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
         SGN = DSIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(SNGL(FNU))
         INUH = INU/2
         IR = INU - 2*INUH
         ARG = (FNU-DBLE(FLOAT(INU-IR)))*SGN
         RHPI = 1.0D0/SGN
C     ZNI = RHPI*DCOS(ARG)
C     ZNR = -RHPI*DSIN(ARG)
         CSGNI = RHPI*DCOS(ARG)
         CSGNR = -RHPI*DSIN(ARG)
         IF (MOD(INUH,2).EQ.0) GO TO 120
C     ZNR = -ZNR
C     ZNI = -ZNI
         CSGNR = -CSGNR
         CSGNI = -CSGNI
  120    CONTINUE
         ZTI = -FMM
         RTOL = 1.0D0/TOL
         ASCLE = UFL*RTOL
         DO 130 I=1,NN
C       STR = CYR(I)*ZNR - CYI(I)*ZNI
C       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
C       CYR(I) = STR
C       STR = -ZNI*ZTI
C       ZNI = ZNR*ZTI
C       ZNR = STR
            AA = CYR(I)
            BB = CYI(I)
            ATOL = 1.0D0
            IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 135
            AA = AA*RTOL
            BB = BB*RTOL
            ATOL = TOL
  135       CONTINUE
            STR = AA*CSGNR - BB*CSGNI
            STI = AA*CSGNI + BB*CSGNR
            CYR(I) = STR*ATOL
            CYI(I) = STI*ATOL
            STR = -CSGNI*ZTI
            CSGNI = CSGNR*ZTI
            CSGNR = STR
  130    CONTINUE
         RETURN
  140    CONTINUE
         IF (ZNR.LT.0.0D0) GO TO 230
         RETURN
  230    CONTINUE
         NZ=0
         IERR=2
         RETURN
  240    CONTINUE
         IF(NW.EQ.(-1)) GO TO 230
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE ZBESI(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)

c*********************************************************************72
c
cc ZBESI computes a sequence of complex Bessel I functions.
c
C***BEGIN PROLOGUE  ZBESI
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C                    ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESI RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)
C
C         WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL I FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=I(FNU+J-1,Z), J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(J)=I(FNU+J-1,Z)  OR
C                    CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
C                    DEPENDING ON KODE, X=REAL(Z)
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0)
C                              J = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
C                            LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
C         SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
C         THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
C         NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
C         UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
C         FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
C         SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
C
C         THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
C         CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
C
C         I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z).GT.0.0
C                       M = +I OR -I,  I**2=-1
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
C         NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBINU,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESI
C     COMPLEX CONE,CSGN,CW,CY,CZERO,Z,ZN
         EXTERNAL ZABS
         DOUBLE PRECISION AA, ALIM, ARG, CONEI, CONER, CSGNI, CSGNR, CYI,
     *    CYR, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR,
     *    ZR, D1MACH, AZ, BB, FN, ZABS, ASCLE, RTOL, ATOL, STI
         INTEGER I, IERR, INU, K, KODE, K1,K2,N,NZ,NN, I1MACH
         DIMENSION CYR(N), CYI(N)
         DATA PI /3.14159265358979324D0/
         DATA CONER, CONEI /1.0D0,0.0D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESI
         IERR = 0
         NZ=0
         IF (FNU.LT.0.0D0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
         TOL = DMAX1(D1MACH(4),1.0D-18)
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         R1M5 = D1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
         K1 = I1MACH(14) - 1
         AA = R1M5*DBLE(FLOAT(K1))
         DIG = DMIN1(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + DMAX1(-AA,-41.45D0)
         RL = 1.2D0*DIG + 3.0D0
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
         AZ = ZABS(ZR,ZI)
         FN = FNU+DBLE(FLOAT(N-1))
         AA = 0.5D0/TOL
         BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
         AA = DMIN1(AA,BB)
         IF (AZ.GT.AA) GO TO 260
         IF (FN.GT.AA) GO TO 260
         AA = DSQRT(AA)
         IF (AZ.GT.AA) IERR=3
         IF (FN.GT.AA) IERR=3
         ZNR = ZR
         ZNI = ZI
         CSGNR = CONER
         CSGNI = CONEI
         IF (ZR.GE.0.0D0) GO TO 40
         ZNR = -ZR
         ZNI = -ZI
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         INU = INT(SNGL(FNU))
         ARG = (FNU-DBLE(FLOAT(INU)))*PI
         IF (ZI.LT.0.0D0) ARG = -ARG
         CSGNR = DCOS(ARG)
         CSGNI = DSIN(ARG)
         IF (MOD(INU,2).EQ.0) GO TO 40
         CSGNR = -CSGNR
         CSGNI = -CSGNI
   40    CONTINUE
         CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     *    ELIM, ALIM)
         IF (NZ.LT.0) GO TO 120
         IF (ZR.GE.0.0D0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
C-----------------------------------------------------------------------
         NN = N - NZ
         IF (NN.EQ.0) RETURN
         RTOL = 1.0D0/TOL
         ASCLE = D1MACH(1)*RTOL*1.0D+3
         DO 50 I=1,NN
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
            AA = CYR(I)
            BB = CYI(I)
            ATOL = 1.0D0
            IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 55
            AA = AA*RTOL
            BB = BB*RTOL
            ATOL = TOL
   55       CONTINUE
            STR = AA*CSGNR - BB*CSGNI
            STI = AA*CSGNI + BB*CSGNR
            CYR(I) = STR*ATOL
            CYI(I) = STI*ATOL
            CSGNR = -CSGNR
            CSGNI = -CSGNI
   50    CONTINUE
         RETURN
  120    CONTINUE
         IF(NZ.EQ.(-2)) GO TO 130
         NZ = 0
         IERR=2
         RETURN
  130    CONTINUE
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)

c*********************************************************************72
c
cc ZBESJ computes a sequence of complex Bessel J functions.
c
C***BEGIN PROLOGUE  ZBESJ
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF FIRST KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESJ RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=J(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=J(FNU+I-1,Z)  OR
C                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE, Y=AIMAG(Z).
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET  ZERO DUE
C                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
C                              I = N-NZ+1,...,N
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0
C
C         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0
C
C         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
C         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBINU,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESJ
C
C     COMPLEX CI,CSGN,CY,Z,ZN
         EXTERNAL ZABS
         DOUBLE PRECISION AA, ALIM, ARG, CII, CSGNI, CSGNR, CYI, CYR, DIG,
     *    ELIM, FNU, FNUL, HPI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR, ZR,
     *    D1MACH, BB, FN, AZ, ZABS, ASCLE, RTOL, ATOL, STI
         INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, N, NL, NZ, I1MACH
         DIMENSION CYR(N), CYI(N)
         DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESJ
         IERR = 0
         NZ=0
         IF (FNU.LT.0.0D0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
         TOL = DMAX1(D1MACH(4),1.0D-18)
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         R1M5 = D1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
         K1 = I1MACH(14) - 1
         AA = R1M5*DBLE(FLOAT(K1))
         DIG = DMIN1(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + DMAX1(-AA,-41.45D0)
         RL = 1.2D0*DIG + 3.0D0
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
         AZ = ZABS(ZR,ZI)
         FN = FNU+DBLE(FLOAT(N-1))
         AA = 0.5D0/TOL
         BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
         AA = DMIN1(AA,BB)
         IF (AZ.GT.AA) GO TO 260
         IF (FN.GT.AA) GO TO 260
         AA = DSQRT(AA)
         IF (AZ.GT.AA) IERR=3
         IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
         CII = 1.0D0
         INU = INT(SNGL(FNU))
         INUH = INU/2
         IR = INU - 2*INUH
         ARG = (FNU-DBLE(FLOAT(INU-IR)))*HPI
         CSGNR = DCOS(ARG)
         CSGNI = DSIN(ARG)
         IF (MOD(INUH,2).EQ.0) GO TO 40
         CSGNR = -CSGNR
         CSGNI = -CSGNI
   40    CONTINUE
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
         ZNR = ZI
         ZNI = -ZR
         IF (ZI.GE.0.0D0) GO TO 50
         ZNR = -ZNR
         ZNI = -ZNI
         CSGNI = -CSGNI
         CII = -CII
   50    CONTINUE
         CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     *    ELIM, ALIM)
         IF (NZ.LT.0) GO TO 130
         NL = N - NZ
         IF (NL.EQ.0) RETURN
         RTOL = 1.0D0/TOL
         ASCLE = D1MACH(1)*RTOL*1.0D+3
         DO 60 I=1,NL
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
            AA = CYR(I)
            BB = CYI(I)
            ATOL = 1.0D0
            IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 55
            AA = AA*RTOL
            BB = BB*RTOL
            ATOL = TOL
   55       CONTINUE
            STR = AA*CSGNR - BB*CSGNI
            STI = AA*CSGNI + BB*CSGNR
            CYR(I) = STR*ATOL
            CYI(I) = STI*ATOL
            STR = -CSGNI*CII
            CSGNI = CSGNR*CII
            CSGNR = STR
   60    CONTINUE
         RETURN
  130    CONTINUE
         IF(NZ.EQ.(-2)) GO TO 140
         NZ = 0
         IERR = 2
         RETURN
  140    CONTINUE
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE ZBESK(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)

c*********************************************************************72
c
cc ZBESK computes a sequence of complex Bessel K functions.
c
C***BEGIN PROLOGUE  ZBESK
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
C             BESSEL FUNCTION OF THE THIRD KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C
C         ON KODE=1, ZBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
C         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESK
C         RETURNS THE SCALED K FUNCTIONS,
C
C         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
C
C         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0D0
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=K(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
C                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C                    DEPENDING ON KODE
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
C                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
C                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
C                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
C                              IN THE SEQUENCE.
C
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
C         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
C         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
C         HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
C         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
C
C         FOR NEGATIVE ORDERS, THE FORMULA
C
C                       K(-FNU,Z) = K(FNU,Z)
C
C         CAN BE USED.
C
C         ZBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
C         AVAILABLE.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESK
C
C     COMPLEX CY,Z
         EXTERNAL ZABS
         DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM, FN,
     *    FNU, FNUL, RL, R1M5, TOL, UFL, ZI, ZR, D1MACH, ZABS, BB
         INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, I1MACH
         DIMENSION CYR(N), CYI(N)
C***FIRST EXECUTABLE STATEMENT  ZBESK
         IERR = 0
         NZ=0
         IF (ZI.EQ.0.0E0 .AND. ZR.EQ.0.0E0) IERR=1
         IF (FNU.LT.0.0D0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
         TOL = DMAX1(D1MACH(4),1.0D-18)
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         R1M5 = D1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
         K1 = I1MACH(14) - 1
         AA = R1M5*DBLE(FLOAT(K1))
         DIG = DMIN1(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + DMAX1(-AA,-41.45D0)
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
         RL = 1.2D0*DIG + 3.0D0
C-----------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
         AZ = ZABS(ZR,ZI)
         FN = FNU + DBLE(FLOAT(NN-1))
         AA = 0.5D0/TOL
         BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
         AA = DMIN1(AA,BB)
         IF (AZ.GT.AA) GO TO 260
         IF (FN.GT.AA) GO TO 260
         AA = DSQRT(AA)
         IF (AZ.GT.AA) IERR=3
         IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
C     UFL = DEXP(-ELIM)
         UFL = D1MACH(1)*1.0D+3
         IF (AZ.LT.UFL) GO TO 180
         IF (FNU.GT.FNUL) GO TO 80
         IF (FN.LE.1.0D0) GO TO 60
         IF (FN.GT.2.0D0) GO TO 50
         IF (AZ.GT.TOL) GO TO 60
         ARG = 0.5D0*AZ
         ALN = -FN*DLOG(ARG)
         IF (ALN.GT.ELIM) GO TO 180
         GO TO 60
   50    CONTINUE
         CALL ZUOIK(ZR, ZI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     *    ALIM)
         IF (NUF.LT.0) GO TO 180
         NZ = NZ + NUF
         NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
         IF (NN.EQ.0) GO TO 100
   60    CONTINUE
         IF (ZR.LT.0.0D0) GO TO 70
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
C-----------------------------------------------------------------------
         CALL ZBKNU(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 200
         NZ=NW
         RETURN
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
C-----------------------------------------------------------------------
   70    CONTINUE
         IF (NZ.NE.0) GO TO 180
         MR = 1
         IF (ZI.LT.0.0D0) MR = -1
         CALL ZACON(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     *    TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 200
         NZ=NW
         RETURN
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
   80    CONTINUE
         MR = 0
         IF (ZR.GE.0.0D0) GO TO 90
         MR = 1
         IF (ZI.LT.0.0D0) MR = -1
   90    CONTINUE
         CALL ZBUNK(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 200
         NZ = NZ + NW
         RETURN
  100    CONTINUE
         IF (ZR.LT.0.0D0) GO TO 180
         RETURN
  180    CONTINUE
         NZ = 0
         IERR=2
         RETURN
  200    CONTINUE
         IF(NW.EQ.(-1)) GO TO 180
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         NZ=0
         IERR=4
         RETURN
      END
      SUBROUTINE ZBESY(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     *           CWRKI, IERR)

c*********************************************************************72
c
cc ZBESY computes a sequence of complex Bessel Y functions.
c
C***BEGIN PROLOGUE  ZBESY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
C             BESSEL FUNCTION OF SECOND KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C
C         ON KODE=1, ZBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESY RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT
C           CWRKI    AT LEAST N
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON KODE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON KODE=2)
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT IN TERMS OF THE I(FNU,Z) AND
C         K(FNU,Z) BESSEL FUNCTIONS IN THE RIGHT HALF PLANE BY
C
C             Y(FNU,Z) = I*CC*I(FNU,ARG) - (2/PI)*CONJG(CC)*K(FNU,ARG)
C
C             Y(FNU,Z) = CONJG(Y(FNU,CONJG(Z)))
C
C         FOR AIMAG(Z).GE.0 AND AIMAG(Z).LT.0 RESPECTIVELY, WHERE
C         CC=EXP(I*PI*FNU/2), ARG=Z*EXP(-I*PI/2) AND I**2=-1.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBESI,ZBESK,I1MACH,D1MACH
C***END PROLOGUE  ZBESY
C
C     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
         DOUBLE PRECISION ARG, ASCLE, CIPI, CIPR, CSGNI, CSGNR, CSPNI,
     *    CSPNR, CWRKI, CWRKR, CYI, CYR, D1M5, D1MACH, ELIM, EXI, EXR, EY,
     *    FNU, FFNU, HPI, RHPI, STR, STI, TAY, TOL, ATOL, RTOL, ZI, ZR,
     *    ZNI, ZNR, ZUI, ZUR, ZVI, ZVR, ZZI, ZZR
         INTEGER I, IERR, IFNU, I4, K, KODE, K1, K2, N, NZ, NZ1, NZ2,
     *    I1MACH
         DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N), CIPR(4), CIPI(4)
         DATA CIPR(1),CIPR(2),CIPR(3),CIPR(4)/1.0D0, 0.0D0, -1.0D0, 0.0D0/
         DATA CIPI(1),CIPI(2),CIPI(3),CIPI(4)/0.0D0, 1.0D0, 0.0D0, -1.0D0/
         DATA HPI / 1.57079632679489662D0 /
C***FIRST EXECUTABLE STATEMENT  ZBESY
         IERR = 0
         NZ=0
         IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
         IF (FNU.LT.0.0D0) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (N.LT.1) IERR=1
         IF (IERR.NE.0) RETURN
         ZZR = ZR
         ZZI = ZI
         IF (ZI.LT.0.0D0) ZZI = -ZZI
         ZNR = ZZI
         ZNI = -ZZR
         CALL ZBESI(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ1, IERR)
         IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
         CALL ZBESK(ZNR, ZNI, FNU, KODE, N, CWRKR, CWRKI, NZ2, IERR)
         IF (IERR.NE.0.AND.IERR.NE.3) GO TO 90
         NZ = MIN(NZ1,NZ2)
         IFNU = INT(SNGL(FNU))
         FFNU = FNU - DBLE(FLOAT(IFNU))
         ARG = HPI*FFNU
         CSGNR = COS(ARG)
         CSGNI = SIN(ARG)
         I4 = MOD(IFNU,4) + 1
         STR = CSGNR*CIPR(I4) - CSGNI*CIPI(I4)
         CSGNI = CSGNR*CIPI(I4) + CSGNI*CIPR(I4)
         CSGNR = STR
         RHPI = 1.0D0/HPI
         CSPNR = CSGNR*RHPI
         CSPNI = -CSGNI*RHPI
         STR = -CSGNI
         CSGNI = CSGNR
         CSGNR = STR
         IF (KODE.EQ.2) GO TO 60
         DO 50 I=1,N
C       CY(I) = CSGN*CY(I)-CSPN*CWRK(I)
            STR = CSGNR*CYR(I) - CSGNI*CYI(I)
            STR = STR - (CSPNR*CWRKR(I) - CSPNI*CWRKI(I))
            STI = CSGNR*CYI(I) + CSGNI*CYR(I)
            STI = STI - (CSPNR*CWRKI(I) + CSPNI*CWRKR(I))
            CYR(I) = STR
            CYI(I) = STI
            STR = - CSGNI
            CSGNI = CSGNR
            CSGNR = STR
            STR = CSPNI
            CSPNI = -CSPNR
            CSPNR = STR
   50    CONTINUE
         IF (ZI.LT.0.0D0) THEN
            DO 55 I=1,N
               CYI(I) = -CYI(I)
   55       CONTINUE
         ENDIF
         RETURN
   60    CONTINUE
         EXR = COS(ZR)
         EXI = SIN(ZR)
         TOL = MAX(D1MACH(4),1.0D-18)
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         K = MIN(IABS(K1),IABS(K2))
         D1M5 = D1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
         ELIM = 2.303D0*(DBLE(FLOAT(K))*D1M5-3.0D0)
         EY = 0.0D0
         TAY = ABS(ZI+ZI)
         IF (TAY.LT.ELIM) EY = EXP(-TAY)
         STR = (EXR*CSPNR - EXI*CSPNI)*EY
         CSPNI = (EXR*CSPNI + EXI*CSPNR)*EY
         CSPNR = STR
         NZ = 0
         RTOL = 1.0D0/TOL
         ASCLE = D1MACH(1)*RTOL*1.0D+3
         DO 80 I=1,N
C----------------------------------------------------------------------
C       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN
C       SCALED MODE IF CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO
C       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.
C----------------------------------------------------------------------
            ZVR = CWRKR(I)
            ZVI = CWRKI(I)
            ATOL=1.0D0
            IF (MAX(ABS(ZVR),ABS(ZVI)).GT.ASCLE) GO TO 75
            ZVR = ZVR*RTOL
            ZVI = ZVI*RTOL
            ATOL = TOL
   75       CONTINUE
            STR = (ZVR*CSPNR - ZVI*CSPNI)*ATOL
            ZVI = (ZVR*CSPNI + ZVI*CSPNR)*ATOL
            ZVR = STR
            ZUR = CYR(I)
            ZUI = CYI(I)
            ATOL=1.0D0
            IF (MAX(ABS(ZUR),ABS(ZUI)).GT.ASCLE) GO TO 85
            ZUR = ZUR*RTOL
            ZUI = ZUI*RTOL
            ATOL = TOL
   85       CONTINUE
            STR = (ZUR*CSGNR - ZUI*CSGNI)*ATOL
            ZUI = (ZUR*CSGNI + ZUI*CSGNR)*ATOL
            ZUR = STR
            CYR(I) = ZUR - ZVR
            CYI(I) = ZUI - ZVI
            IF (ZI.LT.0.0D0) CYI(I) = -CYI(I)
            IF (CYR(I).EQ.0.0D0 .AND. CYI(I).EQ.0.0D0 .AND. EY.EQ.0.0D0)
     $          NZ = NZ + 1
            STR = -CSGNI
            CSGNI = CSGNR
            CSGNR = STR
            STR = CSPNI
            CSPNI = -CSPNR
            CSPNR = STR
   80    CONTINUE
         RETURN
   90    CONTINUE
         NZ = 0
         RETURN
      END
      SUBROUTINE ZBINU(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL,
     * TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZBINU
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY
C
C     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  ZABS,ZASYI,ZBUNI,ZMLRI,ZSERI,ZUOIK,ZWRSK
C***END PROLOGUE  ZBINU
         EXTERNAL ZABS
         DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU,
     *    FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR, ZABS
         INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
         DIMENSION CYR(N), CYI(N), CWR(2), CWI(2)
         DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
C
         NZ = 0
         AZ = ZABS(ZR,ZI)
         NN = N
         DFNU = FNU + DBLE(FLOAT(N-1))
         IF (AZ.LE.2.0D0) GO TO 10
         IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
         CALL ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
         INW = IABS(NW)
         NZ = NZ + INW
         NN = NN - INW
         IF (NN.EQ.0) RETURN
         IF (NW.GE.0) GO TO 120
         DFNU = FNU + DBLE(FLOAT(NN-1))
   20    CONTINUE
         IF (AZ.LT.RL) GO TO 40
         IF (DFNU.LE.1.0D0) GO TO 30
         IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30    CONTINUE
         CALL ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 130
         GO TO 120
   40    CONTINUE
         IF (DFNU.LE.1.0D0) GO TO 70
   50    CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
         CALL ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM,
     *    ALIM)
         IF (NW.LT.0) GO TO 130
         NZ = NZ + NW
         NN = NN - NW
         IF (NN.EQ.0) RETURN
         DFNU = FNU+DBLE(FLOAT(NN-1))
         IF (DFNU.GT.FNUL) GO TO 110
         IF (AZ.GT.FNUL) GO TO 110
   60    CONTINUE
         IF (AZ.GT.RL) GO TO 80
   70    CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
         CALL ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
         IF(NW.LT.0) GO TO 130
         GO TO 120
   80    CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
         CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM,
     *    ALIM)
         IF (NW.GE.0) GO TO 100
         NZ = NN
         DO 90 I=1,NN
            CYR(I) = ZEROR
            CYI(I) = ZEROI
   90    CONTINUE
         RETURN
  100    CONTINUE
         IF (NW.GT.0) GO TO 130
         CALL ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL,
     *    ELIM, ALIM)
         IF (NW.LT.0) GO TO 130
         GO TO 120
  110    CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
         NUI = INT(SNGL(FNUL-DFNU)) + 1
         NUI = MAX0(NUI,0)
         CALL ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL,
     *    TOL, ELIM, ALIM)
         IF (NW.LT.0) GO TO 130
         NZ = NZ + NW
         IF (NLAST.EQ.0) GO TO 120
         NN = NLAST
         GO TO 60
  120    CONTINUE
         RETURN
  130    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
      SUBROUTINE ZBIRY(ZR, ZI, ID, KODE, BIR, BII, IERR)

c*********************************************************************72
c
cc ZBIRY computes a sequence of complex Airy Bi functions.
c
C***BEGIN PROLOGUE  ZBIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801, 930101   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
C         ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
C         DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
C         BOTH THE LEFT AND RIGHT HALF PLANES WHERE
C         ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
C         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             BI=BI(Z)                 ON ID=0 OR
C                             BI=DBI(Z)/DZ             ON ID=1
C                        = 2  RETURNS
C                             BI=CEXP(-AXZTA)*BI(Z)     ON ID=0 OR
C                             BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
C                             AND AXZTA=ABS(XZTA)
C
C         OUTPUT     BIR,BII ARE DOUBLE PRECISION
C           BIR,BII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         BI AND DBI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE I BESSEL
C         FUNCTIONS BY
C
C                BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
C               DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
C                               C=1.0/SQRT(3.0)
C                             ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
C                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
C                 PP 265-273.
C
C***ROUTINES CALLED  ZBINU,ZABS,ZDIV,ZSQRT,D1MACH,I1MACH
C***END PROLOGUE  ZBIRY
C     COMPLEX BI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
         EXTERNAL ZABS
         DOUBLE PRECISION AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BII, BIR,
     *    BK, CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2,
     *    DIG, DK, D1, D2, EAA, ELIM, FID, FMR, FNU, FNUL, PI, RL, R1M5,
     *    SFAC, STI, STR, S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I,
     *    TRM2R, TTH, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS
         INTEGER ID, IERR, K, KODE, K1, K2, NZ, I1MACH
         DIMENSION CYR(2), CYI(2)
         DATA TTH, C1, C2, COEF, PI /6.66666666666666667D-01,
     *    6.14926627446000736D-01,4.48288357353826359D-01,
     *    5.77350269189625765D-01,3.14159265358979324D+00/
         DATA CONER, CONEI /1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZBIRY
         IERR = 0
         NZ=0
         IF (ID.LT.0 .OR. ID.GT.1) IERR=1
         IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
         IF (IERR.NE.0) RETURN
         AZ = ZABS(ZR,ZI)
         TOL = DMAX1(D1MACH(4),1.0D-18)
         FID = DBLE(FLOAT(ID))
         IF (AZ.GT.1.0E0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
         S1R = CONER
         S1I = CONEI
         S2R = CONER
         S2I = CONEI
         IF (AZ.LT.TOL) GO TO 130
         AA = AZ*AZ
         IF (AA.LT.TOL/AZ) GO TO 40
         TRM1R = CONER
         TRM1I = CONEI
         TRM2R = CONER
         TRM2I = CONEI
         ATRM = 1.0D0
         STR = ZR*ZR - ZI*ZI
         STI = ZR*ZI + ZI*ZR
         Z3R = STR*ZR - STI*ZI
         Z3I = STR*ZI + STI*ZR
         AZ3 = AZ*AA
         AK = 2.0D0 + FID
         BK = 3.0D0 - FID - FID
         CK = 4.0D0 - FID
         DK = 3.0D0 + FID + FID
         D1 = AK*DK
         D2 = BK*CK
         AD = DMIN1(D1,D2)
         AK = 24.0D0 + 9.0D0*FID
         BK = 30.0D0 - 9.0D0*FID
         DO 30 K=1,25
            STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
            TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
            TRM1R = STR
            S1R = S1R + TRM1R
            S1I = S1I + TRM1I
            STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
            TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
            TRM2R = STR
            S2R = S2R + TRM2R
            S2I = S2I + TRM2I
            ATRM = ATRM*AZ3/AD
            D1 = D1 + AK
            D2 = D2 + BK
            AD = DMIN1(D1,D2)
            IF (ATRM.LT.TOL*AD) GO TO 40
            AK = AK + 18.0D0
            BK = BK + 18.0D0
   30    CONTINUE
   40    CONTINUE
         IF (ID.EQ.1) GO TO 50
         BIR = C1*S1R + C2*(ZR*S2R-ZI*S2I)
         BII = C1*S1I + C2*(ZR*S2I+ZI*S2R)
         IF (KODE.EQ.1) RETURN
         CALL ZSQRT(ZR, ZI, STR, STI)
         ZTAR = TTH*(ZR*STR-ZI*STI)
         ZTAI = TTH*(ZR*STI+ZI*STR)
         AA = ZTAR
         AA = -DABS(AA)
         EAA = DEXP(AA)
         BIR = BIR*EAA
         BII = BII*EAA
         RETURN
   50    CONTINUE
         BIR = S2R*C2
         BII = S2I*C2
         IF (AZ.LE.TOL) GO TO 60
         CC = C1/(1.0D0+FID)
         STR = S1R*ZR - S1I*ZI
         STI = S1R*ZI + S1I*ZR
         BIR = BIR + CC*(STR*ZR-STI*ZI)
         BII = BII + CC*(STR*ZI+STI*ZR)
   60    CONTINUE
         IF (KODE.EQ.1) RETURN
         CALL ZSQRT(ZR, ZI, STR, STI)
         ZTAR = TTH*(ZR*STR-ZI*STI)
         ZTAI = TTH*(ZR*STI+ZI*STR)
         AA = ZTAR
         AA = -DABS(AA)
         EAA = DEXP(AA)
         BIR = BIR*EAA
         BII = BII*EAA
         RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70    CONTINUE
         FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
         K1 = I1MACH(15)
         K2 = I1MACH(16)
         R1M5 = D1MACH(5)
         K = MIN0(IABS(K1),IABS(K2))
         ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
         K1 = I1MACH(14) - 1
         AA = R1M5*DBLE(FLOAT(K1))
         DIG = DMIN1(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + DMAX1(-AA,-41.45D0)
         RL = 1.2D0*DIG + 3.0D0
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
         AA=0.5D0/TOL
         BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
         AA=DMIN1(AA,BB)
         AA=AA**TTH
         IF (AZ.GT.AA) GO TO 260
         AA=DSQRT(AA)
         IF (AZ.GT.AA) IERR=3
         CALL ZSQRT(ZR, ZI, CSQR, CSQI)
         ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
         ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
         SFAC = 1.0D0
         AK = ZTAI
         IF (ZR.GE.0.0D0) GO TO 80
         BK = ZTAR
         CK = -DABS(BK)
         ZTAR = CK
         ZTAI = AK
   80    CONTINUE
         IF (ZI.NE.0.0D0 .OR. ZR.GT.0.0D0) GO TO 90
         ZTAR = 0.0D0
         ZTAI = AK
   90    CONTINUE
         AA = ZTAR
         IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         BB = DABS(AA)
         IF (BB.LT.ALIM) GO TO 100
         BB = BB + 0.25D0*DLOG(AZ)
         SFAC = TOL
         IF (BB.GT.ELIM) GO TO 190
  100    CONTINUE
         FMR = 0.0D0
         IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
         FMR = PI
         IF (ZI.LT.0.0D0) FMR = -PI
         ZTAR = -ZTAR
         ZTAI = -ZTAI
  110    CONTINUE
C-----------------------------------------------------------------------
C     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
C     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM ZBESI
C-----------------------------------------------------------------------
         CALL ZBINU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, RL, FNUL, TOL,
     *    ELIM, ALIM)
         IF (NZ.LT.0) GO TO 200
         AA = FMR*FNU
         Z3R = SFAC
         STR = DCOS(AA)
         STI = DSIN(AA)
         S1R = (STR*CYR(1)-STI*CYI(1))*Z3R
         S1I = (STR*CYI(1)+STI*CYR(1))*Z3R
         FNU = (2.0D0-FID)/3.0D0
         CALL ZBINU(ZTAR, ZTAI, FNU, KODE, 2, CYR, CYI, NZ, RL, FNUL, TOL,
     *    ELIM, ALIM)
         CYR(1) = CYR(1)*Z3R
         CYI(1) = CYI(1)*Z3R
         CYR(2) = CYR(2)*Z3R
         CYI(2) = CYI(2)*Z3R
C-----------------------------------------------------------------------
C     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
C-----------------------------------------------------------------------
         CALL ZDIV(CYR(1), CYI(1), ZTAR, ZTAI, STR, STI)
         S2R = (FNU+FNU)*STR + CYR(2)
         S2I = (FNU+FNU)*STI + CYI(2)
         AA = FMR*(FNU-1.0D0)
         STR = DCOS(AA)
         STI = DSIN(AA)
         S1R = COEF*(S1R+S2R*STR-S2I*STI)
         S1I = COEF*(S1I+S2R*STI+S2I*STR)
         IF (ID.EQ.1) GO TO 120
         STR = CSQR*S1R - CSQI*S1I
         S1I = CSQR*S1I + CSQI*S1R
         S1R = STR
         BIR = S1R/SFAC
         BII = S1I/SFAC
         RETURN
  120    CONTINUE
         STR = ZR*S1R - ZI*S1I
         S1I = ZR*S1I + ZI*S1R
         S1R = STR
         BIR = S1R/SFAC
         BII = S1I/SFAC
         RETURN
  130    CONTINUE
         AA = C1*(1.0D0-FID) + FID*C2
         BIR = AA
         BII = 0.0D0
         RETURN
  190    CONTINUE
         IERR=2
         NZ=0
         RETURN
  200    CONTINUE
         IF(NZ.EQ.(-1)) GO TO 190
         NZ=0
         IERR=5
         RETURN
  260    CONTINUE
         IERR=4
         NZ=0
         RETURN
      END
      SUBROUTINE ZBKNU(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
cc ZBKNU computes the K Bessel function in the right half Z plane.
c
C***BEGIN PROLOGUE  ZBKNU
C***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH
C
C     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
C
C***ROUTINES CALLED  DGAMLN,I1MACH,D1MACH,ZKSCL,ZSHCH,ZUCHK,ZABS,ZDIV,
C                    ZEXP,ZLOG,ZMLT,ZSQRT
C***END PROLOGUE  ZBKNU
C
         EXTERNAL ZABS
         DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ,
     *    CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER,
     *    CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOR,
     *    CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS,
     *    FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, HPI, PI, PR, PTI,
     *    PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI,
     *    RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,
     *    TOL, TTH, T1, T2, YI, YR, ZI, ZR, DGAMLN, D1MACH, ZABS, ELM,
     *    CELMR, ZDR, ZDI, AS, ALAS, HELIM, CYR, CYI
         INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ,
     *    IDUM, I1MACH, J, IC, INUB, NW
         DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2),
     *    CYI(2)
C     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
C     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
C
         DATA KMAX / 30 /
         DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,R1/
     1     0.0D0 , 0.0D0 , 1.0D0 , 0.0D0 , 2.0D0 , 2.0D0 /
         DATA DPI, RTHPI, SPI ,HPI, FPI, TTH /
     1        3.14159265358979324D0,       1.25331413731550025D0,
     2        1.90985931710274403D0,       1.57079632679489662D0,
     3        1.89769999331517738D0,       6.66666666666666666D-01/
         DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1        5.77215664901532861D-01,    -4.20026350340952355D-02,
     2       -4.21977345555443367D-02,     7.21894324666309954D-03,
     3       -2.15241674114950973D-04,    -2.01348547807882387D-05,
     4        1.13302723198169588D-06,     6.11609510448141582D-09/
C
         CAZ = ZABS(ZR,ZI)
         CSCLR = 1.0D0/TOL
         CRSCR = TOL
         CSSR(1) = CSCLR
         CSSR(2) = 1.0D0
         CSSR(3) = CRSCR
         CSRR(1) = CRSCR
         CSRR(2) = 1.0D0
         CSRR(3) = CSCLR
         BRY(1) = 1.0D+3*D1MACH(1)/TOL
         BRY(2) = 1.0D0/BRY(1)
         BRY(3) = D1MACH(2)
         NZ = 0
         IFLAG = 0
         KODED = KODE
         RCAZ = 1.0D0/CAZ
         STR = ZR*RCAZ
         STI = -ZI*RCAZ
         RZR = (STR+STR)*RCAZ
         RZI = (STI+STI)*RCAZ
         INU = INT(FNU+0.5D0)
         DNU = FNU - DBLE(FLOAT(INU))
         IF (DABS(DNU).EQ.0.5D0) GO TO 110
         DNU2 = 0.0D0
         IF (DABS(DNU).GT.TOL) DNU2 = DNU*DNU
         IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR CABS(Z).LE.R1
C-----------------------------------------------------------------------
         FC = 1.0D0
         CALL ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
         FMUR = SMUR*DNU
         FMUI = SMUI*DNU
         CALL ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
         IF (DNU.EQ.0.0D0) GO TO 10
         FC = DNU*DPI
         FC = FC/DSIN(FC)
         SMUR = CSHR/DNU
         SMUI = CSHI/DNU
   10    CONTINUE
         A2 = 1.0D0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
         T2 = DEXP(-DGAMLN(A2,IDUM))
         T1 = 1.0D0/(T2*FC)
         IF (DABS(DNU).GT.0.1D0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
         AK = 1.0D0
         S = CC(1)
         DO 20 K=2,8
            AK = AK*DNU2
            TM = CC(K)*AK
            S = S + TM
            IF (DABS(TM).LT.TOL) GO TO 30
   20    CONTINUE
   30    G1 = -S
         GO TO 50
   40    CONTINUE
         G1 = (T1-T2)/(DNU+DNU)
   50    CONTINUE
         G2 = (T1+T2)*0.5D0
         FR = FC*(CCHR*G1+SMUR*G2)
         FI = FC*(CCHI*G1+SMUI*G2)
         CALL ZEXP(FMUR, FMUI, STR, STI)
         PR = 0.5D0*STR/T2
         PI = 0.5D0*STI/T2
         CALL ZDIV(0.5D0, 0.0D0, STR, STI, PTR, PTI)
         QR = PTR/T1
         QI = PTI/T1
         S1R = FR
         S1I = FI
         S2R = PR
         S2I = PI
         AK = 1.0D0
         A1 = 1.0D0
         CKR = CONER
         CKI = CONEI
         BK = 1.0D0 - DNU2
         IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
         IF (CAZ.LT.TOL) GO TO 70
         CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
         CZR = 0.25D0*CZR
         CZI = 0.25D0*CZI
         T1 = 0.25D0*CAZ*CAZ
   60    CONTINUE
         FR = (FR*AK+PR+QR)/BK
         FI = (FI*AK+PI+QI)/BK
         STR = 1.0D0/(AK-DNU)
         PR = PR*STR
         PI = PI*STR
         STR = 1.0D0/(AK+DNU)
         QR = QR*STR
         QI = QI*STR
         STR = CKR*CZR - CKI*CZI
         RAK = 1.0D0/AK
         CKI = (CKR*CZI+CKI*CZR)*RAK
         CKR = STR*RAK
         S1R = CKR*FR - CKI*FI + S1R
         S1I = CKR*FI + CKI*FR + S1I
         A1 = A1*T1*RAK
         BK = BK + AK + AK + 1.0D0
         AK = AK + 1.0D0
         IF (A1.GT.TOL) GO TO 60
   70    CONTINUE
         YR(1) = S1R
         YI(1) = S1I
         IF (KODED.EQ.1) RETURN
         CALL ZEXP(ZR, ZI, STR, STI)
         CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
         RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80    CONTINUE
         IF (CAZ.LT.TOL) GO TO 100
         CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
         CZR = 0.25D0*CZR
         CZI = 0.25D0*CZI
         T1 = 0.25D0*CAZ*CAZ
   90    CONTINUE
         FR = (FR*AK+PR+QR)/BK
         FI = (FI*AK+PI+QI)/BK
         STR = 1.0D0/(AK-DNU)
         PR = PR*STR
         PI = PI*STR
         STR = 1.0D0/(AK+DNU)
         QR = QR*STR
         QI = QI*STR
         STR = CKR*CZR - CKI*CZI
         RAK = 1.0D0/AK
         CKI = (CKR*CZI+CKI*CZR)*RAK
         CKR = STR*RAK
         S1R = CKR*FR - CKI*FI + S1R
         S1I = CKR*FI + CKI*FR + S1I
         STR = PR - FR*AK
         STI = PI - FI*AK
         S2R = CKR*STR - CKI*STI + S2R
         S2I = CKR*STI + CKI*STR + S2I
         A1 = A1*T1*RAK
         BK = BK + AK + AK + 1.0D0
         AK = AK + 1.0D0
         IF (A1.GT.TOL) GO TO 90
  100    CONTINUE
         KFLAG = 2
         A1 = FNU + 1.0D0
         AK = A1*DABS(SMUR)
         IF (AK.GT.ALIM) KFLAG = 3
         STR = CSSR(KFLAG)
         P2R = S2R*STR
         P2I = S2I*STR
         CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
         S1R = S1R*STR
         S1I = S1I*STR
         IF (KODED.EQ.1) GO TO 210
         CALL ZEXP(ZR, ZI, FR, FI)
         CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I)
         CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I)
         GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110    CONTINUE
         CALL ZSQRT(ZR, ZI, STR, STI)
         CALL ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI)
         KFLAG = 2
         IF (KODED.EQ.2) GO TO 120
         IF (ZR.GT.ALIM) GO TO 290
C     BLANK LINE
         STR = DEXP(-ZR)*CSSR(KFLAG)
         STI = -STR*DSIN(ZI)
         STR = STR*DCOS(ZI)
         CALL ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120    CONTINUE
         IF (DABS(DNU).EQ.0.5D0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR CABS(Z).GT.R1
C-----------------------------------------------------------------------
         AK = DCOS(DPI*DNU)
         AK = DABS(AK)
         IF (AK.EQ.CZEROR) GO TO 300
         FHS = DABS(0.25D0-DNU2)
         IF (FHS.EQ.CZEROR) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
         T1 = DBLE(FLOAT(I1MACH(14)-1))
         T1 = T1*D1MACH(5)*3.321928094D0
         T1 = DMAX1(T1,12.0D0)
         T1 = DMIN1(T1,60.0D0)
         T2 = TTH*T1 - 6.0D0
         IF (ZR.NE.0.0D0) GO TO 130
         T1 = HPI
         GO TO 140
  130    CONTINUE
         T1 = DATAN(ZI/ZR)
         T1 = DABS(T1)
  140    CONTINUE
         IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
C-----------------------------------------------------------------------
         ETEST = AK/(DPI*CAZ*TOL)
         FK = CONER
         IF (ETEST.LT.CONER) GO TO 180
         FKS = CTWOR
         CKR = CAZ + CAZ + CTWOR
         P1R = CZEROR
         P2R = CONER
         DO 150 I=1,KMAX
            AK = FHS/FKS
            CBR = CKR/(FK+CONER)
            PTR = P2R
            P2R = CBR*P2R - P1R*AK
            P1R = PTR
            CKR = CKR + CTWOR
            FKS = FKS + FK + FK + CTWOR
            FHS = FHS + FK + FK
            FK = FK + CONER
            STR = DABS(P2R)*FK
            IF (ETEST.LT.STR) GO TO 160
  150    CONTINUE
         GO TO 310
  160    CONTINUE
         FK = FK + SPI*T1*DSQRT(T2/CAZ)
         FHS = DABS(0.25D0-DNU2)
         GO TO 180
  170    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
C-----------------------------------------------------------------------
         A2 = DSQRT(CAZ)
         AK = FPI*AK/(TOL*DSQRT(A2))
         AA = 3.0D0*T1/(1.0D0+CAZ)
         BB = 14.7D0*T1/(28.0D0+CAZ)
         AK = (DLOG(AK)+CAZ*DCOS(AA)/(1.0D0+0.008D0*CAZ))/DCOS(BB)
         FK = 0.12125D0*AK*AK/CAZ + 1.5D0
  180    CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
         K = INT(SNGL(FK))
         FK = DBLE(FLOAT(K))
         FKS = FK*FK
         P1R = CZEROR
         P1I = CZEROI
         P2R = TOL
         P2I = CZEROI
         CSR = P2R
         CSI = P2I
         DO 190 I=1,K
            A1 = FKS - FK
            AK = (FKS+FK)/(A1+FHS)
            RAK = 2.0D0/(FK+CONER)
            CBR = (FK+ZR)*RAK
            CBI = ZI*RAK
            PTR = P2R
            PTI = P2I
            P2R = (PTR*CBR-PTI*CBI-P1R)*AK
            P2I = (PTI*CBR+PTR*CBI-P1I)*AK
            P1R = PTR
            P1I = PTI
            CSR = CSR + P2R
            CSI = CSI + P2I
            FKS = A1 - FK + CONER
            FK = FK - CONER
  190    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
         TM = ZABS(CSR,CSI)
         PTR = 1.0D0/TM
         S1R = P2R*PTR
         S1I = P2I*PTR
         CSR = CSR*PTR
         CSI = -CSI*PTR
         CALL ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
         CALL ZMLT(STR, STI, CSR, CSI, S1R, S1I)
         IF (INU.GT.0 .OR. N.GT.1) GO TO 200
         ZDR = ZR
         ZDI = ZI
         IF(IFLAG.EQ.1) GO TO 270
         GO TO 240
  200    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
C-----------------------------------------------------------------------
         TM = ZABS(P2R,P2I)
         PTR = 1.0D0/TM
         P1R = P1R*PTR
         P1I = P1I*PTR
         P2R = P2R*PTR
         P2I = -P2I*PTR
         CALL ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
         STR = DNU + 0.5D0 - PTR
         STI = -PTI
         CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
         STR = STR + 1.0D0
         CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210    CONTINUE
         STR = DNU + 1.0D0
         CKR = STR*RZR
         CKI = STR*RZI
         IF (N.EQ.1) INU = INU - 1
         IF (INU.GT.0) GO TO 220
         IF (N.GT.1) GO TO 215
         S1R = S2R
         S1I = S2I
  215    CONTINUE
         ZDR = ZR
         ZDI = ZI
         IF(IFLAG.EQ.1) GO TO 270
         GO TO 240
  220    CONTINUE
         INUB = 1
         IF(IFLAG.EQ.1) GO TO 261
  225    CONTINUE
         P1R = CSRR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 230 I=INUB,INU
            STR = S2R
            STI = S2I
            S2R = CKR*STR - CKI*STI + S1R
            S2I = CKR*STI + CKI*STR + S1I
            S1R = STR
            S1I = STI
            CKR = CKR + RZR
            CKI = CKI + RZI
            IF (KFLAG.GE.3) GO TO 230
            P2R = S2R*P1R
            P2I = S2I*P1R
            STR = DABS(P2R)
            STI = DABS(P2I)
            P2M = DMAX1(STR,STI)
            IF (P2M.LE.ASCLE) GO TO 230
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1R = S1R*P1R
            S1I = S1I*P1R
            S2R = P2R
            S2I = P2I
            STR = CSSR(KFLAG)
            S1R = S1R*STR
            S1I = S1I*STR
            S2R = S2R*STR
            S2I = S2I*STR
            P1R = CSRR(KFLAG)
  230    CONTINUE
         IF (N.NE.1) GO TO 240
         S1R = S2R
         S1I = S2I
  240    CONTINUE
         STR = CSRR(KFLAG)
         YR(1) = S1R*STR
         YI(1) = S1I*STR
         IF (N.EQ.1) RETURN
         YR(2) = S2R*STR
         YI(2) = S2I*STR
         IF (N.EQ.2) RETURN
         KK = 2
  250    CONTINUE
         KK = KK + 1
         IF (KK.GT.N) RETURN
         P1R = CSRR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 260 I=KK,N
            P2R = S2R
            P2I = S2I
            S2R = CKR*P2R - CKI*P2I + S1R
            S2I = CKI*P2R + CKR*P2I + S1I
            S1R = P2R
            S1I = P2I
            CKR = CKR + RZR
            CKI = CKI + RZI
            P2R = S2R*P1R
            P2I = S2I*P1R
            YR(I) = P2R
            YI(I) = P2I
            IF (KFLAG.GE.3) GO TO 260
            STR = DABS(P2R)
            STI = DABS(P2I)
            P2M = DMAX1(STR,STI)
            IF (P2M.LE.ASCLE) GO TO 260
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1R = S1R*P1R
            S1I = S1I*P1R
            S2R = P2R
            S2I = P2I
            STR = CSSR(KFLAG)
            S1R = S1R*STR
            S1I = S1I*STR
            S2R = S2R*STR
            S2I = S2I*STR
            P1R = CSRR(KFLAG)
  260    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261    CONTINUE
         HELIM = 0.5D0*ELIM
         ELM = DEXP(-ELIM)
         CELMR = ELM
         ASCLE = BRY(1)
         ZDR = ZR
         ZDI = ZI
         IC = -1
         J = 2
         DO 262 I=1,INU
            STR = S2R
            STI = S2I
            S2R = STR*CKR-STI*CKI+S1R
            S2I = STI*CKR+STR*CKI+S1I
            S1R = STR
            S1I = STI
            CKR = CKR+RZR
            CKI = CKI+RZI
            AS = ZABS(S2R,S2I)
            ALAS = DLOG(AS)
            P2R = -ZDR+ALAS
            IF(P2R.LT.(-ELIM)) GO TO 263
            CALL ZLOG(S2R,S2I,STR,STI,IDUM)
            P2R = -ZDR+STR
            P2I = -ZDI+STI
            P2M = DEXP(P2R)/TOL
            P1R = P2M*DCOS(P2I)
            P1I = P2M*DSIN(P2I)
            CALL ZUCHK(P1R,P1I,NW,ASCLE,TOL)
            IF(NW.NE.0) GO TO 263
            J = 3 - J
            CYR(J) = P1R
            CYI(J) = P1I
            IF(IC.EQ.(I-1)) GO TO 264
            IC = I
            GO TO 262
  263       CONTINUE
            IF(ALAS.LT.HELIM) GO TO 262
            ZDR = ZDR-ELIM
            S1R = S1R*CELMR
            S1I = S1I*CELMR
            S2R = S2R*CELMR
            S2I = S2I*CELMR
  262    CONTINUE
         IF(N.NE.1) GO TO 270
         S1R = S2R
         S1I = S2I
         GO TO 270
  264    CONTINUE
         KFLAG = 1
         INUB = I+1
         S2R = CYR(J)
         S2I = CYI(J)
         J = 3 - J
         S1R = CYR(J)
         S1I = CYI(J)
         IF(INUB.LE.INU) GO TO 225
         IF(N.NE.1) GO TO 240
         S1R = S2R
         S1I = S2I
         GO TO 240
  270    CONTINUE
         YR(1) = S1R
         YI(1) = S1I
         IF(N.EQ.1) GO TO 280
         YR(2) = S2R
         YI(2) = S2I
  280    CONTINUE
         ASCLE = BRY(1)
         CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
         INU = N - NZ
         IF (INU.LE.0) RETURN
         KK = NZ + 1
         S1R = YR(KK)
         S1I = YI(KK)
         YR(KK) = S1R*CSRR(1)
         YI(KK) = S1I*CSRR(1)
         IF (INU.EQ.1) RETURN
         KK = NZ + 2
         S2R = YR(KK)
         S2I = YI(KK)
         YR(KK) = S2R*CSRR(1)
         YI(KK) = S2I*CSRR(1)
         IF (INU.EQ.2) RETURN
         T2 = FNU + DBLE(FLOAT(KK-1))
         CKR = T2*RZR
         CKI = T2*RZI
         KFLAG = 1
         GO TO 250
  290    CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY DEXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
         KODED = 2
         IFLAG = 1
         KFLAG = 2
         GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300    CONTINUE
         S1R = COEFR
         S1I = COEFI
         S2R = COEFR
         S2I = COEFI
         GO TO 210
C
C
  310    CONTINUE
         NZ=-2
         RETURN
      END
      SUBROUTINE ZBUNI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST,
     * FNUL, TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZBUNI
C***REFER TO  ZBESI,ZBESK
C
C     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***ROUTINES CALLED  ZUNI1,ZUNI2,ZABS,D1MACH
C***END PROLOGUE  ZBUNI
C     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
         EXTERNAL ZABS
         DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,
     *    ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,
     *    S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M,
     *    D1MACH
         INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
         DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
         NZ = 0
         AX = DABS(ZR)*1.7321D0
         AY = DABS(ZI)
         IFORM = 1
         IF (AY.GT.AX) IFORM = 2
         IF (NUI.EQ.0) GO TO 60
         FNUI = DBLE(FLOAT(NUI))
         DFNU = FNU + DBLE(FLOAT(N-1))
         GNU = DFNU + FNUI
         IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
         CALL ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     *    ELIM, ALIM)
         GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
         CALL ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     *    ELIM, ALIM)
   20    CONTINUE
         IF (NW.LT.0) GO TO 50
         IF (NW.NE.0) GO TO 90
         STR = ZABS(CYR(1),CYI(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
         BRY(1)=1.0D+3*D1MACH(1)/TOL
         BRY(2) = 1.0D0/BRY(1)
         BRY(3) = BRY(2)
         IFLAG = 2
         ASCLE = BRY(2)
         CSCLR = 1.0D0
         IF (STR.GT.BRY(1)) GO TO 21
         IFLAG = 1
         ASCLE = BRY(1)
         CSCLR = 1.0D0/TOL
         GO TO 25
   21    CONTINUE
         IF (STR.LT.BRY(2)) GO TO 25
         IFLAG = 3
         ASCLE=BRY(3)
         CSCLR = TOL
   25    CONTINUE
         CSCRR = 1.0D0/CSCLR
         S1R = CYR(2)*CSCLR
         S1I = CYI(2)*CSCLR
         S2R = CYR(1)*CSCLR
         S2I = CYI(1)*CSCLR
         RAZ = 1.0D0/ZABS(ZR,ZI)
         STR = ZR*RAZ
         STI = -ZI*RAZ
         RZR = (STR+STR)*RAZ
         RZI = (STI+STI)*RAZ
         DO 30 I=1,NUI
            STR = S2R
            STI = S2I
            S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
            S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
            S1R = STR
            S1I = STI
            FNUI = FNUI - 1.0D0
            IF (IFLAG.GE.3) GO TO 30
            STR = S2R*CSCRR
            STI = S2I*CSCRR
            C1R = DABS(STR)
            C1I = DABS(STI)
            C1M = DMAX1(C1R,C1I)
            IF (C1M.LE.ASCLE) GO TO 30
            IFLAG = IFLAG+1
            ASCLE = BRY(IFLAG)
            S1R = S1R*CSCRR
            S1I = S1I*CSCRR
            S2R = STR
            S2I = STI
            CSCLR = CSCLR*TOL
            CSCRR = 1.0D0/CSCLR
            S1R = S1R*CSCLR
            S1I = S1I*CSCLR
            S2R = S2R*CSCLR
            S2I = S2I*CSCLR
   30    CONTINUE
         YR(N) = S2R*CSCRR
         YI(N) = S2I*CSCRR
         IF (N.EQ.1) RETURN
         NL = N - 1
         FNUI = DBLE(FLOAT(NL))
         K = NL
         DO 40 I=1,NL
            STR = S2R
            STI = S2I
            S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
            S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
            S1R = STR
            S1I = STI
            STR = S2R*CSCRR
            STI = S2I*CSCRR
            YR(K) = STR
            YI(K) = STI
            FNUI = FNUI - 1.0D0
            K = K - 1
            IF (IFLAG.GE.3) GO TO 40
            C1R = DABS(STR)
            C1I = DABS(STI)
            C1M = DMAX1(C1R,C1I)
            IF (C1M.LE.ASCLE) GO TO 40
            IFLAG = IFLAG+1
            ASCLE = BRY(IFLAG)
            S1R = S1R*CSCRR
            S1I = S1I*CSCRR
            S2R = STR
            S2I = STI
            CSCLR = CSCLR*TOL
            CSCRR = 1.0D0/CSCLR
            S1R = S1R*CSCLR
            S1I = S1I*CSCLR
            S2R = S2R*CSCLR
            S2I = S2I*CSCLR
   40    CONTINUE
         RETURN
   50    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
   60    CONTINUE
         IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
         CALL ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     *    ELIM, ALIM)
         GO TO 80
   70    CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
         CALL ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     *    ELIM, ALIM)
   80    CONTINUE
         IF (NW.LT.0) GO TO 50
         NZ = NW
         RETURN
   90    CONTINUE
         NLAST = N
         RETURN
      END
      SUBROUTINE ZBUNK(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZBUNK
C***REFER TO  ZBESK,ZBESH
C
C     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
C
C***ROUTINES CALLED  ZUNK1,ZUNK2
C***END PROLOGUE  ZBUNK
C     COMPLEX Y,Z
         DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
         INTEGER KODE, MR, N, NZ
         DIMENSION YR(N), YI(N)
         NZ = 0
         AX = DABS(ZR)*1.7321D0
         AY = DABS(ZI)
         IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
         CALL ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
         GO TO 20
   10    CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
         CALL ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20    CONTINUE
         RETURN
      END
      SUBROUTINE ZDIV(AR, AI, BR, BI, CR, CI)

c*********************************************************************72
c
cc ZDIV carries out double precision complex division.
c
C***BEGIN PROLOGUE  ZDIV
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
C
C***ROUTINES CALLED  ZABS
C***END PROLOGUE  ZDIV
         EXTERNAL ZABS
         DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
         DOUBLE PRECISION ZABS
         BM = 1.0D0/ZABS(BR,BI)
         CC = BR*BM
         CD = BI*BM
         CA = (AR*CC+AI*CD)*BM
         CB = (AI*CC-AR*CD)*BM
         CR = CA
         CI = CB
         RETURN
      END
      SUBROUTINE ZEXP(AR, AI, BR, BI)

c*********************************************************************72
c
cc ZEXP carries out double precision complex exponentiation.
c
C***BEGIN PROLOGUE  ZEXP
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZEXP
         DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
         ZM = DEXP(AR)
         CA = ZM*DCOS(AI)
         CB = ZM*DSIN(AI)
         BR = CA
         BI = CB
         RETURN
      END
      SUBROUTINE ZKSCL(ZRR,ZRI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZKSCL
C***REFER TO  ZBESK
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***ROUTINES CALLED  ZUCHK,ZABS,ZLOG
C***END PROLOGUE  ZKSCL
C     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
         EXTERNAL ZABS
         DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
     *    CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
     *    S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS,
     *    ZDR, ZDI, CELMR, ELM, HELIM, ALAS
         INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
         DIMENSION YR(N), YI(N), CYR(2), CYI(2)
         DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /
C
         NZ = 0
         IC = 0
         NN = MIN0(2,N)
         DO 10 I=1,NN
            S1R = YR(I)
            S1I = YI(I)
            CYR(I) = S1R
            CYI(I) = S1I
            AS = ZABS(S1R,S1I)
            ACS = -ZRR + DLOG(AS)
            NZ = NZ + 1
            YR(I) = ZEROR
            YI(I) = ZEROI
            IF (ACS.LT.(-ELIM)) GO TO 10
            CALL ZLOG(S1R, S1I, CSR, CSI, IDUM)
            CSR = CSR - ZRR
            CSI = CSI - ZRI
            STR = DEXP(CSR)/TOL
            CSR = STR*DCOS(CSI)
            CSI = STR*DSIN(CSI)
            CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
            IF (NW.NE.0) GO TO 10
            YR(I) = CSR
            YI(I) = CSI
            IC = I
            NZ = NZ - 1
   10    CONTINUE
         IF (N.EQ.1) RETURN
         IF (IC.GT.1) GO TO 20
         YR(1) = ZEROR
         YI(1) = ZEROI
         NZ = 2
   20    CONTINUE
         IF (N.EQ.2) RETURN
         IF (NZ.EQ.0) RETURN
         FN = FNU + 1.0D0
         CKR = FN*RZR
         CKI = FN*RZI
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         HELIM = 0.5D0*ELIM
         ELM = DEXP(-ELIM)
         CELMR = ELM
         ZDR = ZRR
         ZDI = ZRI
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
         DO 30 I=3,N
            KK = I
            CSR = S2R
            CSI = S2I
            S2R = CKR*CSR - CKI*CSI + S1R
            S2I = CKI*CSR + CKR*CSI + S1I
            S1R = CSR
            S1I = CSI
            CKR = CKR + RZR
            CKI = CKI + RZI
            AS = ZABS(S2R,S2I)
            ALAS = DLOG(AS)
            ACS = -ZDR + ALAS
            NZ = NZ + 1
            YR(I) = ZEROR
            YI(I) = ZEROI
            IF (ACS.LT.(-ELIM)) GO TO 25
            CALL ZLOG(S2R, S2I, CSR, CSI, IDUM)
            CSR = CSR - ZDR
            CSI = CSI - ZDI
            STR = DEXP(CSR)/TOL
            CSR = STR*DCOS(CSI)
            CSI = STR*DSIN(CSI)
            CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
            IF (NW.NE.0) GO TO 25
            YR(I) = CSR
            YI(I) = CSI
            NZ = NZ - 1
            IF (IC.EQ.KK-1) GO TO 40
            IC = KK
            GO TO 30
   25       CONTINUE
            IF(ALAS.LT.HELIM) GO TO 30
            ZDR = ZDR - ELIM
            S1R = S1R*CELMR
            S1I = S1I*CELMR
            S2R = S2R*CELMR
            S2I = S2I*CELMR
   30    CONTINUE
         NZ = N
         IF(IC.EQ.N) NZ=N-1
         GO TO 45
   40    CONTINUE
         NZ = KK - 2
   45    CONTINUE
         DO 50 I=1,NZ
            YR(I) = ZEROR
            YI(I) = ZEROI
   50    CONTINUE
         RETURN
      END
      SUBROUTINE ZLOG(AR, AI, BR, BI, IERR)

c*********************************************************************72
c
cc ZLOG carries out double precision complex logarithms
c
C***BEGIN PROLOGUE  ZLOG
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
C     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
C***ROUTINES CALLED  ZABS
C***END PROLOGUE  ZLOG
         EXTERNAL ZABS
         DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
         DOUBLE PRECISION ZABS
         INTEGER IERR
         DATA DPI , DHPI  / 3.141592653589793238462643383D+0,
     1                      1.570796326794896619231321696D+0/
C
         IERR=0
         IF (AR.EQ.0.0D+0) GO TO 10
         IF (AI.EQ.0.0D+0) GO TO 20
         DTHETA = DATAN(AI/AR)
         IF (DTHETA.LE.0.0D+0) GO TO 40
         IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
         GO TO 50
   10    IF (AI.EQ.0.0D+0) GO TO 60
         BI = DHPI
         BR = DLOG(DABS(AI))
         IF (AI.LT.0.0D+0) BI = -BI
         RETURN
   20    IF (AR.GT.0.0D+0) GO TO 30
         BR = DLOG(DABS(AR))
         BI = DPI
         RETURN
   30    BR = DLOG(AR)
         BI = 0.0D+0
         RETURN
   40    IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50    ZM = ZABS(AR,AI)
         BR = DLOG(ZM)
         BI = DTHETA
         RETURN
   60    CONTINUE
         IERR=1
         RETURN
      END
      SUBROUTINE ZMLRI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZMLRI
C***REFER TO  ZBESI,ZBESK
C
C     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***ROUTINES CALLED  DGAMLN,D1MACH,ZABS,ZEXP,ZLOG,ZMLT
C***END PROLOGUE  ZMLRI
C     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
         EXTERNAL ZABS
         DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,
     *    CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I,
     *    P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,
     *    SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, DGAMLN,
     *    D1MACH, ZABS
         INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
         DIMENSION YR(N), YI(N)
         DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
         SCLE = D1MACH(1)/TOL
         NZ=0
         AZ = ZABS(ZR,ZI)
         IAZ = INT(SNGL(AZ))
         IFNU = INT(SNGL(FNU))
         INU = IFNU + N - 1
         AT = DBLE(FLOAT(IAZ)) + 1.0D0
         RAZ = 1.0D0/AZ
         STR = ZR*RAZ
         STI = -ZI*RAZ
         CKR = STR*AT*RAZ
         CKI = STI*AT*RAZ
         RZR = (STR+STR)*RAZ
         RZI = (STI+STI)*RAZ
         P1R = ZEROR
         P1I = ZEROI
         P2R = CONER
         P2I = CONEI
         ACK = (AT+1.0D0)*RAZ
         RHO = ACK + DSQRT(ACK*ACK-1.0D0)
         RHO2 = RHO*RHO
         TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
         TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
         AK = AT
         DO 10 I=1,80
            PTR = P2R
            PTI = P2I
            P2R = P1R - (CKR*PTR-CKI*PTI)
            P2I = P1I - (CKI*PTR+CKR*PTI)
            P1R = PTR
            P1I = PTI
            CKR = CKR + RZR
            CKI = CKI + RZI
            AP = ZABS(P2R,P2I)
            IF (AP.GT.TST*AK*AK) GO TO 20
            AK = AK + 1.0D0
   10    CONTINUE
         GO TO 110
   20    CONTINUE
         I = I + 1
         K = 0
         IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
         P1R = ZEROR
         P1I = ZEROI
         P2R = CONER
         P2I = CONEI
         AT = DBLE(FLOAT(INU)) + 1.0D0
         STR = ZR*RAZ
         STI = -ZI*RAZ
         CKR = STR*AT*RAZ
         CKI = STI*AT*RAZ
         ACK = AT*RAZ
         TST = DSQRT(ACK/TOL)
         ITIME = 1
         DO 30 K=1,80
            PTR = P2R
            PTI = P2I
            P2R = P1R - (CKR*PTR-CKI*PTI)
            P2I = P1I - (CKR*PTI+CKI*PTR)
            P1R = PTR
            P1I = PTI
            CKR = CKR + RZR
            CKI = CKI + RZI
            AP = ZABS(P2R,P2I)
            IF (AP.LT.TST) GO TO 30
            IF (ITIME.EQ.2) GO TO 40
            ACK = ZABS(CKR,CKI)
            FLAM = ACK + DSQRT(ACK*ACK-1.0D0)
            FKAP = AP/ZABS(P1R,P1I)
            RHO = DMIN1(FLAM,FKAP)
            TST = TST*DSQRT(RHO/(RHO*RHO-1.0D0))
            ITIME = 2
   30    CONTINUE
         GO TO 110
   40    CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
         K = K + 1
         KK = MAX0(I+IAZ,K+INU)
         FKK = DBLE(FLOAT(KK))
         P1R = ZEROR
         P1I = ZEROI
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
         P2R = SCLE
         P2I = ZEROI
         FNF = FNU - DBLE(FLOAT(IFNU))
         TFNF = FNF + FNF
         BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) -
     *    DGAMLN(TFNF+1.0D0,IDUM)
         BK = DEXP(BK)
         SUMR = ZEROR
         SUMI = ZEROI
         KM = KK - INU
         DO 50 I=1,KM
            PTR = P2R
            PTI = P2I
            P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
            P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
            P1R = PTR
            P1I = PTI
            AK = 1.0D0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUMR = SUMR + (ACK+BK)*P1R
            SUMI = SUMI + (ACK+BK)*P1I
            BK = ACK
            FKK = FKK - 1.0D0
   50    CONTINUE
         YR(N) = P2R
         YI(N) = P2I
         IF (N.EQ.1) GO TO 70
         DO 60 I=2,N
            PTR = P2R
            PTI = P2I
            P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
            P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
            P1R = PTR
            P1I = PTI
            AK = 1.0D0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUMR = SUMR + (ACK+BK)*P1R
            SUMI = SUMI + (ACK+BK)*P1I
            BK = ACK
            FKK = FKK - 1.0D0
            M = N - I + 1
            YR(M) = P2R
            YI(M) = P2I
   60    CONTINUE
   70    CONTINUE
         IF (IFNU.LE.0) GO TO 90
         DO 80 I=1,IFNU
            PTR = P2R
            PTI = P2I
            P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
            P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
            P1R = PTR
            P1I = PTI
            AK = 1.0D0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUMR = SUMR + (ACK+BK)*P1R
            SUMI = SUMI + (ACK+BK)*P1I
            BK = ACK
            FKK = FKK - 1.0D0
   80    CONTINUE
   90    CONTINUE
         PTR = ZR
         PTI = ZI
         IF (KODE.EQ.2) PTR = ZEROR
         CALL ZLOG(RZR, RZI, STR, STI, IDUM)
         P1R = -FNF*STR + PTR
         P1I = -FNF*STI + PTI
         AP = DGAMLN(1.0D0+FNF,IDUM)
         PTR = P1R - AP
         PTI = P1I
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
         P2R = P2R + SUMR
         P2I = P2I + SUMI
         AP = ZABS(P2R,P2I)
         P1R = 1.0D0/AP
         CALL ZEXP(PTR, PTI, STR, STI)
         CKR = STR*P1R
         CKI = STI*P1R
         PTR = P2R*P1R
         PTI = -P2I*P1R
         CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
         DO 100 I=1,N
            STR = YR(I)*CNORMR - YI(I)*CNORMI
            YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
            YR(I) = STR
  100    CONTINUE
         RETURN
  110    CONTINUE
         NZ=-2
         RETURN
      END
      SUBROUTINE ZMLT(AR, AI, BR, BI, CR, CI)

c*********************************************************************72
c
cc ZMLT carries out double precision complex multiplication.
c
C***BEGIN PROLOGUE  ZMLT
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZMLT
         DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
         CA = AR*BR - AI*BI
         CB = AR*BI + AI*BR
         CR = CA
         CI = CB
         RETURN
      END
      SUBROUTINE ZRATI(ZR, ZI, FNU, N, CYR, CYI, TOL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZRATI
C***REFER TO  ZBESI,ZBESK,ZBESH
C
C     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***ROUTINES CALLED  ZABS,ZDIV
C***END PROLOGUE  ZRATI
C     COMPLEX Z,CY(1),CONE,CZERO,P1,P2,T1,RZ,PT,CDFNU
         EXTERNAL ZABS
         DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     *    CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     *    FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     *    RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZABS
         INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
         DIMENSION CYR(N), CYI(N)
         DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1    0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.41421356237309505D0 /
         AZ = ZABS(ZR,ZI)
         INU = INT(SNGL(FNU))
         IDNU = INU + N - 1
         MAGZ = INT(SNGL(AZ))
         AMAGZ = DBLE(FLOAT(MAGZ+1))
         FDNU = DBLE(FLOAT(IDNU))
         FNUP = DMAX1(AMAGZ,FDNU)
         ID = IDNU - MAGZ - 1
         ITIME = 1
         K = 1
         PTR = 1.0D0/AZ
         RZR = PTR*(ZR+ZR)*PTR
         RZI = -PTR*(ZI+ZI)*PTR
         T1R = RZR*FNUP
         T1I = RZI*FNUP
         P2R = -T1R
         P2I = -T1I
         P1R = CONER
         P1I = CONEI
         T1R = T1R + RZR
         T1I = T1I + RZI
         IF (ID.GT.0) ID = 0
         AP2 = ZABS(P2R,P2I)
         AP1 = ZABS(P1R,P1I)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
         ARG = (AP2+AP2)/(AP1*TOL)
         TEST1 = DSQRT(ARG)
         TEST = TEST1
         RAP1 = 1.0D0/AP1
         P1R = P1R*RAP1
         P1I = P1I*RAP1
         P2R = P2R*RAP1
         P2I = P2I*RAP1
         AP2 = AP2*RAP1
   10    CONTINUE
         K = K + 1
         AP1 = AP2
         PTR = P2R
         PTI = P2I
         P2R = P1R - (T1R*PTR-T1I*PTI)
         P2I = P1I - (T1R*PTI+T1I*PTR)
         P1R = PTR
         P1I = PTI
         T1R = T1R + RZR
         T1I = T1I + RZI
         AP2 = ZABS(P2R,P2I)
         IF (AP1.LE.TEST) GO TO 10
         IF (ITIME.EQ.2) GO TO 20
         AK = ZABS(T1R,T1I)*0.5D0
         FLAM = AK + DSQRT(AK*AK-1.0D0)
         RHO = DMIN1(AP2/AP1,FLAM)
         TEST = TEST1*DSQRT(RHO/(RHO*RHO-1.0D0))
         ITIME = 2
         GO TO 10
   20    CONTINUE
         KK = K + 1 - ID
         AK = DBLE(FLOAT(KK))
         T1R = AK
         T1I = CZEROI
         DFNU = FNU + DBLE(FLOAT(N-1))
         P1R = 1.0D0/AP2
         P1I = CZEROI
         P2R = CZEROR
         P2I = CZEROI
         DO 30 I=1,KK
            PTR = P1R
            PTI = P1I
            RAP1 = DFNU + T1R
            TTR = RZR*RAP1
            TTI = RZI*RAP1
            P1R = (PTR*TTR-PTI*TTI) + P2R
            P1I = (PTR*TTI+PTI*TTR) + P2I
            P2R = PTR
            P2I = PTI
            T1R = T1R - CONER
   30    CONTINUE
         IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
         P1R = TOL
         P1I = TOL
   40    CONTINUE
         CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
         IF (N.EQ.1) RETURN
         K = N - 1
         AK = DBLE(FLOAT(K))
         T1R = AK
         T1I = CZEROI
         CDFNUR = FNU*RZR
         CDFNUI = FNU*RZI
         DO 60 I=2,N
            PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
            PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
            AK = ZABS(PTR,PTI)
            IF (AK.NE.CZEROR) GO TO 50
            PTR = TOL
            PTI = TOL
            AK = TOL*RT2
   50       CONTINUE
            RAK = CONER/AK
            CYR(K) = RAK*PTR*RAK
            CYI(K) = -RAK*PTI*RAK
            T1R = T1R - CONER
            K = K - 1
   60    CONTINUE
         RETURN
      END
      SUBROUTINE ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     * IUF)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZS1S2
C***REFER TO  ZBESK,ZAIRY
C
C     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***ROUTINES CALLED  ZABS,ZEXP,ZLOG
C***END PROLOGUE  ZS1S2
C     COMPLEX CZERO,C1,S1,S1D,S2,ZR
         EXTERNAL ZABS
         DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
     *    S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
         INTEGER IUF, IDUM, NZ
         DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
         NZ = 0
         AS1 = ZABS(S1R,S1I)
         AS2 = ZABS(S2R,S2I)
         IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
         IF (AS1.EQ.0.0D0) GO TO 10
         ALN = -ZRR - ZRR + DLOG(AS1)
         S1DR = S1R
         S1DI = S1I
         S1R = ZEROR
         S1I = ZEROI
         AS1 = ZEROR
         IF (ALN.LT.(-ALIM)) GO TO 10
         CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
         C1R = C1R - ZRR - ZRR
         C1I = C1I - ZRI - ZRI
         CALL ZEXP(C1R, C1I, S1R, S1I)
         AS1 = ZABS(S1R,S1I)
         IUF = IUF + 1
   10    CONTINUE
         AA = DMAX1(AS1,AS2)
         IF (AA.GT.ASCLE) RETURN
         S1R = ZEROR
         S1I = ZEROI
         S2R = ZEROR
         S2I = ZEROI
         NZ = 1
         IUF = 0
         RETURN
      END
      SUBROUTINE ZSERI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZSERI
C***REFER TO  ZBESI,ZBESK
C
C     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***ROUTINES CALLED  DGAMLN,D1MACH,ZUCHK,ZABS,ZDIV,ZLOG,ZMLT
C***END PROLOGUE  ZSERI
C     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
         EXTERNAL ZABS
         DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL,
     *    AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU,
     *    ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,
     *    STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,
     *    ZR, DGAMLN, D1MACH, ZABS
         INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
         DIMENSION YR(N), YI(N), WR(2), WI(2)
         DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
         NZ = 0
         AZ = ZABS(ZR,ZI)
         IF (AZ.EQ.0.0D0) GO TO 160
         ARM = 1.0D+3*D1MACH(1)
         RTR1 = DSQRT(ARM)
         CRSCR = 1.0D0
         IFLAG = 0
         IF (AZ.LT.ARM) GO TO 150
         HZR = 0.5D0*ZR
         HZI = 0.5D0*ZI
         CZR = ZEROR
         CZI = ZEROI
         IF (AZ.LE.RTR1) GO TO 10
         CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10    CONTINUE
         ACZ = ZABS(CZR,CZI)
         NN = N
         CALL ZLOG(HZR, HZI, CKR, CKI, IDUM)
   20    CONTINUE
         DFNU = FNU + DBLE(FLOAT(NN-1))
         FNUP = DFNU + 1.0D0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
         AK1R = CKR*DFNU
         AK1I = CKI*DFNU
         AK = DGAMLN(FNUP,IDUM)
         AK1R = AK1R - AK
         IF (KODE.EQ.2) AK1R = AK1R - ZR
         IF (AK1R.GT.(-ELIM)) GO TO 40
   30    CONTINUE
         NZ = NZ + 1
         YR(NN) = ZEROR
         YI(NN) = ZEROI
         IF (ACZ.GT.DFNU) GO TO 190
         NN = NN - 1
         IF (NN.EQ.0) RETURN
         GO TO 20
   40    CONTINUE
         IF (AK1R.GT.(-ALIM)) GO TO 50
         IFLAG = 1
         SS = 1.0D0/TOL
         CRSCR = TOL
         ASCLE = ARM*SS
   50    CONTINUE
         AA = DEXP(AK1R)
         IF (IFLAG.EQ.1) AA = AA*SS
         COEFR = AA*DCOS(AK1I)
         COEFI = AA*DSIN(AK1I)
         ATOL = TOL*ACZ/FNUP
         IL = MIN0(2,NN)
         DO 90 I=1,IL
            DFNU = FNU + DBLE(FLOAT(NN-I))
            FNUP = DFNU + 1.0D0
            S1R = CONER
            S1I = CONEI
            IF (ACZ.LT.TOL*FNUP) GO TO 70
            AK1R = CONER
            AK1I = CONEI
            AK = FNUP + 2.0D0
            S = FNUP
            AA = 2.0D0
   60       CONTINUE
            RS = 1.0D0/S
            STR = AK1R*CZR - AK1I*CZI
            STI = AK1R*CZI + AK1I*CZR
            AK1R = STR*RS
            AK1I = STI*RS
            S1R = S1R + AK1R
            S1I = S1I + AK1I
            S = S + AK
            AK = AK + 2.0D0
            AA = AA*ACZ*RS
            IF (AA.GT.ATOL) GO TO 60
   70       CONTINUE
            S2R = S1R*COEFR - S1I*COEFI
            S2I = S1R*COEFI + S1I*COEFR
            WR(I) = S2R
            WI(I) = S2I
            IF (IFLAG.EQ.0) GO TO 80
            CALL ZUCHK(S2R, S2I, NW, ASCLE, TOL)
            IF (NW.NE.0) GO TO 30
   80       CONTINUE
            M = NN - I + 1
            YR(M) = S2R*CRSCR
            YI(M) = S2I*CRSCR
            IF (I.EQ.IL) GO TO 90
            CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
            COEFR = STR*DFNU
            COEFI = STI*DFNU
   90    CONTINUE
         IF (NN.LE.2) RETURN
         K = NN - 2
         AK = DBLE(FLOAT(K))
         RAZ = 1.0D0/AZ
         STR = ZR*RAZ
         STI = -ZI*RAZ
         RZR = (STR+STR)*RAZ
         RZI = (STI+STI)*RAZ
         IF (IFLAG.EQ.1) GO TO 120
         IB = 3
  100    CONTINUE
         DO 110 I=IB,NN
            YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
            YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
            AK = AK - 1.0D0
            K = K - 1
  110    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  120    CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
C-----------------------------------------------------------------------
         S1R = WR(1)
         S1I = WI(1)
         S2R = WR(2)
         S2I = WI(2)
         DO 130 L=3,NN
            CKR = S2R
            CKI = S2I
            S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
            S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
            S1R = CKR
            S1I = CKI
            CKR = S2R*CRSCR
            CKI = S2I*CRSCR
            YR(K) = CKR
            YI(K) = CKI
            AK = AK - 1.0D0
            K = K - 1
            IF (ZABS(CKR,CKI).GT.ASCLE) GO TO 140
  130    CONTINUE
         RETURN
  140    CONTINUE
         IB = L + 1
         IF (IB.GT.NN) RETURN
         GO TO 100
  150    CONTINUE
         NZ = N
         IF (FNU.EQ.0.0D0) NZ = NZ - 1
  160    CONTINUE
         YR(1) = ZEROR
         YI(1) = ZEROI
         IF (FNU.NE.0.0D0) GO TO 170
         YR(1) = CONER
         YI(1) = CONEI
  170    CONTINUE
         IF (N.EQ.1) RETURN
         DO 180 I=2,N
            YR(I) = ZEROR
            YI(I) = ZEROI
  180    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
C-----------------------------------------------------------------------
  190    CONTINUE
         NZ = -NZ
         RETURN
      END
      SUBROUTINE ZSHCH(ZR, ZI, CSHR, CSHI, CCHR, CCHI)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZSHCH
C***REFER TO  ZBESK,ZBESH
C
C     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZSHCH
C
         DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR,
     *    DCOSH, DSINH
         SH = DSINH(ZR)
         CH = DCOSH(ZR)
         SN = DSIN(ZI)
         CN = DCOS(ZI)
         CSHR = SH*CN
         CSHI = CH*SN
         CCHR = CH*CN
         CCHI = SH*SN
         RETURN
      END
      SUBROUTINE ZSQRT(AR, AI, BR, BI)

c*********************************************************************72
c
cc ZSQRT carries out double precision complex square roots.
c
C***BEGIN PROLOGUE  ZSQRT
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
C
C***ROUTINES CALLED  ZABS
C***END PROLOGUE  ZSQRT
         EXTERNAL ZABS
         DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
         DOUBLE PRECISION ZABS
         DATA DRT , DPI / 7.071067811865475244008443621D-1,
     1                    3.141592653589793238462643383D+0/
         ZM = ZABS(AR,AI)
         ZM = DSQRT(ZM)
         IF (AR.EQ.0.0D+0) GO TO 10
         IF (AI.EQ.0.0D+0) GO TO 20
         DTHETA = DATAN(AI/AR)
         IF (DTHETA.LE.0.0D+0) GO TO 40
         IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
         GO TO 50
   10    IF (AI.GT.0.0D+0) GO TO 60
         IF (AI.LT.0.0D+0) GO TO 70
         BR = 0.0D+0
         BI = 0.0D+0
         RETURN
   20    IF (AR.GT.0.0D+0) GO TO 30
         BR = 0.0D+0
         BI = DSQRT(DABS(AR))
         RETURN
   30    BR = DSQRT(AR)
         BI = 0.0D+0
         RETURN
   40    IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50    DTHETA = DTHETA*0.5D+0
         BR = ZM*DCOS(DTHETA)
         BI = ZM*DSIN(DTHETA)
         RETURN
   60    BR = ZM*DRT
         BI = ZM*DRT
         RETURN
   70    BR = ZM*DRT
         BI = -ZM*DRT
         RETURN
      END
      SUBROUTINE ZUCHK(YR, YI, NZ, ASCLE, TOL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUCHK
C***REFER TO ZSERI,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ZUCHK
C
C     COMPLEX Y
         DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
         INTEGER NZ
         NZ = 0
         WR = DABS(YR)
         WI = DABS(YI)
         ST = DMIN1(WR,WI)
         IF (ST.GT.ASCLE) RETURN
         SS = DMAX1(WR,WI)
         ST = ST/TOL
         IF (SS.LT.ST) NZ = 1
         RETURN
      END
      SUBROUTINE ZUNHJ(ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUNHJ
C***REFER TO  ZBESI,ZBESK
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***ROUTINES CALLED  ZABS,ZDIV,ZLOG,ZSQRT,D1MACH
C***END PROLOGUE  ZUNHJ
C     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
C    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
C    *ZETA2,ZTH
         EXTERNAL ZABS
         DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,
     *    ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,
     *    CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, HPI,
     *    PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,
     *    RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
     *    SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, THPI, TOL, TZAI,
     *    TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
     *    ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,
     *    ZETA2R, ZI, ZR, ZTHI, ZTHR, ZABS, AC, D1MACH
         INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     *    LRP1, L1, L2, M, IDUM
         DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     *    AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),
     *    DRR(14), DRI(14)
         DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1        AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2        1.00000000000000000D+00,     1.04166666666666667D-01,
     3        8.35503472222222222D-02,     1.28226574556327160D-01,
     4        2.91849026464140464D-01,     8.81627267443757652D-01,
     5        3.32140828186276754D+00,     1.49957629868625547D+01,
     6        7.89230130115865181D+01,     4.74451538868264323D+02,
     7        3.20749009089066193D+03,     2.40865496408740049D+04,
     8        1.98923119169509794D+05,     1.79190200777534383D+06/
         DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1        BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2        1.00000000000000000D+00,    -1.45833333333333333D-01,
     3       -9.87413194444444444D-02,    -1.43312053915895062D-01,
     4       -3.17227202678413548D-01,    -9.42429147957120249D-01,
     5       -3.51120304082635426D+00,    -1.57272636203680451D+01,
     6       -8.22814390971859444D+01,    -4.92355370523670524D+02,
     7       -3.31621856854797251D+03,    -2.48276742452085896D+04,
     8       -2.04526587315129788D+05,    -1.83844491706820990D+06/
         DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1        C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2        C(19), C(20), C(21), C(22), C(23), C(24)/
     3        1.00000000000000000D+00,    -2.08333333333333333D-01,
     4        1.25000000000000000D-01,     3.34201388888888889D-01,
     5       -4.01041666666666667D-01,     7.03125000000000000D-02,
     6       -1.02581259645061728D+00,     1.84646267361111111D+00,
     7       -8.91210937500000000D-01,     7.32421875000000000D-02,
     8        4.66958442342624743D+00,    -1.12070026162229938D+01,
     9        8.78912353515625000D+00,    -2.36408691406250000D+00,
     A        1.12152099609375000D-01,    -2.82120725582002449D+01,
     B        8.46362176746007346D+01,    -9.18182415432400174D+01,
     C        4.25349987453884549D+01,    -7.36879435947963170D+00,
     D        2.27108001708984375D-01,     2.12570130039217123D+02,
     E       -7.65252468141181642D+02,     1.05999045252799988D+03/
         DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1        C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2        C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3       -6.99579627376132541D+02,     2.18190511744211590D+02,
     4       -2.64914304869515555D+01,     5.72501420974731445D-01,
     5       -1.91945766231840700D+03,     8.06172218173730938D+03,
     6       -1.35865500064341374D+04,     1.16553933368645332D+04,
     7       -5.30564697861340311D+03,     1.20090291321635246D+03,
     8       -1.08090919788394656D+02,     1.72772750258445740D+00,
     9        2.02042913309661486D+04,    -9.69805983886375135D+04,
     A        1.92547001232531532D+05,    -2.03400177280415534D+05,
     B        1.22200464983017460D+05,    -4.11926549688975513D+04,
     C        7.10951430248936372D+03,    -4.93915304773088012D+02,
     D        6.07404200127348304D+00,    -2.42919187900551333D+05,
     E        1.31176361466297720D+06,    -2.99801591853810675D+06/
         DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1        C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2        C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3        3.76327129765640400D+06,    -2.81356322658653411D+06,
     4        1.26836527332162478D+06,    -3.31645172484563578D+05,
     5        4.52187689813627263D+04,    -2.49983048181120962D+03,
     6        2.43805296995560639D+01,     3.28446985307203782D+06,
     7       -1.97068191184322269D+07,     5.09526024926646422D+07,
     8       -7.41051482115326577D+07,     6.63445122747290267D+07,
     9       -3.75671766607633513D+07,     1.32887671664218183D+07,
     A       -2.78561812808645469D+06,     3.08186404612662398D+05,
     B       -1.38860897537170405D+04,     1.10017140269246738D+02,
     C       -4.93292536645099620D+07,     3.25573074185765749D+08,
     D       -9.39462359681578403D+08,     1.55359689957058006D+09,
     E       -1.62108055210833708D+09,     1.10684281682301447D+09/
         DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1        C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2        C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3       -4.95889784275030309D+08,     1.42062907797533095D+08,
     4       -2.44740627257387285D+07,     2.24376817792244943D+06,
     5       -8.40054336030240853D+04,     5.51335896122020586D+02,
     6        8.14789096118312115D+08,    -5.86648149205184723D+09,
     7        1.86882075092958249D+10,    -3.46320433881587779D+10,
     8        4.12801855797539740D+10,    -3.30265997498007231D+10,
     9        1.79542137311556001D+10,    -6.56329379261928433D+09,
     A        1.55927986487925751D+09,    -2.25105661889415278D+08,
     B        1.73951075539781645D+07,    -5.49842327572288687D+05,
     C        3.03809051092238427D+03,    -1.46792612476956167D+10,
     D        1.14498237732025810D+11,    -3.99096175224466498D+11,
     E        8.19218669548577329D+11,    -1.09837515608122331D+12/
         DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1        C(105)/
     2        1.00815810686538209D+12,    -6.45364869245376503D+11,
     3        2.87900649906150589D+11,    -8.78670721780232657D+10,
     4        1.76347306068349694D+10,    -2.16716498322379509D+09,
     5        1.43157876718888981D+08,    -3.87183344257261262D+06,
     6        1.82577554742931747D+04/
         DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1        ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2        ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3        ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4       -4.44444444444444444D-03,    -9.22077922077922078D-04,
     5       -8.84892884892884893D-05,     1.65927687832449737D-04,
     6        2.46691372741792910D-04,     2.65995589346254780D-04,
     7        2.61824297061500945D-04,     2.48730437344655609D-04,
     8        2.32721040083232098D-04,     2.16362485712365082D-04,
     9        2.00738858762752355D-04,     1.86267636637545172D-04,
     A        1.73060775917876493D-04,     1.61091705929015752D-04,
     B        1.50274774160908134D-04,     1.40503497391269794D-04,
     C        1.31668816545922806D-04,     1.23667445598253261D-04,
     D        1.16405271474737902D-04,     1.09798298372713369D-04,
     E        1.03772410422992823D-04,     9.82626078369363448D-05/
         DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1        ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2        ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3        ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4        9.32120517249503256D-05,     8.85710852478711718D-05,
     5        8.42963105715700223D-05,     8.03497548407791151D-05,
     6        7.66981345359207388D-05,     7.33122157481777809D-05,
     7        7.01662625163141333D-05,     6.72375633790160292D-05,
     8        6.93735541354588974D-04,     2.32241745182921654D-04,
     9       -1.41986273556691197D-05,    -1.16444931672048640D-04,
     A       -1.50803558053048762D-04,    -1.55121924918096223D-04,
     B       -1.46809756646465549D-04,    -1.33815503867491367D-04,
     C       -1.19744975684254051D-04,    -1.06184319207974020D-04,
     D       -9.37699549891194492D-05,    -8.26923045588193274D-05,
     E       -7.29374348155221211D-05,    -6.44042357721016283D-05/
         DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1        ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2        ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3        ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4       -5.69611566009369048D-05,    -5.04731044303561628D-05,
     5       -4.48134868008882786D-05,    -3.98688727717598864D-05,
     6       -3.55400532972042498D-05,    -3.17414256609022480D-05,
     7       -2.83996793904174811D-05,    -2.54522720634870566D-05,
     8       -2.28459297164724555D-05,    -2.05352753106480604D-05,
     9       -1.84816217627666085D-05,    -1.66519330021393806D-05,
     A       -1.50179412980119482D-05,    -1.35554031379040526D-05,
     B       -1.22434746473858131D-05,    -1.10641884811308169D-05,
     C       -3.54211971457743841D-04,    -1.56161263945159416D-04,
     D        3.04465503594936410D-05,     1.30198655773242693D-04,
     E        1.67471106699712269D-04,     1.70222587683592569D-04/
         DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1        ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2        ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3        ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4        1.56501427608594704D-04,     1.36339170977445120D-04,
     5        1.14886692029825128D-04,     9.45869093034688111D-05,
     6        7.64498419250898258D-05,     6.07570334965197354D-05,
     7        4.74394299290508799D-05,     3.62757512005344297D-05,
     8        2.69939714979224901D-05,     1.93210938247939253D-05,
     9        1.30056674793963203D-05,     7.82620866744496661D-06,
     A        3.59257485819351583D-06,     1.44040049814251817D-07,
     B       -2.65396769697939116D-06,    -4.91346867098485910D-06,
     C       -6.72739296091248287D-06,    -8.17269379678657923D-06,
     D       -9.31304715093561232D-06,    -1.02011418798016441D-05,
     E       -1.08805962510592880D-05,    -1.13875481509603555D-05/
         DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1        ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2        ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3        ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4       -1.17519675674556414D-05,    -1.19987364870944141D-05,
     5        3.78194199201772914D-04,     2.02471952761816167D-04,
     6       -6.37938506318862408D-05,    -2.38598230603005903D-04,
     7       -3.10916256027361568D-04,    -3.13680115247576316D-04,
     8       -2.78950273791323387D-04,    -2.28564082619141374D-04,
     9       -1.75245280340846749D-04,    -1.25544063060690348D-04,
     A       -8.22982872820208365D-05,    -4.62860730588116458D-05,
     B       -1.72334302366962267D-05,     5.60690482304602267D-06,
     C        2.31395443148286800D-05,     3.62642745856793957D-05,
     D        4.58006124490188752D-05,     5.24595294959114050D-05,
     E        5.68396208545815266D-05,     5.94349820393104052D-05/
         DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1        ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2        ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3        ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4        6.06478527578421742D-05,     6.08023907788436497D-05,
     5        6.01577894539460388D-05,     5.89199657344698500D-05,
     6        5.72515823777593053D-05,     5.52804375585852577D-05,
     7        5.31063773802880170D-05,     5.08069302012325706D-05,
     8        4.84418647620094842D-05,     4.60568581607475370D-05,
     9       -6.91141397288294174D-04,    -4.29976633058871912D-04,
     A        1.83067735980039018D-04,     6.60088147542014144D-04,
     B        8.75964969951185931D-04,     8.77335235958235514D-04,
     C        7.49369585378990637D-04,     5.63832329756980918D-04,
     D        3.68059319971443156D-04,     1.88464535514455599D-04/
         DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1        ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2        ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3        ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4        3.70663057664904149D-05,    -8.28520220232137023D-05,
     5       -1.72751952869172998D-04,    -2.36314873605872983D-04,
     6       -2.77966150694906658D-04,    -3.02079514155456919D-04,
     7       -3.12594712643820127D-04,    -3.12872558758067163D-04,
     8       -3.05678038466324377D-04,    -2.93226470614557331D-04,
     9       -2.77255655582934777D-04,    -2.59103928467031709D-04,
     A       -2.39784014396480342D-04,    -2.20048260045422848D-04,
     B       -2.00443911094971498D-04,    -1.81358692210970687D-04,
     C       -1.63057674478657464D-04,    -1.45712672175205844D-04,
     D       -1.29425421983924587D-04,    -1.14245691942445952D-04/
         DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1        ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2        ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3        ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4        1.92821964248775885D-03,     1.35592576302022234D-03,
     5       -7.17858090421302995D-04,    -2.58084802575270346D-03,
     6       -3.49271130826168475D-03,    -3.46986299340960628D-03,
     7       -2.82285233351310182D-03,    -1.88103076404891354D-03,
     8       -8.89531718383947600D-04,     3.87912102631035228D-06,
     9        7.28688540119691412D-04,     1.26566373053457758D-03,
     A        1.62518158372674427D-03,     1.83203153216373172D-03,
     B        1.91588388990527909D-03,     1.90588846755546138D-03,
     C        1.82798982421825727D-03,     1.70389506421121530D-03,
     D        1.55097127171097686D-03,     1.38261421852276159D-03/
         DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1        ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2        1.20881424230064774D-03,     1.03676532638344962D-03,
     3        8.71437918068619115D-04,     7.16080155297701002D-04,
     4        5.72637002558129372D-04,     4.42089819465802277D-04,
     5        3.24724948503090564D-04,     2.20342042730246599D-04,
     6        1.28412898401353882D-04,     4.82005924552095464D-05/
         DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1        BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2        BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3        BETA(19), BETA(20), BETA(21), BETA(22)/
     4        1.79988721413553309D-02,     5.59964911064388073D-03,
     5        2.88501402231132779D-03,     1.80096606761053941D-03,
     6        1.24753110589199202D-03,     9.22878876572938311D-04,
     7        7.14430421727287357D-04,     5.71787281789704872D-04,
     8        4.69431007606481533D-04,     3.93232835462916638D-04,
     9        3.34818889318297664D-04,     2.88952148495751517D-04,
     A        2.52211615549573284D-04,     2.22280580798883327D-04,
     B        1.97541838033062524D-04,     1.76836855019718004D-04,
     C        1.59316899661821081D-04,     1.44347930197333986D-04,
     D        1.31448068119965379D-04,     1.20245444949302884D-04,
     E        1.10449144504599392D-04,     1.01828770740567258D-04/
         DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1        BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2        BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3        BETA(41), BETA(42), BETA(43), BETA(44)/
     4        9.41998224204237509D-05,     8.74130545753834437D-05,
     5        8.13466262162801467D-05,     7.59002269646219339D-05,
     6        7.09906300634153481D-05,     6.65482874842468183D-05,
     7        6.25146958969275078D-05,     5.88403394426251749D-05,
     8       -1.49282953213429172D-03,    -8.78204709546389328D-04,
     9       -5.02916549572034614D-04,    -2.94822138512746025D-04,
     A       -1.75463996970782828D-04,    -1.04008550460816434D-04,
     B       -5.96141953046457895D-05,    -3.12038929076098340D-05,
     C       -1.26089735980230047D-05,    -2.42892608575730389D-07,
     D        8.05996165414273571D-06,     1.36507009262147391D-05,
     E        1.73964125472926261D-05,     1.98672978842133780D-05/
         DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1        BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2        BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3        BETA(63), BETA(64), BETA(65), BETA(66)/
     4        2.14463263790822639D-05,     2.23954659232456514D-05,
     5        2.28967783814712629D-05,     2.30785389811177817D-05,
     6        2.30321976080909144D-05,     2.28236073720348722D-05,
     7        2.25005881105292418D-05,     2.20981015361991429D-05,
     8        2.16418427448103905D-05,     2.11507649256220843D-05,
     9        2.06388749782170737D-05,     2.01165241997081666D-05,
     A        1.95913450141179244D-05,     1.90689367910436740D-05,
     B        1.85533719641636667D-05,     1.80475722259674218D-05,
     C        5.52213076721292790D-04,     4.47932581552384646D-04,
     D        2.79520653992020589D-04,     1.52468156198446602D-04,
     E        6.93271105657043598D-05,     1.76258683069991397D-05/
         DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1        BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2        BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3        BETA(85), BETA(86), BETA(87), BETA(88)/
     4       -1.35744996343269136D-05,    -3.17972413350427135D-05,
     5       -4.18861861696693365D-05,    -4.69004889379141029D-05,
     6       -4.87665447413787352D-05,    -4.87010031186735069D-05,
     7       -4.74755620890086638D-05,    -4.55813058138628452D-05,
     8       -4.33309644511266036D-05,    -4.09230193157750364D-05,
     9       -3.84822638603221274D-05,    -3.60857167535410501D-05,
     A       -3.37793306123367417D-05,    -3.15888560772109621D-05,
     B       -2.95269561750807315D-05,    -2.75978914828335759D-05,
     C       -2.58006174666883713D-05,    -2.41308356761280200D-05,
     D       -2.25823509518346033D-05,    -2.11479656768912971D-05,
     E       -1.98200638885294927D-05,    -1.85909870801065077D-05/
         DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1        BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2        BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3        BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4       -1.74532699844210224D-05,    -1.63997823854497997D-05,
     5       -4.74617796559959808D-04,    -4.77864567147321487D-04,
     6       -3.20390228067037603D-04,    -1.61105016119962282D-04,
     7       -4.25778101285435204D-05,     3.44571294294967503D-05,
     8        7.97092684075674924D-05,     1.03138236708272200D-04,
     9        1.12466775262204158D-04,     1.13103642108481389D-04,
     A        1.08651634848774268D-04,     1.01437951597661973D-04,
     B        9.29298396593363896D-05,     8.40293133016089978D-05,
     C        7.52727991349134062D-05,     6.69632521975730872D-05,
     D        5.92564547323194704D-05,     5.22169308826975567D-05,
     E        4.58539485165360646D-05,     4.01445513891486808D-05/
         DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1        BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2        BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3        BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4        3.50481730031328081D-05,     3.05157995034346659D-05,
     5        2.64956119950516039D-05,     2.29363633690998152D-05,
     6        1.97893056664021636D-05,     1.70091984636412623D-05,
     7        1.45547428261524004D-05,     1.23886640995878413D-05,
     8        1.04775876076583236D-05,     8.79179954978479373D-06,
     9        7.36465810572578444D-04,     8.72790805146193976D-04,
     A        6.22614862573135066D-04,     2.85998154194304147D-04,
     B        3.84737672879366102D-06,    -1.87906003636971558D-04,
     C       -2.97603646594554535D-04,    -3.45998126832656348D-04,
     D       -3.53382470916037712D-04,    -3.35715635775048757D-04/
         DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1        BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2        BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3        BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4       -3.04321124789039809D-04,    -2.66722723047612821D-04,
     5       -2.27654214122819527D-04,    -1.89922611854562356D-04,
     6       -1.55058918599093870D-04,    -1.23778240761873630D-04,
     7       -9.62926147717644187D-05,    -7.25178327714425337D-05,
     8       -5.22070028895633801D-05,    -3.50347750511900522D-05,
     9       -2.06489761035551757D-05,    -8.70106096849767054D-06,
     A        1.13698686675100290D-06,     9.16426474122778849D-06,
     B        1.56477785428872620D-05,     2.08223629482466847D-05,
     C        2.48923381004595156D-05,     2.80340509574146325D-05,
     D        3.03987774629861915D-05,     3.21156731406700616D-05/
         DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1        BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2        BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3        BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4       -1.80182191963885708D-03,    -2.43402962938042533D-03,
     5       -1.83422663549856802D-03,    -7.62204596354009765D-04,
     6        2.39079475256927218D-04,     9.49266117176881141D-04,
     7        1.34467449701540359D-03,     1.48457495259449178D-03,
     8        1.44732339830617591D-03,     1.30268261285657186D-03,
     9        1.10351597375642682D-03,     8.86047440419791759D-04,
     A        6.73073208165665473D-04,     4.77603872856582378D-04,
     B        3.05991926358789362D-04,     1.60315694594721630D-04,
     C        4.00749555270613286D-05,    -5.66607461635251611D-05,
     D       -1.32506186772982638D-04,    -1.90296187989614057D-04/
         DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1        BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2        BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3        BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4       -2.32811450376937408D-04,    -2.62628811464668841D-04,
     5       -2.82050469867598672D-04,    -2.93081563192861167D-04,
     6       -2.97435962176316616D-04,    -2.96557334239348078D-04,
     7       -2.91647363312090861D-04,    -2.83696203837734166D-04,
     8       -2.73512317095673346D-04,    -2.61750155806768580D-04,
     9        6.38585891212050914D-03,     9.62374215806377941D-03,
     A        7.61878061207001043D-03,     2.83219055545628054D-03,
     B       -2.09841352012720090D-03,    -5.73826764216626498D-03,
     C       -7.70804244495414620D-03,    -8.21011692264844401D-03,
     D       -7.65824520346905413D-03,    -6.47209729391045177D-03/
         DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1        BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2        BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3        BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4       -4.99132412004966473D-03,    -3.45612289713133280D-03,
     5       -2.01785580014170775D-03,    -7.59430686781961401D-04,
     6        2.84173631523859138D-04,     1.10891667586337403D-03,
     7        1.72901493872728771D-03,     2.16812590802684701D-03,
     8        2.45357710494539735D-03,     2.61281821058334862D-03,
     9        2.67141039656276912D-03,     2.65203073395980430D-03,
     A        2.57411652877287315D-03,     2.45389126236094427D-03,
     B        2.30460058071795494D-03,     2.13684837686712662D-03,
     C        1.95896528478870911D-03,     1.77737008679454412D-03,
     D        1.59690280765839059D-03,     1.42111975664438546D-03/
         DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1        GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2        GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3        GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4        6.29960524947436582D-01,     2.51984209978974633D-01,
     5        1.54790300415655846D-01,     1.10713062416159013D-01,
     6        8.57309395527394825D-02,     6.97161316958684292D-02,
     7        5.86085671893713576D-02,     5.04698873536310685D-02,
     8        4.42600580689154809D-02,     3.93720661543509966D-02,
     9        3.54283195924455368D-02,     3.21818857502098231D-02,
     A        2.94646240791157679D-02,     2.71581677112934479D-02,
     B        2.51768272973861779D-02,     2.34570755306078891D-02,
     C        2.19508390134907203D-02,     2.06210828235646240D-02,
     D        1.94388240897880846D-02,     1.83810633800683158D-02,
     E        1.74293213231963172D-02,     1.65685837786612353D-02/
         DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1        GAMA(29), GAMA(30)/
     2        1.57865285987918445D-02,     1.50729501494095594D-02,
     3        1.44193250839954639D-02,     1.38184805735341786D-02,
     4        1.32643378994276568D-02,     1.27517121970498651D-02,
     5        1.22761545318762767D-02,     1.18338262398482403D-02/
         DATA EX1, EX2, HPI, GPI, THPI /
     1        3.33333333333333333D-01,     6.66666666666666667D-01,
     2        1.57079632679489662D+00,     3.14159265358979324D+00,
     3        4.71238898038468986D+00/
         DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C
         RFNU = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
         TEST = D1MACH(1)*1.0D+3
         AC = FNU*TEST
         IF (DABS(ZR).GT.AC .OR. DABS(ZI).GT.AC) GO TO 15
         ZETA1R = 2.0D0*DABS(DLOG(TEST))+FNU
         ZETA1I = 0.0D0
         ZETA2R = FNU
         ZETA2I = 0.0D0
         PHIR = 1.0D0
         PHII = 0.0D0
         ARGR = 1.0D0
         ARGI = 0.0D0
         RETURN
   15    CONTINUE
         ZBR = ZR*RFNU
         ZBI = ZI*RFNU
         RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
         FN13 = FNU**EX1
         FN23 = FN13*FN13
         RFN13 = 1.0D0/FN13
         W2R = CONER - ZBR*ZBR + ZBI*ZBI
         W2I = CONEI - ZBR*ZBI - ZBR*ZBI
         AW2 = ZABS(W2R,W2I)
         IF (AW2.GT.0.25D0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(W2).LE.0.25D0
C-----------------------------------------------------------------------
         K = 1
         PR(1) = CONER
         PI(1) = CONEI
         SUMAR = GAMA(1)
         SUMAI = ZEROI
         AP(1) = 1.0D0
         IF (AW2.LT.TOL) GO TO 20
         DO 10 K=2,30
            PR(K) = PR(K-1)*W2R - PI(K-1)*W2I
            PI(K) = PR(K-1)*W2I + PI(K-1)*W2R
            SUMAR = SUMAR + PR(K)*GAMA(K)
            SUMAI = SUMAI + PI(K)*GAMA(K)
            AP(K) = AP(K-1)*AW2
            IF (AP(K).LT.TOL) GO TO 20
   10    CONTINUE
         K = 30
   20    CONTINUE
         KMAX = K
         ZETAR = W2R*SUMAR - W2I*SUMAI
         ZETAI = W2R*SUMAI + W2I*SUMAR
         ARGR = ZETAR*FN23
         ARGI = ZETAI*FN23
         CALL ZSQRT(SUMAR, SUMAI, ZAR, ZAI)
         CALL ZSQRT(W2R, W2I, STR, STI)
         ZETA2R = STR*FNU
         ZETA2I = STI*FNU
         STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
         STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
         ZETA1R = STR*ZETA2R - STI*ZETA2I
         ZETA1I = STR*ZETA2I + STI*ZETA2R
         ZAR = ZAR + ZAR
         ZAI = ZAI + ZAI
         CALL ZSQRT(ZAR, ZAI, STR, STI)
         PHIR = STR*RFN13
         PHII = STI*RFN13
         IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
         SUMBR = ZEROR
         SUMBI = ZEROI
         DO 30 K=1,KMAX
            SUMBR = SUMBR + PR(K)*BETA(K)
            SUMBI = SUMBI + PI(K)*BETA(K)
   30    CONTINUE
         ASUMR = ZEROR
         ASUMI = ZEROI
         BSUMR = SUMBR
         BSUMI = SUMBI
         L1 = 0
         L2 = 30
         BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
         ATOL = TOL
         PP = 1.0D0
         IAS = 0
         IBS = 0
         IF (RFNU2.LT.TOL) GO TO 110
         DO 100 IS=2,7
            ATOL = ATOL/RFNU2
            PP = PP*RFNU2
            IF (IAS.EQ.1) GO TO 60
            SUMAR = ZEROR
            SUMAI = ZEROI
            DO 40 K=1,KMAX
               M = L1 + K
               SUMAR = SUMAR + PR(K)*ALFA(M)
               SUMAI = SUMAI + PI(K)*ALFA(M)
               IF (AP(K).LT.ATOL) GO TO 50
   40       CONTINUE
   50       CONTINUE
            ASUMR = ASUMR + SUMAR*PP
            ASUMI = ASUMI + SUMAI*PP
            IF (PP.LT.TOL) IAS = 1
   60       CONTINUE
            IF (IBS.EQ.1) GO TO 90
            SUMBR = ZEROR
            SUMBI = ZEROI
            DO 70 K=1,KMAX
               M = L2 + K
               SUMBR = SUMBR + PR(K)*BETA(M)
               SUMBI = SUMBI + PI(K)*BETA(M)
               IF (AP(K).LT.ATOL) GO TO 80
   70       CONTINUE
   80       CONTINUE
            BSUMR = BSUMR + SUMBR*PP
            BSUMI = BSUMI + SUMBI*PP
            IF (PP.LT.BTOL) IBS = 1
   90       CONTINUE
            IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
            L1 = L1 + 30
            L2 = L2 + 30
  100    CONTINUE
  110    CONTINUE
         ASUMR = ASUMR + CONER
         PP = RFNU*RFN13
         BSUMR = BSUMR*PP
         BSUMI = BSUMI*PP
  120    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     CABS(W2).GT.0.25D0
C-----------------------------------------------------------------------
  130    CONTINUE
         CALL ZSQRT(W2R, W2I, WR, WI)
         IF (WR.LT.0.0D0) WR = 0.0D0
         IF (WI.LT.0.0D0) WI = 0.0D0
         STR = CONER + WR
         STI = WI
         CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI)
         CALL ZLOG(ZAR, ZAI, ZCR, ZCI, IDUM)
         IF (ZCI.LT.0.0D0) ZCI = 0.0D0
         IF (ZCI.GT.HPI) ZCI = HPI
         IF (ZCR.LT.0.0D0) ZCR = 0.0D0
         ZTHR = (ZCR-WR)*1.5D0
         ZTHI = (ZCI-WI)*1.5D0
         ZETA1R = ZCR*FNU
         ZETA1I = ZCI*FNU
         ZETA2R = WR*FNU
         ZETA2I = WI*FNU
         AZTH = ZABS(ZTHR,ZTHI)
         ANG = THPI
         IF (ZTHR.GE.0.0D0 .AND. ZTHI.LT.0.0D0) GO TO 140
         ANG = HPI
         IF (ZTHR.EQ.0.0D0) GO TO 140
         ANG = DATAN(ZTHI/ZTHR)
         IF (ZTHR.LT.0.0D0) ANG = ANG + GPI
  140    CONTINUE
         PP = AZTH**EX2
         ANG = ANG*EX2
         ZETAR = PP*DCOS(ANG)
         ZETAI = PP*DSIN(ANG)
         IF (ZETAI.LT.0.0D0) ZETAI = 0.0D0
         ARGR = ZETAR*FN23
         ARGI = ZETAI*FN23
         CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI)
         CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI)
         TZAR = ZAR + ZAR
         TZAI = ZAI + ZAI
         CALL ZSQRT(TZAR, TZAI, STR, STI)
         PHIR = STR*RFN13
         PHII = STI*RFN13
         IF (IPMTR.EQ.1) GO TO 120
         RAW = 1.0D0/DSQRT(AW2)
         STR = WR*RAW
         STI = -WI*RAW
         TFNR = STR*RFNU*RAW
         TFNI = STI*RFNU*RAW
         RAZTH = 1.0D0/AZTH
         STR = ZTHR*RAZTH
         STI = -ZTHI*RAZTH
         RZTHR = STR*RAZTH*RFNU
         RZTHI = STI*RAZTH*RFNU
         ZCR = RZTHR*AR(2)
         ZCI = RZTHI*AR(2)
         RAW2 = 1.0D0/AW2
         STR = W2R*RAW2
         STI = -W2I*RAW2
         T2R = STR*RAW2
         T2I = STI*RAW2
         STR = T2R*C(2) + C(3)
         STI = T2I*C(2)
         UPR(2) = STR*TFNR - STI*TFNI
         UPI(2) = STR*TFNI + STI*TFNR
         BSUMR = UPR(2) + ZCR
         BSUMI = UPI(2) + ZCI
         ASUMR = ZEROR
         ASUMI = ZEROI
         IF (RFNU.LT.TOL) GO TO 220
         PRZTHR = RZTHR
         PRZTHI = RZTHI
         PTFNR = TFNR
         PTFNI = TFNI
         UPR(1) = CONER
         UPI(1) = CONEI
         PP = 1.0D0
         BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
         KS = 0
         KP1 = 2
         L = 3
         IAS = 0
         IBS = 0
         DO 210 LR=2,12,2
            LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
            DO 160 K=LR,LRP1
               KS = KS + 1
               KP1 = KP1 + 1
               L = L + 1
               ZAR = C(L)
               ZAI = ZEROI
               DO 150 J=2,KP1
                  L = L + 1
                  STR = ZAR*T2R - T2I*ZAI + C(L)
                  ZAI = ZAR*T2I + ZAI*T2R
                  ZAR = STR
  150          CONTINUE
               STR = PTFNR*TFNR - PTFNI*TFNI
               PTFNI = PTFNR*TFNI + PTFNI*TFNR
               PTFNR = STR
               UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI
               UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI
               CRR(KS) = PRZTHR*BR(KS+1)
               CRI(KS) = PRZTHI*BR(KS+1)
               STR = PRZTHR*RZTHR - PRZTHI*RZTHI
               PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
               PRZTHR = STR
               DRR(KS) = PRZTHR*AR(KS+2)
               DRI(KS) = PRZTHI*AR(KS+2)
  160       CONTINUE
            PP = PP*RFNU2
            IF (IAS.EQ.1) GO TO 180
            SUMAR = UPR(LRP1)
            SUMAI = UPI(LRP1)
            JU = LRP1
            DO 170 JR=1,LR
               JU = JU - 1
               SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU)
               SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU)
  170       CONTINUE
            ASUMR = ASUMR + SUMAR
            ASUMI = ASUMI + SUMAI
            TEST = DABS(SUMAR) + DABS(SUMAI)
            IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180       CONTINUE
            IF (IBS.EQ.1) GO TO 200
            SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI
            SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR
            JU = LRP1
            DO 190 JR=1,LR
               JU = JU - 1
               SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU)
               SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU)
  190       CONTINUE
            BSUMR = BSUMR + SUMBR
            BSUMI = BSUMI + SUMBI
            TEST = DABS(SUMBR) + DABS(SUMBI)
            IF (PP.LT.BTOL .AND. TEST.LT.BTOL) IBS = 1
  200       CONTINUE
            IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210    CONTINUE
  220    CONTINUE
         ASUMR = ASUMR + CONER
         STR = -BSUMR*RFN13
         STI = -BSUMI*RFN13
         CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI)
         GO TO 120
      END
      SUBROUTINE ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     * TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUNI1
C***REFER TO  ZBESI,ZBESK
C
C     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  ZUCHK,ZUNIK,ZUOIK,D1MACH,ZABS
C***END PROLOGUE  ZUNI1
C     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
C    *S2,Y,Z,ZETA1,ZETA2
         EXTERNAL ZABS
         DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC,
     *    CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN,
     *    FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
     *    SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,
     *    ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, D1MACH, ZABS
         INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
         DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3),
     *    CSRR(3), CYR(2), CYI(2)
         DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
C
         NZ = 0
         ND = N
         NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
         CSCL = 1.0D0/TOL
         CRSC = TOL
         CSSR(1) = CSCL
         CSSR(2) = CONER
         CSSR(3) = CRSC
         CSRR(1) = CRSC
         CSRR(2) = CONER
         CSRR(3) = CSCL
         BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
         FN = DMAX1(FNU,1.0D0)
         INIT = 0
         CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,
     *    ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
         IF (KODE.EQ.1) GO TO 10
         STR = ZR + ZETA2R
         STI = ZI + ZETA2I
         RAST = FN/ZABS(STR,STI)
         STR = STR*RAST*RAST
         STI = -STI*RAST*RAST
         S1R = -ZETA1R + STR
         S1I = -ZETA1I + STI
         GO TO 20
   10    CONTINUE
         S1R = -ZETA1R + ZETA2R
         S1I = -ZETA1I + ZETA2I
   20    CONTINUE
         RS1 = S1R
         IF (DABS(RS1).GT.ELIM) GO TO 130
   30    CONTINUE
         NN = MIN0(2,ND)
         DO 80 I=1,NN
            FN = FNU + DBLE(FLOAT(ND-I))
            INIT = 0
            CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,
     *       ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
            IF (KODE.EQ.1) GO TO 40
            STR = ZR + ZETA2R
            STI = ZI + ZETA2I
            RAST = FN/ZABS(STR,STI)
            STR = STR*RAST*RAST
            STI = -STI*RAST*RAST
            S1R = -ZETA1R + STR
            S1I = -ZETA1I + STI + ZI
            GO TO 50
   40       CONTINUE
            S1R = -ZETA1R + ZETA2R
            S1I = -ZETA1I + ZETA2I
   50       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = S1R
            IF (DABS(RS1).GT.ELIM) GO TO 110
            IF (I.EQ.1) IFLAG = 2
            IF (DABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = ZABS(PHIR,PHII)
            RS1 = RS1 + DLOG(APHI)
            IF (DABS(RS1).GT.ELIM) GO TO 110
            IF (I.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0D0) GO TO 60
            IF (I.EQ.1) IFLAG = 3
   60       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF CABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
            S2R = PHIR*SUMR - PHII*SUMI
            S2I = PHIR*SUMI + PHII*SUMR
            STR = DEXP(S1R)*CSSR(IFLAG)
            S1R = STR*DCOS(S1I)
            S1I = STR*DSIN(S1I)
            STR = S2R*S1R - S2I*S1I
            S2I = S2R*S1I + S2I*S1R
            S2R = STR
            IF (IFLAG.NE.1) GO TO 70
            CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 110
   70       CONTINUE
            CYR(I) = S2R
            CYI(I) = S2I
            M = ND - I + 1
            YR(M) = S2R*CSRR(IFLAG)
            YI(M) = S2I*CSRR(IFLAG)
   80    CONTINUE
         IF (ND.LE.2) GO TO 100
         RAST = 1.0D0/ZABS(ZR,ZI)
         STR = ZR*RAST
         STI = -ZI*RAST
         RZR = (STR+STR)*RAST
         RZI = (STI+STI)*RAST
         BRY(2) = 1.0D0/BRY(1)
         BRY(3) = D1MACH(2)
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         C1R = CSRR(IFLAG)
         ASCLE = BRY(IFLAG)
         K = ND - 2
         FN = DBLE(FLOAT(K))
         DO 90 I=3,ND
            C2R = S2R
            C2I = S2I
            S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
            S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
            S1R = C2R
            S1I = C2I
            C2R = S2R*C1R
            C2I = S2I*C1R
            YR(K) = C2R
            YI(K) = C2I
            K = K - 1
            FN = FN - 1.0D0
            IF (IFLAG.GE.3) GO TO 90
            STR = DABS(C2R)
            STI = DABS(C2I)
            C2M = DMAX1(STR,STI)
            IF (C2M.LE.ASCLE) GO TO 90
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1R = S1R*C1R
            S1I = S1I*C1R
            S2R = C2R
            S2I = C2I
            S1R = S1R*CSSR(IFLAG)
            S1I = S1I*CSSR(IFLAG)
            S2R = S2R*CSSR(IFLAG)
            S2I = S2I*CSSR(IFLAG)
            C1R = CSRR(IFLAG)
   90    CONTINUE
  100    CONTINUE
         RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110    CONTINUE
         IF (RS1.GT.0.0D0) GO TO 120
         YR(ND) = ZEROR
         YI(ND) = ZEROI
         NZ = NZ + 1
         ND = ND - 1
         IF (ND.EQ.0) GO TO 100
         CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
         IF (NUF.LT.0) GO TO 120
         ND = ND - NUF
         NZ = NZ + NUF
         IF (ND.EQ.0) GO TO 100
         FN = FNU + DBLE(FLOAT(ND-1))
         IF (FN.GE.FNUL) GO TO 30
         NLAST = ND
         RETURN
  120    CONTINUE
         NZ = -1
         RETURN
  130    CONTINUE
         IF (RS1.GT.0.0D0) GO TO 120
         NZ = N
         DO 140 I=1,N
            YR(I) = ZEROR
            YI(I) = ZEROI
  140    CONTINUE
         RETURN
      END
      SUBROUTINE ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     * TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUNI2
C***REFER TO  ZBESI,ZBESK
C
C     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  ZAIRY,ZUCHK,ZUNHJ,ZUOIK,D1MACH,ZABS
C***END PROLOGUE  ZUNI2
C     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
C    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
         EXTERNAL ZABS
         DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,
     *    ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,
     *    CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,
     *    DAIR, ELIM, FN, FNU, FNUL, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI,
     *    RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI,
     *    ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR,
     *    CYI, D1MACH, ZABS, CAR, SAR
         INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     *    NN, NUF, NW, NZ, IDUM
         DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3),
     *    CSRR(3), CYR(2), CYI(2)
         DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
         DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     *    CIPI(4)/ 1.0D0,0.0D0, 0.0D0,1.0D0, -1.0D0,0.0D0, 0.0D0,-1.0D0/
         DATA HPI, AIC  /
     1         1.57079632679489662D+00,     1.265512123484645396D+00/
C
         NZ = 0
         ND = N
         NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
         CSCL = 1.0D0/TOL
         CRSC = TOL
         CSSR(1) = CSCL
         CSSR(2) = CONER
         CSSR(3) = CRSC
         CSRR(1) = CRSC
         CSRR(2) = CONER
         CSRR(3) = CSCL
         BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
         ZNR = ZI
         ZNI = -ZR
         ZBR = ZR
         ZBI = ZI
         CIDI = -CONER
         INU = INT(SNGL(FNU))
         ANG = HPI*(FNU-DBLE(FLOAT(INU)))
         C2R = DCOS(ANG)
         C2I = DSIN(ANG)
         CAR = C2R
         SAR = C2I
         IN = INU + N - 1
         IN = MOD(IN,4) + 1
         STR = C2R*CIPR(IN) - C2I*CIPI(IN)
         C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
         C2R = STR
         IF (ZI.GT.0.0D0) GO TO 10
         ZNR = -ZNR
         ZBI = -ZBI
         CIDI = -CIDI
         C2I = -C2I
   10    CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
         FN = DMAX1(FNU,1.0D0)
         CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     *    ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
         IF (KODE.EQ.1) GO TO 20
         STR = ZBR + ZETA2R
         STI = ZBI + ZETA2I
         RAST = FN/ZABS(STR,STI)
         STR = STR*RAST*RAST
         STI = -STI*RAST*RAST
         S1R = -ZETA1R + STR
         S1I = -ZETA1I + STI
         GO TO 30
   20    CONTINUE
         S1R = -ZETA1R + ZETA2R
         S1I = -ZETA1I + ZETA2I
   30    CONTINUE
         RS1 = S1R
         IF (DABS(RS1).GT.ELIM) GO TO 150
   40    CONTINUE
         NN = MIN0(2,ND)
         DO 90 I=1,NN
            FN = FNU + DBLE(FLOAT(ND-I))
            CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,
     *       ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
            IF (KODE.EQ.1) GO TO 50
            STR = ZBR + ZETA2R
            STI = ZBI + ZETA2I
            RAST = FN/ZABS(STR,STI)
            STR = STR*RAST*RAST
            STI = -STI*RAST*RAST
            S1R = -ZETA1R + STR
            S1I = -ZETA1I + STI + DABS(ZI)
            GO TO 60
   50       CONTINUE
            S1R = -ZETA1R + ZETA2R
            S1I = -ZETA1I + ZETA2I
   60       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = S1R
            IF (DABS(RS1).GT.ELIM) GO TO 120
            IF (I.EQ.1) IFLAG = 2
            IF (DABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
            APHI = ZABS(PHIR,PHII)
            AARG = ZABS(ARGR,ARGI)
            RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
            IF (DABS(RS1).GT.ELIM) GO TO 120
            IF (I.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0D0) GO TO 70
            IF (I.EQ.1) IFLAG = 3
   70       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
            CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
            CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
            STR = DAIR*BSUMR - DAII*BSUMI
            STI = DAIR*BSUMI + DAII*BSUMR
            STR = STR + (AIR*ASUMR-AII*ASUMI)
            STI = STI + (AIR*ASUMI+AII*ASUMR)
            S2R = PHIR*STR - PHII*STI
            S2I = PHIR*STI + PHII*STR
            STR = DEXP(S1R)*CSSR(IFLAG)
            S1R = STR*DCOS(S1I)
            S1I = STR*DSIN(S1I)
            STR = S2R*S1R - S2I*S1I
            S2I = S2R*S1I + S2I*S1R
            S2R = STR
            IF (IFLAG.NE.1) GO TO 80
            CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 120
   80       CONTINUE
            IF (ZI.LE.0.0D0) S2I = -S2I
            STR = S2R*C2R - S2I*C2I
            S2I = S2R*C2I + S2I*C2R
            S2R = STR
            CYR(I) = S2R
            CYI(I) = S2I
            J = ND - I + 1
            YR(J) = S2R*CSRR(IFLAG)
            YI(J) = S2I*CSRR(IFLAG)
            STR = -C2I*CIDI
            C2I = C2R*CIDI
            C2R = STR
   90    CONTINUE
         IF (ND.LE.2) GO TO 110
         RAZ = 1.0D0/ZABS(ZR,ZI)
         STR = ZR*RAZ
         STI = -ZI*RAZ
         RZR = (STR+STR)*RAZ
         RZI = (STI+STI)*RAZ
         BRY(2) = 1.0D0/BRY(1)
         BRY(3) = D1MACH(2)
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         C1R = CSRR(IFLAG)
         ASCLE = BRY(IFLAG)
         K = ND - 2
         FN = DBLE(FLOAT(K))
         DO 100 I=3,ND
            C2R = S2R
            C2I = S2I
            S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
            S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
            S1R = C2R
            S1I = C2I
            C2R = S2R*C1R
            C2I = S2I*C1R
            YR(K) = C2R
            YI(K) = C2I
            K = K - 1
            FN = FN - 1.0D0
            IF (IFLAG.GE.3) GO TO 100
            STR = DABS(C2R)
            STI = DABS(C2I)
            C2M = DMAX1(STR,STI)
            IF (C2M.LE.ASCLE) GO TO 100
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1R = S1R*C1R
            S1I = S1I*C1R
            S2R = C2R
            S2I = C2I
            S1R = S1R*CSSR(IFLAG)
            S1I = S1I*CSSR(IFLAG)
            S2R = S2R*CSSR(IFLAG)
            S2I = S2I*CSSR(IFLAG)
            C1R = CSRR(IFLAG)
  100    CONTINUE
  110    CONTINUE
         RETURN
  120    CONTINUE
         IF (RS1.GT.0.0D0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
         YR(ND) = ZEROR
         YI(ND) = ZEROI
         NZ = NZ + 1
         ND = ND - 1
         IF (ND.EQ.0) GO TO 110
         CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
         IF (NUF.LT.0) GO TO 140
         ND = ND - NUF
         NZ = NZ + NUF
         IF (ND.EQ.0) GO TO 110
         FN = FNU + DBLE(FLOAT(ND-1))
         IF (FN.LT.FNUL) GO TO 130
C      FN = CIDI
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1R = CIPR(K)
C      S1I = CIPI(K)
C      IF (FN.LT.0.0D0) S1I = -S1I
C      STR = C2R*S1R - C2I*S1I
C      C2I = C2R*S1I + C2I*S1R
C      C2R = STR
         IN = INU + ND - 1
         IN = MOD(IN,4) + 1
         C2R = CAR*CIPR(IN) - SAR*CIPI(IN)
         C2I = CAR*CIPI(IN) + SAR*CIPR(IN)
         IF (ZI.LE.0.0D0) C2I = -C2I
         GO TO 40
  130    CONTINUE
         NLAST = ND
         RETURN
  140    CONTINUE
         NZ = -1
         RETURN
  150    CONTINUE
         IF (RS1.GT.0.0D0) GO TO 140
         NZ = N
         DO 160 I=1,N
            YR(I) = ZEROR
            YI(I) = ZEROI
  160    CONTINUE
         RETURN
      END
      SUBROUTINE ZUNIK(ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR,
     * PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUNIK
C***REFER TO  ZBESI,ZBESK
C
C        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***ROUTINES CALLED  ZDIV,ZLOG,ZSQRT,D1MACH
C***END PROLOGUE  ZUNIK
C     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
C    *ZETA2,ZN,ZR
         DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI,
     *    CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI,
     *    SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R,
     *    ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR, D1MACH
         INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
         DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
         DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
         DATA CON(1), CON(2)  /
     1    3.98942280401432678D-01,  1.25331413731550025D+00 /
         DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1        C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2        C(19), C(20), C(21), C(22), C(23), C(24)/
     3        1.00000000000000000D+00,    -2.08333333333333333D-01,
     4        1.25000000000000000D-01,     3.34201388888888889D-01,
     5       -4.01041666666666667D-01,     7.03125000000000000D-02,
     6       -1.02581259645061728D+00,     1.84646267361111111D+00,
     7       -8.91210937500000000D-01,     7.32421875000000000D-02,
     8        4.66958442342624743D+00,    -1.12070026162229938D+01,
     9        8.78912353515625000D+00,    -2.36408691406250000D+00,
     A        1.12152099609375000D-01,    -2.82120725582002449D+01,
     B        8.46362176746007346D+01,    -9.18182415432400174D+01,
     C        4.25349987453884549D+01,    -7.36879435947963170D+00,
     D        2.27108001708984375D-01,     2.12570130039217123D+02,
     E       -7.65252468141181642D+02,     1.05999045252799988D+03/
         DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1        C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2        C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3       -6.99579627376132541D+02,     2.18190511744211590D+02,
     4       -2.64914304869515555D+01,     5.72501420974731445D-01,
     5       -1.91945766231840700D+03,     8.06172218173730938D+03,
     6       -1.35865500064341374D+04,     1.16553933368645332D+04,
     7       -5.30564697861340311D+03,     1.20090291321635246D+03,
     8       -1.08090919788394656D+02,     1.72772750258445740D+00,
     9        2.02042913309661486D+04,    -9.69805983886375135D+04,
     A        1.92547001232531532D+05,    -2.03400177280415534D+05,
     B        1.22200464983017460D+05,    -4.11926549688975513D+04,
     C        7.10951430248936372D+03,    -4.93915304773088012D+02,
     D        6.07404200127348304D+00,    -2.42919187900551333D+05,
     E        1.31176361466297720D+06,    -2.99801591853810675D+06/
         DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1        C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2        C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3        3.76327129765640400D+06,    -2.81356322658653411D+06,
     4        1.26836527332162478D+06,    -3.31645172484563578D+05,
     5        4.52187689813627263D+04,    -2.49983048181120962D+03,
     6        2.43805296995560639D+01,     3.28446985307203782D+06,
     7       -1.97068191184322269D+07,     5.09526024926646422D+07,
     8       -7.41051482115326577D+07,     6.63445122747290267D+07,
     9       -3.75671766607633513D+07,     1.32887671664218183D+07,
     A       -2.78561812808645469D+06,     3.08186404612662398D+05,
     B       -1.38860897537170405D+04,     1.10017140269246738D+02,
     C       -4.93292536645099620D+07,     3.25573074185765749D+08,
     D       -9.39462359681578403D+08,     1.55359689957058006D+09,
     E       -1.62108055210833708D+09,     1.10684281682301447D+09/
         DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1        C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2        C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3       -4.95889784275030309D+08,     1.42062907797533095D+08,
     4       -2.44740627257387285D+07,     2.24376817792244943D+06,
     5       -8.40054336030240853D+04,     5.51335896122020586D+02,
     6        8.14789096118312115D+08,    -5.86648149205184723D+09,
     7        1.86882075092958249D+10,    -3.46320433881587779D+10,
     8        4.12801855797539740D+10,    -3.30265997498007231D+10,
     9        1.79542137311556001D+10,    -6.56329379261928433D+09,
     A        1.55927986487925751D+09,    -2.25105661889415278D+08,
     B        1.73951075539781645D+07,    -5.49842327572288687D+05,
     C        3.03809051092238427D+03,    -1.46792612476956167D+10,
     D        1.14498237732025810D+11,    -3.99096175224466498D+11,
     E        8.19218669548577329D+11,    -1.09837515608122331D+12/
         DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1        C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2        C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3        1.00815810686538209D+12,    -6.45364869245376503D+11,
     4        2.87900649906150589D+11,    -8.78670721780232657D+10,
     5        1.76347306068349694D+10,    -2.16716498322379509D+09,
     6        1.43157876718888981D+08,    -3.87183344257261262D+06,
     7        1.82577554742931747D+04,     2.86464035717679043D+11,
     8       -2.40629790002850396D+12,     9.10934118523989896D+12,
     9       -2.05168994109344374D+13,     3.05651255199353206D+13,
     A       -3.16670885847851584D+13,     2.33483640445818409D+13,
     B       -1.23204913055982872D+13,     4.61272578084913197D+12,
     C       -1.19655288019618160D+12,     2.05914503232410016D+11,
     D       -2.18229277575292237D+10,     1.24700929351271032D+09/
         DATA C(119), C(120)/
     1       -2.91883881222208134D+07,     1.18838426256783253D+05/
C
         IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
         RFN = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
         TEST = D1MACH(1)*1.0D+3
         AC = FNU*TEST
         IF (DABS(ZRR).GT.AC .OR. DABS(ZRI).GT.AC) GO TO 15
         ZETA1R = 2.0D0*DABS(DLOG(TEST))+FNU
         ZETA1I = 0.0D0
         ZETA2R = FNU
         ZETA2I = 0.0D0
         PHIR = 1.0D0
         PHII = 0.0D0
         RETURN
   15    CONTINUE
         TR = ZRR*RFN
         TI = ZRI*RFN
         SR = CONER + (TR*TR-TI*TI)
         SI = CONEI + (TR*TI+TI*TR)
         CALL ZSQRT(SR, SI, SRR, SRI)
         STR = CONER + SRR
         STI = CONEI + SRI
         CALL ZDIV(STR, STI, TR, TI, ZNR, ZNI)
         CALL ZLOG(ZNR, ZNI, STR, STI, IDUM)
         ZETA1R = FNU*STR
         ZETA1I = FNU*STI
         ZETA2R = FNU*SRR
         ZETA2I = FNU*SRI
         CALL ZDIV(CONER, CONEI, SRR, SRI, TR, TI)
         SRR = TR*RFN
         SRI = TI*RFN
         CALL ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16))
         PHIR = CWRKR(16)*CON(IKFLG)
         PHII = CWRKI(16)*CON(IKFLG)
         IF (IPMTR.NE.0) RETURN
         CALL ZDIV(CONER, CONEI, SR, SI, T2R, T2I)
         CWRKR(1) = CONER
         CWRKI(1) = CONEI
         CRFNR = CONER
         CRFNI = CONEI
         AC = 1.0D0
         L = 1
         DO 20 K=2,15
            SR = ZEROR
            SI = ZEROI
            DO 10 J=1,K
               L = L + 1
               STR = SR*T2R - SI*T2I + C(L)
               SI = SR*T2I + SI*T2R
               SR = STR
   10       CONTINUE
            STR = CRFNR*SRR - CRFNI*SRI
            CRFNI = CRFNR*SRI + CRFNI*SRR
            CRFNR = STR
            CWRKR(K) = CRFNR*SR - CRFNI*SI
            CWRKI(K) = CRFNR*SI + CRFNI*SR
            AC = AC*RFN
            TEST = DABS(CWRKR(K)) + DABS(CWRKI(K))
            IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20    CONTINUE
         K = 15
   30    CONTINUE
         INIT = K
   40    CONTINUE
         IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
         SR = ZEROR
         SI = ZEROI
         DO 50 I=1,INIT
            SR = SR + CWRKR(I)
            SI = SI + CWRKI(I)
   50    CONTINUE
         SUMR = SR
         SUMI = SI
         PHIR = CWRKR(16)*CON(1)
         PHII = CWRKI(16)*CON(1)
         RETURN
   60    CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
         SR = ZEROR
         SI = ZEROI
         TR = CONER
         DO 70 I=1,INIT
            SR = SR + TR*CWRKR(I)
            SI = SI + TR*CWRKI(I)
            TR = -TR
   70    CONTINUE
         SUMR = SR
         SUMI = SI
         PHIR = CWRKR(16)*CON(2)
         PHII = CWRKI(16)*CON(2)
         RETURN
      END
      SUBROUTINE ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUNK1
C***REFER TO  ZBESK
C
C     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  ZKSCL,ZS1S2,ZUCHK,ZUNIK,D1MACH,ZABS
C***END PROLOGUE  ZUNK1
C     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
C    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
         EXTERNAL ZABS
         DOUBLE PRECISION ALIM, ANG, APHI, ASC, ASCLE, BRY, CKI, CKR,
     *    CONER, CRSC, CSCL, CSGNI, CSPNI, CSPNR, CSR, CSRR, CSSR,
     *    CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2M, C2R, ELIM, FMR, FN,
     *    FNF, FNU, PHIDI, PHIDR, PHII, PHIR, PI, RAST, RAZR, RS1, RZI,
     *    RZR, SGN, STI, STR, SUMDI, SUMDR, SUMI, SUMR, S1I, S1R, S2I,
     *    S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R,
     *    ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZI, ZR, ZRI, ZRR, D1MACH, ZABS
         INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     *    KK, KODE, MR, N, NW, NZ, INITD, IC, IPARD, J, M
         DIMENSION BRY(3), INIT(2), YR(N), YI(N), SUMR(2), SUMI(2),
     *    ZETA1R(2), ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2),
     *    CWRKR(16,3), CWRKI(16,3), CSSR(3), CSRR(3), PHIR(2), PHII(2)
         DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
         DATA PI / 3.14159265358979324D0 /
C
         KDFLG = 1
         NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
         CSCL = 1.0D0/TOL
         CRSC = TOL
         CSSR(1) = CSCL
         CSSR(2) = CONER
         CSSR(3) = CRSC
         CSRR(1) = CRSC
         CSRR(2) = CONER
         CSRR(3) = CSCL
         BRY(1) = 1.0D+3*D1MACH(1)/TOL
         BRY(2) = 1.0D0/BRY(1)
         BRY(3) = D1MACH(2)
         ZRR = ZR
         ZRI = ZI
         IF (ZR.GE.0.0D0) GO TO 10
         ZRR = -ZR
         ZRI = -ZI
   10    CONTINUE
         J = 2
         DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
            J = 3 - J
            FN = FNU + DBLE(FLOAT(I-1))
            INIT(J) = 0
            CALL ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT(J), PHIR(J), PHII(J),
     *       ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), SUMR(J), SUMI(J),
     *       CWRKR(1,J), CWRKI(1,J))
            IF (KODE.EQ.1) GO TO 20
            STR = ZRR + ZETA2R(J)
            STI = ZRI + ZETA2I(J)
            RAST = FN/ZABS(STR,STI)
            STR = STR*RAST*RAST
            STI = -STI*RAST*RAST
            S1R = ZETA1R(J) - STR
            S1I = ZETA1I(J) - STI
            GO TO 30
   20       CONTINUE
            S1R = ZETA1R(J) - ZETA2R(J)
            S1I = ZETA1I(J) - ZETA2I(J)
   30       CONTINUE
            RS1 = S1R
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            IF (DABS(RS1).GT.ELIM) GO TO 60
            IF (KDFLG.EQ.1) KFLAG = 2
            IF (DABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = ZABS(PHIR(J),PHII(J))
            RS1 = RS1 + DLOG(APHI)
            IF (DABS(RS1).GT.ELIM) GO TO 60
            IF (KDFLG.EQ.1) KFLAG = 1
            IF (RS1.LT.0.0D0) GO TO 40
            IF (KDFLG.EQ.1) KFLAG = 3
   40       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
            S2R = PHIR(J)*SUMR(J) - PHII(J)*SUMI(J)
            S2I = PHIR(J)*SUMI(J) + PHII(J)*SUMR(J)
            STR = DEXP(S1R)*CSSR(KFLAG)
            S1R = STR*DCOS(S1I)
            S1I = STR*DSIN(S1I)
            STR = S2R*S1R - S2I*S1I
            S2I = S1R*S2I + S2R*S1I
            S2R = STR
            IF (KFLAG.NE.1) GO TO 50
            CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 60
   50       CONTINUE
            CYR(KDFLG) = S2R
            CYI(KDFLG) = S2I
            YR(I) = S2R*CSRR(KFLAG)
            YI(I) = S2I*CSRR(KFLAG)
            IF (KDFLG.EQ.2) GO TO 75
            KDFLG = 2
            GO TO 70
   60       CONTINUE
            IF (RS1.GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
            IF (ZR.LT.0.0D0) GO TO 300
            KDFLG = 1
            YR(I)=ZEROR
            YI(I)=ZEROI
            NZ=NZ+1
            IF (I.EQ.1) GO TO 70
            IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 70
            YR(I-1)=ZEROR
            YI(I-1)=ZEROI
            NZ=NZ+1
   70    CONTINUE
         I = N
   75    CONTINUE
         RAZR = 1.0D0/ZABS(ZRR,ZRI)
         STR = ZRR*RAZR
         STI = -ZRI*RAZR
         RZR = (STR+STR)*RAZR
         RZI = (STI+STI)*RAZR
         CKR = FN*RZR
         CKI = FN*RZI
         IB = I + 1
         IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
         FN = FNU + DBLE(FLOAT(N-1))
         IPARD = 1
         IF (MR.NE.0) IPARD = 0
         INITD = 0
         CALL ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI,
     *    ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, CWRKR(1,3),
     *    CWRKI(1,3))
         IF (KODE.EQ.1) GO TO 80
         STR = ZRR + ZET2DR
         STI = ZRI + ZET2DI
         RAST = FN/ZABS(STR,STI)
         STR = STR*RAST*RAST
         STI = -STI*RAST*RAST
         S1R = ZET1DR - STR
         S1I = ZET1DI - STI
         GO TO 90
   80    CONTINUE
         S1R = ZET1DR - ZET2DR
         S1I = ZET1DI - ZET2DI
   90    CONTINUE
         RS1 = S1R
         IF (DABS(RS1).GT.ELIM) GO TO 95
         IF (DABS(RS1).LT.ALIM) GO TO 100
C----------------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-------------------------------------------------------------------------
         APHI = ZABS(PHIDR,PHIDI)
         RS1 = RS1+DLOG(APHI)
         IF (DABS(RS1).LT.ELIM) GO TO 100
   95    CONTINUE
         IF (DABS(RS1).GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
         IF (ZR.LT.0.0D0) GO TO 300
         NZ = N
         DO 96 I=1,N
            YR(I) = ZEROR
            YI(I) = ZEROI
   96    CONTINUE
         RETURN
C---------------------------------------------------------------------------
C     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
C----------------------------------------------------------------------------
  100    CONTINUE
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         C1R = CSRR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 120 I=IB,N
            C2R = S2R
            C2I = S2I
            S2R = CKR*C2R - CKI*C2I + S1R
            S2I = CKR*C2I + CKI*C2R + S1I
            S1R = C2R
            S1I = C2I
            CKR = CKR + RZR
            CKI = CKI + RZI
            C2R = S2R*C1R
            C2I = S2I*C1R
            YR(I) = C2R
            YI(I) = C2I
            IF (KFLAG.GE.3) GO TO 120
            STR = DABS(C2R)
            STI = DABS(C2I)
            C2M = DMAX1(STR,STI)
            IF (C2M.LE.ASCLE) GO TO 120
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1R = S1R*C1R
            S1I = S1I*C1R
            S2R = C2R
            S2I = C2I
            S1R = S1R*CSSR(KFLAG)
            S1I = S1I*CSSR(KFLAG)
            S2R = S2R*CSSR(KFLAG)
            S2I = S2I*CSSR(KFLAG)
            C1R = CSRR(KFLAG)
  120    CONTINUE
  160    CONTINUE
         IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
         NZ = 0
         FMR = DBLE(FLOAT(MR))
         SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
         CSGNI = SGN
         INU = INT(SNGL(FNU))
         FNF = FNU - DBLE(FLOAT(INU))
         IFN = INU + N - 1
         ANG = FNF*SGN
         CSPNR = DCOS(ANG)
         CSPNI = DSIN(ANG)
         IF (MOD(IFN,2).EQ.0) GO TO 170
         CSPNR = -CSPNR
         CSPNI = -CSPNI
  170    CONTINUE
         ASC = BRY(1)
         IUF = 0
         KK = N
         KDFLG = 1
         IB = IB - 1
         IC = IB - 1
         DO 270 K=1,N
            FN = FNU + DBLE(FLOAT(KK-1))
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
            M=3
            IF (N.GT.2) GO TO 175
  172       CONTINUE
            INITD = INIT(J)
            PHIDR = PHIR(J)
            PHIDI = PHII(J)
            ZET1DR = ZETA1R(J)
            ZET1DI = ZETA1I(J)
            ZET2DR = ZETA2R(J)
            ZET2DI = ZETA2I(J)
            SUMDR = SUMR(J)
            SUMDI = SUMI(J)
            M = J
            J = 3 - J
            GO TO 180
  175       CONTINUE
            IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
            IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
            INITD = 0
  180       CONTINUE
            CALL ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI,
     *       ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI,
     *       CWRKR(1,M), CWRKI(1,M))
            IF (KODE.EQ.1) GO TO 200
            STR = ZRR + ZET2DR
            STI = ZRI + ZET2DI
            RAST = FN/ZABS(STR,STI)
            STR = STR*RAST*RAST
            STI = -STI*RAST*RAST
            S1R = -ZET1DR + STR
            S1I = -ZET1DI + STI
            GO TO 210
  200       CONTINUE
            S1R = -ZET1DR + ZET2DR
            S1I = -ZET1DI + ZET2DI
  210       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = S1R
            IF (DABS(RS1).GT.ELIM) GO TO 260
            IF (KDFLG.EQ.1) IFLAG = 2
            IF (DABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = ZABS(PHIDR,PHIDI)
            RS1 = RS1 + DLOG(APHI)
            IF (DABS(RS1).GT.ELIM) GO TO 260
            IF (KDFLG.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0D0) GO TO 220
            IF (KDFLG.EQ.1) IFLAG = 3
  220       CONTINUE
            STR = PHIDR*SUMDR - PHIDI*SUMDI
            STI = PHIDR*SUMDI + PHIDI*SUMDR
            S2R = -CSGNI*STI
            S2I = CSGNI*STR
            STR = DEXP(S1R)*CSSR(IFLAG)
            S1R = STR*DCOS(S1I)
            S1I = STR*DSIN(S1I)
            STR = S2R*S1R - S2I*S1I
            S2I = S2R*S1I + S2I*S1R
            S2R = STR
            IF (IFLAG.NE.1) GO TO 230
            CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
            IF (NW.EQ.0) GO TO 230
            S2R = ZEROR
            S2I = ZEROI
  230       CONTINUE
            CYR(KDFLG) = S2R
            CYI(KDFLG) = S2I
            C2R = S2R
            C2I = S2I
            S2R = S2R*CSRR(IFLAG)
            S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
            S1R = YR(KK)
            S1I = YI(KK)
            IF (KODE.EQ.1) GO TO 250
            CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  250       CONTINUE
            YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
            YI(KK) = CSPNR*S1I + CSPNI*S1R + S2I
            KK = KK - 1
            CSPNR = -CSPNR
            CSPNI = -CSPNI
            IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
            KDFLG = 1
            GO TO 270
  255       CONTINUE
            IF (KDFLG.EQ.2) GO TO 275
            KDFLG = 2
            GO TO 270
  260       CONTINUE
            IF (RS1.GT.0.0D0) GO TO 300
            S2R = ZEROR
            S2I = ZEROI
            GO TO 230
  270    CONTINUE
         K = N
  275    CONTINUE
         IL = N - K
         IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         CSR = CSRR(IFLAG)
         ASCLE = BRY(IFLAG)
         FN = DBLE(FLOAT(INU+IL))
         DO 290 I=1,IL
            C2R = S2R
            C2I = S2I
            S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
            S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
            S1R = C2R
            S1I = C2I
            FN = FN - 1.0D0
            C2R = S2R*CSR
            C2I = S2I*CSR
            CKR = C2R
            CKI = C2I
            C1R = YR(KK)
            C1I = YI(KK)
            IF (KODE.EQ.1) GO TO 280
            CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  280       CONTINUE
            YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
            YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
            KK = KK - 1
            CSPNR = -CSPNR
            CSPNI = -CSPNI
            IF (IFLAG.GE.3) GO TO 290
            C2R = DABS(CKR)
            C2I = DABS(CKI)
            C2M = DMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 290
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1R = S1R*CSR
            S1I = S1I*CSR
            S2R = CKR
            S2I = CKI
            S1R = S1R*CSSR(IFLAG)
            S1I = S1I*CSSR(IFLAG)
            S2R = S2R*CSSR(IFLAG)
            S2I = S2I*CSSR(IFLAG)
            CSR = CSRR(IFLAG)
  290    CONTINUE
         RETURN
  300    CONTINUE
         NZ = -1
         RETURN
      END
      SUBROUTINE ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUNK2
C***REFER TO  ZBESK
C
C     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  ZAIRY,ZKSCL,ZS1S2,ZUCHK,ZUNHJ,D1MACH,ZABS
C***END PROLOGUE  ZUNK2
C     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
C    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
C    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR
         EXTERNAL ZABS
         DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGDI,
     *    ARGDR, ARGI, ARGR, ASC, ASCLE, ASUMDI, ASUMDR, ASUMI, ASUMR,
     *    BRY, BSUMDI, BSUMDR, BSUMI, BSUMR, CAR, CIPI, CIPR, CKI, CKR,
     *    CONER, CRSC, CR1I, CR1R, CR2I, CR2R, CSCL, CSGNI, CSI,
     *    CSPNI, CSPNR, CSR, CSRR, CSSR, CYI, CYR, C1I, C1R, C2I, C2M,
     *    C2R, DAII, DAIR, ELIM, FMR, FN, FNF, FNU, HPI, PHIDI, PHIDR,
     *    PHII, PHIR, PI, PTI, PTR, RAST, RAZR, RS1, RZI, RZR, SAR, SGN,
     *    STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, YY, ZBI, ZBR, ZEROI,
     *    ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZET1DI, ZET1DR, ZET2DI,
     *    ZET2DR, ZI, ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
         INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     *    KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
         DIMENSION BRY(3), YR(N), YI(N), ASUMR(2), ASUMI(2), BSUMR(2),
     *    BSUMI(2), PHIR(2), PHII(2), ARGR(2), ARGI(2), ZETA1R(2),
     *    ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), CIPR(4),
     *    CIPI(4), CSSR(3), CSRR(3)
         DATA ZEROR,ZEROI,CONER,CR1R,CR1I,CR2R,CR2I /
     1            0.0D0, 0.0D0, 1.0D0,
     1    1.0D0,1.73205080756887729D0 , -0.5D0,-8.66025403784438647D-01 /
         DATA HPI, PI, AIC /
     1        1.57079632679489662D+00,     3.14159265358979324D+00,
     1        1.26551212348464539D+00/
         DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     *    CIPI(4) /
     1     1.0D0,0.0D0 ,  0.0D0,-1.0D0 ,  -1.0D0,0.0D0 ,  0.0D0,1.0D0 /
C
         KDFLG = 1
         NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
         CSCL = 1.0D0/TOL
         CRSC = TOL
         CSSR(1) = CSCL
         CSSR(2) = CONER
         CSSR(3) = CRSC
         CSRR(1) = CRSC
         CSRR(2) = CONER
         CSRR(3) = CSCL
         BRY(1) = 1.0D+3*D1MACH(1)/TOL
         BRY(2) = 1.0D0/BRY(1)
         BRY(3) = D1MACH(2)
         ZRR = ZR
         ZRI = ZI
         IF (ZR.GE.0.0D0) GO TO 10
         ZRR = -ZR
         ZRI = -ZI
   10    CONTINUE
         YY = ZRI
         ZNR = ZRI
         ZNI = -ZRR
         ZBR = ZRR
         ZBI = ZRI
         INU = INT(SNGL(FNU))
         FNF = FNU - DBLE(FLOAT(INU))
         ANG = -HPI*FNF
         CAR = DCOS(ANG)
         SAR = DSIN(ANG)
         C2R = HPI*SAR
         C2I = -HPI*CAR
         KK = MOD(INU,4) + 1
         STR = C2R*CIPR(KK) - C2I*CIPI(KK)
         STI = C2R*CIPI(KK) + C2I*CIPR(KK)
         CSR = CR1R*STR - CR1I*STI
         CSI = CR1R*STI + CR1I*STR
         IF (YY.GT.0.0D0) GO TO 20
         ZNR = -ZNR
         ZBI = -ZBI
   20    CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
         J = 2
         DO 80 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
            J = 3 - J
            FN = FNU + DBLE(FLOAT(I-1))
            CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR(J), PHII(J), ARGR(J),
     *       ARGI(J), ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), ASUMR(J),
     *       ASUMI(J), BSUMR(J), BSUMI(J))
            IF (KODE.EQ.1) GO TO 30
            STR = ZBR + ZETA2R(J)
            STI = ZBI + ZETA2I(J)
            RAST = FN/ZABS(STR,STI)
            STR = STR*RAST*RAST
            STI = -STI*RAST*RAST
            S1R = ZETA1R(J) - STR
            S1I = ZETA1I(J) - STI
            GO TO 40
   30       CONTINUE
            S1R = ZETA1R(J) - ZETA2R(J)
            S1I = ZETA1I(J) - ZETA2I(J)
   40       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = S1R
            IF (DABS(RS1).GT.ELIM) GO TO 70
            IF (KDFLG.EQ.1) KFLAG = 2
            IF (DABS(RS1).LT.ALIM) GO TO 50
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = ZABS(PHIR(J),PHII(J))
            AARG = ZABS(ARGR(J),ARGI(J))
            RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
            IF (DABS(RS1).GT.ELIM) GO TO 70
            IF (KDFLG.EQ.1) KFLAG = 1
            IF (RS1.LT.0.0D0) GO TO 50
            IF (KDFLG.EQ.1) KFLAG = 3
   50       CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
            C2R = ARGR(J)*CR2R - ARGI(J)*CR2I
            C2I = ARGR(J)*CR2I + ARGI(J)*CR2R
            CALL ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM)
            CALL ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM)
            STR = DAIR*BSUMR(J) - DAII*BSUMI(J)
            STI = DAIR*BSUMI(J) + DAII*BSUMR(J)
            PTR = STR*CR2R - STI*CR2I
            PTI = STR*CR2I + STI*CR2R
            STR = PTR + (AIR*ASUMR(J)-AII*ASUMI(J))
            STI = PTI + (AIR*ASUMI(J)+AII*ASUMR(J))
            PTR = STR*PHIR(J) - STI*PHII(J)
            PTI = STR*PHII(J) + STI*PHIR(J)
            S2R = PTR*CSR - PTI*CSI
            S2I = PTR*CSI + PTI*CSR
            STR = DEXP(S1R)*CSSR(KFLAG)
            S1R = STR*DCOS(S1I)
            S1I = STR*DSIN(S1I)
            STR = S2R*S1R - S2I*S1I
            S2I = S1R*S2I + S2R*S1I
            S2R = STR
            IF (KFLAG.NE.1) GO TO 60
            CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
            IF (NW.NE.0) GO TO 70
   60       CONTINUE
            IF (YY.LE.0.0D0) S2I = -S2I
            CYR(KDFLG) = S2R
            CYI(KDFLG) = S2I
            YR(I) = S2R*CSRR(KFLAG)
            YI(I) = S2I*CSRR(KFLAG)
            STR = CSI
            CSI = -CSR
            CSR = STR
            IF (KDFLG.EQ.2) GO TO 85
            KDFLG = 2
            GO TO 80
   70       CONTINUE
            IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
            IF (ZR.LT.0.0D0) GO TO 320
            KDFLG = 1
            YR(I)=ZEROR
            YI(I)=ZEROI
            NZ=NZ+1
            STR = CSI
            CSI =-CSR
            CSR = STR
            IF (I.EQ.1) GO TO 80
            IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 80
            YR(I-1)=ZEROR
            YI(I-1)=ZEROI
            NZ=NZ+1
   80    CONTINUE
         I = N
   85    CONTINUE
         RAZR = 1.0D0/ZABS(ZRR,ZRI)
         STR = ZRR*RAZR
         STI = -ZRI*RAZR
         RZR = (STR+STR)*RAZR
         RZI = (STI+STI)*RAZR
         CKR = FN*RZR
         CKI = FN*RZI
         IB = I + 1
         IF (N.LT.IB) GO TO 180
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
         FN = FNU + DBLE(FLOAT(N-1))
         IPARD = 1
         IF (MR.NE.0) IPARD = 0
         CALL ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,
     *    ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI)
         IF (KODE.EQ.1) GO TO 90
         STR = ZBR + ZET2DR
         STI = ZBI + ZET2DI
         RAST = FN/ZABS(STR,STI)
         STR = STR*RAST*RAST
         STI = -STI*RAST*RAST
         S1R = ZET1DR - STR
         S1I = ZET1DI - STI
         GO TO 100
   90    CONTINUE
         S1R = ZET1DR - ZET2DR
         S1I = ZET1DI - ZET2DI
  100    CONTINUE
         RS1 = S1R
         IF (DABS(RS1).GT.ELIM) GO TO 105
         IF (DABS(RS1).LT.ALIM) GO TO 120
C----------------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-------------------------------------------------------------------------
         APHI = ZABS(PHIDR,PHIDI)
         RS1 = RS1+DLOG(APHI)
         IF (DABS(RS1).LT.ELIM) GO TO 120
  105    CONTINUE
         IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
         IF (ZR.LT.0.0D0) GO TO 320
         NZ = N
         DO 106 I=1,N
            YR(I) = ZEROR
            YI(I) = ZEROI
  106    CONTINUE
         RETURN
  120    CONTINUE
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         C1R = CSRR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 130 I=IB,N
            C2R = S2R
            C2I = S2I
            S2R = CKR*C2R - CKI*C2I + S1R
            S2I = CKR*C2I + CKI*C2R + S1I
            S1R = C2R
            S1I = C2I
            CKR = CKR + RZR
            CKI = CKI + RZI
            C2R = S2R*C1R
            C2I = S2I*C1R
            YR(I) = C2R
            YI(I) = C2I
            IF (KFLAG.GE.3) GO TO 130
            STR = DABS(C2R)
            STI = DABS(C2I)
            C2M = DMAX1(STR,STI)
            IF (C2M.LE.ASCLE) GO TO 130
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1R = S1R*C1R
            S1I = S1I*C1R
            S2R = C2R
            S2I = C2I
            S1R = S1R*CSSR(KFLAG)
            S1I = S1I*CSSR(KFLAG)
            S2R = S2R*CSSR(KFLAG)
            S2I = S2I*CSSR(KFLAG)
            C1R = CSRR(KFLAG)
  130    CONTINUE
  180    CONTINUE
         IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
         NZ = 0
         FMR = DBLE(FLOAT(MR))
         SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
C-----------------------------------------------------------------------
         CSGNI = SGN
         IF (YY.LE.0.0D0) CSGNI = -CSGNI
         IFN = INU + N - 1
         ANG = FNF*SGN
         CSPNR = DCOS(ANG)
         CSPNI = DSIN(ANG)
         IF (MOD(IFN,2).EQ.0) GO TO 190
         CSPNR = -CSPNR
         CSPNI = -CSPNI
  190    CONTINUE
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
         CSR = SAR*CSGNI
         CSI = CAR*CSGNI
         IN = MOD(IFN,4) + 1
         C2R = CIPR(IN)
         C2I = CIPI(IN)
         STR = CSR*C2R + CSI*C2I
         CSI = -CSR*C2I + CSI*C2R
         CSR = STR
         ASC = BRY(1)
         IUF = 0
         KK = N
         KDFLG = 1
         IB = IB - 1
         IC = IB - 1
         DO 290 K=1,N
            FN = FNU + DBLE(FLOAT(KK-1))
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
            IF (N.GT.2) GO TO 175
  172       CONTINUE
            PHIDR = PHIR(J)
            PHIDI = PHII(J)
            ARGDR = ARGR(J)
            ARGDI = ARGI(J)
            ZET1DR = ZETA1R(J)
            ZET1DI = ZETA1I(J)
            ZET2DR = ZETA2R(J)
            ZET2DI = ZETA2I(J)
            ASUMDR = ASUMR(J)
            ASUMDI = ASUMI(J)
            BSUMDR = BSUMR(J)
            BSUMDI = BSUMI(J)
            J = 3 - J
            GO TO 210
  175       CONTINUE
            IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 210
            IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
            CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,
     *       ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,
     *       ASUMDI, BSUMDR, BSUMDI)
  210       CONTINUE
            IF (KODE.EQ.1) GO TO 220
            STR = ZBR + ZET2DR
            STI = ZBI + ZET2DI
            RAST = FN/ZABS(STR,STI)
            STR = STR*RAST*RAST
            STI = -STI*RAST*RAST
            S1R = -ZET1DR + STR
            S1I = -ZET1DI + STI
            GO TO 230
  220       CONTINUE
            S1R = -ZET1DR + ZET2DR
            S1I = -ZET1DI + ZET2DI
  230       CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
            RS1 = S1R
            IF (DABS(RS1).GT.ELIM) GO TO 280
            IF (KDFLG.EQ.1) IFLAG = 2
            IF (DABS(RS1).LT.ALIM) GO TO 240
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
            APHI = ZABS(PHIDR,PHIDI)
            AARG = ZABS(ARGDR,ARGDI)
            RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
            IF (DABS(RS1).GT.ELIM) GO TO 280
            IF (KDFLG.EQ.1) IFLAG = 1
            IF (RS1.LT.0.0D0) GO TO 240
            IF (KDFLG.EQ.1) IFLAG = 3
  240       CONTINUE
            CALL ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM)
            CALL ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM)
            STR = DAIR*BSUMDR - DAII*BSUMDI
            STI = DAIR*BSUMDI + DAII*BSUMDR
            STR = STR + (AIR*ASUMDR-AII*ASUMDI)
            STI = STI + (AIR*ASUMDI+AII*ASUMDR)
            PTR = STR*PHIDR - STI*PHIDI
            PTI = STR*PHIDI + STI*PHIDR
            S2R = PTR*CSR - PTI*CSI
            S2I = PTR*CSI + PTI*CSR
            STR = DEXP(S1R)*CSSR(IFLAG)
            S1R = STR*DCOS(S1I)
            S1I = STR*DSIN(S1I)
            STR = S2R*S1R - S2I*S1I
            S2I = S2R*S1I + S2I*S1R
            S2R = STR
            IF (IFLAG.NE.1) GO TO 250
            CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
            IF (NW.EQ.0) GO TO 250
            S2R = ZEROR
            S2I = ZEROI
  250       CONTINUE
            IF (YY.LE.0.0D0) S2I = -S2I
            CYR(KDFLG) = S2R
            CYI(KDFLG) = S2I
            C2R = S2R
            C2I = S2I
            S2R = S2R*CSRR(IFLAG)
            S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
            S1R = YR(KK)
            S1I = YI(KK)
            IF (KODE.EQ.1) GO TO 270
            CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  270       CONTINUE
            YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
            YI(KK) = S1R*CSPNI + S1I*CSPNR + S2I
            KK = KK - 1
            CSPNR = -CSPNR
            CSPNI = -CSPNI
            STR = CSI
            CSI = -CSR
            CSR = STR
            IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
            KDFLG = 1
            GO TO 290
  255       CONTINUE
            IF (KDFLG.EQ.2) GO TO 295
            KDFLG = 2
            GO TO 290
  280       CONTINUE
            IF (RS1.GT.0.0D0) GO TO 320
            S2R = ZEROR
            S2I = ZEROI
            GO TO 250
  290    CONTINUE
         K = N
  295    CONTINUE
         IL = N - K
         IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
         S1R = CYR(1)
         S1I = CYI(1)
         S2R = CYR(2)
         S2I = CYI(2)
         CSR = CSRR(IFLAG)
         ASCLE = BRY(IFLAG)
         FN = DBLE(FLOAT(INU+IL))
         DO 310 I=1,IL
            C2R = S2R
            C2I = S2I
            S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
            S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
            S1R = C2R
            S1I = C2I
            FN = FN - 1.0D0
            C2R = S2R*CSR
            C2I = S2I*CSR
            CKR = C2R
            CKI = C2I
            C1R = YR(KK)
            C1I = YI(KK)
            IF (KODE.EQ.1) GO TO 300
            CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
            NZ = NZ + NW
  300       CONTINUE
            YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
            YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
            KK = KK - 1
            CSPNR = -CSPNR
            CSPNI = -CSPNI
            IF (IFLAG.GE.3) GO TO 310
            C2R = DABS(CKR)
            C2I = DABS(CKI)
            C2M = DMAX1(C2R,C2I)
            IF (C2M.LE.ASCLE) GO TO 310
            IFLAG = IFLAG + 1
            ASCLE = BRY(IFLAG)
            S1R = S1R*CSR
            S1I = S1I*CSR
            S2R = CKR
            S2I = CKI
            S1R = S1R*CSSR(IFLAG)
            S1I = S1I*CSSR(IFLAG)
            S2R = S2R*CSSR(IFLAG)
            S2I = S2I*CSSR(IFLAG)
            CSR = CSRR(IFLAG)
  310    CONTINUE
         RETURN
  320    CONTINUE
         NZ = -1
         RETURN
      END
      SUBROUTINE ZUOIK(ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL,
     * ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZUOIK
C***REFER TO  ZBESI,ZBESK,ZBESH
C
C     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***ROUTINES CALLED  ZUCHK,ZUNHJ,ZUNIK,D1MACH,ZABS,ZLOG
C***END PROLOGUE  ZUOIK
C     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
C    *ZR
         EXTERNAL ZABS
         DOUBLE PRECISION AARG, AIC, ALIM, APHI, ARGI, ARGR, ASUMI, ASUMR,
     *    ASCLE, AX, AY, BSUMI, BSUMR, CWRKI, CWRKR, CZI, CZR, ELIM, FNN,
     *    FNU, GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, TOL, YI,
     *    YR, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI,
     *    ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
         INTEGER I, IDUM, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
         DIMENSION YR(N), YI(N), CWRKR(16), CWRKI(16)
         DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
         DATA AIC / 1.265512123484645396D+00 /
         NUF = 0
         NN = N
         ZRR = ZR
         ZRI = ZI
         IF (ZR.GE.0.0D0) GO TO 10
         ZRR = -ZR
         ZRI = -ZI
   10    CONTINUE
         ZBR = ZRR
         ZBI = ZRI
         AX = DABS(ZR)*1.7321D0
         AY = DABS(ZI)
         IFORM = 1
         IF (AY.GT.AX) IFORM = 2
         GNU = DMAX1(FNU,1.0D0)
         IF (IKFLG.EQ.1) GO TO 20
         FNN = DBLE(FLOAT(NN))
         GNN = FNU + FNN - 1.0D0
         GNU = DMAX1(GNN,FNN)
   20    CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
         IF (IFORM.EQ.2) GO TO 30
         INIT = 0
         CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     *    ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
         CZR = -ZETA1R + ZETA2R
         CZI = -ZETA1I + ZETA2I
         GO TO 50
   30    CONTINUE
         ZNR = ZRI
         ZNI = -ZRR
         IF (ZI.GT.0.0D0) GO TO 40
         ZNR = -ZNR
   40    CONTINUE
         CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     *    ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
         CZR = -ZETA1R + ZETA2R
         CZI = -ZETA1I + ZETA2I
         AARG = ZABS(ARGR,ARGI)
   50    CONTINUE
         IF (KODE.EQ.1) GO TO 60
         CZR = CZR - ZBR
         CZI = CZI - ZBI
   60    CONTINUE
         IF (IKFLG.EQ.1) GO TO 70
         CZR = -CZR
         CZI = -CZI
   70    CONTINUE
         APHI = ZABS(PHIR,PHII)
         RCZ = CZR
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
         IF (RCZ.GT.ELIM) GO TO 210
         IF (RCZ.LT.ALIM) GO TO 80
         RCZ = RCZ + DLOG(APHI)
         IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
         IF (RCZ.GT.ELIM) GO TO 210
         GO TO 130
   80    CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
         IF (RCZ.LT.(-ELIM)) GO TO 90
         IF (RCZ.GT.(-ALIM)) GO TO 130
         RCZ = RCZ + DLOG(APHI)
         IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
         IF (RCZ.GT.(-ELIM)) GO TO 110
   90    CONTINUE
         DO 100 I=1,NN
            YR(I) = ZEROR
            YI(I) = ZEROI
  100    CONTINUE
         NUF = NN
         RETURN
  110    CONTINUE
         ASCLE = 1.0D+3*D1MACH(1)/TOL
         CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
         CZR = CZR + STR
         CZI = CZI + STI
         IF (IFORM.EQ.1) GO TO 120
         CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
         CZR = CZR - 0.25D0*STR - AIC
         CZI = CZI - 0.25D0*STI
  120    CONTINUE
         AX = DEXP(RCZ)/TOL
         AY = CZI
         CZR = AX*DCOS(AY)
         CZI = AX*DSIN(AY)
         CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
         IF (NW.NE.0) GO TO 90
  130    CONTINUE
         IF (IKFLG.EQ.2) RETURN
         IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  140    CONTINUE
         GNU = FNU + DBLE(FLOAT(NN-1))
         IF (IFORM.EQ.2) GO TO 150
         INIT = 0
         CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     *    ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
         CZR = -ZETA1R + ZETA2R
         CZI = -ZETA1I + ZETA2I
         GO TO 160
  150    CONTINUE
         CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     *    ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
         CZR = -ZETA1R + ZETA2R
         CZI = -ZETA1I + ZETA2I
         AARG = ZABS(ARGR,ARGI)
  160    CONTINUE
         IF (KODE.EQ.1) GO TO 170
         CZR = CZR - ZBR
         CZI = CZI - ZBI
  170    CONTINUE
         APHI = ZABS(PHIR,PHII)
         RCZ = CZR
         IF (RCZ.LT.(-ELIM)) GO TO 180
         IF (RCZ.GT.(-ALIM)) RETURN
         RCZ = RCZ + DLOG(APHI)
         IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
         IF (RCZ.GT.(-ELIM)) GO TO 190
  180    CONTINUE
         YR(NN) = ZEROR
         YI(NN) = ZEROI
         NN = NN - 1
         NUF = NUF + 1
         IF (NN.EQ.0) RETURN
         GO TO 140
  190    CONTINUE
         ASCLE = 1.0D+3*D1MACH(1)/TOL
         CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
         CZR = CZR + STR
         CZI = CZI + STI
         IF (IFORM.EQ.1) GO TO 200
         CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
         CZR = CZR - 0.25D0*STR - AIC
         CZI = CZI - 0.25D0*STI
  200    CONTINUE
         AX = DEXP(RCZ)/TOL
         AY = CZI
         CZR = AX*DCOS(AY)
         CZI = AX*DSIN(AY)
         CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
         IF (NW.NE.0) GO TO 180
         RETURN
  210    CONTINUE
         NUF = -1
         RETURN
      END
      SUBROUTINE ZWRSK(ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI,
     * TOL, ELIM, ALIM)

c*********************************************************************72
c
C***BEGIN PROLOGUE  ZWRSK
C***REFER TO  ZBESI,ZBESK
C
C     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
C
C***ROUTINES CALLED  D1MACH,ZBKNU,ZRATI,ZABS
C***END PROLOGUE  ZWRSK
C     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
         EXTERNAL ZABS
         DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI,
     *    CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT,
     *    STI, STR, TOL, YI, YR, ZRI, ZRR, ZABS, D1MACH
         INTEGER I, KODE, N, NW, NZ
         DIMENSION YR(N), YI(N), CWR(2), CWI(2)
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
         NZ = 0
         CALL ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
         IF (NW.NE.0) GO TO 50
         CALL ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
         CINUR = 1.0D0
         CINUI = 0.0D0
         IF (KODE.EQ.1) GO TO 10
         CINUR = DCOS(ZRI)
         CINUI = DSIN(ZRI)
   10    CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
         ACW = ZABS(CWR(2),CWI(2))
         ASCLE = 1.0D+3*D1MACH(1)/TOL
         CSCLR = 1.0D0
         IF (ACW.GT.ASCLE) GO TO 20
         CSCLR = 1.0D0/TOL
         GO TO 30
   20    CONTINUE
         ASCLE = 1.0D0/ASCLE
         IF (ACW.LT.ASCLE) GO TO 30
         CSCLR = TOL
   30    CONTINUE
         C1R = CWR(1)*CSCLR
         C1I = CWI(1)*CSCLR
         C2R = CWR(2)*CSCLR
         C2I = CWI(2)*CSCLR
         STR = YR(1)
         STI = YI(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0D0/CABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
C-----------------------------------------------------------------------
         PTR = STR*C1R - STI*C1I
         PTI = STR*C1I + STI*C1R
         PTR = PTR + C2R
         PTI = PTI + C2I
         CTR = ZRR*PTR - ZRI*PTI
         CTI = ZRR*PTI + ZRI*PTR
         ACT = ZABS(CTR,CTI)
         RACT = 1.0D0/ACT
         CTR = CTR*RACT
         CTI = -CTI*RACT
         PTR = CINUR*RACT
         PTI = CINUI*RACT
         CINUR = PTR*CTR - PTI*CTI
         CINUI = PTR*CTI + PTI*CTR
         YR(1) = CINUR*CSCLR
         YI(1) = CINUI*CSCLR
         IF (N.EQ.1) RETURN
         DO 40 I=2,N
            PTR = STR*CINUR - STI*CINUI
            CINUI = STR*CINUI + STI*CINUR
            CINUR = PTR
            STR = YR(I)
            STI = YI(I)
            YR(I) = CINUR*CSCLR
            YI(I) = CINUI*CSCLR
   40    CONTINUE
         RETURN
   50    CONTINUE
         NZ = -1
         IF(NW.EQ.(-2)) NZ=-2
         RETURN
      END
