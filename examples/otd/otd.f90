!
!     OTD implementation by Simon Kern (skern@mech.kth.se)
!     Based on the method proposed in Babaee,H. and Sapsis,T.P., "A minimization principle for the
!      description of modes associated with finite-time instabilities", Proc. R. Soc. A 474: 20150779, https://dx.doi.org/10.1098/rspa.2015.0779
!
      subroutine otd_init
         implicit none
         include 'SIZE'
         include 'INPUT' ! param(59), initc
         include 'TSTEP' ! time
         include 'GEOM'  ! [xyz]m1
         include 'OTD'

         ! integer itmp
         ! real rtmp, ltim
         ! logical ltmp
         ! character*20 ctmp
         ! character*2 str1, str2
         ! character*200 lstring
         real xtmp(lx1, ly1, lz1, lelv)
         real ytmp(lx1, ly1, lz1, lelv)
         real ztmp(lx1, ly1, lz1, lelv)

         call white_noise_IC ! Initialise all IC fields to white noise
         if (otd_nusric < 0) then !  Read IC fields
            ! c call mntr_abort(otd_id, 'Choose a valid number of IC files.')
            ! c else if(otd_nusric == 0) then
            ! c call mntr_log(otd_id, lp_inf,
            ! c$'No IC files read. ICs will be white noise.')
            ! c else
            ! c if(otd_nusric > npert) then
            ! c call mntr_abort(otd_id,
            ! c$'otd_nusric > npert! Increase the number of perturbations.')
            ! c endif
            ! save mesh
            call opcopy(xtmp, ytmp, ztmp, xm1, ym1, zm1)
            call read_OTDIC() ! read ICs
            call opcopy(xm1, ym1, zm1, xtmp, ytmp, ztmp) ! restore mesh
         end if

         !   Sanity check
         !   if (otd_ifFTLE) then
         !     if (otd_FTLEpd.eq.0.0) then
         !     ! JUST PRINT THIS
         !     !   call mntr_warn(otd_id,'The FLTEs will be computed continuously. Set otd_FTLEpd in the par-file to define a finite time horizon.')
         !     endif
         !     pcount = 0
         !   endif

         !   Ensure right variable size
         !   if (npert.gt.lpert) then
         !     call mntr_abort(otd_id, 'npert > lpert! Set lpert=npert for the tool to work.')
         !   endif

         call blank(initc, 132) ! set initial conditions 
         call setics

!    Set default timestep at which to start OTD computation we use the first timestep to make sure that the perturbations
!     1. are divergence free
!     2. satisfy the orthonormality constraint
!     3. satisfy the boundary conditions of the problem
!
         startOTD = 1
         return
      end subroutine otd_init
!
      subroutine run_OTD ! Main interface for the OTD module
         implicit none
         include 'SIZE'
         include 'TSTEP' ! istep,iostep
         include 'OTD'
         include 'MASS'

         ! real op_glsc2_wt
         integer ntot!, i
         real iz(lpert)

         if (istep == 0) then ! OTD initialisation
!       Perform orthonormalisation of the basis vectors
!       Initialise fundamental solution matrix for FTLE computation
!       Compute one time step
!
            gsstep_override = .true.        ! initialisation
            call OTD_ON()
!           call outpost_OTDbasis()
            if (otd_ifFTLE) then
               call reset_FTLE
            end if
!
         elseif (istep >= startOTD) then   ! standard OTD

!       Orthonormalize perturbations (if necessary)
            gsstep_override = .false.
            call OTD_ON()

!       Build the action of the linearized NS-operator and compute the reduced operator Lr_{ij}
            call gen_Lu()

!       Compute additional forcing for perturbation equation

            call gen_OTD_forcing()

!       Compute the OTD modes
            if (mod(istep, otd_prntsp) == 0) then
               call proj_OTDmodes()
            end if

!       Output OTD modes
            if (mod(istep, otd_iostep) == 0) then
               call outpost_projmodes()
            end if

!       Output perturbation fields (OTD basis) for restart
            if (mod(istep, otd_rststp) == 0 .or. istep == nsteps .or. lastep == 1) then
               call outpost_OTDbasis()
            end if

!       FTLEs
            if (otd_ifFTLE) then
               call get_FTLE()
            end if
!
         end if ! istep 0 --> init
         return
      end subroutine run_OTD

!=======================================================================

      subroutine gen_Lu() ! Construct action of linearized NS operator on perturbation field
         implicit none
         include 'SIZE'
         include 'INPUT'           ! if3d
         include 'MASS'            ! BM1
         include 'SOLN'            ! V[XYZ]P
         include 'TSTEP'           ! istep
         include 'OTD'             ! conv[xyz],diff[xyz],gradp[xyz]
         ! Lu[xyz], Lr

         ! local variables
         integer ipert, jpert, ntot, i
         real op_ip

!     Build the elements of the linearized NS-operator L_{NS} (u_j)
!
!     L_{NS} (u_j) = 1/Re (grad^2 u)_j - (grad p)_j - (Ub.grad) u_j - (u_j.grad) Ub
!
         do ipert = 1, npert

            ! Convective terms
            call Lu_op_conv(vxp(1, ipert), vyp(1, ipert), vzp(1, ipert), ipert)

            ! Perturbation pressure gradient
            call Lu_op_gradp(prp(1, ipert), ipert)

            ! Diffusive term
            call Lu_op_diff(vxp(1, ipert), vyp(1, ipert), vzp(1, ipert), ipert)

         end do

!     Assemble the action of the operator
!
!     L_{NS} (u_j) = 1/Re grad^2 u_j - grad p - (Ub.grad) u_j - (u_j.grad) Ub
!
         ntot = lx1*ly1*lz1*nelv
         do jpert = 1, npert
            do i = 1, ntot
               Lux(i, jpert) = diffx(i, jpert) - gradpx(i, jpert) - convx(i, jpert)
               Luy(i, jpert) = diffy(i, jpert) - gradpy(i, jpert) - convy(i, jpert)
               if (if3d) then
                  Luz(i, jpert) = diffz(i, jpert) - gradpz(i, jpert) - convz(i, jpert)
               end if
            end do
         end do

!     Compute the innner product < L_{NS}(u_i),u_j > with i,j = 1,...,r
         call rzero(Lr, lpert*lpert)
         do ipert = 1, npert
            do jpert = 1, npert
               Lr(ipert, jpert) = op_ip(ipert, jpert, 2)
            end do
         end do
!
! Debug output: Print reduced operator Lr_{ij}
         if (nid == 0 .and. otd_debug) then
            if (mod(istep, otd_prntsp) == 0) then
               ! write out Lr
               call outmat(Lr, npert, npert, 'Lrmat  ', istep)
            end if    ! otd_prntsp
         end if      ! otd_debug.and.nid.eq.0

!     Here you could add internal rotations phi_rot into the method.
!     This does NOT change the subspace.
!
!       dU / dt = L_{NS}U - U (Lr - phi_rot)
!
         call rzero(phi_rot, lpert*lpert)
!
!     The rotation matrix Phi_ij must be skew-symmetric (but is otherwise
!     arbitrary
!
!       phi_rot(i,j) = -phi_rot(j,i)
!
!     e.g. to obtain an evolution that corresponds to continuously
!     performing Gram-Schmidt on the basis (i.e. turning Lr into a lower
!     triangular matrix), set the rotation matrix to
!
!                  / -<Lu_j,u_i>     j < i
!       phi_rot = {   0              j = i
!                  \  <Lu_j,u_i>     j > i
!
         if (npert > 1) then
            do jpert = 1, npert
               do ipert = jpert + 1, npert
                  phi_rot(ipert, jpert) = op_ip(ipert, jpert, 2)
                  phi_rot(jpert, ipert) = -phi_rot(ipert, jpert)
               end do
            end do
         end if

         if (nid == 0 .and. otd_debug) then
            if (mod(istep, otd_prntsp) == 0) then
               call outmat(Lr, npert, npert, 'Lr-mat', istep)
               call outmat(phi_rot, npert, npert, 'phimat', istep)
            end if
         end if

         call sub2(Lr, phi_rot, lpert*lpert) ! add internal rotation if defined
         if (nid == 0 .and. otd_debug) then
            if (mod(istep, otd_prntsp) == 0) then
               call outmat(Lr, npert, npert, 'Lrpmat', istep)
            end if
         end if

         return
      end subroutine gen_Lu

!=======================================================================
! Compute eigenspectrum of the reduced operator Lr_{ij} and project the velocity
! perturbations onto the eigendirections to obtain the most unstable modes
!
      subroutine proj_OTDmodes()
         implicit none

         include 'SIZE'
         include 'INPUT' ! if3d
         include 'TSTEP' ! istep
         include 'SOLN'  ! v[xyz]p
         include 'OTD'   ! EIG[RI], EV[RL], Lr, OTDmr[xyz]

         ! Local variables
         real :: tmp(lpert, lpert)
         integer :: ntot, i, j
         character(len=32) :: fmtr, fmti

         ! Compute the eigenspectrum of the reduced operator, sort the eigenvalues in decreasing order
         ! LAPACK: compute eigenvalues of Lsym = (Lr + Lr^T) / 2
         call copy(tmp, Lr, lpert*lpert)  ! Save Lr
         Lr = 0.5*(tmp + transpose(tmp))  ! Compute Lsym
         call eig_wrapper(npert, 'r')       ! Compute lambdas
         call sorteigs('Ls ')

         ! ! Output
         ! if (mod(istep, otd_prntsp) == 0 .and. nid == 0) then
         !    write (fmtr, '(A, I0, A)') '(', npert, '(E15.7,1X))'
         !    write (6, 100, ADVANCE='NO') istep, time, 'Ls | Re'
         !    write (6, trim(adjustl(fmtr))) EIGR(1:npert)
         ! end if

         ! Compute eigenvalues of Lr
         call copy(Lr, tmp, lpert*lpert)  ! Restore Lr
         call eig_wrapper(npert, 'r')       ! Compute lambdas
         call sorteigs('Lr ')
         call copy(Lr, tmp, lpert*lpert)  ! Restore Lr

100      format('  [OTD] ', I7, 1x, E14.7, 1x, A8)

         ! Output
         if (mod(istep, otd_prntsp) == 0 .and. nid == 0) then
            write (fmtr, '(A, I0, A)') '(', npert, '(E15.7,1X))'
            write (fmti, '(A, I0, A)') '(', npert, '(I4,1X))'

            write (6, 100, ADVANCE='NO') istep, time, 'Lr | Re '
            ! write (6, trim(adjustl(fmtr))) EIGR(1:npert)
            write (6, 100, ADVANCE='NO') istep, time, 'Lr | Im '
            ! write (6, trim(adjustl(fmtr))) EIGI(1:npert)

            ! Print order for reference
            write (6, 100, ADVANCE='NO') istep, time, 's-idx   '
            ! write (6, trim(adjustl(fmti))) idx(1:npert)

            ! Print out non-zero elements of rotated Lr
            write (fmtr, '(A, I0, A)') '(', npert*(npert + 1)/2, '(E15.7,1X))'
            write (6, 100, ADVANCE='NO') istep, time, 'Lrmat   '
            ! write (6, trim(adjustl(fmtr))) ((Lr(i, j), j=i, npert), i=1, npert)
         end if

         ! Project the perturbation velocity field (OTD basis) onto the eigendirections
         ! of the reduced operator to obtain the most unstable directions
         ntot = lx1*ly1*lz1*lelv
         call mxm(vxp, ntot, EVRr, lpert, OTDmrx, npert)
         call mxm(vxp, ntot, EVRi, lpert, OTDmix, npert)
         call mxm(vyp, ntot, EVRr, lpert, OTDmry, npert)
         call mxm(vyp, ntot, EVRi, lpert, OTDmiy, npert)
         if (if3d) then
            call mxm(vzp, ntot, EVRr, lpert, OTDmrz, npert)
            call mxm(vzp, ntot, EVRi, lpert, OTDmiz, npert)
         end if

         return
      end subroutine proj_OTDmodes

!=======================================================================
! Output the projection of the OTD basis on the eigendirection as fields
!
      subroutine outpost_projmodes()
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'TSTEP'
         include 'OTD'

         ! local variables
         integer ipert
         character*2 str
         character*3 oname

         if (mod(istep, otd_prntsp) /= 0) call proj_OTDmodes()
         do ipert = 1, npert
            if (lpert >= 10) then
               write (str, '(I2.2)') ipert
               oname = 'o'//trim(str)
            else
               write (str, '(I1)') ipert
               oname = 'ot'//trim(str)
            end if
            call outpost(OTDmrx(1, ipert), OTDmry(1, ipert), OTDmrz(1, ipert), prp(1, ipert), t, oname)
         end do

         return
      end subroutine outpost_projmodes

!=======================================================================
! Output the OTD basis directly to restart.
!     We could alternatively reconstruct the OTD basis from the modes
!     but for this we would need both real and imaginary part. Since we
!     currently only outpost the real part, it's cheaper to just outpost
!     the OTD basis directly when we also outpost the baseflow.
!
      subroutine outpost_OTDbasis()
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'TSTEP'
         include 'OTD'

         ! local variables
         integer ipert
         character*2 str
         character*3 oname

         ! orthonormalize
         gsstep_override = .true.
         call OTD_ON()

         do ipert = 1, npert
            write (str, '(I2.2)') ipert
            oname = 'r'//trim(str)
            call outpost(vxp(1, ipert), vyp(1, ipert), vzp(1, ipert), prp(1, ipert), t, oname)
         end do

         return
      end subroutine outpost_OTDbasis

!=======================================================================
! Initialise all perturbation fields to random noise
!
      subroutine white_noise_IC()
         implicit none
         include 'SIZE'
         include 'INPUT'           ! if3d
         include 'NEKUSE'
         include 'GEOM'            ! [xyz]m1
         include 'PARALLEL'        ! lglel
         include 'OTD'
         integer ix, iy, iz, ie, ieg, i, ijke
         real mth_rand, xl(LDIM), fcoeff(3)
         do i = 1, npert
            do ix = 1, lx1
            do iy = 1, ly1
            do iz = 1, lz1
               do ie = 1, lelv
                  xl(1) = XM1(ix, iy, iz, ie)
                  xl(2) = YM1(ix, iy, iz, ie)
                  if (if3d) then
                     xl(NDIM) = ZM1(ix, iy, iz, ie)
                  end if
                  ijke = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(ie - 1)))
                  ieg = lglel(ie)
                  fcoeff(1) = sin(real(i))**2*3.0e4
                  fcoeff(2) = sin(real(i))**2*(-1.5e3)
                  fcoeff(3) = sin(real(i))**2*0.5e5
                  vxpic(ijke, i) = mth_rand(ix, iy, iz, ieg, xl, fcoeff)
                  fcoeff(1) = sin(real(i))**2*2.3e4
                  fcoeff(2) = sin(real(i))**2*2.3e3
                  fcoeff(3) = sin(real(i))**2*(-2.0e5)
                  vypic(ijke, i) = mth_rand(ix, iy, iz, ieg, xl, fcoeff)
                  if (if3d) then
                     fcoeff(1) = sin(real(i))**2*2.e4
                     fcoeff(2) = sin(real(i))**2*1.e3
                     fcoeff(3) = sin(real(i))**2*1.e5
                     vzpic(ijke, i) = mth_rand(ix, iy, iz, ieg, xl, fcoeff)
                  end if
               end do
            end do
            end do
            end do
         end do
         return
      end subroutine white_noise_IC

!=======================================================================
! Read OTD IC fields and run setics
!
      subroutine read_OTDIC()
         implicit none
         include 'SIZE'
         include 'INPUT'
         include 'SOLN'
         include 'OTD'
         include 'TSTEP'

         ! local variables
         logical exist_IC
         integer i
         character*2 istr
         character*132 ifile

         !   call mntr_log(otd_id,lp_inf,'INIT - read ICs')
         do i = 1, otd_nusric
            write (istr, '(I0.2)') i
            ifile = 'OTDIC_'//trim(istr)//'.fld'
            inquire (file=ifile, exist=exist_IC)
            if (exist_IC) then
               !   call mntr_log(OTD_id,lp_inf,'INIT - Reading IC file '//trim(ifile))
               call load_fld(ifile)
               ! Copy the initial conditions into the fields v[xyz]pic
               call opcopy(vxpic(1, i), vypic(1, i), vzpic(1, i), vx, vy, vz)
               OTDrsttime = time
               !   call mntr_logr(otd_id,lp_inf,'IC file '//trim(ifile)//' time:',time)
               ! else
               !   call mntr_abort(otd_id,'Cannot open '//trim(ifile)//' !')
            end if
         end do
         !   call mntr_log(otd_id,lp_inf,'INIT - read ICs :: done')

         return
      end subroutine read_OTDIC

!=======================================================================
! Create forcing for OTD evolution equation
!
      subroutine gen_OTD_forcing()
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'OTD'
         include 'INPUT'   ! if3d
         include 'TSTEP'

         integer ntot
         ntot = lx1*ly1*lz1*lelv
         call mxm(VXP, ntot, Lr, lpert, OTDfx, npert)
         call mxm(VYP, ntot, Lr, lpert, OTDfy, npert)
         call mxm(VZP, ntot, Lr, lpert, OTDfz, npert)

         return
      end subroutine gen_OTD_forcing

!=======================================================================
! Set forcing for OTD evolution equation
!     This function is called in userf for each GLL point. Therefore we
!     need a manual switch for when to start setting the forcing.
!
      subroutine set_OTD_forcing(FFX, FFY, FFZ, ix, iy, iz, ieg)
         implicit none
         include 'SIZE'
         include 'SOLN'           ! jp
         include 'TSTEP'          ! istep
         include 'PARALLEL'       ! gllel
         include 'OTD'

         integer ix, iy, iz, ieg
         real ffx, ffy, ffz
         integer ijke, e
         e = gllel(ieg)
         ijke = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
         if (jp /= 0) then ! only for the perturbations
            FFX = FFX - OTDfx(ijke, jp)
            FFY = FFY - OTDfy(ijke, jp)
            FFZ = FFZ - OTDfz(ijke, jp)
         end if

         return
      end subroutine set_OTD_forcing

!=======================================================================

      subroutine OTD_ON() ! Orthonormalize OTD basis
         implicit none
         include 'SIZE'
         include 'OTD'
         include 'SOLN'    ! v[xyz]p
         include 'TSTEP'   ! istep

         ! local variables
         real N, O
         logical runON!
         call compute_NO(N, O, 'pre   ', .true.)

         runON = .false.
         if (otd_gsstep /= 0) then                 ! gsstep=0 => no GS
            if (mod(istep, otd_gsstep) == 0) then
               runON = .true.
            end if
         end if
         if (gsstep_override) runON = .true.       ! override when needed

!     runON
         if (runON) then
!           call CGS()        ! Classical Gram-Schmidt
            call MGS()        ! Modified Gram-Schmidt
         end if
         call compute_NO(N, O, 'post  ', .false.)

         return
      end subroutine OTD_ON

!=======================================================================

      subroutine reset_FTLE() ! Reset FTLE computation
         implicit none
         include 'SIZE'
         include 'OTD'

         call rzero(FTLEv, lpert)
         ! Phi
         !call ident(Phi,lpert)
         ! Blanchard
         call rzero(trapz, lpert)

      end subroutine reset_FTLE

!=======================================================================
! Compute FTLEs (Based on Babaee et al., 2017)
!
!     1. Advect fundamental solution matrix Phi
!
!       dPhi/dt = Lr    , with Phi(t=t0) = I(rxr)
!
!     2. Compute FTLEs
!
!       FTLE(i) = 1/T * log(svd(Phi(t)))
!
      subroutine compute_FTLE_Phi(Lrp, Lrc, deltat, ftledt)
         implicit none

         include 'SIZE'
         include 'TSTEP'
         include 'OTD'

         ! argument list
         real deltat                       ! time interval for FTLE computation
         real ftledt                       ! dt for Phi advection
         real Lrp(lpert, lpert)          ! Lr matrix from previous step
         real Lrc(lpert, lpert)          ! Lr matrix for current step

         ! local variables
         real lhs(lpert, lpert)          ! lhs of Phi advection equation
     $   , invlhs(lpert, lpert)          ! inverse of lhs
     $   , tmp(lpert, lpert)
     $   , rhs(lpert, lpert)          ! rhs of Phi advection equation
         real fact
         integer i

         ! Advect the fundamental solution matrix (using the implicit CN scheme)
         fact = 0.5*ftledt
         ! build LHS
         call ident(lhs, lpert)
         call add2s2(lhs, Lrc, -fact, lpert*lpert)
         ! build RHS
         call ident(tmp, npert)
         call add2s2(tmp, Lrp, fact, lpert*lpert)
         call mxm(tmp, lpert, Phi, lpert, rhs, lpert)
         ! invert LHS and solve system
         call invmt(lhs, invlhs, tmp, npert)
         call mxm(invlhs, lpert, rhs, lpert, Phi, lpert)

         ! compute SVD
         call copy(VMATX, Phi, lpert*lpert)
         call svd_wrapper(lpert, lpert, 'A')

         ! compute FTLEs
         do i = 1, npert
            FTLEv(i) = log(OSIGMA(i))/deltat
         end do

      end subroutine compute_FTLE_Phi

!=======================================================================
! Compute FTLEs (Based on Blanchard & Sapsis, 2019)
!
!       FTLE(i) = 1/T * ( int_(t_0)^t <Lu_i,u_i> d tau )
!
      subroutine compute_FTLE_Blanchard(Lrp, Lrc, deltat, ftledt)
         implicit none

         include 'SIZE'
         include 'TSTEP'
         include 'OTD'

         ! argument list
         real deltat                       ! time interval for FTLE computation
         real ftledt                       ! dt for Phi advection
         real Lrp(lpert, lpert)           ! Lr matrix from previous step
         real Lrc(lpert, lpert)           ! Lr matrix for current step

         ! local variables
         integer i

         ! compute FTLEs
         do i = 1, npert !! MOVE TO EXT3 and add lag fields
            ! Trapezoid rule, this could be done better using the AB scheme
            trapz(i) = trapz(i) + 0.5*ftledt*(Lrp(i, i) + Lrc(i, i))
            FTLEv(i) = trapz(i)/deltat
         end do

      end subroutine compute_FTLE_Blanchard

!=======================================================================
      subroutine get_FTLE()
         include 'SIZE'
         include 'TSTEP'
         include 'OTD' ! Lr, FTLEv

         ! local variables
         real pfrac                ! current fraction of the FTLE comp. period
         real ftledt               ! dt for FTLE computation
         real Lrp(lpert, lpert)     ! Lr from previous step
         real Lrc(lpert, lpert)     ! (linear approx.) of Lr at end of int. interval
         real fact
         real t0
         integer icalld
         ! save variables
         save Lrp
         data t0/0.0d0/
         save t0
         data icalld/0/
         save icalld

         ! determine FTLE horizon
         if (otd_FTLEpd == 0.0d0) then
            if (icalld == 0) then
               t0 = time
               icalld = 1
            end if
            period = time - t0
            pfrac = time - t0
         else
            period = otd_FTLEpd
            pfrac = mod(time, period)
         end if
         ! initialisation
         if (istep == startOTD) then
            call copy(Lrp, Lr, lpert*lpert)
         end if

         if (pfrac < dt) then
            ftledt = dt - pfrac
            ! compute approximation of Lrc at end of period (linear interp.)
            call copy(Lrc, Lrp, lpert*lpert)
            fact = ftledt/dt
            call add2s2(Lrc, Lrp, -fact, lpert*lpert)
            call add2s2(Lrc, Lr, fact, lpert*lpert)
            !call compute_FTLE_Phi(Lrp,Lrc,period,ftledt)
            call compute_FTLE_Blanchard(Lrp, Lrc, period, ftledt)
            pcount = pcount + 1
            if (nid == 0) then
               write (6, 100, ADVANCE='NO') pcount, istep, time
               do i = 1, npert
                  write (6, 102, ADVANCE='NO') FTLEv(i)
               end do
               write (6, *)
            end if
            call reset_FTLE
            call copy(Lrp, Lrc, lpert*lpert)
            call copy(Lrc, Lr, lpert*lpert)
            !call compute_FTLE_Phi(Lrp,Lrc,period,pfrac)
            call compute_FTLE_Blanchard(Lrp, Lrc, period, pfrac)
         else
            call copy(Lrc, Lr, lpert*lpert)
            !call compute_FTLE_Phi(Lrp,Lrc,pfrac,dt)
            call compute_FTLE_Blanchard(Lrp, Lrc, pfrac, dt)
         end if
         ! update Lrp
         call copy(Lrp, Lr, lpert*lpert)

         if (nid == 0 .and. mod(istep, otd_prntsp) == 0) then
            write (6, 101, ADVANCE='NO') istep, time, pfrac
            do i = 1, npert
               write (6, 102, ADVANCE='NO') FTLEv(i)
            end do
            write (6, *)
         end if
100      format(' [OTD] FTLE PRD', 1x, I5, 1x, I7, ' t=', 1x, E15.7)
101      format(' [OTD] FTLE (t)', 1x, I7, ' t=', 2(1x, E15.7))
102      format(1x, E15.7)

      end subroutine get_FTLE

!======================================================================
!Construct the convective terms for Lu
!
!     Lu_conv = (u.grad) Ub + (Ub.grad) u
!
!     Using the convop routine takes care of the dealiasing.
!
      subroutine Lu_op_conv(uxp, uyp, uzp, ipert)
         implicit none

         include 'SIZE'
         include 'INPUT'           ! if3d
         include 'SOLN'            ! v[xyz]
         include 'OTD'             ! conv[xyz]

         ! argument list
         real uxp(lx1*ly1*lz1*lelv)   ! perturbation velocity components
     $   , uyp(lx1*ly1*lz1*lelv)
     $   , uzp(lx1*ly1*lz1*lelv)
         integer ipert                 ! index of the considered perturbation

         ! local variables
         integer ntot, i
         real TA1(LX1, LY1, LZ1, LELV)
     $   , TA2(LX1, LY1, LZ1, LELV)
     $   , TA3(LX1, LY1, LZ1, LELV)
     $   , TB1(LX1, LY1, LZ1, LELV)
     $   , TB2(LX1, LY1, LZ1, LELV)
     $   , TB3(LX1, LY1, LZ1, LELV)

         ntot = lx1*ly1*lz1*lelv
!
         if (if3d) then
            call opcopy(tb1, tb2, tb3, vx, vy, vz)         ! Save velocity
            call opcopy(vx, vy, vz, uxp, uyp, uzp)         ! U <-- u
!       convop(conv,fld): builds the convective term for the scalar field fld
!       conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
            call convop(ta1, tb1)                      ! (u.grad) Ub
            call convop(ta2, tb2)
            call convop(ta3, tb3)
            ! Copy fields into the correct variables
            call opcopy(convx(1, ipert), convy(1, ipert), convz(1, ipert), ta1, ta2, ta3)
            call opcopy(vx, vy, vz, tb1, tb2, tb3)         ! Restore velocity
!       conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
            call convop(tb1, uxp)                      ! (Ub.grad) u
            call convop(tb2, uyp)
            call convop(tb3, uzp)
            ! Add fields to the convective term
            call opadd2(convx(1, ipert), convy(1, ipert), convz(1, ipert), tb1, tb2, tb3)

!   DIAGNOSTICS
!        call opcopy  (convx1(1,ipert),convy1(1,ipert),convz1(1,ipert)
!     $          ,              ta1,ta2,ta3)
!        call opcopy  (convx2(1,ipert),convy2(1,ipert),convz2(1,ipert)
!     $          ,              tb1,tb2,tb3)
         else ! 2D
            call opcopy(tb1, tb2, tb3, vx, vy, vz)         ! Save velocity
            call opcopy(vx, vy, vz, uxp, uyp, uzp)         ! U <-- u
! convop(conv,fld): builds the convective term for the scalar field fld
!       conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
            call convop(ta1, tb1)                      ! (u.grad) Ub
            call convop(ta2, tb2)
            ! Copy fields into the correct variables
            call opcopy(convx(1, ipert), convy(1, ipert), ta3, ta1, ta2, ta3)
            call opcopy(vx, vy, vz, tb1, tb2, tb3)         ! Restore velocity
!       conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
            call convop(tb1, uxp)                      ! (Ub.grad) u
            call convop(tb2, uyp)
            ! Add fields to the convective term
            call opadd2(convx(1, ipert), convy(1, ipert), tb3, tb1, tb2, tb3)
!   DIAGNOSTICS
!        call opcopy  (convx1(1,ipert),convy1(1,ipert),ta3,ta1,ta2,ta3)
!        call opcopy  (convx2(1,ipert),convy2(1,ipert),tb1,tb1,tb2,tb3)
         end if ! if3d

         return
      end subroutine Lu_op_conv

!======================================================================
!Construct the pressure gradient term for Lu
!
!     Lu_gradp = grad p
!
!     Note: The pressure gradient is computed directly on the v-mesh!
!
      subroutine Lu_op_gradp(prpert, ipert)
         implicit none

         include 'SIZE'
         include 'OTD'             ! gradp[xyz]

         ! argument list
         real prpert(lx2*ly2*lz2*lelv, 1) ! perturbation pressure field
         integer ipert                    ! number of the considered pert.

         ! local variables
         integer ntot
         real ta1(lx1, ly1, lz1, lelv), ta2(lx1, ly1, lz1, lelv), wrk(lx1, ly1, lz1, lelv)

         ntot = lx2*ly2*lz2*lelv
!     Map the perturbation pressure to the velocity mesh
         call mappr(wrk, prpert, ta1, ta2)
!     compute the gradient on the velocity mesh directly
         call gradm1(gradpx(1, ipert), gradpy(1, ipert), gradpz(1, ipert), wrk)

         return
      end subroutine Lu_op_gradp

!======================================================================
!Construct the diffusive term for Lu
!
!     Lu_op_diff = 1/Re*grad^2 u
!
      subroutine Lu_op_diff(uxp, uyp, uzp, ipert)
         implicit none

         include 'SIZE'
         include 'INPUT'           ! if3d
         include 'SOLN'            ! vdiff
         include 'OTD'             ! diff[xyz]

         ! argument list
         real uxp(lx1*ly1*lz1*lelv, 1) ! perturbation velocity components
     $   , uyp(lx1*ly1*lz1*lelv, 1)
     $   , uzp(lx1*ly1*lz1*lelv, 1)
         integer ipert, ntot! number of the considered pert.
         ntot = lx1*ly1*lz1*lelv
         ! compute laplacian
         call laplacian(diffx(1, ipert), uxp)
         call laplacian(diffy(1, ipert), uyp)
         if (if3d) call laplacian(diffz(1, ipert), uzp)
         ! multiply by 1/Re                > remove for operator diagnostics
         call col2(diffx(1, ipert), vdiff, ntot)
         call col2(diffy(1, ipert), vdiff, ntot)
         if (if3d) call col2(diffz(1, ipert), vdiff, ntot)

         return
      end subroutine Lu_op_diff

!======================================================================
!Construct the diffusion term (laplacian of u) for direction i
!
      subroutine laplacian(lapu, up)
         implicit none
!
         include 'SIZE'
         include 'INPUT'           ! if3d
         include 'DXYZ'            ! dxm1,d[xy]tm1
         include 'GEOM'            ! r[xy]m1,s[xy]m1,t[xy]m1,jacmi
!
         ! argument list
         real up(lx1*ly1*lz1*lelv, 1)       ! perturbation velocity component
!
         ! output
         real lapu(lx1*ly1*lz1, lelv)
!
         ! local variables
         real ux(lx1*ly1*lz1, lelv)
     $   , uy(lx1*ly1*lz1, lelv)
     $   , uz(lx1*ly1*lz1, lelv)
     $   , ur(lx1*ly1*lz1)
     $   , us(lx1*ly1*lz1)
     $   , ut(lx1*ly1*lz1)
!
         common/ctmp1/ur, us, ut
!
         integer e, i, lxyz, nel, N

         lxyz = lx1*ly1*lz1
         nel = nx1 - 1
         call gradm1(ux, uy, uz, up)
         do e = 1, lelt
            if (if3d) then
               call local_grad3(ur, us, ut, ux, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  lapu(i, e) = jacmi(i, e)*(ur(i)*rxm1(i, 1, 1, e)
     $   +us(i)*sxm1(i, 1, 1, e)
     $   +ut(i)*txm1(i, 1, 1, e))
               end do
               call local_grad3(ur, us, ut, uy, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  lapu(i, e) = lapu(i, e) + jacmi(i, e)*(ur(i)*rym1(i, 1, 1, e)
     $   +us(i)*sym1(i, 1, 1, e)
     $   +ut(i)*tym1(i, 1, 1, e))
               end do
               call local_grad3(ur, us, ut, uz, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  lapu(i, e) = lapu(i, e) + jacmi(i, e)*(ur(i)*rzm1(i, 1, 1, e)
     $   +us(i)*szm1(i, 1, 1, e)
     $   +ut(i)*tzm1(i, 1, 1, e))
               end do
            else ! 2D
               call local_grad2(ur, us, ux, nel, e, dxm1, dytm1)
               do i = 1, lxyz
                  lapu(i, e) = jacmi(i, e)*(ur(i)*rxm1(i, 1, 1, e)
     $   +us(i)*sxm1(i, 1, 1, e))
               end do
               call local_grad2(ur, us, uy, nel, e, dxm1, dytm1)
               do i = 1, lxyz
                  lapu(i, e) = lapu(i, e)
     $   +jacmi(i, e)*(ur(i)*rym1(i, 1, 1, e)
     $   +us(i)*sym1(i, 1, 1, e))
               end do
            end if ! if3d
         end do
!
         return
      end subroutine laplacian
!
!======================================================================
!Compute the global inner product with the pert. velocities
!
!     < vc_i, v[xyz]p_j >
!
      real function op_ip_vp(vcx, vcy, vcz, jpert)
         implicit none

         include 'SIZE'
         include 'SOLN'            ! v[xyz]p, jp
         include 'TSTEP'           ! ifield
         include 'MASS'            ! bm1

         ! argument list
         real vcx(lx1*ly1*lz1*lelv)
         real vcy(lx1*ly1*lz1*lelv)
         real vcz(lx1*ly1*lz1*lelv)
         integer jpert
         ! functions and local variables
         real op_glsc2_wt

         ifield = 1
         op_ip_vp = 0.5*op_glsc2_wt(VCX, VCY, VCZ,
     $   VXP(1, jpert), VYP(1, jpert), VZP(1, jpert),
     $   bm1)

         return
      end function op_ip_vp
!
!======================================================================
!Compute the global inner product for the perturbations
!
!     < v[xyz]p_i , v[xyz]p_j >  or  < LU[xyz]_i , v[xyz]p_j >
!
      real function op_ip(ipert, jpert, iflag)
         implicit none

         include 'SIZE'
         include 'SOLN'            ! v[xyz]p, jp
         include 'OTD'             ! LU[xyz], otd_id

         ! argument list
         integer ipert, jpert       ! [IN]  perturbation indices
         integer iflag             ! [IN]  iflag = 1: < U_i, U_j >
         !                                 iflag = 2: < LU_i, U_j >
         ! functions and local variables
         real op_ip_vp

         if (iflag == 1) then
            op_ip = op_ip_vp(VXP(1, ipert), VYP(1, ipert), VZP(1, ipert), jpert)
         elseif (iflag == 2) then
            op_ip = op_ip_vp(LUx(1, ipert), LUy(1, ipert), LUz(1, ipert), jpert)
            ! else
            !    call mntr_abort(otd_id, 'Error in op_ip!')
         end if

         return
      end function op_ip

!======================================================================

      subroutine op_norm(uxp, uyp, uzp) ! Normalize vector field u_i = v_i/||v_i||
         implicit none

         include 'SIZE'
         include 'TSTEP'           ! ifield
         include 'MASS'            ! bm1
         include 'INPUT'           ! if3d
         include 'OTD'

         ! argument list
         real uxp(lx1*ly1*lz1*lelv, 1) ! perturbation velocity components
     $   , uyp(lx1*ly1*lz1*lelv, 1)
     $   , uzp(lx1*ly1*lz1*lelv, 1)

         ! functions and local variables
         real invnorm, n2, ntot
         real op_glsc2_wt

         ifield = 1
         ntot = lx1*ly1*lz1*lelv
         n2 = 0.5*op_glsc2_wt(uxp, uyp, uzp, uxp, uyp, uzp, bm1)
         ! if (n2 <= 0.0) then
         !    call mntr_abort(otd_id, 'Error in op_norm!')
         ! end if
         invnorm = 1./sqrt(n2)
         call cmult(uxp, invnorm, ntot)
         call cmult(uyp, invnorm, ntot)
         if (if3d) call cmult(uzp, invnorm, ntot)

         return
      end subroutine op_norm
!
!======================================================================
!Perform Classical Gram-Schmidt (CGS) orthonormalization on the
!         perturbation velocity field
!
!     u_k = v_k - sum_{j = 1 -> k-1} proj_{u_j} (u_k)
!
!        with proj_{u_j} (u_k) = < u_k , u_j >/||u_j|| * u_j
!
      subroutine CGS()
         implicit none
!
         include 'SIZE'
         include 'INPUT'   ! if3d
         include 'TSTEP'   ! istep
         include 'SOLN'    ! V[XYZ]P
         include 'OTD'     ! otd_id
!
         ! local variables
         integer i, j, ntot
         real invnorm, proj

         ! functions
         real op_ip

         ntot = lx1*ly1*lz1*lelv
         ! orthonormalize
         do i = 1, npert
            do j = 1, i - 1
               proj = -op_ip(i, j, 1)/op_ip(j, j, 1)
               call add2s2(vxp(1, i), vxp(1, j), proj, ntot)
               call add2s2(vyp(1, i), vyp(1, j), proj, ntot)
               if (if3d) then
                  call add2s2(vzp(1, i), vzp(1, j), proj, ntot)
               end if
            end do
            invnorm = 1/sqrt(op_ip(i, i, 1))
            call cmult(vxp(1, i), invnorm, ntot)
            call cmult(vyp(1, i), invnorm, ntot)
            if (if3d) then
               call cmult(vzp(1, i), invnorm, ntot)
            end if
         end do

         return
      end subroutine CGS
!
!======================================================================
!Perform Modified Gram-Schmidt (MGS) orthonormalization on the
!         perturbation velocity field for improved numerical stability
!
!     do i=1,npert
!       u_i = v_i/||v_i||
!       do j=i+1,npert
!         u_j = v_j - proj_{u_i} (v_j)
!       enddo
!     enddo
!
!        with proj_{u_i} (v_j) = < v_j , u_i >/||u_i|| * v_j
!                              = < v_j , u_i > * v_j   since ||u_i|| = 1
!
      subroutine MGS()
         implicit none
!
         include 'SIZE'
         include 'INPUT'   ! if3d
         include 'TSTEP'   ! istep
         include 'SOLN'    ! V[XYZ]P
         include 'OTD'     ! otd_id
!
         ! local variables
         integer i, j, ntot
         real invnorm, proj

         ! functions
         real op_ip

         ntot = lx1*ly1*lz1*lelv
         ! orthonormalize
         do i = 1, npert
            invnorm = 1/sqrt(op_ip(i, i, 1))
            call cmult(vxp(1, i), invnorm, ntot)
            call cmult(vyp(1, i), invnorm, ntot)
            if (if3d) then
               call cmult(vzp(1, i), invnorm, ntot)
            end if
            do j = i + 1, npert
               proj = -op_ip(i, j, 1)
               call add2s2(vxp(1, j), vxp(1, i), proj, ntot)
               call add2s2(vyp(1, j), vyp(1, i), proj, ntot)
               if (if3d) then
                  call add2s2(vzp(1, j), vzp(1, i), proj, ntot)
               end if
            end do
         end do

         return
      end subroutine MGS
!
!======================================================================
!Compute measures of orthogonality and normality of the
!         perturbations
!
      subroutine compute_NO(N, O, info, flag)
         implicit none

         include 'SIZE'
         include 'TSTEP'
         include 'OTD'

         ! argument list
         real N, O
         logical flag
         character*6 info

         ! local variables
         integer i, j, n9
         real ip(npert, npert)

         ! function
         real op_ip

         do i = 1, npert
            do j = 1, npert
               ip(i, j) = op_ip(i, j, 1)
            end do
         end do

         if (nid == 0) then
            ! compute measure for basis vector normality
            N = 0
            do i = 1, npert
               N = N + ip(i, i)**2
            end do
            N = sqrt(N/npert)
            ! compute measure for basis vector orthogonality
            if (npert > 1) then
               O = 0
               do i = 1, npert
                  do j = i + 1, npert
                     O = O + ip(i, j)**2
                  end do
               end do
               O = sqrt(2*O)/(npert*(npert - 1))
            end if
            if (otd_debug) then ! debug output pre MGS
               n9 = min(npert, 9)
               write (6, *) ' [OTD] ipout istp ipert < u_i, u_j > ON '//trim(info)
               do i = 1, npert
                  write (6, 100) istep, i, (ip(i, j), j=1, n9)
               end do
            end if
            if (flag) write (6, 101) istep, N - 1.0, O
         end if
100      format('  [OTD] ipout', I7, 1x, I5, 2x, 9(1x, F9.6))
101      format('  [OTD] NOout', I7, 1x, 'N-1.0', 1x, E15.7, 1x, 'O', 1x, E15.7)

         return
      end subroutine compute_NO

!======================================================================
! 1. Sort the eigenvalues l_i such that their real parts are ranked in decreasing order
!
!    Re(l_1) .ge. Re(l_i) .ge. Re(l_r), i = 1,...,r
!
!  2. Apply the same sorting to the columns of the right
!     eigenvector matrix and separate real and imaginary parts
!
!    EVR => EVRR + i*EVRI
!
      subroutine sorteigs(str)
         implicit none
!
         include 'SIZE'
         include 'OTD'     ! EIG[IR], EVR, EVR[RI]
         include 'TSTEP'

         ! argument list
         character*3 str
!
         ! local variables
         integer i, j, id
         logical mk(lpert)
         real wrk1(lpert)
         real wrk2(lpert)

         ! zero out indices, output and mask
         call izero(idx, lpert)
         call rzero(EVRR, lpert*lpert)
         call rzero(EVRI, lpert*lpert)
         do i = 1, npert
            mk(i) = .true.
         end do
         call copy(wrk1, EIGR, lpert)
         call copy(wrk2, EIGI, lpert)

         ! we need to exclude the trailing zeros for the sorting to work
         if (lpert > npert) then
            do i = npert + 1, lpert
               mk(i) = .false.
            end do
         end if
         if (otd_debug) then
            if (nid == 0) then
               call print_eigenvalues(' [OTD debug] e-vals: '//trim(str), npert, EIGR, EIGI)
               call print_eigenvectors(' [OTD debug] right e-vecs: '//trim(str), npert, EIGI, EVR, npert)
            end if
         end if
         ! sorting
         j = 1
         do while (j <= npert)
            EIGR(j) = maxval(wrk1, mask=mk)       ! find largest real eigenvalue in remaining list
            id = maxloc(wrk1, 1, mk)          ! find its index
            EIGI(j) = wrk2(id)                   ! extract corresponding imaginary part
            mk(id) = .false.                    ! update mask
            idx(j) = id
            if (EIGI(j) == 0.0) then
               do i = 1, npert
                  EVRR(i, j) = EVR(i, id)
                  ! EVRI(i,j) = 0.0
               end do
               j = j + 1
            else  ! complex conjugate eigenvectors!
               EIGI(j + 1) = -EIGI(j)
               EIGR(j + 1) = EIGR(j)
               mk(id + 1) = .false.
               idx(j + 1) = id + 1
               do i = 1, npert
                  EVRR(i, j) = EVR(i, id)
                  EVRR(i, j + 1) = EVR(i, id)
                  EVRI(i, j) = EVR(i, id + 1)
                  EVRI(i, j + 1) = -EVR(i, id + 1)
               end do
               j = j + 2
            end if
         end do

         ! debugging output
         if (otd_debug) then
            if (nid == 0) then
               call print_eigenvalues(' [OTD debug] sorted e-vals: '//trim(str), npert, EIGR, EIGI)
               write (6, *)
               write (6, *) ' [OTD debug] sorted right e-vecs: '//trim(str)
               do i = 1, npert
                  do j = 1, npert
                     write (6, 200, ADVANCE='NO') EVRR(i, j), EVRI(i, j)
                  end do
                  write (6, *)
               end do
               write (6, *)
            end if
         end if
200      format(9(:, 3x, F6.2, ' + i*', F6.2))

         return
      end subroutine sorteigs


