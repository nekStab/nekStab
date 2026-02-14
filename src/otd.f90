      !-----------------------------------------------------------------------
      ! otd.f90 — Optimally Time-Dependent (OTD) subspace tracking
      !
      ! Purpose:
      !   Implements the OTD method for tracking dominant instability
      !   directions in time-varying flows, including FTLE computation,
      !   basis orthonormalization, and eigenvalue decomposition.
      !
      ! Public interface:
      !   otd, otd_construct_linear_operator, otd_compute_OTD_modes,
      !   otd_outpost_OTD_modes, otd_outpost_orthonormal_basis,
      !   otd_white_noise, otd_generate_forces, otd_orthonormalize_basis,
      !   otd_zero_FTLE, otd_compute_FTLE, otd_construct_convective_terms,
      !   otd_construct_pressure_gradient_terms,
      !   otd_construct_diffusive_terms, otd_compute_laplacian,
      !   otd_op_glsc2, otd_inner_product, otd_normalize_vector_field,
      !   otd_mod_Gram_Schmidt, otd_compute_orthonormality_measures,
      !   otd_sort_eigenvalues, otd_compute_eig_wrapper
      !
      ! Dependencies:
      !   SIZE, TOTAL, INPUT, MASS, SOLN, TSTEP, NEKUSE, GEOM,
      !   PARALLEL, DXYZ
      !-----------------------------------------------------------------------

      module nekstab_otd
         use krylov_subspace
         use nekstab_io
         use nekstab_noise
         implicit none
         include 'SIZE'
         private

      ! OTD variables come from NEKSTAB.inc (via include 'SIZE').
      ! Common blocks: OTD_modes, OTD_params, OTD_oper, OTD_Lu,
      !   OTD_FTLE, OTD_dgeev, rlpvec, ilpvec, clpvec.
      ! Variables are NOT re-exported; use include 'SIZE' to access.

      ! ── Public subroutines ──
         public :: otd, otd_construct_linear_operator,
     $             otd_compute_OTD_modes,
     $             otd_outpost_OTD_modes,
     $             otd_outpost_orthonormal_basis,
     $             otd_white_noise, otd_generate_forces,
     $             otd_orthonormalize_basis,
     $             otd_zero_FTLE, otd_compute_FTLE,
     $             otd_construct_convective_terms,
     $             otd_construct_pressure_gradient_terms,
     $             otd_construct_diffusive_terms,
     $             otd_compute_laplacian,
     $             otd_op_glsc2, otd_inner_product,
     $             otd_normalize_vector_field,
     $             otd_mod_Gram_Schmidt,
     $             otd_compute_orthonormality_measures,
     $             otd_gram_matrix,
     $             otd_sort_eigenvalues,
     $             otd_compute_eig_wrapper
      contains

      !-----------------------------------------------------------------------
      ! otd — main OTD driver (initialization and time-stepping)
      !
      ! Implementation by Simon Kern (skern@mech.kth.se)
      ! Based on Babaee & Sapsis (2016) https://dx.doi.org/10.1098/rspa.2015.0779
      !-----------------------------------------------------------------------
      subroutine otd
         implicit none
         include 'SIZE'
         include 'TOTAL'
         logical, save :: init
         logical :: exist_IC
         data init/.false./
         integer i, j, maxnum
         character(len=3) :: istr
         character(len=30) :: filename

         if (.not. init) then

            ifpert = .true.; call bcast(ifpert, lsize)
            param(31) = lpert; npert = int(param(31))

            ifotd = .true.; call bcast(ifotd, lsize)

            if (nid == 0) then
               open (unit=457, file='otd_growth_rates.dat', status='replace', form='formatted'); close (457)
               open (unit=458, file='otd_eigenvalues.dat', status='replace', form='formatted'); close (458)
               open (unit=460, file='otd_residuals.dat', status='replace', form='formatted'); close (460)
            end if

            write (filename, '(A,A,A)') 'BF_', trim(SESSION), '0.f00001'
            inquire(file=filename, exist=exist_IC)
            if (exist_IC) then
               if (nid == 0) write (6, *) 'OTD: loading base flow from ', trim(filename)
               call load_fld(filename)
               call opcopy(ubic, vbic, wbic, vx, vy, vz)
            else
               if (nid == 0) write (6, *) 'OTD: no BF_ file found, using useric as base flow'
               call opcopy(ubic, vbic, wbic, vx, vy, vz)
            end if
            call outpost(ubic, vbic, wbic, pr, t, 'bf0')

            call otd_white_noise ! Initialise all IC fields to white noise

            ! Search for existing IC files
            maxnum = 0
            do j = 1, 99
               write (filename, '(A,A,A,I5.5)') 'r01', trim(SESSION), '0.f', j
               inquire (file=filename, exist=exist_IC)
               if (exist_IC) then
                  maxnum = j
               end if
            end do
            if (maxnum > 0) then
               if (nid == 0) write (6, *) 'Found ', maxnum, ' IC files'
               do i = 1, npert
                  write (istr, '(A1,I2.2)') 'r', i
                  write (filename, '(2A,I5.5)') trim(istr), trim(SESSION)//'0.f', maxnum
                  if (nid == 0) write (6, *) 'Looking for ICs in ', trim(filename)
                  inquire (file=filename, exist=exist_IC)
                  if (exist_IC) then
                     call load_fld(filename)
                     call opcopy(vxpic(1, i), vypic(1, i), vzpic(1, i), vx, vy, vz)
                  end if
               end do
            else
               if (nid == 0) write (6, *) 'No IC files found.. using white noise'
            end if

            call outpost(ubic, vbic, wbic, pr, t, 'bf0')
            call opcopy(vx, vy, vz, ubic, vbic, wbic) ! restore baseflow

            call blank(initc, 132) ! set initial conditions
            call setics

            time = 0.0d0 ! set time to zero

            call outpost(ubic, vbic, wbic, pr, t, 'ip0')
            do i = 1, npert
               write (istr, '(I1)') i
               call outpost(upic(1, i), vpic(1, i), wpic(1, i), pr, t, 'ip'//trim(istr))
            end do

            if (uparam(1) > 5) then
               ifbase = .true.; call bcast(ifbase, lsize)
            elseif (uparam(1) == 5) then
               if (nid == 0) write (6, *) 'OTD in Frozen baseflow mode!'
               ifbase = .false.; call bcast(ifbase, lsize)
            end if

            gsstep_override = .true.
            call otd_orthonormalize_basis
            call otd_zero_FTLE ! zero Phi matrix

            init = .true.
         end if ! init

      !  Time-stepping: always runs (including istep=0 after init,
      !  which initialises t0=0 in otd_compute_FTLE).
         gsstep_override = .false.
         call otd_orthonormalize_basis
         call otd_construct_linear_operator
         call otd_generate_forces
         call otd_compute_OTD_modes
         if (ifoutfld) then
      ! call otd_outpost_OTD_modes ! M*
            call otd_outpost_orthonormal_basis ! r*
         end if
         call otd_compute_FTLE

      end subroutine otd

      !-----------------------------------------------------------------------
      ! otd_construct_linear_operator — build action of linearized NS
      !   operator on perturbation field
      !-----------------------------------------------------------------------
      subroutine otd_construct_linear_operator
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d
         include 'MASS' ! BM1
         include 'SOLN' ! V[XYZ]P
         include 'TSTEP' ! istep
         integer ipert, jpert, nv, ldv, i

      !     Build the elements of the linearized NS-operator L_{NS} (u_j)
      !         L_{NS} (u_j) = 1/Re (grad^2 u)_j - (grad p)_j - (Ub.grad) u_j - (u_j.grad) Ub

         do ipert = 1, npert
            call otd_construct_convective_terms(vxp(1, ipert), vyp(1, ipert), vzp(1, ipert), ipert)
            call otd_construct_pressure_gradient_terms(prp(1, ipert), ipert)
            call otd_construct_diffusive_terms(vxp(1, ipert), vyp(1, ipert), vzp(1, ipert), ipert)
         end do

      !     Assemble the action of the operator
      !
      !     L_{NS} (u_j) = 1/Re grad^2 u_j - grad p - (Ub.grad) u_j - (u_j.grad) Ub
      !
         nv = lx1*ly1*lz1*nelv
         ldv = lx1*ly1*lz1*lelv
         do jpert = 1, npert
            do i = 1, nv
               otd_Lux(i, jpert) = diffx(i, jpert) - gradpx(i, jpert) - convx(i, jpert)
               otd_Luy(i, jpert) = diffy(i, jpert) - gradpy(i, jpert) - convy(i, jpert)
               if (if3d) otd_Luz(i, jpert) = diffz(i, jpert) - gradpz(i, jpert) - convz(i, jpert)
            end do
         end do

      !     Compute Lr(i,j) = <L_NS(u_i), u_j> via batched dgemm.
      !     Weight otd_Lu by mass matrix in-place (safe: not reused).
         do jpert = 1, npert
            call col2(otd_Lux(1, jpert), bm1, nv)
            call col2(otd_Luy(1, jpert), bm1, nv)
            if (if3d) call col2(otd_Luz(1, jpert), bm1, nv)
         end do
         call otd_gram_matrix(otd_Lr,
     $        otd_Lux, otd_Luy, otd_Luz,
     $        vxp, vyp, vzp, nv, ldv)

      !     --- Internal rotation matrix phi_rot ---
      !     Skew-symmetric: phi_rot(i,j) = -phi_rot(j,i)
      !     Extracted from Lr (same inner products, already computed).
         call rzero(phi_rot, lpert*lpert)
         if (npert > 1) then
            do jpert = 1, npert
               do ipert = jpert + 1, npert
                  phi_rot(ipert, jpert) =  otd_Lr(ipert, jpert)
                  phi_rot(jpert, ipert) = -otd_Lr(ipert, jpert)
               end do
            end do
         end if
         call sub2(otd_Lr, phi_rot, lpert*lpert) ! add internal rotation if defined

      end subroutine otd_construct_linear_operator

      !-----------------------------------------------------------------------
      ! otd_compute_OTD_modes — eigenspectrum of reduced operator and
      !   projection onto eigendirections
      !-----------------------------------------------------------------------
      subroutine otd_compute_OTD_modes

      ! Compute eigenspectrum of the reduced operator Lr_{ij} and project the velocity
      !       perturbations onto the eigendirections to obtain the most unstable modes

         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d
         include 'TSTEP' ! istep
         include 'SOLN' ! v[xyz]p

         real tmp(lpert, lpert)
         integer nv, i, j
         character(len=20) fmtr, fmti

         call copy(tmp, otd_Lr, lpert*lpert) ! Save otd_Lr in tmp
         do i = 1, npert
            do j = 1, npert
               otd_Lr(i, j) = 0.50d0*(tmp(i, j) + tmp(j, i)) ! Compute Lsym = (otd_Lr+otd_Lr^T)/2
            end do
         end do
         call otd_compute_eig_wrapper(npert, 'r') ! Compute lambdas from otd_Lr
         call otd_sort_eigenvalues('otd_Ls') ! Sort lambdas

         if (nid == 0) then ! save to file
            open (unit=457, file='otd_growth_rates.dat', position='append', status='unknown', form='formatted')
            write (fmtr, '("(",I0,"(E15.7,1X))")') npert + 1
            write (457, fmtr) time, (EIGR(i), i=1, npert)
            close (457)
         end if

         if (mod(istep, otd_printStep) == 0 .and. nid == 0) then ! print out
            write (fmtr, '("(",I0,"(E15.7,1X))")') npert
            if (nid == 0) write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Ls | Re '
            write (6, fmtr) (EIGR(i), i=1, npert)
         end if

         call copy(otd_Lr, tmp, lpert*lpert) ! Restore otd_Lr
         call otd_compute_eig_wrapper(npert, 'r') ! Compute lambdas eigenvalues of otd_Lr
         call otd_sort_eigenvalues('otd_Lr ')
         call copy(otd_Lr, tmp, lpert*lpert) ! Restore otd_Lr

         if (nid == 0) then ! save to file
            open (unit=458, file='otd_eigenvalues.dat', position='append', status='unknown', form='formatted')
            write (fmtr, '("(",I0,"(E15.7,1X))")') npert + 1
            write (458, fmtr) time, (EIGR(i), i=1, npert)
            close (458)
         end if

         if (mod(istep, otd_printStep) == 0 .and. nid == 0) then ! print out
            write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Lr | Re '
            write (6, fmtr) (EIGR(i), i=1, npert)

            write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Lr | Im '
            write (6, fmtr) (EIGI(i), i=1, npert) ! print order for reference
            write (fmti, '("(",I0,"(I4,1X))")') npert

            write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 's-otd_idx   '
            write (6, fmti) (otd_idx(i), i=1, npert) ! print out non-zero elements of rotated otd_Lr

            write (fmtr, '("(",I0,"(E15.7,1X))")') npert*(npert + 1)/2
            write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Lrmat   '
            write (6, fmtr) ((otd_Lr(i, j), j=i, npert), i=1, npert)

         end if

      !   Project the perturbation velocity field (OTD basis) onto the eigendirections of
      !   the reduced operator to obtain the most unstable directions

         nv = lx1*ly1*lz1*lelv
         call mxm(vxp, nv, EVRr, lpert, OTDmrx, npert)
         call mxm(vxp, nv, EVRi, lpert, OTDmix, npert)
         call mxm(vyp, nv, EVRr, lpert, OTDmry, npert)
         call mxm(vyp, nv, EVRi, lpert, OTDmiy, npert)
         if (if3d) then
            call mxm(vzp, nv, EVRr, lpert, OTDmrz, npert)
            call mxm(vzp, nv, EVRi, lpert, OTDmiz, npert)
         end if

      end subroutine otd_compute_OTD_modes

      !-----------------------------------------------------------------------
      ! otd_outpost_OTD_modes — output projection of orthonormal basis
      !   on eigendirections
      !-----------------------------------------------------------------------
      subroutine otd_outpost_OTD_modes
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'TSTEP'

         integer ipert
         character(len=2) :: str
         character(len=3) :: oname

         call otd_compute_OTD_modes
         do ipert = 1, npert
            if (lpert >= 10) then
               write (str, '(I2.2)') ipert
               oname = 'M'//trim(str)
            else
               write (str, '(I1)') ipert
               oname = 'M'//trim(str)//'_'
            end if
            call outpost(OTDmrx(1, ipert), OTDmry(1, ipert), OTDmrz(1, ipert), prp(1, ipert), t, oname)

         end do

      end subroutine otd_outpost_OTD_modes

      !-----------------------------------------------------------------------
      ! otd_outpost_orthonormal_basis — output the OTD basis directly
      !   for restart
      !-----------------------------------------------------------------------
      subroutine otd_outpost_orthonormal_basis
      ! Output the OTD basis directly to restart.
      !     We could alternatively reconstruct the OTD basis from the modes
      !     but for this we would need both real and imaginary part. Since we
      !     currently only outpost the real part, it's cheaper to just outpost
      !     the OTD basis directly when we also outpost the baseflow.
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'TSTEP'

         integer ipert
         character(len=2) :: str
         character(len=3) :: oname

         gsstep_override = .true.
         call otd_orthonormalize_basis ! orthonormalize

         do ipert = 1, npert
            write (str, '(I2.2)') ipert
            oname = 'r'//trim(str)
            call outpost(vxp(1, ipert), vyp(1, ipert), vzp(1, ipert), prp(1, ipert), t, oname)
         end do

      end subroutine otd_outpost_orthonormal_basis

      !-----------------------------------------------------------------------
      ! otd_white_noise — initialize all IC fields to white noise
      !
      ! Purpose:
      !   Fills vxpic/vypic/vzpic with deterministic pseudo-random noise.
      !   Each mode gets a DISTINCT pattern because sin(i)^2 scales the
      !   frequency coefficients BEFORE hashing (not after), so the hash
      !   function produces genuinely different fields per mode.
      !   DSS averaging and velocity BCs are applied after generation.
      !-----------------------------------------------------------------------
      subroutine otd_white_noise
         use nekstab_noise, only: mth_rand
         implicit none
         include 'SIZE'
         include 'INPUT'
         include 'GEOM'
         include 'PARALLEL'
         include 'SOLN'

         integer :: i, ix, iy, iz, ie, ieg, nv, ijke
         real :: xl(ldim), fc(3), sin2
         real :: glmin, glmax, nmin, nmax

         nv = nx1*ny1*nz1*nelv

         do i = 1, npert
            sin2 = sin(real(i))**2
            call rzero(vxpic(1, i), nv)
            call rzero(vypic(1, i), nv)
            if (if3d) call rzero(vzpic(1, i), nv)

            do ie = 1, nelv
               do iz = 1, nz1
                  do iy = 1, ny1
                     do ix = 1, nx1
                        xl(1) = xm1(ix, iy, iz, ie)
                        xl(2) = ym1(ix, iy, iz, ie)
                        if (if3d) xl(ndim) = zm1(ix, iy, iz, ie)
                        ijke = ix + lx1*((iy - 1)
     $                       + ly1*((iz - 1) + lz1*(ie - 1)))
                        ieg = lglel(ie)

                        fc(1) = sin2*3.0e4
                        fc(2) = sin2*(-1.5e3)
                        fc(3) = sin2*0.5e5
                        vxpic(ijke, i) =
     $                       mth_rand(ix, iy, iz, ieg, xl, fc)

                        fc(1) = sin2*2.3e4
                        fc(2) = sin2*2.3e3
                        fc(3) = sin2*(-2.0e5)
                        vypic(ijke, i) =
     $                       mth_rand(ix, iy, iz, ieg, xl, fc)

                        if (if3d) then
                           fc(1) = sin2*2.0e4
                           fc(2) = sin2*1.0e3
                           fc(3) = sin2*1.0e5
                           vzpic(ijke, i) =
     $                          mth_rand(ix, iy, iz, ieg, xl, fc)
                        end if
                     end do
                  end do
               end do
            end do

      !     DSS averaging + velocity BCs
            call opdssum(vxpic(1, i), vypic(1, i), vzpic(1, i))
            call opcolv(vxpic(1, i), vypic(1, i),
     $           vzpic(1, i), vmult)
            call dsavg(vxpic(1, i))
            call dsavg(vypic(1, i))
            if (if3d) call dsavg(vzpic(1, i))
            call bcdirVC(vxpic(1, i), vypic(1, i),
     $           vzpic(1, i), v1mask, v2mask, v3mask)
         end do

         nmin = glmin(vxpic, nv)
         nmax = glmax(vxpic, nv)
         if (nid == 0) write (6, *) 'OTD noise min,max', nmin, nmax

      end subroutine otd_white_noise

      !-----------------------------------------------------------------------
      ! otd_generate_forces — create forcing for OTD evolution equation
      !-----------------------------------------------------------------------
      subroutine otd_generate_forces
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'INPUT'
         include 'TSTEP'

         integer nv; nv = lx1*ly1*lz1*lelv
         call mxm(VXP, nv, otd_Lr, lpert, OTDfx, npert)
         call mxm(VYP, nv, otd_Lr, lpert, OTDfy, npert)
         if (if3d) call mxm(VZP, nv, otd_Lr, lpert, OTDfz, npert)

      end subroutine otd_generate_forces

      !-----------------------------------------------------------------------
      ! otd_orthonormalize_basis — orthonormalize perturbation basis
      !-----------------------------------------------------------------------
      subroutine otd_orthonormalize_basis
         implicit none
         include 'SIZE'
         include 'SOLN' ! v[xyz]p
         include 'TSTEP' ! istep

         real :: N, O
         logical :: runON

         call otd_compute_orthonormality_measures(N, O, 'pre   ', .true.)

         runON = .false.
         if (otd_gsStep /= 0) then ! gsstep=0 => no GS
            if (mod(istep, otd_gsStep) == 0) then
               runON = .true.
            end if
         end if
         if (gsstep_override) runON = .true. ! override when needed

         if (runON) then
      !call otd_classic_Gram_Schmidt    ! Classical Gram-Schmidt
            call otd_mod_Gram_Schmidt ! Modified Gram-Schmidt
         end if

      end subroutine otd_orthonormalize_basis

      !-----------------------------------------------------------------------
      ! otd_zero_FTLE — reset FTLE computation
      !-----------------------------------------------------------------------
      subroutine otd_zero_FTLE
         implicit none
         include 'SIZE'

         call rzero(FTLEv, lpert)
         call rzero(LEintegral, lpert) ! Blanchard

      end subroutine otd_zero_FTLE

      !-----------------------------------------------------------------------
      ! otd_compute_FTLE — compute finite-time Lyapunov exponents
      !-----------------------------------------------------------------------
      subroutine otd_compute_FTLE
      !
      !  Compute finite-time Lyapunov exponents via trapezoidal
      !  quadrature of the diagonal of L_r:
      !
      !     lambda_i(T) = (1/T) * integral_0^T  L_r(i,i) dt
      !
      !  When otd_FTLEPeriod > 0, the integral resets every period.
      !  Otherwise it accumulates from the first call onward.
      !
         implicit none
         include 'SIZE'
         include 'TSTEP'

         integer i
         real pfrac ! time past the most recent period boundary
         real ftledt ! sub-step width for integration
         real Lrc(lpert, lpert) ! interpolated L_r at period boundary
         real fact, period
         integer, save :: icalld
         data icalld/0/
         real, save :: t0, Lrp(lpert, lpert)
         data t0/0.0d0/
         logical, save :: init
         data init/.false./
         character(len=20) fmte

      !  --- Determine FTLE horizon ---
         if (otd_FTLEPeriod > 0.0d0) then
            period = otd_FTLEPeriod
            pfrac = mod(time, period)
         else
            if (icalld == 0) then
               t0 = time
               icalld = 1
            end if
            period = time - t0
            pfrac = time - t0
         end if

      !  --- First call: initialise ---
         if (.not. init) then
            call copy(Lrp, otd_Lr, lpert*lpert)
            init = .true.
            if (nid == 0) then
               open (unit=459, file='otd_ftle.dat',
     $              status='replace', form='formatted')
               close (459)
            end if
         end if

      !  --- Integrate L_r(i,i) with trapezoidal rule ---
         if (pfrac < dt) then
      !     We just crossed a period boundary.  Split the step:
      !       [t_prev .. t_boundary]  then  [t_boundary .. t_curr]
      !     Linearly interpolate L_r at the boundary.
            ftledt = dt - pfrac
            call copy(Lrc, Lrp, lpert*lpert)
            fact = ftledt/dt
            call add2s2(Lrc, Lrp, -fact, lpert*lpert)
            call add2s2(Lrc, otd_Lr, fact, lpert*lpert)

      !     Finish previous period
            call integrate_trap(Lrp, Lrc, period, ftledt)

            if (nid == 0) then
               write (6, *) '[OTD] FTLE PRD', istep,
     $              't=', time, (FTLEv(i), i=1, npert)
            end if

      !     Reset and start new period
            call otd_zero_FTLE
            call copy(Lrp, Lrc, lpert*lpert)
            call copy(Lrc, otd_Lr, lpert*lpert)
            call integrate_trap(Lrp, Lrc, period, pfrac)

         else
      !     Normal step (no boundary crossing)
            call integrate_trap(Lrp, otd_Lr, pfrac, dt)

         end if

      !  --- Shift history ---
         call copy(Lrp, otd_Lr, lpert*lpert)

      !  --- Output (skip istep=0 where period=0 produces zeros) ---
         if (nid == 0 .and. istep > 0) then
            open (unit=459, file='otd_ftle.dat',
     $           position='append', status='unknown',
     $           form='formatted')
            write (fmte, '("(",I0,"(E15.7,1X))")') npert + 1
            write (459, fmte) time, (FTLEv(i), i=1, npert)
            close (459)
         end if

         if (nid == 0 .and.
     $       mod(istep, otd_printStep) == 0) then
            write (6, *) '[OTD] FTLE', istep, 't=', time,
     $           pfrac, (FTLEv(i), i=1, npert)
         end if

      !  --- Convergence check (at printStep frequency) ---
         if (istep > 0 .and.
     $       mod(istep, otd_printStep) == 0
     $       .and. otd_convTol > 0.0d0) then
            call otd_check_FTLE_convergence
         end if

      contains

         subroutine integrate_trap(Lrprev, Lrcurr,
     $                             deltat, hstep)
      !     Trapezoidal quadrature: accumulate LEintegral and
      !     update FTLEv = LEintegral / deltat.
            implicit none
            include 'SIZE'
            include 'TSTEP'

            real, intent(in) :: deltat ! FTLE averaging window
            real, intent(in) :: hstep ! integration sub-step
            real, intent(in), dimension(lpert, lpert) :: Lrprev
            real, intent(in), dimension(lpert, lpert) :: Lrcurr
            integer :: i

            do i = 1, npert
               LEintegral(i) = LEintegral(i)
     $              + 0.50d0*hstep*(Lrprev(i, i) + Lrcurr(i, i))
               if (deltat > 0.0d0) then
                  FTLEv(i) = LEintegral(i)/deltat
               else
                  FTLEv(i) = 0.0d0
               end if
            end do

         end subroutine integrate_trap

         subroutine otd_check_FTLE_convergence
      !     Compute step-to-step FTLE residuals, write to file,
      !     and stop early if all modes have converged.
            implicit none
            include 'SIZE'
            include 'TSTEP'

            real :: resid(lpert), rmax
            integer :: i
            character(len=20) :: fmtr

      !     Compute absolute change per mode
            rmax = 0.0d0
            do i = 1, npert
               resid(i) = abs(FTLEv(i) - FTLEv_prev(i))
               if (resid(i) /= resid(i)) resid(i) = 0.0d0
               if (resid(i) > rmax) rmax = resid(i)
            end do

      !     Write residual to file
            if (nid == 0) then
               open (unit=460, file='otd_residuals.dat',
     $              position='append', status='unknown',
     $              form='formatted')
               write (fmtr,
     $              '("(",I0,"(E15.7,1X))")')
     $              npert + 1
               write (460, fmtr)
     $              time, (resid(i), i=1, npert)
               close (460)
            end if

      !     Screen output
            if (nid == 0) then
               write (6, '(A,I7,1x,E14.7,A,E10.3)')
     $              '  [OTD] FTLE resid ',
     $              istep, time,
     $              '  max|dFTLE|= ', rmax
            end if

      !     Check convergence
            if (istep > otd_minSteps
     $           .and. rmax < otd_convTol) then
               if (nid == 0) then
                  write (6, *) ' '
                  write (6, *) '=============================='
     $                 //'========================'
                  write (6, '(A,E10.3,A,I7)')
     $                 '  [OTD] FTLEs converged!'
     $                 //' max|dFTLE|= ',
     $                 rmax, '  at step ', istep
                  write (6, '(A,E10.3)')
     $                 '  [OTD] Tolerance = ',
     $                 otd_convTol
                  write (fmtr,
     $                 '("(",I0,"(E15.7,1X))")')
     $                 npert
                  write (6, '(A)', ADVANCE='NO')
     $                 '  [OTD] Final FTLEs: '
                  write (6, fmtr)
     $                 (FTLEv(i), i=1, npert)
                  write (6, *) '=============================='
     $                 //'========================'
               end if
               lastep = 1
            end if

      !     Update history
            call copy(FTLEv_prev, FTLEv, lpert)

         end subroutine otd_check_FTLE_convergence

      end subroutine otd_compute_FTLE

      !-----------------------------------------------------------------------
      ! otd_construct_convective_terms — Lu_conv = (u.grad) Ub + (Ub.grad) u
      !-----------------------------------------------------------------------
      subroutine otd_construct_convective_terms(uxp, uyp, uzp, ipert)
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d
         include 'SOLN' ! v[xyz]

         real, dimension(lx1*ly1*lz1*lelv), intent(in) :: uxp, uyp, uzp
         real, dimension(lx1, ly1, lz1, lelv), save :: ta1, ta2, ta3, tb1, tb2, tb3
         integer, intent(in) :: ipert
         integer nv

         nv = lx1*ly1*lz1*lelv

         if (if3d) then
            call opcopy(tb1, tb2, tb3, vx, vy, vz) ! Save velocity
            call opcopy(vx, vy, vz, uxp, uyp, uzp) ! U <-- u
      !     convop(conv,fld): builds the convective term for the scalar field fld
      !     conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
            call convop(ta1, tb1) ! (u.grad) Ub (takes care of dealiasing)
            call convop(ta2, tb2)
            call convop(ta3, tb3)
      !     Copy fields into the correct variables
            call opcopy(convx(1, ipert), convy(1, ipert), convz(1, ipert), ta1, ta2, ta3)
            call opcopy(vx, vy, vz, tb1, tb2, tb3) ! Restore velocity
      !     conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
            call convop(tb1, uxp) ! (Ub.grad) u
            call convop(tb2, uyp)
            call convop(tb3, uzp)
      !     Add fields to the convective term
            call opadd2(convx(1, ipert), convy(1, ipert), convz(1, ipert), tb1, tb2, tb3)
         else ! 2D
            call opcopy(tb1, tb2, tb3, vx, vy, vz) ! Save velocity
            call opcopy(vx, vy, vz, uxp, uyp, uzp) ! U <-- u
      !     convop(conv,fld): builds the convective term for the scalar field fld
      !     conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
            call convop(ta1, tb1) ! (u.grad) Ub
            call convop(ta2, tb2)
      !     Copy fields into the correct variables
            call opcopy(convx(1, ipert), convy(1, ipert), ta3, ta1, ta2, ta3)
            call opcopy(vx, vy, vz, tb1, tb2, tb3) ! Restore velocity
      !     conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
            call convop(tb1, uxp) ! (Ub.grad) u
            call convop(tb2, uyp)
      !     Add fields to the convective term
            call opadd2(convx(1, ipert), convy(1, ipert), tb3, tb1, tb2, tb3)
         end if ! if3d

      end subroutine otd_construct_convective_terms

      !-----------------------------------------------------------------------
      ! otd_construct_pressure_gradient_terms — pressure gradient of
      !   perturbation field
      !-----------------------------------------------------------------------
      subroutine otd_construct_pressure_gradient_terms(prpert, ipert)
         implicit none
         include 'SIZE'

         real, intent(in) :: prpert(lx2*ly2*lz2*lelv, 1) ! perturbation pressure field
         integer, intent(in) :: ipert ! number of the considered pert.
         real, dimension(lx1, ly1, lz1, lelv), save :: ta1, ta2, wrk

         call mappr(wrk, prpert, ta1, ta2) ! Map the perturbation pressure to the velocity mesh
         call gradm1(gradpx(1, ipert), gradpy(1, ipert), gradpz(1, ipert), wrk) ! gradient on the velocity mesh directly

      end subroutine otd_construct_pressure_gradient_terms

      !-----------------------------------------------------------------------
      ! otd_construct_diffusive_terms — 1/Re * grad^2 u
      !-----------------------------------------------------------------------
      subroutine otd_construct_diffusive_terms(uxp, uyp, uzp, ipert)
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d
         include 'SOLN' ! vdiff

         real, dimension(lx1*ly1*lz1*lelv), intent(in) :: uxp, uyp, uzp
         integer, intent(in) :: ipert
         integer nv

         nv = lx1*ly1*lz1*nelv
         call otd_compute_laplacian(diffx(1, ipert), uxp)
         call otd_compute_laplacian(diffy(1, ipert), uyp)
         if (if3d) call otd_compute_laplacian(diffz(1, ipert), uzp)
      ! multiply by 1/Re > remove for operator diagnostics
         call col2(diffx(1, ipert), vdiff, nv)
         call col2(diffy(1, ipert), vdiff, nv)
         if (if3d) call col2(diffz(1, ipert), vdiff, nv)

      end subroutine otd_construct_diffusive_terms

      !-----------------------------------------------------------------------
      ! otd_compute_laplacian — construct diffusion term for one
      !   velocity component
      !-----------------------------------------------------------------------
      subroutine otd_compute_laplacian(lapu, up)
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d
         include 'DXYZ' ! dxm1,d[xy]tm1
         include 'GEOM' ! r[xy]m1,s[xy]m1,t[xy]m1,jacmi

         real up(lx1*ly1*lz1*lelv, 1) ! perturbation velocity component
         real, dimension(lx1*ly1*lz1, lelv) :: lapu
         real, dimension(lx1*ly1*lz1, lelv), save :: ux, uy, uz
         real, dimension(lx1*ly1*lz1) :: otd_ur, otd_us, otd_ut
      ! common/ctmp1/otd_ur, otd_us, otd_ut
         integer e, i, lxyz, nel

         lxyz = lx1*ly1*lz1
         nel = nx1 - 1
         call gradm1(ux, uy, uz, up)
         do e = 1, nelv
            if (if3d) then
               call local_grad3(otd_ur, otd_us, otd_ut, ux, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  lapu(i, e) = jacmi(i, e)*(otd_ur(i)*rxm1(i, 1, 1, e) +
     $   otd_us(i)*sxm1(i, 1, 1, e) + otd_ut(i)*txm1(i, 1, 1, e))
               end do
               call local_grad3(otd_ur, otd_us, otd_ut, uy, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  lapu(i, e) = lapu(i, e) + jacmi(i, e)*(otd_ur(i)*rym1(i, 1, 1, e) +
     $   otd_us(i)*sym1(i, 1, 1, e) + otd_ut(i)*tym1(i, 1, 1, e))
               end do
               call local_grad3(otd_ur, otd_us, otd_ut, uz, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  lapu(i, e) = lapu(i, e) + jacmi(i, e)*(otd_ur(i)*rzm1(i, 1, 1, e) +
     $   otd_us(i)*szm1(i, 1, 1, e) + otd_ut(i)*tzm1(i, 1, 1, e))
               end do
            else ! 2D
               call local_grad2(otd_ur, otd_us, ux, nel, e, dxm1, dytm1)
               do i = 1, lxyz
                  lapu(i, e) = jacmi(i, e)*(otd_ur(i)*rxm1(i, 1, 1, e) + otd_us(i)*sxm1(i, 1, 1, e))
               end do
               call local_grad2(otd_ur, otd_us, uy, nel, e, dxm1, dytm1)
               do i = 1, lxyz
                  lapu(i, e) = lapu(i, e) + jacmi(i, e)*(otd_ur(i)*rym1(i, 1, 1, e) + otd_us(i)*sym1(i, 1, 1, e))
               end do
            end if ! if3d
         end do

      end subroutine otd_compute_laplacian

      !-----------------------------------------------------------------------
      ! otd_op_glsc2 — weighted inner product of velocity field with
      !   perturbation j
      !-----------------------------------------------------------------------
      real function otd_op_glsc2(vcx, vcy, vcz, jpert)
         implicit none
         include 'SIZE'
         include 'SOLN' ! v[xyz]p, jp
         include 'TSTEP' ! ifield
         include 'MASS' ! bm1

         real, dimension(lx1*ly1*lz1*lelv), intent(in) :: vcx, vcy, vcz
         integer, intent(in) :: jpert
         real :: op_glsc2_wt

         ifield = 1
         otd_op_glsc2 = 0.50d0*op_glsc2_wt(vcx, vcy, vcz, VXP(1, jpert), VYP(1, jpert), VZP(1, jpert), bm1)

      end function otd_op_glsc2

      !-----------------------------------------------------------------------
      ! otd_inner_product — inner product using velocity or operator
      !   fields (iflag=1: velocity, iflag=2: L_NS)
      !-----------------------------------------------------------------------
      real function otd_inner_product(ipert, jpert, iflag)
         implicit none
         include 'SIZE'
         include 'SOLN'
         include 'PARALLEL'

         integer, intent(in) :: ipert, jpert
         integer, intent(in) :: iflag

         otd_inner_product = 0.0d0

         if (iflag == 1) then
            otd_inner_product = otd_op_glsc2(VXP(1, ipert), VYP(1, ipert), VZP(1, ipert), jpert)
         elseif (iflag == 2) then
            otd_inner_product = otd_op_glsc2(otd_Lux(1, ipert), otd_Luy(1, ipert), otd_Luz(1, ipert), jpert)
         else
            if (nid == 0) write (6, *) 'Error: Invalid iflag in otd_inner_product'
            call exitt
         end if

      end function otd_inner_product

      !-----------------------------------------------------------------------
      ! otd_normalize_vector_field — u_i = v_i / ||v_i||
      !-----------------------------------------------------------------------
      subroutine otd_normalize_vector_field(uxp, uyp, uzp)
         implicit none
         include 'SIZE'
         include 'TSTEP' ! ifield
         include 'MASS' ! bm1
         include 'INPUT' ! if3d

         real, dimension(lx1*ly1*lz1*lelv, 1), intent(in) :: uxp, uyp, uzp
         integer :: nv
         real :: invnorm, n2, op_glsc2_wt

         ifield = 1
         nv = lx1*ly1*lz1*lelv
         n2 = 0.50d0*op_glsc2_wt(uxp, uyp, uzp, uxp, uyp, uzp, bm1)
         if (n2 <= 0.0d0) then
            if (nid == 0) write (6, *) 'Error in otd_normalize_vector_field!'; call nek_end
         end if
         invnorm = 1.0d0/sqrt(n2)
         call cmult(uxp(1:nv, :), invnorm, nv)
         call cmult(uyp(1:nv, :), invnorm, nv)
         if (if3d) call cmult(uzp(1:nv, :), invnorm, nv)

      end subroutine otd_normalize_vector_field

      !-----------------------------------------------------------------------
      ! otd_mod_Gram_Schmidt — Modified Gram-Schmidt orthonormalization
      !   on the perturbation velocity field
      !-----------------------------------------------------------------------
      subroutine otd_mod_Gram_Schmidt
      !Perform Modified Gram-Schmidt orthonormalization on the
      !         perturbation velocity field for improved numerical stability
      !
      !     do i=1,npert
      !       u_i = v_i/||v_i||
      !       do j=i+1,npert
      !         u_j = v_j - proj_{u_i} (v_j)
      !       enddo
      !     enddo
      !
      !        with proj_{u_i} (v_j) = < v_j , u_i >/||u_i|| * u_i
      !                              = < v_j , u_i > * u_i   since ||u_i|| = 1
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d
         include 'TSTEP' ! istep
         include 'SOLN' ! V[XYZ]P

         integer i, j, nv
         real invnorm, proj

         nv = lx1*ly1*lz1*nelv
         do i = 1, npert ! orthonormalize
            invnorm = 1/sqrt(otd_inner_product(i, i, 1))
            call cmult(vxp(1, i), invnorm, nv)
            call cmult(vyp(1, i), invnorm, nv)
            if (if3d) call cmult(vzp(1, i), invnorm, nv)
            do j = i + 1, npert
               proj = -otd_inner_product(i, j, 1)
               call add2s2(vxp(1, j), vxp(1, i), proj, nv)
               call add2s2(vyp(1, j), vyp(1, i), proj, nv)
               if (if3d) call add2s2(vzp(1, j), vzp(1, i), proj, nv)
            end do
         end do

         if (nid == 0) write (6, *) 'V[XYZ]P orthonormalized (otd_mod_Gram_Schmidt)', istep

      end subroutine otd_mod_Gram_Schmidt

      !-----------------------------------------------------------------------
      ! otd_gram_matrix — batch Gram matrix via BLAS dgemm + single gop
      !
      ! G(i,j) = 0.5 * sum_k( ax(k,i)*bx(k,j)
      !                      + ay(k,i)*by(k,j)
      !                      + az(k,i)*bz(k,j) )   (global)
      !
      ! ax/ay/az must be pre-weighted by the mass matrix bm1.
      ! bx/by/bz are unweighted perturbation fields.
      ! npts = nx1*ny1*nz1*nelv (active points), ldv = lx1*ly1*lz1*lelv
      !-----------------------------------------------------------------------
      subroutine otd_gram_matrix(G, ax, ay, az,
     $     bx, by, bz, npts, ldv)
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3d

         integer, intent(in) :: npts, ldv
         real, intent(out) :: G(lpert, lpert)
         real, intent(in), dimension(ldv, lpert) ::
     $        ax, ay, az, bx, by, bz
         real :: wk_gop(lpert*lpert)

      !     X-component (beta=0 initialises G)
         call dgemm('T', 'N', npert, npert, npts,
     $        0.5d0, ax, ldv, bx, ldv,
     $        0.0d0, G, lpert)
      !     Y-component (accumulate)
         call dgemm('T', 'N', npert, npert, npts,
     $        0.5d0, ay, ldv, by, ldv,
     $        1.0d0, G, lpert)
      !     Z-component (3D only)
         if (if3d) then
            call dgemm('T', 'N', npert, npert, npts,
     $           0.5d0, az, ldv, bz, ldv,
     $           1.0d0, G, lpert)
         end if

      !     Single global reduction
         call gop(G, wk_gop, '+  ', lpert*lpert)

      end subroutine otd_gram_matrix

      !-----------------------------------------------------------------------
      ! otd_compute_orthonormality_measures — normality and orthogonality
      !   diagnostics for the perturbation basis
      !-----------------------------------------------------------------------
      subroutine otd_compute_orthonormality_measures(
     $     normality, orthogonality, info, flag)
         implicit none
         include 'SIZE'
         include 'INPUT'
         include 'MASS'
         include 'SOLN'
         include 'TSTEP'

         real, intent(out) :: normality, orthogonality
         logical, intent(in) :: flag
         character(len=6), intent(in) :: info
         integer :: pert_i, pert_j, nv, ldv, j
         real :: G(lpert, lpert)

         nv = lx1*ly1*lz1*nelv
         ldv = lx1*ly1*lz1*lelv

      !     Weight vxp/vyp/vzp into diffx/diffy/diffz workspace
         do j = 1, npert
            call copy(diffx(1, j), vxp(1, j), nv)
            call col2(diffx(1, j), bm1, nv)
            call copy(diffy(1, j), vyp(1, j), nv)
            call col2(diffy(1, j), bm1, nv)
         end do
         if (if3d) then
            do j = 1, npert
               call copy(diffz(1, j), vzp(1, j), nv)
               call col2(diffz(1, j), bm1, nv)
            end do
         end if

      !     Gram matrix via shared helper
         call otd_gram_matrix(G, diffx, diffy, diffz,
     $        vxp, vyp, vzp, nv, ldv)

         if (nid == 0) then
            normality = 0.0d0
            do pert_i = 1, npert
               normality = normality
     $              + G(pert_i, pert_i)**2
            end do
            normality = sqrt(normality/npert)

            if (npert > 1) then
               orthogonality = 0.0d0
               do pert_i = 1, npert
                  do pert_j = pert_i + 1, npert
                     orthogonality = orthogonality
     $                    + G(pert_i, pert_j)**2
                  end do
               end do
               orthogonality =
     $              sqrt(2.0d0*orthogonality)
     $              /(npert*(npert - 1))
            end if

            if (flag) then
               write (6, *) '  [OTD] NOout',
     $              istep, info,
     $              normality - 1.0d0, orthogonality
            end if
         end if

      end subroutine otd_compute_orthonormality_measures

      !-----------------------------------------------------------------------
      ! otd_sort_eigenvalues — sort eigenvalues by decreasing real part
      !   and reorder eigenvector columns accordingly
      !-----------------------------------------------------------------------
      subroutine otd_sort_eigenvalues(str)
      ! 1. Sort the eigenvalues l_i such that their real parts are ranked in decreasing order
      !       Re(l_1) .ge. Re(l_i) .ge. Re(l_r), i = 1,...,r
      !
      !  2. Apply the same sorting to the columns of the right
      !     eigenvector matrix and separate real and imaginary parts
      !       EVR => EVRR + i*EVRI
      !
         implicit none
         include 'SIZE'
         include 'TSTEP'

         character(len=*), intent(in) :: str
         integer i, j, id
         logical mk(lpert)
         real, dimension(lpert) :: wrk1, wrk2

         if (nid == 0) write (6, *) 'OTD: Sorting eigs of ', str

      !   zero out indices, output and mask
         call izero(otd_idx, lpert)
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

      ! sorting
         j = 1
         do while (j <= npert)
            EIGR(j) = maxval(wrk1, mask=mk) ! find largest real eigenvalue in remaining list
            id = maxloc(wrk1, 1, mk) ! find its index
            EIGI(j) = wrk2(id) ! extract corresponding imaginary part
            mk(id) = .false. ! update mask
            otd_idx(j) = id
            if (abs(EIGI(j)) < 1.0d-12) then
               do i = 1, npert
                  EVRR(i, j) = EVR(i, id)
               end do
               j = j + 1
            else ! complex conjugate eigenvectors!
               EIGI(j + 1) = -EIGI(j)
               EIGR(j + 1) = EIGR(j)
               mk(id + 1) = .false.
               otd_idx(j + 1) = id + 1
               do i = 1, npert
                  EVRR(i, j) = EVR(i, id)
                  EVRR(i, j + 1) = EVR(i, id)
                  EVRI(i, j) = EVR(i, id + 1)
                  EVRI(i, j + 1) = -EVR(i, id + 1)
               end do
               j = j + 2
            end if
         end do

      end subroutine otd_sort_eigenvalues

      !-----------------------------------------------------------------------
      ! otd_compute_eig_wrapper — LAPACK dgeev interface for
      !   non-symmetric eigenvalue problem on the reduced operator
      !-----------------------------------------------------------------------
      subroutine otd_compute_eig_wrapper(n, kind)
      ! LAPACK interface for the non-symmetric eigenvalue solver.
      ! Upon finishing, RITZR and RITZI contain the real and imaginary parts
      ! of the computed eigenvalues. Complex conjugate pairs of the
      ! eigenvalues appear with the eigenvalue having the positive
      ! imaginary part first.
      ! The corresponding eigenvectors are stored in EVEC. If the j:th and
      ! (j+1):th eigenvalue form a complex conjugate pair, then:
      ! v(j) = EVEC(:,j)+i*EVEC(:,j+1), v(j+1) = EVEC(:,j)-i*EVEC(:,j+1)
         implicit none
         include 'SIZE'

         character(len=1) :: jobvl, jobvr, kind
         integer :: n, lda, ldvl, ldvr, info, i, i0

      ! Input parameter 'kind' determines whether left and/or right eigenvectors should be computed
         if (kind == 'r') then
            jobvl = 'N'; jobvr = 'V'
         elseif (kind == 'l') then
            jobvl = 'V'; jobvr = 'N'
         elseif (kind == 'b') then
            jobvl = 'V'; jobvr = 'V'
         else
            if (nid == 0) write (6, *) 'ERROR: choose left/right/both eigenvectors', kind
         end if

         lda = npert; ldvl = npert; ldvr = npert

      ! Compute the eigenvalues/eigenvectors in double precision
         call dgeev(jobvl, jobvr, n, otd_Lr, lda, EIGR, EIGI, EVL, ldvl, EVR, ldvr, RWORK, LWORKR, info)

      ! Error-check
         if (info < 0) then
            if (nid == 0) write (6, *) 'ERROR: the i:th argument had an illegal value.', abs(info)
            call exitt
         elseif (info > 0) then
            if (nid == 0) then
               write (6, *) 'ERROR: the QR algorithm failed.', info
               write (6, *) '         Converged eigenvalues:'
               i0 = info + 1
               do i = i0, n
                  write (6, *) EIGR(i), EIGI(i)
               end do
            end if
            call exitt
      ! else
      !    if (nid == 0) then
      !       write (6, *) 'DGEEV: successful exit!'
      !       write (6, *) '        Optimal LWORKR=', int(RWORK(1)), LWORKR
      !    end if ! nid.eq.0
         end if ! info < 0

      end subroutine otd_compute_eig_wrapper

      end module nekstab_otd
