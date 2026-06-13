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
!   nekstab_nek_bridge, krylov_subspace, nekstab_io, nekstab_noise
!-----------------------------------------------------------------------

module nekstab_otd
   use nekstab_nek_bridge
   use krylov_subspace
   use nekstab_io
   use nekstab_noise
   implicit none
   private

   ! ── Module-level persistent variables (replacing DATA statements) ──
   !  NOTE: OTDfx, OTDfy, OTDfz, otd_Lr, etc. are NOT declared here.
   !  They live in the OTD_modes / OTD_params common blocks in NEKSTAB.inc,
   !  accessible via `use nekstab_nek_bridge`.  Declaring them here as
   !  module variables would create a duplicate-declaration compiler error.
   logical :: OTD_init = .false.
   integer :: OTD_ftle_icalld = 0
   real :: OTD_ftle_t0 = 0.0d0
   logical :: OTD_ftle_init = .false.

   ! ── Public subroutines ──
   public :: otd, otd_construct_linear_operator, &
      otd_compute_OTD_modes, &
      otd_outpost_OTD_modes, &
      otd_outpost_orthonormal_basis, &
      otd_white_noise, otd_generate_forces, &
      otd_orthonormalize_basis, &
      otd_zero_FTLE, otd_compute_FTLE, &
      otd_construct_convective_terms, &
      otd_construct_pressure_gradient_terms, &
      otd_construct_diffusive_terms, &
      otd_compute_laplacian, &
      otd_op_glsc2, otd_inner_product, &
      otd_normalize_vector_field, &
      otd_mod_Gram_Schmidt, &
      otd_compute_orthonormality_measures, &
      otd_gram_matrix, &
      otd_sort_eigenvalues, &
      otd_compute_eig_wrapper
contains

!-----------------------------------------------------------------------
! otd — main OTD driver (initialization and time-stepping)
!
! Implementation based on Babaee & Sapsis (2016)
! https://dx.doi.org/10.1098/rspa.2015.0779
!-----------------------------------------------------------------------
subroutine otd
   logical :: file_found
   integer :: mode, restart_id, latest_restart
   character(len=3) :: mode_tag
   character(len=30) :: filename

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.1.
   if (.not. OTD_init) then

      ifpert = .true.
      call bcast(ifpert, lsize)
      param(31) = lpert
      npert = int(param(31))

      ifotd = .true.
      call bcast(ifotd, lsize)

      if (nid == 0) then
         open (unit=457, file='otd_growth_rates.dat', status='replace', form='formatted'); close (457)
         open (unit=458, file='otd_eigenvalues.dat', status='replace', form='formatted'); close (458)
         open (unit=460, file='otd_residuals.dat', status='replace', form='formatted'); close (460)
      end if

      write (filename, '(A,A,A)') 'BF_', trim(SESSION), '0.f00001'
      inquire(file=filename, exist=file_found)
      if (file_found) then
         if (nid == 0) write (6, *) 'OTD: loading base flow from ', trim(filename)
         call load_fld(filename)
         call opcopy(ubic, vbic, wbic, vx, vy, vz)
      else
         if (nid == 0) write (6, *) 'OTD: no BF_ file found, using useric as base flow'
         call opcopy(ubic, vbic, wbic, vx, vy, vz)
      end if
      call outpost(ubic, vbic, wbic, pr, t, 'bf0')

      call otd_white_noise ! Initialise all IC fields to white noise

      latest_restart = 0
      do restart_id = 1, 99
         write (filename, '(A,A,A,I5.5)') &
            'r01', trim(SESSION), '0.f', restart_id
         inquire (file=filename, exist=file_found)
         if (file_found) then
            latest_restart = restart_id
         end if
      end do

      if (latest_restart > 0) then
         if (nid == 0) write (6, *) 'Found ', latest_restart, ' IC files'
         do mode = 1, npert
            write (mode_tag, '(A1,I2.2)') 'r', mode
            write (filename, '(2A,I5.5)') &
               trim(mode_tag), trim(SESSION)//'0.f', latest_restart
            if (nid == 0) write (6, *) 'Looking for ICs in ', trim(filename)
            inquire (file=filename, exist=file_found)
            if (file_found) then
               call load_fld(filename)
               call opcopy(vxpic(1, mode), &
                  vypic(1, mode), vzpic(1, mode), vx, vy, vz)
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
      do mode = 1, npert
         write (mode_tag, '(I1)') mode
         call outpost(upic(1, mode), vpic(1, mode), &
            wpic(1, mode), pr, t, 'ip'//trim(mode_tag))
      end do

      if (uparam(1) > 5) then
         ifbase = .true.
         call bcast(ifbase, lsize)
      elseif (uparam(1) == 5) then
         if (nid == 0) write (6, *) 'OTD in Frozen baseflow mode!'
         ifbase = .false.
         call bcast(ifbase, lsize)
      end if

      gsstep_override = .true.
      call otd_orthonormalize_basis
      call otd_zero_FTLE ! zero Phi matrix

      OTD_init = .true.
   end if ! init

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.1.
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
   integer :: mode, row, active_count, storage_count

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.2.
   ! Build the linearized Navier-Stokes action on each OTD basis vector:
   !   L_NS(u_j) = (1/Re) grad^2 u_j - grad p_j
   !             - (Ub.grad) u_j - (u_j.grad) Ub.
   do mode = 1, npert
      call otd_construct_convective_terms( &
         vxp(1, mode), vyp(1, mode), vzp(1, mode), mode)
      call otd_construct_pressure_gradient_terms(prp(1, mode), mode)
      call otd_construct_diffusive_terms( &
         vxp(1, mode), vyp(1, mode), vzp(1, mode), mode)
   end do

   active_count = lx1*ly1*lz1*nelv
   storage_count = lx1*ly1*lz1*lelv

   do mode = 1, npert
      do row = 1, active_count
         otd_Lux(row, mode) = &
            diffx(row, mode) - gradpx(row, mode) - convx(row, mode)
         otd_Luy(row, mode) = &
            diffy(row, mode) - gradpy(row, mode) - convy(row, mode)
         if (if3d) then
            otd_Luz(row, mode) = &
               diffz(row, mode) - gradpz(row, mode) - convz(row, mode)
         end if
      end do
   end do

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.2-§1.3.
   ! Apply the mass weight before forming <L_NS(u_i), u_j>.
   do mode = 1, npert
      call col2(otd_Lux(1, mode), bm1, active_count)
      call col2(otd_Luy(1, mode), bm1, active_count)
      if (if3d) call col2(otd_Luz(1, mode), bm1, active_count)
   end do

   call otd_gram_matrix(otd_Lr, &
      otd_Lux, otd_Luy, otd_Luz, &
      vxp, vyp, vzp, active_count, storage_count)

   ! Internal rotation is skew-symmetric:
   !   phi_rot(i,j) = -phi_rot(j,i).
   ! It is extracted from the just-computed reduced operator Lr.
   call rzero(phi_rot, lpert*lpert)
   if (npert > 1) then
      do mode = 1, npert
         do row = mode + 1, npert
            phi_rot(row, mode) =  otd_Lr(row, mode)
            phi_rot(mode, row) = -otd_Lr(row, mode)
         end do
      end do
   end if
   call sub2(otd_Lr, phi_rot, lpert*lpert)

end subroutine otd_construct_linear_operator

!-----------------------------------------------------------------------
! otd_compute_OTD_modes — eigenspectrum of reduced operator and
!   projection onto eigendirections
!-----------------------------------------------------------------------
subroutine otd_compute_OTD_modes

   real :: saved_Lr(lpert, lpert)
   integer :: active_count, row, col
   character(len=20) :: real_fmt, int_fmt

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.5.
   ! Symmetric growth-rate operator:
   !   Lsym = (Lr + Lr^T)/2.
   ! The full Lr is saved because the nonlinear OTD forcing still uses it.
   call copy(saved_Lr, otd_Lr, lpert*lpert)
   do row = 1, npert
      do col = 1, npert
         otd_Lr(row, col) = &
            0.50d0*(saved_Lr(row, col) + saved_Lr(col, row))
      end do
   end do

   call otd_compute_eig_wrapper(npert, 'r')
   call otd_sort_eigenvalues('otd_Ls')

   if (nid == 0) then ! save to file
      open (unit=457, file='otd_growth_rates.dat', position='append', status='unknown', form='formatted')
      write (real_fmt, '("(",I0,"(E15.7,1X))")') npert + 1
      write (457, real_fmt) time, (EIGR(row), row=1, npert)
      close (457)
   end if

   if (mod(istep, otd_printStep) == 0 .and. nid == 0) then ! print out
      write (real_fmt, '("(",I0,"(E15.7,1X))")') npert
      if (nid == 0) write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Ls | Re '
      write (6, real_fmt) (EIGR(row), row=1, npert)
   end if

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.5.
   ! Restore the unsymmetrized reduced operator and eigendecompose it to
   ! obtain the OTD-mode directions.
   call copy(otd_Lr, saved_Lr, lpert*lpert)
   call otd_compute_eig_wrapper(npert, 'r')
   call otd_sort_eigenvalues('otd_Lr ')
   call copy(otd_Lr, saved_Lr, lpert*lpert)

   if (nid == 0) then ! save to file
      open (unit=458, file='otd_eigenvalues.dat', position='append', status='unknown', form='formatted')
      write (real_fmt, '("(",I0,"(E15.7,1X))")') npert + 1
      write (458, real_fmt) time, (EIGR(row), row=1, npert)
      close (458)
   end if

   if (mod(istep, otd_printStep) == 0 .and. nid == 0) then ! print out
      write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Lr | Re '
      write (6, real_fmt) (EIGR(row), row=1, npert)

      write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Lr | Im '
      write (6, real_fmt) (EIGI(row), row=1, npert) ! print order for reference
      write (int_fmt, '("(",I0,"(I4,1X))")') npert

      write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 's-otd_idx   '
      write (6, int_fmt) (otd_idx(row), row=1, npert) ! print out non-zero elements of rotated otd_Lr

      write (real_fmt, '("(",I0,"(E15.7,1X))")') npert*(npert + 1)/2
      write (6, '(A,I7,1x,E14.7,1x,A8)', ADVANCE='NO') '  [OTD] ', istep, time, 'Lrmat   '
      write (6, real_fmt) ((otd_Lr(row, col), col=row, npert), row=1, npert)

   end if

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.5.
   ! Project the physical OTD basis fields onto the sorted eigendirections:
   ! real parts use EVRR, imaginary parts use EVRI.
   active_count = lx1*ly1*lz1*nelv
   call mxm(vxp, active_count, EVRr, lpert, OTDmrx, npert)
   call mxm(vxp, active_count, EVRi, lpert, OTDmix, npert)
   call mxm(vyp, active_count, EVRr, lpert, OTDmry, npert)
   call mxm(vyp, active_count, EVRi, lpert, OTDmiy, npert)
   if (if3d) then
      call mxm(vzp, active_count, EVRr, lpert, OTDmrz, npert)
      call mxm(vzp, active_count, EVRi, lpert, OTDmiz, npert)
   end if

end subroutine otd_compute_OTD_modes

!-----------------------------------------------------------------------
! otd_outpost_OTD_modes — output projection of orthonormal basis
!   on eigendirections
!-----------------------------------------------------------------------
subroutine otd_outpost_OTD_modes

   integer :: mode
   character(len=2) :: mode_tag
   character(len=3) :: field_prefix

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §3.3.
   ! docs/otd-spec.md §3.4: OTDm* fields are NEKSTAB.inc common-block
   ! storage; keep the public handoff and field-prefix layout unchanged.
   ! These M* snapshots are diagnostic OTD modes, not restart fields.
   call otd_compute_OTD_modes
   do mode = 1, npert
      if (lpert >= 10) then
         write (mode_tag, '(I2.2)') mode
         field_prefix = 'M'//trim(mode_tag)
      else
         write (mode_tag, '(I1)') mode
         field_prefix = 'M'//trim(mode_tag)//'_'
      end if
      call outpost(OTDmrx(1, mode), OTDmry(1, mode), &
         OTDmrz(1, mode), prp(1, mode), t, field_prefix)

   end do

end subroutine otd_outpost_OTD_modes

!-----------------------------------------------------------------------
! otd_outpost_orthonormal_basis — output the OTD basis directly
!   for restart
!-----------------------------------------------------------------------
subroutine otd_outpost_orthonormal_basis

   integer :: mode
   character(len=2) :: mode_tag
   character(len=3) :: field_prefix

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §3.3.
   ! Force the same orthonormalization path used by the time-step driver, then
   ! write restart-compatible rNN fields from the shared perturbation basis.
   ! A later OTD run discovers these rNN files during initialization and uses
   ! them as the perturbation restart basis after the white-noise fallback.
   gsstep_override = .true.
   call otd_orthonormalize_basis

   do mode = 1, npert
      write (mode_tag, '(I2.2)') mode
      field_prefix = 'r'//trim(mode_tag)
      call outpost(vxp(1, mode), vyp(1, mode), &
         vzp(1, mode), prp(1, mode), t, field_prefix)
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

   integer :: mode, ix, iy, iz, ie, ieg, nv, ijke
   real :: xl(ldim), fc(3), sin2
   real :: glmin, glmax, nmin, nmax

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.6.
   ! The mode-dependent sin(mode)^2 factor scales the frequency coefficients
   ! before each mth_rand call, so every perturbation mode gets a distinct
   ! deterministic hash field while keeping restart/IC common-block layout.
   ! Scaling before hashing matters: changing the frequency tuple changes the
   ! hash trajectory itself, instead of merely rescaling one shared noise field.
   nv = nx1*ny1*nz1*nelv

   do mode = 1, npert
      sin2 = sin(real(mode))**2
      call rzero(vxpic(1, mode), nv)
      call rzero(vypic(1, mode), nv)
      if (if3d) call rzero(vzpic(1, mode), nv)

      do ie = 1, nelv
         do iz = 1, nz1
            do iy = 1, ny1
               do ix = 1, nx1
                  xl(1) = xm1(ix, iy, iz, ie)
                  xl(2) = ym1(ix, iy, iz, ie)
                  if (if3d) xl(ndim) = zm1(ix, iy, iz, ie)
                  ijke = ix + lx1*((iy - 1) &
                     + ly1*((iz - 1) + lz1*(ie - 1)))
                  ieg = lglel(ie)

                  fc(1) = sin2*3.0e4
                  fc(2) = sin2*(-1.5e3)
                  fc(3) = sin2*0.5e5
                  vxpic(ijke, mode) = &
                     mth_rand(ix, iy, iz, ieg, xl, fc)

                  fc(1) = sin2*2.3e4
                  fc(2) = sin2*2.3e3
                  fc(3) = sin2*(-2.0e5)
                  vypic(ijke, mode) = &
                     mth_rand(ix, iy, iz, ieg, xl, fc)

                  if (if3d) then
                     fc(1) = sin2*2.0e4
                     fc(2) = sin2*1.0e3
                     fc(3) = sin2*1.0e5
                     vzpic(ijke, mode) = &
                        mth_rand(ix, iy, iz, ieg, xl, fc)
                  end if
               end do
            end do
         end do
      end do

      ! docs/otd-spec.md §2.6: enforce continuity, mass weighting, averaging,
      ! and velocity Dirichlet masks after all components are generated.
      call opdssum(vxpic(1, mode), vypic(1, mode), vzpic(1, mode))
      call opcolv(vxpic(1, mode), vypic(1, mode), &
         vzpic(1, mode), vmult)
      call dsavg(vxpic(1, mode))
      call dsavg(vypic(1, mode))
      if (if3d) call dsavg(vzpic(1, mode))
      call bcdirVC(vxpic(1, mode), vypic(1, mode), &
         vzpic(1, mode), v1mask, v2mask, v3mask)
   end do

   nmin = glmin(vxpic, nv)
   nmax = glmax(vxpic, nv)
   if (nid == 0) write (6, *) 'OTD noise min,max', nmin, nmax

end subroutine otd_white_noise

!-----------------------------------------------------------------------
! otd_generate_forces — create forcing for OTD evolution equation
!-----------------------------------------------------------------------
subroutine otd_generate_forces

   integer :: active_count

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.4.
   active_count = lx1*ly1*lz1*nelv
   call mxm(VXP, active_count, otd_Lr, lpert, OTDfx, npert)
   call mxm(VYP, active_count, otd_Lr, lpert, OTDfy, npert)
   if (if3d) call mxm(VZP, active_count, otd_Lr, lpert, OTDfz, npert)

end subroutine otd_generate_forces

!-----------------------------------------------------------------------
! otd_orthonormalize_basis — orthonormalize perturbation basis
!-----------------------------------------------------------------------
subroutine otd_orthonormalize_basis

   real :: N, O
   logical :: needs_gs

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.3.
   call otd_compute_orthonormality_measures(N, O, 'pre   ', .true.)

   needs_gs = .false.
   if (otd_gsStep /= 0) then ! gsstep=0 => no GS
      if (mod(istep, otd_gsStep) == 0) then
         needs_gs = .true.
      end if
   end if
   if (gsstep_override) needs_gs = .true. ! override when needed

   if (needs_gs) then
      call otd_mod_Gram_Schmidt ! Modified Gram-Schmidt
   end if

end subroutine otd_orthonormalize_basis

!-----------------------------------------------------------------------
! otd_zero_FTLE — reset FTLE computation
!-----------------------------------------------------------------------
subroutine otd_zero_FTLE

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.6.
   ! This clears the running integral behind
   !   lambda_i(T) = (1/T) integral_0^T Lr(i,i) dt.
   call rzero(FTLEv, lpert)
   call rzero(LEintegral, lpert)

end subroutine otd_zero_FTLE

!-----------------------------------------------------------------------
! otd_compute_FTLE — compute finite-time Lyapunov exponents
!-----------------------------------------------------------------------
subroutine otd_compute_FTLE
   integer :: mode
   real :: pfrac
   real :: ftledt
   real :: Lrc(lpert, lpert)
   real :: fact, period
   real, save :: Lrp(lpert, lpert)
   character(len=20) :: fmte

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §1.6.
   ! Accumulate finite-time Lyapunov exponents by trapezoidal integration of
   ! diag(Lr), using the exact period-boundary split-step order in the spec.
   ! The exposed value is lambda_i(T) = (1/T) integral_0^T Lr(i,i) dt,
   ! where T is either the configured FTLE period or elapsed time since reset.
   if (otd_FTLEPeriod > 0.0d0) then
      period = otd_FTLEPeriod
      pfrac = mod(time, period)
   else
      if (OTD_ftle_icalld == 0) then
         OTD_ftle_t0 = time
         OTD_ftle_icalld = 1
      end if
      period = time - OTD_ftle_t0
      pfrac = time - OTD_ftle_t0
   end if

   if (.not. OTD_ftle_init) then
      call copy(Lrp, otd_Lr, lpert*lpert)
      OTD_ftle_init = .true.
      if (nid == 0) then
         open (unit=459, file='otd_ftle.dat', &
            status='replace', form='formatted')
         close (459)
      end if
   end if

   if (pfrac < dt) then
      ! docs/otd-spec.md §1.6: split a boundary-crossing step and linearly
      ! interpolate Lr at the period boundary before resetting the integral.
      ! First integrate the tail of the old period, print/reset at the
      ! boundary, then integrate the head of the new period with the same
      ! trapezoid rule.  This preserves the boundary value in both halves.
      ftledt = dt - pfrac
      call copy(Lrc, Lrp, lpert*lpert)
      fact = ftledt/dt
      call add2s2(Lrc, Lrp, -fact, lpert*lpert)
      call add2s2(Lrc, otd_Lr, fact, lpert*lpert)

      call integrate_trap(Lrp, Lrc, period, ftledt)

      if (nid == 0) then
         write (6, *) '[OTD] FTLE PRD', istep, &
            't=', time, (FTLEv(mode), mode=1, npert)
      end if

      call otd_zero_FTLE
      call copy(Lrp, Lrc, lpert*lpert)
      call copy(Lrc, otd_Lr, lpert*lpert)
      call integrate_trap(Lrp, Lrc, period, pfrac)

   else
      call integrate_trap(Lrp, otd_Lr, pfrac, dt)

   end if

   call copy(Lrp, otd_Lr, lpert*lpert)

   ! docs/otd-spec.md §3.2: preserve otd_ftle.dat layout and precision.
   if (nid == 0 .and. istep > 0) then
      open (unit=459, file='otd_ftle.dat', &
         position='append', status='unknown', &
         form='formatted')
      write (fmte, '("(",I0,"(E15.7,1X))")') npert + 1
      write (459, fmte) time, (FTLEv(mode), mode=1, npert)
      close (459)
   end if

   if (nid == 0 .and. &
      mod(istep, otd_printStep) == 0) then
      write (6, *) '[OTD] FTLE', istep, 't=', time, &
         pfrac, (FTLEv(mode), mode=1, npert)
   end if

   if (istep > 0 .and. &
      mod(istep, otd_printStep) == 0 &
      .and. otd_convTol > 0.0d0) then
      call otd_check_FTLE_convergence
   end if

contains

   subroutine integrate_trap(Lrprev, Lrcurr, &
      deltat, hstep)
      real, intent(in) :: deltat
      real, intent(in) :: hstep
      real, intent(in), dimension(lpert, lpert) :: Lrprev
      real, intent(in), dimension(lpert, lpert) :: Lrcurr
      integer :: mode

      ! docs/otd-spec.md §1.6: use only diagonal entries and update FTLEv
      ! immediately from the current averaging horizon.
      ! Trapezoid update for lambda_i:
      !   integral_i += h/2 * (Lrprev(i,i) + Lrcurr(i,i)).
      do mode = 1, npert
         LEintegral(mode) = LEintegral(mode) &
            + 0.50d0*hstep*(Lrprev(mode, mode) + Lrcurr(mode, mode))
         if (deltat > 0.0d0) then
            FTLEv(mode) = LEintegral(mode)/deltat
         else
            FTLEv(mode) = 0.0d0
         end if
      end do

   end subroutine integrate_trap

   subroutine otd_check_FTLE_convergence
      real :: resid(lpert), rmax
      integer :: mode
      character(len=20) :: fmtr

      ! docs/otd-spec.md §3.4: residuals compare against the NEKSTAB.inc
      ! FTLEv_prev common-block handoff, then refresh it.
      ! The convergence diagnostic is max_i |lambda_i^n - lambda_i^{n-1}|.
      rmax = 0.0d0
      do mode = 1, npert
         resid(mode) = abs(FTLEv(mode) - FTLEv_prev(mode))
         if (resid(mode) /= resid(mode)) resid(mode) = 0.0d0
         if (resid(mode) > rmax) rmax = resid(mode)
      end do

      if (nid == 0) then
         open (unit=460, file='otd_residuals.dat', &
            position='append', status='unknown', &
            form='formatted')
         write (fmtr, &
            '("(",I0,"(E15.7,1X))")') &
            npert + 1
         write (460, fmtr) &
            time, (resid(mode), mode=1, npert)
         close (460)
      end if

      if (nid == 0) then
         write (6, '(A,I7,1x,E14.7,A,E10.3)') &
            '  [OTD] FTLE resid ', &
            istep, time, &
            '  max|dFTLE|= ', rmax
      end if

      !     Check convergence
      if (istep > otd_minSteps &
         .and. rmax < otd_convTol) then
         if (nid == 0) then
            write (6, *) ' '
            write (6, *) '==============================' &
               //'========================'
            write (6, '(A,E10.3,A,I7)') &
               '  [OTD] FTLEs converged!' &
               //' max|dFTLE|= ', &
               rmax, '  at step ', istep
            write (6, '(A,E10.3)') &
               '  [OTD] Tolerance = ', &
               otd_convTol
            write (fmtr, &
               '("(",I0,"(E15.7,1X))")') &
               npert
            write (6, '(A)', ADVANCE='NO') &
               '  [OTD] Final FTLEs: '
            write (6, fmtr) &
               (FTLEv(mode), mode=1, npert)
            write (6, *) '==============================' &
               //'========================'
         end if
         lastep = 1
      end if

      call copy(FTLEv_prev, FTLEv, lpert)

   end subroutine otd_check_FTLE_convergence

end subroutine otd_compute_FTLE

!-----------------------------------------------------------------------
! otd_construct_convective_terms — Lu_conv = (u.grad) Ub + (Ub.grad) u
!-----------------------------------------------------------------------
subroutine otd_construct_convective_terms(uxp, uyp, uzp, ipert)

   real, dimension(lx1*ly1*lz1*lelv), intent(in) :: uxp, uyp, uzp
   integer, intent(in) :: ipert
   real, dimension(lx1, ly1, lz1, lelv), save :: work_x, work_y, work_z
   real, dimension(lx1, ly1, lz1, lelv), save :: saved_x, saved_y, saved_z

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.5.
   ! convop differentiates its second argument along the current velocity
   ! field and handles Nek's dealiasing path internally.
   ! First pass: set the current velocity to u and compute (u.grad) Ub.
   ! Second pass: restore Ub as the current velocity and compute (Ub.grad) u.
   ! In 2D, the third component is only a dummy slot required by opcopy/opadd2.
   if (if3d) then
      call opcopy(saved_x, saved_y, saved_z, vx, vy, vz)
      call opcopy(vx, vy, vz, uxp, uyp, uzp)
      call convop(work_x, saved_x)
      call convop(work_y, saved_y)
      call convop(work_z, saved_z)
      call opcopy(convx(1, ipert), convy(1, ipert), &
         convz(1, ipert), work_x, work_y, work_z)

      call opcopy(vx, vy, vz, saved_x, saved_y, saved_z)
      call convop(saved_x, uxp)
      call convop(saved_y, uyp)
      call convop(saved_z, uzp)
      call opadd2(convx(1, ipert), convy(1, ipert), &
         convz(1, ipert), saved_x, saved_y, saved_z)
   else ! 2D
      call opcopy(saved_x, saved_y, saved_z, vx, vy, vz)
      call opcopy(vx, vy, vz, uxp, uyp, uzp)
      call convop(work_x, saved_x)
      call convop(work_y, saved_y)
      call opcopy(convx(1, ipert), convy(1, ipert), &
         work_z, work_x, work_y, work_z)

      call opcopy(vx, vy, vz, saved_x, saved_y, saved_z)
      call convop(saved_x, uxp)
      call convop(saved_y, uyp)
      call opadd2(convx(1, ipert), convy(1, ipert), &
         saved_z, saved_x, saved_y, saved_z)
   end if ! if3d

end subroutine otd_construct_convective_terms

!-----------------------------------------------------------------------
! otd_construct_pressure_gradient_terms — pressure gradient of
!   perturbation field
!-----------------------------------------------------------------------
subroutine otd_construct_pressure_gradient_terms(prpert, ipert)

   real, intent(in) :: prpert(lx2*ly2*lz2*lelv, 1)
   integer, intent(in) :: ipert
   real, dimension(lx1, ly1, lz1, lelv), save :: map_x, map_y, mapped_p

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.5.
   call mappr(mapped_p, prpert, map_x, map_y)
   call gradm1(gradpx(1, ipert), gradpy(1, ipert), &
      gradpz(1, ipert), mapped_p)

end subroutine otd_construct_pressure_gradient_terms

!-----------------------------------------------------------------------
! otd_construct_diffusive_terms — 1/Re * grad^2 u
!-----------------------------------------------------------------------
subroutine otd_construct_diffusive_terms(uxp, uyp, uzp, ipert)

   real, dimension(lx1*ly1*lz1*lelv), intent(in) :: uxp, uyp, uzp
   integer, intent(in) :: ipert
   integer :: active_count

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.5.
   ! otd_compute_laplacian returns grad^2 u; multiplying by vdiff applies
   ! the viscous coefficient, i.e. 1/Re for constant-viscosity cases.
   active_count = lx1*ly1*lz1*nelv
   call otd_compute_laplacian(diffx(1, ipert), uxp)
   call otd_compute_laplacian(diffy(1, ipert), uyp)
   if (if3d) call otd_compute_laplacian(diffz(1, ipert), uzp)
   call col2(diffx(1, ipert), vdiff, active_count)
   call col2(diffy(1, ipert), vdiff, active_count)
   if (if3d) call col2(diffz(1, ipert), vdiff, active_count)

end subroutine otd_construct_diffusive_terms

!-----------------------------------------------------------------------
! otd_compute_laplacian — construct diffusion term for one
!   velocity component
!-----------------------------------------------------------------------
subroutine otd_compute_laplacian(lapu, up)

   real, intent(in) :: up(lx1*ly1*lz1*lelv, 1)
   real, intent(out), dimension(lx1*ly1*lz1, lelv) :: lapu
   real, dimension(lx1*ly1*lz1, lelv), save :: grad_x, grad_y, grad_z
   real, dimension(lx1*ly1*lz1) :: ref_r, ref_s, ref_t
   integer :: elem, point, lxyz, polynomial_degree

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.5.
   lxyz = lx1*ly1*lz1
   polynomial_degree = nx1 - 1
   call gradm1(grad_x, grad_y, grad_z, up)
   do elem = 1, nelv
      if (if3d) then
         call local_grad3(ref_r, ref_s, ref_t, grad_x, &
            polynomial_degree, elem, dxm1, dxtm1)
         do point = 1, lxyz
            lapu(point, elem) = jacmi(point, elem) &
               *(ref_r(point)*rxm1(point, 1, 1, elem) &
               + ref_s(point)*sxm1(point, 1, 1, elem) &
               + ref_t(point)*txm1(point, 1, 1, elem))
         end do
         call local_grad3(ref_r, ref_s, ref_t, grad_y, &
            polynomial_degree, elem, dxm1, dxtm1)
         do point = 1, lxyz
            lapu(point, elem) = lapu(point, elem) &
               + jacmi(point, elem) &
               *(ref_r(point)*rym1(point, 1, 1, elem) &
               + ref_s(point)*sym1(point, 1, 1, elem) &
               + ref_t(point)*tym1(point, 1, 1, elem))
         end do
         call local_grad3(ref_r, ref_s, ref_t, grad_z, &
            polynomial_degree, elem, dxm1, dxtm1)
         do point = 1, lxyz
            lapu(point, elem) = lapu(point, elem) &
               + jacmi(point, elem) &
               *(ref_r(point)*rzm1(point, 1, 1, elem) &
               + ref_s(point)*szm1(point, 1, 1, elem) &
               + ref_t(point)*tzm1(point, 1, 1, elem))
         end do
      else ! 2D
         call local_grad2(ref_r, ref_s, grad_x, &
            polynomial_degree, elem, dxm1, dytm1)
         do point = 1, lxyz
            lapu(point, elem) = jacmi(point, elem) &
               *(ref_r(point)*rxm1(point, 1, 1, elem) &
               + ref_s(point)*sxm1(point, 1, 1, elem))
         end do
         call local_grad2(ref_r, ref_s, grad_y, &
            polynomial_degree, elem, dxm1, dytm1)
         do point = 1, lxyz
            lapu(point, elem) = lapu(point, elem) &
               + jacmi(point, elem) &
               *(ref_r(point)*rym1(point, 1, 1, elem) &
               + ref_s(point)*sym1(point, 1, 1, elem))
         end do
      end if ! if3d
   end do

end subroutine otd_compute_laplacian

!-----------------------------------------------------------------------
! otd_op_glsc2 — weighted inner product of velocity field with
!   perturbation j
!-----------------------------------------------------------------------
real function otd_op_glsc2(vcx, vcy, vcz, jpert)

   real, dimension(lx1*ly1*lz1*lelv), intent(in) :: vcx, vcy, vcz
   integer, intent(in) :: jpert
   real :: op_glsc2_wt

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.2.
   ifield = 1
   otd_op_glsc2 = 0.50d0*op_glsc2_wt( &
      vcx, vcy, vcz, VXP(1, jpert), VYP(1, jpert), &
      VZP(1, jpert), bm1)

end function otd_op_glsc2

!-----------------------------------------------------------------------
! otd_inner_product — inner product using velocity or operator
!   fields (iflag=1: velocity, iflag=2: L_NS)
!-----------------------------------------------------------------------
real function otd_inner_product(ipert, jpert, iflag)

   integer, intent(in) :: ipert, jpert
   integer, intent(in) :: iflag

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.2.
   otd_inner_product = 0.0d0

   if (iflag == 1) then
      otd_inner_product = otd_op_glsc2( &
         VXP(1, ipert), VYP(1, ipert), VZP(1, ipert), jpert)
   elseif (iflag == 2) then
      otd_inner_product = otd_op_glsc2( &
         otd_Lux(1, ipert), otd_Luy(1, ipert), &
         otd_Luz(1, ipert), jpert)
   else
      if (nid == 0) write (6, *) 'Error: Invalid iflag in otd_inner_product'
      call exitt
   end if

end function otd_inner_product

!-----------------------------------------------------------------------
! otd_normalize_vector_field — u_i = v_i / ||v_i||
!-----------------------------------------------------------------------
subroutine otd_normalize_vector_field(uxp, uyp, uzp)

   real, dimension(lx1*ly1*lz1*lelv, 1), intent(inout) :: uxp, uyp, uzp
   integer :: active_count
   real :: invnorm, norm_sq, op_glsc2_wt

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.3.
   ifield = 1
   active_count = lx1*ly1*lz1*nelv
   norm_sq = 0.50d0*op_glsc2_wt(uxp, uyp, uzp, uxp, uyp, uzp, bm1)
   if (norm_sq <= 0.0d0) then
      if (nid == 0) write (6, *) 'Error in otd_normalize_vector_field!'
      call nek_end
   end if
   invnorm = 1.0d0/sqrt(norm_sq)
   call cmult(uxp(1:active_count, :), invnorm, active_count)
   call cmult(uyp(1:active_count, :), invnorm, active_count)
   if (if3d) call cmult(uzp(1:active_count, :), invnorm, active_count)

end subroutine otd_normalize_vector_field

!-----------------------------------------------------------------------
! otd_mod_Gram_Schmidt — Modified Gram-Schmidt orthonormalization
!   on the perturbation velocity field
!-----------------------------------------------------------------------
subroutine otd_mod_Gram_Schmidt
   integer :: basis_col, target_col, active_count
   real :: invnorm, projection

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.3.
   ! Modified Gram-Schmidt, written in OTD basis notation:
   !   u_i = v_i / ||v_i||
   !   v_j <- v_j - <v_j, u_i> u_i,  j > i.
   ! The stored projection is negated because add2s2 performs an addition.
   active_count = lx1*ly1*lz1*nelv
   do basis_col = 1, npert
      invnorm = otd_inner_product(basis_col, basis_col, 1)
      if (invnorm > 0.0d0) then
         invnorm = 1.0d0/sqrt(invnorm)
      else
         invnorm = 0.0d0
      end if

      call cmult(vxp(1, basis_col), invnorm, active_count)
      call cmult(vyp(1, basis_col), invnorm, active_count)
      if (if3d) call cmult(vzp(1, basis_col), invnorm, active_count)

      do target_col = basis_col + 1, npert
         projection = -otd_inner_product(basis_col, target_col, 1)
         call add2s2(vxp(1, target_col), &
            vxp(1, basis_col), projection, active_count)
         call add2s2(vyp(1, target_col), &
            vyp(1, basis_col), projection, active_count)
         if (if3d) then
            call add2s2(vzp(1, target_col), &
               vzp(1, basis_col), projection, active_count)
         end if
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
! nv = nx1*ny1*nz1*nelv (active points), ldv = lx1*ly1*lz1*lelv
!-----------------------------------------------------------------------
subroutine otd_gram_matrix(G, ax, ay, az, &
   bx, by, bz, nv, ldv)

   integer, intent(in) :: nv, ldv
   real, intent(out) :: G(lpert, lpert)
   real, intent(in), dimension(ldv, lpert) :: &
      ax, ay, az, bx, by, bz
   real :: wk_gop(lpert*lpert)

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.1.
   ! Accumulate component contributions locally with BLAS:
   !   X initializes G with beta=0,
   !   Y accumulates into G with beta=1,
   !   Z accumulates likewise in 3D.
   ! One gop then performs the global reduction for the full matrix.
   call dgemm('T', 'N', npert, npert, nv, &
      0.5d0, ax, ldv, bx, ldv, &
      0.0d0, G, lpert)
   call dgemm('T', 'N', npert, npert, nv, &
      0.5d0, ay, ldv, by, ldv, &
      1.0d0, G, lpert)
   if (if3d) then
      call dgemm('T', 'N', npert, npert, nv, &
         0.5d0, az, ldv, bz, ldv, &
         1.0d0, G, lpert)
   end if

   call gop(G, wk_gop, '+  ', lpert*lpert)

end subroutine otd_gram_matrix

!-----------------------------------------------------------------------
! otd_compute_orthonormality_measures — normality and orthogonality
!   diagnostics for the perturbation basis
!-----------------------------------------------------------------------
subroutine otd_compute_orthonormality_measures( &
   normality, orthogonality, info, flag)

   real, intent(out) :: normality, orthogonality
   logical, intent(in) :: flag
   character(len=6), intent(in) :: info
   integer :: row, col, nv, ldv, mode
   real :: G(lpert, lpert)

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.4.
   ! Reuse diff* as mass-weighted workspace, then form the shared Gram matrix
   ! so normality and orthogonality use the same global reduction as the core.
   ! With G_ij = <u_i,u_j>, the diagnostics are
   !   normality = sqrt(sum_i G_ii^2 / npert)
   !   orthogonality = sqrt(2 sum_{i<j} G_ij^2) / (npert*(npert - 1)).
   nv = lx1*ly1*lz1*nelv
   ldv = lx1*ly1*lz1*lelv

   do mode = 1, npert
      call copy(diffx(1, mode), vxp(1, mode), nv)
      call col2(diffx(1, mode), bm1, nv)
      call copy(diffy(1, mode), vyp(1, mode), nv)
      call col2(diffy(1, mode), bm1, nv)
   end do
   if (if3d) then
      do mode = 1, npert
         call copy(diffz(1, mode), vzp(1, mode), nv)
         call col2(diffz(1, mode), bm1, nv)
      end do
   end if

   call otd_gram_matrix(G, diffx, diffy, diffz, &
      vxp, vyp, vzp, nv, ldv)

   if (nid == 0) then
      normality = 0.0d0
      do row = 1, npert
         normality = normality &
            + G(row, row)**2
      end do
      normality = sqrt(normality/npert)

      if (npert > 1) then
         orthogonality = 0.0d0
         do row = 1, npert
            do col = row + 1, npert
               orthogonality = orthogonality &
                  + G(row, col)**2
            end do
         end do
         orthogonality = &
            sqrt(2.0d0*orthogonality) &
            /(npert*(npert - 1))
      end if

      if (flag) then
         write (6, *) '  [OTD] NOout', &
            istep, info, &
            normality - 1.0d0, orthogonality
      end if
   end if

end subroutine otd_compute_orthonormality_measures

!-----------------------------------------------------------------------
! otd_sort_eigenvalues — sort eigenvalues by decreasing real part
!   and reorder eigenvector columns accordingly
!-----------------------------------------------------------------------
subroutine otd_sort_eigenvalues(str)
   character(len=*), intent(in) :: str
   integer :: row, sorted_col, source_col
   logical :: active(lpert)
   real, dimension(lpert) :: real_part, imag_part

   if (nid == 0) write (6, *) 'OTD: Sorting eigs of ', str

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.7.
   ! Sorting has two steps: rank eigenvalues by decreasing real part, then
   ! reorder the corresponding dgeev eigenvector columns into EVRR/EVRI.
   ! Real modes occupy one sorted column; complex pairs occupy two columns.
   call izero(otd_idx, lpert)
   call rzero(EVRR, lpert*lpert)
   call rzero(EVRI, lpert*lpert)
   do row = 1, npert
      active(row) = .true.
   end do
   call copy(real_part, EIGR, lpert)
   call copy(imag_part, EIGI, lpert)

   if (lpert > npert) then
      do row = npert + 1, lpert
         active(row) = .false.
      end do
   end if

   sorted_col = 1
   do while (sorted_col <= npert)
      EIGR(sorted_col) = maxval(real_part, mask=active)
      source_col = maxloc(real_part, 1, active)
      EIGI(sorted_col) = imag_part(source_col)
      active(source_col) = .false.
      otd_idx(sorted_col) = source_col

      if (abs(EIGI(sorted_col)) < 1.0d-12) then
         do row = 1, npert
            EVRR(row, sorted_col) = EVR(row, source_col)
         end do
         sorted_col = sorted_col + 1
      else ! complex conjugate eigenvectors!
         EIGI(sorted_col + 1) = -EIGI(sorted_col)
         EIGR(sorted_col + 1) = EIGR(sorted_col)
         active(source_col + 1) = .false.
         otd_idx(sorted_col + 1) = source_col + 1
         do row = 1, npert
            EVRR(row, sorted_col) = EVR(row, source_col)
            EVRR(row, sorted_col + 1) = EVR(row, source_col)
            EVRI(row, sorted_col) = EVR(row, source_col + 1)
            EVRI(row, sorted_col + 1) = -EVR(row, source_col + 1)
         end do
         sorted_col = sorted_col + 2
      end if
   end do

end subroutine otd_sort_eigenvalues

!-----------------------------------------------------------------------
! otd_compute_eig_wrapper — LAPACK dgeev interface for
!   non-symmetric eigenvalue problem on the reduced operator
!-----------------------------------------------------------------------
subroutine otd_compute_eig_wrapper(n, kind)
   integer, intent(in) :: n
   character(len=1) :: jobvl, jobvr
   character(len=1), intent(in) :: kind
   integer :: lda, ldvl, ldvr, info, eig_idx

   ! Implementation based on Babaee & Sapsis (2016); docs/otd-spec.md §2.7.
   ! dgeev stores complex-conjugate eigenpairs as adjacent real columns:
   ! for lambda = a +/- ib, EVR(:,k) is the real part and EVR(:,k+1)
   ! is the imaginary part. otd_sort_eigenvalues converts that contract
   ! into explicit EVRR/EVRI columns for later projection.
   if (kind == 'r') then
      jobvl = 'N'
      jobvr = 'V'
   elseif (kind == 'l') then
      jobvl = 'V'
      jobvr = 'N'
   elseif (kind == 'b') then
      jobvl = 'V'
      jobvr = 'V'
   else
      if (nid == 0) write (6, *) 'ERROR: choose left/right/both eigenvectors', kind
   end if

   lda = npert
   ldvl = npert
   ldvr = npert

   call dgeev(jobvl, jobvr, n, otd_Lr, lda, &
      EIGR, EIGI, EVL, ldvl, EVR, ldvr, RWORK, LWORKR, info)

   if (info < 0) then
      if (nid == 0) write (6, *) 'ERROR: the i:th argument had an illegal value.', abs(info)
      call exitt
   elseif (info > 0) then
      if (nid == 0) then
         write (6, *) 'ERROR: the QR algorithm failed.', info
         write (6, *) '         Converged eigenvalues:'
         do eig_idx = info + 1, n
            write (6, *) EIGR(eig_idx), EIGI(eig_idx)
         end do
      end if
      call exitt
   end if ! info < 0

end subroutine otd_compute_eig_wrapper

end module nekstab_otd
