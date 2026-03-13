      !-----------------------------------------------------------------------
      ! energy_budget.f90 — Perturbation kinetic energy budget analysis
      !
      ! Purpose:
      !   Computes spatial and integrated PKE budget terms (production,
      !   dissipation) for stability analysis. Includes steady-state and
      !   Floquet (time-periodic) formulations.
      !
      ! Public interface:
      !   stability_energy_budget         — PKE budget for steady modes
      !   stability_energy_budget_floquet — PKE budget for Floquet modes
      !   compute_velocity_gradient_tensor — 9-component gradient tensor
      !   compute_dissipation              — viscous dissipation term
      !   compute_production               — Reynolds stress production
      !   compute_gradients                — smooth gradient via dsavg
      !   compute_laplacian                — Laplacian via double gradient
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL, ADJOINT
      !-----------------------------------------------------------------------

   module nekstab_energy_budget
   use nekstab_vectors
   use nekstab_nek_bridge
   use nekstab_io
   use nekstab_matvec
   use nekstab_eigensolvers
   implicit none
   private
   public :: stability_energy_budget, &
      stability_energy_budget_floquet, &
      compute_velocity_gradient_tensor, &
      compute_dissipation, compute_production, &
      compute_gradients, compute_laplacian
    logical, save :: energy_init = .false.

   contains

      !-----------------------------------------------------------------------
      ! stability_energy_budget — PKE budget for steady base flows
      !
      ! Purpose:
      !   Loads converged eigenmodes (dRe*, dIm*), normalizes to unit
      !   norm, computes spatial production and dissipation fields, and
      !   integrates to verify growth rate from the energy equation.
      !
      ! Notes:
      !   - Loops over maxmodes converged eigenvalues
      !   - Outputs spatial fields (K01, K02, K03) and integrated budget
      !   - Production has 9 components (3 gradients × 3 velocity comps)
      !   - Dissipation is a single scalar field
      !-----------------------------------------------------------------------
    subroutine stability_energy_budget
     use krylov_subspace

        !     ----- Arrays to store the instability mode real and imaginary parts.
    !     Allocatable (heap) to avoid stack overflow on large 3D meshes;
    !     each real(lv) array is ~(lx1^3 * nelv * 8) bytes.
    real, allocatable :: vx_dRe(:), vy_dRe(:), vz_dRe(:), t_dRe(:)
    real, allocatable :: vx_dIm(:), vy_dIm(:), vz_dIm(:), t_dIm(:)
    real, allocatable :: pr_dRe(:), pr_dIm(:)

       !     ----- Energy budget terms.
    real, allocatable :: energy_budget(:,:)
    real, allocatable :: integrals(:)

      !     ----- Miscellaneous.
   real :: alpha, beta, glsc2
   integer :: i, k, mode, n
   character(len=80) :: filename
   character(len=6) :: mode_str
   character(len=3) :: mode_str2

   n = nx1*ny1*nz1*nelv

      !     --> Allocate local arrays.
   allocate(vx_dRe(lv), vy_dRe(lv), vz_dRe(lv), t_dRe(lt))
   allocate(vx_dIm(lv), vy_dIm(lv), vz_dIm(lv), t_dIm(lt))
   allocate(pr_dRe(lp), pr_dIm(lp))
   allocate(energy_budget(lv, 10), integrals(10))

      !     #####
      !     #####
      !     #####     PREPROCESSING
      !     #####
      !     #####

      !     --> Load the base flow.
   write (filename, '(a, a, a)') 'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)
   call nopcopy(ubase, vbase, wbase, pbase, tbase, vx, vy, vz, pr, t)

   do mode = 1, maxmodes
   energy_budget(:, :) = 0.0d0
   integrals(:) = 0.0d0

      ! Format the mode number with leading zeros
   write (mode_str, '(i5.5)') mode

      !      --> Load the real part of the mode.
   write (filename, '(a, a, a, a)') 'dRe', trim(SESSION), '0.f', trim(mode_str)
   call load_fld(filename)
   call nopcopy(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, vx, vy, vz, pr, t)

      !     --> Load the imaginary part of the mode.
   write (filename, '(a, a, a, a)') 'dIm', trim(SESSION), '0.f', trim(mode_str)
   call load_fld(filename)
   call nopcopy(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, vx, vy, vz, pr, t)

      !     --> Normalize eigenmode to unit-norm.
      !         norm() returns ||q||, nopcmult(q, c) does q = c*q.
      !         Old bug: alpha = sqrt(alpha**2 + beta**2)
      !         gave ||q_new|| = alpha*||q|| = alpha^2 (wrong).
      !         Fix: invert so nopcmult divides by ||q_total||.
   call norm(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, alpha)
   call norm(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, beta)

   alpha = 1.0d0 / sqrt(alpha**2 + beta**2)

   call nopcmult(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, alpha)
   call nopcmult(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, alpha)

      !     #####
      !     #####
      !     #####     COMPUTE THE ENERGY BUDGET
      !     #####
      !     #####

      !     --> Compute the production terms.
   call compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, 1, &
      energy_budget(:, 1), energy_budget(:, 2), energy_budget(:, 3))
   call compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, 2, &
      energy_budget(:, 4), energy_budget(:, 5), energy_budget(:, 6))
   call compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, 3, &
      energy_budget(:, 7), energy_budget(:, 8), energy_budget(:, 9))

      !     --> Compute the dissipation term.
   call compute_dissipation(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, energy_budget(:, 10))

      !     --> Compute the integrals and the sum.
   do i = 1, 10
   integrals(i) = glsc2(bm1, energy_budget(:, i), n)
   end do

   if (if3d) then
   k = 9
   else
   k = 6
   end if

   do i = 1, k, 3
   call opcopy(vx, vy, vz, energy_budget(:, i), energy_budget(:, i + 1), energy_budget(:, i + 2))

   write (mode_str2, "('K',I2.2)") mode
   call outpost(vx, vy, vz, pr, t, mode_str2)
   end do

   if (nid == 0) then
   write (filename, '(A,A,A,A)') 'PKE_dRe', trim(SESSION), '0.f', trim(mode_str)
   open (101, file=filename, form='formatted')
   do i = 1, 10
   write (101, '(1E15.7)') integrals(i)
   write (*, *) 'Integral ', i, ' = ', integrals(i)
   end do
   write (101, '(1E15.7)') sum(integrals(1:9))
   write (101, '(1E15.7)') sum(integrals(1:9)) - integrals(10)
   close (101)
   end if ! nid == 0

   end do

   deallocate(vx_dRe, vy_dRe, vz_dRe, t_dRe)
   deallocate(vx_dIm, vy_dIm, vz_dIm, t_dIm)
   deallocate(pr_dRe, pr_dIm)
   deallocate(energy_budget, integrals)

   end subroutine stability_energy_budget

      !-----------------------------------------------------------------------
      ! stability_energy_budget_floquet — Orbit-averaged PKE budget
      !
      ! Purpose:
      !   Computes orbit-integrated PKE budget for Floquet modes on
      !   time-periodic base flows. Corrects for exponential growth,
      !   verifies energy equation over one period.
      !
      ! Mathematical formulation:
      !   For Floquet mode u'(x,t) = exp(sigma_r * t) * u_hat(x,t),
      !   orbit-averaging the Reynolds-Orr equation gives:
      !     2 * sigma_r = (P_bar - D_bar) / E_bar
      !   where bars denote orbit averages and growth correction is
      !   applied via weight = exp(-2*sigma_r*t) * dt / T.
      !
      ! Notes:
      !   - Real multipliers only (Mode A, Mode B)
      !   - Base flow orbit stored on first mode, reused for subsequent
      !   - Factor-of-1/2 convention in compute_* helpers (see docstring)
      !   - Right-endpoint rectangle rule for orbit integral (first-order)
      !-----------------------------------------------------------------------
    subroutine stability_energy_budget_floquet
     use krylov_subspace

        !     ----- Perturbation arrays -----
    real, allocatable :: vx_d(:), vy_d(:), vz_d(:)
    real, allocatable :: vx_zero(:), vy_zero(:), vz_zero(:)

       !     ----- Budget accumulators (allocatable to avoid stack overflow) --
    real, allocatable :: budget_avg(:,:)
    real, allocatable :: energy_avg(:)
    real, dimension(10) :: integrals
    real, allocatable :: prod_x(:), prod_y(:), prod_z(:), diss_tmp(:)

      !     ----- Eigenvalue data -----
   real :: sigma_r, omega, period, growth_corr, weight
   real :: E_bar, sigma_check, rel_error

      !     ----- Miscellaneous -----
   real :: glsc2
   integer :: i, j, k, n, mode, m, col
   character(len=80) :: filename
   character(len=6) :: mode_str
   character(len=3) :: mode_str2

   logical ifverbose
   logical :: mode_reloaded

   n = nx1*ny1*nz1*nelv
   nt = nx1*ny1*nz1*nelt

      !     --> Allocate local arrays.
   allocate(vx_d(lv), vy_d(lv), vz_d(lv))
   allocate(vx_zero(lv), vy_zero(lv), vz_zero(lv))
   allocate(prod_x(lv), prod_y(lv), prod_z(lv), diss_tmp(lv))

      !     #####
      !     #####     READ EIGENVALUE
      !     #####

   call read_eigenvalue(sigma_r, omega)
   if (nid == 0) then
   write (6, *) 'Floquet PKE budget: sigma_r =', sigma_r
   write (6, *) 'Floquet PKE budget: omega   =', omega
   end if

      !     #####
      !     #####     LOAD BASE FLOW AND PREPARE SOLVER
      !     #####

   write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)

      !     --> Set up linearized solver (computes nsteps, dt from param(10))
   call prepare_linearized_solver

   period = nsteps * dt
   if (nid == 0) write (6, *) &
      'Orbit period T =', period, ' nsteps =', nsteps

      !     --> Setup linearized solver flags
   ifpert = .true.; ifadj = .false.
   call bcast(ifpert, lsize); call bcast(ifadj, lsize)

      !     --> Zero arrays for imaginary part (real multiplier only)
   call rzero(vx_zero, lv)
   call rzero(vy_zero, lv)
   call rzero(vz_zero, lv)

      !     --> Allocate orbit storage on first call
   if (.not. energy_init) then
   call allocate_orbit(nsteps)
   end if

      !     --> Allocate budget accumulators on heap
   allocate(budget_avg(lv, 10))
   allocate(energy_avg(lv))

      !     #####
      !     #####     LOOP OVER MODES
      !     #####

   do mode = 1, maxmodes

   budget_avg = 0.0d0
   energy_avg = 0.0d0
   integrals  = 0.0d0

      !        --> Format mode number
   write (mode_str, '(i5.5)') mode

      !        --> Load eigenmode (real part only for real multiplier)
   write (filename, '(a,a,a,a)') &
      'dRe', trim(SESSION), '0.f', trim(mode_str)
   call load_fld(filename)

      !        --> Pass eigenmode as perturbation IC
   call opcopy(vxp(:,1), vyp(:,1), vzp(:,1), &
      vx, vy, vz)
   if (ifto) call copy(tp(:,:,1), t, nt)

      !        --> Reload base flow IC into vx,vy,vz for orbit evolution
   if (.not. energy_init) then
   write (filename, '(a,a,a)') &
      'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)
   ifbase = .true.
   else
   ifbase = .false.
      !           --> Restore first orbit step from stored orbit
   call orbit_restore(1)
   end if

      !        ─────────────────────────────────────────
      !        ORBIT INTEGRATION + BUDGET ACCUMULATION
      !        ─────────────────────────────────────────

   time = 0.0d0
   do istep = 1, nsteps

      !           --> Log progress
   if (nid == 0) write (6, &
      "(' PKE_FLOQUET mode',I3,':',I6,'/',I6)") &
      mode, istep, nsteps

      !           --> Advance (BF + perturbation simultaneously)
   call nekstab_usrchk()
   call nek_advance()

      !           --> Store/load orbit
   if (.not. energy_init) then
   call orbit_store(istep)
   else
   call orbit_restore(istep)
   end if

      !           --> Copy current base flow to ubase/vbase/wbase
      !               (compute_production reads from ubase,vbase,wbase)
   call opcopy(ubase, vbase, wbase, vx, vy, vz)

      !           --> Get current perturbation
   call opcopy(vx_d, vy_d, vz_d, &
      vxp(:,1), vyp(:,1), vzp(:,1))

      !           --> Growth correction + trapezoidal weight
   growth_corr = exp(-2.0d0 * sigma_r * time)
   weight = growth_corr * dt / period

      !           --> Production terms (9 = 3 components x 3 gradients)
   do j = 1, 3
   call compute_production(vx_d, vy_d, vz_d, &
      vx_zero, vy_zero, vz_zero, j, &
      prod_x, prod_y, prod_z)
   col = (j - 1) * 3
   do i = 1, n
   budget_avg(i,col+1) = budget_avg(i,col+1) &
      + weight * prod_x(i)
   budget_avg(i,col+2) = budget_avg(i,col+2) &
      + weight * prod_y(i)
   budget_avg(i,col+3) = budget_avg(i,col+3) &
      + weight * prod_z(i)
   end do
   end do

      !           --> Dissipation term (1 scalar field)
   call compute_dissipation(vx_d, vy_d, vz_d, &
      vx_zero, vy_zero, vz_zero, diss_tmp)
   do i = 1, n
   budget_avg(i,10) = budget_avg(i,10) &
      + weight * diss_tmp(i)
   end do

      !           --> Perturbation kinetic energy
   do i = 1, n
   energy_avg(i) = energy_avg(i) + weight * 0.5d0 &
      * (vx_d(i)**2 + vy_d(i)**2 + vz_d(i)**2)
   end do

   end do ! istep

      !        --> Mark orbit as stored after first mode
   if (.not. energy_init) then
   ifbase = .false.
   energy_init = .true.
   end if

      !        ─────────────────────────────────────────
      !        INTEGRATE AND VERIFY
      !        ─────────────────────────────────────────

   do i = 1, 10
   integrals(i) = glsc2(bm1, budget_avg(:,i), n)
   end do
   E_bar = glsc2(bm1, energy_avg, n)

   sigma_check = (sum(integrals(1:9)) - integrals(10)) &
      / (2.0d0 * E_bar)

   rel_error = abs(sigma_check - sigma_r) &
      / max(abs(sigma_r), 1.0d-30)

   if (nid == 0) then
   write (6, *) ''
   write (6, *) &
      '=== Orbit-averaged PKE budget (mode', mode, &
      ') ==='
   write (6, '(A,E15.7)') &
      '  sigma_r (eigenvalue) = ', sigma_r
   write (6, '(A,E15.7)') &
      '  sigma_r (budget)     = ', sigma_check
   write (6, '(A,E15.7)') &
      '  relative error       = ', rel_error
   write (6, '(A,E15.7)') &
      '  E_bar                = ', E_bar
   write (6, '(A,E15.7)') &
      '  P_bar (total)        = ', &
      sum(integrals(1:9))
   write (6, '(A,E15.7)') &
      '  D_bar                = ', integrals(10)
   do i = 1, 10
   write (6, '(A,I2,A,E15.7)') &
      '  Integral ', i, ' = ', integrals(i)
   end do
   write (6, *) ''
   end if

      !        ─────────────────────────────────────────
      !        OUTPUT
      !        ─────────────────────────────────────────

      !        --> Output spatial budget fields (same format as steady)
   if (if3d) then
   k = 9
   else
   k = 6
   end if

   do i = 1, k, 3
   call opcopy(vx, vy, vz, &
      budget_avg(:,i), budget_avg(:,i+1), &
      budget_avg(:,i+2))
   write (mode_str2, "('F',I2.2)") mode
   call outpost(vx, vy, vz, pr, t, mode_str2)
   end do

      !        --> Output scalar integrals to file
   if (nid == 0) then
   write (filename, '(A,A,A,A)') &
      'PKE_floquet_', trim(SESSION), '0.f', &
      trim(mode_str)
   open (101, file=filename, form='formatted')
   do i = 1, 10
   write (101, '(1E15.7)') integrals(i)
   end do
   write (101, '(1E15.7)') sum(integrals(1:9))
   write (101, '(1E15.7)') &
      sum(integrals(1:9)) - integrals(10)
   write (101, '(1E15.7)') E_bar
   write (101, '(1E15.7)') sigma_check
   close (101)
   end if

   end do ! mode

      !     --> Cleanup
   if (allocated(budget_avg)) deallocate(budget_avg)
   if (allocated(energy_avg)) deallocate(energy_avg)
   if (allocated(uor)) deallocate(uor, vor, wor)
   if (allocated(tor)) deallocate(tor)

   end subroutine stability_energy_budget_floquet

      !-----------------------------------------------------------------------
      ! compute_velocity_gradient_tensor — 9-component gradient tensor
      !
      ! Purpose:
      !   Computes all 9 components of velocity gradient tensor using
      !   gradm1 and smooths at element interfaces via dsavg.
      !
      ! Arguments:
      !   vx_in, vy_in, vz_in [in]  — input velocity field
      !   dudx..dwdz [out]          — 9 gradient components
      !-----------------------------------------------------------------------
   subroutine compute_velocity_gradient_tensor( &
      vx_in, vy_in, vz_in, &
      dudx, dudy, dudz, &
      dvdx, dvdy, dvdz, &
      dwdx, dwdy, dwdz)

    real, dimension(lx1*ly1*lz1*lelv), intent(in) :: vx_in, vy_in, vz_in
   real, dimension(lx1*ly1*lz1*lelv), intent(out) :: dudx, dudy, dudz
   real, dimension(lx1*ly1*lz1*lelv), intent(out) :: dvdx, dvdy, dvdz
   real, dimension(lx1*ly1*lz1*lelv), intent(out) :: dwdx, dwdy, dwdz

   call gradm1(dudx, dudy, dudz, vx_in, nelv)
   call gradm1(dvdx, dvdy, dvdz, vy_in, nelv)
   call gradm1(dwdx, dwdy, dwdz, vz_in, nelv)

   call dsavg(dudx); call dsavg(dudy); call dsavg(dudz)
   call dsavg(dvdx); call dsavg(dvdy); call dsavg(dvdz)
   call dsavg(dwdx); call dsavg(dwdy); call dsavg(dwdz)

   end subroutine compute_velocity_gradient_tensor

      !-----------------------------------------------------------------------
      ! compute_dissipation — Viscous dissipation term
      !
      ! Purpose:
      !   Computes 0.5 * (mu/rho) * u'_i * lap(u'_i) for complex mode
      !   (Re^2 + Im^2 formulation).
      !
      ! Arguments:
      !   vx_dRe, vy_dRe, vz_dRe [in] — real part of mode
      !   vx_dIm, vy_dIm, vz_dIm [in] — imaginary part of mode
      !   dissipation [out]           — dissipation field
      !
      ! Notes:
      !   Factor of 0.5 is convention for complex-mode formulation.
      !-----------------------------------------------------------------------
   subroutine compute_dissipation(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, dissipation)
    use krylov_subspace

    real, dimension(lv), intent(in) :: vx_dRe, vy_dRe, vz_dRe
   real, dimension(lv), intent(in) :: vx_dIm, vy_dIm, vz_dIm

   real, allocatable, dimension(:) :: Laplacian_ax, Laplacian_ay, Laplacian_az
   real, allocatable, dimension(:) :: Laplacian_bx, Laplacian_by, Laplacian_bz

   real, dimension(lv), intent(out) :: dissipation
   real, allocatable, dimension(:) :: dummy

   allocate(Laplacian_ax(lv), Laplacian_ay(lv), Laplacian_az(lv))
   allocate(Laplacian_bx(lv), Laplacian_by(lv), Laplacian_bz(lv))
   allocate(dummy(lv))

      !     --> Compute Laplacians.
   call compute_laplacian(vx_dRe, Laplacian_ax)
   call compute_laplacian(vy_dRe, Laplacian_ay)
   call compute_laplacian(vz_dRe, Laplacian_az)

   call compute_laplacian(vx_dIm, Laplacian_bx)
   call compute_laplacian(vy_dIm, Laplacian_by)
   call compute_laplacian(vz_dIm, Laplacian_bz)

      !     --> Compute dissipation term.
   dissipation = 0.0d+00

   dummy = vx_dRe*Laplacian_ax + vx_dIm*Laplacian_bx
   dissipation = dissipation + dummy

   dummy = vy_dRe*Laplacian_ay + vy_dIm*Laplacian_by
   dissipation = dissipation + dummy

   dummy = vz_dRe*Laplacian_az + vz_dIm*Laplacian_bz
   dissipation = dissipation + dummy

   dissipation = 0.5*dissipation*param(2)/param(1)

   end subroutine compute_dissipation

      !-----------------------------------------------------------------------
      ! compute_production — Reynolds stress production term
      !
      ! Purpose:
      !   Computes -0.5 * u'_i u'_j dU_i/dx_j for one gradient direction.
      !
      ! Arguments:
      !   vx_dRe, vy_dRe, vz_dRe [in] — real part of mode
      !   vx_dIm, vy_dIm, vz_dIm [in] — imaginary part of mode
      !   component [in]              — which base flow gradient (1,2,3)
      !   prod_x, prod_y, prod_z [out] — 3 production components
      !
      ! Notes:
      !   Factor of -0.5 is convention for complex-mode formulation.
      !   Reads base flow from ubase, vbase, wbase globals.
      !-----------------------------------------------------------------------
   subroutine compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, component, prod_x, prod_y, prod_z)
    use krylov_subspace

    real, dimension(lv), intent(in) :: vx_dRe, vy_dRe, vz_dRe
   real, dimension(lv), intent(in) :: vx_dIm, vy_dIm, vz_dIm
   real, dimension(lv), intent(out) :: prod_x, prod_y, prod_z
   real, allocatable, dimension(:) :: dcdx, dcdy, dcdz
   integer, intent(in) :: component

   allocate(dcdx(lv), dcdy(lv), dcdz(lv))

   if (component == 1) then
   call gradm1(dcdx, dcdy, dcdz, ubase, nelv)

   prod_x = -0.5*(vx_dRe**2 + vx_dIm**2)*dcdx
   prod_y = -0.5*(vx_dRe*vy_dRe + vy_dIm*vx_dIm)*dcdy
   prod_z = -0.5*(vx_dRe*vz_dRe + vz_dIm*vx_dIm)*dcdz

   else if (component == 2) then
   call gradm1(dcdx, dcdy, dcdz, vbase, nelv)

   prod_x = -0.5*(vx_dRe*vy_dRe + vy_dIm*vx_dIm)*dcdx
   prod_y = -0.5*(vy_dRe**2 + vy_dIm**2)*dcdy
   prod_z = -0.5*(vy_dRe*vz_dRe + vz_dIm*vy_dIm)*dcdz

   else if (component == 3) then
   call gradm1(dcdx, dcdy, dcdz, wbase, nelv)

   prod_x = -0.5*(vx_dRe*vz_dRe + vz_dIm*vx_dIm)*dcdx
   prod_y = -0.5*(vy_dRe*vz_dRe + vz_dIm*vy_dIm)*dcdy
   prod_z = -0.5*(vz_dRe**2 + vz_dIm**2)*dcdz
   end if

   end subroutine compute_production

      !-----------------------------------------------------------------------
      ! compute_gradients — Smoothed gradient via dsavg
      !
      ! Purpose:
      !   Computes gradient of scalar field and smooths at element
      !   interfaces using direct stiffness averaging.
      !
      ! Arguments:
      !   u [in]            — input scalar field
      !   dudx, dudy, dudz [out] — gradient components
      !-----------------------------------------------------------------------
   subroutine compute_gradients(u, dudx, dudy, dudz)
    use krylov_subspace

    real, dimension(lv), intent(in) :: u
   real, dimension(lv), intent(out) :: dudx, dudy, dudz

   call gradm1(dudx, dudy, dudz, u)
   call dsavg(dudx); call dsavg(dudy); call dsavg(dudz)

   end subroutine compute_gradients

      !-----------------------------------------------------------------------
      ! compute_laplacian — Laplacian via double gradient
      !
      ! Purpose:
      !   Computes Laplacian by taking gradient twice (d²u/dx² + ...).
      !
      ! Arguments:
      !   a [in]     — input scalar field
      !   Lap_a [out] — Laplacian field
      !-----------------------------------------------------------------------
   subroutine compute_laplacian(a, Lap_a)
    use krylov_subspace

    real, dimension(lv), intent(in) :: a
   real, allocatable, save, dimension(:) :: dadx, dady, dadz
   real, allocatable, save, dimension(:) :: d2adx2, d2ady2, d2adz2
   real, dimension(lv), intent(out) :: Lap_a
   real, allocatable, save, dimension(:) :: wrk1, wrk2

   if (.not. allocated(dadx)) then
      allocate(dadx(lv), dady(lv), dadz(lv))
      allocate(d2adx2(lv), d2ady2(lv), d2adz2(lv))
      allocate(wrk1(lv), wrk2(lv))
   end if
   call compute_gradients(a, dadx, dady, dadz)
   call compute_gradients(dadx, d2adx2, wrk1, wrk2)
   call compute_gradients(dady, wrk1, d2ady2, wrk2)
   call compute_gradients(dadz, wrk1, wrk2, d2adz2)

   Lap_a = d2adx2 + d2ady2 + d2adz2

   end subroutine compute_laplacian
      !-----------------------------------------------------------------------

   end module nekstab_energy_budget
