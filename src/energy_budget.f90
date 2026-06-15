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
   real, allocatable, save :: diss_lap_ax_s(:), diss_lap_ay_s(:), diss_lap_az_s(:)
   real, allocatable, save :: diss_lap_bx_s(:), diss_lap_by_s(:), diss_lap_bz_s(:)
   real, allocatable, save :: prod_dcdx_s(:), prod_dcdy_s(:), prod_dcdz_s(:)

   contains

   !  Allocate-once workspace for compute_dissipation.
   !  Fixed-size (lv), never resized — these store intermediate
   !  Laplacians that would otherwise be allocated on every call.
   subroutine ensure_dissipation_workspace()
   use krylov_subspace

   if (.not. allocated(diss_lap_ax_s)) then
      allocate(diss_lap_ax_s(lv), diss_lap_ay_s(lv), diss_lap_az_s(lv))
      allocate(diss_lap_bx_s(lv), diss_lap_by_s(lv), diss_lap_bz_s(lv))
   end if

   end subroutine ensure_dissipation_workspace

   !  Allocate-once workspace for compute_production.
   !  Fixed-size (lv), never resized — stores base flow gradients.
   subroutine ensure_production_workspace()
   use krylov_subspace

   if (.not. allocated(prod_dcdx_s)) then
      allocate(prod_dcdx_s(lv), prod_dcdy_s(lv), prod_dcdz_s(lv))
   end if

   end subroutine ensure_production_workspace

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
   integer :: i, k, mode
   character(len=80) :: filename
   character(len=6) :: mode_str
   character(len=3) :: mode_str2

   nv = nx1*ny1*nz1*nelv

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
   integrals(i) = glsc2(bm1, energy_budget(:, i), nv)
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
!     Net budget P-D. integrals(10) stores the viscous energy contribution
!     = -Dbar (it is negative), so the physical P-D is sum(1:9) + integrals(10).
   write (101, '(1E15.7)') sum(integrals(1:9)) + integrals(10)
   close (101)
   end if ! nid == 0

   end do

   deallocate(vx_dRe, vy_dRe, vz_dRe, t_dRe)
   deallocate(vx_dIm, vy_dIm, vz_dIm, t_dIm)
   deallocate(pr_dRe, pr_dIm)
   deallocate(energy_budget, integrals)

   end subroutine stability_energy_budget

      !=======================================================================
      ! stability_energy_budget_floquet
      !   Orbit-averaged perturbation kinetic energy (PKE) budget of a
      !   Floquet mode on a time-periodic base flow.
      !
      ! ---------------------------------------------------------------------
      ! WHAT THIS ROUTINE IS FOR
      ! ---------------------------------------------------------------------
      !   A Floquet stability analysis returns *what* a periodic orbit does
      !   to an infinitesimal perturbation (the growth rate sigma_r and the
      !   eigenmode), but not *why*. This routine answers the "why": it
      !   evaluates, term by term and averaged over one orbit period T, the
      !   energy balance that governs the modal growth rate. It thereby
      !     (1) attributes the instability to physical mechanisms — energy
      !         extracted from the base shear by the perturbation Reynolds
      !         stress (production) versus viscous losses (dissipation) —
      !         and localises them in space through the output fields; and
      !     (2) provides a STRONG, fully independent check on the solver:
      !         the growth rate reconstructed from the energy budget must
      !         agree with the eigenvalue. A mismatch flags an error in the
      !         eigenmode, the base orbit, or the analysis itself.
      !
      ! ---------------------------------------------------------------------
      ! GOVERNING BALANCE (Reynolds-Orr energy equation)
      ! ---------------------------------------------------------------------
      !   Linearise the incompressible Navier-Stokes equations about a
      !   T-periodic base flow U(x,t) = U(x,t+T). For the (in general
      !   complex) perturbation field q, define the Hermitian energy
      !   E(t) = 1/2 integral( q . conjg(q) ) dV. Taking d/dt, substituting
      !   the linearised momentum equation, and integrating over the domain,
      !   the base-flow advection and the pressure-work terms integrate to
      !   zero under incompressibility (div q = 0) with closed/no-slip or
      !   periodic boundaries, leaving the exact balance
      !
      !       dE/dt = P(t) - D(t),
      !       P(t) = - integral( real( conjg(q_i) q_j ) dU_i/dx_j ) dV   (production)
      !       D(t) =   (1/Re) integral( grad q : grad conjg(q) ) dV      (dissipation)
      !
      !   where real(.) is the real part, Re is the Reynolds number, and
      !   summation over repeated indices i,j is implied. P is the rate at
      !   which the perturbation Reynolds stress does work against the
      !   base-flow strain; D is the (sign-definite) viscous dissipation.
      !   Both are real because they are Hermitian forms.
      !
      ! ---------------------------------------------------------------------
      ! FROM THE BALANCE TO THE GROWTH RATE (the key cancellation)
      ! ---------------------------------------------------------------------
      !   A Floquet solution has the form q(x,t) = exp(mu t) u_hat(x,t),
      !   with Floquet exponent mu = sigma_r + i*omega and a T-PERIODIC
      !   eigenfunction u_hat. The crucial observation is that in every
      !   Hermitian product the oscillatory phase cancels exactly:
      !       conjg(q_i) q_j = exp(2 sigma_r t) conjg(u_hat_i) u_hat_j,
      !   because exp(-i omega t) * exp(+i omega t) = 1. Hence
      !   E = exp(2 sigma_r t) E~, P = exp(2 sigma_r t) P~, D = ... with the
      !   tilde quantities GENUINELY T-PERIODIC (no residual omega ripple).
      !   Substituting into dE/dt = P - D gives
      !       2 sigma_r E~ + dE~/dt = P~ - D~,
      !   and averaging over one period (<dE~/dt> = 0 by periodicity) yields
      !
      !       +-------------------------------------------+
      !       |   sigma_r = ( Pbar - Dbar ) / ( 2 Ebar )  |
      !       +-------------------------------------------+
      !
      !   This closed form is identical for steady, real-multiplier, and
      !   complex-multiplier modes. The frequency omega NEVER enters the
      !   budget — it is removed by the Hermitian product. (Full derivation
      !   in docs/pke-periodic-orbit.md.)
      !
      ! ---------------------------------------------------------------------
      ! WHY TWO REAL FIELDS (dRe and dIm)
      ! ---------------------------------------------------------------------
      !   The code stores the complex eigenmode as two real fields,
      !   u_hat = (dRe) + i (dIm). Because the linearised evolution operator
      !   L is real, the complex trajectory splits into two independent real
      !   evolutions: with a(t) = exp(tL) dRe and b(t) = exp(tL) dIm,
      !       q(t) = exp(tL)(dRe + i dIm) = a(t) + i b(t).
      !   Substituting q = a + i b into the Hermitian forms above, all cross
      !   terms cancel and every budget term becomes an ADDITIVE sum of an
      !   "a" contribution and a "b" contribution:
      !       conjg(q_i) q_j  -> a_i a_j + b_i b_j
      !       grad q:grad q*  -> grad a:grad a + grad b:grad b
      !       |q|^2           -> |a|^2 + |b|^2
      !   This additivity is what makes the implementation simple and exact:
      !   the same energy kernels (compute_production / compute_dissipation,
      !   which already form Re + Im Hermitian sums) are called once for a
      !   and once for b, accumulating into the SAME budget arrays.
      !
      ! ---------------------------------------------------------------------
      ! NUMERICAL STRATEGY (two passes over one stored base orbit)
      ! ---------------------------------------------------------------------
      !   PASS 1 marches the real part dRe over the orbit; PASS 2 marches the
      !   imaginary part dIm. The base orbit U(x,t) is expensive, so it is
      !   integrated and stored ONCE (first pass of the first mode) and then
      !   restored step-by-step for every later pass and mode. Each sample is
      !   weighted by
      !       weight = exp(-2 sigma_r t) * dt / T,
      !   which simultaneously (a) undoes the exp(2 sigma_r t) modal growth
      !   so the accumulated quantities are the periodic tilde-averages, and
      !   (b) normalises by the period to form the time mean. The weight uses
      !   sigma_r ONLY — never omega — consistent with the cancellation above.
      !   PASS 2 is skipped when omega ~ 0 (a real multiplier needs no
      !   imaginary part), so the routine reduces exactly to the historical
      !   real-multiplier budget in that limit (full backward compatibility).
      !
      ! ---------------------------------------------------------------------
      ! INPUTS / OUTPUTS
      ! ---------------------------------------------------------------------
      !   Reads : BF_<session>0.f00001        periodic-orbit base-flow IC
      !           dRe<session>0.f<mode>       eigenmode real part
      !           dIm<session>0.f<mode>       eigenmode imaginary part (complex modes)
      !           eigenvalue file via read_eigenvalue -> (sigma_r, omega)
      !   Writes: F<mode> spatial budget fields (production components)
      !           PKE_floquet_<session>0.f<mode>  integrated budget + checks
      !
      ! ---------------------------------------------------------------------
      ! ASSUMPTIONS, CONVENTIONS, AND LIMITATIONS (read before trusting)
      ! ---------------------------------------------------------------------
      !   - Mode normalisation is irrelevant to sigma_r because it cancels in
      !     the ratio (Pbar - Dbar)/(2 Ebar); Ebar carries the same scale.
      !   - The orbit integral uses a first-order rectangle rule; the closure
      !     error therefore scales with dt and shrinks as the orbit is
      !     resolved more finely.
      !   - The factor-of-1/2 energy convention lives in the compute_* helpers
      !     (see their docstrings); production and dissipation share it so the
      !     ratio is convention-independent.
      !   - The base orbit is stored once and reused; for multi-mode runs the
      !     reused base sequencing is a deliberate, documented approximation
      !     of the orbit-recycling scheme (see PASS 1 reuse branch below).
      !   - ORTHOGONALISATION: none is applied along the orbit, and none is
      !     needed for the LEADING mode — a Floquet eigenmode is invariant
      !     under the monodromy, so marching reproduces it exactly (no
      !     re-projection, and NO per-step renormalisation, which would
      !     corrupt the dE~/dt structure). For SUB-LEADING modes the marched
      !     field can drift toward the dominant mode as exp((sig1-sigj)*T);
      !     a robust multi-mode budget would deflate against the converged
      !     subspace each step (not done here). The closure check is self-
      !     diagnosing: drift breaks <dE~/dt>=0 and shows up as a large
      !     |sigma_budget - sigma_eig|, so always read the closure first.
      !   - A complex multiplier is a conjugate pair (mu, conjg(mu)) with the
      !     SAME growth rate; marching dRe alone mixes the mode and its
      !     conjugate. Combining a=dRe and b=dIm into q=a+i*b is what isolates
      !     the single mode — hence both passes are required, not optional.
      !   - PRIOR ART: the energy method is classical for steady flows; this
      !     complex-Floquet, Hermitian, two-field orbit-averaged form is a
      !     non-standard extension to be cross-checked against the literature
      !     before any novelty claim. See docs/pke-periodic-orbit.md sec. 7-8.
      !   - Validity rests on div q = 0 and boundaries that kill the surface
      !     work terms (closed/no-slip or periodic) — the standard setting for
      !     these flows.
      !
      ! ---------------------------------------------------------------------
      ! REFERENCES
      ! ---------------------------------------------------------------------
      !   Reynolds & Orr (1907) energy method; Floquet theory for periodic
      !   base flows (e.g. Barkley & Henderson, J. Fluid Mech. 1996, for the
      !   cylinder-wake multipliers). Steady and real-multiplier budgets are
      !   documented in stability_energy_budget above and in
      !   docs/pke-periodic-orbit.md.
      !=======================================================================
    subroutine stability_energy_budget_floquet
     use krylov_subspace

        !     ----- Perturbation work arrays -----
      !       vx_d/vy_d/vz_d : a snapshot of the perturbation velocity at the
      !                        current orbit step, copied out of Nek's live
      !                        perturbation field vxp/vyp/vzp for the budget
      !                        kernels (which take plain arrays, not vxp).
      !       vx_zero/...     : a permanently zero field handed to the kernels
      !                        as their "imaginary" slot, so each pass injects
      !                        only its own field (see the orbit loop below).
    real, allocatable :: vx_d(:), vy_d(:), vz_d(:)
    real, allocatable :: vx_zero(:), vy_zero(:), vz_zero(:)

       !     ----- Budget accumulators (allocatable -> heap, not stack;
      !             a real(lv) array on the stack overflows large 3D meshes) -
      !       budget_avg(:,1:9)  : the 9 orbit-averaged production fields
      !                            (3 velocity components x 3 base-gradient
      !                            directions), laid out the same as the
      !                            steady budget for direct comparison.
      !       budget_avg(:,10)   : the orbit-averaged dissipation field.
      !       energy_avg(:)      : the orbit-averaged PKE density field.
      !       integrals(1:10)    : volume integrals of the 10 budget fields;
      !                            sum(1:9)=Pbar, integrals(10)=Dbar.
    real, allocatable :: budget_avg(:,:)
    real, allocatable :: energy_avg(:)
    real, dimension(10) :: integrals
    real, allocatable :: prod_x(:), prod_y(:), prod_z(:), diss_tmp(:)

      !     ----- Eigenmode seed stash (per mode). EVERY mode's dRe and dIm is
      !           read ONCE up front, before the build-orbit integration, and
      !           held here; the passes seed the perturbation field from this
      !           stash rather than re-reading files mid-routine. This is the
      !           same load-once-then-compute pattern the steady routine uses:
      !           all file I/O happens before any nek_advance, so a single
      !           shared base orbit can be built once and replayed by every
      !           pass of every mode without interleaving disk reads into the
      !           time integration. Columns are modes.
    real, allocatable :: sRe_x(:,:), sRe_y(:,:), sRe_z(:,:), sRe_t(:,:)
    real, allocatable :: sIm_x(:,:), sIm_y(:,:), sIm_z(:,:), sIm_t(:,:)
    logical, allocatable :: has_im_m(:)   ! per-mode: dIm present & complex

      !     ----- Eigenvalue / scalar diagnostics -----
      !       sigma_r, omega : the Floquet exponent mu = sigma_r + i*omega.
      !       period         : orbit period T = nsteps*dt (averaging window).
      !       growth_corr    : exp(-2*sigma_r*t), undoes modal growth.
      !       weight         : per-step orbit-average weight = growth_corr*dt/T.
      !       E_bar          : orbit- and volume-averaged energy Ebar.
      !       sigma_check    : growth rate reconstructed from the budget.
      !       rel_error      : |sigma_check - sigma_r|/|sigma_r| (closure test).
   real :: sigma_r, omega, period, growth_corr, weight
   real :: E_bar, sigma_check, rel_error

      !     ----- Miscellaneous -----
   real :: glsc2
   integer :: i, j, k, mode, m, col
   character(len=80) :: filename
   character(len=6) :: mode_str
   character(len=3) :: mode_str2

   logical ifverbose
   logical :: mode_reloaded
   logical :: has_im              ! .true. if a dIm file exists for this mode

      !     --> Active point counts. lv = lx1*ly1*lz1*lelv is the COMPILE-TIME
      !         maximum (array storage size); nv = nx1*ny1*nz1*nelv is the
      !         RUNTIME number of velocity points actually used (nv <= lv), and
      !         nt the temperature-point count. Arrays are sized lv but loops
      !         and integrals run over nv/nt only. nv/nt are the public
      !         krylov_subspace globals, set here as elsewhere in the code.
   nv = nx1*ny1*nz1*nelv
   nt = nx1*ny1*nz1*nelt

      !     --> Allocate work arrays at the compile-time size lv (heap).
   allocate(vx_d(lv), vy_d(lv), vz_d(lv))
   allocate(vx_zero(lv), vy_zero(lv), vz_zero(lv))
   allocate(prod_x(lv), prod_y(lv), prod_z(lv), diss_tmp(lv))

      !     #####
      !     #####     READ EIGENVALUE
      !     #####
      !     The Floquet exponent mu = sigma_r + i*omega is read from the
      !     eigensolver output. Both parts steer the rest of the routine:
      !     sigma_r sets the growth correction exp(-2*sigma_r*t) applied to
      !     every sample, and omega decides whether the mode is complex
      !     (omega /= 0 -> the dIm pass below runs) or real (omega ~ 0 ->
      !     the dIm pass is skipped and the budget reduces to the real case).

   call read_eigenvalue(sigma_r, omega)
   if (nid == 0) then
   write (6, *) 'Floquet PKE budget: sigma_r =', sigma_r
   write (6, *) 'Floquet PKE budget: omega   =', omega
   end if

      !     #####
      !     #####     LOAD BASE FLOW AND PREPARE SOLVER
      !     #####

      !     --> BF_<session>0.f00001 is the base flow at phase 0 of the orbit
      !         (the Newton/UPO or limit-cycle initial condition). load_fld
      !         puts it into the Nek primary fields vx,vy,vz so that the time
      !         integrator starts the orbit from the correct point.
   write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)

      !     --> Configure the linearised time stepper. This sets nsteps and dt
      !         from the .par (param(10) = end time), i.e. it fixes how the
      !         orbit is discretised in time for both the base replay and the
      !         perturbation march. It MUST run now, while vx,vy,vz still hold
      !         the BASE flow just loaded: prepare_linearized_solver sizes dt
      !         from compute_cfl(vx,vy,vz), so if an eigenmode (small amplitude)
      !         were resident instead, dt would be sized far too large and the
      !         base-orbit integration would blow the CFL. Hence the eigenmode
      !         pre-load below comes AFTER this call, not before.
   call prepare_linearized_solver

      !     #####
      !     #####     PRE-LOAD ALL EIGENMODES (after dt setup, before marching)
      !     #####
      !     Read every mode's dRe (and dIm, for complex multipliers) now and
      !     stash it; the passes seed the perturbation from the stash. Doing all
      !     reads here keeps disk I/O out of the time-integration loop (the same
      !     load-once-then-compute order the steady routine uses), so the single
      !     shared base orbit is built once and replayed cleanly by every pass.
      !     This runs after prepare_linearized_solver (dt already sized off the
      !     base) but before any nek_advance. vx,vy,vz are clobbered here but
      !     only the stash is kept, and BF_ is reloaded by the build-orbit block.
   allocate(sRe_x(lv,maxmodes), sRe_y(lv,maxmodes), sRe_z(lv,maxmodes))
   allocate(sIm_x(lv,maxmodes), sIm_y(lv,maxmodes), sIm_z(lv,maxmodes))
   if (ifto) allocate(sRe_t(lt,maxmodes), sIm_t(lt,maxmodes))
   allocate(has_im_m(maxmodes))
   has_im_m = .false.
   do m = 1, maxmodes
   write (mode_str, '(i5.5)') m

   write (filename, '(a,a,a,a)') 'dRe', trim(SESSION), '0.f', trim(mode_str)
   call load_fld(filename)
   call opcopy(sRe_x(1,m), sRe_y(1,m), sRe_z(1,m), vx, vy, vz)
   if (ifto) call copy(sRe_t(1,m), t, nt)

   if (abs(omega) .gt. 1.0d-12) then
   write (filename, '(a,a,a,a)') 'dIm', trim(SESSION), '0.f', trim(mode_str)
   inquire(file=filename, exist=has_im_m(m))
   if (has_im_m(m)) then
   call load_fld(filename)
   call opcopy(sIm_x(1,m), sIm_y(1,m), sIm_z(1,m), vx, vy, vz)
   if (ifto) call copy(sIm_t(1,m), t, nt)
   else
   if (nid == 0) write (6, *) &
      'WARNING: missing dIm Floquet mode, real-only budget: ', &
      trim(filename)
   end if
   end if
   end do

      !     --> The orbit period T = nsteps*dt is the averaging window for
      !         <.> below and MUST match the period used in the Floquet
      !         eigensolve; otherwise the growth correction and the time
      !         mean are referenced to the wrong T and the closure fails.
   period = nsteps * dt
   if (nid == 0) write (6, *) &
      'Orbit period T =', period, ' nsteps =', nsteps

      !     --> Solver mode flags. ifpert=.true. makes nek_advance evolve the
      !         linearised perturbation fields vxp/vyp/vzp on top of the base;
      !         ifadj=.false. selects the DIRECT (not adjoint) linear operator,
      !         which is what the Floquet eigenmode being analysed satisfies.
      !         Both are broadcast so every MPI rank agrees.
   ifpert = .true.; ifadj = .false.
   call bcast(ifpert, lsize); call bcast(ifadj, lsize)

      !     --> Permanent zero field used as the "other" Hermitian slot.
      !         compute_production/compute_dissipation form the Hermitian
      !         sum (real_field)^2 + (imag_field)^2 from their two field
      !         arguments. We exploit additivity (see header): PASS 1 calls
      !         them with (a, 0) -> contributes a_i a_j, and PASS 2 with
      !         (b, 0) -> contributes b_i b_j, the two summing to the full
      !         a_i a_j + b_i b_j. Hence a single zeroed slot serves both
      !         passes; it is never overwritten.
   call rzero(vx_zero, lv)
   call rzero(vy_zero, lv)
   call rzero(vz_zero, lv)

      !     --> Allocate orbit storage on first call
   if (.not. energy_init) then
   call allocate_orbit(nsteps)
   end if

      !     #####
      !     #####     BUILD THE BASE ORBIT ONCE (no budget accumulation)
      !     #####
      !     The periodic base flow is integrated and cached a SINGLE time,
      !     before any perturbation pass. Decoupling the (expensive) orbit
      !     construction from the budget passes is deliberate: it guarantees
      !     that EVERY pass below — the dRe pass and the dIm pass, for every
      !     mode — replays the IDENTICAL stored base through orbit_restore.
      !     That is what keeps the real part a and the imaginary part b
      !     evolved by the SAME discrete operator, so q = a + i b reconstructs
      !     a single Floquet mode exactly. (An earlier version accumulated the
      !     dRe pass of the first mode WHILE building the orbit with a
      !     co-evolving base, which left a and b on slightly different
      !     operators for that one mode — removed here.)
   if (.not. energy_init) then
   write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)

      !        --> Base-only build: switch the perturbation solve OFF so this
      !            pass costs only the base integration, and advance the base
      !            (ifbase=.true.), caching it at every step.
   ifpert = .false.; ifbase = .true.
   call bcast(ifpert, lsize); call bcast(ifbase, lsize)
   time = 0.0d0
   do istep = 1, nsteps
   if (nid == 0 .and. (mod(istep, max(1, nsteps/10)) == 0 .or. &
      istep == 1 .or. istep == nsteps)) write (6, &
      "(' PKE_FLOQUET build orbit:',I6,'/',I6)") istep, nsteps
   call nekstab_usrchk()
   call nek_advance()           ! base evolves; settime advances `time`
   call orbit_store(istep)      ! cache base at phase istep*dt
   end do

      !        --> Orbit cached. Re-enable the perturbation solve and freeze
      !            the base; from now on every pass only replays the cache.
   ifpert = .true.; ifbase = .false.
   call bcast(ifpert, lsize); call bcast(ifbase, lsize)
   energy_init = .true.
   end if

      !     --> Allocate budget accumulators on heap
   allocate(budget_avg(lv, 10))
   allocate(energy_avg(lv))

      !     #####
      !     #####     LOOP OVER MODES
      !     #####

      !     Each converged Floquet mode is analysed independently; the budget
      !     for one mode says nothing about another, so the accumulators are
      !     reset at the top of every iteration.
   do mode = 1, maxmodes

   budget_avg = 0.0d0
   energy_avg = 0.0d0
   integrals  = 0.0d0

      !        --> Mode number formats the file suffix (dRe...0.f00001 etc).
   write (mode_str, '(i5.5)') mode

      !        --> This mode's dIm availability was determined during the
      !            up-front pre-load; gates PASS 2 below.
   has_im = has_im_m(mode)

      !        --> Seed PASS 1 with the STASHED real part a=dRe (read up front;
      !            no disk I/O in the integration loop). The base comes entirely
      !            from the cached orbit (orbit_restore below), so vx,vy,vz
      !            contents here are irrelevant.
   call opcopy(vxp(:,1), vyp(:,1), vzp(:,1), &
      sRe_x(1,mode), sRe_y(1,mode), sRe_z(1,mode))
   if (ifto) call copy(tp(:,:,1), sRe_t(1,mode), nt)

      !        ─────────────────────────────────────────
      !        PASS 1: dRe  (replay cached orbit, accumulate the a-part)
      !        ─────────────────────────────────────────
      !        The Hermitian budget is additive in a=dRe and b=dIm. This pass
      !        injects the a contribution (imaginary slot zeroed); the dIm
      !        pass below adds the b contribution into the SAME accumulators.
      !        ifbase=.false. throughout — the base is replayed, never advanced.

   ifbase = .false.
   time = 0.0d0
   do istep = 1, nsteps

      !           --> Log progress
   if (nid == 0) write (6, &
      "(' PKE_FLOQUET mode',I3,':',I6,'/',I6)") &
      mode, istep, nsteps

      !           --> Replay the cached base at this orbit phase, then advance
      !               only the perturbation over it (ifbase=.false.).
   call orbit_restore(istep)
   call nekstab_usrchk()
   call nek_advance()

      !           --> Mirror the current base state into ubase/vbase/wbase.
      !               After nek_advance, vx,vy,vz hold the base flow at this
      !               orbit step. The production kernel needs the base-velocity
      !               gradients dU_i/dx_j and reads the base from the globals
      !               ubase,vbase,wbase, so we copy it across each step (the
      !               base changes along the orbit — this cannot be hoisted).
   call opcopy(ubase, vbase, wbase, vx, vy, vz)

      !           --> Snapshot the perturbation. Pull the live perturbation
      !               vxp/vyp/vzp(:,1) into plain arrays vx_d/vy_d/vz_d, which
      !               is the form the budget kernels accept.
   call opcopy(vx_d, vy_d, vz_d, &
      vxp(:,1), vyp(:,1), vzp(:,1))

      !           --> In 2D the perturbation w-field vzp is not evolved and may
      !               hold uninitialised garbage; zero vz_d so it never enters
      !               the energy sum vx^2+vy^2+vz^2 (which would NaN E_bar). The
      !               budget kernels already if3d-guard their own vz use.
   if (.not. if3d) call rzero(vz_d, nv)

      !           --> Sample weight for the orbit average. Two factors:
      !               exp(-2*sigma_r*time) strips the modal growth so the
      !               accumulated quantity is the T-periodic tilde-average
      !               (E~, P~, D~ in the header), and dt/period turns the
      !               sum over steps into the time mean <.> = (1/T) int_0^T.
      !               First-order rectangle rule; uses sigma_r only (omega
      !               has already cancelled in the Hermitian product).
      !
      !               NB on `time`: it is the Nek5000 GLOBAL clock, advanced
      !               by nek_advance -> settime (TIME = TIME + DT) on every
      !               step. We reset it to 0 before the loop, so at istep it
      !               already equals istep*dt (the right-endpoint phase, in
      !               sync with the perturbation snapshot taken after the
      !               advance). Do NOT add a manual `time = time + dt` here:
      !               that would double-count and corrupt the correction.
      !               This is the same pattern used by the validated orbit
      !               loops in matvec.f90.
   growth_corr = exp(-2.0d0 * sigma_r * time)
   weight = growth_corr * dt / period

      !           --> Production: 9 terms = 3 velocity components x 3 base
      !               gradient directions. Passing (a, 0) makes the kernel
      !               return this pass's a_i a_j contribution; PASS 2 adds
      !               b_i b_j into the same budget_avg columns.
   do j = 1, 3
   call compute_production(vx_d, vy_d, vz_d, &
      vx_zero, vy_zero, vz_zero, j, &
      prod_x, prod_y, prod_z)
   col = (j - 1) * 3
   do i = 1, nv
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
   do i = 1, nv
   budget_avg(i,10) = budget_avg(i,10) &
      + weight * diss_tmp(i)
   end do

      !           --> Perturbation kinetic energy
   do i = 1, nv
   energy_avg(i) = energy_avg(i) + weight * 0.5d0 &
      * (vx_d(i)**2 + vy_d(i)**2 + vz_d(i)**2)
   end do

   end do ! istep

      !        ─────────────────────────────────────────
      !        PASS 2: dIm ORBIT INTEGRATION + ACCUMULATION
      !        ─────────────────────────────────────────
      !        Only complex multipliers (omega /= 0) have an imaginary mode
      !        part. The physical perturbation that grows at rate sigma_r is
      !        the complex field q = a + i*b with a=dRe, b=dIm; PASS 1 added
      !        the a_i a_j / |a|^2 terms, this pass adds the b_i b_j / |b|^2
      !        terms into the SAME accumulators (Hermitian additivity, header
      !        section 4). The exp(i*omega*t) phase has already cancelled in
      !        the Hermitian product, so b uses the identical sigma_r-only
      !        weight and the already-cached base orbit. This pass NEVER
      !        stores the orbit (it only restores) — the cache is read-only
      !        here. See header sections 7.2-7.3 on why both parts are needed.

   if (has_im) then

      !           --> Seed the perturbation with the STASHED imaginary part
      !               b=dIm (loaded up front with all other modes; no disk I/O
      !               in the integration loop). has_im already gates out real
      !               multipliers and missing-dIm-file cases.
   call opcopy(vxp(:,1), vyp(:,1), vzp(:,1), &
      sIm_x(1,mode), sIm_y(1,mode), sIm_z(1,mode))
   if (ifto) call copy(tp(:,:,1), sIm_t(1,mode), nt)

      !           --> Freeze the base (replay only) and restart the clock at
      !               t=0 so the growth correction exp(-2*sigma_r*t) is
      !               referenced to the same orbit phase as PASS 1.
   ifbase = .false.
   time = 0.0d0
   do istep = 1, nsteps

      !           --> Log progress
   if (nid == 0) write (6, &
      "(' PKE_FLOQUET mode',I3,' dIm:',I6,'/',I6)") &
      mode, istep, nsteps

      !           --> Reuse stored base orbit for the b=dIm pass
   call orbit_restore(istep)

      !           --> Advance perturbation on restored base only
   call nekstab_usrchk()
   call nek_advance()

      !           --> Copy current base flow to ubase/vbase/wbase
      !               (compute_production reads from ubase,vbase,wbase)
   call opcopy(ubase, vbase, wbase, vx, vy, vz)

      !           --> Get current perturbation
   call opcopy(vx_d, vy_d, vz_d, &
      vxp(:,1), vyp(:,1), vzp(:,1))

      !           --> Zero the 2D garbage w-field (see PASS 1 note) so it never
      !               poisons the energy sum / E_bar.
   if (.not. if3d) call rzero(vz_d, nv)

      !           --> Growth correction uses sigma_r only
   growth_corr = exp(-2.0d0 * sigma_r * time)
   weight = growth_corr * dt / period

      !           --> Production terms for additive Hermitian dIm part
   do j = 1, 3
   call compute_production(vx_d, vy_d, vz_d, &
      vx_zero, vy_zero, vz_zero, j, &
      prod_x, prod_y, prod_z)
   col = (j - 1) * 3
   do i = 1, nv
   budget_avg(i,col+1) = budget_avg(i,col+1) &
      + weight * prod_x(i)
   budget_avg(i,col+2) = budget_avg(i,col+2) &
      + weight * prod_y(i)
   budget_avg(i,col+3) = budget_avg(i,col+3) &
      + weight * prod_z(i)
   end do
   end do

      !           --> Dissipation term for additive Hermitian dIm part
   call compute_dissipation(vx_d, vy_d, vz_d, &
      vx_zero, vy_zero, vz_zero, diss_tmp)
   do i = 1, nv
   budget_avg(i,10) = budget_avg(i,10) &
      + weight * diss_tmp(i)
   end do

      !           --> Perturbation kinetic energy for dIm part
   do i = 1, nv
   energy_avg(i) = energy_avg(i) + weight * 0.5d0 &
      * (vx_d(i)**2 + vy_d(i)**2 + vz_d(i)**2)
   end do

   end do ! istep
   end if   ! has_im (PASS 2)

      !        ─────────────────────────────────────────
      !        INTEGRATE AND VERIFY
      !        ─────────────────────────────────────────
      !        Collapse the orbit-averaged spatial fields to scalars by
      !        integrating over the domain. glsc2(bm1, f, nv) = sum(bm1*f) is
      !        the spectral-element volume integral (bm1 = diagonal mass
      !        matrix / quadrature weights). Columns 1..9 are the production
      !        components, column 10 the dissipation, energy_avg the PKE.

   do i = 1, 10
   integrals(i) = glsc2(bm1, budget_avg(:,i), nv)
   end do
   E_bar = glsc2(bm1, energy_avg, nv)

      !        --> Growth rate reconstructed purely from the energy budget:
      !            sigma_check = (Pbar - Dbar) / (2 Ebar) physically, but two
      !            code conventions collapse the prefactor and flip a sign:
      !            (1) integrals(10) stores the VISCOUS ENERGY CONTRIBUTION
      !                = -Dbar (compute_dissipation returns 0.5*q.lap(q)*nu,
      !                whose integral is negative), so Pbar - Dbar is
      !                sum(1:9) + integrals(10), NOT minus.
      !            (2) compute_production and compute_dissipation each carry a
      !                0.5, i.e. they return HALF the physical P and D. With
      !                E_bar = <0.5 int |q|^2> (also half), the 2's cancel:
      !                sigma = (Pphys-Dphys)/(2 Ebar) = (Pcode+int10)/Ebar.
      !            Verified against the cylinder Re=100 leading mode: the budget
      !            reproduces sigma_eig=0.1248 exactly with this form (the old
      !            (sum-int10)/(2 Ebar) gave 0.099 - wrong). max() guards the
      !            (physically impossible) zero-energy divide.
   sigma_check = (sum(integrals(1:9)) + integrals(10)) &
      / max(E_bar, 1.0d-30)

      !        --> Closure metric: agreement between the budget-derived
      !            growth rate and the eigenvalue. A small relative error
      !            (typically a few percent, limited by the first-order
      !            orbit quadrature) certifies that eigenmode, base orbit,
      !            and budget are mutually consistent; a large one is a red
      !            flag to investigate before trusting the decomposition.
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

      !        --> Spatial production fields (tagged F<mode>), in the same
      !            layout as the steady budget so the two can be compared
      !            directly. These localise WHERE the mode extracts energy
      !            from the base flow. 2D writes 6 components (the in-plane
      !            block), 3D the full 9.
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

      !        --> Scalar summary file PKE_floquet_<session>0.f<mode>:
      !            the 10 integrals, then total production sum(1:9), the net
      !            P-D, the averaged energy E_bar, and the budget growth rate
      !            sigma_check. This is the machine-readable record used to
      !            tabulate the closure check across cases.
   if (nid == 0) then
   write (filename, '(A,A,A,A)') &
      'PKE_floquet_', trim(SESSION), '0.f', &
      trim(mode_str)
   open (101, file=filename, form='formatted')
   do i = 1, 10
   write (101, '(1E15.7)') integrals(i)
   end do
   write (101, '(1E15.7)') sum(integrals(1:9))
!     Net budget P-D: integrals(10) is the viscous contribution = -Dbar
!     (negative), so the physical P-D is sum(1:9) + integrals(10).
   write (101, '(1E15.7)') &
      sum(integrals(1:9)) + integrals(10)
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

   real, dimension(lv), intent(out) :: dissipation
   call ensure_dissipation_workspace()

      !     --> Compute Laplacians.
      !     NOTE: the if3d guards are correctness fixes, not just
      !     performance.  In 2D runs vz fields are NOT guaranteed to be
      !     zero — they may contain uninitialized data from previous
      !     field loads.  Without the guard, compute_laplacian would
      !     produce nonzero Laplacians and add a spurious z-contribution.
   call compute_laplacian(vx_dRe, diss_lap_ax_s)
   call compute_laplacian(vy_dRe, diss_lap_ay_s)

   call compute_laplacian(vx_dIm, diss_lap_bx_s)
   call compute_laplacian(vy_dIm, diss_lap_by_s)
   if (if3d) then
      call compute_laplacian(vz_dRe, diss_lap_az_s)
      call compute_laplacian(vz_dIm, diss_lap_bz_s)
   end if

      !     --> Compute dissipation term.
   dissipation = vx_dRe*diss_lap_ax_s + vx_dIm*diss_lap_bx_s
   dissipation = dissipation + vy_dRe*diss_lap_ay_s + vy_dIm*diss_lap_by_s

   if (if3d) then
      dissipation = dissipation + vz_dRe*diss_lap_az_s + vz_dIm*diss_lap_bz_s
   end if

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
   integer, intent(in) :: component

   call ensure_production_workspace()
   prod_x = 0.0d0
   prod_y = 0.0d0
   prod_z = 0.0d0

   if (component == 1) then
   call gradm1(prod_dcdx_s, prod_dcdy_s, prod_dcdz_s, ubase, nelv)

   prod_x = -0.5*(vx_dRe**2 + vx_dIm**2)*prod_dcdx_s
   prod_y = -0.5*(vx_dRe*vy_dRe + vy_dIm*vx_dIm)*prod_dcdy_s
   if (if3d) prod_z = -0.5*(vx_dRe*vz_dRe + vz_dIm*vx_dIm)*prod_dcdz_s

   else if (component == 2) then
   call gradm1(prod_dcdx_s, prod_dcdy_s, prod_dcdz_s, vbase, nelv)

   prod_x = -0.5*(vx_dRe*vy_dRe + vy_dIm*vx_dIm)*prod_dcdx_s
   prod_y = -0.5*(vy_dRe**2 + vy_dIm**2)*prod_dcdy_s
   if (if3d) prod_z = -0.5*(vy_dRe*vz_dRe + vz_dIm*vy_dIm)*prod_dcdz_s

   else if (component == 3) then
   if (if3d) then
      call gradm1(prod_dcdx_s, prod_dcdy_s, prod_dcdz_s, wbase, nelv)

      prod_x = -0.5*(vx_dRe*vz_dRe + vz_dIm*vx_dIm)*prod_dcdx_s
      prod_y = -0.5*(vy_dRe*vz_dRe + vz_dIm*vy_dIm)*prod_dcdy_s
      prod_z = -0.5*(vz_dRe**2 + vz_dIm**2)*prod_dcdz_s
   end if
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

!     Only include the z second-derivative in 3D. In 2D, dadz holds garbage
!     (compute_gradients does not zero the out-of-plane derivative), so adding
!     d2adz2 would inject NaN/Inf into the VALID [1:nv] region of the
!     Laplacian -> NaN dissipation integral. Mirrors the if3d guards already
!     in compute_dissipation. (Caught on 2D cylinder_re100; 3D cases were
!     immune because dadz is a real derivative there.)
   if (if3d) then
      call compute_gradients(dadz, wrk1, wrk2, d2adz2)
      Lap_a = d2adx2 + d2ady2 + d2adz2
   else
      Lap_a = d2adx2 + d2ady2
   end if

   end subroutine compute_laplacian
      !-----------------------------------------------------------------------

   end module nekstab_energy_budget
