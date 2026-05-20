      !-----------------------------------------------------------------------
      ! main.f90 — nekStab entry points and dispatcher
      !
      ! Purpose:
      !   Provides the three main entry points called by Nek5000:
      !   nekStab_setDefault (parameter initialization), nekStab_init
      !   (framework startup), and nekStab (mode dispatcher). These
      !   remain bare subroutines (not in a module) because Nek5000
      !   calls them without use statements.
      !
      ! Public interface:
      !   nekStab_setDefault — initialize all defaults
      !   nekStab_init       — framework initialization
      !   nekStab            — main mode dispatcher
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL
      !
      ! See also:
      !   mode_config.f90 — mode resolution and validation
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! nekStab_setDefault — Initialize all default parameters
      !
      ! Purpose:
      !   Sets default values for all nekStab parameters, mode flags,
      !   and configuration variables. Called once at startup before
      !   user overrides in nekStab_usrchk. Broadcasts all values
      !   to ensure MPI consistency.
      !-----------------------------------------------------------------------
      subroutine nekStab_setDefault
         use krylov_subspace, only: TN_MANUAL
         implicit none
         include 'SIZE'
         include 'TOTAL'

         k_dim = 100 ! standard value, increas  in .usr
         schur_tgt = 2 ! schur target for schur step factorizaiton
         eigen_tol = 1.0e-6 ! tolerance for eigenmodes convergence
         schur_del = 0.10d0 !
         maxmodes = 20 ! max number of converged modes to disk
         glob_skip = 10 ! global energy computation skip frequency
         findiff_order = 2 ! finite difference order for the Frechet derivative
         epsilon_base = 1.0e-6 ! finite difference perturbation scale parameter

         bst_skp = 10 ! boostconv skip iterations
         bst_snp = 10 ! bootsconv residual subspace matrix size

         ifres = .false. ! outpost restart files (KRY*, HES*)
         ifvor = .false. ! outpost vorticity (vor* omega_x,omega_y,omega_z components)
         ifvox = .false. ! outpost vortex (vox*: q,lambda2,omega criterions)
         ifldbf = .true. ! load base flow for stability computations
         ifbf2D = .false. ! force 2D base flow solution
         ifstorebase = .true. ! store base flow for Floquet analysis (dynamic allocated)
         ifdyntol = .false. ! dynamical tolerances for SFD and Newton (potential speed-up)
         ew_tol_cap = 0.0d0 ! EW solver cap (0=uncapped, e.g. 1e-5 for conservative)
         thermal_norm_weight = 1.0d0
         thermal_buoyancy_coeff = 0.0d0
         thermal_norm_min = 1.0d0
         thermal_norm_max = 50.0d0
         thermal_norm_mode = TN_MANUAL
      !  For TN_AUTO/TN_CLIP, set thermal_buoyancy_coeff in nekStab_usrchk
      !  using the case scaling (e.g. Ri for buoyant cylinder, Pr*Ra for
      !  thermosyphon). Default TN_MANUAL preserves the historical norm.

         ifseed_nois = .true. ! noise as initial seed
         ifseed_symm = .false. ! symmetry initial seed
         ifseed_load = .false. ! loading initial seed (e.g. Re_ )
      !  Note: if ifseed_* all are false, 'useric' subroutine prescribes the initial seed

      !  Define here the probe position for zero-crossing vertical velocity analysis !
         xck = 2.0d0
         yck = 0.0d0
         zck = 0.0d0

      ! Sponge zone parameters (modified from KTH Toolbox)
         xLspg = 0.0d0 ! x left
         xRspg = 0.0d0 ! x right
         yLspg = 0.0d0
         yRspg = 0.0d0
         zLspg = 0.0d0
         zRspg = 0.0d0
         acc_spg = 0.333d0 ! acceleration phase fraction (e.g. 1/3)
         spng_st = 0.0d0

         ifotd = .false.
         otd_printStep = 100
         otd_gsStep = 10
         otd_FTLEPeriod = 0.0
         otd_convTol = 1.0d-6
         otd_minSteps = 200
         call rzero(FTLEv_prev, lpert)

         evop = '_' ! initialize output prefix

      !  ───────────────────────────────────────────────────────────────
      !  User-settable mode flags (alternative to uparam encoding)
      !  These can be set in nekStab_usrchk to override uparam(1)
      !  ───────────────────────────────────────────────────────────────
         nekstab_mode = ''  ! Empty = use uparam(1) or flag detection

      !  Mode 0: DNS
         ifDNS = .false.     ! Direct Numerical Simulation
         ifLinDNS = .false.  ! Linearized DNS

      !  Mode 1: Fixed point methods
         ifSFD = .false.       ! Selective Frequency Damping
         ifBoostConv = .false. ! BoostConv acceleration
         ifTDF = .false.       ! Time-Delayed Feedback
         ifDMT = .false.       ! Dynamic Mode Tracking

      !  Mode 2: Newton-Krylov (use existing isNewtonFP, isNewtonPO, isNewtonPO_T)
         isNewtonFP = .false.   ! Newton for fixed points
         isNewtonPO = .false.   ! Newton for periodic orbits
         isNewtonPO_T = .false. ! Newton for forced periodic orbits

      !  Mode 3: Stability analysis
         isDirect = .false.           ! Direct eigenmodes
         isAdjoint = .false.          ! Adjoint eigenmodes
         isTransientGrowth = .false.  ! Transient growth
         ifFloquet = .false.          ! Floquet modifier flag
         isFloquetDirect = .false.    ! Floquet direct
         isFloquetAdjoint = .false.   ! Floquet adjoint
         isFloquetTransientGrowth = .false.  ! Floquet transient growth

      !  Mode 4: Postprocessing
         ifEnergyBudget = .false.    ! Stability energy budget
         ifWavemaker = .false.       ! Wavemaker computation
         ifBFSensitivity = .false.   ! Base flow sensitivity
         ifForceSensReal = .false.   ! Steady force sensitivity (real)
         ifForceSensImag = .false.   ! Steady force sensitivity (imag)
         ifDeltaForcing = .false.    ! Delta forcing
         ifAnimateMode = .false.     ! Animate mode only
         ifAnimateBFDeform = .false. ! Animate + base flow deformation
         ifAnimateFloquet = .false.  ! Animate Floquet mode
         animate_mode_num = 1        ! Default mode number for animation

      !  Mode 6: Modal analysis
         ifpod = .false.      ! POD analysis
         ifdmd = .false.      ! DMD analysis
         ifspod = .false.     ! SPOD analysis
         ifwinamp = .true.    ! Amplitude normalization (PySPOD compatible)
         use_cgs = .true.     ! CGS2 orthogonalization (2 gop calls vs 2k for MGS)
         modal_nsnap = 100    ! Default snapshot count
         modal_nsave = 10     ! Default number of modes to save
         modal_dt = 0.1d0     ! Default time between snapshots
         modal_prefix = 'dns' ! Default file prefix
         dmd_rank = 0         ! DMD rank (0=auto based on energy)
         spod_nfft = 64       ! SPOD FFT block size
         spod_noverlap = 32   ! SPOD block overlap (50%)

      !  Broadcast all defaults to ensure MPI consistency
         call bcast_nStab_defaults()

      end subroutine nekStab_setDefault
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! bcast_nStab_defaults — Broadcast all nekStab defaults to MPI ranks
      !
      ! Purpose:
      !   Groups all bcast calls by common block for maintainability.
      !   Called once from nekStab_setDefault after all assignments.
      !-----------------------------------------------------------------------
      subroutine bcast_nStab_defaults()
         implicit none
         include 'SIZE'
         include 'TOTAL'

c        nStab_real (tolerances, domain bounds)
         call bcast(eigen_tol, wdsize)
         call bcast(schur_del, wdsize)
         call bcast(epsilon_base, wdsize)
         call bcast(ew_tol_cap, wdsize)
         call bcast(thermal_norm_weight, wdsize)
         call bcast(thermal_buoyancy_coeff, wdsize)
         call bcast(thermal_norm_min, wdsize)
         call bcast(thermal_norm_max, wdsize)

c        nStab_sponge
         call bcast(xLspg, wdsize)
         call bcast(xRspg, wdsize)
         call bcast(yLspg, wdsize)
         call bcast(yRspg, wdsize)
         call bcast(zLspg, wdsize)
         call bcast(zRspg, wdsize)
         call bcast(acc_spg, wdsize)
         call bcast(spng_st, wdsize)

c        nStab_fd (probe position, finite difference)
         call bcast(xck, wdsize)
         call bcast(yck, wdsize)
         call bcast(zck, wdsize)
         call bcast(findiff_order, isize)

c        nStab_int
         call bcast(schur_tgt, isize)
         call bcast(maxmodes, isize)
         call bcast(glob_skip, isize)
         call bcast(thermal_norm_mode, isize)

c        nStab_boostconv
         call bcast(k_dim, isize)
         call bcast(bst_skp, isize)
         call bcast(bst_snp, isize)

c        nStab_logical
         call bcast(ifres, lsize)
         call bcast(ifvor, lsize)
         call bcast(ifvox, lsize)
         call bcast(ifseed_nois, lsize)
         call bcast(ifseed_symm, lsize)
         call bcast(ifseed_load, lsize)
         call bcast(ifldbf, lsize)
         call bcast(ifbf2D, lsize)
         call bcast(ifstorebase, lsize)
         call bcast(ifdyntol, lsize)
         call bcast(ifotd, lsize)

c        nStab_mode_flags
         call bcast(ifDNS, lsize)
         call bcast(ifLinDNS, lsize)
         call bcast(ifSFD, lsize)
         call bcast(ifBoostConv, lsize)
         call bcast(ifTDF, lsize)
         call bcast(ifDMT, lsize)
         call bcast(ifFloquet, lsize)
         call bcast(isNewtonFP, lsize)
         call bcast(isNewtonPO, lsize)
         call bcast(isNewtonPO_T, lsize)
         call bcast(isDirect, lsize)
         call bcast(isAdjoint, lsize)
         call bcast(isTransientGrowth, lsize)
         call bcast(isFloquetDirect, lsize)
         call bcast(isFloquetAdjoint, lsize)
         call bcast(isFloquetTransientGrowth, lsize)
         call bcast(ifEnergyBudget, lsize)
         call bcast(ifWavemaker, lsize)
         call bcast(ifBFSensitivity, lsize)
         call bcast(ifForceSensReal, lsize)
         call bcast(ifForceSensImag, lsize)
         call bcast(ifDeltaForcing, lsize)
         call bcast(ifAnimateMode, lsize)
         call bcast(ifAnimateBFDeform, lsize)
         call bcast(ifAnimateFloquet, lsize)
         call bcast(ifpod, lsize)
         call bcast(ifdmd, lsize)
         call bcast(ifspod, lsize)
         call bcast(ifwinamp, lsize)
         call bcast(use_cgs, lsize)

c        nStab_mode_int / nStab_mode_real / nStab_mode_char
         call bcast(animate_mode_num, isize)
         call bcast(modal_nsnap, isize)
         call bcast(modal_nsave, isize)
         call bcast(modal_dt, wdsize)
         call bcast(dmd_rank, isize)
         call bcast(spod_nfft, isize)
         call bcast(spod_noverlap, isize)

c        OTD_params
         call bcast(otd_printStep, isize)
         call bcast(otd_gsStep, isize)
         call bcast(otd_FTLEPeriod, wdsize)
         call bcast(otd_convTol, wdsize)
         call bcast(otd_minSteps, isize)

      end subroutine bcast_nStab_defaults
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! nekStab_init — Framework initialization
      !
      ! Purpose:
      !   Initializes the nekStab framework: sets defaults, calls user
      !   configuration hook, resolves operating mode, prints parameters,
      !   computes domain bounds, and prepares forcing arrays. Called
      !   once at istep=0.
      !-----------------------------------------------------------------------
      subroutine nekStab_init
         use krylov_subspace
         use nekstab_mode_config
         use nekstab_diagnostics
         use nekstab_forcing_mod
         use nekstab_vectors, only: zero_forcing
         implicit none
         include 'SIZE'
         include 'TOTAL'
         logical scal
         real glmin, glmax
         integer i

         nv = nx1*ny1*nz1*nelv

         if (.not. isNekStabinit) then
            call nekStab_setDefault
            call nekStab_usrchk ! where user can change defaults or set mode flags
      !     Broadcast user-set values to all MPI ranks (32 chars * csize bytes)
            call bcast(nekstab_mode, 32*csize)
            call bcast(modal_prefix, 3*csize)
            call nekStab_resolve_mode ! resolve mode from string/flags/uparam
            call nekStab_printNEKParams

            xmn = glmin(xm1, nv); xmx = glmax(xm1, nv)
            ymn = glmin(ym1, nv); ymx = glmax(ym1, nv)
            zmn = glmin(zm1, nv); zmx = glmax(zm1, nv)

            if (nid == 0) then
               print *, '                 __   _____  __          __  '
               print *, '   ____   ___   / /__/ ___/ / /_ ____ _ / /_ '
               print *, '  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \'
               print *, ' / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /'
               print *, '/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/ '
               print *, '(c) 2020-2025 DynFluid Laboratoire Paris ', NSVERSION
               print *, 'Nek5000 ', NVERSION
               print *, ''
            end if

            call copy(bm1s, bm1, nv) ! never comment this !
            ifbfcv = .false.

            if (spng_st > 0) call activate_sponge

            nof = 0
            scal = .false.
            do i = 1, size(ifpsco)
               if (ifpsco(i)) then
                  scal = .true.
                  nof = nof + 1
               end if
            end do
            if (ifto .or. scal) then
               if (nid == 0) write (6, *) 'Scalars found:'
               if (nid == 0) write (6, *) ' ifto=', ifto
               if (nid == 0) write (6, *) ' ifpsco=', ifpsco
                if (ifto) nof = nof + 1
                if (nid == 0) write (6, *) 'number of possible scalars (ldimt)=', ldimt
                if (nid == 0) write (6, *) 'number of scalars (nof)=', nof, npscal
             end if

             if (ifheat) call configure_thermal_norm_weight()

             call zero_forcing

             isNekStabinit = .true.
         elseif (nid == 0) then
            print *, 'NekStab already initialized'
         end if

      end subroutine nekStab_init
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! nekStab — Main mode dispatcher
      !
      ! Purpose:
      !   Routes execution to the appropriate solver based on mode flags
      !   resolved by nekStab_resolve_mode. Called every timestep by
      !   Nek5000's userchk.
      !-----------------------------------------------------------------------
      subroutine nekStab
         use nekstab_diagnostics
         use nekstab_noise
         use nekstab_fixedpoint
         use nekstab_newton
         use nekstab_eigensolvers
         use nekstab_energy_budget
         use nekstab_sensitivity
         use nekstab_otd, only: otd
         use nekstab_dmt, only: dmt
         use nekstab_modal_analysis
         use nekstab_vectors, only: zero_forcing
         implicit none
         include 'SIZE'
         include 'TOTAL'

         if (istep == 0) call nekStab_init

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 0: DNS (Direct Numerical Simulation)
      !  ═══════════════════════════════════════════════════════════════════
         if (ifDNS .or. ifLinDNS) then

            if (ifLinDNS) then
               ifbase = .true.; call bcast(ifbase, lsize)
               ifpert = .true.; call bcast(ifpert, lsize)
               param(31) = lpert; npert = int(param(31))
               if (istep == 0) call op_add_noise(vxp, vyp, vzp)
               if (nid == 0) write (6, *) 'Linearized+DNS: ifbase=,',
     &              ifbase, ' ifpert=', ifpert
               if (ifoutfld) call outpost2(vxp, vyp, vzp, prp, tp, nof,
     &              'pr_')
               call nekStab_energy(vxp(:, 1), vyp(:, 1), vzp(:, 1),
     &              tp(:, :, 1), 'total_ene_p1.dat', glob_skip)
               call nekStab_enstrophy(vxp(:, 1), vyp(:, 1), vzp(:, 1),
     &              tp(:, :, 1), 'total_ens_p1.dat', glob_skip)
            end if

             call nekStab_outpost
             call nekStab_comment
             call zero_forcing

            if (ifbfcv) call nek_end
            return
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 1: Fixed-point methods (SFD, BoostConv, TDF)
      !
      !  WARNING: This block must be separate from Mode 0 (DNS).
      !  nekstab_mode='sfd' sets ifSFD=.true. but NOT ifDNS, so nesting
      !  inside if(ifDNS) would make the fixed-point path unreachable.
      !  Nek5000 handles time-stepping; these routines apply damping/feedback.
      !  ═══════════════════════════════════════════════════════════════════
         if (ifSFD .or. ifBoostConv .or. ifTDF .or. ifDMT) then

            call nekStab_outpost
            call nekStab_comment
            !  Must zero fcx/fcy/fcz/fct before SFD/TDF/BoostConv each timestep.
            !  Without zeroing, old forcing values accumulate and corrupt the flow.
            !  Sponge forcing is computed on-the-fly in nekStab_forcing, not stored
            !  in fcx, so no sponge contribution is lost by zeroing here.
            call zero_forcing

            if (ifSFD) then
               call SFD
               if (uparam(5) == 0) call nekStab_energy(vx, vy, vz, t,
     &              'total_energy.dat', glob_skip)
            elseif (ifBoostConv) then
               if (nid == 0) write (6, *) 'BOOSTCONV'
               call BoostConv
            elseif (ifTDF) then
               if (nid == 0) write (6, *) 'TDF'
               call TDF
            elseif (ifDMT) then
               if (nid == 0) write (6, *) 'DMT'
               call dmt
            end if

            if (ifbfcv) call nek_end
            return
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 2: Newton-Krylov Solver
      !  ═══════════════════════════════════════════════════════════════════
         if (isNewtonFP .or. isNewtonPO .or. isNewtonPO_T) then

            if (nid == 0) then
               if (isNewtonFP) then
                  write (6, *) 'Newton-Krylov for fixed points...'
               elseif (isNewtonPO) then
                  write (6, *) 'Newton-Krylov for UPOs...'
               elseif (isNewtonPO_T) then
                  write (6, *) 'Newton-Krylov for forced UPOs...'
               end if
            end if

            call newton_krylov
            call nek_end
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 3: Eigenvalue Problem (Direct/Adjoint/TransientGrowth)
      !  ═══════════════════════════════════════════════════════════════════
         if (isDirect .or. isFloquetDirect .or.
     &       isAdjoint .or. isFloquetAdjoint .or.
     &       isTransientGrowth .or. isFloquetTransientGrowth) then

            call krylov_schur
            call nek_end
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 4: Postprocessing
      !  ═══════════════════════════════════════════════════════════════════
         if (ifEnergyBudget) then
            if (ifFloquet) then
               call stability_energy_budget_floquet
            else
               call stability_energy_budget
            end if
         end if
         if (ifWavemaker) call wave_maker
         if (ifBFSensitivity) call bf_sensitivity
         if (ifForceSensReal .or. ifForceSensImag)
     &        call ts_steady_force_sensitivity
         if (ifDeltaForcing) call delta_forcing
         if (ifAnimateMode)
     &        call animate_mode_only(animate_mode_num, 'd')
         if (ifAnimateBFDeform)
     &        call animate_mode(animate_mode_num, 'd')
         if (ifAnimateFloquet)
     &        call animate_mode_Floquet(animate_mode_num, 'd')

         if (ifEnergyBudget .or. ifWavemaker .or. ifBFSensitivity .or.
     &       ifForceSensReal .or. ifForceSensImag .or. ifDeltaForcing
     &       .or. ifAnimateMode .or. ifAnimateBFDeform
     &       .or. ifAnimateFloquet) then
            call nek_end
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 5: OTD (Optimally Time-Dependent)
      !  ═══════════════════════════════════════════════════════════════════
         if (ifotd) then
            call otd
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 6: Modal Analysis (POD/DMD/SPOD)
      !  ═══════════════════════════════════════════════════════════════════
         if (ifpod .or. ifdmd .or. ifspod) then
            call modal_analysis
            call nek_end
         end if

      end subroutine nekStab
      !-----------------------------------------------------------------------
