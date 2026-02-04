      !---------------------------------------------------------------------
      subroutine nekStab_setDefault
      !     specifying default values for nekStab

         implicit none
         include 'SIZE'
         include 'TOTAL'

         k_dim = 100 ! standard value, increas  in .usr
         schur_tgt = 2 ! schur target for schur step factorizaiton
         eigen_tol = 1.0e-6 ! tolerance for eigenmodes convergence
         schur_del = 0.10d0 !
         maxmodes = 20 ! max number of converged modes to disk
         glob_skip = 10 ! global energy computation skip frequency
         findiff_order = 1 ! finite difference order for the Frechet derivative
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

         ifseed_nois = .true. ! noise as initial seed
         ifseed_symm = .false. ! symmetry initial seed
         ifseed_load = .false. ! loading initial seed (e.g. Re_ )
      !  Note: if ifseed_* all are false, 'useric' subroutine prescribes the initial seed

      !  Define here the probe position for zero-crossing vertical velocity analysis !
         xck = 2.0d0; call bcast(xck, wdsize)
         yck = 0.0d0; call bcast(yck, wdsize)
         zck = 0.0d0; call bcast(zck, wdsize)

      ! Sponge zone parameters (modified from KTH Toolbox)
         xLspg = 0.0d0; call bcast(xLspg, wdsize) ! x left
         xRspg = 0.0d0; call bcast(xRspg, wdsize) ! x right
         yLspg = 0.0d0; call bcast(yLspg, wdsize)
         yRspg = 0.0d0; call bcast(yRspg, wdsize)
         zLspg = 0.0d0; call bcast(zLspg, wdsize)
         zRspg = 0.0d0; call bcast(zRspg, wdsize)
         acc_spg = 0.333d0; call bcast(acc_spg, wdsize) !percentage for the acceleration phase in the sponge (e.g. 1/3)
         spng_st = 0.0d0; call bcast(spng_st, wdsize)

         ifotd = .false.; call bcast(ifotd, lsize)
         otd_printStep = 100; call bcast(otd_printStep, isize)
         otd_gsStep = 10; call bcast(otd_gsStep, isize)
         otd_FTLEPeriod = 0.0; call bcast(otd_FTLEPeriod, wdsize)

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
         modal_nsnap = 100    ! Default snapshot count
         modal_prefix = 'dns' ! Default file prefix

      !     !Broadcast all defaults !
         call bcast(eigen_tol, wdsize) ! wdsize for real
         call bcast(schur_del, wdsize)
         call bcast(epsilon_base, wdsize)
         ! Note: xck, yck, zck already broadcast above
         call bcast(xLspg, wdsize)
         call bcast(xRspg, wdsize)
         call bcast(yLspg, wdsize)
         call bcast(yRspg, wdsize)
         call bcast(zLspg, wdsize)
         call bcast(zRspg, wdsize)
         call bcast(acc_spg, wdsize)
         call bcast(spng_st, wdsize)

         call bcast(schur_tgt, isize) ! isize for integer
         call bcast(maxmodes, isize)
         call bcast(k_dim, isize)
         call bcast(bst_skp, isize)
         call bcast(bst_snp, isize)
         call bcast(glob_skip, isize)
         call bcast(findiff_order, isize)

         call bcast(ifres, lsize) !lsize for boolean
         call bcast(ifvor, lsize)
         call bcast(ifvox, lsize)
         call bcast(ifseed_nois, lsize)
         call bcast(ifseed_symm, lsize)
         call bcast(ifseed_load, lsize)
         call bcast(ifldbf, lsize)
         call bcast(ifbf2D, lsize)
         call bcast(ifstorebase, lsize)
         call bcast(ifdyntol, lsize)

      !  Broadcast new mode flags
         call bcast(ifDNS, lsize)
         call bcast(ifLinDNS, lsize)
         call bcast(ifSFD, lsize)
         call bcast(ifBoostConv, lsize)
         call bcast(ifTDF, lsize)
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
         call bcast(animate_mode_num, isize)
         call bcast(modal_nsnap, isize)

      end subroutine nekStab_setDefault
      !---------------------------------------------------------------------
      subroutine nekStab_init
         use krylov_subspace
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

            call oprzero(fcx, fcy, fcz) ! never comment this!
            call rzero(fct, nx1*ny1*nz1*nelv)

            isNekStabinit = .true.
         elseif (nid == 0) then
            print *, 'NekStab already initialized'
         end if

      end subroutine nekStab_init
      !---------------------------------------------------------------------
      subroutine nekStab
      !  Main dispatcher - routes to appropriate solver based on resolved mode flags
      !  Flags are set by nekStab_resolve_mode (string > flags > uparam priority)
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
            return
         end if

      !  ═══════════════════════════════════════════════════════════════════
      !  MODE 1: Fixed Point Methods (SFD, BoostConv, TDF)
      !  ═══════════════════════════════════════════════════════════════════
         if (ifSFD .or. ifBoostConv .or. ifTDF) then

            call nekStab_outpost
            call nekStab_comment

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
         if (ifEnergyBudget) call stability_energy_budget
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
      !---------------------------------------------------------------------
      subroutine nekStab_resolve_mode
      !  Determines operating mode from three possible sources:
      !    1. nekstab_mode string (highest priority) - human readable
      !    2. Individual if* flags set by user - flexible
      !    3. uparam(1) decoding (lowest priority) - backward compatible
      !
      !  Called in nekStab_init AFTER nekStab_usrchk (where user sets preferences)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         logical :: any_mode_flag_set

      !  ─────────────────────────────────────────────────────────────────
      !  Check if user explicitly set any mode flag in nekStab_usrchk
      !  NOTE: Modal flags (ifpod, ifdmd, ifspod) excluded because they
      !  may be used alongside other modes for snapshot collection
      !  ─────────────────────────────────────────────────────────────────
         any_mode_flag_set = ifDNS .or. ifLinDNS .or.
     &        ifSFD .or. ifBoostConv .or. ifTDF .or.
     &        isNewtonFP .or. isNewtonPO .or. isNewtonPO_T .or.
     &        isDirect .or. isAdjoint .or. isTransientGrowth .or.
     &        isFloquetDirect .or. isFloquetAdjoint .or.
     &        isFloquetTransientGrowth .or.
     &        ifEnergyBudget .or. ifWavemaker .or. ifBFSensitivity .or.
     &        ifForceSensReal .or. ifForceSensImag .or. ifDeltaForcing
     &        .or. ifAnimateMode .or. ifAnimateBFDeform
     &        .or. ifAnimateFloquet .or. ifotd

      !  ─────────────────────────────────────────────────────────────────
      !  Priority 1: String mode (nekstab_mode) - highest priority
      !  ─────────────────────────────────────────────────────────────────
         if (len_trim(nekstab_mode) > 0) then
            call nekStab_mode_from_string(nekstab_mode)
            if (nid == 0) write (6, *) 'Mode set via nekstab_mode = ',
     &           trim(nekstab_mode)

      !  ─────────────────────────────────────────────────────────────────
      !  Priority 2: Flag mode - user set explicit flags
      !  ─────────────────────────────────────────────────────────────────
         elseif (any_mode_flag_set) then
            call nekStab_mode_from_flags
            if (nid == 0) write (6, *) 'Mode set via if-flags'

      !  ─────────────────────────────────────────────────────────────────
      !  Priority 3: uparam(1) decoding - backward compatible default
      !  ─────────────────────────────────────────────────────────────────
         else
            call nekStab_mode_from_uparam
            if (nid == 0) write (6, *) 'Mode set via uparam(1) =',
     &           uparam(1)
         end if

      !  Validate: check for conflicting modes
         call nekStab_validate_mode

      end subroutine nekStab_resolve_mode
      !---------------------------------------------------------------------
      subroutine nekStab_mode_from_string(mode_str)
      !  Parses nekstab_mode string and sets appropriate flags
      !  Case-insensitive comparison for user convenience
         implicit none
         include 'SIZE'
         include 'TOTAL'
         character(len=*), intent(in) :: mode_str
         character(len=32) :: mode_lower
         integer :: i

      !  Convert to lowercase for case-insensitive comparison
         mode_lower = adjustl(mode_str)
         do i = 1, len_trim(mode_lower)
            if (mode_lower(i:i) >= 'A' .and. mode_lower(i:i) <= 'Z')
     &           then
               mode_lower(i:i) = char(ichar(mode_lower(i:i)) + 32)
            end if
         end do

      !  Match mode string and set corresponding flag
         select case (trim(mode_lower))

      !  Mode 0: DNS
         case ('dns')
            ifDNS = .true.
         case ('linear_dns', 'lindns', 'linearized_dns')
            ifLinDNS = .true.

      !  Mode 1: Fixed point methods
         case ('sfd')
            ifSFD = .true.
         case ('boostconv', 'boost')
            ifBoostConv = .true.
         case ('tdf')
            ifTDF = .true.

      !  Mode 2: Newton-Krylov
         case ('newton_fp', 'newton')
            isNewtonFP = .true.
         case ('newton_po', 'upo')
            isNewtonPO = .true.
         case ('newton_po_t', 'forced_upo')
            isNewtonPO_T = .true.

      !  Mode 3: Stability analysis
         case ('direct')
            isDirect = .true.
         case ('floquet_direct', 'floquetdirect')
            isFloquetDirect = .true.
         case ('adjoint')
            isAdjoint = .true.
         case ('floquet_adjoint', 'floquetadjoint')
            isFloquetAdjoint = .true.
         case ('transient_growth', 'tg')
            isTransientGrowth = .true.
         case ('floquet_tg', 'floquet_transient_growth')
            isFloquetTransientGrowth = .true.

      !  Mode 4: Postprocessing
         case ('energy_budget')
            ifEnergyBudget = .true.
         case ('wavemaker')
            ifWavemaker = .true.
         case ('bf_sensitivity', 'baseflow_sensitivity')
            ifBFSensitivity = .true.
         case ('force_sensitivity_real', 'force_sens_real')
            ifForceSensReal = .true.
         case ('force_sensitivity_imag', 'force_sens_imag')
            ifForceSensImag = .true.
         case ('delta_forcing')
            ifDeltaForcing = .true.
         case ('animate_mode', 'animate')
            ifAnimateMode = .true.
         case ('animate_bf_deform', 'animate_deform')
            ifAnimateBFDeform = .true.
         case ('animate_floquet')
            ifAnimateFloquet = .true.

      !  Mode 5: OTD
         case ('otd')
            ifotd = .true.

      !  Mode 6: Modal analysis
         case ('pod')
            ifpod = .true.
         case ('dmd')
            ifdmd = .true.
         case ('spod')
            ifspod = .true.

         case default
            if (nid == 0) then
               write (6, *) 'ERROR: Unknown nekstab_mode: ',
     &              trim(mode_str)
               write (6, *) 'Valid modes: dns, linear_dns, sfd, ',
     &              'boostconv, tdf,'
               write (6, *) '  newton_fp, newton_po, newton_po_t,'
               write (6, *) '  direct, adjoint, transient_growth,'
               write (6, *) '  floquet_direct, floquet_adjoint, ',
     &              'floquet_tg,'
               write (6, *) '  energy_budget, wavemaker, ',
     &              'bf_sensitivity,'
               write (6, *) '  animate_mode, otd, pod, dmd, spod'
            end if
            call nek_end
         end select

      end subroutine nekStab_mode_from_string
      !---------------------------------------------------------------------
      subroutine nekStab_mode_from_flags
      !  Handles ifFloquet modifier flag transformation
      !  ifFloquet=.true. + isDirect=.true. → isFloquetDirect=.true.
         implicit none
         include 'SIZE'
         include 'TOTAL'

      !  Apply ifFloquet modifier to base stability modes
         if (ifFloquet) then
            if (isDirect) then
               isFloquetDirect = .true.
               isDirect = .false.
            end if
            if (isAdjoint) then
               isFloquetAdjoint = .true.
               isAdjoint = .false.
            end if
            if (isTransientGrowth) then
               isFloquetTransientGrowth = .true.
               isTransientGrowth = .false.
            end if
         end if

      end subroutine nekStab_mode_from_flags
      !---------------------------------------------------------------------
      subroutine nekStab_mode_from_uparam
      !  Decodes uparam(1) into mode flags (backward compatible)
      !  Uses tolerance-based comparison to avoid floating-point issues
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real :: up1
         real, parameter :: tol = 1.0e-4

         up1 = uparam(1)

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 0: DNS
      !  ─────────────────────────────────────────────────────────────────
         if (abs(up1 - 0.0) < tol) then
            ifDNS = .true.
         elseif (abs(up1 - 0.1) < tol) then
            ifLinDNS = .true.

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 1: Fixed point methods
      !  ─────────────────────────────────────────────────────────────────
         elseif (abs(up1 - 1.1) < tol) then
            ifSFD = .true.
         elseif (abs(up1 - 1.2) < tol) then
            ifBoostConv = .true.
         elseif (abs(up1 - 1.4) < tol) then
            ifTDF = .true.

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 2: Newton-Krylov
      !  ─────────────────────────────────────────────────────────────────
         elseif (abs(up1 - 2.0) < tol) then
            isNewtonFP = .true.
         elseif (abs(up1 - 2.1) < tol) then
            isNewtonPO = .true.
         elseif (abs(up1 - 2.2) < tol) then
            isNewtonPO_T = .true.

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 3: Eigenvalue problem (stability analysis)
      !  ─────────────────────────────────────────────────────────────────
         elseif (abs(up1 - 3.1) < tol) then
            isDirect = .true.
         elseif (abs(up1 - 3.11) < tol) then
            isFloquetDirect = .true.
         elseif (abs(up1 - 3.2) < tol) then
            isAdjoint = .true.
         elseif (abs(up1 - 3.21) < tol) then
            isFloquetAdjoint = .true.
         elseif (abs(up1 - 3.3) < tol) then
            isTransientGrowth = .true.
         elseif (abs(up1 - 3.31) < tol) then
            isFloquetTransientGrowth = .true.

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 4: Postprocessing
      !  ─────────────────────────────────────────────────────────────────
         elseif (abs(up1 - 4.0) < tol) then
            ifEnergyBudget = .true.
            ifWavemaker = .true.
            ifBFSensitivity = .true.
         elseif (abs(up1 - 4.1) < tol) then
            ifEnergyBudget = .true.
         elseif (abs(up1 - 4.2) < tol) then
            ifWavemaker = .true.
         elseif (abs(up1 - 4.3) < tol) then
            ifBFSensitivity = .true.
         elseif (abs(up1 - 4.41) < tol) then
            ifForceSensReal = .true.
         elseif (abs(up1 - 4.42) < tol) then
            ifForceSensImag = .true.
         elseif (abs(up1 - 4.43) < tol) then
            ifDeltaForcing = .true.
         elseif (abs(up1 - 4.50) < tol) then
            ifAnimateMode = .true.
            animate_mode_num = int(uparam(7))
         elseif (abs(up1 - 4.51) < tol) then
            ifAnimateBFDeform = .true.
            animate_mode_num = int(uparam(7))
         elseif (abs(up1 - 4.52) < tol) then
            ifAnimateFloquet = .true.
            animate_mode_num = int(uparam(7))

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 5: OTD (tolerance-based, consistent with other modes)
      !  ─────────────────────────────────────────────────────────────────
         elseif (abs(up1 - 5.0) < tol) then
            ifotd = .true.

      !  ─────────────────────────────────────────────────────────────────
      !  Mode 6: Modal analysis
      !  ─────────────────────────────────────────────────────────────────
         elseif (abs(up1 - 6.1) < tol) then
            ifpod = .true.
         elseif (abs(up1 - 6.2) < tol) then
            ifdmd = .true.
         elseif (abs(up1 - 6.3) < tol) then
            ifspod = .true.

         end if

      end subroutine nekStab_mode_from_uparam
      !---------------------------------------------------------------------
      subroutine nekStab_validate_mode
      !  Validates that only one main mode category is active
      !  Prevents conflicting modes (e.g., direct + adjoint)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer :: nmodes

         nmodes = 0

      !  Count active mode categories
         if (ifDNS .or. ifLinDNS) nmodes = nmodes + 1
         if (ifSFD .or. ifBoostConv .or. ifTDF) nmodes = nmodes + 1
         if (isNewtonFP .or. isNewtonPO .or. isNewtonPO_T)
     &        nmodes = nmodes + 1
         if (isDirect .or. isFloquetDirect .or.
     &       isAdjoint .or. isFloquetAdjoint .or.
     &       isTransientGrowth .or. isFloquetTransientGrowth)
     &        nmodes = nmodes + 1
         if (ifEnergyBudget .or. ifWavemaker .or. ifBFSensitivity .or.
     &       ifForceSensReal .or. ifForceSensImag .or. ifDeltaForcing
     &       .or. ifAnimateMode .or. ifAnimateBFDeform
     &       .or. ifAnimateFloquet) nmodes = nmodes + 1
         if (ifotd) nmodes = nmodes + 1
      !  Note: Modal analysis (POD/DMD/SPOD) not counted - can run standalone

         if (nmodes > 1) then
            if (nid == 0) then
               write (6, *) 'ERROR: Multiple conflicting modes active'
               write (6, *) 'Set only ONE mode category at a time'
            end if
            call nek_end
         end if

      !  Check for conflicting stability sub-modes
         if ((isDirect .or. isFloquetDirect) .and.
     &       (isAdjoint .or. isFloquetAdjoint)) then
            if (nid == 0) then
               write (6, *) 'ERROR: Both direct and adjoint modes set'
               write (6, *) 'Choose one: direct OR adjoint'
            end if
            call nek_end
         end if

         if ((isDirect .or. isFloquetDirect .or.
     &        isAdjoint .or. isFloquetAdjoint) .and.
     &       (isTransientGrowth .or. isFloquetTransientGrowth)) then
            if (nid == 0) then
               write (6, *) 'ERROR: Eigenmode and transient growth ',
     &              'both set'
               write (6, *) 'Choose one: eigenmode OR transient growth'
            end if
            call nek_end
         end if

      !  Check for orphaned ifFloquet (set without base mode)
         if (ifFloquet .and. .not. (isFloquetDirect .or.
     &       isFloquetAdjoint .or. isFloquetTransientGrowth)) then
            if (nid == 0) then
               write (6, *) 'ERROR: ifFloquet set without base mode'
               write (6, *) 'Set isDirect, isAdjoint, or ',
     &              'isTransientGrowth with ifFloquet'
            end if
            call nek_end
         end if

      !  Check for no mode selected (nmodes == 0 with no modal analysis)
         if (nmodes == 0 .and. .not. (ifpod .or. ifdmd .or. ifspod))
     &        then
            if (nid == 0) then
               write (6, *) 'WARNING: No operating mode selected'
               write (6, *) 'Defaulting to DNS mode (uparam(1)=0)'
            end if
            ifDNS = .true.
         end if

      end subroutine nekStab_validate_mode
      !---------------------------------------------------------------------
