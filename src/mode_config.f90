!-----------------------------------------------------------------------
! mode_config.f90 — Operating mode resolution and validation
!
! Purpose:
!   Determines the nekStab operating mode from three possible
!   sources (string, flags, uparam) with a defined priority order.
!   Validates that only one mode category is active.
!
! Public interface:
!   nekStab_resolve_mode    — main entry: resolves mode from all sources
!   nekStab_mode_from_string — parse mode string into flags
!   nekStab_mode_from_flags  — apply ifFloquet modifier
!   nekStab_mode_from_uparam — decode uparam(1) into flags
!   nekStab_validate_mode    — check for conflicting modes
!
! Dependencies:
!   nekstab_nek_bridge (for uparam, nid, common block flags)
!   nekstab_mode_codes (for MODE_* integer constants, uparam_to_mode_code)
!-----------------------------------------------------------------------

module nekstab_mode_config
   use nekstab_nek_bridge, only: nekStab_log, nekStab_error, nekstab_mode, &
                                 uparam, nid, animate_mode_num, ifbuoyancy, &
                                 thermal_buoyancy_coeff, ifDNS, ifLinDNS, ifSFD, &
                                 ifBoostConv, ifTDF, ifDMT, ifFloquet, isNewtonFP, &
                                 isNewtonPO, isNewtonPO_T, isDirect, isAdjoint, &
                                 isTransientGrowth, isFloquetDirect, &
                                 isFloquetAdjoint, isFloquetTransientGrowth, &
                                 ifEnergyBudget, ifWavemaker, ifBFSensitivity, &
                                 ifForceSensReal, ifForceSensImag, ifDeltaForcing, &
                                 ifAnimateMode, ifAnimateBFDeform, ifAnimateFloquet, &
                                 ifotd, ifpod, ifdmd, ifspod
   use nekstab_mode_codes, only: MODE_DNS, MODE_LINDNS, MODE_SFD, MODE_BOOSTCONV, &
                                 MODE_DMT, MODE_TDF, MODE_NEWTON_FP, MODE_NEWTON_PO, &
                                 MODE_NEWTON_POT, MODE_DIRECT, MODE_FLOQUET_DIRECT, &
                                 MODE_ADJOINT, MODE_FLOQUET_ADJOINT, MODE_TG, &
                                 MODE_FLOQUET_TG, MODE_ALL_POST, MODE_ENERGY_BUDGET, &
                                 MODE_ENERGY_BUD_FLQ, MODE_WAVEMAKER, &
                                 MODE_BF_SENSITIVITY, MODE_FORCE_SENS_REAL, &
                                 MODE_FORCE_SENS_IMAG, MODE_DELTA_FORCING, &
                                 MODE_ANIMATE, MODE_ANIMATE_DEFORM, MODE_ANIMATE_FLOQUET, &
                                 MODE_OTD, MODE_MODAL_ALL, MODE_POD, MODE_DMD, &
                                 MODE_SPOD, uparam_to_mode_code
   implicit none
   private
   public :: nekStab_resolve_mode, &
             nekStab_mode_from_string, &
             nekStab_mode_from_flags, &
             nekStab_mode_from_uparam, &
             nekStab_validate_mode, &
             nekStab_sync_uparam
contains

!-----------------------------------------------------------------------
! nekStab_resolve_mode — Determine operating mode from all sources
!
! Purpose:
!   Checks three sources in priority order:
!     1. nekstab_mode string (highest) — human readable
!     2. Individual if* flags — flexible
!     3. uparam(1) decoding (lowest) — backward compatible
!   Called in nekStab_init AFTER nekStab_usrchk.
!-----------------------------------------------------------------------
   subroutine nekStab_resolve_mode
      logical :: any_mode_flag_set

      !  Check if user explicitly set any mode flag in nekStab_usrchk
      any_mode_flag_set = ifDNS .or. ifLinDNS .or. &
                          ifSFD .or. ifBoostConv .or. ifTDF .or. &
                          isNewtonFP .or. isNewtonPO .or. isNewtonPO_T .or. &
                          isDirect .or. isAdjoint .or. isTransientGrowth .or. &
                          isFloquetDirect .or. isFloquetAdjoint .or. &
                          isFloquetTransientGrowth .or. &
                          ifEnergyBudget .or. ifWavemaker .or. ifBFSensitivity .or. &
                          ifForceSensReal .or. ifForceSensImag .or. ifDeltaForcing &
                          .or. ifAnimateMode .or. ifAnimateBFDeform &
                          .or. ifAnimateFloquet .or. ifotd &
                          .or. ifDMT &
                          .or. ifpod .or. ifdmd .or. ifspod

      !  Priority 1: String mode (nekstab_mode) - highest priority
      if (len_trim(nekstab_mode) > 0) then
         call nekStab_clear_mode_flags
         call nekStab_mode_from_string(nekstab_mode)
         call nekStab_log('Mode set via nekstab_mode = '//trim(nekstab_mode))

         !  Priority 2: Flag mode - user set explicit flags
      elseif (any_mode_flag_set) then
         call nekStab_mode_from_flags
         call nekStab_log('Mode set via if-flags')

         !  Priority 3: uparam(1) decoding - backward compatible default
      else
         call nekStab_mode_from_uparam
         call nekStab_log('Mode set via uparam(1)')
      end if

      !  Validate: check for conflicting modes
      call nekStab_validate_mode

      !  Sync uparam(1) to match resolved flags (downstream code reads it)
      call nekStab_sync_uparam

      !  Post-resolve exhaustive loud validation (single gate for illegal combinations).
      !  Covers cases the basic nmodes + pairwise checks miss.
      !  Always emits exactly one clear ERROR + nek_end on violation.
      if (ifbuoyancy .and. (thermal_buoyancy_coeff <= 0.0d0) .and. &
          (isDirect .or. isAdjoint .or. ifEnergyBudget .or. isTransientGrowth)) then
         call nekStab_error('thermal stability/energy_budget with ifbuoyancy but thermal_buoyancy_coeff <=0')
      end if

   end subroutine nekStab_resolve_mode

!-----------------------------------------------------------------------
! nekStab_clear_mode_flags — Clear user mode flags before string override
!
! Purpose:
!   String mode has documented priority over flags and uparam(1).  A legacy
!   .usr may still set if* flags before assigning nekstab_mode; clear only the
!   mode-selection flags so the string selector is a true override.
!-----------------------------------------------------------------------
   subroutine nekStab_clear_mode_flags

      ifDNS = .false.
      ifLinDNS = .false.
      ifSFD = .false.
      ifBoostConv = .false.
      ifTDF = .false.
      ifDMT = .false.
      ifFloquet = .false.

      isNewtonFP = .false.
      isNewtonPO = .false.
      isNewtonPO_T = .false.

      isDirect = .false.
      isAdjoint = .false.
      isTransientGrowth = .false.
      isFloquetDirect = .false.
      isFloquetAdjoint = .false.
      isFloquetTransientGrowth = .false.

      ifEnergyBudget = .false.
      ifWavemaker = .false.
      ifBFSensitivity = .false.
      ifForceSensReal = .false.
      ifForceSensImag = .false.
      ifDeltaForcing = .false.
      ifAnimateMode = .false.
      ifAnimateBFDeform = .false.
      ifAnimateFloquet = .false.

      ifotd = .false.
      ifpod = .false.
      ifdmd = .false.
      ifspod = .false.

   end subroutine nekStab_clear_mode_flags

!-----------------------------------------------------------------------
! nekStab_mode_from_string — Parse mode string into flags
!
! Purpose:
!   Converts a human-readable mode string (e.g., 'direct',
!   'floquet_adjoint') into the corresponding boolean flags.
!   Case-insensitive.
!
! Arguments:
!   mode_str [in] — mode string to parse
!-----------------------------------------------------------------------
   subroutine nekStab_mode_from_string(mode_str)
      character(len=*), intent(in) :: mode_str
      character(len=32) :: mode_lower
      integer :: i

      !  Convert to lowercase for case-insensitive comparison
      mode_lower = adjustl(mode_str)
      do i = 1, len_trim(mode_lower)
         if (mode_lower(i:i) >= 'A' .and. mode_lower(i:i) <= 'Z') then
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

      case ('dmt')
         ifDMT = .true.

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
      case ('energy_budget_floquet')
         ifEnergyBudget = .true.
         ifFloquet = .true.
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
         call nekStab_error('Unknown nekstab_mode: '//trim(mode_str))
      end select

   end subroutine nekStab_mode_from_string

!-----------------------------------------------------------------------
! nekStab_mode_from_flags — Apply ifFloquet modifier flag
!
! Purpose:
!   Transforms base stability mode + ifFloquet modifier into
!   the combined Floquet flags. E.g., isDirect + ifFloquet
!   becomes isFloquetDirect.
!-----------------------------------------------------------------------
   subroutine nekStab_mode_from_flags

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

!-----------------------------------------------------------------------
! nekStab_mode_from_uparam — Decode uparam(1) into mode flags
!
! Purpose:
!   Backward-compatible decoder that maps uparam(1) floating-point
!   values to boolean mode flags.
!   Integer select case on nint(uparam(1)*100) using MODE_* constants from
!   nekstab_mode_codes. A value that is not a mode code stops the run with an
!   error; it never falls back to DNS.
!-----------------------------------------------------------------------
   subroutine nekStab_mode_from_uparam
      real :: up1
      real, parameter :: tol = 1.0e-4
      integer :: icode
      character(len=80) :: msg

      up1 = uparam(1)
      icode = uparam_to_mode_code(real(up1, 8))

      !  nint() alone would read 3.104 as 3.1. Accept only a value within tol of a
      !  code (the tolerance the former float fallback used), so a typo cannot
      !  select a neighbouring mode.
      if (abs(real(up1, 8) - 0.01d0*icode) > tol) then
         write (msg, '(a,f9.4)') 'Unknown mode, userParam01 = ', up1
         call nekStab_error(trim(msg)//'. It is not a mode code (see the list in the docs).')
      end if

      !  Primary decoder: integer select case on nint(uparam(1)*100).
      !  Unambiguous and collision-free (see docs/mode-code-constants-design.md).
      select case (icode)

         !  Mode 0: DNS
      case (MODE_DNS)
         ifDNS = .true.
      case (MODE_LINDNS)
         ifLinDNS = .true.

         !  Mode 1: Fixed point methods
      case (MODE_SFD)
         ifSFD = .true.
      case (MODE_BOOSTCONV)
         ifBoostConv = .true.
      case (MODE_DMT)
         ifDMT = .true.
      case (MODE_TDF)
         ifTDF = .true.

         !  Mode 2: Newton-Krylov
      case (MODE_NEWTON_FP)
         isNewtonFP = .true.
      case (MODE_NEWTON_PO)
         isNewtonPO = .true.
      case (MODE_NEWTON_POT)
         isNewtonPO_T = .true.

         !  Mode 3: Stability analysis
      case (MODE_DIRECT)
         isDirect = .true.
      case (MODE_FLOQUET_DIRECT)
         isFloquetDirect = .true.
      case (MODE_ADJOINT)
         isAdjoint = .true.
      case (MODE_FLOQUET_ADJOINT)
         isFloquetAdjoint = .true.
      case (MODE_TG)
         isTransientGrowth = .true.
      case (MODE_FLOQUET_TG)
         isFloquetTransientGrowth = .true.

         !  Mode 4: Post-processing
      case (MODE_ALL_POST)
         ifEnergyBudget = .true.
         ifWavemaker = .true.
         ifBFSensitivity = .true.
      case (MODE_ENERGY_BUDGET)
         ifEnergyBudget = .true.
      case (MODE_ENERGY_BUD_FLQ)
         ifEnergyBudget = .true.
         ifFloquet = .true.
      case (MODE_WAVEMAKER)
         ifWavemaker = .true.
      case (MODE_BF_SENSITIVITY)
         ifBFSensitivity = .true.
      case (MODE_FORCE_SENS_REAL)
         ifForceSensReal = .true.
      case (MODE_FORCE_SENS_IMAG)
         ifForceSensImag = .true.
      case (MODE_DELTA_FORCING)
         ifDeltaForcing = .true.
      case (MODE_ANIMATE)
         ifAnimateMode = .true.
         animate_mode_num = int(uparam(7))
      case (MODE_ANIMATE_DEFORM)
         ifAnimateBFDeform = .true.
         animate_mode_num = int(uparam(7))
      case (MODE_ANIMATE_FLOQUET)
         ifAnimateFloquet = .true.
         animate_mode_num = int(uparam(7))

         !  Mode 5: OTD
      case (MODE_OTD)
         ifotd = .true.

         !  Mode 6: Modal analysis
      case (MODE_MODAL_ALL)
         ifpod = .true.
         ifdmd = .true.
         ifspod = .true.
      case (MODE_POD)
         ifpod = .true.
      case (MODE_DMD)
         ifdmd = .true.
      case (MODE_SPOD)
         ifspod = .true.

      case default
         !  An unknown code used to leave every flag false, and nekStab_validate_mode
         !  then fell back to DNS with a warning: a typo such as userParam01 = 3.15
         !  ran a long DNS instead of the stability analysis that was meant.
         write (msg, '(a,f9.4)') 'Unknown mode, userParam01 = ', up1
         call nekStab_error(trim(msg)//'. Valid: 0, 0.1, 1.1-1.4, 2.0-2.2, '// &
                            '3.1, 3.11, 3.2, 3.21, 3.3, 3.31, 4, 4.1-4.3, 4.11, '// &
                            '4.41-4.43, 4.5-4.52, 5, 6, 6.1-6.3')
      end select

   end subroutine nekStab_mode_from_uparam

!-----------------------------------------------------------------------
! nekStab_validate_mode — Check for conflicting mode selections
!
! Purpose:
!   Ensures only one main mode category is active and checks for
!   conflicting sub-modes (e.g., direct + adjoint). Falls back
!   to DNS if no mode is selected.
!-----------------------------------------------------------------------
   subroutine nekStab_validate_mode
      integer :: nmodes

      nmodes = 0

      !  Count active mode categories
      if (ifDNS .or. ifLinDNS) nmodes = nmodes + 1
      if (ifSFD .or. ifBoostConv .or. ifTDF .or. ifDMT) nmodes = nmodes + 1
      if (isNewtonFP .or. isNewtonPO .or. isNewtonPO_T) &
         nmodes = nmodes + 1
      if (isDirect .or. isFloquetDirect .or. &
          isAdjoint .or. isFloquetAdjoint .or. &
          isTransientGrowth .or. isFloquetTransientGrowth) &
         nmodes = nmodes + 1
      if (ifEnergyBudget .or. ifWavemaker .or. ifBFSensitivity .or. &
          ifForceSensReal .or. ifForceSensImag .or. ifDeltaForcing &
          .or. ifAnimateMode .or. ifAnimateBFDeform &
          .or. ifAnimateFloquet) nmodes = nmodes + 1
      if (ifotd) nmodes = nmodes + 1
      if (ifpod .or. ifdmd .or. ifspod) nmodes = nmodes + 1

      if (nmodes > 1) then
         call nekStab_error('Multiple conflicting modes active - set only ONE mode category')
      end if

      !  Check for conflicting stability sub-modes
      if ((isDirect .or. isFloquetDirect) .and. &
          (isAdjoint .or. isFloquetAdjoint)) then
         call nekStab_error('Both direct and adjoint modes set - choose one')
      end if

      if ((isDirect .or. isFloquetDirect .or. &
           isAdjoint .or. isFloquetAdjoint) .and. &
          (isTransientGrowth .or. isFloquetTransientGrowth)) then
         call nekStab_error('Eigenmode and transient growth both set - choose one')
      end if

      !  Check for orphaned ifFloquet
      if (ifFloquet .and. .not. (isFloquetDirect .or. &
                                 isFloquetAdjoint .or. isFloquetTransientGrowth .or. &
                                 ifEnergyBudget)) then
         call nekStab_error('ifFloquet set without base mode (direct/adjoint/TG)')
      end if

      !  Check for no mode selected
      if (nmodes == 0) then
         if (nid == 0) then
            write (6, *) 'WARNING: No operating mode selected'
            write (6, *) 'Defaulting to DNS mode (uparam(1)=0)'
         end if
         ifDNS = .true.
      end if

   end subroutine nekStab_validate_mode

!-----------------------------------------------------------------------
! nekStab_sync_uparam — Write resolved flags back to uparam(1)
!
! Purpose:
!   Ensures uparam(1) matches the resolved mode flags so that
!   downstream code dispatching on uparam(1) works correctly
!   regardless of which mode-selection method was used.
!-----------------------------------------------------------------------
   subroutine nekStab_sync_uparam

      !  Map resolved flags back to uparam(1)
      !  Mode 0: DNS
      if (ifDNS) then
         uparam(1) = 0.0
      elseif (ifLinDNS) then
         uparam(1) = 0.1

         !  Mode 1: Fixed point
      elseif (ifSFD) then
         uparam(1) = 1.1
      elseif (ifBoostConv) then
         uparam(1) = 1.2
      elseif (ifDMT) then
         uparam(1) = 1.3
      elseif (ifTDF) then
         uparam(1) = 1.4

         !  Mode 2: Newton-Krylov
      elseif (isNewtonFP) then
         uparam(1) = 2.0
      elseif (isNewtonPO) then
         uparam(1) = 2.1
      elseif (isNewtonPO_T) then
         uparam(1) = 2.2

         !  Mode 3: Stability analysis
      elseif (isDirect) then
         uparam(1) = 3.1
      elseif (isFloquetDirect) then
         uparam(1) = 3.11
      elseif (isAdjoint) then
         uparam(1) = 3.2
      elseif (isFloquetAdjoint) then
         uparam(1) = 3.21
      elseif (isTransientGrowth) then
         uparam(1) = 3.3
      elseif (isFloquetTransientGrowth) then
         uparam(1) = 3.31

         !  Mode 4: Postprocessing
      elseif (ifEnergyBudget .and. ifFloquet) then
         uparam(1) = 4.11
      elseif (ifEnergyBudget) then
         uparam(1) = 4.1
      elseif (ifWavemaker) then
         uparam(1) = 4.2
      elseif (ifBFSensitivity) then
         uparam(1) = 4.3
      elseif (ifForceSensReal) then
         uparam(1) = 4.41
      elseif (ifForceSensImag) then
         uparam(1) = 4.42
      elseif (ifDeltaForcing) then
         uparam(1) = 4.43
      elseif (ifAnimateMode) then
         uparam(1) = 4.50
      elseif (ifAnimateBFDeform) then
         uparam(1) = 4.51
      elseif (ifAnimateFloquet) then
         uparam(1) = 4.52

         !  Mode 5: OTD
      elseif (ifotd) then
         uparam(1) = 5.0

         !  Mode 6: Modal analysis
      elseif (ifpod) then
         uparam(1) = 6.1
      elseif (ifdmd) then
         uparam(1) = 6.2
      elseif (ifspod) then
         uparam(1) = 6.3
      end if

   end subroutine nekStab_sync_uparam

end module nekstab_mode_config
