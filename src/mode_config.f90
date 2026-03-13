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
!-----------------------------------------------------------------------

module nekstab_mode_config
   use nekstab_nek_bridge
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
      .or. ifpod .or. ifdmd .or. ifspod

   !  Priority 1: String mode (nekstab_mode) - highest priority
   if (len_trim(nekstab_mode) > 0) then
      call nekStab_mode_from_string(nekstab_mode)
      if (nid == 0) write (6, *) 'Mode set via nekstab_mode = ', &
         trim(nekstab_mode)

   !  Priority 2: Flag mode - user set explicit flags
   elseif (any_mode_flag_set) then
      call nekStab_mode_from_flags
      if (nid == 0) write (6, *) 'Mode set via if-flags'

   !  Priority 3: uparam(1) decoding - backward compatible default
   else
      call nekStab_mode_from_uparam
      if (nid == 0) write (6, *) 'Mode set via uparam(1) =', &
         uparam(1)
   end if

   !  Validate: check for conflicting modes
   call nekStab_validate_mode

   !  Sync uparam(1) to match resolved flags (downstream code reads it)
   call nekStab_sync_uparam

end subroutine nekStab_resolve_mode

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
      if (nid == 0) then
         write (6, *) 'ERROR: Unknown nekstab_mode: ', &
            trim(mode_str)
         write (6, *) 'Valid modes: dns, linear_dns, sfd, ', &
            'boostconv, tdf,'
         write (6, *) '  newton_fp, newton_po, newton_po_t,'
         write (6, *) '  direct, adjoint, transient_growth,'
         write (6, *) '  floquet_direct, floquet_adjoint, ', &
            'floquet_tg,'
         write (6, *) '  energy_budget, ', &
            'energy_budget_floquet,'
         write (6, *) '  wavemaker, ', &
            'bf_sensitivity,'
         write (6, *) '  animate_mode, otd, pod, dmd, spod'
      end if
      call nek_end
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
!   values to boolean mode flags. Uses tolerance-based comparison.
!-----------------------------------------------------------------------
subroutine nekStab_mode_from_uparam
   real :: up1
   real, parameter :: tol = 1.0e-4

   up1 = uparam(1)

   !  Mode 0: DNS
   if (abs(up1 - 0.0) < tol) then
      ifDNS = .true.
   elseif (abs(up1 - 0.1) < tol) then
      ifLinDNS = .true.

   !  Mode 1: Fixed point methods
   elseif (abs(up1 - 1.1) < tol) then
      ifSFD = .true.
   elseif (abs(up1 - 1.2) < tol) then
      ifBoostConv = .true.
   elseif (abs(up1 - 1.4) < tol) then
      ifTDF = .true.

   !  Mode 2: Newton-Krylov
   elseif (abs(up1 - 2.0) < tol) then
      isNewtonFP = .true.
   elseif (abs(up1 - 2.1) < tol) then
      isNewtonPO = .true.
   elseif (abs(up1 - 2.2) < tol) then
      isNewtonPO_T = .true.

   !  Mode 3: Eigenvalue problem (stability analysis)
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

   !  Mode 4: Postprocessing
   elseif (abs(up1 - 4.0) < tol) then
      ifEnergyBudget = .true.
      ifWavemaker = .true.
      ifBFSensitivity = .true.
   elseif (abs(up1 - 4.1) < tol) then
      ifEnergyBudget = .true.
   elseif (abs(up1 - 4.11) < tol) then
      ifEnergyBudget = .true.
      ifFloquet = .true.
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

   !  Mode 5: OTD
   elseif (abs(up1 - 5.0) < tol) then
      ifotd = .true.

   !  Mode 6: Modal analysis
   elseif (abs(up1 - 6.0) < tol) then
      ifpod = .true.
      ifdmd = .true.
      ifspod = .true.
   elseif (abs(up1 - 6.1) < tol) then
      ifpod = .true.
   elseif (abs(up1 - 6.2) < tol) then
      ifdmd = .true.
   elseif (abs(up1 - 6.3) < tol) then
      ifspod = .true.

   end if

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
   if (ifSFD .or. ifBoostConv .or. ifTDF) nmodes = nmodes + 1
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
      if (nid == 0) then
         write (6, *) 'ERROR: Multiple conflicting modes active'
         write (6, *) 'Set only ONE mode category at a time'
      end if
      call nek_end
   end if

   !  Check for conflicting stability sub-modes
   if ((isDirect .or. isFloquetDirect) .and. &
      (isAdjoint .or. isFloquetAdjoint)) then
      if (nid == 0) then
         write (6, *) 'ERROR: Both direct and adjoint modes set'
         write (6, *) 'Choose one: direct OR adjoint'
      end if
      call nek_end
   end if

   if ((isDirect .or. isFloquetDirect .or. &
      isAdjoint .or. isFloquetAdjoint) .and. &
      (isTransientGrowth .or. isFloquetTransientGrowth)) then
      if (nid == 0) then
         write (6, *) 'ERROR: Eigenmode and transient growth ', &
            'both set'
         write (6, *) 'Choose one: eigenmode OR transient growth'
      end if
      call nek_end
   end if

   !  Check for orphaned ifFloquet
   if (ifFloquet .and. .not. (isFloquetDirect .or. &
      isFloquetAdjoint .or. isFloquetTransientGrowth .or. &
      ifEnergyBudget)) then
      if (nid == 0) then
         write (6, *) 'ERROR: ifFloquet set without base mode'
         write (6, *) 'Set isDirect, isAdjoint, or ', &
            'isTransientGrowth with ifFloquet'
      end if
      call nek_end
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
