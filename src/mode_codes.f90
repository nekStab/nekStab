!-----------------------------------------------------------------------
! mode_codes.f90 — Integer mode-code constants for nekStab
!
! Purpose:
!   Named integer parameters for every nekStab operating mode.
!   Each constant equals nint(uparam(1) * 100.0), providing an
!   unambiguous, collision-free identifier for each mode.
!
!   Phase 1 of the mode-code cleanup plan.
!   Used internally by nekStab_mode_from_uparam in mode_config.f90;
!   no external .par or .usr changes required in this phase.
!
! Public interface:
!   MODE_* constants — integer codes for every operating mode
!   uparam_to_mode_code — real(8) -> integer conversion
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_mode_codes
   implicit none
   private

   public :: MODE_DNS, MODE_LINDNS, MODE_SFD, MODE_BOOSTCONV, MODE_DMT, &
             MODE_TDF, MODE_NEWTON_FP, MODE_NEWTON_PO, MODE_NEWTON_POT, &
             MODE_DIRECT, MODE_FLOQUET_DIRECT, MODE_ADJOINT, &
             MODE_FLOQUET_ADJOINT, MODE_TG, MODE_FLOQUET_TG, MODE_ALL_POST, &
             MODE_ENERGY_BUDGET, MODE_ENERGY_BUD_FLQ, MODE_WAVEMAKER, &
             MODE_BF_SENSITIVITY, MODE_FORCE_SENS_REAL, MODE_FORCE_SENS_IMAG, &
             MODE_DELTA_FORCING, MODE_ANIMATE, MODE_ANIMATE_DEFORM, &
             MODE_ANIMATE_FLOQUET, MODE_OTD, MODE_MODAL_ALL, MODE_POD, &
             MODE_DMD, MODE_SPOD, uparam_to_mode_code

   !  Integer mode codes: nint(uparam(1) * 100.0)
   !  Phase 1 of mode-code cleanup design.

   !  Mode 0: DNS
   integer, parameter :: MODE_DNS = 0
   integer, parameter :: MODE_LINDNS = 10

   !  Mode 1: Fixed-point convergence
   integer, parameter :: MODE_SFD = 110
   integer, parameter :: MODE_BOOSTCONV = 120
   integer, parameter :: MODE_DMT = 130
   integer, parameter :: MODE_TDF = 140

   !  Mode 2: Newton-Krylov
   integer, parameter :: MODE_NEWTON_FP = 200
   integer, parameter :: MODE_NEWTON_PO = 210
   integer, parameter :: MODE_NEWTON_POT = 220

   !  Mode 3: Stability analysis
   integer, parameter :: MODE_DIRECT = 310
   integer, parameter :: MODE_FLOQUET_DIRECT = 311
   integer, parameter :: MODE_ADJOINT = 320
   integer, parameter :: MODE_FLOQUET_ADJOINT = 321
   integer, parameter :: MODE_TG = 330
   integer, parameter :: MODE_FLOQUET_TG = 331

   !  Mode 4: Post-processing
   integer, parameter :: MODE_ALL_POST = 400
   integer, parameter :: MODE_ENERGY_BUDGET = 410
   integer, parameter :: MODE_ENERGY_BUD_FLQ = 411
   integer, parameter :: MODE_WAVEMAKER = 420
   integer, parameter :: MODE_BF_SENSITIVITY = 430
   integer, parameter :: MODE_FORCE_SENS_REAL = 441
   integer, parameter :: MODE_FORCE_SENS_IMAG = 442
   integer, parameter :: MODE_DELTA_FORCING = 443
   integer, parameter :: MODE_ANIMATE = 450
   integer, parameter :: MODE_ANIMATE_DEFORM = 451
   integer, parameter :: MODE_ANIMATE_FLOQUET = 452

   !  Mode 5: OTD
   integer, parameter :: MODE_OTD = 500

   !  Mode 6: Modal analysis
   integer, parameter :: MODE_MODAL_ALL = 600
   integer, parameter :: MODE_POD = 610
   integer, parameter :: MODE_DMD = 620
   integer, parameter :: MODE_SPOD = 630

contains

   integer function uparam_to_mode_code(up1)
      !  Convert uparam(1) float to integer mode code.
      !  Uses double precision to avoid rounding errors for
      !  values such as 3.11 and 4.11 (two decimal places).
      real(8), intent(in) :: up1
      uparam_to_mode_code = nint(up1*100.0d0)
   end function uparam_to_mode_code

end module nekstab_mode_codes
