! Stand-in for nekstab_nek_bridge, used only by test_mode_selector.py.
! It declares the mode flags, uparam and nekstab_mode that src/mode_config.f90
! imports, so mode resolution can run without a Nek5000 build.
module nekstab_nek_bridge
   implicit none
   public
   integer :: nid = 0
   integer :: animate_mode_num = 0
   real(8) :: uparam(20) = 0.0d0
   real(8) :: thermal_buoyancy_coeff = 0.0d0
   character(len=32) :: nekstab_mode = ''
   logical :: ifbuoyancy = .false.
   logical :: ifDNS = .false., ifLinDNS = .false., ifSFD = .false.
   logical :: ifBoostConv = .false., ifTDF = .false., ifDMT = .false.
   logical :: ifFloquet = .false.
   logical :: isNewtonFP = .false., isNewtonPO = .false., isNewtonPO_T = .false.
   logical :: isDirect = .false., isAdjoint = .false., isTransientGrowth = .false.
   logical :: isFloquetDirect = .false., isFloquetAdjoint = .false.
   logical :: isFloquetTransientGrowth = .false.
   logical :: ifEnergyBudget = .false., ifWavemaker = .false., ifBFSensitivity = .false.
   logical :: ifForceSensReal = .false., ifForceSensImag = .false., ifDeltaForcing = .false.
   logical :: ifAnimateMode = .false., ifAnimateBFDeform = .false., ifAnimateFloquet = .false.
   logical :: ifotd = .false., ifpod = .false., ifdmd = .false., ifspod = .false.
contains
   subroutine nekStab_error(msg)
      character(len=*), intent(in) :: msg
      write (6, '(a)') 'ERROR: '//trim(msg)
      stop 3
   end subroutine nekStab_error

   subroutine nekStab_log(msg)
      character(len=*), intent(in) :: msg
      ! Silent: the driver prints only the resolved mode.
      if (len(msg) < 0) stop 9
   end subroutine nekStab_log
end module nekstab_nek_bridge
