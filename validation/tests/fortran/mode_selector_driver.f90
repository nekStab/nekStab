! Resolve one operating mode and print it, for test_mode_selector.py.
!
! Usage: mode_selector_driver <uparam1> [mode-string] [flag ...]
!   uparam1      value of userParam01 as a decimal number
!   mode-string  value of nekstab_mode; '-' leaves it empty
!   flag         name of a mode flag set before resolution (for example isDirect)
! Output, one line: the names of the flags that are true, then the value of uparam(1)
! after resolution. An error stops the program with status 3.
program mode_selector_driver
   use nekstab_nek_bridge
   use nekstab_mode_config, only: nekStab_resolve_mode
   implicit none
   character(len=64) :: arg
   character(len=512) :: line
   integer :: i

   call get_command_argument(1, arg)
   read (arg, *) uparam(1)
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      if (trim(arg) /= '-') nekstab_mode = trim(arg)
   end if
   do i = 3, command_argument_count()
      call get_command_argument(i, arg)
      call set_flag(trim(arg))
   end do

   call nekStab_resolve_mode

   line = ''
   call add(ifDNS, 'ifDNS'); call add(ifLinDNS, 'ifLinDNS')
   call add(ifSFD, 'ifSFD'); call add(ifBoostConv, 'ifBoostConv')
   call add(ifTDF, 'ifTDF'); call add(ifDMT, 'ifDMT'); call add(ifFloquet, 'ifFloquet')
   call add(isNewtonFP, 'isNewtonFP'); call add(isNewtonPO, 'isNewtonPO')
   call add(isNewtonPO_T, 'isNewtonPO_T')
   call add(isDirect, 'isDirect'); call add(isAdjoint, 'isAdjoint')
   call add(isTransientGrowth, 'isTransientGrowth')
   call add(isFloquetDirect, 'isFloquetDirect'); call add(isFloquetAdjoint, 'isFloquetAdjoint')
   call add(isFloquetTransientGrowth, 'isFloquetTransientGrowth')
   call add(ifEnergyBudget, 'ifEnergyBudget'); call add(ifWavemaker, 'ifWavemaker')
   call add(ifBFSensitivity, 'ifBFSensitivity'); call add(ifForceSensReal, 'ifForceSensReal')
   call add(ifForceSensImag, 'ifForceSensImag'); call add(ifDeltaForcing, 'ifDeltaForcing')
   call add(ifAnimateMode, 'ifAnimateMode'); call add(ifAnimateBFDeform, 'ifAnimateBFDeform')
   call add(ifAnimateFloquet, 'ifAnimateFloquet')
   call add(ifotd, 'ifotd'); call add(ifpod, 'ifpod'); call add(ifdmd, 'ifdmd')
   call add(ifspod, 'ifspod')
   write (6, '(a,f7.3)') trim(line)//' | uparam1=', uparam(1)

contains

   subroutine add(flag, name)
      logical, intent(in) :: flag
      character(len=*), intent(in) :: name
      if (flag) line = trim(line)//' '//name
   end subroutine add

   subroutine set_flag(name)
      character(len=*), intent(in) :: name
      select case (name)
      case ('ifDNS'); ifDNS = .true.
      case ('ifSFD'); ifSFD = .true.
      case ('isNewtonFP'); isNewtonFP = .true.
      case ('isDirect'); isDirect = .true.
      case ('isAdjoint'); isAdjoint = .true.
      case ('isTransientGrowth'); isTransientGrowth = .true.
      case ('ifFloquet'); ifFloquet = .true.
      case ('ifEnergyBudget'); ifEnergyBudget = .true.
      case ('ifWavemaker'); ifWavemaker = .true.
      case default
         write (6, '(a)') 'driver: unknown flag '//name
         stop 4
      end select
   end subroutine set_flag

end program mode_selector_driver
