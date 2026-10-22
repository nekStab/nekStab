      !-----------------------------------------------------------------------
      ! usr_wrappers.f90 — Wrappers for user .f files to call module procedures
      !
      ! Purpose:
      !   Provides non-module subroutines that user .f files can call
      !   without use statements. They delegate to the actual module
      !   procedures (historical Nek5000 .usr compatibility).
      !
      ! Public interface:
      !   nekStab_torque          — wrapper for torque/drag
      !   nekStab_forcing         — wrapper for velocity forcing
      !   nekStab_forcing_temp    — wrapper for temperature forcing
      !   nekStab_define_obj      — wrapper for object definition
      !
      ! Dependencies:
      !   (delegates to nekstab_torque_mod, nekstab_forcing_mod)
      !-----------------------------------------------------------------------

      subroutine nekStab_torque(fname)
         use nekstab_torque_mod, only: mod_nekStab_torque => nekStab_torque
         implicit none
         character(len=*), intent(in) :: fname
         call mod_nekStab_torque(fname)

      end subroutine nekStab_torque

      subroutine nekStab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg)
         use nekstab_forcing_mod, only: mod_nekStab_forcing => nekStab_forcing
         implicit none
         real, intent(inout) :: ffx, ffy, ffz
         integer, intent(in) :: ix, iy, iz, ieg
         call mod_nekStab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg)

      end subroutine nekStab_forcing

      subroutine nekStab_forcing_temp(temp, ix, iy, iz, ieg, m)
         use nekstab_forcing_mod,
     &      only: mod_nekStab_forcing_temp => nekStab_forcing_temp
         implicit none
         real, intent(inout) :: temp
         integer, intent(in) :: ix, iy, iz, ieg, m
         call mod_nekStab_forcing_temp(temp, ix, iy, iz, ieg, m)

      end subroutine nekStab_forcing_temp

      subroutine nekStab_define_obj()
         use nekstab_torque_mod, only: mod_nekStab_define_obj => nekStab_define_obj
         implicit none
         call mod_nekStab_define_obj()

      end subroutine nekStab_define_obj
