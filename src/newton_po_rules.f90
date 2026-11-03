!-----------------------------------------------------------------------
! newton_po_rules.f90 — Step count and warnings for a Newton periodic orbit
!
! Purpose:
!   dt = T/N for one period, the CFL warning and the explicit-filter warning.
!   The routines have no dependency on Nek5000.
!   A stand-alone test can compile them.
!   Newton and the linearized solver call them.
!
! Public interface:
!   newton_po_set_step          — dt and the step count for one period
!   newton_po_cfl_warn          — CFL number above 1.5 times the target
!   newton_po_filter_warn       — explicit filter on a periodic-orbit Newton run
!   newton_po_kept_steps        — step count kept from iteration 1
!   newton_po_steps_after_stage — kept count after the Newton stage
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_newton_po_rules
   implicit none
   private
   public :: newton_po_set_step, newton_po_cfl_warn, newton_po_filter_warn, &
             newton_po_kept_steps, newton_po_steps_after_stage

contains

!-----------------------------------------------------------------------
! newton_po_set_step — Set dt and the step count for one period
!
!   n_kept > 0 gives nsteps = n_kept and dt = period/n_kept; the value of
!   dt on entry is not used in that case.
!   n_kept = 0 keeps the provisional step dt on entry.
!   The count is ceiling(period/dt). The final step is period/count.
!-----------------------------------------------------------------------
   pure subroutine newton_po_set_step(period, n_kept, dt, nsteps)
      real, intent(in) :: period
      integer, intent(in) :: n_kept
      real, intent(inout) :: dt
      integer, intent(out) :: nsteps

      if (n_kept > 0) then
         ! Not ceiling(period/(period/n_kept)): rounding can give n_kept + 1,
         ! and Newton keeps n_kept (period 1.0 and count 49 give 50).
         nsteps = n_kept
      else
         nsteps = ceiling(period/dt)
      end if
      dt = period/nsteps

   end subroutine newton_po_set_step

!-----------------------------------------------------------------------
! newton_po_cfl_warn — Return true when the CFL number is above the limit
!
!   The kept count must be set and the target must be positive.
!   A CFL number equal to 1.5 times the target does not warn.
!-----------------------------------------------------------------------
   pure logical function newton_po_cfl_warn(n_kept, target_cfl, cfl_now) result(warn)
      integer, intent(in) :: n_kept
      real, intent(in) :: target_cfl, cfl_now

      warn = n_kept > 0 .and. target_cfl > 0.0d0 .and. cfl_now > 1.5d0*target_cfl

   end function newton_po_cfl_warn

!-----------------------------------------------------------------------
! newton_po_filter_warn — Return true for an explicit filter on an orbit
!
!   Nek stores filterType as a real. The value 1 is the explicit filter.
!   A positive weight is required.
!   The relaxation filter is filter_type 2 and does not warn.
!-----------------------------------------------------------------------
   pure logical function newton_po_filter_warn(is_po, is_po_t, filter_type, weight) result(warn)
      logical, intent(in) :: is_po, is_po_t
      real, intent(in) :: filter_type, weight

      warn = (is_po .or. is_po_t) .and. filter_type == 1 .and. weight > 0.0d0

   end function newton_po_filter_warn

!-----------------------------------------------------------------------
! newton_po_kept_steps — Return the step count stored for the next iteration
!
!   Iteration 1 of a periodic orbit stores nsteps.
!   Every other call returns the stored count.
!-----------------------------------------------------------------------
   pure integer function newton_po_kept_steps(iteration, is_po, is_po_t, nsteps, stored) &
      result(kept)
      integer, intent(in) :: iteration, nsteps, stored
      logical, intent(in) :: is_po, is_po_t

      if (iteration == 1 .and. (is_po .or. is_po_t)) then
         kept = nsteps
      else
         kept = stored
      end if

   end function newton_po_kept_steps

!-----------------------------------------------------------------------
! newton_po_steps_after_stage — Return the step count after the Newton stage
!
!   The result is zero. The next stage chooses its own step count.
!-----------------------------------------------------------------------
   pure integer function newton_po_steps_after_stage() result(kept)
      kept = 0
   end function newton_po_steps_after_stage

end module nekstab_newton_po_rules
