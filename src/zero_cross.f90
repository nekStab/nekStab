!-----------------------------------------------------------------------
! zero_cross.f90 — Period of one probe signal from its upward crossings
!
! Purpose:
!   The crossing time, the period and the mean over that period.
!   The routines have no dependency on Nek5000, so a stand-alone
!   test can compile them. The probe loop in probes.f90 calls them.
!
! Public interface:
!   crossing_state  — samples, mean and crossing times kept between steps
!   crossing_reset  — set the state at the first step
!   crossing_update — advance the state by one sample
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_zero_cross
   implicit none
   private
   public :: crossing_state, crossing_reset, crossing_update

   type :: crossing_state
      real :: v_mean = 0.0d0
      real :: v_sum = 0.0d0
      real :: t_cross = 0.0d0
      real :: t_cross_old = 0.0d0
      real :: p_old = 0.0d0
      integer :: n_cross = 0
      logical :: restart = .false.
      logical :: running = .false.
      integer :: n_samples = 1
   end type crossing_state

contains

!-----------------------------------------------------------------------
! crossing_theta — Fraction of the step at which the signal crosses
!
!   v_prev is below the level and v_now is at or above it.
!   The fraction is clamped to [0, 1].
!-----------------------------------------------------------------------
   pure function crossing_theta(v_prev, v_now, v_level) result(theta)
      real, intent(in) :: v_prev, v_now, v_level
      real :: theta

      theta = (v_level - v_prev)/(v_now - v_prev)
      theta = max(0.0d0, min(1.0d0, theta))

   end function crossing_theta

!-----------------------------------------------------------------------
! crossing_time — Time of the crossing on the current step
!
!   The signal is taken linear between t_now-dt and t_now.
!-----------------------------------------------------------------------
   pure function crossing_time(t_now, dt, theta) result(t_hit)
      real, intent(in) :: t_now, dt, theta
      real :: t_hit

      t_hit = (t_now - dt) + theta*dt

   end function crossing_time

!-----------------------------------------------------------------------
! crossing_reset — Set the state at the first step
!
!   v_mean_init is the level used until the first period closes. With running,
!   it is the first sample and the level follows the mean of the samples until the
!   first crossing: a fixed first sample can be an extremum that the signal never
!   rises through again.
!-----------------------------------------------------------------------
   pure subroutine crossing_reset(state, v_mean_init, running)
      type(crossing_state), intent(out) :: state
      real, intent(in) :: v_mean_init
      logical, intent(in), optional :: running

      state%v_mean = v_mean_init
      state%v_sum = 0.0d0
      state%t_cross = 0.0d0
      state%t_cross_old = 0.0d0
      state%p_old = 0.0d0
      state%n_cross = 0
      state%restart = .false.
      state%running = .false.
      if (present(running)) state%running = running
      state%n_samples = 1

   end subroutine crossing_reset

!-----------------------------------------------------------------------
! crossing_update — Advance the period estimate by one sample
!
!   An upward crossing is v_prev below the level and v_now at or
!   above it, with v_now above v_prev. active is false on the first
!   step, where the solver only stores the sample.
!   phase is 0 with no crossing, 1 at the first crossing, 2 when
!   a period closes and 3 when the start of the period moves to a later
!   crossing of the new level. v_logged is the level of that period. The mean
!   stored for the next period is the integral over the one that
!   just closed. period_change is |T - T_old| at the log line.
!-----------------------------------------------------------------------
   pure subroutine crossing_update(state, v_prev, v_now, t_now, dt, active, &
                                   phase, theta, t_hit, period, v_logged, period_change)
      type(crossing_state), intent(inout) :: state
      real, intent(in) :: v_prev, v_now, t_now, dt
      logical, intent(in) :: active
      integer, intent(out) :: phase
      real, intent(out) :: theta, t_hit, period, v_logged, period_change
      logical :: hit

      phase = 0
      theta = 0.0d0
      t_hit = state%t_cross
      period = 0.0d0
      v_logged = state%v_mean
      period_change = 0.0d0

      hit = active .and. v_prev < state%v_mean .and. v_now >= state%v_mean &
            .and. v_now > v_prev
      if (hit .and. state%restart) then
         ! A period closed with a new level that its last step did not cross. This
         ! crossing, on the same edge or one period later, is the first of the new
         ! level. Make it the start of the next period (phase 3) and close nothing.
         ! Testing the time since the last crossing
         ! instead would stop the probe for good when the first gap is much longer
         ! than the period of the settled flow.
         theta = crossing_theta(v_prev, v_now, state%v_mean)
         t_hit = crossing_time(t_now, dt, theta)
         phase = 3
         state%restart = .false.
         state%t_cross = t_hit
         state%v_sum = 0.5d0*(state%v_mean + v_now)*(1.0d0 - theta)*dt
         return
      end if

      if (hit) then
         theta = crossing_theta(v_prev, v_now, state%v_mean)
         t_hit = crossing_time(t_now, dt, theta)
         state%n_cross = state%n_cross + 1

         if (state%n_cross == 1) then
            phase = 1
            state%t_cross = t_hit
            state%t_cross_old = state%t_cross
            state%v_sum = 0.5d0*(state%v_mean + v_now)*(1.0d0 - theta)*dt
         else
            phase = 2
            state%v_sum = state%v_sum + 0.5d0*(v_prev + state%v_mean)*theta*dt
            state%t_cross_old = state%t_cross
            state%t_cross = t_hit
            period = state%t_cross - state%t_cross_old
            if (state%n_cross == 2) state%p_old = period
            period_change = abs(period - state%p_old)
            v_logged = state%v_mean
            if (period > 0.0d0) state%v_mean = state%v_sum/period
            state%p_old = period
            state%v_sum = 0.5d0*(v_logged + v_now)*(1.0d0 - theta)*dt
            ! A period must start and end on one level, or it spans crossings of two
            ! (0.9597 for a sine of period 1 and mean -0.25). When this step crosses the
            ! new level too, as it does once the mean has settled, the period starts
            ! there. Otherwise the next crossing of the new level starts the period
            ! (phase 3) and one line is not logged. Moving the start back along the edge
            ! would extrapolate: 1.15 periods for a sine with mean -0.9.
            if (state%v_mean >= v_prev .and. state%v_mean < v_now) then
               state%t_cross = crossing_time(t_now, dt, crossing_theta(v_prev, v_now, state%v_mean))
            else
               state%restart = .true.
            end if
         end if
      else
         if (state%n_cross >= 1) then
            state%v_sum = state%v_sum + 0.5d0*(v_prev + v_now)*dt
         else if (state%running .and. active) then
            state%n_samples = state%n_samples + 1
            state%v_mean = state%v_mean + (v_now - state%v_mean)/state%n_samples
         end if
      end if

   end subroutine crossing_update

end module nekstab_zero_cross
