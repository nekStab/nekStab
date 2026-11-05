! Advance one sine, or one sine plus a constant, through crossing_update.
! Line 1: sine
! Line 2: n dt period offset v_mean_init
! Or line 1: switch, for a sine that changes period at t_switch without a jump.
! Line 2: n dt period offset v_mean_init t_switch period_after
! Or line 1: cosine, which starts at its maximum, with the level from the samples.
! Line 2: n dt period offset
! One line per closed period: period, level used, level stored for the next period.
! Last line: done, crossing count, stored level.
program zero_cross_driver
   use nekstab_zero_cross, only: crossing_reset, crossing_state, crossing_update
   implicit none
   character(len=16) :: cmd
   integer :: n, i, phase, ios
   real :: dt, period, offset, v_init, t, v_prev, v_now
   real :: theta, t_hit, p, v_logged, dp, t_switch, period_after, angle
   type(crossing_state) :: state
   real, parameter :: pi = acos(-1.0d0)

   read (*, *, iostat=ios) cmd
   if (ios /= 0) stop 1
   t_switch = huge(t_switch)
   period_after = 0.0d0
   select case (cmd)
   case ('sine')
      read (*, *) n, dt, period, offset, v_init
   case ('switch')
      read (*, *) n, dt, period, offset, v_init, t_switch, period_after
   case ('cosine')
      read (*, *) n, dt, period, offset
      v_init = offset + 1.0d0
   case default
      stop 1
   end select
   if (n < 2 .or. dt <= 0.0d0 .or. period <= 0.0d0) stop 1

   if (cmd == 'cosine') then
      call crossing_reset(state, v_init, running=.true.)
      v_prev = v_init
   else
      call crossing_reset(state, v_init)
      v_prev = offset
   end if
   do i = 1, n - 1
      t = real(i, kind(t))*dt
      if (t < t_switch) then
         angle = 2.0d0*pi*t/period
      else
         angle = 2.0d0*pi*(t_switch/period + (t - t_switch)/period_after)
      end if
      v_now = offset + sin(angle)
      if (cmd == 'cosine') v_now = offset + cos(angle)
      call crossing_update(state, v_prev, v_now, t, dt, .true., &
                           phase, theta, t_hit, p, v_logged, dp)
      if (phase == 2) write (6, '(3es24.16)') p, v_logged, state%v_mean
      v_prev = v_now
   end do
   write (6, '(a,1x,i6,1x,es24.16)') 'done', state%n_cross, state%v_mean
end program zero_cross_driver
