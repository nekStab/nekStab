! Advance one sine, or one sine plus a constant, through crossing_update.
! Line 1: sine
! Line 2: n dt period offset v_mean_init
! One line per closed period: period, level used, level stored for the next period.
! Last line: done, crossing count, stored level.
program zero_cross_driver
   use nekstab_zero_cross, only: crossing_reset, crossing_state, crossing_update
   implicit none
   character(len=16) :: cmd
   integer :: n, i, phase, ios
   real :: dt, period, offset, v_init, t, v_prev, v_now
   real :: theta, t_hit, p, v_logged, dp
   type(crossing_state) :: state
   real, parameter :: pi = acos(-1.0d0)

   read (*, *, iostat=ios) cmd
   if (ios /= 0) stop 1
   if (cmd /= 'sine') stop 1
   read (*, *) n, dt, period, offset, v_init
   if (n < 2 .or. dt <= 0.0d0 .or. period <= 0.0d0) stop 1

   call crossing_reset(state, v_init)
   v_prev = offset
   do i = 1, n - 1
      t = real(i, kind(t))*dt
      v_now = offset + sin(2.0d0*pi*t/period)
      call crossing_update(state, v_prev, v_now, t, dt, .true., &
                           phase, theta, t_hit, p, v_logged, dp)
      if (phase == 2) write (6, '(3es24.16)') p, v_logged, state%v_mean
      v_prev = v_now
   end do
   write (6, '(a,1x,i6,1x,es24.16)') 'done', state%n_cross, state%v_mean
end program zero_cross_driver
