! Drive the Newton periodic-orbit rules.
! Line 1 is the command: step, cfl, filter, keep or clear.
! step:   period, n_kept, dt_in  -> dt, nsteps
! cfl:    n_kept, target, cfl    -> 1 when the warning fires
! filter: is_po, is_po_t, type, weight -> 1 when the warning fires
! keep:   iteration, is_po, is_po_t, nsteps, stored -> stored count
! clear:  no arguments -> count after the stage
program newton_po_rules_driver
   use nekstab_newton_po_rules, only: newton_po_cfl_warn, newton_po_filter_warn, &
                                      newton_po_kept_steps, newton_po_set_step, &
                                      newton_po_steps_after_stage
   implicit none
   character(len=16) :: cmd
   integer :: ios, n_kept, nsteps, flag, is_po_i, is_po_t_i
   integer :: iteration, stored, kept
   real :: period, dt_in, dt, target_cfl, cfl_now, weight, filter_type
   logical :: is_po, is_po_t, warn

   read (*, *, iostat=ios) cmd
   if (ios /= 0) stop 1

   if (cmd == 'step') then
      read (*, *) period, n_kept, dt_in
      dt = dt_in
      call newton_po_set_step(period, n_kept, dt, nsteps)
      write (6, '(es24.16,1x,i8)') dt, nsteps
   else if (cmd == 'cfl') then
      read (*, *) n_kept, target_cfl, cfl_now
      warn = newton_po_cfl_warn(n_kept, target_cfl, cfl_now)
      flag = 0
      if (warn) flag = 1
      write (6, '(i1)') flag
   else if (cmd == 'filter') then
      read (*, *) is_po_i, is_po_t_i, filter_type, weight
      is_po = is_po_i /= 0
      is_po_t = is_po_t_i /= 0
      warn = newton_po_filter_warn(is_po, is_po_t, filter_type, weight)
      flag = 0
      if (warn) flag = 1
      write (6, '(i1)') flag
   else if (cmd == 'keep') then
      read (*, *) iteration, is_po_i, is_po_t_i, nsteps, stored
      is_po = is_po_i /= 0
      is_po_t = is_po_t_i /= 0
      kept = newton_po_kept_steps(iteration, is_po, is_po_t, nsteps, stored)
      write (6, '(i8)') kept
   else if (cmd == 'clear') then
      write (6, '(i8)') newton_po_steps_after_stage()
   else
      stop 1
   end if
end program newton_po_rules_driver
