      !----------------------------------------------------------------------
      subroutine whereyouwant(resnam, posfil) !file numbering suffix counter
         character(len=3) resnam
         integer posfil, iprefix
         common/nopenf2/nopen(1000, 2)
      
         iprefix = i_find_prefix(resnam, 99)
         nopen(iprefix, 1) = posfil - 1
      
      end
      !----------------------------------------------------------------------
      subroutine load_files(Q, mstart, kd, fname)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      ! Krylov basis V for the projection M*V = V*H
         integer, intent(in) :: mstart, kd
         type(krylov_vector), dimension(kd) :: Q
         character(len=3), intent(in) :: fname
      
         integer :: i
         character(len=60) filename
      
      ! Upload the snapshots
         do i = 1, mstart
            write (filename, '(A,A,"0.f",I5.5)') trim(fname), trim(SESSION), i
            call load_fld(filename)
            call nopcopy(Q(i)%vx, Q(i)%vy, Q(i)%vz, Q(i)%pr, Q(i)%t, vx, vy, vz, pr, t)
         end do
      end subroutine load_files
      !----------------------------------------------------------------------
      subroutine k_load(Q, fname)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         type(krylov_vector), intent(out) :: Q
         character(len=*), intent(in) :: fname
         character(len=256) :: fname_local ! ifx: trim() on assumed-length args can segfault
         character(len=60) :: filename

         fname_local = fname ! copy to local before trim() for ifx compatibility

         if (index(fname_local, '.f') == 0) then
            write (filename, '(2A)') trim(fname_local), trim(SESSION)//'0.f00001'
         else
            filename = fname_local
         end if
      
         call load_fld(filename)
         call nopcopy(Q%vx, Q%vy, Q%vz, Q%pr, Q%t, vx, vy, vz, pr, t)
      
      end subroutine k_load
