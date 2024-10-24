      !----------------------------------------------------------------------
      subroutine whereyouwant(resnam, posfil) !file numbering suffix counter
         character(len=3) resnam
         integer posfil, iprefix
         common/RES_WANT/nopen(99, 2)
      
      !     change prepost.f line 1094 from "save nopen" to "common /RES_WANT/ nopen"
      
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
