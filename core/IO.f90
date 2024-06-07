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
      !   subroutine load_files(Q, mstart, kd, fname)
      !      use krylov_subspace
      !      implicit none
      !      include 'SIZE'
      !      include 'TOTAL'
      
      !   !     ----- Krylov basis V for the projection M*V = V*H -----
      !      integer, intent(in) :: mstart, kd
      !      type(krylov_vector), dimension(kd) :: Q
      !      character(len=3), intent(in) :: fname
      
      !      integer :: i, j
      !      character(len=7) tl
      !      character(len=20) fmt
      !      character(len=60) filename
      
      !   !----Upload the snapshots - ----
      
      !      do i = 1, mstart
      
      !         if (i < 10) fmt = '(i1.1)' ! 1 to 9
      !         if (i >= 10 .and. i < 100) fmt = '(i2.2)' ! 10 to 99
      !         if (i >= 100 .and. i < 1000) fmt = '(i3.3)' ! 100 to 999
      !         if (i >= 1000 .and. i < 10000) fmt = '(i4.4)' ! 1000 to 9999
      !         if (i >= 10000 .and. i < 100000) fmt = '(i5.5)' ! 10000 to 99999
      
      !         j = i
      !         write (tl, fmt) j
      !         if (i < 10)
      !     $         filename = trim(fname)//trim(SESSION)//'0.f0000'//trim(tl)
      !         if (i >= 10 .and. i < 100)
      !     $         filename = trim(fname)//trim(SESSION)//'0.f000'//trim(tl)
      !         if (i >= 100 .and. i < 1000)
      !     $         filename = trim(fname)//trim(SESSION)//'0.f00'//trim(tl)
      !         if (i >= 1000 .and. i < 10000)
      !     $         filename = trim(fname)//trim(SESSION)//'0.f0'//trim(tl)
      !         if (i >= 10000 .and. i < 100000)
      !     $         filename = trim(fname)//trim(SESSION)//'0.f'//trim(tl)
      
      !         call load_fld(filename)
      !         call opcopy(Q(i)%vx, Q(i)%vy, Q(i)%vz, vx, vy, vz)
      !         if (ifpo) call copy(Q(i)%pr, pr, nx2*ny2*nz2*nelv)
      !         if (ifto) call copy(Q(i)%t, t(1, 1, 1, 1, 1), nx1*ny1*nz1*nelv)
      
      !      end do
      !   end subroutine load_files
      !   !----------------------------------------------------------------------
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
