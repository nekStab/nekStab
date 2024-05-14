      !----------------------------------------------------------------------
      subroutine whereyouwant(resnam, posfil) !file numbering suffix counter
         character*3 resnam
         integer posfil
         common/RES_WANT/nopen(99, 2)
      
      !     change prepost.f line 1094 from "save nopen" to "common /RES_WANT/ nopen"
      
         iprefix = i_find_prefix(resnam, 99)
         nopen(iprefix, 1) = posfil - 1
      
         return
      end
      !----------------------------------------------------------------------
      subroutine load_files(Q, mstart, kd, fname)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Krylov basis V for the projection M*V = V*H -----
         type(krylov_vector), dimension(kd) :: Q
      
         integer :: kd, mstart, i, j
         character*3 fname
         character*60 filename
         character(len=7) :: tl
         character(len=20) :: fmt
      
      !----Upload the snapshots - ----
      
         do i = 1, mstart
      
            if (i < 10) fmt = '(i1.1)' ! 1 to 9
            if (i >= 10 .and. i < 100) fmt = '(i2.2)' ! 10 to 99
            if (i >= 100 .and. i < 1000) fmt = '(i3.3)' ! 100 to 999
            if (i >= 1000 .and. i < 10000) fmt = '(i4.4)' ! 1000 to 9999
            if (i >= 10000 .and. i < 100000) fmt = '(i5.5)' ! 10000 to 99999
      
            j = i
            write (tl, fmt) j
            if (i < 10)
     $   filename = trim(fname)//trim(SESSION)//'0.f0000'//trim(tl)
            if (i >= 10 .and. i < 100)
     $   filename = trim(fname)//trim(SESSION)//'0.f000'//trim(tl)
            if (i >= 100 .and. i < 1000)
     $   filename = trim(fname)//trim(SESSION)//'0.f00'//trim(tl)
            if (i >= 1000 .and. i < 10000)
     $   filename = trim(fname)//trim(SESSION)//'0.f0'//trim(tl)
            if (i >= 10000 .and. i < 100000)
     $   filename = trim(fname)//trim(SESSION)//'0.f'//trim(tl)
      
            call load_fld(filename)
            call opcopy(Q(i)%vx, Q(i)%vy, Q(i)%vz, vx, vy, vz)
            if (ifpo) call copy(Q(i)%pr, pr, nx2*ny2*nz2*nelv)
            if (ifto) call copy(Q(i)%t, t(1, 1, 1, 1, 1), nx1*ny1*nz1*nelv)
      
         end do
         return
      end subroutine load_files
      !----------------------------------------------------------------------
