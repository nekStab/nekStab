      !----------------------------------------------------------------------
      subroutine whereyouwant(resnam, posfil) !file numbering suffix counter
         character(len=3) resnam
         integer posfil, iprefix
         common/nopenf2/nopen(1000, 2)
      
         iprefix = i_find_prefix(resnam, 99)
         nopen(iprefix, 1) = posfil - 1
      
      end
      !----------------------------------------------------------------------
      subroutine read_eigenvalue(sigma, omega)
!
!     Read the leading eigenvalue (sigma, omega) from the convergence
!     file written by the eigensolver, and broadcast to all MPI ranks.
!
!     OUTPUT
!       sigma : real part (growth rate)
!       omega : imaginary part (frequency), returned as abs(omega)
!
         implicit none
         include 'SIZE'
         include 'PARALLEL'  ! nid
         include 'TSTEP'     ! wdsize

         real, intent(out) :: sigma, omega

         if (nid == 0) then
            open (unit=10, file='Spectre_NSd_conv.dat',
     $         status='old', action='read')
            read (10, '(2E15.7)') sigma, omega
            close (10)
            omega = abs(omega)
            write (6, *) 'read_eigenvalue: sigma =', sigma,
     $         ' omega =', omega
         end if
         call bcast(sigma, wdsize)
         call bcast(omega, wdsize)

      end subroutine read_eigenvalue
      !----------------------------------------------------------------------
      subroutine load_mode_pair(mode, Re, Im)
!
!     Load the real and imaginary parts of a direct ('d') or adjoint ('a')
!     eigenmode pair into krylov_vectors Re and Im.
!
         use krylov_subspace
         implicit none
         include 'SIZE'

         character(len=*), intent(in) :: mode
         type(krylov_vector), intent(out) :: Re, Im

         if (mode == 'd') then
            call k_load(Re, 'dRe')
            call k_load(Im, 'dIm')
         else if (mode == 'a') then
            call k_load(Re, 'aRe')
            call k_load(Im, 'aIm')
         end if

      end subroutine load_mode_pair
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
