      !-----------------------------------------------------------------------
      ! IO.f90 -- File I/O utilities for nekStab
      !
      ! Purpose:
      !   Provides routines for file numbering, loading eigenvalues,
      !   loading eigenmodes, and reading snapshot sequences from disk.
      !
      ! Public interface:
      !   whereyouwant    -- set file numbering suffix counter
      !   read_eigenvalue -- read leading eigenvalue from convergence file
      !   load_mode_pair  -- load real/imaginary eigenmode pair
      !   load_files      -- load a sequence of Krylov snapshots
      !   k_load          -- load a single field into a krylov_vector
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL, PARALLEL, TSTEP
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! whereyouwant -- Set file numbering suffix counter
      !
      ! Arguments:
      !   resnam [in] -- 3-character file prefix
      !   posfil [in] -- starting file number
      !-----------------------------------------------------------------------
      subroutine whereyouwant(resnam, posfil)
         character(len=3), intent(in) :: resnam
         integer, intent(in) :: posfil
         integer :: iprefix
         common/nopenf2/nopen(1000, 2)
      
         iprefix = i_find_prefix(resnam, 99)
         nopen(iprefix, 1) = posfil - 1
      
      end subroutine whereyouwant

      !-----------------------------------------------------------------------
      ! read_eigenvalue -- Read leading eigenvalue from convergence file
      !
      ! Purpose:
      !   Reads (sigma, omega) from Spectre_NSd_conv.dat on rank 0
      !   and broadcasts to all MPI ranks.
      !
      ! Arguments:
      !   sigma [out] -- real part (growth rate)
      !   omega [out] -- imaginary part (frequency), returned as abs
      !-----------------------------------------------------------------------
      subroutine read_eigenvalue(sigma, omega)
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

      !-----------------------------------------------------------------------
      ! load_mode_pair -- Load real/imaginary eigenmode pair from disk
      !
      ! Purpose:
      !   Loads the real and imaginary parts of a direct ('d') or
      !   adjoint ('a') eigenmode pair into krylov_vectors.
      !
      ! Arguments:
      !   mode [in]  -- 'd' for direct, 'a' for adjoint
      !   Re   [out] -- real part of the eigenmode
      !   Im   [out] -- imaginary part of the eigenmode
      !-----------------------------------------------------------------------
      subroutine load_mode_pair(mode, Re, Im)
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

      !-----------------------------------------------------------------------
      ! load_files -- Load a sequence of snapshots into Krylov vectors
      !
      ! Arguments:
      !   Q      [out]   -- array of krylov_vectors to fill
      !   mstart [in]    -- number of snapshots to load
      !   kd     [in]    -- declared dimension of Q
      !   fname  [in]    -- 3-character file prefix
      !-----------------------------------------------------------------------
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

      !-----------------------------------------------------------------------
      ! k_load -- Load a single field file into a krylov_vector
      !
      ! Arguments:
      !   Q     [out] -- krylov_vector to fill
      !   fname [in]  -- file prefix or full filename
      !-----------------------------------------------------------------------
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
