      subroutine prepare_linearized_solver
      implicit none
      include 'SIZE'
      include 'TOTAL'

! --> Force only single perturbation mode.
      if (param(31) .gt. 1) then
         write(6, *) 'nekStab not ready for npert > 1 -- jp loops not yet implemented. Stoppiing.'
         call nek_end
      endif
      param(31) = 1 ; npert = param(31)

! --> Force deactivate OIFS.
      if (ifchar) write(6, *) 'OIFS not working with linearized solver. Turning it off.'
      ifchar = .false. ; call bcast(ifchar, lsize)

! --> Encore CFL target for EXTk
      if (param(26) .gt. 0.5) then
         if (nid .eq. 0) write(6, *) "Force the maximum CFL to 0.5 for practical reasons."
         param(26) = 0.5D+00 ; ctarg = param(26)
      endif

! --> Set nsteps/endTime accordingly.
      if (param(10) .gt. 0) then
         call compute_cfl(dt, vx, vy, vz, 1.0)
         dt = ctarg / dt
         nsteps = ceiling(param(10) / dt)
         dt = param(10) / nsteps
         if(nid.eq.0)write(6,*)'endTime specified! computing CFL from BASE FLOW!'
         if(nid.eq.0)write(6,*)' computing timeStep dt=',dt
         if(nid.eq.0)write(6,*)' computing numSteps=',nsteps
         if(nid.eq.0)write(6,*)' sampling period =',nsteps*dt
         param(12) = dt
      endif

! --> Force constant time step.
      param(12) = -abs(param(12))

! --> Broadcast parameters.
      call bcast(param, 200*wdsize)

      return
      end subroutine prepare_linearized_solver





      subroutine matvec(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      real, dimension(lt) :: fx, fy, fz, ft
      real, dimension(lt2) :: fp

      real :: umax
      integer :: n

! --> Pass the baseflow to vx, vy, vz
      call opcopy(vx, vy, vz, ubase, vbase, wbase)
      if (ifheat) call copy(t, tbase, n)

! --> Standard setup for the linearized solver.
      call prepare_linearized_solver()

! --> Direct solver only.
      if (uparam(01) .eq. 3) then
          evop = 'd'
          call forward_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      endif

! --> Adjoint solver only.
      if (uparam(01) .eq. 3.2) then
          evop = 'a'
          call adjoint_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      endif

! --> Direct-Adjoint for optimal transient growth.

! --> Linearized forward map for the Newton-Krylov solver.
      if (floor(uparam(01)) .eq. 1) then
          call newton_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      endif

! --> Adjoint solver for the steady force sensitivity analysis.

      return
      end subroutine matvec





      subroutine forward_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      real, dimension(lt) :: fx, fy, fz, ft
      real, dimension(lt2) :: fp

      real :: umax
      integer :: n

      n = nx1*ny1*nz1*nelt

! --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .false.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

! --> Turning-off the base flow side-by-side computation. Need to change for Floquet.
      ifbase = .false.

! --> Pass the initial condition for the perturbation.
      call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1), qx, qy, qz, qp, qt)

      time = 0.0D+00
      do istep = 1, nsteps
! --> Output current info to logfile.
         if(nid .eq. 0) write(6,"(' DIRECT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

         ! --> Integrate forward in time.
         call nekstab_usrchk()
         call nek_advance()
      end do

! --> Copy the solution.
      call nopcopy(fx, fy, fz, fp, ft, vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1))

      return
      end subroutine forward_linearized_map





      subroutine adjoint_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      real, dimension(lt) :: fx, fy, fz, ft
      real, dimension(lt2) :: fp

      real :: umax
      integer :: n

      n = nx1*ny1*nz1*nelt

! --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .true.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

! --> Turning-off base flow computation.
      ifbase = .false.

! --> Pass the initial condition for the perturbation.
      call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1), qx, qy, qz, qp, qt)

      time = 0.0D+00
      do istep = 1, nsteps
! --> Output current info to logfile.
         if(nid .eq. 0) write(6,"(' ADJOINT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

         ! --> Integrate forward in time.
         call nekstab_usrchk()
         call nek_advance()
      end do

! --> Copy the solution.
      call nopcopy(fx, fy, fz, fp, ft, vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1))

      return
      end subroutine adjoint_linearized_map





      subroutine newton_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      real, dimension(lt) :: fx, fy, fz, ft
      real, dimension(lt2) :: fp

      real :: umax
      integer :: n

! --> Evaluate exp(t*L) * q0.
      call forward_linearized_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)

! --> Evaluate (I - exp(t*L)) * q0.
      call nopsub2(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
      call nopchsign(fx, fy, fz, fp, ft)

      return
      end subroutine newton_linearized_map
