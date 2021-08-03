      subroutine prepare_linearized_solver

      implicit none
      include 'SIZE'
      include 'TOTAL'

!     --> Force only single perturbation mode.
      if (param(31) .gt. 1) then
         write(6, *) 'nekStab not ready for npert > 1 -- jp loops not yet implemented. Stoppiing.'
         call nek_end
      endif
      param(31) = 1 ; npert = param(31)

!     --> Force deactivate OIFS.
      if (ifchar) write(6, *) 'OIFS not working with linearized solver. Turning it off.'
      ifchar = .false. ; call bcast(ifchar, lsize)

!     --> Encore CFL target for EXTk
      if (param(26) .gt. 0.5) then
         if (nid .eq. 0) write(6, *) "Force the maximum CFL to 0.5 for practical reasons."
         param(26) = 0.5D+00 ; ctarg = param(26)
      endif

!     --> Set nsteps/endTime accordingly.
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
         FINTIM = nsteps*dt
         !NSTEPS = PARAM(11)
      endif

!     --> Force constant time step.
      param(12) = -abs(param(12))

!     --> Broadcast parameters.
      call bcast(param, 200*wdsize)

      return
      end subroutine prepare_linearized_solver





!-----------------------------------------------------------------------





      subroutine matvec(f, q)

!     Dispatch the correct matrix-vector product to the Arnoldi factorization.
!     All subroutines need to have the same interface.
!     
!     NOTE : The baseflow needs to be pass to (ubase, vbase, wbase, tbase)
!     before this function is called.
!     
!     INPUTS
!     ------
!     
!     qx, qy, qz, qt : nek-arrays of size lt
!     Initial velocity and temperature components.
!     
!     qp : nek-array of size lt2
!     Initial pressure component.
!     
!     OUTPUTS
!     -------
!     
!     fx, fy, fz, ft : nek-arrays of size lt
!     Final velocity and temperature components.
!     
!     fp : nek-array of size lt2
!     Final pressure component.
!     

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      type(krylov_vector) :: q, f

      integer :: n

      n = nx1*ny1*nz1*nelt

!     --> Pass the baseflow to vx, vy, vz
      call opcopy(vx, vy, vz, ubase, vbase, wbase)
      if (ifheat) call copy(t, tbase, n)

!     --> Standard setup for the linearized solver.
      call prepare_linearized_solver()

      !     --> Linearized forward map for the Newton-Krylov solver.
      if (floor(uparam(01)) .eq. 2) then
            call newton_linearized_map(f, q)
         endif

!     --> Direct solver only steady and periodic!
      if (uparam(01) .ge. 3.0 .and. uparam(01) .lt. 3.2 ) then
         evop = 'd'
         call forward_linearized_map(f, q)
      endif

!     --> Adjoint solver only steady and periodic!
      if (uparam(01) .ge. 3.2 .and. uparam(01) .lt. 3.3 ) then
         evop = 'a'
         call adjoint_linearized_map(f, q)
      endif

!     --> Direct-Adjoint for optimal transient growth steady only.
      if (uparam(01) .eq. 3.3) then
         evop="p"
         call transient_growth_map(f, q)
      endif

      !     --> Adjoint solver for the steady force sensitivity analysis.
      if (floor(uparam(01)) .eq. 4) then
         call ts_force_sensitivity_map(f, q)
      end if

      return
      end subroutine matvec





!-----------------------------------------------------------------------





      subroutine forward_linearized_map(f, q)

!     Integrate forward in time the linearized Navier-Stokes equations.
!     Denoting by L the Jacobian of the Navier-Stokes equations, the corresponding
!     matrix vector product is thus
!     
!     x(t) = exp(t * L) * x(0)
!     
!     where x(0) is the initial condition (qx, qy, qz, qp, qt) and x(t) the final
!     one (fx, fy, fz, fp, ft).

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      type(krylov_vector) :: q, f

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      real, save, allocatable, dimension(:, :) :: uor,vor,wor

      logical, save :: init
      data             init /.false./

!     --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .false.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

!     --> Turning-off the base flow side-by-side computation. Need to change for Floquet.
      ifbase = .false.
      if(uparam(01) .eq. 3.11)ifbase=.true. ! activate floquet
      if(ifstorebase.and.init)ifbase=.false.

      if(ifstorebase .and.ifbase.and..not.init)then
       if(nid .eq. 0) write(6,*) 'ALLOCATING UOR WITH NSTEPS:',nsteps
       allocate(uor(lt, nsteps), vor(lt, nsteps), wor(lt,nsteps))
      endif

!     --> Pass the initial condition for the perturbation.
      call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1), q%vx, q%vy, q%vz, q%pr, q%theta)

      time = 0.0D+00
      do istep = 1, nsteps
!     --> Output current info to logfile.
         if(nid .eq. 0) write(6,"(' DIRECT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

         ! --> Integrate forward in time.
         call nekstab_usrchk()
         call nek_advance()
         if    (ifstorebase .and.ifbase .and..not.  init)then !storing first time
           if(nid.eq.0)write(6,*)'storing first series:',istep,'/',nsteps
           call opcopy(uor(:,istep),vor(:,istep),wor(:,istep),vx,vy,vz)
         elseif(ifstorebase .and.init   .and..not.ifbase)then !just moving in memory
           call opcopy(vx,vy,vz,uor(:,istep),vor(:,istep),wor(:,istep))
         endif
      end do
      if(ifstorebase .and..not.init.and.ifbase)then
        ifbase=.false.;init=.true.
      endif

!     --> Copy the solution.
      call nopcopy(f%vx, f%vy, f%vz, f%pr, f%theta, vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1))

      return
      end subroutine forward_linearized_map





!-----------------------------------------------------------------------





      subroutine adjoint_linearized_map(f, q)

!     Integrate forward in time the adjoint Navier-Stokes equations.
!     Denoting by L adjoint Navier-Stokes operator, the corresponding
!     matrix vector product is thus
!     
!     x(t) = exp(t * L) * x(0)
!     
!     where x(0) is the initial condition (qx, qy, qz, qp, qt) and x(t) the final
!     one (fx, fy, fz, fp, ft).

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      type(krylov_vector) :: q, f

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      real, save, allocatable, dimension(:, :) :: uor,vor,wor

      logical, save :: init
      data             init /.false./
!     --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .true.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

!     --> Turning-off the base flow side-by-side computation. Need to change for Floquet.
      ifbase = .false.
      if(uparam(01) .eq. 3.21)ifbase=.true. ! activate floquet
      if(ifstorebase.and.init)ifbase=.false.

      if(ifstorebase .and.ifbase.and..not.init)then
       if(nid .eq. 0) write(6,*) 'ALLOCATING UOR WITH NSTEPS:',nsteps
       allocate(uor(lt, nsteps), vor(lt, nsteps), wor(lt,nsteps))
      endif

!     --> Pass the initial condition for the perturbation.
      call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1), q%vx, q%vy, q%vz, q%pr, q%theta)

      time = 0.0D+00
      do istep = 1, nsteps
!     --> Output current info to logfile.
         if(nid .eq. 0) write(6,"(' ADJOINT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

         ! --> Integrate forward in time.
         call nekstab_usrchk()
         call nek_advance()
         if    (ifstorebase .and.ifbase .and..not.  init)then !storing first time
           if(nid.eq.0)write(6,*)'storing first series:',istep,'/',nsteps
           call opcopy(uor(:,istep),vor(:,istep),wor(:,istep),vx,vy,vz)
         elseif(ifstorebase .and.init   .and..not.ifbase)then !just moving in memory
           call opcopy(vx,vy,vz,uor(:,istep),vor(:,istep),wor(:,istep))
         endif
      end do
      if(ifstorebase .and..not.init.and.ifbase)then
        ifbase=.false.;init=.true.
      endif

!     --> Copy the solution.
      call nopcopy(f%vx, f%vy, f%vz, f%pr, f%theta, vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, 1, 1))

      return
      end subroutine adjoint_linearized_map





!-----------------------------------------------------------------------





      subroutine newton_linearized_map(f, q)

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      type(krylov_vector) :: f, q
      type(krylov_vector) :: bvec, btvec

!     ----------------------------------
!     -----     REGULAR NEWTON     -----
!     ----------------------------------

!     --> Evaluate exp(t*L) * q0.
      call forward_linearized_map(f, q)

!     --> Evaluate (exp(t*L) - I) * q0.
      call krylov_sub2(f, q)

!     ----------------------------------
!     -----     NEWTON FOR UPO     -----
!     ----------------------------------

      if ( uparam(3) .eq. 3.1 ) then

         call krylov_zero(bvec)
         call krylov_zero(btvec)

         call compute_bvec(bvec, fc_nwt)
         call krylov_cmult(bvec, q%time)
         call krylov_add2(f, bvec)

         call compute_bvec(btvec, ic_nwt)
         call krylov_inner_product(f%time, btvec, q)
      else
         f%time = 0.0D+00
      end if

      return
      end subroutine newton_linearized_map


      subroutine compute_bvec(bvec, qbase)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      type(krylov_vector) :: bvec, qbase
      type(krylov_vector) :: wrk1, wrk2

!     --> Setup the paramtemers for the solver.
      ifpert = .false. ; ifadj = .false.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

      ifbase = .false.
      call bcast(ifbase, lsize)

!     --> Pass the initial condition.
      call nopcopy(vx, vy, vz, pr, t, qbase%vx, qbase%vy, qbase%vz, qbase%pr, qbase%theta)
      param(10) = qbase%time
      call newton_krylov_prepare()
      call krylov_copy(wrk1, qbase)

!     --> Single time-step to approximate the time-derivative.
      time = 0.0D+00
      do istep = 1, 1
         call nekStab_usrchk()
         call nek_advance()
      enddo
      call nopcopy(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%theta, vx, vy, vz, pr, t)

!     --> Approximate the time-derivative.
      call krylov_sub2(wrk2, wrk1)
      call krylov_copy(bvec, wrk2)
      call krylov_cmult(bvec, 1.0/dt)
      bvec%time = 0.0D+00

      return
      end subroutine compute_bvec



!-----------------------------------------------------------------------





      subroutine ts_force_sensitivity_map(f, q)

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      type(krylov_vector) :: f, q

!     --> Evaluate exp(t*L) * q0.
      call adjoint_linearized_map(f, q)

!     --> Evaluate (I - exp(t*L)) * q0.
      call krylov_sub2(f, q)
      call krylov_cmult(f, -1.0D+00)

      return
      end subroutine ts_force_sensitivity_map




!-----------------------------------------------------------------------




      subroutine transient_growth_map(f, q)

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      type(krylov_vector) :: f, q, wrk

!     --> Evaluate the forward map.
      call forward_linearized_map(wrk, q)

!     --> Evaluate the adjoint map.
      call adjoint_linearized_map(f, wrk)

      return
      end subroutine transient_growth_map
