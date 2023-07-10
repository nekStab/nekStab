      subroutine prepare_linearized_solver

      implicit none
      include 'SIZE'
      include 'TOTAL'

      call nekgsync

!     --> Force only single perturbation mode.
      if (param(31) .gt. 1) then
         write(6, *) 'nekStab not ready for npert > 1 -- jp loops not yet implemented. Stoppiing.'
         call nek_end
      endif
      param(31) = 1 ; npert = param(31)

!     --> Force deactivate OIFS.
      if (ifchar) write(6, *) 'OIFS not working with linearized solver. Turning it off.'
      ifchar = .false. ; call bcast(ifchar, lsize)

!     --> Enforce CFL target for EXTk
      if (param(26) .gt. 1.0) then
         write(6, *) "Forcing target CFL to 0.5!"
         param(26) = 0.5d0
      endif

!     --> Set nsteps/endTime accordingly.
      if (param(10) .gt. 0) then
         if(nid.eq.0)write(6,*)'param(10),time=',param(10),time
         if(nid.eq.0)write(6,*)'endTime specified! Recomputing dt and nsteps to match endTime'
         call compute_cfl(ctarg, vx, vy, vz, 1.0d0) ! ctarg contains the sum ( ux_i / dx_i )
         if(nid.eq.0)write(6,*)'max spatial restriction:',ctarg
         dt = param(26) / ctarg ! dt given target CFL
         nsteps = ceiling(param(10) / dt) ! computing a safe value of nsteps
         dt = param(10) / nsteps ! reducing dt to match param(10)
         if(nid.eq.0)write(6,*)' new timeStep dt=',dt
         if(nid.eq.0)write(6,*)' new numSteps nsteps=',nsteps
         if(nid.eq.0)write(6,*)' resulting sampling period =',nsteps*dt
         param(12) = dt
         call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
         if(nid.eq.0)write(6,*)' current CFL and target=',ctarg,param(26)
         lastep = 0             ! subs1.f:279
         fintim = nsteps*dt
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
!     qx, qy, qz, qt : nek-arrays of size lv
!     Initial velocity and temperature components.
!
!     qp : nek-array of size lp
!     Initial pressure component.
!
!     OUTPUTS
!     -------
!
!     fx, fy, fz, ft : nek-arrays of size lv
!     Final velocity and temperature components.
!
!     fp : nek-array of size lp
!     Final pressure component.
!

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      type(krylov_vector) :: q, f

      logical, save :: init
      data init /.false./
      integer m

!     --> Pass the baseflow to vx, vy, vz
      call opcopy(vx, vy, vz, ubase, vbase, wbase)
      if (ifheat) call copy(t(1,1,1,1,1), tbase(1,1,1,1,1),  nx1*ny1*nz1*nelv)
      if (ldimt.gt.1) then
         do m = 2,ldimt
            call copy(t(1,1,1,1,m), tbase(1,1,1,1,m),  nx1*ny1*nz1*nelv)
         enddo
      endif
      if(ifbf2d .and. if3d)then
         call rzero(vz,nx1*ny1*nz1*nelv);if(nid.eq.0)write(6,*)'Forcing vz=0'
      endif

!     --> Standard setup for the linearized solver.
      if(.not.init) then
         call prepare_linearized_solver
         init = .true.
      endif

      lastep = 0
      fintim = param(10)

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

!     --> Direct-Adjoint for optimal transient growth.
      if (uparam(01) .ge. 3.3 .and. uparam(01) .lt. 3.4) then
         evop="p"
         call transient_growth_map(f, q)
      endif

!     --> Adjoint solver for the steady force sensitivity analysis.
      if (floor(uparam(01)) .eq. 4) then
         call ts_force_sensitivity_map(f, q)
      end if

!     --> Linearized forward map for the Newton-Krylov solver.
      if (floor(uparam(01)) .eq. 2) then
         evop = 'n'
         call newton_linearized_map(f, q)
         init = .false.
      endif

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

      logical, save :: init
      data             init /.false./
      n = nx1*ny1*nz1*nelv

!     --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .false.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

!     --> Turning-off the base flow side-by-side computation.
      ifbase=.false.
      if(uparam(01) .eq. 3.11)ifbase=.true. ! activate Floquet
      if(uparam(01) .eq. 3.31)ifbase=.true. ! activate Floquet for intracycle transient growth
      if(uparam(01) .eq. 2.1 .or. uparam(01) .eq. 2.2)then
         init=.true.            ! use stored baseflow if ifstorebase
         ifbase=.true.          ! activate baseflow evolution for UPO
      endif
      if(ifstorebase.and.init)ifbase=.false. ! deactivte ifbase if baseflow stored

      if(ifstorebase .and.ifbase.and..not.init)then
         if(nid .eq. 0) write(6,*) 'ALLOCATING ORBIT WITH NSTEPS:',nsteps
         allocate(uor(lv, nsteps), vor(lv, nsteps))
         if(if3d)then
            allocate(wor(lv,nsteps))
         else                   ! 2D
            allocate(wor(1,1))
         endif
         if(ifheat)allocate(tor(lv, nsteps))
      endif

!     --> Pass the initial condition for the perturbation.
      call nopcopy(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1),
     $     q%vx, q%vy, q%vz, q%pr, q%theta)

      time = 0.0D+00
      do istep = 1, nsteps
!     --> Output current info to logfile.
         if(nid .eq. 0) write(6,"(' DIRECT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

         ! --> Integrate forward in time.
         call nekstab_usrchk()
         call nek_advance()

         if    (ifstorebase .and.ifbase .and..not.  init)then !storing first time
            if(nid.eq.0)write(6,*)'storing first series:',istep,'/',nsteps
            call opcopy(uor(1,istep),vor(1,istep),wor(1,istep),vx,vy,vz)
            if(ifheat)call copy(tor(1,istep),t(1,1,1,1,1),n)
         elseif(ifstorebase .and.init .and..not.ifbase)then !just moving in memory
            if(nid.eq.0)write(6,*)'using stored baseflow'
            call opcopy(vx,vy,vz,uor(1,istep),vor(1,istep),wor(1,istep))
            if(ifheat)call copy(t(1,1,1,1,1),tor(1,istep),n)
         endif
      end do
      if(ifstorebase .and..not.init.and.ifbase)then
         ifbase=.false.;init=.true.
      endif

!     --> Copy the solution.
      call nopcopy(f%vx,f%vy,f%vz,f%pr,f%theta,
     $     vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:, :,1))

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

      logical, save :: init
      data             init /.false./
      n = nx1*ny1*nz1*nelv

!     --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .true.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

!     --> Turning-off the base flow side-by-side computation. Need to change for Floquet.
      ifbase=.false.
      if(uparam(01) .eq. 3.21)ifbase=.true. ! activate floquet

      if(ifstorebase.and.init)ifbase=.false.

      if(ifstorebase .and.ifbase.and..not.init)then
         if(nid .eq. 0) write(6,*) 'ALLOCATING ORBIT WITH NSTEPS:',nsteps
         allocate(uor(lv, nsteps), vor(lv, nsteps))
         if(if3d)then
            allocate(wor(lv,nsteps))
         else                   ! 2D
            allocate(wor(1,1))
         endif
         if(ifheat)allocate(tor(lv, nsteps))
      endif
      if(uparam(01) .eq. 3.31)init=.true. ! activate Floquet for intracycle transient growth (base flow already computed)

!     --> Pass the initial condition for the perturbation.
      call nopcopy(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1),
     $     q%vx,    q%vy,    q%vz,    q%pr,    q%theta)

      time = 0.0D+00
      do istep = 1, nsteps
!     --> Output current info to logfile.
         if(nid .eq. 0) write(6,"(' ADJOINT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt

         ! --> Integrate backward in time.
         call nekstab_usrchk()
         call nek_advance()

         if    (ifstorebase .and.ifbase .and..not.  init)then !storing first time
            if(nid.eq.0)write(6,*)'storing first series:',istep,'/',nsteps
            call opcopy(uor(1,istep),vor(1,istep),wor(1,istep),vx,vy,vz)
            if(ifheat)call copy(tor(1,istep),t(1,1,1,1,1),n)
         elseif(ifstorebase .and.init .and..not.ifbase)then !just moving in memory
            if(nid.eq.0)write(6,*)'using stored baseflow'
            call opcopy(vx,vy,vz,uor(1,istep),vor(1,istep),wor(1,istep))
            if(ifheat)call copy(t(1,1,1,1,1),tor(1,istep),n)
         endif
      end do
      if(ifstorebase .and..not.init.and.ifbase)then
         ifbase=.false.;init=.true.
      endif

!     --> Copy the solution.
      call nopcopy(f%vx,f%vy,f%vz,f%pr,f%theta,
     $     vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1))

      return
      end subroutine adjoint_linearized_map


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


      subroutine newton_linearized_map(f, q)


      use krylov_subspace

      implicit none
      include 'SIZE'
      include 'TOTAL'

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

      if ( uparam(1) .eq. 2.1 ) then

         call krylov_zero(bvec)
         call krylov_zero(btvec)

         call compute_bvec(bvec, fc_nwt)
         call krylov_cmult(bvec, q%time)
         call krylov_add2(f, bvec)

         call compute_bvec(btvec, ic_nwt)
         f%time = krylov_inner_product(btvec, q)

         if(nid .eq. 0) write(6,*)'Newton period correction:',f%time

      else

         f%time = 0.0D+00

      end if

      return
      end subroutine newton_linearized_map



!-----------------------------------------------------------------------


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
      call prepare_linearized_solver
      wrk1 = qbase

!     --> Single time-step to approximate the time-derivative.
      time = 0.0D+00
      do istep = 1, 1
         call nekStab_usrchk()
         call nek_advance()
      enddo
      call nopcopy(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%theta,
     $     vx,      vy,      vz,      pr,      t(1,1,1,1,1))

!     --> Approximate the time-derivative.
      call krylov_sub2(wrk2, wrk1)
      bvec = wrk2
      call krylov_cmult(bvec, 1.0/dt)
      bvec%time = 0.0D+00


      return
      end subroutine compute_bvec
