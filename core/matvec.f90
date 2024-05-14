      subroutine prepare_linearized_solver
      
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         call nekgsync
      
      !     --> Force only single perturbation mode.
         if (param(31) > 1) then
            write (6, *) 'not ready for npert >1 : Stopping.'
            call nek_end
         end if
         param(31) = 1; npert = param(31)
      
      !     --> Force deactivate OIFS.
         if (ifchar) write (6, *) 'OIFS not working with linearized solver.'
         ifchar = .false.; call bcast(ifchar, lsize)
      
      !     --> Enforce CFL target for EXTk
         if (param(26) > 1.0) then
            write (6, *) "Forcing target CFL to 0.5!"
            param(26) = 0.5d0
         end if
      
      !     --> Set nsteps/endTime accordingly.
         if (param(10) > 0) then
            if (nid == 0) write (6, *) 'param(10),time=', param(10), time
            if (nid == 0) write (6, *) 'endTime specified! Recomputing dt and nsteps to match endTime'
            call compute_cfl(ctarg, vx, vy, vz, 1.0d0) ! ctarg contains the sum ( ux_i / dx_i )
            if (nid == 0) write (6, *) 'max spatial restriction:', ctarg
            dt = param(26)/ctarg ! dt given target CFL
            nsteps = ceiling(param(10)/dt) ! computing a safe value of nsteps
            dt = param(10)/nsteps ! reducing dt to match param(10)
            if (nid == 0) write (6, *) ' new timeStep dt=', dt
            if (nid == 0) write (6, *) ' new numSteps nsteps=', nsteps
            if (nid == 0) write (6, *) ' resulting sampling period =', nsteps*dt
            param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
            if (nid == 0) write (6, *) ' current CFL and target=', ctarg, param(26)
            lastep = 0             ! subs1.f:279
            fintim = nsteps*dt
         end if
      
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
         data init/.false./
      
      !     --> Pass the baseflow to vx, vy, vz
         call nopcopy(vx, vy, vz, pr, t, ubase, vbase, wbase, pbase, tbase)
      
         if (ifbf2d .and. if3d) then
            call rzero(vz, nx1*ny1*nz1*nelv)
            if (nid == 0) write (6, *) 'Forcing vz=0'
         end if
      
      !     --> Standard setup for the linearized solver.
         if (.not. init) then
            call prepare_linearized_solver
            init = .true.
         end if
      
         lastep = 0
         fintim = param(10)
      
      !     --> Direct solver only steady and periodic!
         if (uparam(01) >= 3.0 .and. uparam(01) < 3.2) then
            evop = 'd'
            if (iffindiff) then
               if (nid == 0) write (*, *) "Using the finite-difference approximation of the Fréchet derivative."
               call forward_finite_difference_map(f, q)
            else
               call forward_linearized_map(f, q)
            end if
         end if
      
      !     --> Adjoint solver only steady and periodic!
         if (uparam(01) >= 3.2 .and. uparam(01) < 3.3) then
            evop = 'a'
            call adjoint_linearized_map(f, q)
         end if
      
      !     --> Direct-Adjoint for optimal transient growth.
         if (uparam(01) >= 3.3 .and. uparam(01) < 3.4) then
            evop = "p"
            call transient_growth_map(f, q)
         end if
      
      !     --> Adjoint solver for the steady force sensitivity analysis.
         if (floor(uparam(01)) == 4) then
            call ts_force_sensitivity_map(f, q)
         end if
      
      !     --> Linearized forward map for the Newton-Krylov solver.
         if (floor(uparam(01)) == 2) then
            evop = 'n'
            call newton_linearized_map(f, q)
            init = .false.
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
      
         logical, save :: init
         data init/.false./
      
         integer m
         nt = nx1*ny1*nz1*nelt
      
      !     --> Setup the parameters for the linearized solver.
         ifpert = .true.; ifadj = .false.
         call bcast(ifpert, lsize); call bcast(ifadj, lsize)
      
      !     --> Turning-off the base flow side-by-side computation.
         ifbase = .false.
         if (uparam(01) == 3.11) ifbase = .true. ! activate Floquet
         if (uparam(01) == 3.31) ifbase = .true. ! activate Floquet for intracycle transient growth
         if (uparam(01) == 2.1 .or. uparam(01) == 2.2) then
            init = .true.            ! use stored baseflow if ifstorebase
            ifbase = .true.          ! activate baseflow evolution for UPO
         end if
         if (ifstorebase .and. init) ifbase = .false. ! deactivte ifbase if baseflow stored
      
         if (ifstorebase .and. ifbase .and. .not. init) then
            if (nid == 0) write (6, *) 'ALLOCATING ORBIT WITH NSTEPS:', nsteps
            allocate (uor(lv, nsteps), vor(lv, nsteps))
            if (if3d) then
               allocate (wor(lv, nsteps))
            else                   ! 2D
               allocate (wor(1, 1))
            end if
            if (ifto .or. ldimt > 1) allocate (tor(lt, nsteps, ldimt))
         end if
      
      !     --> Pass the initial condition for the perturbation.
         call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1),
     $   q%vx, q%vy, q%vz, q%pr, q%t)
      
         time = 0.0d+00
         do istep = 1, nsteps
      !     --> Output current info to logfile.
            if (nid == 0) write (6, "(' DIRECT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") istep, nsteps, mstep, k_dim, schur_cnt
      
      ! --> Integrate forward in time.
            call nekstab_usrchk()
            call nek_advance()
      
            if (ifstorebase .and. ifbase .and. .not. init) then !storing first time
               if (nid == 0) write (6, *) 'storing first series:', istep, '/', nsteps
               call opcopy(uor(:, istep), vor(:, istep), wor(:, istep), vx, vy, vz)
               if (ifto) call copy(tor(:, istep, 1), t(:, :, :, :, 1), nt)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(tor(:, istep, m), t(:, :, :, :, m), nt)
               end do
               end if
            elseif (ifstorebase .and. init .and. .not. ifbase) then !just moving in memory
               if (nid == 0) write (6, *) 'using stored baseflow'
               call opcopy(vx, vy, vz, uor(:, istep), vor(:, istep), wor(:, istep))
               if (ifto) call copy(t(:, :, :, :, 1), tor(:, istep, 1), nt)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(t(:, :, :, :, m), tor(:, istep, m), nt)
               end do
               end if
            end if
         end do
         if (ifstorebase .and. .not. init .and. ifbase) then
            ifbase = .false.; init = .true.
         end if
      
      !     --> Copy the solution.
         call nopcopy(f%vx, f%vy, f%vz, f%pr, f%t,
     $   vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1))
      
         return
      end subroutine forward_linearized_map
      
      !-----------------------------------------------------------------------
      
      subroutine forward_finite_difference_map(f, q)
      
      ! Approximate the linearized forward map by finite-differencing the nonlinear solver.
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         include 'ADJOINT'
      
      !integer, parameter :: findiff_order = 2
         real, dimension(findiff_order) :: coeficients
         real, dimension(findiff_order) :: amplitudes
      
         type(krylov_vector) :: q, f, pert
         type(krylov_vector) :: work
      
         real :: epsilon0, dummy
         integer :: i, m
      
         logical, save :: init
         data init/.false./
         nv = nx1*ny1*nz1*nelv
      
      ! --> Setup the Nek parameters for the finite-differences approximation.
         ifpert = .false.; ifadj = .false.
         call bcast(ifpert, lsize); call bcast(ifadj, lsize)
      
         call k_zero(f)
         call k_zero(work)
         call nopcopy(work%vx, work%vy, work%vz, work%pr, work%t, ubase, vbase, wbase, pbase, tbase)
      
         call k_norm(dummy, work)
         epsilon0 = 1e-6*dummy
      
         if (findiff_order == 2) then
            amplitudes(1) = 1; amplitudes(2) = -1
            coeficients(1) = 1; coeficients(2) = -1
            coeficients = coeficients/2.0d+00
         else if (findiff_order == 4) then
            amplitudes(1) = 1; amplitudes(2) = -1
            amplitudes(3) = 2; amplitudes(4) = -2
            coeficients(1) = 8; coeficients(2) = -8
            coeficients(3) = -1; coeficients(4) = 1
            coeficients = coeficients/12.0d+00
         end if
         amplitudes = amplitudes*epsilon0
      
      ! --> Turning-off the base flow side-by-side computation.
         ifbase = .false.
         if (uparam(01) == 3.11) ifbase = .true. ! activate Floquet
         if (uparam(01) == 3.31) ifbase = .true. ! activate Floquet for intracycle transient growth
         if (uparam(01) == 2.1 .or. uparam(01) == 2.2) then
            init = .true.          ! Use stored baseflow if ifstorebase.
            ifbase = .true.        ! activate baseflow evolution for UPO.
         end if
         if (ifstorebase .and. init) ifbase = .false. ! deactivate ifbase if baseflow stored.
      
         if (ifstorebase .and. ifbase .and. .not. init) then
            if (nid == 0) write (6, *) 'ALLOCATING ORBIT WITH NSTEPS:', nsteps
            allocate (uor(lv, nsteps), vor(lv, nsteps))
            if (if3d) then
               allocate (wor(lv, nsteps))
            else                   ! 2D
               allocate (wor(1, 1))
            end if
            if (ifto .or. ldimt > 1) allocate (tor(lt, nsteps, ldimt))
         end if
      
      !-----------------------------------------------------------------------
      !-----                                                             -----
      !-----     FINITE-DIFFERENCE APPROX. OF THE FRECHET DERIVATIVE     -----
      !-----                                                             -----
      !-----------------------------------------------------------------------
      
         do i = 1, findiff_order
      
      ! --> Scale the perturbation.
            call k_copy(pert, q)
            call k_cmult(pert, amplitudes(i))
      
      ! --> Initial condition for the each evaluation.
            call nopcopy(vx, vy, vz, pr, t, ubase, vbase, wbase, pbase, tbase)
            call nopadd2(vx, vy, vz, pr, t, pert%vx, pert%vy, pert%vz, pert%pr, pert%t)
            if (ifbf2d .and. if3d) then
               call rzero(vz, nv); if (nid == 0) write (6, *) 'Forcing vz=0'
            end if
      
      ! --> Time-integration of the nonlinear Nek5000 equations.
            time = 0.0d+00
            do istep = 1, nsteps
      !     --> Output current info to logfile.
               if (nid == 0) write (6, "(' DIRECT FD [',I1,'/',I1,']:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") i,
     $   findiff_order, istep, nsteps, mstep, k_dim, schur_cnt
      
      ! --> Nek5000 computational core.
               call nekstab_usrchk()
               call nek_advance()
      
               if (i == 1 .and. ifstorebase .and. ifbase .and. .not. init) then !storing first time
                  if (nid == 0) write (6, *) 'storing first series:', istep, '/', nsteps
                  call opcopy(uor(:, istep), vor(:, istep), wor(:, istep), vx, vy, vz)
                  if (ifto) call copy(tor(:, istep, 1), t(:, :, :, :, 1), nt)
                  if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) call copy(tor(:, istep, m), t(:, :, :, :, m), nt)
                  end do
                  end if
               elseif (i > 1 .and. ifstorebase .and. init .and. .not. ifbase) then !just moving in memory
                  if (nid == 0) write (6, *) 'using stored baseflow'
                  call opcopy(vx, vy, vz, uor(:, istep), vor(:, istep), wor(:, istep))
                  if (ifto) call copy(t(:, :, :, :, 1), tor(:, istep, 1), nt)
                  if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) call copy(t(:, :, :, :, m), tor(:, istep, m), nt)
                  end do
                  end if
               end if
            end do
            if (ifstorebase .and. .not. init .and. ifbase) then
               ifbase = .false.; init = .true.
            end if
      
      !     --> Copy the solution and compute the approximation of the Fréchet dérivative.
            call nopcopy(work%vx, work%vy, work%vz, work%pr, work%t, vx, vy, vz, pr, t)
            call k_cmult(work, coeficients(i))
            call k_add2(f, work)
      
         end do
      
      ! --> Rescale the approximate Fréchet derivative with the step size.
         call k_cmult(f, 1.0d+00/epsilon0)
      
         return
      end subroutine forward_finite_difference_map
      
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
         data init/.false./
      
         integer m
         nt = nx1*ny1*nz1*nelt
      
      !     --> Setup the parameters for the linearized solver.
         ifpert = .true.; ifadj = .true.
         call bcast(ifpert, lsize); call bcast(ifadj, lsize)
      
      !     --> Turning-off the base flow side-by-side computation. Need to change for Floquet.
         ifbase = .false.
         if (uparam(01) == 3.21) ifbase = .true. ! activate floquet
      
         if (ifstorebase .and. init) ifbase = .false.
      
         if (ifstorebase .and. ifbase .and. .not. init) then
            if (nid == 0) write (6, *) 'ALLOCATING ORBIT WITH NSTEPS:', nsteps
            allocate (uor(lv, nsteps), vor(lv, nsteps))
            if (if3d) then
               allocate (wor(lv, nsteps))
            else                   ! 2D
               allocate (wor(1, 1))
            end if
            if (ifto .or. ldimt > 1) allocate (tor(lt, nsteps, ldimt))
         end if
         if (uparam(01) == 3.31) init = .true. ! activate Floquet for intracycle transient growth (base flow already computed)
      
      !     --> Pass the initial condition for the perturbation.
         call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1),
     $   q%vx, q%vy, q%vz, q%pr, q%t)
      
         time = 0.0d+00
         do istep = 1, nsteps
      !     --> Output current info to logfile.
            if (nid == 0) write (6, "(' ADJOINT:',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')")
     $   istep, nsteps, mstep, k_dim, schur_cnt
      
      ! --> Integrate backward in time.
            call nekstab_usrchk()
            call nek_advance()
      
            if (ifstorebase .and. ifbase .and. .not. init) then !storing first time
               if (nid == 0) write (6, *) 'storing first series:', istep, '/', nsteps
               call opcopy(uor(:, istep), vor(:, istep), wor(:, istep), vx, vy, vz)
               if (ifto) call copy(tor(:, istep, 1), t(:, :, :, :, 1), nt)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(tor(:, istep, m), t(:, :, :, :, m), nt)
               end do
               end if
      
            elseif (ifstorebase .and. init .and. .not. ifbase) then !just moving in memory
               if (nid == 0) write (6, *) 'using stored baseflow'
               call opcopy(vx, vy, vz, uor(:, istep), vor(:, istep), wor(:, istep))
               if (ifto) call copy(t(:, :, :, :, 1), tor(:, istep, 1), nt)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(t(:, :, :, :, m), tor(:, istep, m), nt)
               end do
               end if
            end if
         end do
         if (ifstorebase .and. .not. init .and. ifbase) then
            ifbase = .false.; init = .true.
         end if
      
      !     --> Copy the solution.
         call nopcopy(f%vx, f%vy, f%vz, f%pr, f%t,
     $   vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1))
      
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
         call k_sub2(f, q)
         call k_cmult(f, -1.0d+00)
      
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
         if (iffindiff) then
            if (nid == 0) write (*, *) "Using the finite-difference approximation of the Fréchet derivative."
            call forward_finite_difference_map(f, q)
         else
            call forward_linearized_map(f, q)
         end if
      
      !     --> Evaluate (exp(t*L) - I) * q0.
         call k_sub2(f, q)
      
      !     ----------------------------------
      !     -----     NEWTON FOR UPO     -----
      !     ----------------------------------
      
         if (uparam(1) == 2.1) then
      
            call k_zero(bvec)
            call k_zero(btvec)
      
            call compute_bvec(bvec, fc_nwt)
            call k_cmult(bvec, q%time)
            call k_add2(f, bvec)
      
            call compute_bvec(btvec, ic_nwt)
            call k_dot(f%time, btvec, q)
      
            if (nid == 0) write (6, *) 'Newton period correction:', f%time
      
         else
      
            f%time = 0.0d+00
      
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
         ifpert = .false.; ifadj = .false.
         call bcast(ifpert, lsize); call bcast(ifadj, lsize)
      
         ifbase = .false.
         call bcast(ifbase, lsize)
      
      !     --> Pass the initial condition.
         call nopcopy(vx, vy, vz, pr, t, qbase%vx, qbase%vy, qbase%vz, qbase%pr, qbase%t)
         param(10) = qbase%time
         call prepare_linearized_solver
         call k_copy(wrk1, qbase)
      
      !     --> Single time-step to approximate the time-derivative.
         time = 0.0d+00
         do istep = 1, 1
            call nekStab_usrchk()
            call nek_advance()
         end do
         call nopcopy(wrk2%vx, wrk2%vy, wrk2%vz, wrk2%pr, wrk2%t, vx, vy, vz, pr, t)
      
      !     --> Approximate the time-derivative.
         call k_sub2(wrk2, wrk1)
         call k_copy(bvec, wrk2)
         call k_cmult(bvec, 1.0/dt)
         bvec%time = 0.0d+00
      
         return
      end subroutine compute_bvec
