!-----------------------------------------------------------------------





      subroutine newton_krylov(qx, qy, qz, qp, qt)

!     Implementation of a simple Newton-Krylov solver for fixed point
!     computation using a time-stepper formulation. The resolution of
!     the linear system for each Newton iteration is sovled by means
!     of GMRES.
!     
!     INPUTS
!     ------
!     
!     OUTPUTS
!     -------
!     
!     
!     Last edited : March 26th by JC Loiseau.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

!     ----- Right-hand side vector for Newton solver.
      real, dimension(lt) :: f_x, f_y, f_z, f_t
      real, dimension(lt2) :: f_p

!     ----- Current estimate of the solution.
      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

!     ----- Newton correction obtained from GMRES.
      real, dimension(lt) :: dqx, dqy, dqz, dqt
      real, dimension(lt2) :: dqp

!     ----- Miscellaneous
      real :: tol, residual
      integer :: maxiter_newton, maxiter_gmres
      integer :: i, j, k, n     !, n2

      n = nx1*ny1*nz1*nelt      !; n2 = nx2*ny2*nz2*nelt
      maxiter_newton = 10 ; maxiter_gmres = 10
      tol = max(param(21), param(22))

!     --> Initialize arrays.
      call noprzero(f_x, f_y, f_z, f_p, f_t)

!     --> Newton iteration.
      newton : do i = 1, maxiter_newton

!     --> Setup/Update the nek-parameters for the Newton solver.
      call newton_krylov_prepare

!     --> Outpost current estimate of the solution
      time = i ; call outpost(qx, qy, qz, qp, qt, "nwt")

!     --> Compute rhs of Newton iteration f(q).
      call forward_map(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)

!     --> Check residual || f(q) ||
      call norm(f_x, f_y, f_z, f_p, f_t, residual) ; residual = residual ** 2

      if (nid .eq. 0) write(6, *) "NEWTON --- Iteration :", i, " residual :", residual
      if (residual .lt. tol) exit newton

      if (nid.eq.0) write(*, *) "LINEAR SOLVER"
!     --> Solve the linear system.
      call ts_gmres(f_x, f_y, f_z, f_p, f_t, dqx, dqy, dqz, dqp, dqt, maxiter_gmres, k_dim)

!     --> Update Newton solution.
      call nopsub2(qx, qy, qz, qp, qt, dqx, dqy, dqz, dqp, dqt)

      enddo newton

      return
      end subroutine newton_krylov


!-----------------------------------------------------------------------



      subroutine ts_gmres(rhs_x, rhs_y, rhs_z, rhs_p, rhs_t, sol_x, sol_y, sol_z, sol_p, sol_t, maxiter, ksize)

!     Implementation of simple GMRES to be part of the Newton-Krylov solver
!     for fixed point computation. The rank of the Krylov subspace is set as the user parameter k_dim.
!     
!     INPUT
!     -----
!     
!     rhs_x, rhs_y, rhs_z, rhs_t : nek arrays of size (lt).
!     Arrays containing the right-hand side of the linear problem to be solved.
!     
!     rhs_p : nek array of size (lt2)
!     Array containing the right-hand side of the linear problem to be solved (pressure component).
!     
!     maxiter : integer
!     Maximum number of restarts for the GMRES computation.
!     
!     ksize : integer
!     Dimension of the Krylov subspace.
!     
!     RETURNS
!     -------
!     
!     sol_x, sol_y, sol_z, sol_t : nek arrays of size (lt).
!     Arrays containing the solution of the linear problem.
!     
!     sol_p : nek array of size (lt2).
!     Array containing the solution of the linear problem (pressure component).
!     
!     
!     NOTE : This is a plain implementation of GMRES following the algorithm given in
!     Y. Saad. Iterative methods for sparse linear systems. Section 6.5 GMRES, alg. 6.9
!     
!     Last Edit : March 26th 2021 by JC Loiseau.
!     

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      integer :: ksize

!     ----- Right-hand side vector of A x = b -----
      real, dimension(lt) :: rhs_x, rhs_y, rhs_z, rhs_t
      real, dimension(lt2) :: rhs_p

!     ----- GMRES solution vector.
      real, dimension(lt) :: sol_x, sol_y, sol_z, sol_t
      real, dimension(lt2) :: sol_p

!     ----- Krylov basis for the Arnoldi factorization.
      real, allocatable, dimension(:, :) :: qx, qy, qz, qt
      real, allocatable, dimension(:, :) :: qp

!     ----- Upper Hessenberg matrix.
      real, allocatable, dimension(:, :) :: H
      real, allocatable, dimension(:) :: yvec, evec

!     ----- Miscellaneous.
      integer :: i, j, k, maxiter !, n, n2
      real :: beta, tol

!     Initialize arrays.
      tol = max(param(21), param(22))

!     ----- Allocate arrays -----
      allocate(qx(lt, ksize+1), qy(lt, ksize+1), qz(lt, ksize+1), qt(lt, ksize+1), qp(lt2, ksize+1))
      allocate(H(ksize+1, ksize), yvec(ksize), evec(ksize+1))

      qx = 0.0D+00 ; qy = 0.0D+00 ; qz = 0.0D+00 ; qp = 0.0D+00 ; qt = 0.0D+00
      sol_x = 0.0D+00 ; sol_y = 0.0D+00 ; sol_z = 0.0D+00 ; sol_p = 0.0D+00 ; sol_t = 0.0D+00
      H = 0.0D+00 ; yvec = 0.0D+00 ; evec = 0.0D+00

      call nopcopy(qx(:, 1),qy(:, 1),qz(:, 1),qp(:, 1),qt(:, 1), rhs_x,rhs_y,rhs_z,rhs_p,rhs_t)
      call normalize(qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), beta)

      gmres : do i = 1, maxiter

!     --> Zero-out stuff.
      H = 0.0D+00 ; yvec = 0.0D+00
      evec = 0.0D+00 ; evec(1) = beta

      qx(:, 2:k_dim+1) = 0.0D+00 ; qy(:, 2:k_dim+1) = 0.0D+00 ; qz(:, 2:k_dim+1) = 0.0D+00
      qp(:, 2:k_dim+1) = 0.0D+00 ; qt(:, 2:k_dim+1) = 0.0D+00

      arnoldi : do k = 1, k_dim

!     --> Arnoldi factorization.
      call arnoldi_factorization(qx, qy, qz, qp, qt, H, k, k, ksize)

!     --> Least-squares problem.
      call lstsq(H(1:k+1, 1:k), evec(1:k+1), yvec(1:k), k+1, k)

!     --> Compute residual.
      beta = norm2(evec(1:k+1) - matmul(H(1:k+1, 1:k), yvec(1:k)))
      if (nid.eq.0) write(6, *) "Inner iteration residual : ", beta**2
      if (beta**2 .lt. tol) exit arnoldi

      enddo arnoldi

!     --> Update solution.
      sol_x = sol_x + matmul(qx(:, 1:k), yvec(1:k))
      sol_y = sol_y + matmul(qy(:, 1:k), yvec(1:k))
      if (if3d) sol_z = sol_z + matmul(qz(:, 1:k), yvec(1:k))
      if (ifpo) sol_p = sol_p + matmul(qp(:, 1:k), yvec(1:k))
      if (ifto) sol_t = sol_t + matmul(qt(:, 1:k), yvec(1:k))

!     --> Recompute residual for sanity check and initialize new Krylov seed if needed.

      call nopcopy(qx(:, 1),qy(:, 1),qz(:, 1),qp(:, 1),qt(:, 1), sol_x,sol_y,sol_z,sol_p,sol_t)
      call initialize_gmres_vector(beta, qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

      if (nid.EQ.0) write(6, *) "GMRES --- Iteration : ", i, " residual : ", beta**2
      if (beta**2 .lt. tol) exit gmres

      enddo gmres

!     ----- Deallocate arrays -----
      deallocate(qx, qy, qz, qt, qp)
      deallocate(H, yvec, evec)

      return
      end subroutine ts_gmres




!-----------------------------------------------------------------------





      subroutine initialize_gmres_vector(beta, qx, qy, qz, qp, qt, rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv

!     ----- Right-hand side vector of A x = b (input) -----
      real, dimension(lt) :: rhs_x, rhs_y, rhs_z, rhs_t
      real, dimension(lt2) :: rhs_p

!     ----- Seed for the GMRES Krylov subspace (output) -----
      real, dimension(lt) :: qx, qy, qz, qt, f_x, f_y, f_z, f_t
      real, dimension(lt2) :: qp, f_p
      real :: beta

!     --> Initial Krylov vector.
      call matvec(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)
      call nopsub2(f_x,f_y,f_z,f_p,f_t, rhs_x,rhs_y,rhs_z,rhs_p,rhs_t)

!     --> Normalize the starting vector.
      call normalize(f_x, f_y, f_z, f_p, f_t, beta)
      call nopcopy(qx,qy,qz,qp,qt, f_x,f_y,f_z,f_p,f_t)

      return
      end subroutine initialize_gmres_vector





!-----------------------------------------------------------------------
      subroutine forward_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv

!     ----- Initial condition for the forward simulation.
      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

!     ----- Right-hand side of the Newton.
      real, dimension(lt) :: fx, fy, fz, ft
      real, dimension(lt2) :: fp

!     --> Fixed point/UPO computation.
      if (uparam(01) .eq. 1) call nonlinear_forward_map(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)

      return
      end subroutine forward_map



      subroutine nonlinear_forward_map(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv

!     ----- Initial condition for the forward simulation.
      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

!     ----- Right-hand side of the Newton.
      real, dimension(lt) :: f_x, f_y, f_z, f_t
      real, dimension(lt2) :: f_p

      integer n

      n = nx1*ny1*nz1*nelt

!     --> Copy the initial condition to Nek.
      call nopcopy(vx, vy, vz, pr, t(1,1,1,1,1), qx, qy, qz, qp, qt)

!     --> Turn-off the linearized solver.
      ifpert = .false. ; call bcast(ifpert, lsize)

!     --> Run the simulation forward.
      time = 0.0d0
      do istep = 1, nsteps
         call nekStab_usrchk()
         call nek_advance()
      enddo

!     --> Compute the right hand side of the time-stepper Newton.
      call nopcopy(f_x,f_y,f_z,f_p,f_t, vx,vy,vz,pr,t(1,1,1,1,1))
      call nopsub2(f_x,f_y,f_z,f_p,f_t, qx,qy,qz,qp,qt)
      call nopchsign(f_x,f_y,f_z,f_p,f_t)

!     --> Pass current guess as base flow for the linearized calculation.
      call opcopy(ubase, vbase, wbase, qx, qy, qz) ; call copy(tbase, qt, n)

      return
      end subroutine nonlinear_forward_map


      subroutine newton_krylov_prepare

!     This
!     
!     INPUT
!     -----
!     
!     RETURNS
!     -------
!     
!     Last edit : March 26th 2021 by RAS Frantz.

      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

!     forcing npert to unity
      if(param(31).gt.1)then
         write(6,*)'nekStab not ready for npert>1 -- jp loops NOT implemented! STOPPING!'
         call nek_end
      endif
      param(31) = 1 ; npert = param(31)

!     force OIFS deactivation
      if(ifchar.and.nid.eq.0)write(6,*)'OIFS not working with linearized and adjoint -> turning OFF!'
      ifchar = .false.
      call bcast(ifchar, lsize)

!     enforce CFL target for EXTk
      if( param(26).gt.0.5 )then
         if(nid.eq.0)write(6,*)'reducing target CFL to 0.5!'
         param(26)=0.50d0 ; ctarg = param(26)
      endif

      if(param(10).gt.0)then
!     if param(10)=endTime=0 -> param(11) = numSteps
         call compute_cfl(dt,vx,vy,vz,1.0) ! dt=1 ! vx at this point is base flow
         dt = ctarg/dt
         nsteps = ceiling(param(10)/dt)
         dt = param(10)/nsteps
         if(nid.eq.0)write(6,*)'endTime specified! computing CFL from BASE FLOW!'
         if(nid.eq.0)write(6,*)' computing timeStep dt=',dt
         if(nid.eq.0)write(6,*)' computing numSteps=',nsteps
         if(nid.eq.0)write(6,*)' sampling period =',nsteps*dt
         param(12) = dt
      endif

!     deactivate variable time step! !freeze dt
      param(12) = -abs(param(12))

!     broadcast all parameters to processors
      call bcast(param,200*wdsize)

      return
      end subroutine newton_krylov_prepare
