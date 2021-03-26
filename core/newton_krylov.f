!-----------------------------------------------------------------------





      subroutine newton_krylov()

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

! ----- Right-hand side vector for Newton solver.
      real, dimension(lt) :: f_x, f_y, f_z, f_t
      real, dimension(lt2) :: f_p

! ----- Current estimate of the solution.
      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      ! ----- Newton correction obtained from GMRES.
      real, dimension(lt) :: dqx, dqy, dqz, dqt
      real, dimension(lt2) :: dqp

! ----- Miscellaneous
      real :: tol, residual
      integer :: maxiter_newton, maxiter_gmres
      integer :: i, j, k, n, n2

      n = nx1*ny1*nz1*nelt ; n2 = nx2*ny2*nz2*nelt
      maxiter_newton = 10 ; maxiter_gmres = 10
      tol = max(param(21), param(22))

! --> Initialize arrays.
      call oprzero(f_x, f_y, f_z)
      call rzero(f_p, n2) ; call rzero(f_t, n)

! --> Copy initial guess for the Newton solver.
      call opcopy(qx, qy, qz, vx, vy, vz)
      call copy(qt, t, n) ; call copy(qp, pr, n2)

! --> Newton iteration.
      newton : do i = 1, maxiter_newton

        ! --> Outpost current estimate of the solution
         time = i ; call outpost(qx, qy, qz, qp, qt, "nwt")

         ! --> Compute rhs of Newton iteration f(q).
         call forward_map(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)

         ! --> Check residual || f(q) ||
         call norm(f_x, f_y, f_z, f_p, f_t, residual) ; residual = residual ** 2

         if (nid .eq. 0) write(6, *) "NEWTON --- Iteration :", i, " residual :", residual
         if (residual .lt. tol) exit newton

         if (nid.eq.0) write(*, *) "LINEAR SOLVER"
         ! --> Copy the current guess into the base flow.
         call opcopy(ubase, vbase, wbase, qx, qy, qz) ; call copy(tbase, t, n)

         ! --> Solve the linear system.
         call time_stepper_gmres(f_x, f_y, f_z, f_p, f_t, dqx, dqy, dqz, dqp, dqt, maxiter_gmres)

         ! --> Update Newton solution.
         call opsub2(qx, qy, qz, dqx, dqy, dqz) ; call sub2(qp, dqp, n2) ; call sub2(qt, dqt, n)

       enddo newton

       ! --> Outpost solution.
      call outpost(qx, qy, qz, qp, qt, "NBF")

      return
      end subroutine newton_krylov


!-----------------------------------------------------------------------



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
      ! if param(10)=endTime=0 -> param(11) = numSteps
        call compute_cfl(dt,vx,vy,vz,1.0) ! dt=1
        dt = ctarg/dt
        nsteps = ceiling(param(10)/dt)
        dt = param(10)/nsteps
                  if(nid.eq.0)write(6,*)'endTime specified! computing CFL from initial condition!'
        if(nid.eq.0)write(6,*)' computing timeStep dt=',dt
        if(nid.eq.0)write(6,*)' computing numSteps=',nsteps
        if(nid.eq.0)write(6,*)' sampling period =',nsteps*dt
        param(12) = dt
      endif

      ! deactivate variable time step! !freeze dt
      param(12) = -abs(param(12))

      ! broadcast all parameters to processors
      call bcast(param,200*wdsize) 

      return 
      end



!-----------------------------------------------------------------------




      subroutine time_stepper_gmres(rhs_x, rhs_y, rhs_z, rhs_p, rhs_t, sol_x, sol_y, sol_z, sol_p, sol_t, maxiter)

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

        !     ----- Right-hand side vector of A x = b -----
        real, dimension(lt) :: rhs_x, rhs_y, rhs_z, rhs_t
        real, dimension(lt2) :: rhs_p

        !     ----- GMRES solution vector.
        real, dimension(lt) :: sol_x, sol_y, sol_z, sol_t
        real, dimension(lt2) :: sol_p

        !     ----- Krylov basis for the Arnoldi factorization.
        real, dimension(lt, k_dim+1) :: qx, qy, qz, qt
        real, dimension(lt2, k_dim+1) :: qp

        !     ----- Upper Hessenberg matrix.
        real, dimension(k_dim+1, k_dim) :: H
        real, dimension(k_dim) :: yvec
        real, dimension(k_dim+1) :: evec

        !     ----- Miscellaneous.
        integer :: i, j, k, maxiter, n, n2
        real :: beta, tol

        n = nx1*ny1*nz1*nelt ; n2 = nx2*ny2*nz2*nelt
        uparam(01) = 3 ; ifpert = .true. ; call bcast(ifpert, lsize)

        ! Initialize arrays.
        tol = max(param(21), param(22))

        qx = 0.0D+00 ; qy = 0.0D+00 ; qz = 0.0D+00 ; qp = 0.0D+00 ; qt = 0.0D+00
        sol_x = 0.0D+00 ; sol_y = 0.0D+00 ; sol_z = 0.0D+00 ; sol_p = 0.0D+00 ; sol_t = 0.0D+00
        H = 0.0D+00 ; yvec = 0.0D+00 ; evec = 0.0D+00

        call opcopy(qx(:, 1), qy(:, 1), qz(:, 1), rhs_x, rhs_y, rhs_z)
        call copy(qp(:, 1), rhs_p, n2) ; call copy(qt(:, 1), rhs_t, n)
        call normalize(qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), beta)

        gmres : do i = 1, maxiter

          ! --> Zero-out stuff.
          H = 0.0D+00 ; yvec = 0.0D+00
          evec = 0.0D+00 ; evec(1) = beta

          qx(:, 2:k_dim+1) = 0.0D+00 ; qy(:, 2:k_dim+1) = 0.0D+00 ; qz(:, 2:k_dim+1) = 0.0D+00
          qp(:, 2:k_dim+1) = 0.0D+00 ; qt(:, 2:k_dim+1) = 0.0D+00

          arnoldi : do k = 1, k_dim
            ! --> Arnoldi factorization.
            call arnoldi_factorization(qx, qy, qz, qp, qt, H, k, k)

            ! --> Least-squares problem.
            call lstsq(H(1:k+1, 1:k), evec(1:k+1), yvec(1:k), k+1, k)

            ! --> Compute residual.
            beta = norm2(evec(1:k+1) - matmul(H(1:k+1, 1:k), yvec(1:k)))
            if (nid.eq.0) write(*, *) "Inner iteration residual : ", beta**2
            if (beta**2 .lt. tol) exit arnoldi

          enddo arnoldi

          ! --> Update solution.
          sol_x = sol_x + matmul(qx(:, 1:k), yvec(1:k))
          sol_y = sol_y + matmul(qy(:, 1:k), yvec(1:k))
          sol_z = sol_z + matmul(qz(:, 1:k), yvec(1:k))
          sol_p = sol_p + matmul(qp(:, 1:k), yvec(1:k))
          sol_t = sol_t + matmul(qt(:, 1:k), yvec(1:k))

          ! --> Recompute residual for sanity check and initialize new Krylov seed if needed.
          call opcopy(qx(:, 1), qy(:, 1), qz(:, 1), sol_x, sol_y, sol_z)
          call copy(qp(:, 1), sol_p, n2) ; call copy(qt(:, 1), sol_t, n)
          call initialize_gmres_vector(beta, qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

          if (nid.EQ.0) write(*, *) "GMRES --- Iteration : ", i, " residual : ", beta**2
          if (beta**2 .lt. tol) exit gmres

        enddo gmres

        uparam(01) = 1 ; ifpert = .false. ; call bcast(ifpert, lsize)

        return
      end subroutine time_stepper_gmres




!-----------------------------------------------------------------------





      subroutine initialize_gmres_vector(beta, qx, qy, qz, qp, qt, rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv

! ----- Right-hand side vector of A x = b (input) -----
      real, dimension(lt) :: rhs_x, rhs_y, rhs_z, rhs_t
      real, dimension(lt2) :: rhs_p

! ----- Seed for the GMRES Krylov subspace (output) -----
      real, dimension(lt) :: qx, qy, qz, qt, f_x, f_y, f_z, f_t
      real, dimension(lt2) :: qp, f_p
      real :: beta
      integer :: n, n2

      n = nx1*ny1*nz1*nelt ; n2 = nx2*ny2*nz2*nelt

! --> Initial Krylov vector.
      call matrix_vector_product(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)
      f_x = rhs_x - f_x ; f_y = rhs_y - f_y ; f_z = rhs_z - f_z
      f_p = rhs_p - f_p ; f_t = rhs_t - f_t

! --> Normalize the starting vector.
      call normalize(f_x, f_y, f_z, f_p, f_t, beta)
      call opcopy(qx, qy, qz, f_x, f_y, f_z)
      call copy(qp, f_p, n2) ; call copy(qt, f_t, n)

      return
      end subroutine initialize_gmres_vector





!-----------------------------------------------------------------------




      subroutine forward_map(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv

! ----- Initial condition for the forward simulation.
      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

! ----- Right-hand side of the Newton.
      real, dimension(lt) :: f_x, f_y, f_z, f_t
      real, dimension(lt2) :: f_p

      integer :: n, n2

      n = nx1*ny1*nz1*nelt ; n2 = nx2*ny2*nz2*nelt

! --> Copy the initial condition to Nek.
      call opcopy(vx, vy, vz, qx, qy, qz)
      call copy(t, qt, n) ; call copy(pr, qp, n2)

! --> Turn-off the linearized solver.
      ifpert = .false. ; call bcast(ifpert, lsize)

! --> Ensuring sampling-period consistency.
      call compute_cfl(dt,vx,vy,vz,1.0); dt = ctarg/dt
      nsteps = ceiling(param(10)/dt); dt = param(10)/nsteps
      param(12) = -abs(dt)
      call bcast(param,200*wdsize) 

! --> Run the simulation forward.
      do istep = 1, nsteps
         call nekStab_usrchk()
         call nek_advance()
      enddo

! --> Compute the right hand side of the time-stepper Newton.
      call opcopy(f_x, f_y, f_z, vx, vy, vz)
      call copy(f_p, pr, n2) ; call copy(f_t, t, n)

      call opsub2(f_x, f_y, f_z, qx, qy, qz)
      call sub2(f_p, qp, n2) ;  call sub2(f_t, qt, n)

      call chsign(f_x, n) ; call chsign(f_y, n) ; call chsign(f_z, n)
      call chsign(f_p, n2) ; call chsign(f_t, n)

      return
      end subroutine forward_map
