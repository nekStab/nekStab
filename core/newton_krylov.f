!-----------------------------------------------------------------------





      subroutine newton_krylov()

!     Implementation of a simple Newton-Krylov solver for fixed point
!     computation using a time-stepper formulation. The resolution of
!     the linear system for each Newton iteration is sovled by means
!     of GMRES.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

! ----- Right-hand side vector for Newton solver.
      real, dimension(lt) :: f_x, f_y, f_z, f_t
      real, dimension(lt2) :: f_p

! -----
      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      real, dimension(lt) :: dqx, dqy, dqz, dqt
      real, dimension(lt2) :: dqp

! ----- Miscellaneous
      real :: tol, residual
      integer :: maxiter_newton, maxiter_gmres
      integer :: i, j, k, n, n2
      logical :: converged

      n = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

      maxiter_newton = 10
      maxiter_gmres = 10
      converged = .false.
      tol = 1e-8

! --> Initialize arrays.
      call oprzero(f_x, f_y, f_z)
      call rzero(f_p, n2)
      call rzero(f_t, n)

! --> Copy initial guess for the Newton solver.
      call opcopy(qx, qy, qz, vx, vy, vz)
      call copy(qt, t, n)
      call copy(qp, pr, n2)

! --> Newton iteration.
      i = 0
      do while ((.not. converged) .and. (i < maxiter_newton))

         time = i
         call outpost(qx, qy, qz, qp, qt, "nwt")

! --> Compute rhs of Newton iteration f(q).
         call forward_map(f_x, f_y, f_z, f_p, f_t, qx, qy, qz, qp, qt)

! --> Check residual || f(q) ||
         call norm(f_x, f_y, f_z, f_p, f_t, residual)
         residual = residual ** 2

         if (nid .eq. 0) then
            write(6, *) "NEWTON --- Iteration :", i, " residual :", residual
         endif

         if (residual .lt. max(param(21),param(22))) then
            converged = .true.
         else
            write(*, *) "LINEAR SOLVER"
! --> Copy the current guess into the base flow.
            call opcopy(ubase, vbase, wbase, qx, qy, qz)
            call copy(tbase, t, n)

! --> Zero-initial guess for the solution of the linear system.
            call oprzero(dqx, dqy, dqz)
            call rzero(dqp, n2)
            call rzero(dqt, n)

! --> Solve the linear system.
            param(12) = -abs(param(12)) !freeze dt
            param(31) = 1 ; npert = param(31)
            call bcast(param,200*wdsize) !broadcast all parameters to processors

            call time_stepper_gmres(f_x, f_y, f_z, f_p, f_t, dqx, dqy, dqz, dqp, dqt, maxiter_gmres)
            ifpert = .false.
            call bcast(ifpert, lsize)

! --> Update Newton solution.
            call opsub2(qx, qy, qz, dqx, dqy, dqz)
            call sub2(qp, dqp, n2)
            call sub2(qt, dqt, n)

         end if

! -->
         i = i + 1
      enddo

      call outpost(qx, qy, qz, qp, qt, "NBF")

      return
      end subroutine newton_krylov




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
        !     Last Edit : March 22nd 2021 by JC Loiseau.
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
        logical :: converged
        real :: beta, residual, tol, rhs_norm, ls_residual

        n = nx1*ny1*nz1*nelt
        n2 = nx2*ny2*nz2*nelt

        uparam(01) = 3
        ifpert = .true.
        call bcast(ifpert, lsize)

        ! Initialize arrays.
        tol = max(param(21),param(22))
        converged = .false.

        qx = 0.0D+00 ; qy = 0.0D+00 ; qz = 0.0D+00 ; qp = 0.0D+00 ; qt = 0.0D+00
        H = 0.0D+00 ; yvec = 0.0D+00 ; evec = 0.0D+00

        call opcopy(qx(:, 1), qy(:, 1), qz(:, 1), rhs_x, rhs_y, rhs_z)
        call copy(qp(:, 1), rhs_p, n2) ; call copy(qt(:, 1), rhs_t, n)

        call normalize(qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), beta)

        i = 0
        do while ((.not. converged) .and. (i < maxiter))

          ! -->
          H = 0.0D+00 ; yvec = 0.0D+00
          evec = 0.0D+00 ; evec(1) = beta

          arnoldi : do k = 1, k_dim
            ! --> Arnoldi factorization.
            call arnoldi_factorization(qx, qy, qz, qp, qt, H, k, k)
            ! --> Solve the least-squares problem.
            call lstsq(H(1:k+1, 1:k), evec(1:k+1), yvec(1:k), k+1, k)

            ! --> Compute residual.
            beta = norm2(evec(1:k+1) - matmul(H(1:k+1, 1:k), yvec(1:k)))
            if (nid.eq.0) write(*, *) "Inner iteration residual : ", beta**2

            if (beta**2 .lt. tol) then
              converged = .true.
              exit arnoldi
            endif

          enddo arnoldi

          ! --> Update solution.
          sol_x = sol_x + matmul(qx(:, 1:k_dim), yvec)
          sol_y = sol_y + matmul(qy(:, 1:k_dim), yvec)
          sol_z = sol_z + matmul(qz(:, 1:k_dim), yvec)
          sol_p = sol_p + matmul(qp(:, 1:k_dim), yvec)
          sol_t = sol_t + matmul(qt(:, 1:k_dim), yvec)

          ! --> Recompute residual.
          call opcopy(qx(:, 1), qy(:, 1), qz(:, 1), sol_x, sol_y, sol_z)
          call copy(qp(:, 1), sol_p, n2) ; call copy(qt(:, 1), sol_t, n)
          call initialize_gmres_vector(beta, qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

          if (nid.EQ.0) write(*, *) "GMRES --- Iteration : ", i, " residual : ", beta**2
          if (beta**2 .lt. tol) converged = .true.
          i = i+1

        enddo

        uparam(01) = 1
        ifpert = .false.

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

! Initialize stuff.
      n = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

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


      n = nx1*ny1*nz1*nelt
      n2 = nx2*ny2*nz2*nelt

! --> Copy the initial condition to Nek.
      call opcopy(vx, vy, vz, qx, qy, qz)
      call copy(t, qt, n)
      call copy(pr, qp, n2)

! --> Turn-off the linearized solver.
      ifpert = .false.
      call bcast(ifpert, lsize)

! --> Run the simulation forward.
      do istep = 1, nsteps
         call nekStab_chk()
         call nek_advance()
      enddo

! --> Compute the right hand side of the time-stepper Newton.
      call opcopy(f_x, f_y, f_z, vx, vy, vz)
      call copy(f_p, pr, n2)
      call copy(f_t, t, n)

      call opsub2(f_x, f_y, f_z, qx, qy, qz)
      call sub2(f_p, qp, n2)
      call sub2(f_t, qt, n)

      call chsign(f_x, n)
      call chsign(f_y, n)
      call chsign(f_z, n)
      call chsign(f_p, n2)
      call chsign(f_t, n)

      return
      end subroutine forward_map
