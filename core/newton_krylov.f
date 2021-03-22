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

      return
      end subroutine newton_krylov




!-----------------------------------------------------------------------




      subroutine gmres(rhs_x, rhs_y, rhs_z, rhs_p, rhs_t, sol_x, sol_y, sol_z, sol_p, sol_t, maxiter)

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
        real dimension(lt2, k_dim+1) :: qp

        !     ----- Upper Hessenberg matrix.
        real, dimension(k_dim+1, k_dim) :: H
        real, dimension(k_dim) :: yvec
        real, dimension(k_dim+1) :: evec

        !     ----- Miscellaneous.
        integer :: i, j, k, maxiter, n
        logical :: converged
        real :: beta, residual, tol

        n = nx1*ny1*nz1*nelv
        n2 = nx2*ny2*nz2*nelv

        ! Initialize arrays.
        ifres = .false. ! Turns off the checkpointing in Arnoldi factorization.
        converged = .false.
        tol = 1e-8

        do i = 1, k_dim+1
          call oprzero(qx(:, i), qy(:, i), qz(:, i))
          call rzero(qp(:, i), n2)
          call rzero(qt(:, i), n)
        enddo

        i = 0
        do while ((.not. converged) .and. (i < maxiter))

          ! --> Initialize arrays.
          H(:, :) = 0.0D+00
          yvec = 0.0D+00
          evec = 0.0D+00


          !         --> Generate seed for the Krylov subspace.
          call opcopy(qx(:, 1), qy(:, 1), qz(:, 1), sol_x, sol_y, solz)
          call copy(qp, sol_p, n2)
          call copy(qt, sol_t, n)
          call initialize_gmres_vector(beta, qx(:, 1), qy(:, 1), qz(:, 1), qp(:, 1), qt(:, 1), rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

          evec(1) = beta

          !           --> Fill-in the Krylov subspace.
          call arnoldi_factorization(qx, qy, qz, qp, qt, H, 1, k_dim)

          !        --> Solve the least-squares problem.
          call lstsq(H, evec, yvec, k_dim+1, k_dim)

          !         --> Compute the residual.
          residual = sum(abs(matmul(H, yvec) - evec)**2)

          ! --> Check convergence.
          if (residual .lt. tol) converged = .true.

          ! --> Update the solution vector.
          sol_x = sol_x + matmul(qx(:, 1:k_dim), yvec)
          sol_y = sol_y + matmul(qy(:, 1:k_dim), yvec)
          sol_z = sol_z + matmul(qz(:, 1:k_dim), yvec)
          sol_p = sol_p + matmul(qp(:, 1:k_dim), yvec)
          sol_t = sol_t + matmul(qt(:, 1:k_dim), yvec)


          ! --> Update counter of GMRES restarts.
          i = i + 1

        enddo

        return
      end subroutine gmres




      !-----------------------------------------------------------------------





      subroutine initialize_gmres_vector(alpha, qx, qy, qz, qp, qt, rhs_x, rhs_y, rhs_z, rhs_p, rhs_t)

        implicit none
        include 'SIZE'
        include 'TOTAL'

        integer, parameter :: lt = lx1*ly1*lz1*lelv
        integer, parameter :: lt2 = lx2*ly2*lz2*lelv

        ! ----- Right-hand side vector of A x = b (input) -----
        real, dimension(lt) :: rhs_x, rhs_y, rhs_z, rhs_t
        real, dimension(lt2) :: rhs_p

        ! ----- Seed for the GMRES Krylov subspace (output) -----
        real, dimension(lt) :: qx, qy, qz, qp
        real, dimension(lt2) :: qp
        real :: alpha

        ! Initialize stuff.
        n = nx1*ny1*nz1*nelv
        n2 = nx2*ny2*nz2*nelv

        ! --> Initial Krylov vector.
        call matrix_vector_product(qx, qy, qz, qp, qt, qx, qy, qz, qp, qt)
        call opsub2(qx, qy, qz, rhs_x, rhs_y, rhs_z)
        call sub2(qp, rhs_p, n2)
        call sub2(qt, rhs_t, n)

        ! --> Change the sign to r = b - Ax.
        call chsign(qx, n)
        call chsign(qy, n)
        call chsign(qz, n)
        call chsign(qt, n)
        call chsign(qp, n2)

        ! --> Normalize the starting vector.
        call normalize(qx, qy, qz, qp, qt, alpha)

        return
      end subroutine initialize_gmres_vector
