      module krylov_subspace
      implicit none
      include 'SIZE'

      private
      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      type, public :: krylov_vector
      real, dimension(lt) :: vx, vy, vz
      real, dimension(lt2) :: pr
      real, dimension(lt) :: theta
      real :: time
      end type krylov_vector

      type(krylov_vector), public :: ic_nwt, fc_nwt

      contains
      end module krylov_subspace


      subroutine krylov_inner_product(alpha, p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q
      real :: alpha
      real :: glsc3
      integer :: n

      n = nx1*ny1*nz1*nelt

!     --> Kinetic energy.
      alpha = glsc3(p%vx, bm1s, q%vx, n) + glsc3(p%vy, bm1s, q%vy, n)
      if (if3d) alpha = alpha + glsc3(p%vz, bm1s, q%vz, n)

!     --> Potential energy.
      if (ifheat) alpha = alpha + glsc3(p%theta, bm1s, q%theta, n)

!     --> Time component.
      if ( uparam(3) .eq. 3.1 ) then
         alpha = alpha + p%time * q%time
      end if

      return
      end subroutine krylov_inner_product

      subroutine krylov_norm(alpha, p)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q
      real :: alpha

      call krylov_inner_product(alpha, p, p)
      alpha = dsqrt(alpha)
      return
      end subroutine krylov_norm

      subroutine krylov_normalize(p, alpha)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p
      real :: alpha, inv_alpha

!     --> Compute the user-defined norm.
      call krylov_norm(alpha, p)
      inv_alpha = 1.0D+00 / alpha

!     --> Normalize the vector.
      call krylov_cmult(p, inv_alpha)

      return
      end subroutine krylov_normalize

      subroutine krylov_cmult(p, alpha)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p
      real alpha

      p%vx = p%vx * alpha
      p%vy = p%vy * alpha
      p%pr = p%pr * alpha
      if (if3d) p%vz = p%vz * alpha
      if (ifheat) p%theta = p%theta * alpha
      p%time = p%time * alpha

      return
      end subroutine krylov_cmult

      subroutine krylov_add2(p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q

      p%vx = p%vx + q%vx
      p%vy = p%vy + q%vy
      p%pr = p%pr + q%pr

      if (if3d) p%vz = p%vz + q%vz
      if (ifheat) p%theta = p%theta + q%theta

      p%time = p%time + q%time

      return
      end subroutine krylov_add2


      subroutine krylov_sub2(p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q

      p%vx = p%vx - q%vx
      p%vy = p%vy - q%vy
      p%pr = p%pr - q%pr

      if (if3d) p%vz = p%vz - q%vz
      if (ifheat) p%theta = p%theta - q%theta

      p%time = p%time - q%time

      return
      end subroutine krylov_sub2

      subroutine krylov_zero(p)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p

      p%vx = 0.0D+00
      p%vy = 0.0D+00
      p%pr = 0.0D+00

      if (if3D) p%vz = 0.0D+00
      if (ifheat) p%theta = 0.0D+00

      p%time = 0.0D+00

      return
      end subroutine krylov_zero

      subroutine krylov_copy(p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q

      p%vx = q%vx
      p%vy = q%vy
      p%pr = q%pr

      if (if3D) p%vz = q%vz
      if (ifheat) p%theta = q%theta

      p%time = q%time

      return
      end subroutine krylov_copy

      subroutine krylov_matmul(dq, Q, yvec, k)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      integer :: i, j, k
      type(krylov_vector) :: dq
      type(krylov_vector), dimension(k) :: Q
      real, dimension(k) :: yvec

      real, dimension(lt, k) :: qx, qy, qz, qt
      real, dimension(lt2, k) :: qp
      real, dimension(k) :: time_comp

      do i = 1, k
         qx(:, i) = Q(i)%vx
         qy(:, i) = Q(i)%vy
         if (if3D) qz(:, i) = Q(i)%vz
         qp(:, i) = Q(i)%pr
         if (ifheat) qt(:, i) = Q(i)%theta
         time_comp(i) = Q(i)%time
      enddo

      dq%vx = matmul(qx, yvec)
      dq%vy = matmul(qy, yvec)
      dq%vz = matmul(qy, yvec)
      dq%pr = matmul(qp, yvec)
      dq%theta = matmul(qt, yvec)
      dq%time = dot_product(time_comp, yvec)

      return
      end subroutine krylov_matmul
