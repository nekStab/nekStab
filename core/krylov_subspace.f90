      module krylov_subspace
         implicit none
         include 'SIZE'
      
         private
      
         integer, public, parameter :: lv = lx1*ly1*lz1*lelv
         integer, public, parameter :: lt = lx1*ly1*lz1*lelt
         integer, public, parameter :: lp = lx2*ly2*lz2*lelv ! lp is used for norm -> lp
         integer, save, public :: nv, nt, n2 ! np conflits with mpi variable -> n2
      ! we use lp and n2 -> for pressure fields
         type, public :: krylov_vector
            real, dimension(lv) :: vx, vy, vz
            real, dimension(lp) :: pr
            real, dimension(lv, ldimt) :: t
            real :: time
         end type krylov_vector
      
         type(krylov_vector), save, public :: ic_nwt, fc_nwt
         real, save, allocatable, dimension(:, :), public :: uor, vor, wor
         real, save, allocatable, dimension(:, :, :), public :: tor
      
      contains
      end module krylov_subspace
      
      subroutine k_dot(alpha, p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(in) :: p, q
         real, intent(out) :: alpha
         real :: glsc3
         integer m
      
         nv = nx1*ny1*nz1*nelv
         nt = nx1*ny1*nz1*nelt
      
      !     --> Kinetic energy.
         alpha = glsc3(p%vx, bm1s, q%vx, nv) + glsc3(p%vy, bm1s, q%vy, nv)
         if (if3d) alpha = alpha + glsc3(p%vz, bm1s, q%vz, nv)
      
      !     --> Potential energy.
         if (ifto) alpha = alpha + glsc3(p%t(:, 1), bm1s, q%t(:, 1), nt)
         if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) alpha = alpha + glsc3(p%t(:, m), bm1s, q%t(:, m), nt)
         end do
         end if
      
      !     --> Time component.
         if (uparam(1) == 2.1) then
            alpha = alpha + p%time*q%time
         end if
      
      !     --> Check integrity.
      !  if (isnan(alpha)) call nek_end
      
         return
      end subroutine k_dot
      
      subroutine k_norm(alpha, p)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(in) :: p
         real, intent(out) :: alpha
      
         call k_dot(alpha, p, p)
         alpha = sqrt(alpha)
         return
      end subroutine k_norm
      
      subroutine k_normalize(p, alpha)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         real, intent(out) :: alpha
         real :: inv_alpha
      
      !     --> Compute the user-defined norm.
         call k_norm(alpha, p)
         inv_alpha = 1.0d+00/alpha
      
      !     --> Normalize the vector.
         call k_cmult(p, inv_alpha)
      
         return
      end subroutine k_normalize
      
      subroutine k_cmult(p, c)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector) :: p
         real c
         call nopcmult(p%vx, p%vy, p%vz, p%pr, p%t, c)
         p%time = p%time*c
         return
      end subroutine k_cmult
      
      subroutine k_add2(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector) :: p, q
         call nopadd2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
         p%time = p%time + q%time
         return
      end subroutine k_add2
      
      subroutine k_sub2(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector) :: p, q
         call nopsub2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
         p%time = p%time - q%time
      
         return
      end subroutine k_sub2
      
      subroutine k_sub3(p, q, r)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector) :: p, q, r
         call nopsub3(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz,
     $   q%pr, q%t, r%vx, r%vy, r%vz, r%pr, r%t)
         p%time = q%time - r%time
         return
      end subroutine k_sub3
      
      subroutine k_zero(p)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector) :: p
         call noprzero(p%vx, p%vy, p%vz, p%pr, p%t)
         p%time = 0.0d+00
         return
      end subroutine k_zero
      
      subroutine k_copy(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector) :: p, q
         call nopcopy(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
         p%time = q%time
         return
      end subroutine k_copy
      
      subroutine k_matmul(dq, Q, yvec, k)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         integer :: i, m, k
         type(krylov_vector) :: dq
         type(krylov_vector), dimension(k) :: Q
         real, dimension(k) :: yvec
      
         real, dimension(lv, k) :: qx, qy, qz
         real, dimension(lp, k) :: qp
         real, dimension(lv, k, ldimt) :: qt
      
         real, dimension(k) :: time_comp
      
         do i = 1, k
            qx(:, i) = Q(i)%vx(:)
            qy(:, i) = Q(i)%vy(:)
            qp(:, i) = Q(i)%pr(:)
            if (if3d) qz(:, i) = Q(i)%vz(:)
            if (ifto) qt(:, i, 1) = Q(i)%t(:, 1)
            if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) qt(:, i, m) = Q(i)%t(:, m)
            end do
            end if
            time_comp(i) = Q(i)%time
         end do
      
         call k_zero(dq)
      
         dq%vx(:) = matmul(qx(:, :), yvec(:))
         dq%vy(:) = matmul(qy(:, :), yvec(:))
         dq%pr(:) = matmul(qp(:, :), yvec(:))
         if (if3d) dq%vz(:) = matmul(qz(:, :), yvec(:))
         if (ifto) dq%t(:, 1) = matmul(qt(:, :, 1), yvec(:))
         if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) dq%t(:, m) = matmul(qt(:, :, m), yvec(:))
         end do
         end if
         dq%time = dot_product(time_comp(:), yvec(:))
      
         return
      end subroutine k_matmul
