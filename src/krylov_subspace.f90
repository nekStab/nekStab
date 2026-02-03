      !-----------------------------------------------------------------------
      ! krylov_subspace: High-performance Krylov operations for spectral elements
      !
      ! Purpose:
      ! - Optimized vector operations for Krylov subspace methods
      ! - MPI-compatible for parallel execution on supercomputers
      ! - Minimizes memory usage and FLOPs in critical operations
      !
      ! Data Layout:
      ! - Velocity/Temperature: lv = lx1*ly1*lz1*lelv (first-order elements)
      ! - Pressure: lp = lx2*ly2*lz2*lelv (second-order elements)
      ! - Temperature scalars: lt = lx1*ly1*lz1*lelt
      !
      ! Memory Management:
      ! - Static allocation for vectors to avoid runtime overhead
      ! - Minimal temporary storage in vector operations
      ! - Optional orbit storage for UPO computations (deallocated when unused)
      !-----------------------------------------------------------------------
       module krylov_subspace
         implicit none
         include 'SIZE'
      
         private
      
         ! Spectral element dimensions (fixed at compile time)
         integer, public, parameter :: lv = lx1*ly1*lz1*lelv  ! Velocity/temp points
         integer, public, parameter :: lt = lx1*ly1*lz1*lelt  ! Temperature points
         integer, public, parameter :: lp = lx2*ly2*lz2*lelv  ! Pressure points
         
         ! Runtime dimensions (may be smaller than max)
         integer, save, public :: nv    ! Active velocity/temp points
         integer, save, public :: nt    ! Active temperature points  
         integer, save, public :: n2    ! Active pressure points (renamed from np for MPI)
      
         ! Core data structure for Krylov operations
         type, public :: krylov_vector
            real, dimension(lv) :: vx, vy, vz     ! Velocity components
            real, dimension(lp) :: pr             ! Pressure field
            real, dimension(lv, ldimt) :: t       ! Temperature/passive scalars
            real :: time                          ! Time value
         end type krylov_vector
      
         ! Global storage (minimized to essential data only)
         type(krylov_vector), save, public :: ic_nwt, fc_nwt  ! Newton conditions
         real, save, allocatable, dimension(:, :), public :: uor, vor, wor  ! Velocity orbits
         real, save, allocatable, dimension(:, :, :), public :: tor        ! Temperature orbits
      
      contains
      end module krylov_subspace

      !====================================================================
      !     Module: krylov_subspace
      !     
      !     Purpose: Implements Krylov subspace operations for numerical computations
      !              with support for velocity, pressure, and temperature fields.
      !              Designed for both serial and MPI parallel processing.
      !
      !     Operations:
      !     - k_dot(alpha,p,q):      alpha = <p,q>     ! Inner product
      !     - k_norm(alpha,p):       alpha = ||p||      ! L2 norm
      !     - k_normalize(p,alpha):  p = p/||p||        ! Returns norm in alpha
      !     - k_cmult(p,c):         p = c*p            ! Scalar multiplication
      !     - k_add2(p,q):          p = p + q          ! Vector addition
      !     - k_sub2(p,q):          p = p - q          ! Vector subtraction
      !     - k_sub3(p,q,r):        p = q - r          ! Three vector operation
      !     - k_zero(p):            p = 0              ! Zero vector
      !     - k_copy(p,q):          p = q              ! Copy vector
      !     - k_matmul(dq,Q,y,k):   dq = Q*y          ! Matrix-vector product
      !     
      !     Note: All operations preserve MPI compatibility
      !====================================================================
      
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
      
      !     --> Kinetic energy. ! Note : glsc3(a,b,mult,n)
         alpha = glsc3(p%vx, q%vx, bm1s, nv) + glsc3(p%vy, q%vy, bm1s, nv)
         if (if3d) alpha = alpha + glsc3(p%vz, q%vz, bm1s, nv)
      
      !     --> Potential energy.
         if (ifto) alpha = alpha + glsc3(p%t(:, 1), q%t(:, 1), bm1s, nt)
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) alpha = alpha + glsc3(p%t(:, m), q%t(:, m), bm1s, nt)
            end do
         end if
      
      !     --> Time component.
         if (uparam(1) == 2.1) then
            alpha = alpha + p%time*q%time
         end if
      
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
         inv_alpha = 1.0d0/alpha
      
      !     --> Normalize the vector.
         call k_cmult(p, inv_alpha)
      
         return
      end subroutine k_normalize
      
      subroutine k_cmult(p, c)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         real, intent(in) :: c
         call nopcmult(p%vx, p%vy, p%vz, p%pr, p%t, c)
         p%time = p%time*c
         return
      end subroutine k_cmult
      
      subroutine k_add2(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         type(krylov_vector), intent(in) :: q
         call nopadd2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
         p%time = p%time + q%time
         return
      end subroutine k_add2
      
      subroutine k_sub2(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         type(krylov_vector), intent(in) :: q
         call nopsub2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
         p%time = p%time - q%time
      
         return
      end subroutine k_sub2
      
      subroutine k_sub3(p, q, r)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         type(krylov_vector), intent(in) :: q, r
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
         type(krylov_vector), intent(inout) :: p
         call noprzero(p%vx, p%vy, p%vz, p%pr, p%t)
         p%time = 0.0d0
         return
      end subroutine k_zero
      
      subroutine k_copy(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(out) :: p
         type(krylov_vector), intent(in) :: q
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
         type(krylov_vector), intent(out) :: dq
         type(krylov_vector), dimension(k), intent(in) :: Q
         real, dimension(k), intent(in) :: yvec
      
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