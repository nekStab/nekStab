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
            real, dimension(lt, ldimt) :: t       ! Temperature/passive scalars
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
      !     - k_add2s2(p,q,c):      p = p + c*q        ! Scaled addition (AXPY)
      !     - k_axpby(p,a,q,b):    p = a*p + b*q      ! Scaled combination (AXPBY)
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
         if (isNewtonPO) then
            alpha = alpha + p%time*q%time
         end if

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
      end subroutine k_norm
      
      subroutine k_normalize(p, alpha)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         real, intent(out) :: alpha
         real :: inv_alpha

      !     --> Minimum norm threshold for warning about potential breakdown.
      !         A very small norm indicates the Arnoldi process may have found an
      !         invariant subspace (lucky breakdown), or severe loss of orthogonality.
      !         NOTE: We do NOT clamp alpha because:
      !         1. The caller uses alpha for H(k+1,k) in the Hessenberg matrix
      !         2. Modifying alpha would corrupt the Arnoldi factorization
      !         3. The lucky breakdown detection in update_hessenberg_matrix handles this
         real, parameter :: NORM_WARN_TOL = 1.0d-14

      !     --> Compute the user-defined norm.
         call k_norm(alpha, p)

      !     --> Warn if norm is dangerously small (but don't modify alpha).
         if (alpha < NORM_WARN_TOL) then
            if (nid == 0) then
               write(6,*) 'WARNING [k_normalize]: Near-zero norm =', alpha
               write(6,*) '         Possible invariant subspace or loss of orthogonality.'
            end if
         end if

      !     --> Guard against division by exactly zero (should be extremely rare).
      !         IMPORTANT: Do NOT modify alpha - the caller needs the true value for H(k+1,k).
      !         If alpha is zero, this indicates exact breakdown (invariant subspace).
      !         We zero the vector to avoid NaN but preserve alpha=0 for proper detection.
         if (alpha > 0.0d0) then
            inv_alpha = 1.0d0/alpha
            call k_cmult(p, inv_alpha)
         else
            if (nid == 0) then
               write(6,*) 'INFO [k_normalize]: Exact zero norm detected.'
               write(6,*) '     This indicates an invariant subspace (exact breakdown).'
               write(6,*) '     Vector zeroed, alpha preserved as 0 for H(k+1,k).'
            end if
            call k_zero(p)  ! Zero the vector to avoid undefined state
            ! alpha remains 0 - DO NOT modify it
         end if

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
      end subroutine k_add2

      subroutine k_add2s2(p, q, c)
!        p = p + c*q  (BLAS-style AXPY operation)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         type(krylov_vector), intent(in) :: q
         real, intent(in) :: c
         call nopadd2s2(p%vx, p%vy, p%vz, p%pr, p%t,
     $                  q%vx, q%vy, q%vz, q%pr, q%t, c)
         p%time = p%time + c*q%time
      end subroutine k_add2s2

      subroutine k_axpby(p, alpha, q, beta)
!        p = alpha*p + beta*q  (BLAS-style AXPBY operation)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         real, intent(in) :: alpha
         type(krylov_vector), intent(in) :: q
         real, intent(in) :: beta
         call nopaxpby(p%vx, p%vy, p%vz, p%pr, p%t,
     $                 alpha, q%vx, q%vy, q%vz, q%pr, q%t, beta)
         p%time = alpha*p%time + beta*q%time
      end subroutine k_axpby

      subroutine k_sub2(p, q)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         type(krylov_vector), intent(in) :: q
         call nopsub2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
         p%time = p%time - q%time
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
      end subroutine k_sub3
      
      subroutine k_zero(p)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), intent(inout) :: p
         call noprzero(p%vx, p%vy, p%vz, p%pr, p%t)
         p%time = 0.0d0
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
      end subroutine k_copy
      
      subroutine k_matmul(dq, Q, yvec, k)
!        dq = Q(:) * yvec(:)  (matrix-vector product, zero temporary storage)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         integer, intent(in) :: k
         type(krylov_vector), intent(out) :: dq
         type(krylov_vector), dimension(k), intent(in) :: Q
         real, dimension(k), intent(in) :: yvec
         integer :: i

         call k_zero(dq)
         do i = 1, k
            call k_add2s2(dq, Q(i), yvec(i))
         end do
      end subroutine k_matmul

      subroutine allocate_orbit(nsteps_in)
!
!     Allocate orbit storage arrays (uor, vor, wor, tor) on the heap.
!     Guards on if3d and ifto/ldimt ensure minimal memory usage.
!     Call once before the first orbit integration.
!
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'PARALLEL' ! nid
         include 'INPUT'    ! if3d, ifto, ifpsco

         integer, intent(in) :: nsteps_in

         if (nid == 0) write (6, *)
     $      'ALLOCATING ORBIT WITH NSTEPS:', nsteps_in

         allocate(uor(lv, nsteps_in), vor(lv, nsteps_in))
         if (if3d) then
            allocate(wor(lv, nsteps_in))
         else
            allocate(wor(1, 1))
         end if
         if (ifto .or. ldimt > 1)
     $      allocate(tor(lt, nsteps_in, ldimt))

      end subroutine allocate_orbit

      subroutine orbit_store(istep_in)
!
!     Store the current velocity (vx,vy,vz) and temperature (t) fields
!     into the orbit arrays at the given time-step index.
!
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'SOLN'    ! vx, vy, vz, t
         include 'INPUT'   ! if3d, ifto, ifpsco

         integer, intent(in) :: istep_in
         integer :: m

         nt = nx1*ny1*nz1*nelt

         call opcopy(uor(:,istep_in), vor(:,istep_in),
     $      wor(:,istep_in), vx, vy, vz)
         if (ifto) call copy(tor(:,istep_in,1),
     $      t(:,:,:,:,1), nt)
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m-1)) call copy(
     $            tor(:,istep_in,m), t(:,:,:,:,m), nt)
            end do
         end if

      end subroutine orbit_store

      subroutine orbit_restore(istep_in)
!
!     Restore the velocity (vx,vy,vz) and temperature (t) fields
!     from the orbit arrays at the given time-step index.
!
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'SOLN'    ! vx, vy, vz, t
         include 'INPUT'   ! if3d, ifto, ifpsco

         integer, intent(in) :: istep_in
         integer :: m

         nt = nx1*ny1*nz1*nelt

         call opcopy(vx, vy, vz,
     $      uor(:,istep_in), vor(:,istep_in),
     $      wor(:,istep_in))
         if (ifto) call copy(t(:,:,:,:,1),
     $      tor(:,istep_in,1), nt)
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m-1)) call copy(
     $            t(:,:,:,:,m), tor(:,istep_in,m), nt)
            end do
         end if

      end subroutine orbit_restore