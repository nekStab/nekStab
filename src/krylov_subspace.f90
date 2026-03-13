!-----------------------------------------------------------------------
! krylov_subspace.f90 — Krylov vector type and subspace operations
!
! Purpose:
!   Defines the krylov_vector derived type and provides optimized
!   vector operations for Krylov subspace methods in spectral
!   element discretizations. MPI-compatible for parallel execution.
!
! Public interface:
!   k_dot, k_norm, k_normalize, k_cmult, k_add2, k_add2s2,
!   k_axpby, k_sub2, k_sub3, k_zero, k_copy, k_matmul,
!   allocate_orbit, orbit_store, orbit_restore
!
! Dependencies:
!   nekstab_vectors, nekstab_nek_bridge, PARALLEL, SOLN, INPUT
!-----------------------------------------------------------------------
module krylov_subspace
   use nekstab_vectors
   use nekstab_nek_bridge
   implicit none

   private

   ! Mathematical constants
   real, public, parameter :: NEKSTAB_PI = 4.0d0*atan(1.0d0)

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

   public :: inner_product, norm, dnekclock, &
      k_dot, k_norm, k_normalize, k_cmult, &
      k_add2, k_add2s2, k_axpby, k_sub2, k_sub3, &
      k_zero, k_copy, k_matmul, &
      allocate_orbit, orbit_store, orbit_restore

   ! Interface for Nek5000 wall-clock timer
    interface
       function dnekclock() result(t)
          real :: t
       end function
    end interface

contains

!-----------------------------------------------------------------------
! Krylov vector operations
!
! Operations:
!   k_dot(alpha,p,q):      alpha = <p,q>      Inner product
!   k_norm(alpha,p):       alpha = ||p||       L2 norm
!   k_normalize(p,alpha):  p = p/||p||         Returns norm in alpha
!   k_cmult(p,c):          p = c*p             Scalar multiplication
!   k_add2(p,q):           p = p + q           Vector addition
!   k_add2s2(p,q,c):       p = p + c*q         Scaled addition (AXPY)
!   k_axpby(p,a,q,b):      p = a*p + b*q       Scaled combination
!   k_sub2(p,q):           p = p - q           Vector subtraction
!   k_sub3(p,q,r):         p = q - r           Three vector operation
!   k_zero(p):             p = 0               Zero vector
!   k_copy(p,q):           p = q               Copy vector
!   k_matmul(dq,Q,y,k):    dq = Q*y            Matrix-vector product
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! inner_product — Raw-array weighted L2 inner product (canonical primitive)
!
!   alpha = <p, q>_bm1s  (sponge-masked mass matrix)
!   Computes kinetic energy (vx,vy,[vz]) + thermal energy ([t,scalars])
!   Pressure dummy args exist for interface consistency but are unused.
!-----------------------------------------------------------------------
subroutine inner_product(alpha, &
     px, py, pz, pp, pt, qx, qy, qz, qp, qt)
   real, dimension(lv), intent(in) :: px, py, pz, qx, qy, qz
   real, dimension(lp), intent(in) :: pp, qp  ! not used
   real, dimension(lt, ldimt), intent(in) :: pt, qt
   real, intent(out) :: alpha
   real :: glsc3
   integer :: m

   nv = nx1*ny1*nz1*nelv
   nt = nx1*ny1*nz1*nelt

   alpha = glsc3(px, qx, bm1s, nv) + glsc3(py, qy, bm1s, nv)
   if (if3D) alpha = alpha + glsc3(pz, qz, bm1s, nv)
   if (ifto) alpha = alpha &
        + glsc3(pt(:,1), qt(:,1), bm1s, nt)
   if (ldimt > 1) then
      do m = 2, ldimt
         if (ifpsco(m-1)) alpha = alpha &
              + glsc3(pt(:,m), qt(:,m), bm1s, nt)
      end do
   end if

end subroutine inner_product

!-----------------------------------------------------------------------
! norm — Raw-array norm: alpha = sqrt(<q, q>_bm1s)
!-----------------------------------------------------------------------
subroutine norm(qx, qy, qz, qp, qt, alpha)
   real, intent(in), dimension(lv) :: qx, qy, qz
   real, intent(in), dimension(lp) :: qp
   real, intent(in), dimension(lt, ldimt) :: qt
   real, intent(out) :: alpha

   call inner_product(alpha, &
        qx, qy, qz, qp, qt, qx, qy, qz, qp, qt)
   alpha = sqrt(alpha)

end subroutine norm

!-----------------------------------------------------------------------
! k_dot — Weighted inner product of two Krylov vectors
!   Delegates to inner_product, then adds Newton PO time component.
!-----------------------------------------------------------------------
subroutine k_dot(alpha, p, q)
   type(krylov_vector), intent(in) :: p, q
   real, intent(out) :: alpha

   call inner_product(alpha, &
        p%vx, p%vy, p%vz, p%pr, p%t, &
        q%vx, q%vy, q%vz, q%pr, q%t)

   if (isNewtonPO) alpha = alpha + p%time*q%time

end subroutine k_dot

!-----------------------------------------------------------------------
! k_norm — L2 norm of a Krylov vector
!-----------------------------------------------------------------------
subroutine k_norm(alpha, p)
   type(krylov_vector), intent(in) :: p
   real, intent(out) :: alpha

   call k_dot(alpha, p, p)
   alpha = sqrt(alpha)

end subroutine k_norm

!-----------------------------------------------------------------------
! k_normalize — Normalize a Krylov vector to unit norm
!-----------------------------------------------------------------------
subroutine k_normalize(p, alpha)
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

!-----------------------------------------------------------------------
! k_cmult — Scalar multiplication p = c*p
!-----------------------------------------------------------------------
subroutine k_cmult(p, c)
   type(krylov_vector), intent(inout) :: p
   real, intent(in) :: c

   call nopcmult(p%vx, p%vy, p%vz, p%pr, p%t, c)
   p%time = p%time*c

end subroutine k_cmult

!-----------------------------------------------------------------------
! k_add2 — Vector addition p = p + q
!-----------------------------------------------------------------------
subroutine k_add2(p, q)
   type(krylov_vector), intent(inout) :: p
   type(krylov_vector), intent(in) :: q

   call nopadd2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
   p%time = p%time + q%time

end subroutine k_add2

!-----------------------------------------------------------------------
! k_add2s2 — Scaled addition p = p + c*q (AXPY)
!-----------------------------------------------------------------------
subroutine k_add2s2(p, q, c)
   type(krylov_vector), intent(inout) :: p
   type(krylov_vector), intent(in) :: q
   real, intent(in) :: c

   call nopadd2s2(p%vx, p%vy, p%vz, p%pr, p%t, &
        q%vx, q%vy, q%vz, q%pr, q%t, c)
   p%time = p%time + c*q%time

end subroutine k_add2s2

!-----------------------------------------------------------------------
! k_axpby — Scaled combination p = alpha*p + beta*q (AXPBY)
!-----------------------------------------------------------------------
subroutine k_axpby(p, alpha, q, beta)
   type(krylov_vector), intent(inout) :: p
   real, intent(in) :: alpha
   type(krylov_vector), intent(in) :: q
   real, intent(in) :: beta

   call nopaxpby(p%vx, p%vy, p%vz, p%pr, p%t, &
        alpha, q%vx, q%vy, q%vz, q%pr, q%t, beta)
   p%time = alpha*p%time + beta*q%time

end subroutine k_axpby

!-----------------------------------------------------------------------
! k_sub2 — Vector subtraction p = p - q
!-----------------------------------------------------------------------
subroutine k_sub2(p, q)
   type(krylov_vector), intent(inout) :: p
   type(krylov_vector), intent(in) :: q

   call nopsub2(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
   p%time = p%time - q%time

end subroutine k_sub2

!-----------------------------------------------------------------------
! k_sub3 — Three-vector subtraction p = q - r
!-----------------------------------------------------------------------
subroutine k_sub3(p, q, r)
   type(krylov_vector), intent(inout) :: p
   type(krylov_vector), intent(in) :: q, r

   call nopsub3(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, &
        q%pr, q%t, r%vx, r%vy, r%vz, r%pr, r%t)
   p%time = q%time - r%time

end subroutine k_sub3

!-----------------------------------------------------------------------
! k_zero — Zero all fields of a Krylov vector
!-----------------------------------------------------------------------
subroutine k_zero(p)
   type(krylov_vector), intent(inout) :: p

   call noprzero(p%vx, p%vy, p%vz, p%pr, p%t)
   p%time = 0.0d0

end subroutine k_zero

!-----------------------------------------------------------------------
! k_copy — Copy Krylov vector q into p
!-----------------------------------------------------------------------
subroutine k_copy(p, q)
   type(krylov_vector), intent(out) :: p
   type(krylov_vector), intent(in) :: q

   call nopcopy(p%vx, p%vy, p%vz, p%pr, p%t, q%vx, q%vy, q%vz, q%pr, q%t)
   p%time = q%time

end subroutine k_copy

!-----------------------------------------------------------------------
! k_matmul — Matrix-vector product dq = Q * yvec
!-----------------------------------------------------------------------
subroutine k_matmul(dq, Q, yvec, k)
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

!-----------------------------------------------------------------------
! allocate_orbit — Allocate orbit storage arrays on the heap
!-----------------------------------------------------------------------
subroutine allocate_orbit(nsteps_in)

   integer, intent(in) :: nsteps_in

   if (nid == 0) write (6, *) &
        'ALLOCATING ORBIT WITH NSTEPS:', nsteps_in

   allocate(uor(lv, nsteps_in), vor(lv, nsteps_in))
   if (if3d) then
      allocate(wor(lv, nsteps_in))
   else
      allocate(wor(1, 1))
   end if
   if (ifto .or. ldimt > 1) &
        allocate(tor(lt, nsteps_in, ldimt))

end subroutine allocate_orbit

!-----------------------------------------------------------------------
! orbit_store — Store current fields into orbit arrays
!-----------------------------------------------------------------------
subroutine orbit_store(istep_in)

   integer, intent(in) :: istep_in
   integer :: m

   nt = nx1*ny1*nz1*nelt

   call opcopy(uor(:,istep_in), vor(:,istep_in), &
        wor(:,istep_in), vx, vy, vz)
   if (ifto) call copy(tor(:,istep_in,1), &
        t(:,:,:,:,1), nt)
   if (ldimt > 1) then
      do m = 2, ldimt
         if (ifpsco(m-1)) call copy( &
              tor(:,istep_in,m), t(:,:,:,:,m), nt)
      end do
   end if

end subroutine orbit_store

!-----------------------------------------------------------------------
! orbit_restore — Restore fields from orbit arrays
!-----------------------------------------------------------------------
subroutine orbit_restore(istep_in)

   integer, intent(in) :: istep_in
   integer :: m

   nt = nx1*ny1*nz1*nelt

   call opcopy(vx, vy, vz, &
        uor(:,istep_in), vor(:,istep_in), &
        wor(:,istep_in))
   if (ifto) call copy(t(:,:,:,:,1), &
        tor(:,istep_in,1), nt)
   if (ldimt > 1) then
      do m = 2, ldimt
         if (ifpsco(m-1)) call copy( &
              t(:,:,:,:,m), tor(:,istep_in,m), nt)
      end do
   end if

end subroutine orbit_restore

end module krylov_subspace
