!-----------------------------------------------------------------------
! krylov_subspace.f90 — Krylov vector type and subspace operations
!
! Purpose:
!   Defines the krylov_vector derived type and provides optimized
!   vector operations for Krylov subspace methods in spectral
!   element discretizations. MPI-compatible for parallel execution.
!
! Public interface:
!   k_dot      — weighted inner product of two Krylov vectors
!   k_norm     — L2 norm of a Krylov vector
!   k_normalize — normalize to unit norm (returns norm in alpha)
!   k_cmult    — scalar multiplication p = c*p
!   k_add2     — vector addition p = p + q
!   k_add2s2   — scaled addition p = p + c*q (AXPY)
!   k_axpby    — scaled combination p = alpha*p + beta*q
!   k_sub2     — vector subtraction p = p - q
!   k_sub3     — three-vector subtraction p = q - r
!   k_zero     — zero all fields of a Krylov vector
!   k_copy     — copy Krylov vector q into p
!   k_matmul   — matrix-vector product dq = Q * yvec
!   allocate_orbit, orbit_store, orbit_restore — Newton orbit storage
!
! Dependencies:
!   nekstab_vectors, nekstab_nek_bridge, PARALLEL, SOLN, INPUT
!-----------------------------------------------------------------------
module nekstab_krylov_subspace
   use nekstab_vectors, only: noprzero, nopcmult, nopaxpby, nopcopy, &
                              nopsub2, nopsub3, nopadd2, nopadd2s2
   use nekstab_nek_bridge, only: lx1, ly1, lz1, lx2, ly2, lz2, lelv, lelt, &
                                 ldimt, if3D, if3d, ifto, ifpsco, ifheat, nid, &
                                 nekStab_error, bm1s, vx, vy, vz, t, &
                                 thermal_norm_mode, thermal_norm_weight, &
                                 thermal_norm_min, thermal_norm_max, &
                                 thermal_buoyancy_coeff, isNewtonPO
   implicit none

   private

   ! Mathematical constants
   real, public, parameter :: NEKSTAB_PI = 4.0d0*atan(1.0d0)

   ! Spectral element dimensions (fixed at compile time)
   integer, public, parameter :: lv = lx1*ly1*lz1*lelv ! Velocity/temp points
   integer, public, parameter :: lt = lx1*ly1*lz1*lelt ! Temperature points
   integer, public, parameter :: lp = lx2*ly2*lz2*lelv ! Pressure points

   ! Runtime dimensions (may be smaller than max)
   integer, save, public :: nv ! Active velocity/temp points
   integer, save, public :: nt ! Active temperature points
   integer, save, public :: n2 ! Active pressure points (renamed from np for MPI)

   ! Core data structure for Krylov operations
   !  Temperature uses lt (=lx1*ly1*lz1*lelt), NOT lv (=lx1*ly1*lz1*lelv).
   !  In CHT cases lelt > lelv so lt > lv; using lv would silently truncate
   !  the solid-domain temperature and corrupt inner products.  All code that
   !  passes temperature arrays to nopcopy/dgemm must use lt as the leading
   !  dimension (not lv) to match this stride.
   type, public :: krylov_vector
      real, dimension(lv) :: vx, vy, vz ! Velocity components
      real, dimension(lp) :: pr ! Pressure field
      real, dimension(lt, ldimt) :: t ! Temperature/passive scalars
      real :: time ! Time value
   end type krylov_vector

   ! Global storage (minimized to essential data only)
   type(krylov_vector), save, public :: ic_nwt, fc_nwt ! Newton conditions
   real, save, allocatable, dimension(:, :), public :: uor, vor, wor ! Velocity orbits
   real, save, allocatable, dimension(:, :, :), public :: tor ! Temperature orbits

   public :: inner_product, norm, dnekclock, &
             k_dot, k_norm, k_normalize, k_cmult, &
             k_add2, k_add2s2, k_axpby, k_sub2, k_sub3, &
             k_zero, k_copy, k_matmul, &
             configure_thermal_norm_weight, &
             TN_MANUAL, TN_AUTO, TN_CLIP, &
             allocate_orbit, orbit_store, orbit_restore

   ! thermal_norm_mode selector values
   integer, parameter :: TN_MANUAL = 0, TN_AUTO = 1, TN_CLIP = 2

   ! Interface for Nek5000 wall-clock timer
   interface
      function dnekclock() result(t)
         real :: t
      end function dnekclock
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
      real, dimension(lp), intent(in) :: pp, qp ! not used
      real, dimension(lt, ldimt), intent(in) :: pt, qt
      real, intent(out) :: alpha
      real :: glsc3, temp_inner
      integer :: m

      ! nv, nt from globals (this module defines them)

      alpha = glsc3(px, qx, bm1s, nv) + glsc3(py, qy, bm1s, nv)
      if (if3D) alpha = alpha + glsc3(pz, qz, bm1s, nv)
      if (ifto) then
         ! Temperature is scalar slot 1.
         !
         ! Keep the implementation local to the inner product so that:
         ! 1. All Newton/Krylov norms inherit the same weighting automatically.
         ! 2. Vector update routines remain pure algebra on the stored fields.
         ! 3. Pressure stays excluded from the norm, which is appropriate for
         !    incompressible flow where pressure acts as a Lagrange multiplier.
         !
         ! Default thermal_norm_weight = 1.0 reproduces the historical behavior.
         ! The explicit multiplier exists so thermal continuation can rebalance
         ! velocity and temperature contributions without changing the state itself.
         temp_inner = glsc3(pt(:, 1), qt(:, 1), bm1s, nt)
         alpha = alpha + thermal_norm_weight*temp_inner
      end if
      if (ldimt > 1) then
         do m = 2, ldimt
            ! Additional passive scalars currently keep unit weight.
            ! This preserves existing behavior while leaving a clean extension
            ! path for future RANS-specific per-scalar weights.
            if (ifpsco(m - 1)) alpha = alpha &
                                       + glsc3(pt(:, m), qt(:, m), bm1s, nt)
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
! configure_thermal_norm_weight — Set the temperature weight in the norm
!
! Modes (thermal_norm_mode):
!   TN_MANUAL : leave thermal_norm_weight as set by the user
!   TN_AUTO   : weight = thermal_buoyancy_coeff * (E_u / E_T)
!   TN_CLIP   : TN_AUTO clipped to [thermal_norm_min, thermal_norm_max]
!
! Rationale (TN_AUTO):
!   In Boussinesq flow the momentum equation carries beta*T as forcing,
!   so scaling the temperature term by beta*(E_u/E_T) equalizes the
!   dynamical contribution of the two fields to the residual rather
!   than their raw L2 magnitude.
!
! glsc3 uses MPI_Allreduce, so every rank already holds identical
! results — no explicit bcast is needed on the outputs.
!-----------------------------------------------------------------------
   subroutine configure_thermal_norm_weight()
      real :: glsc3
      real :: kinetic_energy, thermal_energy

      if (.not. ifheat) return

      if (thermal_norm_mode == TN_MANUAL) then
         if (nid == 0) write (6, '(A,1PE12.4)') &
            'thermal_norm_weight (manual) = ', thermal_norm_weight
         return
      end if

      if (thermal_buoyancy_coeff <= 0.0d0) then
         call thermal_norm_fallback( &
            'thermal_norm_mode > 0 but thermal_buoyancy_coeff <= 0')
         return
      end if

      ! nv, nt from globals

      kinetic_energy = glsc3(vx, vx, bm1s, nv) + glsc3(vy, vy, bm1s, nv)
      if (if3D) kinetic_energy = kinetic_energy + glsc3(vz, vz, bm1s, nv)
      thermal_energy = glsc3(t(1, 1, 1, 1, 1), t(1, 1, 1, 1, 1), bm1s, nt)

      if (thermal_energy <= 1.0d-30) then
         call thermal_norm_fallback( &
            'thermal_energy is ~0 during automatic thermal norm setup')
         return
      end if

      thermal_norm_weight = thermal_buoyancy_coeff*(kinetic_energy/thermal_energy)
      if (thermal_norm_mode == TN_CLIP) &
         thermal_norm_weight = max(thermal_norm_min, &
                                   min(thermal_norm_max, thermal_norm_weight))

      if (nid == 0) then
         write (6, '(A,1PE12.4)') 'thermal_buoyancy_coeff = ', thermal_buoyancy_coeff
         write (6, '(A,1PE12.4)') 'thermal kinetic_energy = ', kinetic_energy
         write (6, '(A,1PE12.4)') 'thermal scalar_energy  = ', thermal_energy
         write (6, '(A,1PE12.4)') 'thermal_norm_weight    = ', thermal_norm_weight
      end if

   end subroutine configure_thermal_norm_weight

!-----------------------------------------------------------------------
! thermal_norm_fallback — Warn and reset to unit weight on bad input
!-----------------------------------------------------------------------
   subroutine thermal_norm_fallback(msg)
      character(len=*), intent(in) :: msg
      if (nid == 0) then
         write (6, '(A,A)') 'WARNING: ', trim(msg)
         write (6, '(A)') '         Falling back to thermal_norm_weight = 1.0.'
      end if
      thermal_norm_weight = 1.0d0

   end subroutine thermal_norm_fallback

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

      if (alpha /= alpha) then
         if (nid == 0) write (6, *) &
            'ERROR [k_normalize]: NaN norm; refusing to zero Krylov vector.'
         call exitti('k_normalize NaN norm$', 1)
      end if

      !     --> Warn if norm is dangerously small (but don't modify alpha).
      if (alpha < NORM_WARN_TOL) then
         if (nid == 0) then
            write (6, *) 'WARNING [k_normalize]: Near-zero norm =', alpha
            write (6, *) '         Possible invariant subspace or loss of orthogonality.'
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
            write (6, *) 'INFO [k_normalize]: Exact zero norm detected.'
            write (6, *) '     This indicates an invariant subspace (exact breakdown).'
            write (6, *) '     Vector zeroed, alpha preserved as 0 for H(k+1,k).'
         end if
         call k_zero(p) ! Zero the vector to avoid undefined state
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
      integer :: alloc_stat
      real :: request_gb

      if (nid == 0) write (6, *) &
         'ALLOCATING ORBIT WITH NSTEPS:', nsteps_in

      if (allocated(uor)) deallocate (uor, vor, wor)
      if (allocated(tor)) deallocate (tor)

      request_gb = 2.0d0*real(lv)*real(nsteps_in)*8.0d0/1.0d9
      allocate (uor(lv, nsteps_in), vor(lv, nsteps_in), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call nekStab_error('orbit allocation failed for uor/vor')
         call exitti('orbit allocation failed$', alloc_stat)
      end if
      if (if3d) then
         request_gb = real(lv)*real(nsteps_in)*8.0d0/1.0d9
         allocate (wor(lv, nsteps_in), stat=alloc_stat)
      else
         request_gb = 8.0d0/1.0d9
         allocate (wor(1, 1), stat=alloc_stat)
      end if
      if (alloc_stat /= 0) then
         call nekStab_error('orbit allocation failed for wor')
         call exitti('orbit allocation failed$', alloc_stat)
      end if
      if (ifto .or. ldimt > 1) then
         request_gb = real(lt)*real(nsteps_in)*real(ldimt)*8.0d0/1.0d9
         allocate (tor(lt, nsteps_in, ldimt), stat=alloc_stat)
         if (alloc_stat /= 0) then
            call nekStab_error('orbit allocation failed for tor')
            call exitti('orbit allocation failed$', alloc_stat)
         end if
      end if

   end subroutine allocate_orbit

!-----------------------------------------------------------------------
! orbit_store — Store current fields into orbit arrays
!-----------------------------------------------------------------------
   subroutine orbit_store(istep_in)

      integer, intent(in) :: istep_in
      integer :: m

      ! nt from globals

      call opcopy(uor(:, istep_in), vor(:, istep_in), &
                  wor(:, istep_in), vx, vy, vz)
      if (ifto) call copy(tor(:, istep_in, 1), &
                          t(:, :, :, :, 1), nt)
      if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) call copy( &
               tor(:, istep_in, m), t(:, :, :, :, m), nt)
         end do
      end if

   end subroutine orbit_store

!-----------------------------------------------------------------------
! orbit_restore — Restore fields from orbit arrays
!-----------------------------------------------------------------------
   subroutine orbit_restore(istep_in)

      integer, intent(in) :: istep_in
      integer :: m

      ! nt from globals

      call opcopy(vx, vy, vz, &
                  uor(:, istep_in), vor(:, istep_in), &
                  wor(:, istep_in))
      if (ifto) call copy(t(:, :, :, :, 1), &
                          tor(:, istep_in, 1), nt)
      if (ldimt > 1) then
         do m = 2, ldimt
            if (ifpsco(m - 1)) call copy( &
               t(:, :, :, :, m), tor(:, istep_in, m), nt)
         end do
      end if

   end subroutine orbit_restore

end module nekstab_krylov_subspace
