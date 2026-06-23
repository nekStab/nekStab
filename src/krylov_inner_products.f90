!-----------------------------------------------------------------------
! krylov_inner_products.f90 — Batch inner product operations
!
! Purpose:
!   Provides batched, BLAS-accelerated inner products and projections
!   (k_gram_matrix, k_project) plus complex Gram matrices for Krylov
!   solvers (Arnoldi, GMRES) and modal analysis (POD, DMD, SPOD).
!   Replaces O(k^2) per-pair reductions with dense dgemm + one gop.
!
! WHY THIS FILE EXISTS
! --------------------
! A krylov_vector stores one snapshot of the Nek5000 solution:
! velocity (vx,vy,vz), pressure (pr), and temperature/scalars (t).
! Many algorithms (Arnoldi, GMRES, POD, DMD, SPOD) need the
! weighted inner product  <a, b> = integral(a * b * bm1s)  between
! many pairs of krylov_vectors.
!
! The naive approach — calling glsc3 once per pair — triggers one
! MPI_Allreduce per call.  For k vectors that is O(k^2) reductions.
! This module replaces that with:
!   1. Gather the field data into dense contiguous arrays  (copy)
!   2. Weight the left array by the mass matrix             (col2)
!   3. Compute ALL dot products at once via BLAS dgemm/dgemv
!   4. Do ONE MPI_Allreduce for the entire result matrix    (gop)
!
! HOW IT WORKS (per field component)
! -----------------------------------
! For velocity-x, given m left vectors P and n right vectors Q:
!
!   Aw(:,i) = P(i)%vx * bm1s    — mass-weighted left vectors
!   B(:,j)  = Q(j)%vx            — unweighted right vectors
!   G      += Aw^T * B            — local dot products via dgemm
!
! This is repeated for vy, vz (if 3D), and temperature/scalars.
! Each field ACCUMULATES into G (beta=1 after the first field).
! After all fields, a single gop sums G across MPI ranks.
!
! BLAS STRIDE CONVENTION
! ----------------------
! Arrays are allocated with leading dimension lv (compile-time max),
! but only nv rows (runtime actual grid points) are filled and used.
! dgemm is called with M=nv but LDA=lv so the stride matches the
! Fortran column layout.  Rows nv+1:lv are never read.
! Temperature uses lt/nt instead of lv/nv (different grid for CHT).
!
! WORKSPACE STRATEGY
! ------------------
! The gather arrays (Aw, B, wk_gop, etc.) are allocated once as
! module-level persistent buffers and reused across calls.  This
! avoids allocator overhead on every Arnoldi/GMRES step.  Each
! buffer tracks its capacity; reallocation happens only if a
! caller requests more columns than any previous call.
!
! Public interface:
!   k_gram_matrix  — G(i,j) = <P(i), Q(j)>_bm1s     (matrix)
!   k_project      — h(i)   = <Q(i), f>_bm1s         (vector)
!   k_gram_complex — CSD(i,j) for complex krylov_vectors (SPOD)
!
! Dependencies:
!   nekstab_krylov_subspace, nekstab_nek_bridge
!-----------------------------------------------------------------------

module nekstab_krylov_inner_products
   use nekstab_krylov_subspace, only: lv, lt, nv, nt, krylov_vector
   use nekstab_nek_bridge, only: bm1s, if3d, ifto, ldimt, ifpsco, isNewtonPO, &
                                 nekStab_dp, ifKEnorm
   implicit none
   private

   public :: k_gram_matrix, k_project, k_gram_complex

   ! ── Persistent workspace (module-level, allocated once) ──────────
   !
   ! These arrays replace local allocatables that were allocated and
   ! deallocated on every call to k_gram_matrix / k_project /
   ! k_gram_complex.  Since those routines sit on the inner Krylov
   ! path (Arnoldi, GMRES, POD, DMD, SPOD), the per-call allocator
   ! churn was a measurable overhead once the numerical work itself
   ! was already batched into dgemm/dgemv + single gop.
   !
   ! Each buffer is sized to the high-water mark seen so far.  The
   ! ensure_* helpers below reallocate only when the requested
   ! dimension exceeds the current capacity, so repeated calls with
   ! the same (or smaller) dimensions cost nothing.
   !
   ! Velocity workspace uses leading dimension lv (compile-time);
   ! temperature workspace uses lt.  They are separate because
   ! lv != lt when lelt > lelv (CHT cases).
   !
   ! The _cap variables track the current second-dimension capacity
   ! so we know when reallocation is needed.  We cannot use size()
   ! instead because Fortran does NOT guarantee short-circuit
   ! evaluation of .or. — writing (.not. allocated(X) .or. size(X,2) < n)
   ! can evaluate size() on an unallocated array (undefined behavior).

   ! k_gram_matrix workspace
   real, allocatable, save :: gram_vel_aw_s(:, :), gram_vel_b_s(:, :)
   real, allocatable, save :: gram_temp_aw_s(:, :), gram_temp_b_s(:, :)

   ! k_gram_complex workspace (4 real Gram sub-matrices + pack buffer)
   real, allocatable, save :: gram_complex_rr_s(:, :), gram_complex_ii_s(:, :)
   real, allocatable, save :: gram_complex_ri_s(:, :), gram_complex_ir_s(:, :)
   real, allocatable, save :: gram_complex_all_s(:)
   real, allocatable, save :: gram_complex_vel_aw_s(:, :), gram_complex_vel_b_s(:, :)
   real, allocatable, save :: gram_complex_temp_aw_s(:, :), gram_complex_temp_b_s(:, :)

   ! k_project workspace
   real, allocatable, save :: proj_vel_aw_s(:, :), proj_vel_fwrk_s(:)
   real, allocatable, save :: proj_temp_aw_s(:, :), proj_temp_fwrk_s(:)

   ! Shared gop scratch (used by all 3 public routines — never concurrent)
   real, allocatable, save :: wk_gop_s(:)

   ! Capacity tracking (second dimension of each 2D buffer)
   integer, save :: gram_vel_m_cap = 0, gram_vel_n_cap = 0
   integer, save :: gram_temp_m_cap = 0, gram_temp_n_cap = 0
   integer, save :: gram_complex_nblk_cap = 0
   integer, save :: gram_complex_vel_nblk_cap = 0
   integer, save :: gram_complex_temp_nblk_cap = 0
   integer, save :: gram_complex_all_cap = 0
   integer, save :: proj_vel_k_cap = 0, proj_temp_k_cap = 0
   integer, save :: wk_gop_cap = 0

contains

! ── Workspace helpers ───────────────────────────────────────────────
! Each helper guarantees that the corresponding module-level buffer
! has at least the requested capacity.  If the buffer is already
! large enough it is reused as-is (no allocation, no zero-fill).
! Reallocation only happens when a caller needs MORE columns than
! any previous call — which is rare after the first Krylov cycle.

   subroutine ensure_gram_velocity_workspace(m, n)
      integer, intent(in) :: m, n

      if (.not. allocated(gram_vel_aw_s) .or. gram_vel_m_cap < m) then
         if (allocated(gram_vel_aw_s)) deallocate (gram_vel_aw_s)
         allocate (gram_vel_aw_s(lv, m))
         gram_vel_m_cap = m
      end if

      if (.not. allocated(gram_vel_b_s) .or. gram_vel_n_cap < n) then
         if (allocated(gram_vel_b_s)) deallocate (gram_vel_b_s)
         allocate (gram_vel_b_s(lv, n))
         gram_vel_n_cap = n
      end if

   end subroutine ensure_gram_velocity_workspace

   subroutine ensure_gram_temperature_workspace(m, n)
      integer, intent(in) :: m, n

      if (.not. allocated(gram_temp_aw_s) .or. gram_temp_m_cap < m) then
         if (allocated(gram_temp_aw_s)) deallocate (gram_temp_aw_s)
         allocate (gram_temp_aw_s(lt, m))
         gram_temp_m_cap = m
      end if

      if (.not. allocated(gram_temp_b_s) .or. gram_temp_n_cap < n) then
         if (allocated(gram_temp_b_s)) deallocate (gram_temp_b_s)
         allocate (gram_temp_b_s(lt, n))
         gram_temp_n_cap = n
      end if

   end subroutine ensure_gram_temperature_workspace

   subroutine ensure_gop_workspace(n)
      integer, intent(in) :: n

      if (.not. allocated(wk_gop_s) .or. wk_gop_cap < n) then
         if (allocated(wk_gop_s)) deallocate (wk_gop_s)
         allocate (wk_gop_s(n))
         wk_gop_cap = n
      end if

   end subroutine ensure_gop_workspace

   subroutine ensure_complex_gram_matrix_workspace(nblk)
      integer, intent(in) :: nblk
      integer :: nn

!  Sub-matrices: exact-match (not grow-only) because both dimensions
!  equal nblk — a larger allocation would give the wrong LDA when
!  passed to gram4_field_vel/temp.
      if (.not. allocated(gram_complex_rr_s) .or. gram_complex_nblk_cap /= nblk) then
         if (allocated(gram_complex_rr_s)) then
            deallocate (gram_complex_rr_s, gram_complex_ii_s, &
                        gram_complex_ri_s, gram_complex_ir_s)
         end if
         allocate (gram_complex_rr_s(nblk, nblk), gram_complex_ii_s(nblk, nblk))
         allocate (gram_complex_ri_s(nblk, nblk), gram_complex_ir_s(nblk, nblk))
         gram_complex_nblk_cap = nblk
      end if

!  Pack buffer (1D, no LDA issue — grow-only is safe)
      nn = 4*nblk*nblk
      if (.not. allocated(gram_complex_all_s) .or. gram_complex_all_cap < nn) then
         if (allocated(gram_complex_all_s)) deallocate (gram_complex_all_s)
         allocate (gram_complex_all_s(nn))
         gram_complex_all_cap = nn
      end if

   end subroutine ensure_complex_gram_matrix_workspace

   subroutine ensure_complex_velocity_workspace(nblk)
      integer, intent(in) :: nblk

      if (.not. allocated(gram_complex_vel_aw_s) .or. &
          gram_complex_vel_nblk_cap < nblk) then
         if (allocated(gram_complex_vel_aw_s)) then
            deallocate (gram_complex_vel_aw_s, gram_complex_vel_b_s)
         end if
         allocate (gram_complex_vel_aw_s(lv, nblk), gram_complex_vel_b_s(lv, nblk))
         gram_complex_vel_nblk_cap = nblk
      end if

   end subroutine ensure_complex_velocity_workspace

   subroutine ensure_complex_temperature_workspace(nblk)
      integer, intent(in) :: nblk

      if (.not. allocated(gram_complex_temp_aw_s) .or. &
          gram_complex_temp_nblk_cap < nblk) then
         if (allocated(gram_complex_temp_aw_s)) then
            deallocate (gram_complex_temp_aw_s, gram_complex_temp_b_s)
         end if
         allocate (gram_complex_temp_aw_s(lt, nblk), gram_complex_temp_b_s(lt, nblk))
         gram_complex_temp_nblk_cap = nblk
      end if

   end subroutine ensure_complex_temperature_workspace

   subroutine ensure_project_velocity_workspace(k)
      integer, intent(in) :: k

      if (.not. allocated(proj_vel_aw_s) .or. proj_vel_k_cap < k) then
         if (allocated(proj_vel_aw_s)) deallocate (proj_vel_aw_s)
         allocate (proj_vel_aw_s(lv, k))
         proj_vel_k_cap = k
      end if

      if (.not. allocated(proj_vel_fwrk_s)) then
         allocate (proj_vel_fwrk_s(lv))
      end if

   end subroutine ensure_project_velocity_workspace

   subroutine ensure_project_temperature_workspace(k)
      integer, intent(in) :: k

      if (.not. allocated(proj_temp_aw_s) .or. proj_temp_k_cap < k) then
         if (allocated(proj_temp_aw_s)) deallocate (proj_temp_aw_s)
         allocate (proj_temp_aw_s(lt, k))
         proj_temp_k_cap = k
      end if

      if (.not. allocated(proj_temp_fwrk_s)) then
         allocate (proj_temp_fwrk_s(lt))
      end if

   end subroutine ensure_project_temperature_workspace

!-----------------------------------------------------------------------
! k_gram_matrix — Batch Gram matrix: G(i,j) = <P(i), Q(j)>_bm1s
!
! Computes an m-by-n matrix of mass-weighted inner products between
! two sets of krylov_vectors.  Used by POD (correlation matrix),
! DMD (snapshot Gram matrix), and Arnoldi (orthogonalization).
!
! Algorithm:
!   For each field component (vx, vy, [vz], [temperature], [scalars]):
!     1. Gather left vectors into Aw(:,i) = P(i)%field * bm1s
!     2. Gather right vectors into B(:,j) = Q(j)%field
!     3. dgemm: G += Aw^T * B  (beta=0 on first field, beta=1 after)
!   Single gop over all m*n entries at the end.
!
! The beta=0/1 trick avoids zeroing G explicitly: the first dgemm
! writes G = Aw_vx^T * B_vx, and subsequent calls accumulate into
! it.  After all fields, each G(i,j) holds the LOCAL sum of
! dot products across all components.  The final gop sums across
! MPI ranks.
!
! ldg: leading dimension of G (may exceed m if caller has a larger
!      array, e.g., POD passes ldg=nsnap for an nsnap-by-nsnap matrix).
!-----------------------------------------------------------------------
   subroutine k_gram_matrix(G, P, Q, m, n, ldg)
      integer, intent(in) :: m, n, ldg
      real, intent(out) :: G(ldg, n)
      type(krylov_vector), intent(in) :: P(m), Q(n)
      integer :: i, j, mm

!  Runtime grid sizes (nv/nt) come from the globals in nekstab_krylov_subspace
!  (set once in nekStab_init). See rationale in nekstab_krylov_subspace.f90.

!  Ensure persistent workspace is large enough for this call.
!  gop needs a scratch buffer of the same length as G's data (m*n).
      call ensure_gram_velocity_workspace(m, n)
      call ensure_gop_workspace(m*n)

!  ── vx: first field, beta=0 overwrites G (no explicit zero needed) ──
      do i = 1, m
         call copy(gram_vel_aw_s(1, i), P(i)%vx, nv)
         call col2(gram_vel_aw_s(1, i), bm1s, nv)
      end do
      do j = 1, n
         call copy(gram_vel_b_s(1, j), Q(j)%vx, nv)
      end do
      call dgemm('T', 'N', m, n, nv, &
                 1.0d0, gram_vel_aw_s, lv, gram_vel_b_s, lv, 0.0d0, G, ldg)

!  ── vy: beta=1 accumulates into G (G += Aw_vy^T * B_vy) ──
      do i = 1, m
         call copy(gram_vel_aw_s(1, i), P(i)%vy, nv)
         call col2(gram_vel_aw_s(1, i), bm1s, nv)
      end do
      do j = 1, n
         call copy(gram_vel_b_s(1, j), Q(j)%vy, nv)
      end do
      call dgemm('T', 'N', m, n, nv, &
                 1.0d0, gram_vel_aw_s, lv, gram_vel_b_s, lv, 1.0d0, G, ldg)

      ! ── vz (3D only, accumulate) ──
      if (if3d) then
         do i = 1, m
            call copy(gram_vel_aw_s(1, i), P(i)%vz, nv)
            call col2(gram_vel_aw_s(1, i), bm1s, nv)
         end do
         do j = 1, n
            call copy(gram_vel_b_s(1, j), Q(j)%vz, nv)
         end do
         call dgemm('T', 'N', m, n, nv, &
                    1.0d0, gram_vel_aw_s, lv, gram_vel_b_s, lv, 1.0d0, G, ldg)
      end if

!  ── Temperature / passive scalars ──
!  Temperature uses separate workspace with leading dimension lt
!  (not lv) because lelt may differ from lelv in CHT setups.
!  ifto: temperature is active.  ldimt>1: additional passive scalars
!  are present; ifpsco(mm-1) checks if scalar mm is solved.
      if ((ifto .or. (ldimt > 1)) .and. .not. ifKEnorm) then  ! School B
         call ensure_gram_temperature_workspace(m, n)

         if (ifto) then
            do i = 1, m
               call copy(gram_temp_aw_s(1, i), P(i)%t(1, 1), nt)
               call col2(gram_temp_aw_s(1, i), bm1s, nt)
            end do
            do j = 1, n
               call copy(gram_temp_b_s(1, j), Q(j)%t(1, 1), nt)
            end do
!        LDA=lt (not lv) — temperature grid may have different size
            call dgemm('T', 'N', m, n, nt, &
                       1.0d0, gram_temp_aw_s, lt, gram_temp_b_s, lt, 1.0d0, G, ldg)
         end if

!     Additional passive scalars (species, turbulence quantities, etc.)
         if (ldimt > 1) then
            do mm = 2, ldimt
               if (ifpsco(mm - 1)) then
                  do i = 1, m
                     call copy(gram_temp_aw_s(1, i), P(i)%t(1, mm), nt)
                     call col2(gram_temp_aw_s(1, i), bm1s, nt)
                  end do
                  do j = 1, n
                     call copy(gram_temp_b_s(1, j), Q(j)%t(1, mm), nt)
                  end do
                  call dgemm('T', 'N', m, n, nt, &
                             1.0d0, gram_temp_aw_s, lt, gram_temp_b_s, lt, 1.0d0, G, ldg)
               end if
            end do
         end if
      end if

!  ── Global reduction: sum local dot products across MPI ranks ──
!  G currently holds rank-local partial sums.  gop performs a single
!  MPI_Allreduce over all m*n entries at once — this is the payoff
!  for batching (one allreduce replaces O(m*n) individual ones).
      call gop(G, wk_gop_s, '+  ', m*n)

!  ── Newton periodic orbit: add time-component inner product ──
!  Newton-PO augments the state vector with a scalar time variable.
!  Since 'time' is a single scalar (identical on all ranks), it does
!  not need MPI reduction — so it is added AFTER the gop.
      if (isNewtonPO) then
         do j = 1, n
            do i = 1, m
               G(i, j) = G(i, j) + P(i)%time*Q(j)%time
            end do
         end do
      end if

   end subroutine k_gram_matrix

!-----------------------------------------------------------------------
! k_project — Batch projection: h(i) = <Q(i), f>_bm1s for i=1..k
!
! Computes k inner products of a single vector f against a set of k
! basis vectors Q.  This is the inner loop of Arnoldi/GMRES
! orthogonalization (h = V^T * f in Gram-Schmidt) and a building
! block for DMD shifted Gram rows.
!
! Identical pattern to k_gram_matrix but uses dgemv (matrix-vector)
! instead of dgemm because f is a single vector, not a matrix.
! The result h is a vector of length k (not a matrix).
!
! The left vectors Q are mass-weighted (gathered into Aw * bm1s);
! f is gathered into a plain work vector fwrk.  Each field
! contributes via  h += Aw^T * fwrk  (beta=0 for vx, 1 after).
!-----------------------------------------------------------------------
   subroutine k_project(h, Q, f, k)
      integer, intent(in) :: k
      real, intent(out) :: h(k)
      type(krylov_vector), dimension(k), intent(in) :: Q
      type(krylov_vector), intent(in) :: f
      integer :: i, mm

      ! nv/nt from globals (set at init)

      call ensure_project_velocity_workspace(k)
      call ensure_gop_workspace(k)

!  ── vx: beta=0 overwrites h ──
      do i = 1, k
         call copy(proj_vel_aw_s(1, i), Q(i)%vx, nv)
         call col2(proj_vel_aw_s(1, i), bm1s, nv)
      end do
      call copy(proj_vel_fwrk_s, f%vx, nv)
      call dgemv('T', nv, k, 1.0d0, &
                 proj_vel_aw_s, lv, proj_vel_fwrk_s, 1, 0.0d0, h, 1)

      ! ── vy (accumulate) ──
      do i = 1, k
         call copy(proj_vel_aw_s(1, i), Q(i)%vy, nv)
         call col2(proj_vel_aw_s(1, i), bm1s, nv)
      end do
      call copy(proj_vel_fwrk_s, f%vy, nv)
      call dgemv('T', nv, k, 1.0d0, &
                 proj_vel_aw_s, lv, proj_vel_fwrk_s, 1, 1.0d0, h, 1)

      ! ── vz (3D only) ──
      if (if3d) then
         do i = 1, k
            call copy(proj_vel_aw_s(1, i), Q(i)%vz, nv)
            call col2(proj_vel_aw_s(1, i), bm1s, nv)
         end do
         call copy(proj_vel_fwrk_s, f%vz, nv)
         call dgemv('T', nv, k, 1.0d0, &
                    proj_vel_aw_s, lv, proj_vel_fwrk_s, 1, 1.0d0, h, 1)
      end if

!  ── Temperature / passive scalars (LDA=lt, M=nt) ──
      if ((ifto .or. (ldimt > 1)) .and. .not. ifKEnorm) then  ! School B
         call ensure_project_temperature_workspace(k)

         if (ifto) then
            do i = 1, k
               call copy(proj_temp_aw_s(1, i), Q(i)%t(1, 1), nt)
               call col2(proj_temp_aw_s(1, i), bm1s, nt)
            end do
            call copy(proj_temp_fwrk_s, f%t(1, 1), nt)
            call dgemv('T', nt, k, 1.0d0, &
                       proj_temp_aw_s, lt, proj_temp_fwrk_s, 1, 1.0d0, h, 1)
         end if

         if (ldimt > 1) then
            do mm = 2, ldimt
               if (ifpsco(mm - 1)) then
                  do i = 1, k
                     call copy(proj_temp_aw_s(1, i), Q(i)%t(1, mm), nt)
                     call col2(proj_temp_aw_s(1, i), bm1s, nt)
                  end do
                  call copy(proj_temp_fwrk_s, f%t(1, mm), nt)
                  call dgemv('T', nt, k, 1.0d0, &
                             proj_temp_aw_s, lt, proj_temp_fwrk_s, 1, 1.0d0, h, 1)
               end if
            end do
         end if
      end if

!  ── Single global reduction (one allreduce for all k entries) ──
      call gop(h, wk_gop_s, '+  ', k)

!  ── Newton PO time component (rank-local, added after gop) ──
      if (isNewtonPO) then
         do i = 1, k
            h(i) = h(i) + Q(i)%time*f%time
         end do
      end if

   end subroutine k_project

!-----------------------------------------------------------------------
! k_gram_complex — Hermitian CSD matrix for complex krylov_vectors
!
! Used by SPOD: each frequency bin has complex-valued Fourier
! coefficients stored as separate real and imaginary krylov_vectors.
! The cross-spectral density matrix requires the complex inner
! product:
!
!   CSD(i,j) = <a_re + i*a_im,  b_re + i*b_im>
!
! Expanding this product gives four real-valued Gram matrices:
!
!   Re(CSD) = <a_re, b_re> + <a_im, b_im>  =  G_rr + G_ii
!   Im(CSD) = <a_re, b_im> - <a_im, b_re>  =  G_ri - G_ir
!
! We compute all four real Grams (G_rr, G_ii, G_ri, G_ir) using the
! same batch-dgemm pattern as k_gram_matrix, pack them into one
! contiguous array, do a SINGLE gop over all 4*nblk^2 entries, then
! unpack and assemble into complex CSD.  This gives one MPI_Allreduce
! instead of 4*nblk^2 individual ones.
!-----------------------------------------------------------------------
   subroutine k_gram_complex(CSD, blk_re, blk_im, nblk, ldcsd)
      integer, intent(in) :: nblk, ldcsd
      complex(nekStab_dp), intent(out) :: CSD(ldcsd, nblk)
      type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)

      integer :: i, j, mm, off

      ! nv/nt from globals

      call ensure_complex_gram_matrix_workspace(nblk)
      call ensure_gop_workspace(4*nblk*nblk)
      call ensure_complex_velocity_workspace(nblk)

!  ── Velocity fields: each call computes all 4 sub-matrices ──
!  vx: beta_flag=0 overwrites the sub-matrices
      call gram4_field_vel(gram_complex_rr_s, gram_complex_ii_s, &
                           gram_complex_ri_s, gram_complex_ir_s, blk_re, blk_im, nblk, &
                           gram_complex_vel_aw_s, gram_complex_vel_b_s, 'x', 0)
!  vy: beta_flag=1 accumulates
      call gram4_field_vel(gram_complex_rr_s, gram_complex_ii_s, &
                           gram_complex_ri_s, gram_complex_ir_s, blk_re, blk_im, nblk, &
                           gram_complex_vel_aw_s, gram_complex_vel_b_s, 'y', 1)
!  vz (3D only)
      if (if3d) then
         call gram4_field_vel(gram_complex_rr_s, gram_complex_ii_s, &
                              gram_complex_ri_s, gram_complex_ir_s, blk_re, blk_im, nblk, &
                              gram_complex_vel_aw_s, gram_complex_vel_b_s, 'z', 1)
      end if

!  ── Temperature / passive scalars ──
      if ((ifto .or. (ldimt > 1)) .and. .not. ifKEnorm) then  ! School B
         call ensure_complex_temperature_workspace(nblk)

         if (ifto) then
            call gram4_field_temp(gram_complex_rr_s, gram_complex_ii_s, &
                                  gram_complex_ri_s, gram_complex_ir_s, blk_re, blk_im, nblk, &
                                  gram_complex_temp_aw_s, gram_complex_temp_b_s, 1)
         end if
         if (ldimt > 1) then
            do mm = 2, ldimt
               if (ifpsco(mm - 1)) then
                  call gram4_field_temp(gram_complex_rr_s, gram_complex_ii_s, &
                                        gram_complex_ri_s, gram_complex_ir_s, blk_re, blk_im, &
                                        nblk, gram_complex_temp_aw_s, gram_complex_temp_b_s, mm)
               end if
            end do
         end if
      end if

!  ── Pack 4 sub-matrices into one contiguous 1D array ──
!  gop requires a contiguous buffer.  We lay out [G_rr | G_ii | G_ri | G_ir]
!  end-to-end so a single gop call reduces all 4*nblk^2 entries.
      off = 0
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_all_s(off) = gram_complex_rr_s(i, j)
         end do
      end do
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_all_s(off) = gram_complex_ii_s(i, j)
         end do
      end do
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_all_s(off) = gram_complex_ri_s(i, j)
         end do
      end do
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_all_s(off) = gram_complex_ir_s(i, j)
         end do
      end do

!  ── Single global reduction over all 4 matrices at once ──
      call gop(gram_complex_all_s, wk_gop_s, '+  ', 4*nblk*nblk)

!  ── Unpack back into named sub-matrices ──
      off = 0
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_rr_s(i, j) = gram_complex_all_s(off)
         end do
      end do
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_ii_s(i, j) = gram_complex_all_s(off)
         end do
      end do
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_ri_s(i, j) = gram_complex_all_s(off)
         end do
      end do
      do j = 1, nblk
         do i = 1, nblk
            off = off + 1
            gram_complex_ir_s(i, j) = gram_complex_all_s(off)
         end do
      end do

!  ── Assemble complex CSD from the 4 real sub-matrices ──
!  Re(CSD) = G_rr + G_ii,  Im(CSD) = G_ri - G_ir
      do j = 1, nblk
         do i = 1, nblk
            CSD(i, j) = dcmplx(gram_complex_rr_s(i, j) + &
                               gram_complex_ii_s(i, j), &
                               gram_complex_ri_s(i, j) - gram_complex_ir_s(i, j))
         end do
      end do

   end subroutine k_gram_complex

!-----------------------------------------------------------------------
! gram4_field_vel — Accumulate 4 Gram sub-matrices for one velocity
!                   field component (vx, vy, or vz)
!
! For each velocity field we need four combinations:
!   G_rr(i,j) += <re(i), re(j)>     (real-real)
!   G_ii(i,j) += <im(i), im(j)>     (imag-imag)
!   G_ri(i,j) += <re(i), im(j)>     (real-imag, cross term)
!   G_ir(i,j) += <im(i), re(j)>     (imag-real, cross term)
!
! Each requires gathering nblk vectors into Aw (mass-weighted left)
! and B (plain right), then one dgemm.  That is 4 dgemm calls per
! field, totalling 8 gather+dgemm for vx+vy (2D) or 12 for 3D.
!
! The two workspace arrays Aw and B are shared across all four
! sub-matrices — each gather overwrites the previous contents.
! Where possible we reuse data that is still in the buffers:
!   - G_ri reuses B (still holds im from G_ii), only regathers Aw
!   - G_ir regathers both (im for Aw, re for B)
!
! beta_flag: 0 = first field (overwrite sub-matrices),
!            1 = subsequent fields (accumulate into sub-matrices).
!-----------------------------------------------------------------------
   subroutine gram4_field_vel(G_rr, G_ii, G_ri, G_ir, &
                              blk_re, blk_im, nblk, Aw, B, field, beta_flag)
      integer, intent(in) :: nblk, beta_flag
      real, intent(inout) :: G_rr(nblk, nblk), G_ii(nblk, nblk)
      real, intent(inout) :: G_ri(nblk, nblk), G_ir(nblk, nblk)
      type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)
!  Caller-provided workspace arrays (lv x nblk each).
!  Only rows 1:nv are written by copy+col2; dgemm uses M=nv, LDA=lv
!  so the stride matches the Fortran column layout.  Rows nv+1:lv
!  are never read — they are padding from compile-time allocation.
      real, intent(out) :: Aw(lv, nblk), B(lv, nblk)
      character(len=1), intent(in) :: field
      real :: dbeta
      integer :: i

      dbeta = dble(beta_flag)

!  ── G_rr: <re, re> — gather real parts into both Aw and B ──
      do i = 1, nblk
         if (field == 'x') then
            call copy(Aw(1, i), blk_re(i)%vx, nv)
            call col2(Aw(1, i), bm1s, nv)
         elseif (field == 'y') then
            call copy(Aw(1, i), blk_re(i)%vy, nv)
            call col2(Aw(1, i), bm1s, nv)
         else
            call copy(Aw(1, i), blk_re(i)%vz, nv)
            call col2(Aw(1, i), bm1s, nv)
         end if
      end do
      do i = 1, nblk
         if (field == 'x') then
            call copy(B(1, i), blk_re(i)%vx, nv)
         elseif (field == 'y') then
            call copy(B(1, i), blk_re(i)%vy, nv)
         else
            call copy(B(1, i), blk_re(i)%vz, nv)
         end if
      end do
      call dgemm('T', 'N', nblk, nblk, nv, &
                 1.0d0, Aw, lv, B, lv, dbeta, G_rr, nblk)

!  ── G_ii: <im, im> — gather imaginary parts into Aw and B ──
!  This overwrites the real-part data that was in Aw and B.
      do i = 1, nblk
         if (field == 'x') then
            call copy(Aw(1, i), blk_im(i)%vx, nv)
            call col2(Aw(1, i), bm1s, nv)
         elseif (field == 'y') then
            call copy(Aw(1, i), blk_im(i)%vy, nv)
            call col2(Aw(1, i), bm1s, nv)
         else
            call copy(Aw(1, i), blk_im(i)%vz, nv)
            call col2(Aw(1, i), bm1s, nv)
         end if
      end do
      do i = 1, nblk
         if (field == 'x') then
            call copy(B(1, i), blk_im(i)%vx, nv)
         elseif (field == 'y') then
            call copy(B(1, i), blk_im(i)%vy, nv)
         else
            call copy(B(1, i), blk_im(i)%vz, nv)
         end if
      end do
      call dgemm('T', 'N', nblk, nblk, nv, &
                 1.0d0, Aw, lv, B, lv, dbeta, G_ii, nblk)

!  ── G_ri: <re, im> ──
!  B still holds im from G_ii above, so only Aw needs regathering
!  with the real parts.  (This is why the order rr→ii→ri→ir matters.)
      do i = 1, nblk
         if (field == 'x') then
            call copy(Aw(1, i), blk_re(i)%vx, nv)
            call col2(Aw(1, i), bm1s, nv)
         elseif (field == 'y') then
            call copy(Aw(1, i), blk_re(i)%vy, nv)
            call col2(Aw(1, i), bm1s, nv)
         else
            call copy(Aw(1, i), blk_re(i)%vz, nv)
            call col2(Aw(1, i), bm1s, nv)
         end if
      end do
!  B still has im from G_ii — no regather needed
      call dgemm('T', 'N', nblk, nblk, nv, &
                 1.0d0, Aw, lv, B, lv, dbeta, G_ri, nblk)

!  ── G_ir: <im, re> — regather both Aw (im) and B (re) ──
!  Neither buffer has the right data: Aw was just overwritten with re
!  for G_ri, and B has im.  Both must be regathered.
      do i = 1, nblk
         if (field == 'x') then
            call copy(Aw(1, i), blk_im(i)%vx, nv)
            call col2(Aw(1, i), bm1s, nv)
         elseif (field == 'y') then
            call copy(Aw(1, i), blk_im(i)%vy, nv)
            call col2(Aw(1, i), bm1s, nv)
         else
            call copy(Aw(1, i), blk_im(i)%vz, nv)
            call col2(Aw(1, i), bm1s, nv)
         end if
      end do
      do i = 1, nblk
         if (field == 'x') then
            call copy(B(1, i), blk_re(i)%vx, nv)
         elseif (field == 'y') then
            call copy(B(1, i), blk_re(i)%vy, nv)
         else
            call copy(B(1, i), blk_re(i)%vz, nv)
         end if
      end do
      call dgemm('T', 'N', nblk, nblk, nv, &
                 1.0d0, Aw, lv, B, lv, dbeta, G_ir, nblk)

   end subroutine gram4_field_vel

!-----------------------------------------------------------------------
! gram4_field_temp — Accumulate 4 Gram sub-matrices for one
!                    temperature or passive scalar field
!
! Same logic as gram4_field_vel but uses lt/nt instead of lv/nv,
! and always accumulates (beta=1) because temperature is never the
! first field — velocity was already processed.
!
! iscal: scalar index into krylov_vector%t(:,iscal).
!        iscal=1 is temperature, iscal>1 are passive scalars.
!-----------------------------------------------------------------------
   subroutine gram4_field_temp(G_rr, G_ii, G_ri, G_ir, &
                               blk_re, blk_im, nblk, Aw, B, iscal)
      integer, intent(in) :: nblk, iscal
      real, intent(inout) :: G_rr(nblk, nblk), G_ii(nblk, nblk)
      real, intent(inout) :: G_ri(nblk, nblk), G_ir(nblk, nblk)
      type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)
!  Workspace — leading dimension lt (not lv), rows 1:nt used.
      real, intent(out) :: Aw(lt, nblk), B(lt, nblk)
      integer :: i

!  All dgemm calls below use beta=1.0 (accumulate) because velocity
!  components were already written into G_rr/G_ii/G_ri/G_ir.
!  Same gather order as gram4_field_vel: rr → ii → ri → ir.

!  G_rr += <re, re>
      do i = 1, nblk
         call copy(Aw(1, i), blk_re(i)%t(1, iscal), nt)
         call col2(Aw(1, i), bm1s, nt)
      end do
      do i = 1, nblk
         call copy(B(1, i), blk_re(i)%t(1, iscal), nt)
      end do
      call dgemm('T', 'N', nblk, nblk, nt, &
                 1.0d0, Aw, lt, B, lt, 1.0d0, G_rr, nblk)

      ! G_ii += <im, im>
      do i = 1, nblk
         call copy(Aw(1, i), blk_im(i)%t(1, iscal), nt)
         call col2(Aw(1, i), bm1s, nt)
      end do
      do i = 1, nblk
         call copy(B(1, i), blk_im(i)%t(1, iscal), nt)
      end do
      call dgemm('T', 'N', nblk, nblk, nt, &
                 1.0d0, Aw, lt, B, lt, 1.0d0, G_ii, nblk)

!  G_ri += <re, im> — regather Aw with re; B still has im from G_ii
      do i = 1, nblk
         call copy(Aw(1, i), blk_re(i)%t(1, iscal), nt)
         call col2(Aw(1, i), bm1s, nt)
      end do
      call dgemm('T', 'N', nblk, nblk, nt, &
                 1.0d0, Aw, lt, B, lt, 1.0d0, G_ri, nblk)

!  G_ir += <im, re> — regather both Aw (im) and B (re)
      do i = 1, nblk
         call copy(Aw(1, i), blk_im(i)%t(1, iscal), nt)
         call col2(Aw(1, i), bm1s, nt)
      end do
      do i = 1, nblk
         call copy(B(1, i), blk_re(i)%t(1, iscal), nt)
      end do
      call dgemm('T', 'N', nblk, nblk, nt, &
                 1.0d0, Aw, lt, B, lt, 1.0d0, G_ir, nblk)

   end subroutine gram4_field_temp

end module nekstab_krylov_inner_products
