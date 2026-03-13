!-----------------------------------------------------------------------
! krylov_inner_products.f90 — Batch inner product operations
!
! Purpose:
!   Provides batch inner product operations on arrays of
!   krylov_vector using BLAS dgemm/dgemv with a single MPI
!   reduction (gop), replacing O(n^2) individual glsc3 calls
!   that each trigger a separate MPI_Allreduce.
!
! Public interface:
!   k_gram_matrix  — G(i,j) = <P(i), Q(j)>_bm1s
!   k_project      — h(i)   = <Q(i), f>_bm1s
!   k_gram_complex — CSD(i,j) = <P_re(i)+iP_im(i), Q_re(j)+iQ_im(j)>
!
! Dependencies:
!   krylov_subspace
!-----------------------------------------------------------------------

module krylov_inner_products
   use krylov_subspace
   use nekstab_nek_bridge
   implicit none
   private

   public :: k_gram_matrix, k_project, k_gram_complex

contains

!-----------------------------------------------------------------------
! k_gram_matrix — Batch Gram matrix: G(i,j) = <P(i), Q(j)>_bm1s
!
! Algorithm:
!   For each field (vx, vy, [vz], [t]):
!     1. Gather left vectors into Aw(:,i) = P(i)%field * bm1s
!     2. Gather right vectors into B(:,j) = Q(j)%field
!     3. dgemm: G += Aw^T * B  (beta=0 for vx, beta=1 after)
!   Single gop at the end.
!-----------------------------------------------------------------------
subroutine k_gram_matrix(G, P, Q, m, n, ldg)
   integer, intent(in) :: m, n, ldg
   real, intent(out) :: G(ldg, n)
   type(krylov_vector), intent(in) :: P(m), Q(n)
   real, allocatable :: Aw(:,:), B(:,:), wk_gop(:)
   integer :: i, j, mm

   nv = nx1*ny1*nz1*nelv
   nt = nx1*ny1*nz1*nelt

   ! Allocate velocity-sized workspace
   allocate(Aw(lv, m), B(lv, n), wk_gop(m*n))

   ! ── vx (beta=0 initializes G) ──
   do i = 1, m
      call copy(Aw(1,i), P(i)%vx, nv)
      call col2(Aw(1,i), bm1s, nv)
   end do
   do j = 1, n
      call copy(B(1,j), Q(j)%vx, nv)
   end do
   call dgemm('T','N', m, n, nv, &
        1.0d0, Aw, lv, B, lv, 0.0d0, G, ldg)

   ! ── vy (accumulate, beta=1) ──
   do i = 1, m
      call copy(Aw(1,i), P(i)%vy, nv)
      call col2(Aw(1,i), bm1s, nv)
   end do
   do j = 1, n
      call copy(B(1,j), Q(j)%vy, nv)
   end do
   call dgemm('T','N', m, n, nv, &
        1.0d0, Aw, lv, B, lv, 1.0d0, G, ldg)

   ! ── vz (3D only, accumulate) ──
   if (if3d) then
      do i = 1, m
         call copy(Aw(1,i), P(i)%vz, nv)
         call col2(Aw(1,i), bm1s, nv)
      end do
      do j = 1, n
         call copy(B(1,j), Q(j)%vz, nv)
      end do
      call dgemm('T','N', m, n, nv, &
           1.0d0, Aw, lv, B, lv, 1.0d0, G, ldg)
   end if

   deallocate(Aw, B)

   ! ── Temperature fields (lt-sized workspace) ──
   if (ifto .or. (ldimt > 1)) then
      allocate(Aw(lt, m), B(lt, n))

      if (ifto) then
         do i = 1, m
            call copy(Aw(1,i), P(i)%t(1,1), nt)
            call col2(Aw(1,i), bm1s, nt)
         end do
         do j = 1, n
            call copy(B(1,j), Q(j)%t(1,1), nt)
         end do
         call dgemm('T','N', m, n, nt, &
              1.0d0, Aw, lt, B, lt, 1.0d0, G, ldg)
      end if

      if (ldimt > 1) then
         do mm = 2, ldimt
            if (ifpsco(mm-1)) then
               do i = 1, m
                  call copy(Aw(1,i), P(i)%t(1,mm), nt)
                  call col2(Aw(1,i), bm1s, nt)
               end do
               do j = 1, n
                  call copy(B(1,j), Q(j)%t(1,mm), nt)
               end do
               call dgemm('T','N', m, n, nt, &
                    1.0d0, Aw, lt, B, lt, 1.0d0, G, ldg)
            end if
         end do
      end if

      deallocate(Aw, B)
   end if

   ! ── Single global reduction ──
   call gop(G, wk_gop, '+  ', m*n)
   deallocate(wk_gop)

   ! ── Newton PO time component (scalar, no MPI — added after gop) ──
   if (isNewtonPO) then
      do j = 1, n
         do i = 1, m
            G(i,j) = G(i,j) + P(i)%time * Q(j)%time
         end do
      end do
   end if

end subroutine k_gram_matrix

!-----------------------------------------------------------------------
! k_project — Batch projection: h(i) = <Q(i), f>_bm1s
!
! Same as one column of k_gram_matrix but uses dgemv.
!-----------------------------------------------------------------------
subroutine k_project(h, Q, f, k)
   integer, intent(in) :: k
   real, intent(out) :: h(k)
   type(krylov_vector), dimension(k), intent(in) :: Q
   type(krylov_vector), intent(in) :: f
   real, allocatable :: Aw(:,:), fwrk(:), wk_gop(:)
   integer :: i, mm

   nv = nx1*ny1*nz1*nelv
   nt = nx1*ny1*nz1*nelt

   ! Allocate workspace
   allocate(Aw(lv, k), fwrk(lv), wk_gop(k))

   ! ── vx (beta=0 initializes h) ──
   do i = 1, k
      call copy(Aw(1,i), Q(i)%vx, nv)
      call col2(Aw(1,i), bm1s, nv)
   end do
   call copy(fwrk, f%vx, nv)
   call dgemv('T', nv, k, 1.0d0, &
        Aw, lv, fwrk, 1, 0.0d0, h, 1)

   ! ── vy (accumulate) ──
   do i = 1, k
      call copy(Aw(1,i), Q(i)%vy, nv)
      call col2(Aw(1,i), bm1s, nv)
   end do
   call copy(fwrk, f%vy, nv)
   call dgemv('T', nv, k, 1.0d0, &
        Aw, lv, fwrk, 1, 1.0d0, h, 1)

   ! ── vz (3D only) ──
   if (if3d) then
      do i = 1, k
         call copy(Aw(1,i), Q(i)%vz, nv)
         call col2(Aw(1,i), bm1s, nv)
      end do
      call copy(fwrk, f%vz, nv)
      call dgemv('T', nv, k, 1.0d0, &
           Aw, lv, fwrk, 1, 1.0d0, h, 1)
   end if

   deallocate(Aw, fwrk)

   ! ── Temperature fields ──
   if (ifto .or. (ldimt > 1)) then
      allocate(Aw(lt, k), fwrk(lt))

      if (ifto) then
         do i = 1, k
            call copy(Aw(1,i), Q(i)%t(1,1), nt)
            call col2(Aw(1,i), bm1s, nt)
         end do
         call copy(fwrk, f%t(1,1), nt)
         call dgemv('T', nt, k, 1.0d0, &
              Aw, lt, fwrk, 1, 1.0d0, h, 1)
      end if

      if (ldimt > 1) then
         do mm = 2, ldimt
            if (ifpsco(mm-1)) then
               do i = 1, k
                  call copy(Aw(1,i), Q(i)%t(1,mm), nt)
                  call col2(Aw(1,i), bm1s, nt)
               end do
               call copy(fwrk, f%t(1,mm), nt)
               call dgemv('T', nt, k, 1.0d0, &
                    Aw, lt, fwrk, 1, 1.0d0, h, 1)
            end if
         end do
      end if

      deallocate(Aw, fwrk)
   end if

   ! ── Single global reduction ──
   call gop(h, wk_gop, '+  ', k)
   deallocate(wk_gop)

   ! ── Newton PO time component (scalar, no MPI — added after gop) ──
   if (isNewtonPO) then
      do i = 1, k
         h(i) = h(i) + Q(i)%time * f%time
      end do
   end if

end subroutine k_project

!-----------------------------------------------------------------------
! k_gram_complex — Hermitian CSD matrix for complex krylov_vectors
!
!   CSD(i,j) = <P_re(i)+iP_im(i), Q_re(j)+iQ_im(j)>
!            = (G_rr + G_ii) + i(G_ri - G_ir)
!
! Computes 4 real Gram matrices contiguously, single gop over all,
! then assembles complex CSD.
!-----------------------------------------------------------------------
subroutine k_gram_complex(CSD, blk_re, blk_im, nblk, ldcsd)
   integer, intent(in) :: nblk, ldcsd
   complex(nekStab_dp), intent(out) :: CSD(ldcsd, nblk)
   type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)

   real, allocatable :: G_rr(:,:), G_ii(:,:)
   real, allocatable :: G_ri(:,:), G_ir(:,:)
   real, allocatable :: G_all(:), wk_gop(:)
   real, allocatable :: Aw(:,:), B(:,:)
   integer :: i, j, mm, off

   nv = nx1*ny1*nz1*nelv
   nt = nx1*ny1*nz1*nelt

   ! 4 Gram matrices packed contiguously for single gop
   allocate(G_all(4*nblk*nblk), wk_gop(4*nblk*nblk))
   allocate(G_rr(nblk,nblk), G_ii(nblk,nblk))
   allocate(G_ri(nblk,nblk), G_ir(nblk,nblk))

   G_rr = 0.0d0; G_ii = 0.0d0
   G_ri = 0.0d0; G_ir = 0.0d0

   ! Velocity-sized workspace (reused for each component)
   allocate(Aw(lv, nblk), B(lv, nblk))

   ! ── Process each velocity field ──
   ! vx
   call gram4_field_vel(G_rr, G_ii, G_ri, G_ir, &
        blk_re, blk_im, nblk, Aw, B, 'x', 0)
   ! vy
   call gram4_field_vel(G_rr, G_ii, G_ri, G_ir, &
        blk_re, blk_im, nblk, Aw, B, 'y', 1)
   ! vz
   if (if3d) then
      call gram4_field_vel(G_rr, G_ii, G_ri, G_ir, &
           blk_re, blk_im, nblk, Aw, B, 'z', 1)
   end if

   deallocate(Aw, B)

   ! ── Temperature fields ──
   if (ifto .or. (ldimt > 1)) then
      allocate(Aw(lt, nblk), B(lt, nblk))

      if (ifto) then
         call gram4_field_temp(G_rr, G_ii, G_ri, G_ir, &
              blk_re, blk_im, nblk, Aw, B, 1)
      end if
      if (ldimt > 1) then
         do mm = 2, ldimt
            if (ifpsco(mm-1)) then
               call gram4_field_temp(G_rr, G_ii, G_ri, G_ir, &
                    blk_re, blk_im, nblk, Aw, B, mm)
            end if
         end do
      end if

      deallocate(Aw, B)
   end if

   ! ── Pack into contiguous array for single gop ──
   off = 0
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_all(off) = G_rr(i,j)
      end do
   end do
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_all(off) = G_ii(i,j)
      end do
   end do
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_all(off) = G_ri(i,j)
      end do
   end do
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_all(off) = G_ir(i,j)
      end do
   end do

   ! ── Single global reduction ──
   call gop(G_all, wk_gop, '+  ', 4*nblk*nblk)

   ! ── Unpack ──
   off = 0
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_rr(i,j) = G_all(off)
      end do
   end do
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_ii(i,j) = G_all(off)
      end do
   end do
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_ri(i,j) = G_all(off)
      end do
   end do
   do j = 1, nblk
      do i = 1, nblk
         off = off + 1
         G_ir(i,j) = G_all(off)
      end do
   end do

   ! ── Assemble complex CSD ──
   do j = 1, nblk
      do i = 1, nblk
         CSD(i,j) = dcmplx(G_rr(i,j) + G_ii(i,j), &
              G_ri(i,j) - G_ir(i,j))
      end do
   end do

   deallocate(G_all, wk_gop, G_rr, G_ii, G_ri, G_ir)

end subroutine k_gram_complex

!-----------------------------------------------------------------------
! gram4_field_vel — Helper: accumulate 4 Gram sub-matrices for one
!                   velocity field (vx/vy/vz)
!
!   beta_flag: 0 = initialize (beta=0), 1 = accumulate (beta=1)
!-----------------------------------------------------------------------
subroutine gram4_field_vel(G_rr, G_ii, G_ri, G_ir, &
     blk_re, blk_im, nblk, Aw, B, field, beta_flag)
   integer, intent(in) :: nblk, beta_flag
   real, intent(inout) :: G_rr(nblk,nblk), G_ii(nblk,nblk)
   real, intent(inout) :: G_ri(nblk,nblk), G_ir(nblk,nblk)
   type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)
   real :: Aw(lv, nblk), B(lv, nblk)
   character(len=1), intent(in) :: field
   real :: dbeta
   integer :: i

   nv = nx1*ny1*nz1*nelv
   dbeta = dble(beta_flag)

   ! ── G_rr: <re, re> ──
   do i = 1, nblk
      if (field == 'x') then
         call copy(Aw(1,i), blk_re(i)%vx, nv)
         call col2(Aw(1,i), bm1s, nv)
      elseif (field == 'y') then
         call copy(Aw(1,i), blk_re(i)%vy, nv)
         call col2(Aw(1,i), bm1s, nv)
      else
         call copy(Aw(1,i), blk_re(i)%vz, nv)
         call col2(Aw(1,i), bm1s, nv)
      end if
   end do
   do i = 1, nblk
      if (field == 'x') then
         call copy(B(1,i), blk_re(i)%vx, nv)
      elseif (field == 'y') then
         call copy(B(1,i), blk_re(i)%vy, nv)
      else
         call copy(B(1,i), blk_re(i)%vz, nv)
      end if
   end do
   call dgemm('T','N', nblk, nblk, nv, &
        1.0d0, Aw, lv, B, lv, dbeta, G_rr, nblk)

   ! ── G_ii: <im, im> ──
   do i = 1, nblk
      if (field == 'x') then
         call copy(Aw(1,i), blk_im(i)%vx, nv)
         call col2(Aw(1,i), bm1s, nv)
      elseif (field == 'y') then
         call copy(Aw(1,i), blk_im(i)%vy, nv)
         call col2(Aw(1,i), bm1s, nv)
      else
         call copy(Aw(1,i), blk_im(i)%vz, nv)
         call col2(Aw(1,i), bm1s, nv)
      end if
   end do
   do i = 1, nblk
      if (field == 'x') then
         call copy(B(1,i), blk_im(i)%vx, nv)
      elseif (field == 'y') then
         call copy(B(1,i), blk_im(i)%vy, nv)
      else
         call copy(B(1,i), blk_im(i)%vz, nv)
      end if
   end do
   call dgemm('T','N', nblk, nblk, nv, &
        1.0d0, Aw, lv, B, lv, dbeta, G_ii, nblk)

   ! ── G_ri: <re, im> (regather Aw=re*bm1s, overwritten by G_ii) ──
   do i = 1, nblk
      if (field == 'x') then
         call copy(Aw(1,i), blk_re(i)%vx, nv)
         call col2(Aw(1,i), bm1s, nv)
      elseif (field == 'y') then
         call copy(Aw(1,i), blk_re(i)%vy, nv)
         call col2(Aw(1,i), bm1s, nv)
      else
         call copy(Aw(1,i), blk_re(i)%vz, nv)
         call col2(Aw(1,i), bm1s, nv)
      end if
   end do
   ! B already has im from G_ii
   call dgemm('T','N', nblk, nblk, nv, &
        1.0d0, Aw, lv, B, lv, dbeta, G_ri, nblk)

   ! ── G_ir: <im, re> ──
   ! Aw = im*bm1s (already gathered for G_ii)
   do i = 1, nblk
      if (field == 'x') then
         call copy(Aw(1,i), blk_im(i)%vx, nv)
         call col2(Aw(1,i), bm1s, nv)
      elseif (field == 'y') then
         call copy(Aw(1,i), blk_im(i)%vy, nv)
         call col2(Aw(1,i), bm1s, nv)
      else
         call copy(Aw(1,i), blk_im(i)%vz, nv)
         call col2(Aw(1,i), bm1s, nv)
      end if
   end do
   do i = 1, nblk
      if (field == 'x') then
         call copy(B(1,i), blk_re(i)%vx, nv)
      elseif (field == 'y') then
         call copy(B(1,i), blk_re(i)%vy, nv)
      else
         call copy(B(1,i), blk_re(i)%vz, nv)
      end if
   end do
   call dgemm('T','N', nblk, nblk, nv, &
        1.0d0, Aw, lv, B, lv, dbeta, G_ir, nblk)

end subroutine gram4_field_vel

!-----------------------------------------------------------------------
! gram4_field_temp — Helper: accumulate 4 Gram sub-matrices for one
!                    temperature/scalar field
!-----------------------------------------------------------------------
subroutine gram4_field_temp(G_rr, G_ii, G_ri, G_ir, &
     blk_re, blk_im, nblk, Aw, B, iscal)
   integer, intent(in) :: nblk, iscal
   real, intent(inout) :: G_rr(nblk,nblk), G_ii(nblk,nblk)
   real, intent(inout) :: G_ri(nblk,nblk), G_ir(nblk,nblk)
   type(krylov_vector), intent(in) :: blk_re(nblk), blk_im(nblk)
   real :: Aw(lt, nblk), B(lt, nblk)
   integer :: i

   nt = nx1*ny1*nz1*nelt

   ! G_rr += <re, re>
   do i = 1, nblk
      call copy(Aw(1,i), blk_re(i)%t(1,iscal), nt)
      call col2(Aw(1,i), bm1s, nt)
   end do
   do i = 1, nblk
      call copy(B(1,i), blk_re(i)%t(1,iscal), nt)
   end do
   call dgemm('T','N', nblk, nblk, nt, &
        1.0d0, Aw, lt, B, lt, 1.0d0, G_rr, nblk)

   ! G_ii += <im, im>
   do i = 1, nblk
      call copy(Aw(1,i), blk_im(i)%t(1,iscal), nt)
      call col2(Aw(1,i), bm1s, nt)
   end do
   do i = 1, nblk
      call copy(B(1,i), blk_im(i)%t(1,iscal), nt)
   end do
   call dgemm('T','N', nblk, nblk, nt, &
        1.0d0, Aw, lt, B, lt, 1.0d0, G_ii, nblk)

   ! G_ri += <re, im>
   do i = 1, nblk
      call copy(Aw(1,i), blk_re(i)%t(1,iscal), nt)
      call col2(Aw(1,i), bm1s, nt)
   end do
   ! B already has im from G_ii
   call dgemm('T','N', nblk, nblk, nt, &
        1.0d0, Aw, lt, B, lt, 1.0d0, G_ri, nblk)

   ! G_ir += <im, re>
   do i = 1, nblk
      call copy(Aw(1,i), blk_im(i)%t(1,iscal), nt)
      call col2(Aw(1,i), bm1s, nt)
   end do
   do i = 1, nblk
      call copy(B(1,i), blk_re(i)%t(1,iscal), nt)
   end do
   call dgemm('T','N', nblk, nblk, nt, &
        1.0d0, Aw, lt, B, lt, 1.0d0, G_ir, nblk)

end subroutine gram4_field_temp

end module krylov_inner_products
