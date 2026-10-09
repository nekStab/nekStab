!-----------------------------------------------------------------------
! ks_select.f90 — Ritz value selection for the Krylov-Schur restart
!
! Purpose:
!   Chooses which Ritz values the restart keeps. The routines have no
!   dependency on Nek5000, so a stand-alone test can compile them.
!
! Public interface:
!   ks_select_eigenvalues     — adaptive selection with a pair-aware cap
!   ks_ensure_conjugate_pairs — select both members of a conjugate pair
!
! Dependencies:
!   nekstab_argsort
!-----------------------------------------------------------------------

module nekstab_ks_select
   use nekstab_argsort, only: argsort
   implicit none
   private
   public :: ks_select_eigenvalues, ks_ensure_conjugate_pairs

   integer, parameter :: dp = kind(0.0d0)
   real, parameter :: PAIR_REL_TOL = 1.0d-10

contains

!-----------------------------------------------------------------------
! ks_select_eigenvalues — Adaptive selection for the restart
!
!   Three criteria select Ritz values:
!     1. Near the unit circle: |lambda| >= 1 - delta.
!     2. Partly converged: Schur residual < sqrt(tol). Discarding these
!        wastes the matvecs already spent on them.
!     3. At least nev + 2 by magnitude.
!   Conjugate pairs stay together (real Schur form). The result holds at
!   most n - nev values, so the next Arnoldi run adds at least nev steps.
!   Over the limit, whole blocks (a pair, or a real value) with the
!   smallest residual stay. A block is never split.
!
!   An older version removed single values by residual and then completed
!   the pairs. Completing a pair adds a value, so the count returned to
!   n. No Arnoldi step followed, and the residual test read the Schur
!   form, where only the last block has a nonzero last component. The run
!   then reported n - 2 converged eigenvalues.
!
!   delta    — magnitude selection radius (schur_del)
!   nev      — number of wanted eigenvalues (schur_tgt), at least 1
!   tol      — eigenvalue tolerance (eigen_tol)
!   n_circle — number selected by criterion 1
!   n_resid  — number added by criterion 2
!-----------------------------------------------------------------------
   subroutine ks_select_eigenvalues(selected, nsel, vals, residuals, delta, &
                                    nev, n, tol, n_circle, n_resid)
      integer, intent(in) :: nev, n
      complex(dp), dimension(n), intent(in) :: vals
      real, dimension(n), intent(in) :: residuals
      real, intent(in) :: delta, tol
      logical, dimension(n), intent(out) :: selected
      integer, intent(out) :: nsel, n_circle, n_resid

      integer :: i, min_keep, max_keep
      integer, dimension(n) :: idx
      real, dimension(n) :: work_arr
      real :: sqrt_tol

      sqrt_tol = sqrt(tol)
      selected = .false.

      !  Criterion 1: Ritz values near the unit circle.
      do i = 1, n
         if (abs(vals(i)) >= (1.0d0 - delta)) selected(i) = .true.
      end do
      n_circle = count(selected)

      !  Criterion 2: partly converged Ritz values.
      do i = 1, n
         if (residuals(i) < sqrt_tol) selected(i) = .true.
      end do
      n_resid = count(selected) - n_circle

      !  Criterion 3: at least nev + 2 by magnitude. Ascending sort, so the
      !  largest magnitudes are at the end.
      min_keep = nev + 2
      if (count(selected) < min_keep) then
         do i = 1, n
            idx(i) = i
         end do
         work_arr = abs(vals)
         call argsort(n, work_arr, idx)
         do i = n, 1, -1
            if (count(selected) >= min_keep) exit
            selected(idx(i)) = .true.
         end do
      end if

      call ks_ensure_conjugate_pairs(selected, vals, n)

      !  Cap at n - nev. Whole blocks with the smallest residual stay.
      max_keep = n - nev
      if (count(selected) > max_keep) call cap_selection(selected, vals, residuals, n, max_keep)

      nsel = count(selected)

   end subroutine ks_select_eigenvalues

!-----------------------------------------------------------------------
! cap_selection — Keep whole blocks with the smallest residual
!
!   A block is a conjugate pair (two entries) or a real value. The key of
!   a block is the largest residual of its entries. Blocks are taken in
!   order of the key while they fit under max_keep, so the count never
!   exceeds max_keep and no pair is split.
!-----------------------------------------------------------------------
   subroutine cap_selection(selected, vals, residuals, n, max_keep)
      integer, intent(in) :: n, max_keep
      complex(dp), dimension(n), intent(in) :: vals
      real, dimension(n), intent(in) :: residuals
      logical, dimension(n), intent(out) :: selected

      integer :: i, b, nblk, kept, j
      integer, dimension(n) :: first, blen, idx
      real, dimension(n) :: key

      nblk = 0
      i = 1
      do while (i <= n)
         nblk = nblk + 1
         first(nblk) = i
         if (is_conjugate_pair(vals, i, n)) then
            blen(nblk) = 2
            key(nblk) = max(residuals(i), residuals(i + 1))
         else
            blen(nblk) = 1
            key(nblk) = residuals(i)
         end if
         i = i + blen(nblk)
      end do

      do b = 1, nblk
         idx(b) = b
      end do
      call argsort(nblk, key, idx)

      selected = .false.
      kept = 0
      do b = 1, nblk
         j = idx(b)
         if (kept + blen(j) <= max_keep) then
            selected(first(j):first(j) + blen(j) - 1) = .true.
            kept = kept + blen(j)
         end if
      end do

   end subroutine cap_selection

!-----------------------------------------------------------------------
! is_conjugate_pair — True when vals(i) and vals(i+1) are a conjugate pair
!
!   The real Schur form of dgees stores a pair in two consecutive entries
!   with equal real parts and opposite imaginary parts.
!-----------------------------------------------------------------------
   logical function is_conjugate_pair(vals, i, n)
      integer, intent(in) :: i, n
      complex(dp), dimension(n), intent(in) :: vals
      real :: scale, tol

      is_conjugate_pair = .false.
      if (i >= n) return
      if (aimag(vals(i)) == 0.0d0) return
      scale = max(abs(vals(i)), 1.0d0)
      tol = PAIR_REL_TOL*scale
      is_conjugate_pair = abs(real(vals(i)) - real(vals(i + 1))) < tol .and. &
                          abs(aimag(vals(i)) + aimag(vals(i + 1))) < tol

   end function is_conjugate_pair

!-----------------------------------------------------------------------
! ks_ensure_conjugate_pairs — Select both members of a selected pair
!-----------------------------------------------------------------------
   subroutine ks_ensure_conjugate_pairs(selected, vals, n)
      integer, intent(in) :: n
      complex(dp), dimension(n), intent(in) :: vals
      logical, dimension(n), intent(inout) :: selected

      integer :: i

      i = 1
      do while (i < n)
         if (is_conjugate_pair(vals, i, n)) then
            if (selected(i) .or. selected(i + 1)) then
               selected(i) = .true.
               selected(i + 1) = .true.
            end if
            i = i + 2
         else
            i = i + 1
         end if
      end do

   end subroutine ks_ensure_conjugate_pairs

end module nekstab_ks_select
