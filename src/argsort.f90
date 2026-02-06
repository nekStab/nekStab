      !-----------------------------------------------------------------------
      ! argsort.f90 — Insertion sort with index tracking
      !
      ! Purpose:
      !   Sorts an array in ascending order while tracking original indices.
      !   Stable sort: equal elements maintain their relative order.
      !   O(n^2) but optimal for small arrays (n < 100).
      !
      ! Public interface:
      !   argsort(n, arr, idx)
      !
      ! Dependencies:
      !   None
      !-----------------------------------------------------------------------

      ! argsort — Sort array ascending with index tracking
      subroutine argsort(n, arr, idx)
      implicit none
      integer, intent(in) :: n
      real, intent(inout) :: arr(n)
      integer, intent(inout) :: idx(n)
      integer :: i, j
      real :: key
      integer :: key_idx

      do i = 2, n
         key = arr(i)
         key_idx = idx(i)
         j = i - 1
         do while (j >= 1)
            if (arr(j) <= key) exit
            arr(j+1) = arr(j)
            idx(j+1) = idx(j)
            j = j - 1
         end do
         arr(j+1) = key
         idx(j+1) = key_idx
      end do

      end subroutine argsort
