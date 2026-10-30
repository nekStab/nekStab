!-----------------------------------------------------------------------
! bm1s_mask.f90 — Zero a field where its weight is zero
!
! Purpose:
!   The Krylov norm weights each point by bm1s. A zero weight does not
!   enter the norm. The routines have no dependency on Nek5000, so a
!   stand-alone test can compile them.
!
! Public interface:
!   weight_has_exact_zero  — true when any weight entry is zero
!   zero_where_weight_zero — set the field to zero at those points
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_bm1s_mask
   implicit none
   private
   public :: weight_has_exact_zero, zero_where_weight_zero

contains

!-----------------------------------------------------------------------
! weight_has_exact_zero — Return true when any weight entry is zero
!
!   The test is exact equality with zero.
!-----------------------------------------------------------------------
   logical function weight_has_exact_zero(weight, n)
      integer, intent(in) :: n
      real, intent(in) :: weight(*)

      integer :: i

      weight_has_exact_zero = .false.
      if (n <= 0) return
      do i = 1, n
         if (weight(i) == 0.0d0) then
            weight_has_exact_zero = .true.
            return
         end if
      end do

   end function weight_has_exact_zero

!-----------------------------------------------------------------------
! zero_where_weight_zero — Set the field to zero where the weight is zero
!
!   A point with a nonzero weight keeps its value.
!   The stored weight supplies the zero test. Do not rebuild the fringe.
!-----------------------------------------------------------------------
   subroutine zero_where_weight_zero(field, weight, n)
      integer, intent(in) :: n
      real, intent(inout) :: field(*)
      real, intent(in) :: weight(*)

      integer :: i

      if (n <= 0) return
      do i = 1, n
         if (weight(i) == 0.0d0) field(i) = 0.0d0
      end do

   end subroutine zero_where_weight_zero

end module nekstab_bm1s_mask
