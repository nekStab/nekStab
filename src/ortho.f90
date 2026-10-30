!-----------------------------------------------------------------------
! ortho.f90 — Orthogonality metric of a Krylov basis
!
! Purpose:
!   The mean |<Q(k), Q(j)>| printed by GMRES, and the vector count
!   of the last Arnoldi pass. The routines have no dependency on
!   Nek5000, so a stand-alone test can compile them.
!
! Public interface:
!   mean_abs_inner  — mean absolute inner product with earlier vectors
!   last_pass_count — vectors the last Arnoldi pass built
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_ortho
   implicit none
   private
   public :: mean_abs_inner, last_pass_count

contains

!-----------------------------------------------------------------------
! mean_abs_inner — Mean of |dot(j)| for j = 1 .. n
!
!   dot(j) is <Q(k), Q(j)> for the earlier vectors. n is k-1.
!   The mean is 0 when k is 1.
!-----------------------------------------------------------------------
   pure function mean_abs_inner(dot, n) result(metric)
      integer, intent(in) :: n
      real, intent(in) :: dot(*)
      real :: metric

      integer :: j

      metric = 0.0d0
      if (n <= 0) return
      do j = 1, n
         metric = metric + abs(dot(j))
      end do
      metric = metric/n

   end function mean_abs_inner

!-----------------------------------------------------------------------
! last_pass_count — Vector count of the last Arnoldi pass
!
!   kk is that count. Columns after kk belong to an earlier pass.
!-----------------------------------------------------------------------
   pure integer function last_pass_count(kk) result(n)
      integer, intent(in) :: kk

      n = kk

   end function last_pass_count

end module nekstab_ortho
