! Read one field from standard input and print the mask result.
! Line 1: n. Line 2: n weights. Line 3: n field values.
! Line 1 of the output: 1 when any weight is zero, else 0.
! Then one line per point: zero-flag, bits before, bits after.
! A second pass must leave the masked field unchanged.
program bm1s_mask_driver
   use iso_fortran_env, only: int64
   use nekstab_bm1s_mask, only: weight_has_exact_zero, zero_where_weight_zero
   implicit none
   integer :: n, i, ios
   integer(int64) :: bits_before, bits_after, bits_again
   real, allocatable :: weight(:), field(:), original(:), again(:)
   logical :: has_zero

   read (*, *, iostat=ios) n
   if (ios /= 0) stop 1
   if (n <= 0) stop 1
   allocate (weight(n), field(n), original(n), again(n))
   read (*, *) weight
   read (*, *) field
   original = field

   has_zero = weight_has_exact_zero(weight, n)
   call zero_where_weight_zero(field, weight, n)
   again = field
   call zero_where_weight_zero(again, weight, n)

   write (6, '(i1)') merge(1, 0, has_zero)
   do i = 1, n
      bits_before = transfer(original(i), 0_int64)
      bits_after = transfer(field(i), 0_int64)
      bits_again = transfer(again(i), 0_int64)
      if (bits_again /= bits_after) stop 2
      write (6, '(i1,1x,z16.16,1x,z16.16)') merge(1, 0, weight(i) == 0.0d0), &
         bits_before, bits_after
   end do
end program bm1s_mask_driver
