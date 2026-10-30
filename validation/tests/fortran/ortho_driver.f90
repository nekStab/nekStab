! inner: n, then n inner products. Prints the mean absolute value.
! length: kk. Prints the last-pass count.
! basis: m, k_dim, kk, then the basis in column-major order.
! Prints the last-pass column and its mean absolute inner product
! with the earlier columns of that pass.
program ortho_driver
   use nekstab_ortho, only: last_pass_count, mean_abs_inner
   implicit none
   character(len=16) :: cmd
   integer :: n, kk, k_dim, m, col, j, ios
   real :: metric
   real, allocatable :: dot(:), vec(:, :)

   read (*, *, iostat=ios) cmd
   if (ios /= 0) stop 1

   if (cmd == 'inner') then
      read (*, *) n
      if (n < 0) stop 1
      allocate (dot(max(n, 1)))
      if (n > 0) read (*, *) dot(1:n)
      metric = mean_abs_inner(dot, n)
      write (6, '(es24.16)') metric
   else if (cmd == 'length') then
      read (*, *) kk
      write (6, '(i8)') last_pass_count(kk)
   else if (cmd == 'basis') then
      read (*, *) m, k_dim, kk
      if (m <= 0 .or. k_dim <= 0 .or. kk <= 0 .or. kk > k_dim) stop 1
      allocate (vec(m, k_dim), dot(k_dim))
      read (*, *) vec
      col = last_pass_count(kk)
      if (col < 2) then
         metric = 0.0d0
      else
         do j = 1, col - 1
            dot(j) = dot_product(vec(:, col), vec(:, j))
         end do
         metric = mean_abs_inner(dot, col - 1)
      end if
      write (6, '(i8,1x,es24.16)') col, metric
   else
      stop 1
   end if
end program ortho_driver
