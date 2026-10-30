! Drive the scalar Fourier series.
! accum: N K, then N samples. One line per mode: c_re c_im.
! eval: K T ntime, then K+1 real parts, K+1 imag parts, ntime times.
!        One value per line.
! select: Kacc kmax_fixed, then Kacc+1 energies. Line 1: K missing ntrace.
!         Then one line per doubling: K Knext missing.
! kaccum: nsamples kmax_fixed. One integer.
program fourier_series_driver
   use nekstab_fourier_series, only: fourier_dft_add, fourier_eval_series, &
                                     fourier_k_accum, fourier_select_k
   implicit none
   character(len=16) :: cmd
   integer :: ios

   read (*, *, iostat=ios) cmd
   if (ios /= 0) stop 1
   if (cmd == 'accum') then
      call run_accum()
   else if (cmd == 'eval') then
      call run_eval()
   else if (cmd == 'select') then
      call run_select()
   else if (cmd == 'kaccum') then
      call run_kaccum()
   else
      stop 1
   end if

contains

   subroutine run_accum()
      integer :: n, k, j
      real, allocatable :: samples(:), c_re(:), c_im(:)

      read (*, *) n, k
      if (n < 1 .or. k < 0) stop 1
      allocate (samples(n), c_re(0:k), c_im(0:k))
      read (*, *) samples
      c_re = 0.0d0
      c_im = 0.0d0
      do j = 1, n
         call fourier_dft_add(c_re, c_im, k, samples(j), j, n)
      end do
      do j = 0, k
         write (6, '(2es24.16)') c_re(j), c_im(j)
      end do
   end subroutine run_accum

   subroutine run_eval()
      integer :: k, ntime, j
      real :: period, t
      real, allocatable :: c_re(:), c_im(:), times(:)

      read (*, *) k, period, ntime
      if (k < 0 .or. period == 0.0d0 .or. ntime < 1) stop 1
      allocate (c_re(0:k), c_im(0:k), times(ntime))
      read (*, *) c_re
      read (*, *) c_im
      read (*, *) times
      do j = 1, ntime
         t = fourier_eval_series(c_re, c_im, k, times(j), period)
         write (6, '(es24.16)') t
      end do
   end subroutine run_eval

   subroutine run_select()
      integer :: kacc, kmax_fixed, k_sel, ntrace, j
      integer :: trace_k(32), trace_next(32)
      real :: missing
      real :: trace_miss(32)
      real, allocatable :: ek(:)

      read (*, *) kacc, kmax_fixed
      if (kacc < 0) stop 1
      allocate (ek(0:kacc))
      read (*, *) ek
      call fourier_select_k(ek, kacc, kmax_fixed, k_sel, missing, &
                            ntrace, trace_k, trace_next, trace_miss)
      write (6, '(i6,1x,es24.16,1x,i6)') k_sel, missing, ntrace
      do j = 1, ntrace
         write (6, '(2i6,1x,es24.16)') trace_k(j), trace_next(j), trace_miss(j)
      end do
   end subroutine run_select

   subroutine run_kaccum()
      integer :: nsamples, kmax_fixed

      read (*, *) nsamples, kmax_fixed
      write (6, '(i6)') fourier_k_accum(nsamples, kmax_fixed)
   end subroutine run_kaccum

end program fourier_series_driver
