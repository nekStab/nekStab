!-----------------------------------------------------------------------
! fourier_series.f90 — Fourier coefficients of one scalar series
!
! Purpose:
!   The DFT weights, the mode count and the evaluation of
!     u(t) = c_0 + 2 Re sum_{k=1}^{K} c_k exp(2 pi i k t / T)
!   The routines have no dependency on Nek5000, so a stand-alone
!   test can compile them. The orbit in fourier.f90 calls them.
!
! Public interface:
!   fourier_k_accum    — mode count accumulated over one period
!   fourier_dft_weight — weight of sample istep in mode k
!   fourier_dft_add    — add one sample to the coefficients
!   fourier_eval_weight — weight of mode k at time t
!   fourier_eval_series — evaluate the stored series at time t
!   fourier_select_k   — keep modes 0..K by the doubling rule
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_fourier_series
   implicit none
   private
   public :: fourier_k_accum, fourier_dft_weight, fourier_dft_add, &
             fourier_eval_weight, fourier_eval_series, fourier_select_k

   ! Most modes accumulated when the caller does not fix K.
   ! Below the Nyquist mode every c_k has a distinct conjugate,
   ! so the factor 2 in u(t) holds for all stored k.
   integer, parameter :: KCAP_DEFAULT = 128
   ! K starts here and doubles while the missing energy falls by 10.
   integer, parameter :: K_FLOOR = 8
   real, parameter :: MISSING_TARGET = 1.0d-10
   real, parameter :: TWOPI = 8.0d0*atan(1.0d0)

contains

!-----------------------------------------------------------------------
! fourier_k_accum — Mode count accumulated over one period
!
!   nsamples steps cover one period. kmax_fixed > 0 returns
!   that count. The count never includes the Nyquist mode.
!-----------------------------------------------------------------------
   pure integer function fourier_k_accum(nsamples, kmax_fixed) result(k)
      integer, intent(in) :: nsamples, kmax_fixed

      k = min(KCAP_DEFAULT, (nsamples - 1)/2)
      if (kmax_fixed > 0) k = min(kmax_fixed, (nsamples - 1)/2)

   end function fourier_k_accum

!-----------------------------------------------------------------------
! fourier_dft_weight — Weight of sample istep in mode k
!
!   Samples are U(dt), ..., U(N dt), with N dt = T.
!   Mode 0 takes the real weight 1/N.
!   Mode k returns (1/N) exp(-2 pi i k istep / N).
!-----------------------------------------------------------------------
   pure subroutine fourier_dft_weight(istep, nsamples, k, w_re, w_im)
      integer, intent(in) :: istep, nsamples, k
      real, intent(out) :: w_re, w_im
      real :: theta

      if (k == 0) then
         w_re = 1.0d0/real(nsamples)
         w_im = 0.0d0
         return
      end if
      theta = TWOPI*real(k)*real(istep)/real(nsamples)
      w_re = cos(theta)/real(nsamples)
      w_im = -sin(theta)/real(nsamples)

   end subroutine fourier_dft_weight

!-----------------------------------------------------------------------
! fourier_dft_add — Add one sample to the coefficients
!
!   c_re and c_im hold modes 0..kmax. One call is one time step.
!-----------------------------------------------------------------------
   pure subroutine fourier_dft_add(c_re, c_im, kmax, sample, istep, nsamples)
      integer, intent(in) :: kmax, istep, nsamples
      real, intent(inout) :: c_re(0:kmax), c_im(0:kmax)
      real, intent(in) :: sample
      integer :: kk
      real :: w_re, w_im

      do kk = 0, kmax
         call fourier_dft_weight(istep, nsamples, kk, w_re, w_im)
         c_re(kk) = c_re(kk) + sample*w_re
         c_im(kk) = c_im(kk) + sample*w_im
      end do

   end subroutine fourier_dft_add

!-----------------------------------------------------------------------
! fourier_eval_weight — Weight of mode k at time t
!
!   t is taken modulo the period. Mode 0 has weight 1.
!   Mode k adds 2 Re(c_k exp(2 pi i k t / T)).
!-----------------------------------------------------------------------
   pure subroutine fourier_eval_weight(k, t, period, w_re, w_im)
      integer, intent(in) :: k
      real, intent(in) :: t, period
      real, intent(out) :: w_re, w_im
      real :: theta

      if (k == 0) then
         w_re = 1.0d0
         w_im = 0.0d0
         return
      end if
      theta = TWOPI*real(k)*modulo(t, period)/period
      w_re = 2.0d0*cos(theta)
      w_im = -2.0d0*sin(theta)

   end subroutine fourier_eval_weight

!-----------------------------------------------------------------------
! fourier_eval_series — Evaluate the stored series at time t
!
!   c_re(0) is the mean. c_im(0) is not used.
!-----------------------------------------------------------------------
   pure function fourier_eval_series(c_re, c_im, kmax, t, period) result(value)
      integer, intent(in) :: kmax
      real, intent(in) :: c_re(0:kmax), c_im(0:kmax)
      real, intent(in) :: t, period
      real :: value
      integer :: kk
      real :: w_re, w_im

      value = 0.0d0
      do kk = 0, kmax
         call fourier_eval_weight(kk, t, period, w_re, w_im)
         value = value + c_re(kk)*w_re + c_im(kk)*w_im
      end do

   end function fourier_eval_series

!-----------------------------------------------------------------------
! fourier_select_k — Keep modes 0..K by the doubling rule
!
!   ek(0) is the mean. The fluctuation energy is ek(1)..ek(kacc).
!   kmax_fixed > 0 keeps every accumulated mode.
!   Otherwise K starts at 8. K doubles while the next count cuts
!   the missing fluctuation energy by 10 or more, down to 1e-10.
!   trace_* records each doubling that the run prints, accepted or not.
!-----------------------------------------------------------------------
   pure subroutine fourier_select_k(ek, kacc, kmax_fixed, k_sel, missing_out, &
                                    ntrace, trace_k, trace_next, trace_miss)
      integer, intent(in) :: kacc, kmax_fixed
      real, intent(in) :: ek(0:kacc)
      integer, intent(out) :: k_sel, ntrace
      real, intent(out) :: missing_out
      integer, intent(out) :: trace_k(:), trace_next(:)
      real, intent(out) :: trace_miss(:)
      real :: missing(0:kacc), etot
      integer :: kk, knext

      ntrace = 0
      etot = sum(ek(1:kacc))
      do kk = 0, kacc
         missing(kk) = 0.0d0
         if (etot > 0.0d0) missing(kk) = max(0.0d0, 1.0d0 - sum(ek(1:kk))/etot)
      end do

      k_sel = kacc
      if (kmax_fixed <= 0) then
         k_sel = min(K_FLOOR, kacc)
         do while (k_sel < kacc .and. missing(k_sel) >= MISSING_TARGET)
            knext = min(2*k_sel, kacc)
            ntrace = ntrace + 1
            trace_k(ntrace) = k_sel
            trace_next(ntrace) = knext
            trace_miss(ntrace) = missing(knext)
            if (missing(knext) > 0.1d0*missing(k_sel) .and. missing(knext) >= MISSING_TARGET) exit
            k_sel = knext
         end do
      end if
      missing_out = missing(k_sel)

   end subroutine fourier_select_k

end module nekstab_fourier_series
