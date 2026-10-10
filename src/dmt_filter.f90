!-----------------------------------------------------------------------
! dmt_filter.f90 — Oscillator-pair band-pass filter of the DMT bank
!
! Purpose:
!   Coefficients of one filter of the DMT bank (dmt.f90, mode 6.4). The
!   routines have no dependency on Nek5000, so a stand-alone test can
!   compile them.
!
! Filter (u is the input field, one value per grid point):
!   p' = q,   q' = u - beta*q - omega**2*p
!   A damped oscillator of natural frequency omega and damping beta driven
!   by u. The band-pass output b = beta*q has gain 1 and phase 0 at omega
!   and gain 0 at omega = 0, the second-order-in-time form of the DMT
!   filter. For omega = 0 the same equations give the exponential mean:
!   b = beta*q follows u with the time constant 1/beta.
!
! Complex amplitude:
!   The signal b(t) = Re[Z e^{i omega t}] with slowly varying Z has the
!   derivative b' = beta*q' = beta*(u - beta*q - omega**2*p), which the
!   equations give without differentiating data. Then
!     Z = b - i*b'/omega,   uhat = Z e^{-i omega t},   u ~ Re[uhat e^{i omega t}].
!   The combination beta*(q + i*omega*p) has the same value in the band, but
!   it passes a constant mean flow with the gain i*beta/omega. The
!   derivative form removes that gain. For a mean of 1 and a tone of
!   amplitude 0.01, the error of the amplitude, in units of the tone amplitude,
!   was 5.0 with beta*(q + i*omega*p) and 1.9e-3 with Z, at beta/omega = 0.05
!   after 40 periods; the 1.9e-3 is the start-up transient, exp(-beta*t/2).
!
! Time stepping:
!   u is known at the steps only. It is taken linear between two steps and
!   the oscillator is integrated exactly, with the matrix exponential of
!   the augmented system (p, q, u, du), so a step of any length gives the
!   same accuracy and a variable time step needs no special case.
!
! Public interface:
!   dmt_filter_step_coeffs — z(n+1) = phi*z(n) + a*u(n) + b*u(n+1), z = (p, q)
!   dmt_filter_amp_coeffs  — uhat = c(1)*q + c(2)*p + c(3)*u, real and imaginary part
!
! Dependencies:
!   (none)
!-----------------------------------------------------------------------

module nekstab_dmt_filter
   implicit none
   private
   public :: dmt_filter_step_coeffs, dmt_filter_amp_coeffs

   integer, parameter :: dp = kind(0.0d0)

contains

!-----------------------------------------------------------------------
! expm4 — Matrix exponential of a 4 x 4 matrix
!
!   Scaling and squaring with a Taylor series: the matrix is divided by
!   2**s until its infinity norm is below 0.5, 16 terms of the series
!   then leave a truncation error below 1e-19, and e is squared s times.
!-----------------------------------------------------------------------
   subroutine expm4(m, e)
      real(dp), dimension(4, 4), intent(in) :: m
      real(dp), dimension(4, 4), intent(out) :: e

      real(dp), dimension(4, 4) :: a, term
      real(dp) :: nrm
      integer :: i, k, s

      nrm = maxval(sum(abs(m), dim=2))
      s = 0
      if (nrm > 0.5d0) s = ceiling(log(nrm/0.5d0)/log(2.0d0))
      a = m/2.0d0**s

      e = 0.0d0
      do i = 1, 4
         e(i, i) = 1.0d0
      end do
      term = e
      do k = 1, 16
         term = matmul(term, a)/real(k, dp)
         e = e + term
      end do
      do i = 1, s
         e = matmul(e, e)
      end do

   end subroutine expm4

!-----------------------------------------------------------------------
! dmt_filter_step_coeffs — Exact step of the oscillator for a linear input
!
!   z = (p, q). The input goes from u(n) to u(n+1) linearly over the step h:
!     z(n+1) = phi*z(n) + a*u(n) + b*u(n+1)
!   The augmented system (p, q, u, du) has u' = du/h and du' = 0, and its
!   exponential over h holds phi, the response to u(n) and the response to
!   du = u(n+1) - u(n) in its first two rows.
!
!   omega — natural angular frequency (rad per time unit), 0 for the mean
!   beta  — damping (1 per time unit), 0 < beta < omega, or beta > 0 at omega = 0
!   h     — step length
!-----------------------------------------------------------------------
   subroutine dmt_filter_step_coeffs(omega, beta, h, phi, a, b)
      real(dp), intent(in) :: omega, beta, h
      real(dp), dimension(2, 2), intent(out) :: phi
      real(dp), dimension(2), intent(out) :: a, b

      real(dp), dimension(4, 4) :: m, e
      real(dp), dimension(2) :: g0, g1

      m = 0.0d0
      m(1, 2) = 1.0d0
      m(2, 1) = -omega**2
      m(2, 2) = -beta
      m(2, 3) = 1.0d0
      m(3, 4) = 1.0d0/h
      call expm4(m*h, e)

      phi = e(1:2, 1:2)
      g0 = e(1:2, 3)
      g1 = e(1:2, 4)
      a = g0 - g1
      b = g1

   end subroutine dmt_filter_step_coeffs

!-----------------------------------------------------------------------
! dmt_filter_amp_coeffs — Coefficients of the complex amplitude
!
!   With s = sin(omega*t), c = cos(omega*t) and q' = u - beta*q - omega**2*p:
!     Re uhat = beta*(c*q - s*q'/omega)
!     Im uhat = -beta*(s*q + c*q'/omega)
!   written as  Re uhat = cre(1)*q + cre(2)*p + cre(3)*u  and the same for cim.
!   At omega = 0 the amplitude is the mean, beta*q, and has no imaginary part.
!
!   t — time of the state (p, q) and of u
!-----------------------------------------------------------------------
   subroutine dmt_filter_amp_coeffs(omega, beta, t, cre, cim)
      real(dp), intent(in) :: omega, beta, t
      real(dp), dimension(3), intent(out) :: cre, cim

      real(dp) :: s, c

      if (omega == 0.0d0) then
         cre = (/beta, 0.0d0, 0.0d0/)
         cim = 0.0d0
         return
      end if

      s = sin(omega*t)
      c = cos(omega*t)
      cre(1) = beta*(c + beta*s/omega)
      cre(2) = beta*omega*s
      cre(3) = -beta*s/omega
      cim(1) = beta*(beta*c/omega - s)
      cim(2) = beta*omega*c
      cim(3) = -beta*c/omega

   end subroutine dmt_filter_amp_coeffs

end module nekstab_dmt_filter
