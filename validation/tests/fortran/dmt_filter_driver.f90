! Reads commands from standard input and prints results of nekstab_dmt_filter in es24.16.
!   coef omega beta h           -> phi11 phi12 phi21 phi22 a1 a2 b1 b2
!   amp  omega beta t           -> cre(1:3) cim(1:3)
!   tone omega beta h n u0 a ph -> Re uhat, Im uhat, t after n steps of h, for the input
!                                  u = u0 + a*cos(omega*t + ph), filter started at equilibrium
program dmt_filter_driver
   use nekstab_dmt_filter, only: dmt_filter_step_coeffs, dmt_filter_amp_coeffs
   implicit none
   character(len=8) :: cmd
   character(len=256) :: line
   integer :: ios, n, k
   real(kind(0.0d0)) :: omega, beta, h, t, u0, a, ph, un, up, pp, qq, pn, qn
   real(kind(0.0d0)) :: phi(2, 2), ca(2), cb(2), cre(3), cim(3), re, im

   do
      read (*, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      read (line, *) cmd
      select case (trim(cmd))
      case ('coef')
         read (line, *) cmd, omega, beta, h
         call dmt_filter_step_coeffs(omega, beta, h, phi, ca, cb)
         write (6, '(8es24.16)') phi(1, 1), phi(1, 2), phi(2, 1), phi(2, 2), ca, cb
      case ('amp')
         read (line, *) cmd, omega, beta, t
         call dmt_filter_amp_coeffs(omega, beta, t, cre, cim)
         write (6, '(6es24.16)') cre, cim
      case ('tone')
         read (line, *) cmd, omega, beta, h, n, u0, a, ph
         t = 0.0d0
         up = u0 + a*cos(ph)
         if (omega > 0.0d0) then
            pp = up/omega**2
            qq = 0.0d0
         else
            pp = 0.0d0
            qq = up/beta
         end if
         call dmt_filter_step_coeffs(omega, beta, h, phi, ca, cb)
         do k = 1, n
            t = k*h
            un = u0 + a*cos(omega*t + ph)
            pn = phi(1, 1)*pp + phi(1, 2)*qq + ca(1)*up + cb(1)*un
            qn = phi(2, 1)*pp + phi(2, 2)*qq + ca(2)*up + cb(2)*un
            pp = pn
            qq = qn
            up = un
         end do
         call dmt_filter_amp_coeffs(omega, beta, t, cre, cim)
         re = cre(1)*qq + cre(2)*pp + cre(3)*up
         im = cim(1)*qq + cim(2)*pp + cim(3)*up
         write (6, '(3es24.16)') re, im, t
      case default
         write (6, *) 'unknown command ', trim(cmd)
         stop 1
      end select
   end do
end program dmt_filter_driver
