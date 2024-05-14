      !-----------------------------------------------------------------------c
      subroutine tdf            !Time-delayed Feedback
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv) :: do1, do2, do3
         real :: h1, semi, l2, linf, rate, tol
         real, save :: residu0, gain, porbit
         integer, save :: i, norbit, m
      
         logical, save :: init
         data init/.false./
      
         nt = nx1*ny1*nz1*nt
         if (.not. init) then
      
            if (nid == 0) write (6, *) '   Target frequency specified in uparam(5)=', uparam(5)
            porbit = 1.0d0/uparam(5) ! user input from .usr -> also used on the inflow BC
      
            call compute_cfl(ctarg, vx, vy, vz, 1.0d0) ! ctarg contains the sum ( ux_i / dx_i )
            dt = param(26)/ctarg ! dt given target CFL
            norbit = ceiling(porbit/dt) ! computing a safe value of norbit
            if (nid == 0) write (6, *) ' Computing norbit=', norbit
      
            dt = porbit/norbit   ! reducing dt to match forced orbit to machine accuracy
            param(12) = dt
      
            if (nid == 0) write (6, *) ' Adjusting timeStep dt=', dt
            call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
            if (nid == 0) write (6, *) ' veryfing current CFL and target=', ctarg, param(26)
            param(12) = -abs(param(12))
      
            gain = -0.04432d0*8*atan(1.0d0)/porbit ! Theoretical optimal feedback parameter see reference
      
            if (nid == 0) write (6, *) 'Allocating TDF orbit with nsteps:', norbit, norbit*dt
            allocate (uor(lv, norbit), vor(lv, norbit))
            if (if3d) then
               allocate (wor(lv, norbit))
            else
               allocate (wor(1, 1))
            end if
            call oprzero(uor(:, :), vor(:, :), wor(:, :))
      
            if (ifto .or. ldimt > 1) then
               allocate (tor(lv, norbit, ldimt))
               tor(:, :, :) = 0.0d0
            end if
      
            rate = 0.0d0; residu0 = 0.0d0
            open (unit=10, file='residu.dat')
      
            init = .true.
      
         else
      
            if (istep <= norbit) then !t<T->save solutions
      
               if (nid == 0) write (6, *) ' Storing initial solution in memory:', istep, '/', norbit
               call opcopy(uor(:, istep), vor(:, istep), wor(:, istep), vx, vy, vz)
               if (ifto) call copy(tor(1, istep, 1), t(:, :, :, :, 1), nt)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(tor(1, istep, m), t(:, :, :, :, m), nt)
               end do
               end if
      
            else !t>T->compute forcing !f(t)= - \Lambda * 2*pi*St * ( u(t) - u(t-T) )
      
               call opsub3(do1, do2, do3, vx, vy, vz, uor(:, 1), vor(:, 1), wor(:, 1)) !ub=v-vold
               if (istep > norbit + 1) call normvc(h1, semi, l2, linf, do1, do2, do3); rate = (l2 - residu0); residu0 = l2
      
               call opcmult(do1, do2, do3, gain) !f=fc*-chi
               call opadd2(fcx, fcy, fcz, do1, do2, do3) !FORCE HERE DO NOT COPY, ADD!
      
               do i = 1, norbit - 1     !discard the i=1 solution
      
                  uor(:, i) = uor(:, i + 1)
                  vor(:, i) = vor(:, i + 1)
                  if (if3d) wor(:, i) = wor(:, i + 1)
                  if (ifto) tor(:, i, 1) = tor(:, i + 1, 1)
                  if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) tor(:, i, m) = tor(:, i + 1, m)
                  end do
                  end if
      
               end do !store the last one
               call opcopy(uor(1, norbit), vor(1, norbit), wor(1, norbit), vx, vy, vz) !store solution
               if (ifto) call copy(tor(1, norbit, 1), t(:, :, :, :, 1), nt)
               if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(tor(1, norbit, m), t(:, :, :, :, m), nt)
               end do
               end if
      
               if (nid == 0) then
                  write (10, "(3E15.7)") time, l2, rate
                  write (6, "(' TDF residu=',1pE11.4,' rate of change= ',1pE11.4)") l2, rate
                  write (6, *) ' '
               end if
      
               tol = max(param(21), param(22))
               if (l2 > 0.0d0 .and. l2 < tol) then
                  if (nid == 0) write (6, *) ' Converged base flow to:', tol
                  ifbfcv = .true.
                  call bcast(ifbfcv, lsize)
                  param(63) = 1    ! Enforce 64-bit output
                  call bcast(param, 200*wdsize)
                  call outpost(vx, vy, vz, pr, t, 'BF_')
                  param(63) = 0    ! Enforce 32-bit output
                  call bcast(param, 200*wdsize)
                  call outpost_vort(vx, vy, vz, 'BFV')
               end if
      
            end if                  ! else
         end if                     ! not init
      
         return
      end subroutine TDF
      !-----------------------------------------------------------------------
      
      subroutine SFD
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         type(krylov_vector), save :: oldQ, oldV, qa, qb, qc
         type(krylov_vector) :: tempD, tempM
         real adt, bdt, cdt, cutoff, gain, res, h1, l2, semi, linf, oldRes, frq, sig, rate, dtol
         save oldRes, dtol
      
         frq = abs(uparam(04))*8*atan(1.0d0) ! St to omega
         sig = abs(uparam(05))
      
         if (uparam(5) > 0) then
      
            if (uparam(4) > 0) then ! Akervik
               cutoff = 0.5*frq
               gain = -2*sig
            else ! Casacuberta
               cutoff = 0.5*(sqrt(frq**2 + sig**2) - sig)
               gain = -0.5*(sqrt(frq**2 + sig**2) + sig)
            end if
            if (nid == 0) open (unit=10, file='residu.dat')
      
            if (istep == 0) then
               dtol = max(param(21), param(22)); oldRes = 0.0d0
               call k_zero(tempM)
               call k_zero(qa)
               call k_zero(qb)
               call k_zero(qc)
               call nopcopy(oldQ%vx, oldQ%vy, oldQ%vz, oldQ%pr, oldQ%t, vx, vy, vz, pr, t)
               call nopcopy(oldV%vx, oldV%vy, oldV%vz, oldV%pr, oldV%t, vx, vy, vz, pr, t)
            else
               call setab3(adt, bdt, cdt)
               call k_copy(qc, qb)
               call k_copy(qb, qa)
               call k_sub3(qa, oldV, oldQ)
               call k_copy(tempD, qa)
               call k_cmult(tempD, adt)
               call k_copy(tempM, tempD)
               call k_copy(tempD, qb)
               call k_cmult(tempD, bdt)
               call k_add2(tempM, tempD)
               call k_copy(tempD, qc)
               call k_cmult(tempD, cdt)
               call k_add2(tempM, tempD)
               call k_cmult(tempM, cutoff*dt)
               call k_add2(oldQ, tempM)
            end if
      
            call k_sub3(tempD, oldV, oldQ)
            call k_cmult(tempD, gain)
            call nopadd2(fcx, fcy, fcz, fcp, fct, tempD%vx, tempD%vy, tempD%vz, tempD%pr, tempD%t)
      
         elseif (nid == 0) then
      
            open (unit=10, file='residu.dat')
            write (6, *) ' SFD in continuation mode'
      
         end if
      
         if (istep >= 1) then
            call nopsub2(oldV%vx, oldV%vy, oldV%vz, oldV%pr, oldV%t, vx, vy, vz, pr, t)
            call normvc(h1, semi, l2, linf, oldV%vx, oldV%vy, oldV%vz)
            res = l2; rate = (res - oldRes)/dt; oldRes = res
            call nopcopy(oldV%vx, oldV%vy, oldV%vz, oldV%pr, oldV%t, vx, vy, vz, pr, t)
            if (nid == 0) then
               write (10, "(4E15.7)") time, res, rate, param(21)
               write (6, "(A,3E15.7)") '  SFD residual, rate:', res, rate
               if (uparam(4) > 0) then
                  write (6, *) ' Akervik  cutoff, gain:', cutoff, gain
               else
                  write (6, *) ' Casacub. cutoff, gain:', cutoff, gain
               end if
            end if
            if (ifdyntol .and. mod(istep, 20) == 0 .and. res > 0) call set_solv_tole(res/20.0)
      
            if (istep > 100 .and. res < dtol) then
               if (nid == 0) write (6, *) ' Converged base flow to:', res
               ifbfcv = .true.
               call bcast(ifbfcv, lsize)
               param(63) = 1
               call bcast(param, 200*wdsize)
               call outpost2(vx, vy, vz, pr, t, nof, 'BF_')
               param(63) = 0
               call bcast(param, 200*wdsize)
               call outpost_vort(vx, vy, vz, 'BFV')
            end if
      
            if (istep == nsteps) close (10)
      
         end if
      end subroutine SFD
      !-----------------------------------------------------------------------
      subroutine BoostConv
      !     boostconv core subroutine
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv) :: dvx, dvy, dvz
         real :: residu, h1, semi, linf, rate, tol
         real, save :: residu0
         logical, save :: init
         data init/.false./
      
         if (mod(istep, bst_skp) == 0) then
      
            if (.not. init) then
               residu = 0.0d0; rate = 0.0d0; residu0 = 0.0d0
               open (unit=10, file='residu.dat')
               init = .true.
            end if
      
            call opsub3(dvx, dvy, dvz, vx, vy, vz, vxlag(1, 1, 1, 1, 1), vylag(1, 1, 1, 1, 1), vzlag(1, 1, 1, 1, 1)) !dv=v-vold
            call normvc(h1, semi, residu, linf, dvx, dvy, dvz); rate = (residu - residu0)/dt; residu0 = residu
            call boostconv_core(dvx, dvy, dvz)
            call opadd3(vx, vy, vz, vxlag(1, 1, 1, 1, 1), vylag(1, 1, 1, 1, 1), vzlag(1, 1, 1, 1, 1), dvx, dvy, dvz) !v=vold+dv
      
            if (nid == 0) then
               write (10, "(3E15.7)") time, residu, rate
               write (6, "(' BoostConv residu=',1pE11.4,' delta= ',1pE11.4)") residu, rate
               write (6, *) ' '
            end if
      
            tol = max(param(21), param(22))
            if (residu < tol) then
               if (nid == 0) write (6, *) ' Converged base flow to:', tol
               ifbfcv = .true.
               call bcast(ifbfcv, lsize)
               param(63) = 1       ! Enforce 64-bit output
               call bcast(param, 200*wdsize)
               call outpost(vx, vy, vz, pr, t, 'BF_')
               param(63) = 0       ! Enforce 32-bit output
               call bcast(param, 200*wdsize)
            end if
      
         end if
      
         return
      end subroutine BoostConv
      !-----------------------------------------------------------------------
      subroutine boostconv_core(rbx, rby, rbz)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         integer, save :: rot
         integer :: j
      
         logical, save :: init
         data init/.false./
      
         real, allocatable, save, dimension(:) :: cc, ccb
         real, allocatable, save, dimension(:, :) :: dd
         real, allocatable, save, dimension(:, :) :: q_x, q_y, q_z
         real, allocatable, save, dimension(:, :) :: x_x, x_y, x_z, y_x, y_y, y_z
      
         real, dimension(lv) :: rbx, rby, rbz, dumx, dumy, dumz
      
         real :: glsc3
         nv = nx1*ny1*nz1*nelv
      
         if (.not. init) then
      
            allocate (cc(bst_snp), ccb(bst_snp), dd(bst_snp, bst_snp))
            allocate (q_x(lv, bst_snp), q_y(lv, bst_snp), q_z(lv, bst_snp))
            allocate (x_x(lv, bst_snp), x_y(lv, bst_snp), x_z(lv, bst_snp))
            allocate (y_x(lv, bst_snp), y_y(lv, bst_snp), y_z(lv, bst_snp))
      
            if (nid == 0) write (6, *) 'Allocating BoostConv variables with:', bst_snp
            if (nid == 0) write (6, *) '                     skipping every:', bst_skp
      
            call oprzero(x_x(:, :), x_y(:, :), x_z(:, :))
            call oprzero(y_x(:, :), y_y(:, :), y_z(:, :))
            call oprzero(q_x(:, :), q_y(:, :), q_z(:, :))
            call opcopy(y_x(:, 1), y_y(:, 1), y_z(:, 1), rbx, rby, rbz)
            call opcopy(x_x(:, 1), x_y(:, 1), x_z(:, 1), rbx, rby, rbz)
            dd(:, :) = 1.0d0; rot = 1; init = .true.
      
         else
      
            call opsub2(y_x(:, rot), y_y(:, rot), y_z(:, rot), rbx, rby, rbz)
            call opsub2(x_x(:, rot), x_y(:, rot), x_z(:, rot), y_x(:, rot), y_y(:, rot), y_z(:, rot))
            call qr_dec(dd, q_x, q_y, q_z, y_x, y_y, y_z)
      
            do j = 1, bst_snp
               cc(j) = glsc3(rbx, bm1, q_x(:, j), nv) + glsc3(rby, bm1, q_y(:, j), nv)
               if (if3d) cc(j) = cc(j) + glsc3(rbz, bm1, q_z(:, j), nv)
            end do
      
            call linear_system(ccb, cc, dd, bst_snp); rot = mod(rot, bst_snp) + 1
            call opcopy(y_x(:, rot), y_y(:, rot), y_z(:, rot), rbx, rby, rbz)
      
            do j = 1, bst_snp
               call opcopy(dumx, dumy, dumz, x_x(:, j), x_y(:, j), x_z(:, j))
               call opcmult(dumx, dumy, dumz, ccb(j))
               call opadd2(rbx, rby, rbz, dumx, dumy, dumz)
            end do
            call opcopy(x_x(:, rot), x_y(:, rot), x_z(:, rot), rbx, rby, rbz)
      
         end if
         return
      end subroutine boostconv_core
      !-----------------------------------------------------------------------
      subroutine qr_dec(rr, q_x, q_y, q_z, x_x, x_y, x_z)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer i, j
         real, dimension(lv, bst_snp) :: x_x, x_y, x_z, q_x, q_y, q_z !res subspaces
         real, dimension(lv) :: dum_x1, dum_y1, dum_z1, dum_x, dum_y, dum_z
         real, dimension(bst_snp, bst_snp) :: rr
         real norma, glsc3
      
         nv = nx1*ny1*nz1*nelv
         rr = 0.0d0; norma = 0.0d0
         call oprzero(Q_x(:, :), Q_y(:, :), Q_z(:, :))
      
         call opcopy(dum_x, dum_y, dum_z, x_x(:, 1), x_y(:, 1), x_z(:, 1))
      
         norma = glsc3(dum_x, bm1, dum_x, nv) + glsc3(dum_y, bm1, dum_y, nv)
         if (if3d) norma = norma + glsc3(dum_z, bm1, dum_z, nv)
         norma = sqrt(norma)
      
         call opcmult(dum_x, dum_y, dum_z, 1./norma)
         call opcopy(q_x(:, 1), q_y(:, 1), q_z(:, 1), dum_x, dum_y, dum_z)
         rr(1, 1) = norma
      
         do j = 2, bst_snp
            call opcopy(dum_x, dum_y, dum_z, x_x(:, j), x_y(:, j), x_z(:, j))
            do i = 1, j - 1
      
               rr(i, j) = glsc3(dum_x, bm1, q_x(:, i), nv) + glsc3(dum_y, bm1, q_y(:, i), nv)
               if (if3d) rr(i, j) = rr(i, j) + glsc3(dum_z, bm1, q_z(:, i), nv)
      
               call opcopy(dum_x1, dum_y1, dum_z1, q_x(:, i), q_y(:, i), q_z(:, i))
               call opcmult(dum_x1, dum_y1, dum_z1, rr(i, j))
               call opsub2(dum_x, dum_y, dum_z, dum_x1, dum_y1, dum_z1)
      
            end do
      
            norma = glsc3(dum_x, bm1, dum_x, nv) + glsc3(dum_y, bm1, dum_y, nv)
            if (if3d) norma = norma + glsc3(dum_z, bm1, dum_z, nv)
      
            if (norma < 1e-60) then
               norma = 1.0d0
               q_x(:, j) = 0.0d0; q_y(:, j) = 0.0d0
               if (if3d) q_z(:, j) = 0.0d0
            else
               call opcmult(dum_x, dum_y, dum_z, 1.0d0/sqrt(norma))
               call opcopy(q_x(:, j), q_y(:, j), q_z(:, j), dum_x, dum_y, dum_z)
            end if
      
            rr(j, j) = sqrt(norma)
      
         end do
         return
      end subroutine qr_dec
      !----------------------------------------------------------------------
      subroutine linear_system(outp, inp, m, size_m)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer :: size_m, j, k
         real :: m(size_m, size_m), inp(size_m), outp(size_m)
         outp = 0.0d0
         do j = size_m, 1, -1
      !outp(j) = inp(j) - sum(m(j,j+1:size_m)*outp(j+1:size_m))/m(j,j)
            outp(j) = inp(j)
            do k = j + 1, size_m
               outp(j) = outp(j) - m(j, k)*outp(k)
            end do
            outp(j) = outp(j)/m(j, j)
         end do
         return
      end subroutine linear_system
      !----------------------------------------------------------------------
