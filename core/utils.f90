      !-----------------------------------------------------------------------
      !     function to check the velocity value in a node
      !     posiz   = position in the velocity vector of the good grid point
      !     procmin = processor that contains the good grid point
      subroutine pointcheck(posiz, procmin)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer m, n, procmin, posiz
         real chk(lx1*ly1*lz1*lelt), chkmin, glmin
         n = nx1*ny1*nz1*nelt
         procmin = 0
         if (nid == 0) write (6, *) 'Evaluating probe at:', xck, yck, zck
         do m = 1, n
            chk(m) = (xm1(m, 1, 1, 1) - xck)**2 + (ym1(m, 1, 1, 1) - yck)**2
            if (if3D) chk(m) = chk(m) + (zm1(m, 1, 1, 1) - zck)**2
         end do
         chkmin = glmin(chk, n)
         do m = 1, n
         if (chkmin == chk(m)) then
            procmin = 1
            posiz = m
            print *, 'Point found: ', m
         end if
         end do
         return
      end subroutine pointcheck
      !-----------------------------------------------------------------------
      subroutine quicksort2(n, arr, idx)
         implicit none
         integer, intent(in) :: n
         real, dimension(n), intent(inout) :: arr
         integer, dimension(n), intent(inout) :: idx
         integer, parameter :: m = 7, nstack = 50
         integer :: i, ir, j, jstack, k, l, istack(nstack)
         real :: a, temp
         integer :: b, temp_int
      
         jstack = 0
         l = 1
         ir = n
      
         do
         if (ir - l < m) then
         do j = l + 1, ir
            a = arr(j)
            b = idx(j)
            i = j - 1
            do while (i >= l .and. arr(i) > a)
               arr(i + 1) = arr(i)
               idx(i + 1) = idx(i)
               i = i - 1
            end do
            arr(i + 1) = a
            idx(i + 1) = b
         end do
         if (jstack == 0) return
         ir = istack(jstack)
         l = istack(jstack - 1)
         jstack = jstack - 2
         else
         k = (l + ir)/2
         call swap_real(arr(k), arr(l + 1), temp)
         call swap_int(idx(k), idx(l + 1), temp_int)
      
         if (arr(l) > arr(ir)) then
            call swap_real(arr(l), arr(ir), temp)
            call swap_int(idx(l), idx(ir), temp_int)
         end if
      
         if (arr(l + 1) > arr(ir)) then
            call swap_real(arr(l + 1), arr(ir), temp)
            call swap_int(idx(l + 1), idx(ir), temp_int)
         end if
      
         if (arr(l) > arr(l + 1)) then
            call swap_real(arr(l), arr(l + 1), temp)
            call swap_int(idx(l), idx(l + 1), temp_int)
         end if
      
         i = l + 1
         j = ir
         a = arr(l + 1)
         b = idx(l + 1)
      
         do
         do while (arr(i) < a)
            i = i + 1
         end do
         do while (arr(j) > a)
            j = j - 1
         end do
         if (i >= j) exit
         call swap_real(arr(i), arr(j), temp)
         call swap_int(idx(i), idx(j), temp_int)
         i = i + 1
         j = j - 1
         end do
      
         arr(l + 1) = arr(j)
         arr(j) = a
         idx(l + 1) = idx(j)
         idx(j) = b
         jstack = jstack + 2
         if (jstack > nstack) then
            print *, "NSTACK too small in quicksort2."
            return
         end if
      
         if (ir - i + 1 >= j - l) then
            istack(jstack) = ir
            istack(jstack - 1) = i
            ir = j - 1
         else
            istack(jstack) = j - 1
            istack(jstack - 1) = l
            l = i
         end if
         end if
         end do
      
      contains
      
         subroutine swap_real(x, y, temp)
            real, intent(inout) :: x, y, temp
            temp = x
            x = y
            y = temp
         end subroutine swap_real
      
         subroutine swap_int(x, y, temp)
            integer, intent(inout) :: x, y, temp
            temp = x
            x = y
            y = temp
         end subroutine swap_int
      
      end subroutine quicksort2
      !-----------------------------------------------------------------------
      ! subroutine quicksort2(n, arr, idx)
      !    implicit none
      !    integer n, m, nstack
      !    real arr(n)
      !    integer idx(n)
      !    parameter(m=7, nstack=50)
      !    integer i, ir, j, jstack, k, l, istack(nstack)
      !    real a, temp
      !    integer b, temp_int
      
      !    jstack = 0
      !    l = 1
      !    ir = n
      
      ! 1  if (ir - l .lt. m) then
      !       do 12 j = l + 1, ir
      !          a = arr(j)
      !          b = idx(j)
      !          do 11 i = j - 1, l, -1
      !             if (arr(i) .le. a) goto 2
      !             arr(i + 1) = arr(i)
      !             idx(i + 1) = idx(i)
      ! 11          continue
      !             i = l - 1
      ! 2           arr(i + 1) = a
      !             idx(i + 1) = b
      ! 12          continue
      
      !             if (jstack .eq. 0) return
      
      !             ir = istack(jstack)
      !             l = istack(jstack - 1)
      !             jstack = jstack - 2
      !             else
      !             k = (l + ir)/2
      !             temp = arr(k)
      !             arr(k) = arr(l + 1)
      !             arr(l + 1) = temp
      !             temp_int = idx(k)
      !             idx(k) = idx(l + 1)
      !             idx(l + 1) = temp_int
      
      !             if (arr(l) .gt. arr(ir)) then
      !                temp = arr(l)
      !                arr(l) = arr(ir)
      !                arr(ir) = temp
      
      !                temp_int = idx(l)
      !                idx(l) = idx(ir)
      !                idx(ir) = temp_int
      !             end if
      
      !             if (arr(l + 1) .gt. arr(ir)) then
      !                temp = arr(l + 1)
      !                arr(l + 1) = arr(ir)
      !                arr(ir) = temp
      
      !                temp_int = idx(l + 1)
      !                idx(l + 1) = idx(ir)
      !                idx(ir) = temp_int
      !             end if
      
      !             if (arr(l) .gt. arr(l + 1)) then
      !                temp = arr(l)
      !                arr(l) = arr(l + 1)
      !                arr(l + 1) = temp
      !                temp_int = idx(l)
      !                idx(l) = idx(l + 1)
      !                idx(l + 1) = temp_int
      !             end if
      
      !             i = l + 1
      !             j = ir
      !             a = arr(l + 1)
      !             b = idx(l + 1)
      
      ! 3           continue
      
      !             i = i + 1
      
      !             if (arr(i) .lt. a) goto 3
      
      ! 4           continue
      
      !             j = j - 1
      
      !             if (arr(j) .gt. a) goto 4
      !             if (j .lt. i) goto 5
      
      !             temp = arr(i)
      !             arr(i) = arr(j)
      !             arr(j) = temp
      
      !             temp_int = idx(i)
      !             idx(i) = idx(j)
      !             idx(j) = temp_int
      !             goto 3
      
      ! 5           arr(l + 1) = arr(j)
      !             arr(j) = a
      !             idx(l + 1) = idx(j)
      !             idx(j) = b
      !             jstack = jstack + 2
      !             if (jstack .gt. nstack) pause "..." !'NSTACK too small in quicksort2'
      
      !             if (ir - i + 1 .ge. j - 1) then
      !                istack(jstack) = ir
      !                istack(jstack - 1) = i
      !                ir = j - 1
      !             else
      !                istack(jstack) = j - 1
      !                istack(jstack - 1) = l
      !                l = i
      !             end if
      !             end if
      !             goto 1
      !             end subroutine quicksort2
      !-----------------------------------------------------------------------
      subroutine add_noise_scal(qin, fc1, fc2, fc3)
         implicit none
         include 'SIZE'
         include 'TSTEP'           ! TIME, DT
         include 'PARALLEL'        ! LGLEL
         include 'INPUT'           ! if3D
         include 'SOLN'            ! VX, VY, VZ, VMULT
         include 'GEOM'            ! XM1, YM1, ZM1
         real, intent(inout), dimension(lx1*ly1*lz1*lelt) :: qin
         real, intent(in) :: fc1, fc2, fc3
         real, dimension(lx1, ly1, lz1, lelt) :: q
         integer iel, ieg, il, jl, kl, nt
         real xl(ldim), mth_rand, fc(3), nmin, nmax, glmax, glmin
         fc(1) = fc1; fc(2) = fc2; fc(3) = fc3
         nt = nx1*ny1*nz1*nelt
         call copy(q(:, :, :, :), qin(:), nt)
         do iel = 1, nelv
         do kl = 1, nz1
         do jl = 1, ny1
         do il = 1, nx1
            ieg = lglel(iel)
            xl(1) = xm1(il, jl, kl, iel)
            xl(2) = ym1(il, jl, kl, iel)
            if (if3D) xl(ndim) = zm1(il, jl, kl, iel)
            q(il, jl, kl, iel) = q(il, jl, kl, iel) + mth_rand(il, jl, kl, ieg, xl, fc)
         end do
         end do
         end do
         end do
         call dssum(q, lx1, ly1, lz1)
         call col2(q, vmult, nt)
         call dsavg(q)
         call bcdirSC(q)
         call copy(qin(:), q(:, :, :, :), nt) ! RESHAPE ARRAY TO 1D
         nmin = glmin(qin(:), nt); nmax = glmax(qin(:), nt)
         if (nid == 0) write (6, *) 'noise scal min,max', nmin, nmax
         return
      end subroutine add_noise_scal
      !-----------------------------------------------------------------------
      subroutine op_add_noise(qx, qy, qz)
      !     input random number to fields
         implicit none
         include 'SIZE'            ! NX1, NY1, NZ1, NELV, NID
         include 'TSTEP'           ! TIME, DT
         include 'PARALLEL'        ! LGLEL
         include 'INPUT'           ! if3D
         include 'SOLN'            ! VX, VY, VZ, VMULT
         include 'GEOM'            ! XM1, YM1, ZM1
      
         real, dimension(lx1, ly1, lz1, lelv) :: qx, qy, qz
         integer iel, ieg, il, jl, kl, nv
         real xl(LDIM), mth_rand, fc(3), nmin, nmax, glmax, glmin
         nv = nx1*ny1*nz1*nelv
      
         do iel = 1, NELV
         do kl = 1, NZ1
         do jl = 1, NY1
         do il = 1, NX1
            ieg = LGLEL(iel)
            xl(1) = XM1(il, jl, kl, iel)
            xl(2) = YM1(il, jl, kl, iel)
            if (if3D) xl(NDIM) = ZM1(il, jl, kl, iel)
      
            fc(1) = 3.0e4; fc(2) = -1.5e3; fc(3) = 0.5e5
            qx(il, jl, kl, iel) = qx(il, jl, kl, iel) + mth_rand(il, jl, kl, ieg, xl, fc)
      
            fc(1) = 2.3e4; fc(2) = 2.3e3; fc(3) = -2.0e5
            qy(il, jl, kl, iel) = qy(il, jl, kl, iel) + mth_rand(il, jl, kl, ieg, xl, fc)
      
            if (if3D) then
               fc(1) = 2.e4; fc(2) = 1.e3; fc(3) = 1.e5
               qz(il, jl, kl, iel) = qz(il, jl, kl, iel) + mth_rand(il, jl, kl, ieg, xl, fc)
            end if
      
         end do
         end do
         end do
         end do
      
      !     face averaging
         call opdssum(qx(:, :, :, :), qy(:, :, :, :), qz(:, :, :, :))
         call opcolv(qx(:, :, :, :), qy(:, :, :, :), qz(:, :, :, :), VMULT)
      
         call dsavg(qx(:, :, :, :))
         call dsavg(qy(:, :, :, :))
         if (if3D) call dsavg(qz(:, :, :, :))
      
      !Note: v*mask removes points at wall/inflow
         call bcdirVC(qx(:, :, :, :), qy(:, :, :, :), qz(:, :, :, :), v1mask, v2mask, v3mask)
      
         nmin = glmin(qx, nv); nmax = glmax(qx, nv)
         if (nid == 0) write (6, *) 'noise vx min,max', nmin, nmax
      
         nmin = glmin(qy, nv); nmax = glmax(qy, nv)
         if (nid == 0) write (6, *) 'noise vy min,max', nmin, nmax
      
         if (if3D) then
            nmin = glmin(qz, nv); nmax = glmax(qz, nv)
            if (nid == 0) write (6, *) 'noise vz min,max', nmin, nmax
         end if
         return
      end subroutine op_add_noise
      !-----------------------------------------------------------------------
      subroutine add_symmetric_seed(qx, qy, qz, qp) !generate symmetric seed to fields
         implicit none
         include "SIZE"
         include "TOTAL"
         real, dimension(lx1, ly1, lz1, lelv) :: qx, qy, qz, qp
         integer iel, ieg, il, jl, kl, ntot
         real xlx, yly, zlz, alpha, x, y, z
         real glsc3, amp
      
         ntot = NX1*NY1*NZ1*NELV
         xlx = xmx - xmn
         yly = ymx - ymn
         zlz = zmx - zmn
      
         alpha = 2*pi/zlz
      
      !     --> Create the initial velocity perturbation.
      
         do iel = 1, NELV
         do kl = 1, NZ1
         do jl = 1, NY1
         do il = 1, NX1
      
            ieg = LGLEL(iel)
            x = XM1(il, jl, kl, iel)
            y = YM1(il, jl, kl, iel)
            if (if3D) z = ZM1(il, jl, kl, iel)
      
      !     -> Construct the perturbation. ! Note: Spanwise invariant.
            qx(il, jl, kl, iel) = cos(alpha*z)*sin(2.*pi*y)
            qz(il, jl, kl, iel) = -(2.*pi)/(alpha)*cos(alpha*z)*cos(2.*pi*y)
            qp(il, jl, kl, iel) = cos(alpha*z)*cos(2.*pi*y)
      
         end do
         end do
         end do
         end do
      
         amp = glsc3(qx, bm1, qx, ntot) + glsc3(qy, bm1, qy, ntot)
         if (if3D) amp = amp + glsc3(qz, bm1, qz, ntot)
         amp = 1e-6/(0.50d0*amp)
         call opcmult(qx, qy, qz, amp)
         call cmult(qp, amp, ntot)
      
         return
      end subroutine add_symmetric_seed
      !-----------------------------------------------------------------------
      real function mth_rand(ix, iy, iz, ieg, xl, fc) !generate random number
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3D
         integer ix, iy, iz, ieg
         real xl(LDIM), fc(3)
         mth_rand = fc(1)*(ieg + xl(1)*sin(xl(2))) + fc(2)*ix*iy + fc(3)*ix
         if (if3D) mth_rand = fc(1)*(ieg + xl(NDIM)*sin(mth_rand)) + fc(2)*iz*ix + fc(3)*iz
         mth_rand = cos(1.e3*sin(1.e3*sin(mth_rand)))
         return
      end function mth_rand
      !-----------------------------------------------------------------------
      subroutine outpost_vort(ux, uy, uz, name)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, intent(in) :: ux(lv), uy(lv), uz(lv)
         character(len=3) :: name
         real wo1(lv), wo2(lv), wo3(lv), vort(lv, 3)
         logical ifto_sav, ifpo_sav
      
         if (ifvor) then
      
            call comp_vort3(vort, wo1, wo2, ux, uy, uz)
      
            ifto_sav = ifto
            ifpo_sav = ifpo
            ifto = .false.
            ifpo = .false.
            call outpost(vort(1, 1), vort(1, 2), vort(1, 3), pr, t, name)
            ifto = ifto_sav
            ifpo = ifpo_sav
         end if
      
         return
      end subroutine outpost_vort
      !-----------------------------------------------------------------------
      subroutine norm_grad(vx_, vy_, vz_, pr_, t_, norma)
      
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, intent(in), dimension(lv) :: vx_, vy_, vz_
         real, intent(in), dimension(lp) :: pr_
         real, intent(in), dimension(lt, ldimt) :: t_
         real, intent(out) :: norma
      
         real, dimension(lv) :: dudx, dudy, dudz
         real, dimension(lv) :: dvdx, dvdy, dvdz
         real, dimension(lv) :: dwdx, dwdy, dwdz
      
         real :: glsc3
         nv = nx1*ny1*nz1*nelv
      
      ! gradient computation
         call gradm1(dudx, dudy, dudz, vx_, nelv)
         call gradm1(dvdx, dvdy, dvdz, vy_, nelv)
         if (if3D) call gradm1(dwdx, dwdy, dwdz, vz_, nelv)
      
      ! call dsavg(dudx); call dsavg(dudy); call dsavg(dudz)
      ! call dsavg(dvdx); call dsavg(dvdy); call dsavg(dvdz)
      ! call dsavg(dwdx); call dsavg(dwdy); call dsavg(dwdz)
      
         norma = 0.0d0
      
         norma = norma + glsc3(dudx, bm1s, dudx, nv) + glsc3(dudy, bm1s, dudy, nv)
         norma = norma + glsc3(dvdx, bm1s, dvdx, nv) + glsc3(dvdy, bm1s, dvdy, nv)
      
         if (if3D) then
            norma = norma + glsc3(dudz, bm1s, dudz, nv)
            norma = norma + glsc3(dvdz, bm1s, dvdz, nv)
            norma = norma + glsc3(dwdx, bm1s, dwdx, nv)
            norma = norma + glsc3(dwdy, bm1s, dwdy, nv)
            norma = norma + glsc3(dwdz, bm1s, dwdz, nv)
         end if
      end subroutine norm_grad
      
      !-----------------------------------------------------------------------
      subroutine nekStab_outpost
      !     nekStab custom outpost routine
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real vort(lv, 3), wo1(lv), wo2(lv)
         common/ugrad/vort, wo1, wo2
      
         logical ifto_sav, ifpo_sav, ifpsco_sav(ldimt1)
      
         if ((istep == 0) .or. ifoutfld .and. (ifvor .or. ifvox)) then
      
            ifto_sav = ifto; ifpo_sav = ifpo; ifpsco_sav = ifpsco
            ifto = .false.; ifpo = .false.; ifpsco(:) = .false.
      
      !---  > Compute and oupost vorticity.
            if (ifvor) then
               call oprzero(vort(1, 1), vort(1, 2), vort(1, 3))
               call comp_vort3(vort, wo1, wo2, vx, vy, vz)
               call outpost(vort(1, 1), vort(1, 2), vort(1, 3), pr, t, 'vor')
            end if
      
      !---  > Compute and outpost vortex fields.
            if (ifvox .and. (ifaxis .eqv. .false.)) then
               call vortex_core(vort(:, 1), 'lambda2')
               call vortex_core(vort(:, 2), 'q')
               call vortex_core(vort(:, 3), 'omega')
               call outpost(vort(:, 1), vort(:, 2), vort(:, 3), pr, t, 'vox')
      
               if (.not. if3d) then ! 2D case -> no lambda2, no q
                  ifvo = .false.; ifto = .true. ! just outposting omega field to temperature... v's and p ignored
                  call outpost(vx, vy, vz, pr, vort(:, 3), 'vox')
                  ifvo = .true.
               end if
      
            end if
      
            ifto = ifto_sav; ifpo = ifpo_sav; ifpsco(:) = ifpsco_sav(:)
      
         end if
      
      !---  > Outpost initial condition.
         if (istep == 0) call outpost2(vx, vy, vz, pr, t, nof, '   ')
      
         return
      end subroutine nekStab_outpost
      !-----------------------------------------------------------------------
      subroutine nekStab_comment ! verbose output
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, save :: eetime0, eetime1, eetime2, deltatime
         data eetime0, eetime1, eetime2/0.0d0, 0.0d0, 0.0d0/
         real telapsed, tpernondt, tmiss, dnekclock, ttime
         integer ttime_stp
      
      !     if extrapolation is not OIFS -> ifchar = false
      !     if OIFS -> ifchar = .true. and CFL 2-5
      !     cases can have CFL > 1 in initial time steps
         if (courno > 10) then
         if (nio == 0) then
            write (6, *)
            write (6, *) '    CFL > 10 stopping code'
            write (6, *)
         end if
         call nek_end
         end if
         if (nio /= 0) return
      
         eetime2 = dnekclock()
         if (eetime0 == 0.0 .and. istep == 1) then
            eetime0 = eetime2
            deltatime = time
         else
            eetime1 = eetime2
         end if
      
         if (istep > 0 .and. lastep == 0 .and. iftran) then
            ttime_stp = eetime2 - eetime1
            ttime = eetime2 - eetime0
      
            if (istep == 1) then
               ttime_stp = 0.0d0
               ttime = 0.0d0
            end if
      
            if (mod(istep, 10) == 0) then
               telapsed = ttime/3600.0d0
               tpernondt = ttime/(time - deltatime)
               tmiss = (param(10) - time)*tpernondt/3600.0d0
      
               print *, ''
               write (6, "('      Mean time per timestep: ',F8.4,'  dev:',I8,'ms')")
     $   ttime/istep, nint(((ttime/istep) - ttime_stp)*1000)
               write (6, "('      Estimated remaining time: ',I8,' hr ',I2,' min')")
     $   int(tmiss), nint((tmiss - int(tmiss))*60)
               if (tpernondt > 60.) then
                  write (6, "('      Time per nondimensional time: ',F8.2,' sec')") tpernondt
               else
                  write (6, "('      Time per nondimensional time: ',F8.2,' min ')") tpernondt/60.0
               end if
               write (6, "('      Current local time: ',F10.4,'  File:',I8)")
     $   time - deltatime, int((time - deltatime)/param(14)) + 1
               print *, ''
            end if
         end if
      
      end subroutine nekStab_comment
      !-----------------------------------------------------------------------
      subroutine nekStab_printNEKParams ! print initialization for sanity check
         implicit none
         include 'SIZE'
         include 'TOTAL'
         if (nid == 0) then
            write (6, *) 'P01=', param(1), 'density'
            write (6, *) 'P02=', param(2), 'viscosity (1/Re)'
            write (6, *) 'P07=', param(7), 'rhoCp'
            write (6, *) 'P08=', param(8), 'conductivity (1/(Re*Pr))'
            write (6, *) 'P10=', param(10), 'stop at endTime'
            write (6, *) 'P10=', param(11), 'stop at numSteps'
            write (6, *) 'P14=', param(14), 'io step'
            write (6, *) 'P15=', param(15), 'io time'
            write (6, *) 'P21=', param(21), 'pressure sol tol'
            write (6, *) 'P22=', param(22), 'velocity sol tol'
            write (6, *) 'P26=', param(26), 'target CFL number'
            write (6, *) 'P27=', param(27), 'order in time'
            write (6, *) 'P28=', param(28), 'use same torder for mesh solver'
            write (6, *) 'P31=', param(31), 'numberOfPerturbations'
            write (6, *) 'P41=', param(41), '1 for multiplicative SEMG'
            write (6, *) 'P42=', param(42), 'lin solv for the pres equation 0:GMRES,1:CG'
            write (6, *) 'P43=', param(43), '0:additive multilevel scheme 1:orig 2lvl sch'
            write (6, *) 'P44=', param(44), '0=E-based addit Schwarz PnPn-2;1=A-based'
            write (6, *) 'P93=', param(93), 'num vectors for projection'
            write (6, *) 'P94 =', param(94), 'num projection for helmholz solves'
            write (6, *) 'P95=', param(95), 'projection for pressure solver on/off'
            write (6, *) 'P101=', param(101), 'no additional modes'
            write (6, *) 'P103=', param(103), 'filter weight'
            write (6, *)
            write (6, *) 'uparam01=', uparam(1)
            write (6, *) 'uparam02=', uparam(02)
            write (6, *) 'uparam03=', uparam(03)
            write (6, *) 'uparam04=', uparam(04)
            write (6, *) 'uparam05=', uparam(05)
            write (6, *) 'uparam06=', uparam(06)
            write (6, *) 'uparam07=', uparam(07)
            write (6, *) 'uparam08=', uparam(08)
            write (6, *) 'uparam09=', uparam(09)
            write (6, *) 'uparam10=', uparam(10)
            write (6, *)
            write (6, *) 'x min,max,tot=', xmn, xmx, xmx - xmn
            write (6, *) 'y min,max,tot=', ymn, ymx, ymx - xmn
            write (6, *) 'z min,max,tot=', zmn, zmx, zmx - zmn
            write (6, *)
         end if
      end subroutine nekStab_printNEKParams
      !-----------------------------------------------------------------------
      subroutine nekStab_energy(px, py, pz, pt, fname, skip)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, dimension(lv), intent(in) :: px, py, pz
         real, dimension(lv, ldimt), intent(in) :: pt
         integer, intent(in) :: skip
         character(len=16), intent(in) :: fname
         real glsc3, uek, vek, wek, eek, pot
         save eek
         logical, save :: initialized
         data initialized/.false./
         nv = nx1*ny1*nz1*nelv
         nt = nx1*ny1*nz1*nelt
         eek = 0.50d0/volvm1
         uek = 0.0d0; vek = 0.0d0; wek = 0.0d0; pot = 0.0d0
      
         if (.not. initialized) then
            if (nid == 0) open (730, file=fname, action='write', status='replace')
            initialized = .true.
         end if
      
         if (mod(istep, skip) == 0) then
            uek = glsc3(px, bm1s, px, nv)
            vek = glsc3(py, bm1s, py, nv)
            if (if3d) wek = glsc3(pz, bm1s, pz, nv)
            if (ifheat) pot = glsc3(pt(:, 1), bm1s, ym1, nt)
            if (nid == 0) write (730, "(6E15.7)") time, uek*eek, vek*eek, wek*eek, (uek + vek + wek)*eek, pot*eek
         end if
      
         return
      end subroutine nekStab_energy
      !-----------------------------------------------------------------------
      subroutine nekStab_enstrophy(px, py, pz, pt, fname, skip)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, dimension(lv), intent(in) :: px, py, pz
         real, dimension(lv, ldimt), intent(in) :: pt
      !real, dimension(lx1,ly1,lz1,lelt,ldimt), intent(in) :: pt
         integer, intent(in) :: skip
         character(len=19), intent(in) :: fname
         real vort(lv, 3), wo1(lv), wo2(lv)
         common/ugrad/vort, wo1, wo2
         real glsc3, uek, vek, wek, eek
         save eek
         logical, save :: initialized
         data initialized/.false./
         nv = nx1*ny1*nz1*nelv
         eek = 0.50d0/volvm1
         uek = 0.0d0; vek = 0.0d0; wek = 0.0d0
      
         if (.not. initialized) then
            if (nid == 0) open (736, file=fname, action='write', status='replace')
            initialized = .true.
         end if
      
         if (mod(istep, skip) == 0) then
            call comp_vort3(vort, wo1, wo2, px, py, pz)
            uek = glsc3(vort(1, 1), bm1s, vort(1, 1), nv)
            vek = glsc3(vort(1, 2), bm1s, vort(1, 2), nv)
            if (if3d) wek = glsc3(vort(1, 3), bm1s, vort(1, 3), nv)
            if (nid == 0) write (736, "(5E15.7)") time, uek*eek, vek*eek, wek*eek, (uek + vek + wek)*eek
      
         end if
      
         return
      end subroutine nekStab_enstrophy
      !-----------------------------------------------------------------------
      subroutine nekStab_torque(fname, skip)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         integer, intent(in) :: skip
         character(len=20), intent(in) :: fname
      
         logical, save :: initialized
         data initialized/.false./
      
         integer, save :: bIDs(1), iobj_wall(1)
      
         real, save :: x0(3), scale
         data x0/3*0/
      
         integer :: nij, i, iobj, memtot, mem, ieg, ifc, ie
         integer, parameter :: lr = lx1*ly1*lz1
      
         real glmin, glmax, x1min, x2min, x3min, x1max, x2max, x3max, w1(0:maxobj)
         real flow_rate, base_flow, domain_length, xsec, scale_vf(3), sij, pm1, xm0, ym0, zm0
      
         real ur, us, ut, vr, vs, vt, wr, ws, wt
         common/scruz/ur(lr), us(lr), ut(lr), vr(lr), vs(lr), vt(lr), wr(lr), ws(lr), wt(lr)
         common/cvflow_r/flow_rate, base_flow, domain_length, xsec, scale_vf
         common/scrns/sij(lx1*ly1*lz1*6*lelv)
         common/scrcg/pm1(lx1, ly1, lz1, lelv)
         common/scrsf/xm0(lx1, ly1, lz1, lelt), ym0(lx1, ly1, lz1, lelt), zm0(lx1, ly1, lz1, lelt)
         nv = lx1*ly1*lz1*nelv
      
         if (.not. initialized) then
            if (nid == 0) write (6, *) 'Initializing torque routine...'
            if (nid == 0) open (737, file=fname, action='write', status='replace')
            bIDs(1) = 1
            call create_obj(iobj_wall(1), bIDs, 1)
            scale = 2
            initialized = .true.
         end if
      
         if (mod(istep, skip) == 0) then
      
            call mappr(pm1, pr, xm0, ym0) ! map pressure onto Mesh 1
            if (param(55) /= 0) then
               dpdx_mean = -scale_vf(1)
               dpdy_mean = -scale_vf(2)
               dpdz_mean = -scale_vf(3)
            end if
            call add2s2(pm1, xm1, dpdx_mean, nv) ! Doesn't work if object is cut by
            call add2s2(pm1, ym1, dpdy_mean, nv) ! periodicboundary.  In this case,
            call add2s2(pm1, zm1, dpdz_mean, nv) ! set ._mean=0 and compensate in
            nij = 3
            if (if3d .or. ifaxis) nij = 6
            call comp_sij(sij, nij, vx, vy, vz, ur, us, ut, vr, vs, vt, wr, ws, wt)
            if (istep < 1) call cfill(vdiff, param(2), nv)
            call cadd2(xm0, xm1, -x0(1), nv)
            call cadd2(ym0, ym1, -x0(2), nv)
            call cadd2(zm0, zm1, -x0(3), nv)
            x1min = glmin(xm0(1, 1, 1, 1), nv)
            x2min = glmin(ym0(1, 1, 1, 1), nv)
            x3min = glmin(zm0(1, 1, 1, 1), nv)
            x1max = glmax(xm0(1, 1, 1, 1), nv)
            x2max = glmax(ym0(1, 1, 1, 1), nv)
            x3max = glmax(zm0(1, 1, 1, 1), nv)
            do i = 0, maxobj
               dragpx(i) = 0       ! BIG CODE  :}
               dragvx(i) = 0
               dragx(i) = 0
               dragpy(i) = 0
               dragvy(i) = 0
               dragy(i) = 0
               dragpz(i) = 0
               dragvz(i) = 0
               dragz(i) = 0
               torqpx(i) = 0
               torqvx(i) = 0
               torqx(i) = 0
               torqpy(i) = 0
               torqvy(i) = 0
               torqy(i) = 0
               torqpz(i) = 0
               torqvz(i) = 0
               torqz(i) = 0
            end do
            ifield = 1
            do iobj = 1, nobj
               memtot = nmember(iobj)
               do mem = 1, memtot
                  ieg = object(iobj, mem, 1)
                  ifc = object(iobj, mem, 2)
                  if (gllnid(ieg) == nid) then ! this processor has a contribution
                     ie = gllel(ieg)
                     call drgtrq(dgtq, xm0, ym0, zm0, sij, pm1, vdiff, ifc, ie)
                     call cmult(dgtq, scale, 12)
                     dragpx(iobj) = dragpx(iobj) + dgtq(1, 1) ! pressure
                     dragpy(iobj) = dragpy(iobj) + dgtq(2, 1)
                     dragpz(iobj) = dragpz(iobj) + dgtq(3, 1)
                     dragvx(iobj) = dragvx(iobj) + dgtq(1, 2) ! viscous
                     dragvy(iobj) = dragvy(iobj) + dgtq(2, 2)
                     dragvz(iobj) = dragvz(iobj) + dgtq(3, 2)
                     torqpx(iobj) = torqpx(iobj) + dgtq(1, 3) ! pressure
                     torqpy(iobj) = torqpy(iobj) + dgtq(2, 3)
                     torqpz(iobj) = torqpz(iobj) + dgtq(3, 3)
                     torqvx(iobj) = torqvx(iobj) + dgtq(1, 4) ! viscous
                     torqvy(iobj) = torqvy(iobj) + dgtq(2, 4)
                     torqvz(iobj) = torqvz(iobj) + dgtq(3, 4)
                  end if
               end do
            end do
            call gop(dragpx, w1, '+  ', maxobj + 1)
            call gop(dragpy, w1, '+  ', maxobj + 1)
            call gop(dragpz, w1, '+  ', maxobj + 1)
            call gop(dragvx, w1, '+  ', maxobj + 1)
            call gop(dragvy, w1, '+  ', maxobj + 1)
            call gop(dragvz, w1, '+  ', maxobj + 1)
            call gop(torqpx, w1, '+  ', maxobj + 1)
            call gop(torqpy, w1, '+  ', maxobj + 1)
            call gop(torqpz, w1, '+  ', maxobj + 1)
            call gop(torqvx, w1, '+  ', maxobj + 1)
            call gop(torqvy, w1, '+  ', maxobj + 1)
            call gop(torqvz, w1, '+  ', maxobj + 1)
            do i = 1, nobj
               dragx(i) = dragpx(i) + dragvx(i)
               dragy(i) = dragpy(i) + dragvy(i)
               dragz(i) = dragpz(i) + dragvz(i)
               torqx(i) = torqpx(i) + torqvx(i)
               torqy(i) = torqpy(i) + torqvy(i)
               torqz(i) = torqpz(i) + torqvz(i)
               dragpx(0) = dragpx(0) + dragpx(i)
               dragvx(0) = dragvx(0) + dragvx(i)
               dragx(0) = dragx(0) + dragx(i)
               dragpy(0) = dragpy(0) + dragpy(i)
               dragvy(0) = dragvy(0) + dragvy(i)
               dragy(0) = dragy(0) + dragy(i)
               dragpz(0) = dragpz(0) + dragpz(i)
               dragvz(0) = dragvz(0) + dragvz(i)
               dragz(0) = dragz(0) + dragz(i)
               torqpx(0) = torqpx(0) + torqpx(i)
               torqvx(0) = torqvx(0) + torqvx(i)
               torqx(0) = torqx(0) + torqx(i)
               torqpy(0) = torqpy(0) + torqpy(i)
               torqvy(0) = torqvy(0) + torqvy(i)
               torqy(0) = torqy(0) + torqy(i)
               torqpz(0) = torqpz(0) + torqpz(i)
               torqvz(0) = torqvz(0) + torqvz(i)
               torqz(0) = torqz(0) + torqz(i)
            end do
            do i = 1, nobj
            if (nio == 0) then
            if (if3d .or. ifaxis) then
               write (737, "(i8,19E15.7)") istep, time,
     $   dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), dragz(i), dragpz(i), dragvz(i),
     $   torqx(i), torqpx(i), torqvx(i), torqy(i), torqpy(i), torqvy(i), torqz(i), torqpz(i), torqvz(i)
            else
               write (737, "(i8,10E15.7)") istep, time,
     $   dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), torqz(i), torqpz(i), torqvz(i)
            end if
            end if
            end do
         end if
         return
      end subroutine nekStab_torque
      !-----------------------------------------------------------------------
      subroutine nekStab_define_obj
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer iel, ifc
      
         do iel = 1, nelt
         do ifc = 1, 2*ndim
            if (cbc(ifc, iel, 1) == 'W  ') boundaryID(ifc, iel) = 1
         end do
         end do
      
         return
      end subroutine nekStab_define_obj
      !-----------------------------------------------------------------------
      subroutine zero_crossing(v_mean_init)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, intent(in) :: v_mean_init
         real, save :: T_delayed(lv, 3), do1(lv), do2(lv), do3(lv), v_mean
         integer, parameter :: plor = 2
         real h1, l2, semi, linf
         save l2
         real :: glsum, dtime, vdot, vddot
         real, save :: velp(plor), v_sum, time0
         real, save :: t_cross, t_cross_old
         real, save :: v_cross, v_cross_old
         real, save :: p_now, p_sum, p_mean, p_old
         integer, save :: probe_nel, probe_nid, t_cross_count
         integer :: i
      
         if (istep == 0) then
            if (nid == 0) write (6, *) 'Initializing zero-crossing routine...'
            probe_nel = 0; probe_nid = 0; vdot = 0.0d0; vddot = 0.0d0
            velp(:) = 0.0d0; v_sum = 0.0d0; v_mean = v_mean_init
            p_now = 0.0d0; p_sum = 0.0d0; p_mean = 0.0d0
            t_cross = 0.0d0; t_cross_old = 0.0d0; t_cross_count = 0
            time0 = time; p_old = 0.0d0; l2 = 0.0d0
            call pointcheck(probe_nel, probe_nid) !alter xck, yck, zck in usrchck
            if (nid == 0) open (unit=17, file='zc_period.dat')
            if (nid == 0) open (unit=19, file='zc_poincare.dat')
            call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
         end if
      
         velp(plor) = 0.0d0
         if (probe_nid == 1) velp(plor) = vy(probe_nel, 1, 1, 1)
         velp(plor) = glsum(velp(plor), 1) !broadcast
         dtime = time - time0
      
         if (istep > 1) then
            v_sum = v_sum + velp(plor)*dt !(v_mean*(dtime-dt)+v_now*dt)/dtime
            v_mean = v_mean_init + v_sum/dtime
      !     if(nid.eq.0)write(6,*)' probe: v_mean, v_now= ',v_mean,velp(plor)
      !     if(nid.eq.0)write(6,*)'         v_sum, dtime= ',v_sum,dtime
         end if
         if (velp(plor - 1) <= v_mean .and. velp(plor) >= v_mean) then !period found
      
            p_old = p_now          !save old value
            t_cross_old = t_cross  !save old value
            t_cross = dtime        !update new value
            p_now = t_cross - t_cross_old !compute period
      
            call opsub3(do1, do2, do3, vx, vy, vz, T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3)) !ub=v-vold
            call normvc(h1, semi, l2, linf, do1, do2, do3)
            call opcopy(T_delayed(:, 1), T_delayed(:, 2), T_delayed(:, 3), vx, vy, vz)
            if (nid == 0) write (6, *) ' Zero-crossing T=', p_now, abs(p_now - p_old), l2
            v_cross_old = v_cross; v_cross = velp(plor)
            if (nid == 0) write (17, "(5E15.7)") time, p_now, abs(p_now - p_old), v_mean, l2
      
         end if
      !     https://en.wikipedia.org/wiki/Finite_difference_coefficient
         if (istep > 3) then
            vdot = ((11./6.)*velp(plor) - 3*velp(plor - 1) + 1.5*velp(plor - 2) - (1./3.)*velp(plor - 3))*dt**(-1)
            vddot = (2*velp(plor) - 5*velp(plor - 1) + 4.0*velp(plor - 2) - 1.0*velp(plor - 3))*dt**(-2)
         else
            vdot = 0.0d0; vddot = 0.0d0
         end if
         if (nid == 0) write (19, "(4E15.7)") time, velp(plor), vdot, vddot
      
         do i = 1, plor - 1
            velp(i) = velp(i + 1)
         end do
      
         return
      end subroutine zero_crossing
      !-----------------------------------------------------------------------
