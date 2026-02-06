      !-----------------------------------------------------------------------
      subroutine vortex_core(l2, vortex)
         include 'SIZE'
         include 'TOTAL'
         real l2(lx1, ly1, lz1, 1)
         character*(*) vortex
      
         if (vortex == "lambda2") then
            call lambda2(l2)
         elseif (vortex == "q") then
            call compute_q(l2)
         elseif (vortex == "delta") then
            call compute_delta(l2)
         elseif (vortex == "swirling") then
            call compute_swirling(l2)
         elseif (vortex == "omega") then
            call compute_omega_jc(l2)
         elseif (vortex == "symmetric") then
            call compute_symmetricVec(l2)
         elseif (vortex == "assymetric") then
            call compute_assymetricVec(l2)
         else
            if (nid == 0) write (6, *) "ABORT:unknown vortex type:", vortex
            call exitt
         end if
      
         return
      end subroutine vortex_core
      !---------------------------------------------------------
      subroutine compute_omega_jc(omega)
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         integer, parameter :: n = lx1*ly1*lz1*lelv
         integer, parameter :: nxyz = lx1*ly1*lz1
         real, parameter :: eps = 1.0d-5
      
      !     --> Element-wise velocity pseudo-gradient.
         real gije(lx1*ly1*lz1, ldim, ldim)
      
      !     --> Point-wise symmetric and anti-symmetric parts.
         real ss(ldim, ldim), oo(ldim, ldim)
         real omega(lx1, ly1, lz1, lelv)
         real norm_a, norm_b
      
      !     --> Miscellaneous.
         integer ie, l, i, j
      
      !     --> Loop through the elements.
         do ie = 1, nelv
      
      !     --> Compute the velocity pseudo-gradient.
            call comp_gije(gije, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
      
      !     --> Point-wise computations.
            do l = 1, nxyz
      
               do j = 1, ldim
                  do i = 1, ldim
      !     --> Compute the symmetric and antisymmetric component.
                     ss(i, j) = 0.50d0*(gije(l, i, j) + gije(l, j, i))
                     oo(i, j) = 0.50d0*(gije(l, i, j) - gije(l, j, i))
                  end do
               end do
      
      !     --> Compute the Frobenius norm
               norm_a = (norm2(ss)**2)**2
               norm_b = (norm2(oo)**2)**2
      
      !     --> Compute omega.
               omega(l, 1, 1, ie) = norm_b/(norm_a + norm_b + eps)
      
            end do
         end do
         call filter_s0(omega, 0.5, 1, 'vortx') !filtering is necessary here!
         return
      end subroutine compute_omega_jc
      !---------------------------------------------------------
      subroutine compute_omega(l2)
      !     Generate \Omega criterion vortex of
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer, parameter :: lxyz = lx1*ly1*lz1
         real l2(lx1, ly1, lz1, 1)
         real mygi(lxyz, ldim, ldim)
         integer n, ie, l, nxyz
         real a, b
         common/mygrad/mygi
         nxyz = lx1*ly1*lz1
         n = nxyz*nelv
         do ie = 1, nelv ! Compute velocity gradient tensor
            call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
            do l = 1, nxyz
               call compute_symmetric(A, l)
               call compute_antisymmetric(B, l)
               l2(l, 1, 1, ie) = (B**2)/(B**2 + A**2 + 0.0020d0*glmax_qc)
            end do
         end do
         call filter_s0(l2, 0.5, 1, 'vortx')
         return
      end subroutine compute_omega
      !---------------------------------------------------------
      subroutine compute_symmetricVec(l2)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer, parameter :: lxyz = lx1*ly1*lz1
         real l2(lx1, ly1, lz1, 1)
         real mygi(lxyz, ldim, ldim)
         integer n, ie, l, nxyz
         common/mygrad/mygi
         nxyz = lx1*ly1*lz1
         n = nxyz*nelv
         do ie = 1, nelv ! Compute velocity gradient tensor
            call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
            do l = 1, nxyz
               call compute_symmetric(l2(l, 1, 1, ie), l)
            end do
         end do
         call filter_s0(l2, 0.5, 1, 'vortx')
         return
      end subroutine compute_symmetricVec
      !---------------------------------------------------------
      subroutine compute_assymetricVec(l2)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer, parameter :: lxyz = lx1*ly1*lz1
         real l2(lx1, ly1, lz1, 1)
         real mygi(lxyz, ldim, ldim)
         integer n, ie, l, nxyz
         common/mygrad/mygi
         nxyz = lx1*ly1*lz1
         n = nxyz*nelv
         do ie = 1, nelv ! Compute velocity gradient tensor
            call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
            do l = 1, nxyz
               call compute_antisymmetric(l2(l, 1, 1, ie), l)
            end do
         end do
         call filter_s0(l2, 0.5, 1, 'vortx')
         return
      end subroutine compute_assymetricVec
      !---------------------------------------------------------
      subroutine compute_q(l2)
      !
      !     Generate Q criterion vortex of Hunt, Wray & Moin, CTR-S88 1988
      !     positive second invariant of velocity gradient tensor
      !
         include 'SIZE'
         include 'TOTAL'
      
         parameter(lxyz=lx1*ly1*lz1)
         real l2(lx1, ly1, lz1, 1)
         real mygi(lxyz, ldim, ldim)
         real Q1
         common/mygrad/mygi
      
         nxyz = lx1*ly1*lz1
         n = nxyz*nelv
      
         do ie = 1, nelv
            call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
            do l = 1, nxyz
               call compute_secondInv(Q1, l)
               l2(l, 1, 1, ie) = Q1
            end do
         end do
         call filter_s0(l2, 0.5, 1, 'vortx')
      
         return
      end subroutine compute_q
      !---------------------------------------------------------
      subroutine compute_delta(l2)
      !
      !     Generate  Discriminant (DELTA) criterion vortex of Chong, Perry & Cantwell, Phys. Fluids 1990
      !     complex eigenvalues of velocity gradient tensor
      !
         include 'SIZE'
         include 'TOTAL'
      
         parameter(lxyz=lx1*ly1*lz1)
         real l2(lx1, ly1, lz1, 1)
         real mygi(lxyz, ldim, ldim)
         real P1, Q1, R1, R
         common/mygrad/mygi
      
         nxyz = lx1*ly1*lz1
         n = nxyz*nelv
      
         do ie = 1, nelv
      !     Compute velocity gradient tensor
            call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
            do l = 1, nxyz
               call compute_firstInv(P1, l)
               P1 = -P1 !negative sign
               call compute_secondInv(Q1, l)
               call compute_thirdInv(R1, l)
               R1 = -R1 !negative sign
               R = R1 + (P1**3)*2/27 - P1*Q1/3
               l2(l, 1, 1, ie) = (R/2)**2 + (Q1/3)**3
            end do
         end do
         call filter_s0(l2, 0.5, 1, 'vortx')
         return
      end subroutine compute_delta
      !---------------------------------------------------------
      subroutine compute_swirling(l2)
      !
      !     Generate Swirling Strength criterion vortex of
      !     Zhou, Adrian, Balachandar and Kendall, JFM 1999
      !
      !     imaginary part of complex eigenvalues of velocity gradient tensor
      !     presented as lambda_ci^2
      !
      !     for the 3x3 case
      !     |d11, d12, d13|  |du/dx,du/dy,du/dz|
      !     D = [d_ij] = |d21, d22, d23|= |dv/dx,dv/dy,dv/dz|
      !     |d31, d32, d33|  |dw/dx,dw/dy,dw/dz|
      !
      !     |lambda_r,  0,       0      |
      !     [d_ij]= [vr,vcr,vci]|  0,   lambda_cr, lambda_ci|[vr,vcr,vci]^-1
      !     |  0,  -lambda_ci, lambda_cr|
      !
      !     λ^3+P*λ^2+Qλ+R=0
      !
      !     where P = -I_D   = -tr(D)
      !     Q = II_D  =  0.5[P*P - tr(DD)]
      !     R = -III_D = 1/3.[-P+3QP-tr(DDD)]
      !     for the 2x2 case
      !
      !     λ^2+Qλ+R=0
      !
      !     where Q = II_D  =  0.5[P*P - tr(DD)]
      !     R = -III_D = 1/3.[-P+3QP-tr(DDD)]
      !
      !
         include 'SIZE'
         include 'TOTAL'
         parameter(lxyz=lx1*ly1*lz1)
         real l2(lx1, ly1, lz1, 1)
         real gije(lxyz, ldim, ldim), mygi(lxyz, ldim, ldim)
         real P1, Q1, R1, R, lambdaCi
         common/mygrad/mygi
      
         nxyz = lx1*ly1*lz1
         n = nxyz*nelv
         if (if3d) then ! 3D CASE
            do ie = 1, nelv
      !     Compute velocity gradient tensor
               call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
               do l = 1, nxyz
                  call compute_firstInv(P1, l)
                  P1 = -P1 !negative sign
                  call compute_secondInv(Q1, l)
                  call compute_thirdInv(R1, l)
                  R1 = -R1 !negative sign
      !     Q=Q1-(P1**2)/3
                  R = R1 + (P1**3)*2/27 - P1*Q1/3
      !     Delta = (R/2)**2+(Q/3)**3
                  call cubicLambdaCi(P1, Q1, R1, lambdaCi)
      
                  l2(l, 1, 1, ie) = lambdaCi
               end do
            end do
         elseif (ifaxis) then ! AXISYMMETRIC CASE
            if (nid == 0) write (6, *)
     $   'ABORT:no compute_swirling axisymmetric support for now'
            call exitt
         else ! 2D CASE
            do ie = 1, nelv
      !     Compute velocity gradient tensor
               call comp_gije(mygi, vx(1, 1, 1, ie), vy(1, 1, 1, ie), vz(1, 1, 1, ie), ie)
               do l = 1, nxyz
                  call compute_firstInv(P1, l)
                  P1 = -P1 !negative sign
      !     call compute_secondInv(Q1,l)
      !     Q1=0.
                  call compute_thirdInv(R1, l)
                  R1 = -R1 !negative sign
      !     Q=Q1-(P1**2)/3
      !     R=R1+(P1**3)*2/27-P1*Q1/3
      !     Delta = (R/2)**2+(Q/3)**3
                  call quadLambdaCi(P1, R1, lambdaCi)
      
                  l2(l, 1, 1, ie) = lambdaCi
               end do
            end do
      
         end if
      
         call filter_s0(l2, 0.5, 1, 'vortx')
      
         do ie = 1, nelv
            do l = 1, nxyz
               l2(l, 1, 1, ie) = l2(l, 1, 1, ie)*l2(l, 1, 1, ie)
            end do
         end do
      
         return
      end subroutine compute_swirling
      !---------------------------------------------------------
      subroutine compute_antisymmetric(B, l) !B=0.5(\Nabla u - (\Nabla u)^T)
      !     |11(dudx),12(dudy),13(dudz)|
      !     |21(dvdx),22(dvdy),23(dvdz)|
      !     |31(dwdx),32(dwdy),33(dwdz)|
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer, parameter :: lxyz = lx1*ly1*lz1
         real mygi(lxyz, ldim, ldim), B
         integer l
         common/mygrad/mygi
         B = 0.0d0
         if (if3d) then ! 3D CASE
            B = ((mygi(l, 1, 2) + mygi(l, 2, 1))**2 + (mygi(l, 1, 3) + mygi(l, 3, 1))**2 + (mygi(l, 2, 3) + mygi(l, 3, 2))**2)/4
         else ! 2D CASE
            B = ((mygi(l, 1, 2) + mygi(l, 2, 1))**2)/4
         end if
         return
      end subroutine compute_antisymmetric
      !---------------------------------------------------------
      subroutine compute_symmetric(A, l) !A=0.5(\Nabla u + (\Nabla u)^T)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer, parameter :: lxyz = lx1*ly1*lz1
         real mygi(lxyz, ldim, ldim), A, B
         integer l
         common/mygrad/mygi
         A = 0.0d0; B = 0.0d0
         if (if3d) then ! 3D CASE
            B = ((mygi(l, 1, 2) + mygi(l, 2, 1))**2 + (mygi(l, 1, 3) + mygi(l, 3, 1))**2 + (mygi(l, 2, 3) + mygi(l, 3, 2))**2)/4
            A = B + (mygi(l, 1, 1)**2 + mygi(l, 2, 2)**2 + mygi(l, 3, 3)**2)/2
         else ! 2D CASE
            B = ((mygi(l, 1, 2) + mygi(l, 2, 1))**2)/4
            A = B + (mygi(l, 1, 1)**2 + mygi(l, 2, 2)**2)/2
         end if
         return
      end subroutine compute_symmetric
      !---------------------------------------------------------
      subroutine compute_firstInv(a, l)
      
      !     for the 3x3 case
      !     |d11, d12, d13|
      !     D = [d_ij] = |d21, d22, d23|
      !     |d31, d32, d33|
      !     compute_firstInv returns tr(d_ij)=(d11+d22+d33)
      !
         include 'SIZE'
         include 'TOTAL'
         parameter(lxyz=lx1*ly1*lz1)
         real mygi(lxyz, ldim, ldim), a
         integer l
         common/mygrad/mygi
         if (if3d) then ! 3D CASE
            a = (mygi(l, 1, 1) + mygi(l, 2, 2) + mygi(l, 3, 3))
         elseif (ifaxis) then ! AXISYMMETRIC CASE
            if (nid == 0) write (6, *) 'ABORT: compute_firstInv axisymmetric support for now'
            call exitt
         else ! 2D CASE
            a = (mygi(l, 1, 1) + mygi(l, 2, 2))
         end if
         return
      end subroutine compute_firstInv
      !-------------------------------------------------------------------
      subroutine compute_secondInv(a, l)
      
      !     for the 3x3 case
      !     |d11, d12, d13|
      !     D = [d_ij] = |d21, d22, d23|
      !     |d31, d32, d33|
      !     compute_SecondInv returns
      !     0.5*(tr(d_ij)^2+tr(d_ij*d_ij))=-(d22*d33-d23*d32)-(d11*d22-d12*d21)-(d33*d11-d13*d31)
      
         include 'SIZE'
         include 'TOTAL'
         parameter(lxyz=lx1*ly1*lz1)
         real mygi(lxyz, ldim, ldim)
         real a
         integer l
         common/mygrad/mygi
         if (if3d) then ! 3D CASE
            a = (mygi(l, 2, 2)*mygi(l, 3, 3) - mygi(l, 2, 3)*mygi(l, 3, 2))
     $   +(mygi(l, 1, 1)*mygi(l, 2, 2) - mygi(l, 1, 2)*mygi(l, 2, 1))
     $   +(mygi(l, 3, 3)*mygi(l, 1, 1) - mygi(l, 1, 3)*mygi(l, 3, 1))
         elseif (ifaxis) then ! AXISYMMETRIC CASE
            if (nid == 0) write (6, *) 'ABORT: compute_secondInv axisymmetric support for now'
            call exitt
         else ! 2D CASE
            a = (mygi(l, 1, 1) + mygi(l, 2, 2))*(mygi(l, 1, 1) + mygi(l, 2, 2))
     $   -2*mygi(l, 1, 2)*mygi(l, 2, 1) - mygi(l, 1, 1)*mygi(l, 1, 1) - mygi(l, 2, 2)*mygi(l, 2, 2)
         end if
      
         return
      end subroutine compute_secondInv
      !-------------------------------------------------------------------
      subroutine compute_thirdInv(a, l)
      
      !     for the 3x3 case
      !     |d11, d12, d13|
      !     D = [d_ij] = |d21, d22, d23|
      !     |d31, d32, d33|
      
      !     compute_thirdInv returns det(D)
      !     =-d11*(d23*d32-d22*d33)-d12*(d21*d33-d31*d23)-d13*(d31*d22-d21*d32)
      
         include 'SIZE'
         include 'TOTAL'
         parameter(lxyz=lx1*ly1*lz1)
         real mygi(lxyz, ldim, ldim)
         real a
         integer l
         common/mygrad/mygi
         if (if3d) then ! 3D CASE
            a = -mygi(l, 1, 1)*(mygi(l, 2, 3)*mygi(l, 3, 2)
     $   -mygi(l, 2, 2)*mygi(l, 3, 3))
     $   -mygi(l, 1, 2)*(mygi(l, 2, 1)*mygi(l, 3, 3)
     $   -mygi(l, 3, 1)*mygi(l, 2, 3))
     $   -mygi(l, 1, 3)*(mygi(l, 3, 1)*mygi(l, 2, 2)
     $   -mygi(l, 2, 1)*mygi(l, 3, 2))
         elseif (ifaxis) then ! AXISYMMETRIC CASE
            if (nid == 0) write (6, *) 'ABORT: compute_thirdInv axisymmetric support for now'
            call exitt
         else ! 2D CASE
            a = mygi(l, 1, 1)*mygi(l, 2, 2) - mygi(l, 1, 2)*mygi(l, 2, 1)
      
         end if
         return
      end subroutine compute_thirdInv
      !-------------------------------------------------------------------
      subroutine cubicLambdaCi(b, c, d, lci)
         include 'SIZE'
         include 'TOTAL'
         parameter(lxyz=lx1*ly1*lz1)
         real mygi(lxyz, ldim, ldim)
         real b, c, d
         real f, g, h, r, s, t2, u
         real lci
         complex ci
         complex x1
      
         ci = sqrt(cmplx(-1.))
      
         a = 1
         f = c/a - 1/3.*(b/a)**2.
         g = ((2.*(b**3.)/(a**3.)) - (9.*b*c/(a**2.)) + (27.*d/a))/27.
         h = (g/2.)**2.+(f/3.)**3.
      
         if (h <= 0.) then
            lci = 0.
         else
            r = -(g/2.) + sqrt(h)
            if (r <= 0.) then
               s = sign(abs(r)**(1.0/3.0), r)
            else
               s = (r)**(1./3.)
            end if
            t2 = -(g/2.) - sqrt(h)
            if (t2 <= 0.) then
               u = sign(abs(t2)**(1.0/3.0), t2)
            else
               u = ((t2)**(1/3.))
            end if
            x1 = -(s + u)/2.+(b/3./a) + ci*(s - u)*sqrt(3.)/2.
            lci = aimag(x1)
         end if
      
      end subroutine cubicLambdaCi
      !-------------------------------------------------------------------
      subroutine quadLambdaCi(b, c, lci)
         include 'SIZE'
         include 'TOTAL'
         parameter(lxyz=lx1*ly1*lz1)
         real mygi(lxyz, ldim, ldim)
         real b, c
         real f
         real lci
         complex ci
         complex x1
      
         ci = sqrt(cmplx(-1.))
      
         a = 1
         d = b**2.-4.*a*c
      
         if (d >= 0.) then
            lci = 0.
         else
            f = sqrt(abs(d)) !sign(abs(d)**(1.0/2.0), d)
            x1 = -b/2./a + f*ci/2./a
            lci = aimag(x1)
         end if
      
      end subroutine quadLambdaCi
      !-------------------------------------------------------------------
      subroutine nekStab_avg(ifstatis)
         use krylov_subspace
         include 'SIZE'
         include 'TOTAL'
         include 'AVG'
      
         logical, intent(in) :: ifstatis
         real do1(lv), do2(lv), do3(lv)
      
         real, save :: x0(3)
         data x0/0.0, 0.0, 0.0/
      
         integer, save :: icalld
         data icalld/0/
      
         logical ifverbose
      
         if (ax1 /= lx1 .or. ay1 /= ly1 .or. az1 /= lz1) then
            if (nid == 0) write (6, *) 'ABORT: wrong size of ax1,ay1,az1 in avg_all(), check SIZE!'
            call exitt
         end if
         if (ax2 /= lx2 .or. ay2 /= ay2 .or. az2 /= lz2) then
            if (nid == 0) write (6, *) 'ABORT: wrong size of ax2,ay2,az2 in avg_all(), check SIZE!'
            call exitt
         end if
      
         ntot = lx1*ly1*lz1*nelv
         nto2 = lx2*ly2*lz2*nelv
      
      !     initialization
         if (icalld == 0) then
            icalld = icalld + 1
            atime = 0.; timel = time
            call oprzero(uavg, vavg, wavg)
            call rzero(pavg, nto2)
            call oprzero(urms, vrms, wrms)
            call rzero(prms, nto2)
            call oprzero(vwms, wums, uvms)
         end if
      
         dtime = time - timel
         atime = atime + dtime
      
      !     dump freq
         iastep = param(68)
         if (iastep == 0) iastep = param(15) ! same as iostep
         if (iastep == 0) iastep = 500
      
         ifverbose = .false.
         if (istep <= 10) ifverbose = .true.
         if (mod(istep, iastep) == 0) ifverbose = .true.
      
         if (atime /= 0. .and. dtime /= 0.) then
            if (nio == 0) write (6, *) 'Computing statistics ...'
            beta = dtime/atime
            alpha = 1.-beta
      
      !     compute averages E(X) !avg(k) = alpha*avg(k) + beta*f(k)
            call avg1(uavg, vx, alpha, beta, ntot, 'um  ', ifverbose)
            call avg1(vavg, vy, alpha, beta, ntot, 'vm  ', ifverbose)
            call avg1(wavg, vz, alpha, beta, ntot, 'wm  ', ifverbose)
            call avg1(pavg, pr, alpha, beta, nto2, 'prm ', ifverbose)
            call avg1(tavg, t(1, 1, 1, 1, 1), alpha, beta, ntot, 'tm ', ifverbose)
      
            if (ifstatis) then !compute fluctuations
      !     compute averages E(X^2)
               call avg2(urms, vx, alpha, beta, ntot, 'ums ', ifverbose)
               call avg2(vrms, vy, alpha, beta, ntot, 'vms ', ifverbose)
               call avg2(wrms, vz, alpha, beta, ntot, 'wms ', ifverbose)
               call avg2(prms, pr, alpha, beta, nto2, 'prms', ifverbose)
               call avg2(trms, t(1, 1, 1, 1, 1), alpha, beta, ntot, 'trms', ifverbose)
      
      !     compute averages E(X*Y)
      
            end if
      
      
         end if
      
         if ((mod(istep, iastep) == 0 .and. istep > 1) .or. lastep == 1) then
      
            time_temp = time
            time = atime ! Output the duration of this avg
            dtmp = param(63)
            param(63) = 1 ! Enforce 64-bit output
      
      !     mean fluctuation fields
            ifto = .false.; ifpo = .false.
            call opsub3(do1, do2, do3, uavg, vavg, wavg, ubase, vbase, wbase)
            call outpost(do1, do2, do3, pr, t, 'avt')
      
            call outpost2(uavg, vavg, wavg, pavg, tavg, ldimt, 'avg')
            call outpost2(urms, vrms, wrms, prms, trms, ldimt, 'rms')
            call outpost(uvms, vwms, wums, prms, trms, 'rm2')
      
      
            param(63) = dtmp
            atime = 0.
            time = time_temp ! Restore clock
      
         end if
         timel = time
      
         return
      end subroutine nekStab_avg
      !----------------------------------------------------------------------
      
      subroutine stability_energy_budget
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Arrays to store the instability mode real and imaginary parts.
         real, dimension(lv) :: vx_dRe, vy_dRe, vz_dRe, t_dRe
         real, dimension(lv) :: vx_dIm, vy_dIm, vz_dIm, t_dIm
         real, dimension(lp) :: pr_dRe, pr_dIm
      
      !     ----- Energy budget terms.
         real, dimension(lv, 10) :: energy_budget
         real, dimension(10) :: integrals
      
      !     ----- Miscellaneous.
         real :: alpha, beta, glsc2
         integer :: i, k, mode, n
         character(len=80) :: filename
         character(len=6) :: mode_str
         character(len=3) :: mode_str2
      
         n = nx1*ny1*nz1*nelv
      
      !     #####
      !     #####
      !     #####     PREPROCESSING
      !     #####
      !     #####
      
      !     --> Load the base flow.
         write (filename, '(a, a, a)') 'BF_', trim(SESSION), '0.f00001'
         call load_fld(filename)
         call nopcopy(ubase, vbase, wbase, pr, tbase, vx, vy, vz, pr, t)
      
         do mode = 1, maxmodes
            energy_budget(:, :) = 0.0d0
            integrals(:) = 0.0d0
      
      ! Format the mode number with leading zeros
            write (mode_str, '(i5.5)') mode
      
      !      --> Load the real part of the mode.
            write (filename, '(a, a, a, a)') 'dRe', trim(SESSION), '0.f', trim(mode_str)
            call load_fld(filename)
            call nopcopy(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, vx, vy, vz, pr, t)
      
      !     --> Load the imaginary part of the mode.
            write (filename, '(a, a, a, a)') 'dIm', trim(SESSION), '0.f', trim(mode_str)
            call load_fld(filename)
            call nopcopy(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, vx, vy, vz, pr, t)
      
      !     --> Normalize eigenmode to unit-norm.
      !         norm() returns ||q||, nopcmult(q, c) does q = c*q.
      !         Old bug: alpha = sqrt(alpha**2 + beta**2)
      !         gave ||q_new|| = alpha*||q|| = alpha^2 (wrong).
      !         Fix: invert so nopcmult divides by ||q_total||.
            call norm(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, alpha)
            call norm(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, beta)
      
            alpha = 1.0d0 / sqrt(alpha**2 + beta**2)
      
            call nopcmult(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, alpha)
            call nopcmult(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, alpha)
      
      !     #####
      !     #####
      !     #####     COMPUTE THE ENERGY BUDGET
      !     #####
      !     #####
      
      !     --> Compute the production terms.
            call compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, 1,
     $   energy_budget(:, 1), energy_budget(:, 2), energy_budget(:, 3))
            call compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, 2,
     $   energy_budget(:, 4), energy_budget(:, 5), energy_budget(:, 6))
            call compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, 3,
     $   energy_budget(:, 7), energy_budget(:, 8), energy_budget(:, 9))
      
      !     --> Compute the dissipation term.
            call compute_dissipation(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, energy_budget(:, 10))
      
      !     --> Compute the integrals and the sum.
            do i = 1, 10
               integrals(i) = glsc2(bm1, energy_budget(:, i), n)
            end do
      
            if (if3d) then
               k = 9
            else
               k = 6
            end if
      
            do i = 1, k, 3
               call opcopy(vx, vy, vz, energy_budget(:, i), energy_budget(:, i + 1), energy_budget(:, i + 2))
      
               write (mode_str2, "('K',I2.2)") mode
               call outpost(vx, vy, vz, pr, t, mode_str2)
            end do
      
            if (nid == 0) then
               write (filename, '(A,A,A,A)') 'PKE_dRe', trim(SESSION), '0.f', trim(mode_str)
               open (101, file=filename, form='formatted')
               do i = 1, 10
                  write (101, '(1E15.7)') integrals(i)
                  write (*, *) 'Integral ', i, ' = ', integrals(i)
               end do
               write (101, '(1E15.7)') sum(integrals(1:9))
               write (101, '(1E15.7)') sum(integrals(1:9)) - integrals(10)
               close (101)
            end if ! nid == 0
      
         end do
      
         return
      end subroutine stability_energy_budget

!----------------------------------------------------------------------

      subroutine stability_energy_budget_floquet
!
!     Orbit-averaged perturbation kinetic energy (PKE) budget for
!     Floquet stability analysis of time-periodic base flows.
!
!     MATHEMATICAL FORMULATION
!     ========================
!
!     Starting from the Reynolds-Orr energy equation for the
!     perturbation kinetic energy e = 0.5 * u'_i * u'_i :
!
!        de/dt = -u'_i u'_j (dU_i/dx_j) + (1/Re) u'_i lap(u'_i)
!                 \_________P__________/   \________D__________/
!                    production                dissipation
!
!     Integrating over the domain V gives the global budget:
!
!        dE/dt = P - D       where  E = int_V e dV
!
!     For a Floquet mode with real multiplier mu, the perturbation
!     takes the form:
!
!        u'(x,t) = exp(sigma_r * t) * u_hat(x,t)
!
!     where u_hat is T-periodic and sigma_r = ln|mu| / T is the
!     Floquet exponent (growth rate). Substituting into dE/dt and
!     orbit-averaging over one period T:
!
!        2 * sigma_r = (P_bar - D_bar) / E_bar
!
!     where bars denote orbit averages:  X_bar = (1/T) int_0^T X dt.
!     This is an approximation that becomes exact when E_hat(t) is
!     constant over the orbit; for typical converged Floquet modes
!     the error is O(var(E_hat)/mean(E_hat)), which is negligible.
!
!     GROWTH CORRECTION
!     =================
!
!     The time-stepper evolves the raw (growing) perturbation u'(t),
!     not the periodic part u_hat(t). The instantaneous budget terms
!     thus carry a factor exp(2*sigma_r*t) from the growth:
!
!        P_raw(t) = exp(2*sigma_r*t) * P_hat(t)
!
!     We correct this by weighting each timestep contribution with
!     exp(-2*sigma_r*t), so the accumulated average reflects u_hat:
!
!        P_bar = (1/T) sum_k  P_raw(t_k) * exp(-2*sigma_r*t_k) * dt
!
!     The same correction applies to D and E.
!
!     FACTOR-OF-1/2 CONVENTION
!     ========================
!
!     The helper routines compute_production and compute_dissipation
!     include an extra factor of 1/2 relative to the standard
!     Reynolds-Orr definitions:
!
!        compute_production:  P_ij = -0.5 * u'_i u'_j dU_i/dx_j
!                             (Reynolds-Orr has factor -1, not -0.5)
!
!        compute_dissipation: D = 0.5 * (mu/rho) * u'_i lap(u'_i)
!                             (Reynolds-Orr has factor 1, not 0.5)
!
!     This convention was chosen for consistency with the complex-mode
!     case (Re^2+Im^2 formulation). The factor cancels in the
!     verification formula because all three quantities carry it:
!
!        sigma_check = (0.5*P - 0.5*D) / (2 * E)
!                    = (P - D) / (4 * E)
!
!     and E = 0.5*|u'|^2 already has the physical 1/2, so:
!
!        sigma_check = (P - D) / (4 * 0.5 * |u'|^2_bar)
!                    = (P - D) / (2 * |u'|^2_bar)
!                    = sigma_r       (correct)
!
!     The individual spatial fields and integrals output to disk
!     carry this 0.5 factor. To obtain physical Reynolds-Orr values,
!     multiply production and dissipation fields by 2.
!
!     QUADRATURE
!     ==========
!
!     The orbit integral uses a right-endpoint rectangle rule:
!     fields and weight are both evaluated AFTER nek_advance
!     (i.e., at t_{k+1}). This is first-order accurate with
!     error O(dt/T) per orbit, which is negligible for typical
!     Nek5000 Floquet runs (hundreds of steps per period).
!
!     VERIFICATION
!     ============
!
!     sigma_r from the eigenvalue file must match
!     sigma_check = (P_bar - D_bar) / (2 * E_bar) to within
!     discretization error (~1% for typical resolutions).
!
!     SCOPE AND LIMITATIONS
!     =====================
!
!     - Real Floquet multipliers only (Mode A, Mode B).
!     - Complex multipliers (QP modes) require phase-averaging
!       of two linearly independent Floquet modes (deferred).
!     - Pressure work and boundary terms are neglected (standard
!       for incompressible flow with homogeneous BCs on u').
!     - The orbit-averaging identity is approximate when E_hat(t)
!       varies over the period (exact for steady base flows).

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

!     ----- Perturbation arrays -----
      real, dimension(lv) :: vx_d, vy_d, vz_d
      real, dimension(lv) :: vx_zero, vy_zero, vz_zero

!     ----- Budget accumulators (allocatable to avoid stack overflow) --
      real, allocatable, dimension(:,:) :: budget_avg
      real, allocatable, dimension(:) :: energy_avg
      real, dimension(10) :: integrals
      real, dimension(lv) :: prod_x, prod_y, prod_z, diss_tmp

!     ----- Eigenvalue data -----
      real :: sigma_r, omega, period, growth_corr, weight
      real :: E_bar, sigma_check, rel_error

!     ----- Miscellaneous -----
      real :: glsc2
      integer :: i, j, k, n, mode, m, col
      character(len=80) :: filename
      character(len=6) :: mode_str
      character(len=3) :: mode_str2

      logical, save :: init
      data init/.false./

      n = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt

!     #####
!     #####     READ EIGENVALUE
!     #####

      call read_eigenvalue(sigma_r, omega)
      if (nid == 0) then
         write (6, *) 'Floquet PKE budget: sigma_r =', sigma_r
         write (6, *) 'Floquet PKE budget: omega   =', omega
      end if

!     #####
!     #####     LOAD BASE FLOW AND PREPARE SOLVER
!     #####

      write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
      call load_fld(filename)

!     --> Set up linearized solver (computes nsteps, dt from param(10))
      call prepare_linearized_solver

      period = nsteps * dt
      if (nid == 0) write (6, *)
     $   'Orbit period T =', period, ' nsteps =', nsteps

!     --> Setup linearized solver flags
      ifpert = .true.; ifadj = .false.
      call bcast(ifpert, lsize); call bcast(ifadj, lsize)

!     --> Zero arrays for imaginary part (real multiplier only)
      call rzero(vx_zero, lv)
      call rzero(vy_zero, lv)
      call rzero(vz_zero, lv)

!     --> Allocate orbit storage on first call
      if (.not. init) then
         call allocate_orbit(nsteps)
      end if

!     --> Allocate budget accumulators on heap
      allocate(budget_avg(lv, 10))
      allocate(energy_avg(lv))

!     #####
!     #####     LOOP OVER MODES
!     #####

      do mode = 1, maxmodes

         budget_avg = 0.0d0
         energy_avg = 0.0d0
         integrals  = 0.0d0

!        --> Format mode number
         write (mode_str, '(i5.5)') mode

!        --> Load eigenmode (real part only for real multiplier)
         write (filename, '(a,a,a,a)')
     $      'dRe', trim(SESSION), '0.f', trim(mode_str)
         call load_fld(filename)

!        --> Pass eigenmode as perturbation IC
         call opcopy(vxp(:,1), vyp(:,1), vzp(:,1),
     $      vx, vy, vz)
         if (ifto) call copy(tp(:,:,1), t, nt)

!        --> Reload base flow IC into vx,vy,vz for orbit evolution
         if (.not. init) then
            write (filename, '(a,a,a)')
     $         'BF_', trim(SESSION), '0.f00001'
            call load_fld(filename)
            ifbase = .true.
         else
            ifbase = .false.
!           --> Restore first orbit step from stored orbit
            call orbit_restore(1)
         end if

!        ─────────────────────────────────────────
!        ORBIT INTEGRATION + BUDGET ACCUMULATION
!        ─────────────────────────────────────────

         time = 0.0d0
         do istep = 1, nsteps

!           --> Log progress
            if (nid == 0) write (6,
     $         "(' PKE_FLOQUET mode',I3,':',I6,'/',I6)")
     $         mode, istep, nsteps

!           --> Advance (BF + perturbation simultaneously)
            call nekstab_usrchk()
            call nek_advance()

!           --> Store/load orbit
            if (.not. init) then
               call orbit_store(istep)
            else
               call orbit_restore(istep)
            end if

!           --> Copy current base flow to ubase/vbase/wbase
!               (compute_production reads from ubase,vbase,wbase)
            call opcopy(ubase, vbase, wbase, vx, vy, vz)

!           --> Get current perturbation
            call opcopy(vx_d, vy_d, vz_d,
     $         vxp(:,1), vyp(:,1), vzp(:,1))

!           --> Growth correction + trapezoidal weight
            growth_corr = exp(-2.0d0 * sigma_r * time)
            weight = growth_corr * dt / period

!           --> Production terms (9 = 3 components x 3 gradients)
            do j = 1, 3
               call compute_production(vx_d, vy_d, vz_d,
     $            vx_zero, vy_zero, vz_zero, j,
     $            prod_x, prod_y, prod_z)
               col = (j - 1) * 3
               do i = 1, n
                  budget_avg(i,col+1) = budget_avg(i,col+1)
     $               + weight * prod_x(i)
                  budget_avg(i,col+2) = budget_avg(i,col+2)
     $               + weight * prod_y(i)
                  budget_avg(i,col+3) = budget_avg(i,col+3)
     $               + weight * prod_z(i)
               end do
            end do

!           --> Dissipation term (1 scalar field)
            call compute_dissipation(vx_d, vy_d, vz_d,
     $         vx_zero, vy_zero, vz_zero, diss_tmp)
            do i = 1, n
               budget_avg(i,10) = budget_avg(i,10)
     $            + weight * diss_tmp(i)
            end do

!           --> Perturbation kinetic energy
            do i = 1, n
               energy_avg(i) = energy_avg(i) + weight * 0.5d0
     $            * (vx_d(i)**2 + vy_d(i)**2 + vz_d(i)**2)
            end do

         end do ! istep

!        --> Mark orbit as stored after first mode
         if (.not. init) then
            ifbase = .false.
            init = .true.
         end if

!        ─────────────────────────────────────────
!        INTEGRATE AND VERIFY
!        ─────────────────────────────────────────

         do i = 1, 10
            integrals(i) = glsc2(bm1, budget_avg(:,i), n)
         end do
         E_bar = glsc2(bm1, energy_avg, n)

         sigma_check = (sum(integrals(1:9)) - integrals(10))
     $      / (2.0d0 * E_bar)

         rel_error = abs(sigma_check - sigma_r)
     $      / max(abs(sigma_r), 1.0d-30)

         if (nid == 0) then
            write (6, *) ''
            write (6, *)
     $         '=== Orbit-averaged PKE budget (mode', mode,
     $         ') ==='
            write (6, '(A,E15.7)')
     $         '  sigma_r (eigenvalue) = ', sigma_r
            write (6, '(A,E15.7)')
     $         '  sigma_r (budget)     = ', sigma_check
            write (6, '(A,E15.7)')
     $         '  relative error       = ', rel_error
            write (6, '(A,E15.7)')
     $         '  E_bar                = ', E_bar
            write (6, '(A,E15.7)')
     $         '  P_bar (total)        = ',
     $         sum(integrals(1:9))
            write (6, '(A,E15.7)')
     $         '  D_bar                = ', integrals(10)
            do i = 1, 10
               write (6, '(A,I2,A,E15.7)')
     $            '  Integral ', i, ' = ', integrals(i)
            end do
            write (6, *) ''
         end if

!        ─────────────────────────────────────────
!        OUTPUT
!        ─────────────────────────────────────────

!        --> Output spatial budget fields (same format as steady)
         if (if3d) then
            k = 9
         else
            k = 6
         end if

         do i = 1, k, 3
            call opcopy(vx, vy, vz,
     $         budget_avg(:,i), budget_avg(:,i+1),
     $         budget_avg(:,i+2))
            write (mode_str2, "('F',I2.2)") mode
            call outpost(vx, vy, vz, pr, t, mode_str2)
         end do

!        --> Output scalar integrals to file
         if (nid == 0) then
            write (filename, '(A,A,A,A)')
     $         'PKE_floquet_', trim(SESSION), '0.f',
     $         trim(mode_str)
            open (101, file=filename, form='formatted')
            do i = 1, 10
               write (101, '(1E15.7)') integrals(i)
            end do
            write (101, '(1E15.7)') sum(integrals(1:9))
            write (101, '(1E15.7)')
     $         sum(integrals(1:9)) - integrals(10)
            write (101, '(1E15.7)') E_bar
            write (101, '(1E15.7)') sigma_check
            close (101)
         end if

      end do ! mode

!     --> Cleanup
      if (allocated(budget_avg)) deallocate(budget_avg)
      if (allocated(energy_avg)) deallocate(energy_avg)
      if (allocated(uor)) deallocate(uor, vor, wor)
      if (allocated(tor)) deallocate(tor)

      return
      end subroutine stability_energy_budget_floquet

!----------------------------------------------------------------------

      subroutine compute_velocity_gradient_tensor(
     $   vx_in, vy_in, vz_in,
     $   dudx, dudy, dudz,
     $   dvdx, dvdy, dvdz,
     $   dwdx, dwdy, dwdz)
!
!     Compute the full velocity gradient tensor and smooth at element
!     interfaces using dsavg.  Uses the 5-argument form of gradm1.
!
!     INPUT
!       vx_in, vy_in, vz_in : velocity field
!
!     OUTPUT
!       dudx..dwdz : 9 gradient components
!
         implicit none
         include 'SIZE'
         include 'TOTAL'

         real, dimension(lx1*ly1*lz1*lelv) :: vx_in, vy_in, vz_in
         real, dimension(lx1*ly1*lz1*lelv) :: dudx, dudy, dudz
         real, dimension(lx1*ly1*lz1*lelv) :: dvdx, dvdy, dvdz
         real, dimension(lx1*ly1*lz1*lelv) :: dwdx, dwdy, dwdz

         call gradm1(dudx, dudy, dudz, vx_in, nelv)
         call gradm1(dvdx, dvdy, dvdz, vy_in, nelv)
         call gradm1(dwdx, dwdy, dwdz, vz_in, nelv)

         call dsavg(dudx); call dsavg(dudy); call dsavg(dudz)
         call dsavg(dvdx); call dsavg(dvdy); call dsavg(dvdz)
         call dsavg(dwdx); call dsavg(dwdy); call dsavg(dwdz)

      end subroutine compute_velocity_gradient_tensor
!----------------------------------------------------------------------

      subroutine compute_dissipation(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, dissipation)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv) :: vx_dRe, vy_dRe, vz_dRe
         real, dimension(lv) :: vx_dIm, vy_dIm, vz_dIm
      
         real, dimension(lv) :: Laplacian_ax, Laplacian_ay, Laplacian_az
         real, dimension(lv) :: Laplacian_bx, Laplacian_by, Laplacian_bz
      
         real, dimension(lv) :: dissipation, dummy
      
      !     --> Compute Laplacians.
         call compute_laplacian(vx_dRe, Laplacian_ax)
         call compute_laplacian(vy_dRe, Laplacian_ay)
         call compute_laplacian(vz_dRe, Laplacian_az)
      
         call compute_laplacian(vx_dIm, Laplacian_bx)
         call compute_laplacian(vy_dIm, Laplacian_by)
         call compute_laplacian(vz_dIm, Laplacian_bz)
      
      !     --> Compute dissipation term.
         dissipation = 0.0d+00
      
         dummy = vx_dRe*Laplacian_ax + vx_dIm*Laplacian_bx
         dissipation = dissipation + dummy
      
         dummy = vy_dRe*Laplacian_ay + vy_dIm*Laplacian_by
         dissipation = dissipation + dummy
      
         dummy = vz_dRe*Laplacian_az + vz_dIm*Laplacian_bz
         dissipation = dissipation + dummy
      
         dissipation = 0.5*dissipation*param(2)/param(1)
      
         return
      end subroutine compute_dissipation
      
      subroutine compute_production(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm, component, prod_x, prod_y, prod_z)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv) :: vx_dRe, vy_dRe, vz_dRe
         real, dimension(lv) :: vx_dIm, vy_dIm, vz_dIm
         real, dimension(lv) :: prod_x, prod_y, prod_z
         real, dimension(lv) :: dcdx, dcdy, dcdz
         integer :: component
      
         if (component == 1) then
            call gradm1(dcdx, dcdy, dcdz, ubase, nelv)
      
            prod_x = -0.5*(vx_dRe**2 + vx_dIm**2)*dcdx
            prod_y = -0.5*(vx_dRe*vy_dRe + vy_dIm*vx_dIm)*dcdy
            prod_z = -0.5*(vx_dRe*vz_dRe + vz_dIm*vx_dIm)*dcdz
      
         else if (component == 2) then
            call gradm1(dcdx, dcdy, dcdz, vbase, nelv)
      
            prod_x = -0.5*(vx_dRe*vy_dRe + vy_dIm*vx_dIm)*dcdx
            prod_y = -0.5*(vy_dRe**2 + vy_dIm**2)*dcdy
            prod_z = -0.5*(vy_dRe*vz_dRe + vz_dIm*vy_dIm)*dcdz
      
         else if (component == 3) then
            call gradm1(dcdx, dcdy, dcdz, wbase, nelv)
      
            prod_x = -0.5*(vx_dRe*vz_dRe + vz_dIm*vx_dIm)*dcdx
            prod_y = -0.5*(vy_dRe*vz_dRe + vz_dIm*vy_dIm)*dcdy
            prod_z = -0.5*(vz_dRe**2 + vz_dIm**2)*dcdz
         end if
      
         return
      end subroutine compute_production
      
      subroutine compute_gradients(u, dudx, dudy, dudz)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv) :: u
         real, dimension(lv) :: dudx, dudy, dudz
      
         call gradm1(dudx, dudy, dudz, u)
         call dsavg(dudx); call dsavg(dudy); call dsavg(dudz)
      
         return
      end subroutine compute_gradients
      
      subroutine compute_laplacian(a, Lap_a)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real, dimension(lv) :: a
         real, dimension(lv) :: dadx, dady, dadz
         real, dimension(lv) :: d2adx2, d2ady2, d2adz2
         real, dimension(lv) :: Lap_a, wrk1, wrk2
      
         call compute_gradients(a, dadx, dady, dadz)
         call compute_gradients(dadx, d2adx2, wrk1, wrk2)
         call compute_gradients(dady, wrk1, d2ady2, wrk2)
         call compute_gradients(dadz, wrk1, wrk2, d2adz2)
      
         Lap_a = d2adx2 + d2ady2 + d2adz2
      
         return
      end subroutine compute_laplacian
      !----------------------------------------------------------------------
