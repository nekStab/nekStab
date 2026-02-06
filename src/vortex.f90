      !-----------------------------------------------------------------------
      ! vortex.f90 — Vortex identification criteria
      !
      ! Purpose:
      !   Computes various vortex identification fields from velocity
      !   data: lambda2, Q-criterion, delta, swirling strength, omega,
      !   and symmetric/antisymmetric decompositions.
      !
      ! Public interface:
      !   vortex_core             — dispatch to selected vortex criterion
      !   compute_omega_jc        — omega criterion (Frobenius norm)
      !   compute_omega           — omega criterion (legacy)
      !   compute_q               — Q criterion (Hunt et al., 1988)
      !   compute_delta           — delta criterion (Chong et al., 1990)
      !   compute_swirling        — swirling strength (Zhou et al., 1999)
      !   compute_symmetricVec    — symmetric part of velocity gradient
      !   compute_assymetricVec   — antisymmetric part of velocity gradient
      !
      ! Dependencies:
      !   SIZE, TOTAL
      !-----------------------------------------------------------------------

      module nekstab_vortex
         private
         public :: vortex_core, compute_omega_jc,
     $             compute_omega, compute_symmetricVec,
     $             compute_assymetricVec, compute_q,
     $             compute_delta, compute_swirling,
     $             compute_antisymmetric, compute_symmetric,
     $             compute_firstInv, compute_secondInv,
     $             compute_thirdInv, cubicLambdaCi,
     $             quadLambdaCi
      contains

      !-----------------------------------------------------------------------
      ! vortex_core — Dispatch to selected vortex identification method
      !
      ! Arguments:
      !   l2     [out]   — output vortex field
      !   vortex [in]    — name of vortex criterion to compute
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

      !-----------------------------------------------------------------------
      ! compute_omega_jc — Omega criterion using Frobenius norms
      !
      ! Arguments:
      !   omega [out] — omega vortex identification field
      !-----------------------------------------------------------------------
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

      !-----------------------------------------------------------------------
      ! compute_omega — Omega criterion (legacy implementation)
      !
      ! Arguments:
      !   l2 [out] — omega vortex identification field
      !-----------------------------------------------------------------------
      subroutine compute_omega(l2)
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

      !-----------------------------------------------------------------------
      ! compute_symmetricVec — Symmetric part of velocity gradient
      !
      ! Arguments:
      !   l2 [out] — symmetric part field
      !-----------------------------------------------------------------------
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

      !-----------------------------------------------------------------------
      ! compute_assymetricVec — Antisymmetric part of velocity gradient
      !
      ! Arguments:
      !   l2 [out] — antisymmetric part field
      !-----------------------------------------------------------------------
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

      !-----------------------------------------------------------------------
      ! compute_q — Q criterion (Hunt, Wray & Moin, CTR-S88 1988)
      !
      ! Purpose:
      !   Positive second invariant of velocity gradient tensor.
      !
      ! Arguments:
      !   l2 [out] — Q criterion field
      !-----------------------------------------------------------------------
      subroutine compute_q(l2)
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

      !-----------------------------------------------------------------------
      ! compute_delta — Delta criterion (Chong, Perry & Cantwell, 1990)
      !
      ! Purpose:
      !   Discriminant criterion: complex eigenvalues of velocity
      !   gradient tensor.
      !
      ! Arguments:
      !   l2 [out] — delta criterion field
      !-----------------------------------------------------------------------
      subroutine compute_delta(l2)
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

      !-----------------------------------------------------------------------
      ! compute_swirling — Swirling strength (Zhou et al., JFM 1999)
      !
      ! Purpose:
      !   Imaginary part of complex eigenvalues of velocity gradient
      !   tensor, presented as lambda_ci^2.
      !
      ! Arguments:
      !   l2 [out] — swirling strength field
      !-----------------------------------------------------------------------
      subroutine compute_swirling(l2)
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
                  R = R1 + (P1**3)*2/27 - P1*Q1/3
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
                  call compute_thirdInv(R1, l)
                  R1 = -R1 !negative sign
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

      !-----------------------------------------------------------------------
      ! compute_antisymmetric — Antisymmetric part norm of gradient
      !
      ! Arguments:
      !   B [out] — Frobenius norm of antisymmetric part
      !   l [in]  — point index within element
      !-----------------------------------------------------------------------
      subroutine compute_antisymmetric(B, l)
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

      !-----------------------------------------------------------------------
      ! compute_symmetric — Symmetric part norm of gradient
      !
      ! Arguments:
      !   A [out] — Frobenius norm of symmetric part
      !   l [in]  — point index within element
      !-----------------------------------------------------------------------
      subroutine compute_symmetric(A, l)
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

      !-----------------------------------------------------------------------
      ! compute_firstInv — First invariant (trace) of gradient tensor
      !
      ! Arguments:
      !   a [out] — tr(D) = d11 + d22 + d33
      !   l [in]  — point index within element
      !-----------------------------------------------------------------------
      subroutine compute_firstInv(a, l)
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

      !-----------------------------------------------------------------------
      ! compute_secondInv — Second invariant of gradient tensor
      !
      ! Arguments:
      !   a [out] — II_D
      !   l [in]  — point index within element
      !-----------------------------------------------------------------------
      subroutine compute_secondInv(a, l)
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

      !-----------------------------------------------------------------------
      ! compute_thirdInv — Third invariant (determinant) of gradient
      !
      ! Arguments:
      !   a [out] — det(D)
      !   l [in]  — point index within element
      !-----------------------------------------------------------------------
      subroutine compute_thirdInv(a, l)
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

      !-----------------------------------------------------------------------
      ! cubicLambdaCi — Complex eigenvalue of 3x3 gradient tensor
      !
      ! Arguments:
      !   b   [in]  — first invariant coefficient
      !   c   [in]  — second invariant coefficient
      !   d   [in]  — third invariant coefficient
      !   lci [out] — imaginary part of complex eigenvalue
      !-----------------------------------------------------------------------
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

      !-----------------------------------------------------------------------
      ! quadLambdaCi — Complex eigenvalue of 2x2 gradient tensor
      !
      ! Arguments:
      !   b   [in]  — trace coefficient
      !   c   [in]  — determinant coefficient
      !   lci [out] — imaginary part of complex eigenvalue
      !-----------------------------------------------------------------------
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
            f = sqrt(abs(d))
            x1 = -b/2./a + f*ci/2./a
            lci = aimag(x1)
         end if

      end subroutine quadLambdaCi
      !-----------------------------------------------------------------------

      end module nekstab_vortex
