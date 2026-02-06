      !-----------------------------------------------------------------------
      ! noise.f90 — Random perturbation seeding
      !
      ! Purpose:
      !   Generates pseudo-random perturbation fields for initial seeding
      !   of stability computations. Supports scalar, vector, and
      !   symmetric seed modes.
      !
      ! Public interface:
      !   add_noise_scal    — add noise to a scalar field
      !   op_add_noise      — add noise to velocity vector field
      !   add_symmetric_seed — generate symmetric initial perturbation
      !   mth_rand          — deterministic pseudo-random number generator
      !
      ! Dependencies:
      !   SIZE, TOTAL (TSTEP, PARALLEL, INPUT, SOLN, GEOM)
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! add_noise_scal — Add pseudo-random noise to a scalar field
      !
      ! Purpose:
      !   Perturbs a scalar field with deterministic noise based on
      !   element geometry and frequency coefficients. Applies DSS
      !   averaging and boundary conditions.
      !
      ! Arguments:
      !   qin [inout] — scalar field to perturb (1D layout, lx1*ly1*lz1*lelt)
      !   fc1 [in]    — frequency coefficient 1
      !   fc2 [in]    — frequency coefficient 2
      !   fc3 [in]    — frequency coefficient 3
      !-----------------------------------------------------------------------
      subroutine add_noise_scal(qin, fc1, fc2, fc3)
         implicit none
         include 'SIZE'
         include 'TSTEP' ! TIME, DT
         include 'PARALLEL' ! LGLEL
         include 'INPUT' ! if3D
         include 'SOLN' ! VX, VY, VZ, VMULT
         include 'GEOM' ! XM1, YM1, ZM1
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

      !-----------------------------------------------------------------------
      ! op_add_noise — Add pseudo-random noise to velocity vector field
      !
      ! Purpose:
      !   Perturbs velocity components with deterministic noise using
      !   different frequency coefficients per component. Applies DSS
      !   averaging and velocity boundary conditions.
      !
      ! Arguments:
      !   qx [inout] — x-velocity to perturb
      !   qy [inout] — y-velocity to perturb
      !   qz [inout] — z-velocity to perturb
      !-----------------------------------------------------------------------
      subroutine op_add_noise(qx, qy, qz)
         implicit none
         include 'SIZE'
         include 'TSTEP' ! TIME, DT
         include 'PARALLEL' ! LGLEL
         include 'INPUT' ! if3D
         include 'SOLN' ! VX, VY, VZ, VMULT
         include 'GEOM' ! XM1, YM1, ZM1

         real, intent(inout), dimension(lx1, ly1, lz1, lelv) :: qx, qy, qz
         integer iel, ieg, il, jl, kl, nv
         real xl(LDIM), mth_rand, fc(3), nmin, nmax, glmax, glmin

         nv = nx1*ny1*nz1*nelv

         do iel = 1, NELV
            do kl = 1, NZ1
               do jl = 1, NY1
                  do il = 1, NX1
                     xl(1) = XM1(il, jl, kl, iel)
                     xl(2) = YM1(il, jl, kl, iel)
                     if (if3D) xl(NDIM) = ZM1(il, jl, kl, iel)

                     ieg = LGLEL(iel)
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

      !-----------------------------------------------------------------------
      ! add_symmetric_seed — Generate symmetric initial perturbation
      !
      ! Purpose:
      !   Creates a divergence-free symmetric perturbation field suitable
      !   for seeding stability analysis. In 2D, generates y-dependent
      !   modes; in 3D, adds spanwise modulation.
      !
      ! Arguments:
      !   qx [inout] — x-velocity perturbation
      !   qy [inout] — y-velocity perturbation
      !   qz [inout] — z-velocity perturbation
      !   qp [inout] — pressure perturbation
      !-----------------------------------------------------------------------
      subroutine add_symmetric_seed(qx, qy, qz, qp)
         implicit none
         include "SIZE"
         include "TOTAL"
         real, intent(inout), dimension(lx1, ly1, lz1, lelv) :: qx, qy, qz, qp
         integer iel, il, jl, kl, nv
         real xlx, yly, zlz, alpha, x, y, z
         real glsc3, amp

         nv = NX1*NY1*NZ1*NELV
         xlx = xmx - xmn
         yly = ymx - ymn
         zlz = zmx - zmn

      !     --> Create the initial velocity perturbation.

         do iel = 1, NELV
            do kl = 1, NZ1
               do jl = 1, NY1
                  do il = 1, NX1

                     x = XM1(il, jl, kl, iel)
                     y = YM1(il, jl, kl, iel)

                     if (if3D) then
      !     -> 3D: spanwise-modulated perturbation
                        z = ZM1(il, jl, kl, iel)
                        alpha = 2*pi/zlz
                        qx(il, jl, kl, iel) = cos(alpha*z)*sin(2.*pi*y)
                        qz(il, jl, kl, iel) = -(2.*pi)/(alpha)*cos(alpha*z)*cos(2.*pi*y)
                        qp(il, jl, kl, iel) = cos(alpha*z)*cos(2.*pi*y)
                     else
      !     -> 2D: y-dependent perturbation only
                        qx(il, jl, kl, iel) = sin(2.*pi*y/yly)
                        qy(il, jl, kl, iel) = -cos(2.*pi*y/yly)
                        qp(il, jl, kl, iel) = cos(2.*pi*y/yly)
                     end if

                  end do
               end do
            end do
         end do

         amp = glsc3(qx, qx, bm1s, nv) + glsc3(qy, qy, bm1s, nv)
         if (if3D) amp = amp + glsc3(qz, qz, bm1s, nv)
         amp = 1e-6/(0.50d0*amp)
         call opcmult(qx, qy, qz, amp)
         call cmult(qp, amp, nv)

         return
      end subroutine add_symmetric_seed
      !-----------------------------------------------------------------------

      ! Deterministic pseudo-random number generator based on coordinates
      real function mth_rand(ix, iy, iz, ieg, xl, fc)
         implicit none
         include 'SIZE'
         include 'INPUT' ! if3D
         integer, intent(in) :: ix, iy, iz, ieg
         real, intent(in) :: xl(LDIM), fc(3)

         mth_rand = fc(1)*(ieg + xl(1)*sin(xl(2))) + fc(2)*ix*iy + fc(3)*ix
         if (if3D) mth_rand = fc(1)*(ieg + xl(NDIM)*sin(mth_rand)) + fc(2)*iz*ix + fc(3)*iz
         mth_rand = cos(1.e3*sin(1.e3*sin(mth_rand)))

         return
      end function mth_rand
      !-----------------------------------------------------------------------
