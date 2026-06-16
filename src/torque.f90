!-----------------------------------------------------------------------
! torque.f90 — Drag and torque computation
!
! Purpose:
!   Computes aerodynamic drag and torque coefficients on wall
!   boundaries during DNS. Outputs pressure and viscous
!   contributions separately.
!
! Public interface:
!   nekStab_torque     — compute and output drag/torque
!   nekStab_define_obj — define wall boundary objects
!
! Dependencies:
!   nekstab_nek_bridge
!-----------------------------------------------------------------------

module nekstab_torque_mod
   use nekstab_nek_bridge, only: nekStab_log, nekStab_error, nid, &
                                 NEKSTAB_UNIT_TORQUE, nx1, ny1, nz1, nelv, &
                                 nelt, lx1, ly1, lz1, lelv, lelt, vdiff, param, nobj, &
                                 pr, dpdx_mean, dpdy_mean, dpdz_mean, xm1, ym1, &
                                 zm1, if3d, ifaxis, vx, vy, vz, maxobj, dragpx, &
                                 dragvx, dragx, dragpy, dragvy, dragy, dragpz, &
                                 dragvz, dragz, torqpx, torqvx, torqx, torqpy, &
                                 torqvy, torqy, torqpz, torqvz, torqz, ifield, &
                                 ndim, cbc, nio, istep, time, boundaryID
   implicit none
   private
   public :: nekStab_torque, nekStab_define_obj
   logical, save :: torque_initialized = .false.
   integer, save :: torque_bIDs(1), torque_iobj_wall(1)
   real, save :: torque_x0(3), torque_scale

contains

!-----------------------------------------------------------------------
! nekStab_torque -- Compute drag and torque on wall boundaries
!
! Purpose:
!   Computes pressure and viscous contributions to drag and torque on
!   all 'W  ' (wall) faces, then appends one row per step to `fname`.
!   Output columns (3D): istep, time, then per axis
!   drag, drag_pressure, drag_viscous (x,y,z) and torque (x,y,z).
!
! Arguments:
!   fname [in] -- output filename for drag/torque data
!
! ---------------------------------------------------------------------
! WHY THIS ROUTINE DOES NOT USE create_obj / nmember / object
! ---------------------------------------------------------------------
! Nek5000's stock drag pattern is: mark wall faces in `boundaryID`,
! call create_obj() to gather them into an "object" (nobj, nmember,
! object(...)), then loop the object members calling drgtrq().
!
! That pattern is BROKEN in this v2.0.0 build and silently returns
! zero forces. Root cause (debugged 2026-06, cube_5/400):
!
!   The "wrap all Nek state in a bridge module" refactor put
!   `include 'TOTAL'` *inside* a Fortran module (nekstab_nek_bridge),
!   and the link uses `-Wl,--allow-multiple-definition`. Nek core
!   files (e.g. bdry.f) also `include 'TOTAL'`. The link flag lets the
!   duplicate COMMON blocks coexist instead of erroring, and their
!   storage is NOT merged: module-compiled code and Nek-core-compiled
!   code see DIFFERENT copies of some commons.
!
!   For the object machinery this is fatal and deceptive: create_obj
!   (Nek core) correctly scans boundaryID, finds the wall faces, and
!   writes nmember=Nwall into *core's* copy of common /input2/
!   (INPUT: ...,nmember(maxobj),nobj,...). But the integration loop
!   here (module-compiled) reads the *module's* copy, where nmember=0
!   -> zero members -> zero force. `nobj` happens to alias (it sat at
!   a coinciding offset) which made the split maddening to find:
!   create_obj reports 216 wall faces and nobj=1, yet nmember=0.
!
!   We sidestep the split entirely: `cbc` (the BC character array) IS
!   read consistently through the bridge, so we loop wall faces by
!   `cbc(ifc,ie,1)=='W  '` directly and call drgtrq() per face. No
!   object, no /input2/ dependency, robust for any wall geometry.
!
! Second regression fixed here: the scale factor was dropped. nekStab
! 1.1 set `scale = 2` (the 2/(rho U^2 A) drag-coefficient
! normalization); the refactor left `torque_scale` uninitialised, so
! cmult(dgtq,torque_scale,12) multiplied every contribution by 0 --
! by itself enough to zero the forces. Restored to 2 below.
!
! Verified: with both fixes, lift_drag.dat reproduces the production
! 1.1 values bit-for-bit at the restart time (cube_5/400, t=25000).
!-----------------------------------------------------------------------
   subroutine nekStab_torque(fname)
      character(len=*), intent(in) :: fname
      character(len=256) :: fname_local ! Local copy for ifx compatibility

      ! becuase of this block we can not use implicit none in this routine !
      real sij(lx1*ly1*lz1*6*lelv)
      real pm1(lx1, ly1, lz1, lelv)
      real xm0(lx1, ly1, lz1, lelt), ym0(lx1, ly1, lz1, lelt), zm0(lx1, ly1, lz1, lelt)
      integer, parameter :: lr = lx1*ly1*lz1
      real ur(lr), us(lr), ut(lr), vr(lr), vs(lr), vt(lr), wr(lr), ws(lr), wt(lr)
      real torque_scale_vf(3)
      real dgtq(3, 4)
   !!!

      common/scrns/sij
      common/scrcg/pm1
      common/scrsf/xm0, ym0, zm0
      common/scruz/ur, us, ut, vr, vs, vt, wr, ws, wt
      common/cvflow_r/torque_scale_vf

      ! Local variables only
      integer :: nv, i, ie, iobj, memtot, mem, ieg, ifc, nij, ierr
      real :: glmin, glmax, x1min, x2min, x3min, x1max, x2max, x3max, w1(0:maxobj)

      ! nv matches the nekstab_krylov_subspace global set at init; explicit
      ! decl required because this file compiles before that module.
      nv = nx1*ny1*nz1*nelv

      if (.not. torque_initialized) then
         fname_local = fname ! Copy to local variable for ifx compatibility
         call nekStab_log('Initializing torque routine...')
         call nekStab_log('About to open torque file: '//trim(fname_local))
         if (nid == 0) then
            open (NEKSTAB_UNIT_TORQUE, file=trim(fname_local), action='write', &
                  status='replace', iostat=ierr)
            if (ierr /= 0) then
               call nekStab_error('Could not open torque file')
            end if
         end if
         torque_scale = 2.0d0 ! Cd/Cl normalization (2/(rho U^2 A); dropped in v2.0.0 refactor)
         torque_x0(1) = 0.0d0 ! moment reference point
         torque_x0(2) = 0.0d0
         torque_x0(3) = 0.0d0
         nobj = 1 ! single wall object; faces gathered directly from cbc below
         ! NOTE: create_obj/nmember/object are NOT used. The v2.0.0 module/bridge
         ! refactor splits common /input2/ between this module and Nek core, so
         ! create_obj (core) populates an nmember the module cannot read. We loop
         ! wall faces directly via cbc (which the bridge reads correctly).
         call cfill(vdiff, param(2), nv) ! Fill up viscous array with default on initialization
         torque_initialized = .true.
         call nekStab_log('Torque init: scale=2, wall faces from cbc==W.')
      end if

      call mappr(pm1, pr, xm0, ym0) ! map pressure onto Mesh 1

      if (param(55) /= 0) then ! Add mean_pressure_gradient.X to p:
         dpdx_mean = -torque_scale_vf(1)
         dpdy_mean = -torque_scale_vf(2)
         dpdz_mean = -torque_scale_vf(3)
      end if

      call add2s2(pm1, xm1, dpdx_mean, nv) ! Doesn't work if object is cut by periodic boundary.
      call add2s2(pm1, ym1, dpdy_mean, nv) ! In this case, set ._mean=0 and compensate in post.
      call add2s2(pm1, zm1, dpdz_mean, nv)

      nij = 3 ! Compute sij
      if (if3d .or. ifaxis) nij = 6
      call comp_sij(sij, nij, vx, vy, vz, ur, us, ut, vr, vs, vt, wr, ws, wt)

      call cadd2(xm0, xm1, -torque_x0(1), nv)
      call cadd2(ym0, ym1, -torque_x0(2), nv)
      call cadd2(zm0, zm1, -torque_x0(3), nv)

      x1min = glmin(xm0(1, 1, 1, 1), nv)
      x2min = glmin(ym0(1, 1, 1, 1), nv)
      x3min = glmin(zm0(1, 1, 1, 1), nv)
      x1max = glmax(xm0(1, 1, 1, 1), nv)
      x2max = glmax(ym0(1, 1, 1, 1), nv)
      x3max = glmax(zm0(1, 1, 1, 1), nv)

      do i = 0, maxobj
         dragpx(i) = 0 ! Pressure drag x
         dragvx(i) = 0 ! Viscous drag x
         dragx(i) = 0 ! Total drag x
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
      iobj = 1
      do ie = 1, nelv ! loop local elements; gather wall faces from cbc
         do ifc = 1, 2*ndim
            if (cbc(ifc, ie, 1) == 'W  ') then
               call drgtrq(dgtq, xm0, ym0, zm0, sij, pm1, vdiff, ifc, ie)
               call cmult(dgtq, torque_scale, 12)

               dragpx(iobj) = dragpx(iobj) + dgtq(1, 1) ! Pressure
               dragpy(iobj) = dragpy(iobj) + dgtq(2, 1)
               dragpz(iobj) = dragpz(iobj) + dgtq(3, 1)

               dragvx(iobj) = dragvx(iobj) + dgtq(1, 2) ! Viscous
               dragvy(iobj) = dragvy(iobj) + dgtq(2, 2)
               dragvz(iobj) = dragvz(iobj) + dgtq(3, 2)

               torqpx(iobj) = torqpx(iobj) + dgtq(1, 3) ! Pressure
               torqpy(iobj) = torqpy(iobj) + dgtq(2, 3)
               torqpz(iobj) = torqpz(iobj) + dgtq(3, 3)

               torqvx(iobj) = torqvx(iobj) + dgtq(1, 4) ! Viscous
               torqvy(iobj) = torqvy(iobj) + dgtq(2, 4)
               torqvz(iobj) = torqvz(iobj) + dgtq(3, 4)
            end if
         end do
      end do

      ! Sum contributions from all processors
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

      do i = 1, nobj ! Combine results
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

      ! Match output format to htps routine: use 1p in format string for floats
      do i = 1, nobj ! Write results
         if (nio == 0) then
            if (if3D .or. ifaxis) then ! 3D or axisymmetric
               write (NEKSTAB_UNIT_TORQUE, "(i8,1p19E15.7)") istep, time, &
                  dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), dragz(i), dragpz(i), dragvz(i), &
                  torqx(i), torqpx(i), torqvx(i), torqy(i), torqpy(i), torqvy(i), torqz(i), torqpz(i), torqvz(i)
            else ! 2D or not axisymmetric
               write (NEKSTAB_UNIT_TORQUE, "(i8,1p10E15.7)") istep, time, &
                  dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), torqz(i), torqpz(i), torqvz(i)
            end if ! if3D .or. ifaxis
         end if ! nio == 0
      end do ! i = 1, nobj

      if (nio == 0) flush (NEKSTAB_UNIT_TORQUE) ! Ensure data is written to disk

   end subroutine nekStab_torque
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! nekStab_define_obj -- Assign boundaryID=1 to all wall ('W  ') faces.
!
! Kept for API compatibility: legacy cases still `call nekStab_define_obj`
! from usrdat2. nekStab_torque NO LONGER depends on it -- and in fact it
! is effectively a no-op for the drag path in this v2.0.0 build, because
! a boundaryID write from module-compiled code does not reach the copy of
! the common that Nek-core create_obj reads (the same module/core common
! split documented in nekStab_torque above). Harmless to keep; do not
! rely on it for object creation until the bridge common-split is fixed.
!-----------------------------------------------------------------------
   subroutine nekStab_define_obj
      integer :: iel, iface

      do iel = 1, nelt
         do iface = 1, 2*ndim
            if (cbc(iface, iel, 1) == 'W  ') then
               boundaryID(iface, iel) = 1
            end if
         end do
      end do

   end subroutine nekStab_define_obj
!-----------------------------------------------------------------------

end module nekstab_torque_mod
