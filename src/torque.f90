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
      !   SIZE, TOTAL
      !-----------------------------------------------------------------------

      module nekstab_torque_mod
         implicit none
         private
         public :: nekStab_torque, nekStab_define_obj
      contains

      !-----------------------------------------------------------------------
      ! nekStab_torque — Compute drag and torque on wall boundaries
      !
      ! Purpose:
      !   Computes pressure and viscous contributions to drag and torque
      !   on wall-type boundary objects. Initializes on first call,
      !   then appends results to a file each timestep.
      !
      ! Arguments:
      !   fname [in] — output filename for drag/torque data
      !-----------------------------------------------------------------------
      subroutine nekStab_torque(fname)
         include 'SIZE'
         include 'TOTAL'
         character(len=*), intent(in) :: fname
         character(len=256) :: fname_local  ! Local copy for ifx compatibility

         ! becuase of this block we can not use implicit none in this routine !
         real sij(lx1*ly1*lz1*6*lelv)
         real pm1(lx1,ly1,lz1,lelv)
         real xm0(lx1,ly1,lz1,lelt), ym0(lx1,ly1,lz1,lelt), zm0(lx1,ly1,lz1,lelt)
         integer, parameter :: lr = lx1*ly1*lz1
         real ur(lr), us(lr), ut(lr), vr(lr), vs(lr), vt(lr), wr(lr), ws(lr), wt(lr)
         real scale_vf(3)
         !!!

         common /scrns/ sij
         common /scrcg/ pm1
         common /scrsf/ xm0, ym0, zm0
         common /scruz/ ur, us, ut, vr, vs, vt, wr, ws, wt
         common /cvflow_r/ scale_vf

         ! Local variables only
         logical, save :: initialized
         data initialized/.false./

         integer, save :: bIDs(1), iobj_wall(1)
         real, save :: x0(3), scale
         data x0 /3*0.0d0/
         data scale /2.0d0/

         integer :: nv, i, ie, iobj, memtot, mem, ieg, ifc, nij
         real :: glmin, glmax, x1min, x2min, x3min, x1max, x2max, x3max, w1(0:maxobj)

         nv = nx1*ny1*nz1*nelv

         if (.not. initialized) then
            fname_local = fname  ! Copy to local variable for ifx compatibility
            if (nid == 0) write (6, *) 'Initializing torque routine... (nid =', nid, ')'
            if (nid == 0) write (6, *) 'About to open file ', trim(fname_local)
            if (nid == 0) open (737, file=trim(fname_local), action='write', status='replace')
            if (nid == 0) write (6, *) 'File opened. Setting bIDs(1) = 1'
            bIDs(1) = 1
            if (nid == 0) write (6, *) 'Calling create_obj: iobj_wall(1) [before] =', iobj_wall(1), ', bIDs(1) =', bIDs(1)
            call create_obj(iobj_wall(1), bIDs, 1)
            if (nid == 0) write (6, *) 'Returned from create_obj: iobj_wall(1) [after] =', iobj_wall(1)
            call cfill(vdiff, param(2), nv)  ! Fill up viscous array with default on initialization
            initialized = .true.
            if (nid == 0) write (6, *) 'Initialization block completed.'
         end if

         call mappr(pm1, pr, xm0, ym0) ! map pressure onto Mesh 1

         if (param(55) /= 0) then ! Add mean_pressure_gradient.X to p:
            dpdx_mean = -scale_vf(1)
            dpdy_mean = -scale_vf(2)
            dpdz_mean = -scale_vf(3)
         end if

         call add2s2(pm1, xm1, dpdx_mean, nv)  ! Doesn't work if object is cut by periodic boundary.
         call add2s2(pm1, ym1, dpdy_mean, nv)  ! In this case, set ._mean=0 and compensate in post.
         call add2s2(pm1, zm1, dpdz_mean, nv)

         nij = 3 ! Compute sij
         if (if3d .or. ifaxis) nij = 6
         call comp_sij(sij, nij, vx, vy, vz, ur, us, ut, vr, vs, vt, wr, ws, wt)

         call cadd2(xm0, xm1, -x0(1), nv)
         call cadd2(ym0, ym1, -x0(2), nv)
         call cadd2(zm0, zm1, -x0(3), nv)

         x1min = glmin(xm0(1,1,1,1), nv)
         x2min = glmin(ym0(1,1,1,1), nv)
         x3min = glmin(zm0(1,1,1,1), nv)
         x1max = glmax(xm0(1,1,1,1), nv)
         x2max = glmax(ym0(1,1,1,1), nv)
         x3max = glmax(zm0(1,1,1,1), nv)

         do i = 0, maxobj
            dragpx(i) = 0   ! Pressure drag x
            dragvx(i) = 0   ! Viscous drag x
            dragx(i)  = 0   ! Total drag x
            dragpy(i) = 0
            dragvy(i) = 0
            dragy(i)  = 0
            dragpz(i) = 0
            dragvz(i) = 0
            dragz(i)  = 0
            torqpx(i) = 0
            torqvx(i) = 0
            torqx(i)  = 0
            torqpy(i) = 0
            torqvy(i) = 0
            torqy(i)  = 0
            torqpz(i) = 0
            torqvz(i) = 0
            torqz(i)  = 0
         end do

         ifield = 1
         do iobj = 1, nobj
            memtot = nmember(iobj)
            do mem = 1, memtot
                ieg = object(iobj, mem, 1)
                ifc = object(iobj, mem, 2)
                if (gllnid(ieg) == nid) then ! This processor has a contribution
                    ie = gllel(ieg)
                    call drgtrq(dgtq, xm0, ym0, zm0, sij, pm1, vdiff, ifc, ie)
                    call cmult(dgtq, scale, 12)

                    dragpx(iobj) = dragpx(iobj) + dgtq(1,1)  ! Pressure
                    dragpy(iobj) = dragpy(iobj) + dgtq(2,1)
                    dragpz(iobj) = dragpz(iobj) + dgtq(3,1)

                    dragvx(iobj) = dragvx(iobj) + dgtq(1,2)  ! Viscous
                    dragvy(iobj) = dragvy(iobj) + dgtq(2,2)
                    dragvz(iobj) = dragvz(iobj) + dgtq(3,2)

                    torqpx(iobj) = torqpx(iobj) + dgtq(1,3)  ! Pressure
                    torqpy(iobj) = torqpy(iobj) + dgtq(2,3)
                    torqpz(iobj) = torqpz(iobj) + dgtq(3,3)

                    torqvx(iobj) = torqvx(iobj) + dgtq(1,4)  ! Viscous
                    torqvy(iobj) = torqvy(iobj) + dgtq(2,4)
                    torqvz(iobj) = torqvz(iobj) + dgtq(3,4)
                end if
            end do
         end do

         ! Sum contributions from all processors
         call gop(dragpx, w1, '+  ', maxobj+1)
         call gop(dragpy, w1, '+  ', maxobj+1)
         call gop(dragpz, w1, '+  ', maxobj+1)
         call gop(dragvx, w1, '+  ', maxobj+1)
         call gop(dragvy, w1, '+  ', maxobj+1)
         call gop(dragvz, w1, '+  ', maxobj+1)
         call gop(torqpx, w1, '+  ', maxobj+1)
         call gop(torqpy, w1, '+  ', maxobj+1)
         call gop(torqpz, w1, '+  ', maxobj+1)
         call gop(torqvx, w1, '+  ', maxobj+1)
         call gop(torqvy, w1, '+  ', maxobj+1)
         call gop(torqvz, w1, '+  ', maxobj+1)

         do i = 1, nobj ! Combine results
            dragx(i) = dragpx(i) + dragvx(i)
            dragy(i) = dragpy(i) + dragvy(i)
            dragz(i) = dragpz(i) + dragvz(i)

            torqx(i) = torqpx(i) + torqvx(i)
            torqy(i) = torqpy(i) + torqvy(i)
            torqz(i) = torqpz(i) + torqvz(i)

            dragpx(0) = dragpx(0) + dragpx(i)
            dragvx(0) = dragvx(0) + dragvx(i)
            dragx(0)  = dragx(0)  + dragx(i)

            dragpy(0) = dragpy(0) + dragpy(i)
            dragvy(0) = dragvy(0) + dragvy(i)
            dragy(0)  = dragy(0)  + dragy(i)

            dragpz(0) = dragpz(0) + dragpz(i)
            dragvz(0) = dragvz(0) + dragvz(i)
            dragz(0)  = dragz(0)  + dragz(i)

            torqpx(0) = torqpx(0) + torqpx(i)
            torqvx(0) = torqvx(0) + torqvx(i)
            torqx(0)  = torqx(0)  + torqx(i)

            torqpy(0) = torqpy(0) + torqpy(i)
            torqvy(0) = torqvy(0) + torqvy(i)
            torqy(0)  = torqy(0)  + torqy(i)

            torqpz(0) = torqpz(0) + torqpz(i)
            torqvz(0) = torqvz(0) + torqvz(i)
            torqz(0)  = torqz(0)  + torqz(i)
         end do

         ! Match output format to htps routine: use 1p in format string for floats
         do i = 1, nobj ! Write results
            if (nio == 0) then
               if (if3D .or. ifaxis) then ! 3D or axisymmetric
                  write (737, "(i8,1p19E15.7)") istep, time,
     $   dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), dragz(i), dragpz(i), dragvz(i),
     $   torqx(i), torqpx(i), torqvx(i), torqy(i), torqpy(i), torqvy(i), torqz(i), torqpz(i), torqvz(i)
               else ! 2D or not axisymmetric
               write (737, "(i8,1p10E15.7)") istep, time,
     $   dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), torqz(i), torqpz(i), torqvz(i)
              end if ! if3D .or. ifaxis
            end if ! nio == 0
         end do ! i = 1, nobj

         if (nio == 0) flush(737) ! Ensure data is written to disk

      end subroutine nekStab_torque
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! nekStab_define_obj — Define wall boundary objects for torque
      !
      ! Purpose:
      !   Assigns boundaryID=1 to all wall ('W') faces. Must be called
      !   from userbc or usrdat2 before nekStab_torque.
      !-----------------------------------------------------------------------
      subroutine nekStab_define_obj
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer :: iel, iface

         do iel=1,nelt
            do iface = 1, 2*ndim
               if (cbc(iface,iel,1) .eq. 'W  ') then
                  boundaryID(iface,iel) = 1
               endif
            enddo
         enddo

      end subroutine nekStab_define_obj
      !-----------------------------------------------------------------------

      end module nekstab_torque_mod
