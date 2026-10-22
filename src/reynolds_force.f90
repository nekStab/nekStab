!-----------------------------------------------------------------------
! reynolds_force.f90 — Resolved Reynolds stress and its divergence force
!
! Purpose:
!   Forms the resolved Reynolds stress R_ij = E(u_i u_j) - E(u_i)E(u_j)
!   from the moments accumulated by avg_all (the /chkavg/ means and the
!   /chkrms/ second and cross moments), computes the momentum force
!   f = -div(R) that makes the mean a steady solution of the Navier-
!   Stokes equations, writes both to disk, and applies a previously
!   written force as a frozen volume force in userf. This is the
!   fluctuation correlation, not the k-tau RANS closure.
!
!   Steady mean equation: NS(U) = -div(R). Nek userf adds f on the
!   right-hand side, so f = -div(R) with (div R)_i = dR_ij/dx_j.
!
!   The force is off unless the case calls reynolds_enable from
!   nekStab_usrchk. nekStab_init then calls reynolds_load once, which
!   parks and restores the restart state and arms the stored field for
!   the rest of the run. nekStab_forcing applies the force only when
!   jp == 0, so it is present in every finite-difference residual
!   evaluation and cancels in the difference, and an analytic
!   perturbation equation does not gain a constant source.
!
! Public interface:
!   reynolds_commit                  — form R and f from the live averages
!                                      and write RS1/RS2/FRS (does not arm)
!   reynolds_load                    — load FRS<session>0.f00001 and arm
!                                      (called once by nekStab_init
!                                      when enabled; not for case use)
!   reynolds_add(ffx,ffy,ffz,ix,iy,iz,iel) — add the armed force at a point
!   reynolds_enable                  — arm the init-time load (call from
!                                      nekStab_usrchk; not a mode)
!   reynolds_enabled()               — .true. once reynolds_enable ran
!   reynolds_armed()                 — .true. once a force file is loaded
!
! Dependencies:
!   nekstab_nek_bridge (sizes, nelv, if3d, ifield, vmult, SESSION, vx,
!   vy, vz, t, param, time)
!-----------------------------------------------------------------------
module nekstab_reynolds

   use nekstab_nek_bridge, only: nekStab_log, nekStab_error, &
                                 ax1, ay1, az1, ax2, ay2, az2, ldimt, &
                                 lx1, ly1, lz1, lx2, ly2, lz2, lelt, lelv, &
                                 nelv, nelt, nid, if3d, lsize, ifield, &
                                 vmult, SESSION, vx, vy, vz, pr, t, param, &
                                 time, ifvo, ifto, ifpo, ifxyo, ifpsco, ldimt1, wdsize
   implicit none

   private

   ! Compile-time storage extent — matches the (lx1,ly1,lz1,lelv) layout
   ! of the velocity field so the linear index in reynolds_add is
   ! stride-consistent with base_mut_get. lelv is the compile-time cap
   ! on both nelv and nelt.
   integer, parameter :: lrs = lx1*ly1*lz1*lelv

   real, save :: frs_x(lrs), frs_y(lrs), frs_z(lrs)  ! stored force

   logical, save :: enabled = .false.   ! set by reynolds_enable
   logical, save :: armed = .false.

   public :: reynolds_commit, reynolds_load, reynolds_add, &
             reynolds_enable, reynolds_enabled, reynolds_armed

contains

!-----------------------------------------------------------------------
! reynolds_commit — Form R and f = -div(R) from the live averages
!
! Reads /chkavg/ and /chkrms/ exactly as statistics.f90 does, forms
! R_ij = E(u_i u_j) - E(u_i)E(u_j), differentiates each component with
! rs_grad, and accumulates f = -div(R). Writes three velocity-layout
! files and does NOT arm the force: the march that produced the
! averages is unchanged. The R, derivative, and force fields are all
! local allocatable scratch; the saved frs_* storage is written only
! by reynolds_load, so commit never disturbs an armed force.
!
! Output (single precision, param(63)=0 because 64-bit avg files from
! this case do not reload):
!   RS1 — x=Ruu, y=Rvv, z=Rww
!   RS2 — x=Ruv, y=Rvw, z=Rwu
!   FRS — x=fx, y=fy, z=fz
!-----------------------------------------------------------------------
   subroutine reynolds_commit

      integer :: ntot, nout, i, idx, ix, iy, iz, iel, istat
      real :: p63_sav
      logical :: ifvo_sav, ifto_sav, ifpo_sav, ifxyo_sav
      logical :: ifpsco_sav(ldimt1)
      real, allocatable :: r_uu(:), r_vv(:), r_ww(:)
      real, allocatable :: r_uv(:), r_vw(:), r_wu(:)
      real, allocatable :: wk1(:), wk2(:), wk3(:)
      real, allocatable :: f_x(:), f_y(:), f_z(:)

      ! avg_all common blocks (not in TOTAL) — same layout as statistics.f90
      real :: uavg(ax1, ay1, az1, lelt), vavg(ax1, ay1, az1, lelt), &
              wavg(ax1, ay1, az1, lelt), tavg(ax1, ay1, az1, lelt, ldimt), &
              pavg(ax2, ay2, az2, lelt)
      common /chkavg/ uavg, vavg, wavg, tavg, pavg

      real :: urms(ax1, ay1, az1, lelt), vrms(ax1, ay1, az1, lelt), &
              wrms(ax1, ay1, az1, lelt), trms(ax1, ay1, az1, lelt, ldimt), &
              prms(ax2, ay2, az2, lelt), vwms(ax1, ay1, az1, lelt), &
              wums(ax1, ay1, az1, lelt), uvms(ax1, ay1, az1, lelt)
      common /chkrms/ urms, vrms, wrms, trms, prms, vwms, wums, uvms

!     Same SIZE check as statistics.f90: a SIZE built with ax1 = 1
!     would make the common-block dummy arrays index off the real
!     averages.
      if (ax1 /= lx1 .or. ay1 /= ly1 .or. az1 /= lz1) then
         call nekStab_error('ABORT: wrong size of ax1,ay1,az1 in avg_all(), check SIZE!')
      end if
      if (ax2 /= lx2 .or. ay2 /= ly2 .or. az2 /= lz2) then
         call nekStab_error('ABORT: wrong size of ax2,ay2,az2 in avg_all(), check SIZE!')
      end if

!     lelv is the compile-time cap of the velocity field rs_grad
!     differentiates. nelt is rank-local and legitimately exceeds lelv
!     in conjugate heat transfer, so only nelv can abort collectively;
!     the scratch below is sized nout = nxyz*nelt to cover outpost2.
      if (nelv > lelv) then
         call nekStab_error('reynolds: element count exceeds lelv')
         return
      end if

      ntot = lx1*ly1*lz1*nelv
      nout = lx1*ly1*lz1*nelt

      allocate (r_uu(nout), r_vv(nout), r_ww(nout), &
                r_uv(nout), r_vw(nout), r_wu(nout), &
                wk1(nout), wk2(nout), wk3(nout), &
                f_x(nout), f_y(nout), f_z(nout), stat=istat)
      if (istat /= 0) then
!        A failed allocate leaves its completed members allocated.
         if (allocated(r_uu)) deallocate (r_uu)
         if (allocated(r_vv)) deallocate (r_vv)
         if (allocated(r_ww)) deallocate (r_ww)
         if (allocated(r_uv)) deallocate (r_uv)
         if (allocated(r_vw)) deallocate (r_vw)
         if (allocated(r_wu)) deallocate (r_wu)
         if (allocated(wk1)) deallocate (wk1)
         if (allocated(wk2)) deallocate (wk2)
         if (allocated(wk3)) deallocate (wk3)
         if (allocated(f_x)) deallocate (f_x)
         if (allocated(f_y)) deallocate (f_y)
         if (allocated(f_z)) deallocate (f_z)
         call nekStab_error('reynolds: scratch allocation failed')
         return
      end if

!     R_ij = E(u_i u_j) - E(u_i)E(u_j). urms..wrms are E(u^2) (not
!     variances) and uvms,vwms,wums are E(uv), E(vw), E(wu).
      idx = 0
      do iel = 1, nelv
         do iz = 1, lz1
            do iy = 1, ly1
               do ix = 1, lx1
                  idx = idx + 1
                  r_uu(idx) = urms(ix, iy, iz, iel) &
                              - uavg(ix, iy, iz, iel)*uavg(ix, iy, iz, iel)
                  r_vv(idx) = vrms(ix, iy, iz, iel) &
                              - vavg(ix, iy, iz, iel)*vavg(ix, iy, iz, iel)
                  r_ww(idx) = wrms(ix, iy, iz, iel) &
                              - wavg(ix, iy, iz, iel)*wavg(ix, iy, iz, iel)
                  r_uv(idx) = uvms(ix, iy, iz, iel) &
                              - uavg(ix, iy, iz, iel)*vavg(ix, iy, iz, iel)
                  r_vw(idx) = vwms(ix, iy, iz, iel) &
                              - vavg(ix, iy, iz, iel)*wavg(ix, iy, iz, iel)
                  r_wu(idx) = wums(ix, iy, iz, iel) &
                              - wavg(ix, iy, iz, iel)*uavg(ix, iy, iz, iel)
               end do
            end do
         end do
      end do

!     f_i = -dR_ij/dx_j. rs_grad gives the strong-form SEM derivative,
!     single-valued on the velocity mesh. The z derivatives are skipped
!     when .not. if3d.
      f_x = 0.0d0
      f_y = 0.0d0
      f_z = 0.0d0

      call rs_grad(r_uu, wk1, wk2, wk3)
      do i = 1, ntot
         f_x(i) = f_x(i) - wk1(i)  ! -dRuu/dx
      end do

      call rs_grad(r_uv, wk1, wk2, wk3)
      do i = 1, ntot
         f_x(i) = f_x(i) - wk2(i)  ! -dRuv/dy
         f_y(i) = f_y(i) - wk1(i)  ! -dRuv/dx
      end do

      if (if3d) then
         call rs_grad(r_wu, wk1, wk2, wk3)
         do i = 1, ntot
            f_x(i) = f_x(i) - wk3(i)  ! -dRwu/dz
            f_z(i) = f_z(i) - wk1(i)  ! -dRwu/dx
         end do
      end if

      call rs_grad(r_vv, wk1, wk2, wk3)
      do i = 1, ntot
         f_y(i) = f_y(i) - wk2(i)  ! -dRvv/dy
      end do

      if (if3d) then
         call rs_grad(r_vw, wk1, wk2, wk3)
         do i = 1, ntot
            f_y(i) = f_y(i) - wk3(i)  ! -dRvw/dz
            f_z(i) = f_z(i) - wk2(i)  ! -dRvw/dy
         end do

         call rs_grad(r_ww, wk1, wk2, wk3)
         do i = 1, ntot
            f_z(i) = f_z(i) - wk3(i)  ! -dRww/dz
         end do
      end if

!     outpost2 copies nxyz*nelt words from these buffers, so zero the
!     elements past nelv when nelt > nelv.
      do i = ntot + 1, nout
         r_uu(i) = 0.0d0
         r_vv(i) = 0.0d0
         r_ww(i) = 0.0d0
         r_uv(i) = 0.0d0
         r_vw(i) = 0.0d0
         r_wu(i) = 0.0d0
         f_x(i) = 0.0d0
         f_y(i) = 0.0d0
         f_z(i) = 0.0d0
      end do

!     Single-precision velocity-layout files; outpost is collective.
!     ifxyo and every ifpsco off: an X block would make load_fld
!     overwrite xm1,ym1,zm1 on reload, and scalar planes would leak
!     into the FRS file.
      p63_sav = param(63)
      ifvo_sav = ifvo
      ifto_sav = ifto
      ifpo_sav = ifpo
      ifxyo_sav = ifxyo
      ifpsco_sav = ifpsco
      param(63) = 0.0d0 ! 64-bit avg files from this case do not reload
      call bcast(param(63), wdsize)
      ifvo = .true.
      ifto = .false.
      ifpo = .false.
      ifxyo = .false.
      ifpsco(:) = .false.

      call outpost(r_uu, r_vv, r_ww, pr, t, 'RS1')
      call outpost(r_uv, r_vw, r_wu, pr, t, 'RS2')
      call outpost(f_x, f_y, f_z, pr, t, 'FRS')

      param(63) = p63_sav
      call bcast(param(63), wdsize)
      ifvo = ifvo_sav
      ifto = ifto_sav
      ifpo = ifpo_sav
      ifxyo = ifxyo_sav
      ifpsco(:) = ifpsco_sav(:)

      deallocate (r_uu, r_vv, r_ww, r_uv, r_vw, r_wu, wk1, wk2, wk3, &
                  f_x, f_y, f_z)

      call nekStab_log('reynolds: wrote RS1, RS2, FRS (R_ij and f = -div R)')

   end subroutine reynolds_commit

!-----------------------------------------------------------------------
! rs_grad — Strong-form SEM gradient, single-valued on the velocity mesh
!
! gradm11 differentiates one element; its ux,uy,uz arguments are
! lx1*ly1*lz1 scratch, so each call fills one element of the output
! field and the loop runs over e = 1..nelv only. gradm1 would walk
! nelt elements and read past a velocity-length field when nelt > nelv.
! dssum + vmult under ifield = 1 is dsavg's velocity branch without its
! ifflow switch, which could move the averaging to the temperature
! mesh. In 2D gradm11 does not write uz, so gz is never read.
!-----------------------------------------------------------------------
   subroutine rs_grad(fld, gx, gy, gz)
      real, intent(in) :: fld(*)
      real, intent(inout) :: gx(*), gy(*), gz(*)
      integer :: e, nxyz, ifield_sav

      nxyz = lx1*ly1*lz1
      do e = 1, nelv
         call gradm11(gx((e - 1)*nxyz + 1), gy((e - 1)*nxyz + 1), &
                      gz((e - 1)*nxyz + 1), fld, e)
      end do

      ifield_sav = ifield
      ifield = 1  ! velocity-mesh gather-scatter handle
      call dssum(gx, lx1, ly1, lz1)
      call col2(gx, vmult, lx1*ly1*lz1*nelv)
      call dssum(gy, lx1, ly1, lz1)
      call col2(gy, vmult, lx1*ly1*lz1*nelv)
      if (if3d) then
         call dssum(gz, lx1, ly1, lz1)
         call col2(gz, vmult, lx1*ly1*lz1*nelv)
      end if
      ifield = ifield_sav

   end subroutine rs_grad

!-----------------------------------------------------------------------
! reynolds_load — Load FRS<session>0.f00001 and arm the force
!
! outpost starts each run's per-prefix counter at 1, so one
! reynolds_commit per run writes f00001; a second commit in the same
! run writes f00002 and is not loaded, and an older f00002 left by a
! previous run never wins. Only rank 0 inquires the one name; the
! logical is broadcast so all ranks agree before the collective
! load_fld. Missing file: log and leave
! the force unarmed (no abort, no allocation). load_fld fills
! vx,vy,vz,pr,t and time from the file, so the loaded state is parked
! and restored after the velocity is copied into the force storage.
! The t planes are parked one scalar at a time: t is (..,lelt,ldimt),
! and a contiguous nxyz*nelt*ldimt copy would cross the lelt padding
! and corrupt the later scalars.
!-----------------------------------------------------------------------
   subroutine reynolds_load
      character(len=80) :: filename
      logical :: file_found
      integer :: ntot, ntot2, nplan, m
      real :: time_sav
      real, allocatable :: vx_sav(:), vy_sav(:), vz_sav(:)
      real, allocatable :: pr_sav(:), t_sav(:, :)

      if (armed) return

      ntot = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv
      nplan = lx1*ly1*lz1*nelt

!     One commit per run writes FRS<session>0.f00001. Rank 0 inquires
!     that name; every rank uses the broadcast result before the
!     collective load_fld.
      filename = 'FRS'//trim(SESSION)//'0.f00001'
      file_found = .false.
      if (nid == 0) inquire (file=filename, exist=file_found)
      call bcast(file_found, lsize)

      if (.not. file_found) then
         call nekStab_log('reynolds: no '//trim(filename)// &
                          ' file, force not armed')
         return
      end if

      call nekStab_log('reynolds: loading frozen force '//trim(filename))

!     load_fld fills vx,vy,vz,pr,t and time from the file header, so the
!     whole loaded state is parked and restored around the call.
      time_sav = time
      allocate (vx_sav(ntot), vy_sav(ntot), vz_sav(ntot))
      allocate (pr_sav(ntot2), t_sav(nplan, ldimt))
      call copy(vx_sav, vx, ntot)
      call copy(vy_sav, vy, ntot)
      call copy(vz_sav, vz, ntot)
      call copy(pr_sav, pr, ntot2)
      do m = 1, ldimt
         call copy(t_sav(1, m), t(1, 1, 1, 1, m), nplan)
      end do

      call load_fld(filename)

      call copy(frs_x, vx, ntot)
      call copy(frs_y, vy, ntot)
      call copy(frs_z, vz, ntot)

      call copy(vx, vx_sav, ntot)
      call copy(vy, vy_sav, ntot)
      call copy(vz, vz_sav, ntot)
      call copy(pr, pr_sav, ntot2)
      do m = 1, ldimt
         call copy(t(1, 1, 1, 1, m), t_sav(1, m), nplan)
      end do
      time = time_sav
      call bcast(time, wdsize)

      deallocate (vx_sav, vy_sav, vz_sav, pr_sav, t_sav)

      armed = .true.

   end subroutine reynolds_load

!-----------------------------------------------------------------------
! reynolds_add — Add the stored force at a grid point
!
! No-op unless reynolds_load armed the force. Index uses the compile-
! time lx1,ly1,lz1 strides of the (lx1,ly1,lz1,lelv) storage layout.
!-----------------------------------------------------------------------
   subroutine reynolds_add(ffx, ffy, ffz, ix, iy, iz, iel)
      real, intent(inout) :: ffx, ffy, ffz
      integer, intent(in) :: ix, iy, iz, iel
      integer :: idx

      if (.not. armed) return

      idx = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(iel - 1)))
      ffx = ffx + frs_x(idx)
      ffy = ffy + frs_y(idx)
      if (if3d) ffz = ffz + frs_z(idx)

   end subroutine reynolds_add

!-----------------------------------------------------------------------
! reynolds_enable / reynolds_enabled — Opt-in flag for the init-time load
!
! A case calls reynolds_enable from nekStab_usrchk on every rank.
! nekStab_init reads reynolds_enabled after usrchk and calls
! reynolds_load once. The bcast keeps every rank's flag in step.
!-----------------------------------------------------------------------
   subroutine reynolds_enable
      enabled = .true.
      call bcast(enabled, lsize)
   end subroutine reynolds_enable

   logical function reynolds_enabled()
      reynolds_enabled = enabled
   end function reynolds_enabled

   logical function reynolds_armed()
      reynolds_armed = armed
   end function reynolds_armed

end module nekstab_reynolds
