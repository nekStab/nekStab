!-----------------------------------------------------------------------
! forcing.f90 -- Volume forcing and sponge zone routines
!
! Purpose:
!   Provides volume forcing callbacks for Nek5000 (velocity and
!   temperature), sponge zone initialization and evaluation, and
!   the smooth step function for sponge profiles.
!
! Public interface:
!   nekStab_forcing      -- velocity forcing callback (SFD/TDF/sponge)
!   nekStab_forcing_temp -- temperature forcing callback
!   activate_sponge      -- initialize sponge zone
!   spng_init            -- set sponge parameters and reference fields
!   spng_set             -- compute sponge spatial function
!   fringe_step          -- smooth step function for fringe profiles
!
! Dependencies:
!   nekstab_nek_bridge, krylov_subspace
!-----------------------------------------------------------------------

module nekstab_forcing_mod
   use nekstab_nek_bridge
   implicit none
   private
   public :: nekStab_forcing, nekStab_forcing_temp, &
      nekStab_set_buoyancy, nekStab_qvol, &
      activate_sponge, spng_init, spng_set, fringe_step
contains

!-----------------------------------------------------------------------
! nekStab_set_buoyancy -- Configure Boussinesq buoyancy forcing
!
! Purpose:
!   Enables the shared Boussinesq coupling ff_i += coeff*dir_i*T
!   and broadcasts the configuration to all MPI ranks.
!-----------------------------------------------------------------------
subroutine nekStab_set_buoyancy(dx, dy, dz, coeff)
   real, intent(in) :: dx, dy, dz, coeff

   ifbuoyancy = .true.
   buoyancy_qvol_wired = .false.
   buoyancy_dir(1) = dx
   buoyancy_dir(2) = dy
   buoyancy_dir(3) = dz
   thermal_buoyancy_coeff = coeff

   call bcast(ifbuoyancy, lsize)
   call bcast(buoyancy_qvol_wired, lsize)
   call bcast(buoyancy_dir, 3*wdsize)
   call bcast(thermal_buoyancy_coeff, wdsize)
end subroutine nekStab_set_buoyancy

!-----------------------------------------------------------------------
! nekStab_forcing -- Velocity forcing callback
!
! Purpose:
!   Adds nekStab volume forcing (SFD, TDF, sponge, OTD) at each
!   grid point. Called by Nek5000 during the solve.
!
! Arguments:
!   ffx, ffy, ffz [inout] -- force components, accumulated
!   ix, iy, iz    [in]    -- local grid indices
!   ieg            [in]    -- global element number
!-----------------------------------------------------------------------
subroutine nekStab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg)
   real, intent(inout) :: ffx, ffy, ffz
   integer, intent(in) :: ix, iy, iz, ieg
   integer :: iel, ip
   real :: tloc

   iel = gllel(ieg) ! local element number
   ! compute the corresponding index in the pertubation arrays
   ip = ix + nx1*(iy - 1 + ny1*(iz - 1 + nz1*(iel - 1)))

   ! here we add other forcings (SFD, TDF, if any)
   ffx = ffx + fcx(ix, iy, iz, iel)
   ffy = ffy + fcy(ix, iy, iz, iel)
   if (if3D) ffz = ffz + fcz(ix, iy, iz, iel)

   if (ifbuoyancy .and. ifheat) then
      if (jp == 0 .or. .not. ifadj) then
         if (jp == 0) then
            tloc = t(ix, iy, iz, iel, 1)
         else
            tloc = tp(ip, 1, jp)
         end if
         ffx = ffx + buoyancy_dir(1)*thermal_buoyancy_coeff*tloc
         ffy = ffy + buoyancy_dir(2)*thermal_buoyancy_coeff*tloc
         if (if3D) then
            ffz = ffz + buoyancy_dir(3)*thermal_buoyancy_coeff*tloc
         end if
      end if
   end if

   if (spng_st /= 0) then

      ! compute the corresponding index in the pertubation arrays
      ! computed once above for buoyancy and sponge forcing
      if (jp == 0) then ! dns
         ffx = ffx + spng_fn(ip)*(spng_vr(ip, 1) - vx(ix, iy, iz, iel))*spng_st
         ffy = ffy + spng_fn(ip)*(spng_vr(ip, 2) - vy(ix, iy, iz, iel))*spng_st
         if (if3D) ffz = ffz + spng_fn(ip)*(spng_vr(ip, ndim) - vz(ix, iy, iz, iel))*spng_st
      else ! perturbation
         ffx = ffx - spng_fn(ip)*vxp(ip, jp)*spng_st
         ffy = ffy - spng_fn(ip)*vyp(ip, jp)*spng_st
         if (if3D) ffz = ffz - spng_fn(ip)*vzp(ip, jp)*spng_st

      end if

   end if

   if (ifotd .and. jp /= 0) then
      ip = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(iel - 1)))
      ffx = ffx - otdfx(ip, jp)
      ffy = ffy - otdfy(ip, jp)
      if (if3D) ffz = ffz - otdfz(ip, jp)
   end if
end subroutine nekStab_forcing

!-----------------------------------------------------------------------
! nekStab_qvol -- Temperature forcing callback for adjoint buoyancy
!
! Purpose:
!   Adds the transpose of the Boussinesq momentum forcing to the
!   adjoint temperature equation.
!-----------------------------------------------------------------------
subroutine nekStab_qvol(qvol, ix, iy, iz, ieg)
   real, intent(inout) :: qvol
   integer, intent(in) :: ix, iy, iz, ieg
   integer :: iel, ip
   real :: vdot

   if (ifbuoyancy .and. ifheat .and. ifadj .and. jp > 0) then
      iel = gllel(ieg)
      ip = ix + nx1*(iy - 1 + ny1*(iz - 1 + nz1*(iel - 1)))
      vdot = buoyancy_dir(1)*vxp(ip, jp) + buoyancy_dir(2)*vyp(ip, jp)
      if (if3D) vdot = vdot + buoyancy_dir(3)*vzp(ip, jp)
      qvol = qvol + thermal_buoyancy_coeff*vdot
      buoyancy_qvol_wired = .true.
   end if
end subroutine nekStab_qvol

!-----------------------------------------------------------------------
! nekStab_forcing_temp -- Temperature forcing callback
!
! Purpose:
!   Adds nekStab temperature forcing (sponge) at each grid point.
!
! Arguments:
!   temp           [inout] -- temperature force, accumulated
!   ix, iy, iz     [in]    -- local grid indices
!   ieg            [in]    -- global element number
!   m              [in]    -- scalar field index
!-----------------------------------------------------------------------
subroutine nekStab_forcing_temp(temp, ix, iy, iz, ieg, m)
   real, intent(inout) :: temp
   integer, intent(in) :: ix, iy, iz, ieg
   integer, intent(in) :: m
   integer :: iel, ip

   ! local element number
   iel = gllel(ieg)
   ! compute the corresponding index in the pertubation arrays
   ip = ix + nx1*(iy - 1 + ny1*(iz - 1 + nz1*(iel - 1)))
   if (jp == 0) temp = temp + fct(ix, iy, iz, iel, m)

   if (spng_st /= 0) then

      ! compute the corresponding index in the pertubation arrays
      ! computed once above for temperature sponge forcing
      if (jp == 0) then ! dns ! t(1,1,1,1,ifield-1)
         temp = temp + spng_fn(ip)*(spng_vt(ip, m) - t(ix, iy, iz, iel, m))*spng_st
      else ! perturbation   ! tp(lpx1*lpy1*lpz1*lpelt,ldimt,lpert)
         temp = temp - spng_fn(ip)*tp(ip, m, jp)*spng_st
      end if

   end if
end subroutine nekStab_forcing_temp

!-----------------------------------------------------------------------
! activate_sponge -- Initialize sponge zone and mask mass matrix
!
! Purpose:
!   Sets up the sponge zone function and zeros the mass matrix
!   inside the sponge region so eigensolver ignores that region.
!-----------------------------------------------------------------------
subroutine activate_sponge
   use krylov_subspace
   integer i
   nv = nx1*ny1*nz1*nelv
   if (spng_st /= 0) then !sponge on

      if (nid == 0) write (6, *)
      if (nid == 0) write (6, *) ' Initializing sponge...'
      if (nid == 0) write (6, *) ' Sponge strenght:', spng_st
      if (spng_st < 0) then
         spng_st = abs(spng_st)
         if (nid == 0) write (6, *) ' Ensure positive sponge strenght:', spng_st
      end if
      call spng_init

   !     applying sponge function to BM1 matrix to remove the sponge zone from eigensolver
      do i = 1, nv
         if (spng_fn(i) /= 0) bm1s(i, 1, 1, 1) = 0.0d0
      end do

   !  outposting BM1s to disk for check
   !  ifto_sav = ifto; ifpo_sav = ifpo
   !  ifvo=.false.; ifpo = .false.; ifto = .true.
   !  call outpost(vx,vy,vz,pr,bm1s,'BMS')
   !  ifvo=.true.; ifpo = ifpo_sav; ifto = ifto_sav

      if (nid == 0) write (6, *) 'Sponge activated.'
      if (nid == 0) write (6, *)
   end if
end subroutine activate_sponge

!-----------------------------------------------------------------------
! spng_init -- Compute sponge widths and save reference fields
!-----------------------------------------------------------------------
subroutine spng_init

   integer :: nv
   nv = nx1*ny1*nz1*nelv
   acc_spg = abs(acc_spg)

   spng_wl(1) = (1.0d0 - acc_spg)*xLspg ! Left flat/outer fringe width in x
   spng_wl(2) = (1.0d0 - acc_spg)*yLspg
   if (if3D) spng_wl(3) = (1.0d0 - acc_spg)*zLspg

   spng_wr(1) = (1.0d0 - acc_spg)*xRspg ! Right flat/outer fringe width in x
   spng_wr(2) = (1.0d0 - acc_spg)*yRspg
   if (if3D) spng_wr(3) = (1.0d0 - acc_spg)*zRspg

   spng_dl(1) = (acc_spg)*xLspg ! Left smooth-ramp fringe width in x
   spng_dl(2) = (acc_spg)*yLspg
   if (if3D) spng_dl(3) = (acc_spg)*zLspg

   spng_dr(1) = (acc_spg)*xRspg ! Right smooth-ramp fringe width in x
   spng_dr(2) = (acc_spg)*yRspg
   if (if3D) spng_dr(3) = (acc_spg)*zRspg

   if (nid == 0) then
      write (6, *) '  Left  section width x y z:', spng_wl
      write (6, *) '  Right section width x y z:', spng_wr
      write (6, *) '  Left  drop/rise width x y z:', spng_dl
      write (6, *) '  Right drop/rise width x y z:', spng_dr
   end if

   !     save reference field -> sponge value reference
   call opcopy(spng_vr(1, 1), spng_vr(1, 2), spng_vr(1, NDIM), vx, vy, vz) !only DNS
   if (ifto) call copy(spng_vt(:, 1), t(1, 1, 1, 1, 1), nv) !only DNS - temperature
   call spng_set ! -> compute spng_fn

end subroutine spng_init

!-----------------------------------------------------------------------
! spng_set -- compute the spatial fringe (sponge) mask spng_fn on the mesh
!
! Implementation based on the fringe-region mask described by
! Nordstrom, Nordin & Henningson (1999) and Lundbladh et al. (1999).
!-----------------------------------------------------------------------
subroutine spng_set
   integer :: nv, idir
   real :: bmin(ldim), bmax(ldim)
   real :: xxmin, xxmax, xxmin_c, xxmax_c
   logical :: ltmp, ltmp2
   real lcoord(lx1*ly1*lz1*lelv)
   common /SCRUZ/ lcoord

   nv = nx1*ny1*nz1*nelv

   call rzero(spng_fn, nv)

   bmin(1) = xmn
   bmax(1) = xmx
   bmin(2) = ymn
   bmax(2) = ymx
   if (if3D) then
      bmin(3) = zmn
      bmax(3) = zmx
   end if

   do idir = 1, ndim
      if (spng_wl(idir) > 0.0d0 .or. spng_wr(idir) > 0.0d0) then
         if (spng_wl(idir) < spng_dl(idir) .or. spng_wr(idir) < spng_dr(idir)) then
            if (nid == 0) write (6, *) 'Wrong sponge parameters!'
         end if

         xxmax = bmax(idir) - spng_wr(idir)
         xxmin = bmin(idir) + spng_wl(idir)
         xxmax_c = xxmax + spng_dr(idir)
         xxmin_c = xxmin - spng_dl(idir)

         if (xxmax <= xxmin) then
            if (nid == 0) write (6, *) 'Sponge too wide'
         else
            call load_fringe_coordinate(idir, lcoord, nv)
            call accumulate_fringe_dimension(idir, lcoord, nv, xxmin_c, xxmin, xxmax, xxmax_c)
         end if
      end if
   end do

      ! SPG diagnostic dump: write the mask as the scalar field with no
      ! pressure. outpost2 takes nfldt before the 3-char name; the ifto/
      ! ifpo toggle is the pinned SPG-outpost contract (sponge-spec sec 2).
   ltmp = ifto; ltmp2 = ifpo
   ifto = .true.; ifpo = .false.
   call outpost2(spng_vr(1, 1), spng_vr(1, 2), spng_vr(1, ndim), spng_fn, spng_fn, 1, 'SPG')
   ifto = ltmp; ifpo = ltmp2

contains

      ! Coordinate selection for the Nordstrom, Nordin & Henningson (1999)
      ! and Lundbladh et al. (1999) fringe-region mask pass.
   subroutine load_fringe_coordinate(idir_in, coord, ncoord)
      integer, intent(in) :: idir_in, ncoord
      real, intent(out) :: coord(ncoord)

      select case (idir_in)
      case (1)
         call copy(coord, xm1, ncoord)
      case (2)
         call copy(coord, ym1, ncoord)
      case (3)
         call copy(coord, zm1, ncoord)
      end select
   end subroutine load_fringe_coordinate

      ! Max-accumulation pass for the Nordstrom, Nordin & Henningson (1999)
      ! and Lundbladh et al. (1999) one-dimensional fringe ramps.
   subroutine accumulate_fringe_dimension(idir_in, coord, ncoord, left_outer, left_inner, right_inner, right_outer)
      integer, intent(in) :: idir_in, ncoord
      real, intent(in) :: coord(ncoord)
      real, intent(in) :: left_outer, left_inner, right_inner, right_outer
      integer :: ipt
      real :: lambda

      do ipt = 1, ncoord
         lambda = fringe_profile_value(coord(ipt), spng_wl(idir_in), spng_wr(idir_in), &
            left_outer, left_inner, right_inner, right_outer)
         spng_fn(ipt) = max(spng_fn(ipt), lambda)
      end do
   end subroutine accumulate_fringe_dimension

      ! Pure one-dimensional ramp for the Nordstrom, Nordin & Henningson
      ! (1999) and Lundbladh et al. (1999) fringe-region profile.
   pure real function fringe_profile_value(r, left_width, right_width, left_outer, left_inner, right_inner, right_outer)
      real, intent(in) :: r, left_outer, left_inner, right_inner, right_outer
      real, intent(in) :: left_width, right_width

      if (r <= left_outer) then
         fringe_profile_value = 1.0d0
      else if (r < left_inner) then
         fringe_profile_value = fringe_step((left_inner - r)/left_width)
      else if (r <= right_inner) then
         fringe_profile_value = 0.0d0
      else if (r < right_outer) then
         fringe_profile_value = fringe_step((r - right_inner)/right_width)
      else
         fringe_profile_value = 1.0d0
      end if
   end function fringe_profile_value

end subroutine spng_set

!-----------------------------------------------------------------------
! fringe_step -- smooth step ramp primitive for fringe-region profiles
!
! Implementation based on the Nordstrom, Nordin & Henningson (1999)
! and Lundbladh et al. (1999) fringe-region smooth-step specification.
!-----------------------------------------------------------------------
pure real function fringe_step(x)
   real, intent(in) :: x

   if (x <= 0.0010d0) then
      fringe_step = 0.0d0
   else if (x <= 0.9990d0) then
      fringe_step = 1.0d0/(1.0d0 + exp(1.0d0/(x - 1.0d0) + 1.0d0/x))
   else
      fringe_step = 1.0d0
   end if
end function fringe_step

end module nekstab_forcing_mod

subroutine nekStab_set_buoyancy(dx, dy, dz, coeff)
   use nekstab_forcing_mod, &
      only: mod_nekStab_set_buoyancy => nekStab_set_buoyancy
   implicit none
   real, intent(in) :: dx, dy, dz, coeff

   call mod_nekStab_set_buoyancy(dx, dy, dz, coeff)
end subroutine nekStab_set_buoyancy

subroutine nekStab_qvol(qvol, ix, iy, iz, ieg)
   use nekstab_forcing_mod, only: mod_nekStab_qvol => nekStab_qvol
   implicit none
   real, intent(inout) :: qvol
   integer, intent(in) :: ix, iy, iz, ieg

   call mod_nekStab_qvol(qvol, ix, iy, iz, ieg)
end subroutine nekStab_qvol
