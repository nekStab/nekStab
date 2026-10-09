!-----------------------------------------------------------------------
! diagnostics.f90 — Runtime diagnostics and monitoring
!
! Purpose:
!   Provides energy/enstrophy monitoring, vorticity output,
!   gradient norms, field smoothing, verbose timestep commentary,
!   and parameter printing for nekStab simulations.
!
! Public interface:
!   outpost_vort          — output vorticity field
!   norm_grad             — compute H1 gradient semi-norm
!   smooth_field          — smooth field via DSS + filtering
!   nekStab_outpost       — custom outpost with vorticity/vortex
!   nekStab_comment       — verbose timestep progress output
!   nekStab_printNEKParams — print Nek5000 parameters for sanity check
!   nekStab_energy        — write kinetic energy to file
!   nekStab_enstrophy     — write enstrophy to file
!
! Dependencies:
!   nekstab_nek_bridge, krylov_subspace, nekstab_vortex
!-----------------------------------------------------------------------

module nekstab_diagnostics
   use nekstab_nek_bridge, only: nekStab_error, nekStab_log, nid, nio, &
                                 ifvor, ifvox, ifoutfld, ifvo, ifaxis, if3D, &
                                 ifto, ifpo, ifpsco, iftran, ifheat, &
                                 courno, istep, lastep, time, param, uparam, &
                                 lx1, ly1, lz1, lelt, nelv, ldimt, ldimt1, &
                                 vx, vy, vz, pr, t, bm1, bm1s, binvm1, ifield, &
                                 volvm1, ym1, &
                                 xmn, xmx, ymn, ymx, zmn, zmx, &
                                 NEKSTAB_UNIT_ENERGY, NEKSTAB_UNIT_ENSTRO, &
                                 initc, spng_st, isize, &
                                 isDirect, isAdjoint, isTransientGrowth, &
                                 isFloquetDirect, isFloquetAdjoint, &
                                 isFloquetTransientGrowth, ifprojfld, &
                                 isNewtonFP, isNewtonPO, isNewtonPO_T, &
                                 session, lhis
   use nekstab_vortex, only: vortex_core
   implicit none
   private
   public :: outpost_vort, norm_grad, smooth_field, &
             nekStab_outpost, nekStab_comment, &
             nekStab_printNEKParams, nekStab_energy, &
             nekStab_enstrophy, nekStab_error, nekStab_log, &
             nekStab_check_setup
contains

!-----------------------------------------------------------------------
! outpost_vort — Compute and output vorticity field
!
! Purpose:
!   Computes the 3D vorticity from a velocity field and outputs
!   it using Nek5000's outpost routine. Only active if ifvor=.true.
!
! Arguments:
!   ux   [in] — x-velocity field
!   uy   [in] — y-velocity field
!   uz   [in] — z-velocity field
!   name [in] — 3-character output file prefix
!-----------------------------------------------------------------------
   subroutine outpost_vort(ux, uy, uz, name)
      use nekstab_krylov_subspace, only: lv
      implicit none
      real, intent(in) :: ux(lv), uy(lv), uz(lv)
      character(len=3), intent(in) :: name
      real wo1(lv), wo2(lv), vort(lv, 3)
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

   end subroutine outpost_vort

!-----------------------------------------------------------------------
! norm_grad — Compute H1 gradient semi-norm of a velocity field
!
! Purpose:
!   Computes the L2 norm of velocity gradients (H1 semi-norm),
!   weighted by the mass matrix. Used for convergence monitoring.
!
! Arguments:
!   vx_  [in]  — x-velocity field
!   vy_  [in]  — y-velocity field
!   vz_  [in]  — z-velocity field
!   pr_  [in]  — pressure field (unused, kept for interface)
!   t_   [in]  — temperature field (unused, kept for interface)
!   norma [out] — computed gradient norm
!-----------------------------------------------------------------------
   subroutine norm_grad(vx_, vy_, vz_, pr_, t_, norma)

      use nekstab_krylov_subspace, only: lv, lp, lt, nv
      implicit none

      real, intent(in), dimension(lv) :: vx_, vy_, vz_
      real, intent(in), dimension(lp) :: pr_
      real, intent(in), dimension(lt, ldimt) :: t_
      real, intent(out) :: norma

      real, dimension(lv) :: dudx, dudy, dudz
      real, dimension(lv) :: dvdx, dvdy, dvdz
      real, dimension(lv) :: dwdx, dwdy, dwdz

      real :: glsc3

      ! gradient computation
      call gradm1(dudx, dudy, dudz, vx_, nelv)
      call gradm1(dvdx, dvdy, dvdz, vy_, nelv)
      if (if3D) call gradm1(dwdx, dwdy, dwdz, vz_, nelv)

      norma = 0.0d0

      norma = norma + glsc3(dudx, dudx, bm1s, nv) + glsc3(dudy, dudy, bm1s, nv)
      norma = norma + glsc3(dvdx, dvdx, bm1s, nv) + glsc3(dvdy, dvdy, bm1s, nv)

      if (if3D) then
         norma = norma + glsc3(dudz, dudz, bm1s, nv)
         norma = norma + glsc3(dvdz, dvdz, bm1s, nv)
         norma = norma + glsc3(dwdx, dwdx, bm1s, nv)
         norma = norma + glsc3(dwdy, dwdy, bm1s, nv)
         norma = norma + glsc3(dwdz, dwdz, bm1s, nv)
      end if

   end subroutine norm_grad

!-----------------------------------------------------------------------
! smooth_field — Smooth a field via DSS averaging and filtering
!
! Purpose:
!   Applies direct stiffness summation averaging and spectral
!   filtering to produce a continuous, smooth field in H1.
!
! Arguments:
!   u [inout] — field to smooth
!-----------------------------------------------------------------------
   subroutine smooth_field(u)
      use nekstab_krylov_subspace, only: nv
      implicit none
      real, intent(inout), dimension(lx1, ly1, lz1, lelt) :: u
      integer :: ifld

      ! avg boundary - used in comp_vort3
      call col2(u, bm1, nv)
      call dssum(u, lx1, ly1, lz1)
      call col2(u, binvm1, nv)

      ! ensure continuous field that is in H1
      call dsavg(u) ! direct stiffness avg of u

      ! filter field
      ifld = ifield
      ifield = 1
      !weight, ncut, name
      call filter_s0(u, 0.5, 1, 'field')
      ifield = ifld

   end subroutine smooth_field

!-----------------------------------------------------------------------
! nekStab_outpost — Custom output with vorticity and vortex criteria
!
! Purpose:
!   Computes and outputs vorticity (vor*) and/or vortex criterion
!   fields (vox*, omR*) during DNS. Activated by ifvor and ifvox
!   flags respectively.
!-----------------------------------------------------------------------
   subroutine nekStab_outpost
      use nekstab_krylov_subspace, only: lv
      implicit none

      real vort(lv, 3), wo1(lv), wo2(lv)
      common/ugrad/vort, wo1, wo2

      logical ifto_sav, ifpo_sav, ifpsco_sav(ldimt1)

      if (ifoutfld .and. (ifvor .or. ifvox)) then

         ifto_sav = ifto; ifpo_sav = ifpo; ifpsco_sav = ifpsco
         ifto = .false.; ifpo = .false.; ifpsco(:) = .false.

         !---  > Compute and oupost vorticity.
         if (ifvor) then
            call oprzero(vort(1, 1), vort(1, 2), vort(1, 3))
            call comp_vort3(vort, wo1, wo2, vx, vy, vz)

            call smooth_field(vort(:, 1))
            call smooth_field(vort(:, 2))
            call smooth_field(vort(:, 3))

            call outpost(vort(1, 1), vort(1, 2), vort(1, 3), pr, t, 'vor')
         end if

         !---  > Compute and outpost vortex fields.
         if (ifvox .and. (ifaxis .eqv. .false.)) then

            ifvo = .false.; ifto = .true.
            call vortex_core(vort(:, 3), 'omega')
            call smooth_field(vort(:, 3))
            call outpost(vort(:, 1), vort(:, 2), vort(:, 3), pr, vort(:, 3), 'omR')
            ifvo = .true.; ifto = .true.

            if (.not. if3D) then ! 2D case -> no lambda2, no q
               ifvo = .false.; ifto = .true. ! just outposting omega field to temperature... v's and p ignored
               call outpost(vx, vy, vz, pr, vort(:, 3), 'vox')
               ifvo = .true.
            end if

         end if

         ifto = ifto_sav; ifpo = ifpo_sav; ifpsco(:) = ifpsco_sav(:)

      end if

   end subroutine nekStab_outpost

!-----------------------------------------------------------------------
! nekStab_comment — Verbose timestep progress output
!
! Purpose:
!   Prints timing statistics every 10 timesteps: mean time per step,
!   estimated remaining time, and wall-clock per nondimensional time.
!   Also stops the simulation if CFL > 10.
!-----------------------------------------------------------------------
   subroutine nekStab_comment
      use nekstab_krylov_subspace, only: dnekclock
      implicit none
      real, save :: eetime0, eetime1, eetime2, deltatime
      real telapsed, tpernondt, tmiss, ttime
      real ttime_stp

      !     if extrapolation is not OIFS -> ifchar = false
      !     if OIFS -> ifchar = .true. and CFL 2-5
      !     cases can have CFL > 1 in initial time steps
      if (courno > 10) then
         call nekStab_error('CFL > 10 stopping code')
      end if
      if (nio /= 0) return

      eetime2 = dnekclock()
      if (eetime0 == 0.0d0 .and. istep == 1) then
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
            write (6, "('      Mean time per timestep: ',F8.4,'  dev:',I8,'ms')") &
               ttime/istep, nint(((ttime/istep) - ttime_stp)*1000)
            write (6, "('      Estimated remaining time: ',I8,' hr ',I2,' min')") &
               int(tmiss), nint((tmiss - int(tmiss))*60)
            if (tpernondt > 60.0d0) then
               write (6, "('      Time per nondimensional time: ',F8.2,' sec')") tpernondt
            else
               write (6, "('      Time per nondimensional time: ',F8.2,' min ')") tpernondt/60.0d0
            end if
            write (6, "('      Current local time: ',F10.4,'  File:',I8)") &
               time - deltatime, int((time - deltatime)/param(14)) + 1
            print *, ''
         end if
      end if

   end subroutine nekStab_comment

!-----------------------------------------------------------------------
! nekStab_printNEKParams — Print Nek5000 parameters for sanity check
!
! Purpose:
!   Outputs all relevant Nek5000 and nekStab parameters at startup
!   for verification. Only prints on rank 0.
!-----------------------------------------------------------------------
   subroutine nekStab_printNEKParams
      implicit none

      call nekStab_log('Printing NEK parameters for sanity check')
      if (nid == 0) then
         write (6, *) 'P01=', param(1), 'density'
         write (6, *) 'P02=', param(2), 'viscosity (1/Re)'
         write (6, *) 'P07=', param(7), 'rhoCp'
         write (6, *) 'P08=', param(8), 'conductivity (1/(Re*Pr))'
         write (6, *) 'P10=', param(10), 'stop at endTime'
         write (6, *) 'P11=', param(11), 'stop at numSteps'
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
         write (6, *) 'P44=', param(44), '0:E-based addit Schwarz PnPn-2;1:A-based'
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
         write (6, *) 'y min,max,tot=', ymn, ymx, ymx - ymn
         write (6, *) 'z min,max,tot=', zmn, zmx, zmx - zmn
         write (6, *)
      end if

   end subroutine nekStab_printNEKParams

!-----------------------------------------------------------------------
! nekStab_energy — Write kinetic energy time series to file
!
! Purpose:
!   Computes volume-averaged kinetic energy of a velocity field
!   and appends it to a text file at prescribed intervals.
!
! Arguments:
!   px    [in] — x-velocity field
!   py    [in] — y-velocity field
!   pz    [in] — z-velocity field
!   pt    [in] — temperature field (for potential energy)
!   fname [in] — output filename
!   skip  [in] — output every 'skip' timesteps
!-----------------------------------------------------------------------
   subroutine nekStab_energy(px, py, pz, pt, fname, skip)
      use nekstab_krylov_subspace, only: lv, nv, nt
      implicit none
      real, dimension(lv), intent(in) :: px, py, pz
      real, dimension(lv, ldimt), intent(in) :: pt
      integer, intent(in) :: skip
      character(len=*), intent(in) :: fname
      character(len=256) :: fname_local ! Local copy for ifx compatibility
      real glsc3, uek, vek, wek, eek, pot

      eek = 0.50d0/volvm1
      uek = 0.0d0; vek = 0.0d0; wek = 0.0d0; pot = 0.0d0

      if (mod(istep, skip) == 0) then
         fname_local = fname
         uek = glsc3(px, px, bm1s, nv)
         vek = glsc3(py, py, bm1s, nv)
         if (if3D) wek = glsc3(pz, pz, bm1s, nv)
         if (ifheat) pot = glsc3(pt(:, 1), ym1, bm1s, nt)
         if (nid == 0) then
            open (unit=NEKSTAB_UNIT_ENERGY, file=trim(fname_local), action='write', status='unknown', position='append')
            write (NEKSTAB_UNIT_ENERGY, "(6E15.7)") time, uek*eek, vek*eek, wek*eek, (uek + vek + wek)*eek, pot*eek
            close (NEKSTAB_UNIT_ENERGY)
         end if
      end if

   end subroutine nekStab_energy

!-----------------------------------------------------------------------
! nekStab_enstrophy — Write enstrophy time series to file
!
! Purpose:
!   Computes volume-averaged enstrophy (vorticity squared) and
!   appends it to a text file at prescribed intervals.
!
! Arguments:
!   px    [in] — x-velocity field
!   py    [in] — y-velocity field
!   pz    [in] — z-velocity field
!   pt    [in] — temperature field (unused, kept for interface)
!   fname [in] — output filename
!   skip  [in] — output every 'skip' timesteps
!-----------------------------------------------------------------------
   subroutine nekStab_enstrophy(px, py, pz, pt, fname, skip)
      use nekstab_krylov_subspace, only: lv, nv
      implicit none
      real, dimension(lv), intent(in) :: px, py, pz
      real, dimension(lv, ldimt), intent(in) :: pt
      integer, intent(in) :: skip
      character(len=*), intent(in) :: fname
      character(len=256) :: fname_local ! Local copy for ifx compatibility
      real vort(lv, 3), wo1(lv), wo2(lv)
      common/ugrad/vort, wo1, wo2
      real glsc3, uek, vek, wek, eek

      eek = 0.50d0/volvm1
      uek = 0.0d0; vek = 0.0d0; wek = 0.0d0

      if (mod(istep, skip) == 0) then
         fname_local = fname
         call comp_vort3(vort, wo1, wo2, px, py, pz)
         uek = glsc3(vort(1, 1), vort(1, 1), bm1s, nv)
         vek = glsc3(vort(1, 2), vort(1, 2), bm1s, nv)
         if (if3D) wek = glsc3(vort(1, 3), vort(1, 3), bm1s, nv)
         if (nid == 0) then
            open (unit=NEKSTAB_UNIT_ENSTRO, file=trim(fname_local), action='write', status='unknown', position='append')
            write (NEKSTAB_UNIT_ENSTRO, "(5E15.7)") time, uek*eek, vek*eek, wek*eek, (uek + vek + wek)*eek
            close (NEKSTAB_UNIT_ENSTRO)
         end if
      end if

   end subroutine nekStab_enstrophy

!-----------------------------------------------------------------------
! nekStab_check_setup — Stop or warn when the setup can give a wrong answer
!
! Purpose:
!   Runs once at istep = 0, after the mode flags are resolved. It covers
!   the two setup errors that gave a plausible but wrong spectrum in the
!   cylinder_re180 and cubic_cavity_re1914 examples.
!
!   1. Polynomial order. Nek5000 interpolates a restart file written at
!      another order. The restart file of a stability or Newton mode is the
!      base flow or the orbit. Its interpolation is not the steady state or
!      the orbit of this discretization, so a Floquet multiplier of exactly 1
!      is lost. The run stops for the Floquet modes. For the steady
!      stability modes it prints a warning.
!   2. Sponge. nekStab_forcing applies the sponge to the base flow as well as
!      to the perturbation. The base flow must come from a run with the same
!      sponge. The run prints a note, since this cannot be tested here.
!
!   The header of the restart file is plain text. Only rank 0 reads it.
!-----------------------------------------------------------------------
   subroutine nekStab_check_setup
      implicit none
      logical :: needs_base, orbit_mode
      character(len=132) :: fname
      character(len=132) :: hdr
      integer :: i, u, ios, wdsz, nxr, nyr, nzr, iflag, nprobe

      !  Newton modes are not listed: a Newton run may start from a seed of
      !  another order (a DNS at lower order), and its result is then at the
      !  order of SIZE. The stability modes read that result.
      needs_base = isDirect .or. isAdjoint .or. isTransientGrowth .or. &
                   isFloquetDirect .or. isFloquetAdjoint .or. &
                   isFloquetTransientGrowth
      orbit_mode = isFloquetDirect .or. isFloquetAdjoint .or. &
                   isFloquetTransientGrowth

      !  Velocity residual projection. Every matvec of a Newton or Krylov mode
      !  starts from an unrelated vector, so the projection space carries no
      !  useful history across matvecs. On naca0012 (Newton) the velocity
      !  Helmholtz solve with projection stopped with a non-finite norm at 24
      !  or more ranks. With residualProj = no in [VELOCITY] it converged at 8, 16,
      !  24, 32 and 40 ranks, with identical residuals at 8 to 32.
      if (nid == 0 .and. ifprojfld(1) .and. &
          (needs_base .or. isNewtonFP .or. isNewtonPO .or. isNewtonPO_T)) then
         write (6, *) 'nekStab: NOTE: [VELOCITY] residualProj = yes in a Newton or'
         write (6, *) '  Krylov mode. If the run stops with "non-finite norm", set'
         write (6, *) '  residualProj = no in [VELOCITY]; naca0012 needed it at 24+ ranks.'
      end if

      !  History points. hpts sends the probes to ranks in blocks of lhis. When the
      !  number of probes is a multiple of lhis, the last block stays on a rank
      !  that hpts_out does not read, and the output file holds zeros for it.
      !  Measured on cylinder_re100: 2 probes with lhis = 2 wrote zeros, 2 probes
      !  with lhis = 3 wrote the values.
      if (nid == 0) then
         open (newunit=u, file=trim(adjustl(session))//'.his', status='old', &
               action='read', iostat=ios)
         if (ios == 0) then
            read (u, *, iostat=ios) nprobe
            close (u)
            if (ios == 0 .and. nprobe > 0 .and. lhis > 0) then
               if (mod(nprobe, lhis) == 0) then
                  write (6, *) 'nekStab: WARNING: ', nprobe, ' history points and lhis =', lhis
                  write (6, *) '  hpts writes zeros for the last block when the number of'
                  write (6, *) '  points is a multiple of lhis. Set lhis in SIZE to a value'
                  write (6, *) '  that does not divide the number of points (for example 100).'
               end if
            end if
         end if
      end if

      if (.not. needs_base) return

      iflag = 0 ! 0: fine, 1: order differs
      nxr = 0
      if (nid == 0) then
         fname = adjustl(initc(1))
         i = index(fname, ' ')
         if (i > 0) fname(i:) = ' ' ! first word only: the file name
         if (len_trim(fname) > 0) then
            open (newunit=u, file=trim(fname), access='stream', &
                  form='unformatted', status='old', action='read', iostat=ios)
            if (ios == 0) then
               hdr = ' '
               read (u, iostat=ios) hdr
               close (u)
               if (hdr(1:4) == '#std') then
                  read (hdr, '(5x,i1,1x,i2,1x,i2,1x,i2)', iostat=ios) &
                     wdsz, nxr, nyr, nzr
                  if (ios == 0 .and. nxr /= lx1) iflag = 1
               end if
            end if
         end if
      end if
      call bcast(iflag, isize)

      if (iflag == 1) then
         if (nid == 0) then
            write (6, *) 'nekStab: the restart file ', trim(fname)
            write (6, *) '  has lx1 =', nxr, ' but SIZE has lx1 =', lx1
         end if
         if (orbit_mode) then
            call nekStab_error('Floquet and orbit modes need a base flow written at '// &
                               'the same polynomial order as SIZE. Run the Newton '// &
                               'stage and this stage with the same lx1.')
         else if (nid == 0) then
            write (6, *) 'nekStab: WARNING: the base flow is interpolated from another'
            write (6, *) '  polynomial order. It is not a steady state of this'
            write (6, *) '  discretization, so the eigenvalues carry that error.'
         end if
      end if

      if (nid == 0 .and. spng_st /= 0) then
         write (6, *) 'nekStab: NOTE: the sponge is active in a stability mode.'
         write (6, *) '  It also forces the base flow. The base flow must come from'
         write (6, *) '  a run with the same sponge (same strength, length, field).'
      end if

   end subroutine nekStab_check_setup

end module nekstab_diagnostics
