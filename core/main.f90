      !---------------------------------------------------------------------
      subroutine nekStab_setDefault
         !     specifying default values for nekStab

         implicit none
         include 'SIZE'
         include 'TOTAL'

         k_dim = 100               ! standard value, increas  in .usr
         schur_tgt = 2             ! schur target for schur step factorizaiton
         eigen_tol = 1.0e-6        ! tolerance for eigenmodes convergence
         schur_del = 0.10d0        !
         maxmodes = 20             ! max number of converged modes to disk
         glob_skip = 10            ! global energy computation skip frequency
         findiff_order = 1         ! finite difference order for the Frechet derivative
         epsilon_base = 1.0e-6     ! finite difference perturbation scale parameter

         bst_skp = 10              ! boostconv skip iterations
         bst_snp = 10              ! bootsconv residual subspace matrix size

         ifres = .false.          ! outpost restart files (KRY*, HES*)
         ifvor = .false.          ! outpost vorticity (vor* omega_x,omega_y,omega_z components)
         ifvox = .false.          ! outpost vortex (vox*: q,lambda2,omega criterions)
         ifldbf = .true.           ! load base flow for stability computations
         ifbf2D = .false.          ! force 2D base flow solution
         ifstorebase = .true.      ! store base flow for Floquet analysis (dynamic allocated)
         ifdyntol = .false.        ! dynamical tolerances for SFD and Newton (potential speed-up)

         ifseed_nois = .true.      ! noise as initial seed
         ifseed_symm = .false.     ! symmetry initial seed
         ifseed_load = .false.     ! loading initial seed (e.g. Re_ )
         !  Note: if ifseed_* all are false, 'useric' subroutine prescribes the initial seed

         !  Define here the probe position for zero-crossing vertical velocity analysis !
         xck = 2.0d0; call bcast(xck, wdsize)
         yck = 0.0d0; call bcast(yck, wdsize)
         zck = 0.0d0; call bcast(zck, wdsize)

         ! Sponge zone parameters (modified from KTH Toolbox)
         xLspg = 0.0d0; call bcast(xLspg, wdsize) ! x left
         xRspg = 0.0d0; call bcast(xRspg, wdsize) ! x right
         yLspg = 0.0d0; call bcast(yLspg, wdsize)
         yRspg = 0.0d0; call bcast(yRspg, wdsize)
         zLspg = 0.0d0; call bcast(zLspg, wdsize)
         zRspg = 0.0d0; call bcast(zRspg, wdsize)
         acc_spg = 0.333d0; call bcast(acc_spg, wdsize) !percentage for the acceleration phase in the sponge (e.g. 1/3)
         spng_st = 0.0d0; call bcast(spng_st, wdsize)

         evop = '_' ! initialize output prefix

         !     !Broadcast all defaults !
         call bcast(eigen_tol, wdsize) ! wdsize for real
         call bcast(schur_del, wdsize)
         call bcast(epsilon_base, wdsize)

         call bcast(schur_tgt, isize) ! isize for integer
         call bcast(maxmodes, isize)
         call bcast(k_dim, isize)
         call bcast(bst_skp, isize)
         call bcast(bst_snp, isize)
         call bcast(glob_skip, isize)
         call bcast(findiff_order, isize)

         call bcast(ifres, lsize) !lsize for boolean
         call bcast(ifvor, lsize)
         call bcast(ifvox, lsize)
         call bcast(ifseed_nois, lsize)
         call bcast(ifseed_symm, lsize)
         call bcast(ifseed_load, lsize)
         call bcast(ifldbf, lsize)
         call bcast(ifbf2D, lsize)
         call bcast(ifstorebase, lsize)
         call bcast(ifdyntol, lsize)

      end subroutine nekStab_setDefault
      !---------------------------------------------------------------------
      subroutine nekStab_init
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         logical scal
         real glmin, glmax
         integer i
         nv = nx1*ny1*nz1*nelv

         if (.not. isNekStabinit) then
            call nekStab_setDefault
            call nekStab_usrchk ! where user change defaults
            call nekStab_printNEKParams

            xmn = glmin(xm1, nv); xmx = glmax(xm1, nv)
            ymn = glmin(ym1, nv); ymx = glmax(ym1, nv)
            zmn = glmin(zm1, nv); zmx = glmax(zm1, nv)

            if (nid == 0) then
               print *, '                 __   _____  __          __  '
               print *, '   ____   ___   / /__/ ___/ / /_ ____ _ / /_ '
               print *, '  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \'
               print *, ' / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /'
               print *, '/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/ '
               print *, 'COPYRIGHT (c) 2020-2024 DynFluid Laboratoire Paris ', NSVERSION
               print *, 'Nek5000 ', NVERSION
               print *, ''
            end if

            call copy(bm1s, bm1, nv) ! never comment this !
            ifbfcv = .false.

            nof = 0
            scal = .false.
            do i = 1, size(ifpsco)
               if (ifpsco(i) .eqv. .true.) then
                  scal = .true.
                  nof = nof + 1
               end if
            end do
            if (ifto .eqv. .true. .or. scal .eqv. .true.) then
               if (nid == 0) write (6, *) 'Scalars found:'
               if (nid == 0) write (6, *) ' ifto=', ifto
               if (nid == 0) write (6, *) ' ifpsco=', ifpsco
               if (ifto) nof = nof + 1
               if (nid == 0) write (6, *) 'number of possible scalars (ldimt)=', ldimt
               if (nid == 0) write (6, *) 'number of scalars (nof)=', nof, npscal
            end if

            call oprzero(fcx, fcy, fcz) ! never comment this!
            call rzero(fct, nx1*ny1*nz1*nelv)

            isNekStabinit = .true.
         elseif (nid == 0) then
            print *, 'NekStab already initialized'
         end if

      end subroutine nekStab_init
      !---------------------------------------------------------------------
      subroutine nekStab_drive
         implicit none
         include 'SIZE'
         include 'TOTAL'

         if (istep == 0) call nekStab_init

         select case (floor(uparam(1)))

          case (0)                   ! DNS

            call nekStab_outpost   ! outpost vorticity
            call nekStab_comment   ! print comments
            call nekStab_energy(vx, vy, vz, t, 'total_energy.dat', glob_skip)
            call nekStab_enstrophy(vx, vy, vz, t, 'total_enstrophy.dat', glob_skip)
            if (lastep == 1) call nek_end

          case (1)                   ! fixed points computation

            call nekStab_outpost   ! outpost vorticity
            call nekStab_comment   ! print comments

            if (uparam(1) == 1.1) then
               call SFD
               if (uparam(5) == 0) call nekStab_energy(vx, vy, vz, t, 'total_energy.dat', glob_skip)
            elseif (uparam(1) == 1.2) then
               if (nid == 0) write (6, *) 'BOOSTCONV'
               call BoostConv
            elseif (uparam(1) == 1.3) then
               if (nid == 0) write (6, *) 'DMT'
               if (nid == 0) write (6, *) 'stopping ! not yet ported to this version'; call nek_end
               ! call DMT
            elseif (uparam(1) == 1.4) then
               if (nid == 0) write (6, *) 'TDF'
               call TDF
            end if

            if (ifbfcv) call nek_end

          case (2) ! Newton-Krylov solver

            ! Initialize flags based on the value of uparam(1)
            isNewtonFP = (uparam(1) == 2.0)
            isNewtonPO = (uparam(1) == 2.1)
            isNewtonPO_T = (uparam(1) == 2.2)

            ! Conditional statements for each Newton-Krylov case
            if (nid == 0) then
               if (isNewtonFP) then
                  write (6, *) 'Newton-Krylov for fixed points...'
               elseif (isNewtonPO) then
                  write (6, *) 'Newton-Krylov for UPOs...'
               elseif (isNewtonPO_T) then
                  write (6, *) 'Newton-Krylov for forced UPOs...'
               else
                  write (6, *) 'Unrecognized option...'
                  call nek_end
               end if
            end if

            ! Proceed with Newton-Krylov computation
            call newton_krylov
            call nek_end

          case (3)                   ! eigenvalue problem

            isDirect = (uparam(1) == 3.1)
            isFloquetDirect = (uparam(1) == 3.11)

            isAdjoint = (uparam(1) == 3.2)
            isFloquetAdjoint = (uparam(1) == 3.21)

            isTransientGrowth = (uparam(1) == 3.3)
            isFloquetTransientGrowth = (uparam(1) == 3.31)

            if (isDirect .or. isFloquetDirect .or. isAdjoint .or. isFloquetAdjoint) then
               call krylov_schur
            elseif (isTransientGrowth .or. isFloquetTransientGrowth) then
               call krylov_schur
            end if
            call nek_end

          case (4)                   ! in postprocessing.f

            if (uparam(01) == 4.0) then ! all
               call stability_energy_budget
               call wave_maker
               call bf_sensitivity
            end if

            !     -----> Direct mode kinetic energy budget.
            if (uparam(01) == 4.1) call stability_energy_budget

            !     -----> Wavemaker computation.
            if (uparam(01) == 4.2) call wave_maker

            !     -----> Baseflow sensitivity.
            if (uparam(01) == 4.3) call bf_sensitivity

            !     -----> Sensitivity to steady force.
            if (uparam(01) == 4.41 .or. uparam(01) == 4.42) call ts_steady_force_sensitivity
            if (uparam(01) == 4.43) call delta_forcing

            call nek_end

         end select

      end subroutine nekStab_drive
      !---------------------------------------------------------------------
