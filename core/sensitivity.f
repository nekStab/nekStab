!-----------------------------------------------------------------------





      subroutine wave_maker

!     Provided the direct and adjoint modes have already been computed,
!     this function computes the wavemaker following the formulation by
!     Giannetti et al. [1]. Set uparam(01) = 4.1 in the par file to use it.
!
!     OUTPOST
!     -------
!
!     wm_blah0.f000001 : Nek file. The wavemaker is stored in the array
!     for the temperature.
!
!     References
!     ----------
!
!     [1] Giannetti F. & Luchini P.
!     Structural sensitivity of the first instability of the cylinder wake.
!     J. Fluid Mech., vol 581., 2007.
!
!     NOTE : This implementation does not apply to cases involving temperature
!     or any other scalar.

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      type(krylov_vector) :: dRe_eig, dIm_eig
      type(krylov_vector) :: aRe_eig, aIm_eig

      real, dimension(lv) :: wavemaker, work1, work2

      character(len=80)   :: filename

      ifto=.false. ; ifpo=.false.

!     --> Load real part of the direct mode
      write(filename,'(a,a,a)')'dRe',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(dRe_eig%vx, dRe_eig%vy, dRe_eig%vz, vx, vy, vz)

!     --> Load imaginary part of the direct mode
      write(filename,'(a,a,a)')'dIm',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(dIm_eig%vx, dIm_eig%vy, dIm_eig%vz, vx, vy, vz)

!     --> Load real part of the adjoint mode
      write(filename,'(a,a,a)')'aRe',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call opcopy(aRe_eig%vx, aRe_eig%vy, aRe_eig%vz, vx, vy, vz)

!     --> Load imaginary part of the adjoint mode
      write(filename,'(a,a,a)')'aIm',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call opcopy(aIm_eig%vx, aIm_eig%vy, aIm_eig%vz, vx, vy, vz)

!     --> Normalize the adjoint mode.
      call krylov_biorthogonalize(dRe_eig, dIm_eig, aRe_eig, aIm_eig)

!     --> Compute the wavemaker.
      work1 = sqrt(dRe_eig%vx**2 + dIm_eig%vx**2 + dRe_eig%vy**2 + dIm_eig%vy**2 + dRe_eig%vz**2 + dIm_eig%vz**2)
      work2 = sqrt(aRe_eig%vx**2 + aIm_eig%vx**2 + aRe_eig%vy**2 + aIm_eig%vy**2 + aRe_eig%vz**2 + aIm_eig%vz**2)
      wavemaker = work1 * work2

      ifto = .true. ; ifvo = .false.
      call outpost(vx, vy, vz, pr, wavemaker, "wm_")

      return
      end subroutine wave_maker





!-----------------------------------------------------------------------





      subroutine bf_sensitivity

!     Provided the direct and adjoint modes have been computed,
!     this function computes the baseflow sensitivity following
!     the formulation by Marquet et al. [1].
!
!     OUTPOST
!     -------
!
!     References
!     ----------
!
!     [1] Marquet O., Sipp D. and Jacquin L.
!     Sensitivity analysis and passive control of cylinder flow
!     J. Fluid Mech., vol 615, pp. 221-252, 2008.

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      type(krylov_vector) :: dRe_eig, dIm_eig
      type(krylov_vector) :: aRe_eig, aIm_eig

      type(krylov_vector) :: dx_dRe_eig, dy_dRe_eig, dz_dRe_eig
      type(krylov_vector) :: dx_dIm_eig, dy_dIm_eig, dz_dIm_eig
      type(krylov_vector) :: dx_aRe_eig, dy_aRe_eig, dz_aRe_eig
      type(krylov_vector) :: dx_aIm_eig, dy_aIm_eig, dz_aIm_eig

      real, dimension(lv) :: vx_tr, vy_tr, vz_tr
      real, dimension(lv) :: vx_ti, vy_ti, vz_ti
      real, dimension(lv) :: vx_pr, vy_pr, vz_pr
      real, dimension(lv) :: vx_pi, vy_pi, vz_pi
      real, dimension(lv) :: vx_sr, vy_sr, vz_sr

      character(len=80)   :: filename

      ifto=.false.;ifpo=.false.

!     load real part of the direct mode
      write(filename,'(a,a,a)')'dRe',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call nopcopy(dRe_eig%vx, dRe_eig%vy, dRe_eig%vz, dRe_eig%pr, dRe_eig%theta, vx, vy, vz, pr, t)

!     load imaginary part of the direct mode
      write(filename,'(a,a,a)')'dIm',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call nopcopy(dIm_eig%vx, dIm_eig%vy, dIm_eig%vz, dIm_eig%pr, dIm_eig%theta, vx, vy, vz, pr, t)

!     load real part of the adjoint mode
      write(filename,'(a,a,a)')'aRe',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call nopcopy(aRe_eig%vx, aRe_eig%vy, aRe_eig%vz, aRe_eig%pr, aRe_eig%theta, vx, vy, vz, pr, t)

!     load imaginary part of the adjoint mode
      write(filename,'(a,a,a)')'aIm',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call nopcopy(aIm_eig%vx, aIm_eig%vy, aIm_eig%vz, aIm_eig%pr, aIm_eig%theta, vx, vy, vz, pr, t)

!     --> Normalize the adjoint mode.
      call krylov_biorthogonalize(dRe_eig, dIm_eig, aRe_eig, aIm_eig)

      call opcopy(dRe_eig%vx, dRe_eig%vy, dRe_eig%vz, dRe_eig%vx, dRe_eig%vy, dRe_eig%vz)
      call opcopy(dIm_eig%vx, dIm_eig%vy, dIm_eig%vz, dIm_eig%vx, dIm_eig%vy, dIm_eig%vz)
      call opcopy(aRe_eig%vx, aRe_eig%vy, aRe_eig%vz, aRe_eig%vx, aRe_eig%vy, aRe_eig%vz)
      call opcopy(aIm_eig%vx, aIm_eig%vy, aIm_eig%vz, aIm_eig%vx, aIm_eig%vy, aIm_eig%vz)

!     --> Compute the gradient of the direct and adjoint eigenmode.
      call krylov_gradient(dx_dRe_eig, dy_dRe_eig, dz_dRe_eig, dRe_eig)
      call krylov_gradient(dx_dIm_eig, dy_dIm_eig, dz_dIm_eig, dIm_eig)

      call krylov_gradient(dx_aRe_eig, dy_aRe_eig, dz_aRe_eig, aRe_eig)
      call krylov_gradient(dx_aIm_eig, dy_aIm_eig, dz_aIm_eig, aIm_eig)

!     --> Computation of real part of base flow sensitivity term related to downstream transport of perturbations
      call oprzero  (vx_tr, vy_tr, vz_tr)
      call opaddcol3(vx_tr, vy_tr, vz_tr, -aRe_eig%vx, -aRe_eig%vx, -aRe_eig%vx, dx_dRe_eig%vx, dy_dRe_eig%vx, dz_dRe_eig%vx)
      call opaddcol3(vx_tr, vy_tr, vz_tr, -aRe_eig%vy, -aRe_eig%vy, -aRe_eig%vy, dx_dRe_eig%vy, dy_dRe_eig%vy, dz_dRe_eig%vy)
      call opaddcol3(vx_tr, vy_tr, vz_tr, -aRe_eig%vz, -aRe_eig%vz, -aRe_eig%vz, dx_dRe_eig%vz, dy_dRe_eig%vz, dz_dRe_eig%vz)
      call opaddcol3(vx_tr, vy_tr, vz_tr, -aIm_eig%vx, -aIm_eig%vx, -aIm_eig%vx, dx_dIm_eig%vx, dy_dIm_eig%vx, dz_dIm_eig%vx)
      call opaddcol3(vx_tr, vy_tr, vz_tr, -aIm_eig%vy, -aIm_eig%vy, -aIm_eig%vy, dx_dIm_eig%vy, dy_dIm_eig%vy, dz_dIm_eig%vy)
      call opaddcol3(vx_tr, vy_tr, vz_tr, -aIm_eig%vz, -aIm_eig%vz, -aIm_eig%vz, dx_dIm_eig%vz, dy_dIm_eig%vz, dz_dIm_eig%vz)

!     --> Computation of imaginary part of base flow sensitivity term related to downstream transport of perturbations
      call oprzero  (vx_ti, vy_ti, vz_ti)
      call opaddcol3(vx_ti, vy_ti, vz_ti, aRe_eig%vx, aRe_eig%vx, aRe_eig%vx, dx_dIm_eig%vx, dy_dIm_eig%vx, dz_dIm_eig%vx)
      call opaddcol3(vx_ti, vy_ti, vz_ti, aRe_eig%vy, aRe_eig%vy, aRe_eig%vy, dx_dIm_eig%vy, dy_dIm_eig%vy, dz_dIm_eig%vy)
      call opaddcol3(vx_ti, vy_ti, vz_ti, aRe_eig%vz, aRe_eig%vz, aRe_eig%vz, dx_dIm_eig%vz, dy_dIm_eig%vz, dz_dIm_eig%vz)
      call opaddcol3(vx_ti, vy_ti, vz_ti, -aIm_eig%vx, -aIm_eig%vx, -aIm_eig%vx, dx_dRe_eig%vx, dy_dRe_eig%vx, dz_dRe_eig%vx)
      call opaddcol3(vx_ti, vy_ti, vz_ti, -aIm_eig%vy, -aIm_eig%vy, -aIm_eig%vy, dx_dRe_eig%vy, dy_dRe_eig%vy, dz_dRe_eig%vy)
      call opaddcol3(vx_ti, vy_ti, vz_ti, -aIm_eig%vz, -aIm_eig%vz, -aIm_eig%vz, dx_dRe_eig%vz, dy_dRe_eig%vz, dz_dRe_eig%vz)

!     --> Computation of real part of base flow sensitivity term related to perturbations production
      call oprzero  (vx_pr, vy_pr, vz_pr)
      call opaddcol3(vx_pr, vy_pr, vz_pr, dRe_eig%vx, dRe_eig%vx, dRe_eig%vx, dx_aRe_eig%vx, dx_aRe_eig%vy, dx_aRe_eig%vz)
      call opaddcol3(vx_pr, vy_pr, vz_pr, dRe_eig%vy, dRe_eig%vy, dRe_eig%vy, dy_aRe_eig%vx, dy_aRe_eig%vy, dy_aRe_eig%vz)
      call opaddcol3(vx_pr, vy_pr, vz_pr, dRe_eig%vz, dRe_eig%vz, dRe_eig%vz, dz_aRe_eig%vx, dz_aRe_eig%vy, dz_aRe_eig%vz)
      call opaddcol3(vx_pr, vy_pr, vz_pr, dIm_eig%vx, dIm_eig%vx, dIm_eig%vx, dx_aIm_eig%vx, dx_aIm_eig%vy, dx_aIm_eig%vz)
      call opaddcol3(vx_pr, vy_pr, vz_pr, dIm_eig%vy, dIm_eig%vy, dIm_eig%vy, dy_aIm_eig%vx, dy_aIm_eig%vy, dy_aIm_eig%vz)
      call opaddcol3(vx_pr, vy_pr, vz_pr, dIm_eig%vz, dIm_eig%vz, dIm_eig%vz, dz_aIm_eig%vx, dz_aIm_eig%vy, dz_aIm_eig%vz)

!     --> Computation of imaginary part of base flow sensitivity term related to perturbations production
      call oprzero  (vx_pi, vy_pi, vz_pi)
      call opaddcol3(vx_pi, vy_pi, vz_pi, dRe_eig%vx, dRe_eig%vx, dRe_eig%vx, dx_aIm_eig%vx, dx_aIm_eig%vy, dx_aIm_eig%vz)
      call opaddcol3(vx_pi, vy_pi, vz_pi, dRe_eig%vy, dRe_eig%vy, dRe_eig%vy, dy_aIm_eig%vx, dy_aIm_eig%vy, dy_aIm_eig%vz)
      call opaddcol3(vx_pi, vy_pi, vz_pi, dRe_eig%vz, dRe_eig%vz, dRe_eig%vz, dz_aIm_eig%vx, dz_aIm_eig%vy, dz_aIm_eig%vz)
      call opaddcol3(vx_pi, vy_pi, vz_pi, -dIm_eig%vx, -dIm_eig%vx, -dIm_eig%vx, dx_aRe_eig%vx, dx_aRe_eig%vy, dx_aRe_eig%vz)
      call opaddcol3(vx_pi, vy_pi, vz_pi, -dIm_eig%vy, -dIm_eig%vy, -dIm_eig%vy, dy_aRe_eig%vx, dy_aRe_eig%vy, dy_aRe_eig%vz)
      call opaddcol3(vx_pi, vy_pi, vz_pi, -dIm_eig%vz, -dIm_eig%vz, -dIm_eig%vz, dz_aRe_eig%vx, dz_aRe_eig%vy, dz_aRe_eig%vz)

!     --> Outpost the different fields (transport and production).
      ifvo=.true.; ifpo=.false.; ifto=.false.
      call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'tr_')
      call outpost(vx_ti, vy_ti, vz_ti, pr, t, 'ti_')
      call outpost(vx_pr, vy_pr, vz_pr, pr, t, 'pr_')
      call outpost(vx_pi, vy_pi, vz_pi, pr, t, 'pi_')

!     --> Compute/Outpost the sensitivity fields.
      call opadd2(vx_tr, vy_tr, vz_tr, vx_pr, vy_pr, vz_pr)
      call opadd2(vx_ti, vy_ti, vz_ti, vx_pi, vy_pi, vz_pi)
      call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'sr_')
      call outpost(vx_ti, vy_ti, vz_ti, pr, t, 'si_')

      return
      end subroutine bf_sensitivity





!-----------------------------------------------------------------------





      subroutine ts_steady_force_sensitivity

!     Provided the baseflow sensitivity has been computed,
!     this function computes the sensitivity of the flow to
!     a steady force following the formulation by Marquet at al. [1].
!     A time-stepper formulation of the problem is used and
!     the linearized system is solved using GMRES. Set uparam(01) = 4.41
!     to compute the real part and uparam(01) = 4.42 for the imaginary one.
!
!     OUTPOST
!     -------
!
!     fsr_blah0.f00001 / fsi_blah0.f00001 : Sensitivity fields.
!
!     References
!     ----------
!
!     [1] Marquet O., Sipp D. and Jacquin L.
!     Sensitivity analysis and passive control of cylinder flow
!     J. Fluid Mech., vol 615, pp. 221-252, 2008.

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

!     ----- Right-hand side : baseflow sensitivity
      type(krylov_vector) :: rhs

!     ----- Solution of the linear system.
      type(krylov_vector) :: sol

!     ----- Misc.
      character(len=80) :: filename
      character(len=3) :: prefix
      real :: alpha
      integer :: calls

!     --> Load base flow.
      write(filename, '(a, a, a)') 'BF_', trim(session), '0.f00001'
      call load_fld(filename)
      call opcopy(ubase, vbase, wbase, vx, vy, vz)

!     --> Load the forcing term.
      if (uparam(01) .eq. 4.41) then
         write(filename, '(a, a, a)') 'sr_', trim(session), '0.f00001'
         prefix = 'fsr'
      elseif (uparam(01) .eq. 4.42) then
         write(filename, '(a, a, a)') 'si_', trim(session), '0.f00001'
         prefix = 'fsi'
      endif
      call load_fld(filename)
      call opcopy(rhs%vx, rhs%vy, rhs%vz, vx, vy, vz)

!     --> Zero-out initial guess.
      call krylov_zero(sol)

!     --> Recast rhs into time-stepper/discrete-time framework.
      call initialize_rhs_ts_steady_force_sensitivity(rhs)

!     --> Normalize right-hand side for simplicity in gmres.
      call krylov_normalize(rhs, alpha)

!     --> Solve the linear system.
      call ts_gmres(rhs, sol, 10, k_dim, calls)

!     -->
      call krylov_cmult(sol, alpha)

!     --> Outpost solution.
      call outpost(sol%vx, sol%vy, sol%vz, sol%pr, sol%theta, prefix)

      return
      end subroutine ts_steady_force_sensitivity




!-----------------------------------------------------------------------





      subroutine initialize_rhs_ts_steady_force_sensitivity(rhs)

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

      type(krylov_vector) :: rhs

!     --> Setup the parameters for the linearized solver.
      ifpert = .true. ; ifadj = .true.
      call bcast(ifpert, lsize) ; call bcast(ifadj, lsize)

!     -->
      call opcopy(vx, vy, vz, ubase, vbase, wbase)
      if (ifheat) call copy(t, tbase, nx1*ny1*nz1*nelv)

!     --> General initialization of the linear solver.
      call prepare_linearized_solver()

!     -->
      ifbase = .false.

!     --> Zero-out the initial perturbation.
      call oprzero(vxp, vyp, vzp)

      time = 0.0D+00
      do istep = 1, nsteps
!     --> Pass the forcing to nek.
         call opcopy(fcx, fcy, fcz, rhs%vx, rhs%vy, rhs%vz)

!     --> Integrate forward in time.
         call nekstab_usrchk()
         call nek_advance()
      enddo

!     --> Copy the final solution as the new rhs for the time-stepper formulation.
      call opcopy(rhs%vx, rhs%vy, rhs%vz, vxp(1, 1), vyp(1, 1), vzp(1, 1))
      call oprzero(fcx, fcy, fcz)

      return
      end subroutine initialize_rhs_ts_steady_force_sensitivity






c-----------------------------------------------------------------------
      subroutine delta_forcing

!     Provided the base flow and the steady force sensitivity have been computed,
!     this function computes the variations of a leading eigenvalue induced
!     by a steady pointwise force according to eq. (5.1) by Marquet et al. [1].
!
!     OUTPOST
!     -------
!
!     dfr_blah0.f00001 : Eigenvalue variations (x_comp -> delta_lambda/alpha)
!                                              (y_comp -> delta_omega /alpha).
!
!     References
!     ----------
!
!     [1] Marquet O., Sipp D. and Jacquin L.
!     Sensitivity analysis and passive control of cylinder flow
!     J. Fluid Mech., vol 615, pp. 221-252, 2008.
      use krylov_subspace
      implicit none
      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP

      real, dimension(lv)  :: vx_bf, vy_bf, vz_bf
      real, dimension(lv)  :: fsrx, fsry, fsrz
      real, dimension(lv)  :: fsix, fsiy, fsiz
      real, dimension(lv)  :: work, workr, worki
      real, dimension(lv)  :: delta_lambda, delta_omega

!     ----- Misc.
      character(len=80) :: filename
      real :: alpha

      alpha   = 1.0d0
!     load base flow
      write(filename,'(a,a,a)')'BF_',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_bf, vy_bf, vz_bf, vx, vy, vz)

!     load growth rate sensitivity
      write(filename, '(a, a, a)') 'fsr', trim(session), '0.f00001'
      call load_fld(filename)
      call opcopy(fsrx, fsry, fsrz, vx, vy, vz)

!     load frequency sensitivity
      write(filename, '(a, a, a)') 'fsi', trim(session), '0.f00001'
      call load_fld(filename)
      call opcopy(fsix, fsiy, fsiz, vx, vy, vz)

      work  = sqrt(vx_bf**2 + vy_bf**2 + vz_bf**2)
      workr = fsrx * vx_bf + fsry * vy_bf + fsrz * vz_bf
      worki = fsix * vx_bf + fsiy * vy_bf + fsiz * vz_bf

      delta_lambda = - alpha * work * workr
      delta_omega  =   alpha * work * worki

      ifvo=.true.; ifpo=.false.; ifto=.false.
      call outpost(delta_lambda, delta_omega, t, pr, t, 'dfr')

      return
      end subroutine delta_forcing
c-----------------------------------------------------------------------
