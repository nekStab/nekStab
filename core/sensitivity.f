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

      real, dimension(lv) :: vx_dRe, vy_dRe, vz_dRe
      real, dimension(lv) :: vx_dIm, vy_dIm, vz_dIm

      real, dimension(lv) :: vx_aRe, vy_aRe, vz_aRe
      real, dimension(lv) :: vx_aIm, vy_aIm, vz_aIm

      real, dimension(lv) :: wavemaker, work1, work2

      character(len=80)   :: filename

      ifto=.false. ; ifpo=.false.

!     --> Load real part of the direct mode
      write(filename,'(a,a,a)')'dRe',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_dRe, vy_dRe, vz_dRe, vx, vy, vz)

!     --> Load imaginary part of the direct mode
      write(filename,'(a,a,a)')'dIm',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_dIm, vy_dIm, vz_dIm, vx, vy, vz)

!     --> Load real part of the adjoint mode
      write(filename,'(a,a,a)')'aRe',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call opcopy(vx_aRe, vy_aRe, vz_aRe, vx, vy, vz)

!     --> Load imaginary part of the adjoint mode
      write(filename,'(a,a,a)')'aIm',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call opcopy(vx_aIm, vy_aIm, vz_aIm, vx, vy, vz)

!     --> Normalize the adjoint mode.
      call biorthogonalize(vx_dRe, vy_dRe, vz_dRe, pr, t, 
     $     vx_dIm, vy_dIm, vz_dIm, pr, t, 
     $     vx_aRe, vy_aRe, vz_aRe, pr, t, 
     $     vx_aIm, vy_aIm, vz_aIm, pr, t)

!     --> Compute the wavemaker.
      work1 = sqrt(vx_dRe**2 + vx_dIm**2 + vy_dRe**2 + vy_dIm**2 + vz_dRe**2 + vz_dIm**2)
      work2 = sqrt(vx_aRe**2 + vx_aIm**2 + vy_aRe**2 + vy_aIm**2 + vz_aRe**2 + vz_aIm**2)
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
!     [1]

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real, dimension(lv) :: vx_dRe, vy_dRe, vz_dRe
      real, dimension(lv) :: vx_dIm, vy_dIm, vz_dIm
      real, dimension(lv) :: vx_aRe, vy_aRe, vz_aRe
      real, dimension(lv) :: vx_aIm, vy_aIm, vz_aIm
      real, dimension(lv) :: vx_aRet, vy_aRet, vz_aRet
      real, dimension(lv) :: vx_aImt, vy_aImt, vz_aImt

      real, dimension(lv) :: dudx_dRe, dudy_dRe, dudz_dRe
      real, dimension(lv) :: dvdx_dRe, dvdy_dRe, dvdz_dRe
      real, dimension(lv) :: dwdx_dRe, dwdy_dRe, dwdz_dRe

      real, dimension(lv) :: dudx_dIm, dudy_dIm, dudz_dIm
      real, dimension(lv) :: dvdx_dIm, dvdy_dIm, dvdz_dIm
      real, dimension(lv) :: dwdx_dIm, dwdy_dIm, dwdz_dIm

      real, dimension(lv) :: dudx_aRe, dudy_aRe, dudz_aRe
      real, dimension(lv) :: dvdx_aRe, dvdy_aRe, dvdz_aRe
      real, dimension(lv) :: dwdx_aRe, dwdy_aRe, dwdz_aRe

      real, dimension(lv) :: dudx_aIm, dudy_aIm, dudz_aIm
      real, dimension(lv) :: dvdx_aIm, dvdy_aIm, dvdz_aIm
      real, dimension(lv) :: dwdx_aIm, dwdy_aIm, dwdz_aIm

      real, dimension(lv) :: vx_tr,vy_tr,vz_tr
      real, dimension(lv) :: vx_ti,vy_ti,vz_ti
      real, dimension(lv) :: vx_pr,vy_pr,vz_pr
      real, dimension(lv) :: vx_pi,vy_pi,vz_pi
      real, dimension(lv) :: vx_sr,vy_sr,vz_sr

      character(len=80)   :: filename

      real :: alpha, beta, gamma, delta, epsilon, zeta, eta, theta, iota, kappa
      alpha   = 0.0d0
      beta    = 0.0d0
      gamma   = 0.0d0
      delta   = 0.0d0
      epsilon = 0.0d0
      zeta    = 0.0d0
      eta     = 0.0d0
      theta   = 0.0d0
      iota    = 0.0d0
      kappa   =-1.0d0

      ifto=.false.;ifpo=.false.

!     load real part of the direct mode
      write(filename,'(a,a,a)')'dRe',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_dRe,vy_dRe,vz_dRe,vx,vy,vz)

!     load imaginary part of the direct mode
      write(filename,'(a,a,a)')'dIm',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_dIm,vy_dIm,vz_dIm,vx,vy,vz)

!     load real part of the adjoint mode
      write(filename,'(a,a,a)')'aRe',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call opcopy(vx_aRe,vy_aRe,vz_aRe,vx,vy,vz)

!     load imaginary part of the adjoint mode
      write(filename,'(a,a,a)')'aIm',trim(SESSION),'0.f00002'
      call load_fld(filename)
      call opcopy(vx_aIm,vy_aIm,vz_aIm,vx,vy,vz)

!     --> Normalize the adjoint mode.
      call biorthogonalize(vx_dRe, vy_dRe, vz_dRe, pr, t,
     $     vx_dIm, vy_dIm, vz_dIm, pr, t, 
     $     vx_aRe, vy_aRe, vz_aRe, pr, t, 
     $     vx_aIm, vy_aIm, vz_aIm, pr, t)

!     gradient computation
!     real part of the direct mode
      call gradm1(dudx_dRe, dudy_dRe, dudz_dRe, vx_dRe, nelv)
      call gradm1(dvdx_dRe, dvdy_dRe, dvdz_dRe, vy_dRe, nelv)
      call gradm1(dwdx_dRe, dwdy_dRe, dwdz_dRe, vz_dRe, nelv)
      call dsavg(dudx_dRe);call dsavg(dudy_dRe); call dsavg(dudz_dRe)
      call dsavg(dvdx_dRe);call dsavg(dvdy_dRe); call dsavg(dvdz_dRe)
      call dsavg(dwdx_dRe);call dsavg(dwdy_dRe); call dsavg(dwdz_dRe)

!     imaginary part of the direct mode
      call gradm1(dudx_dIm, dudy_dIm, dudz_dIm, vx_dIm, nelv)
      call gradm1(dvdx_dIm, dvdy_dIm, dvdz_dIm, vy_dIm, nelv)
      call gradm1(dwdx_dIm, dwdy_dIm, dwdz_dIm, vz_dIm, nelv)
      call dsavg(dudx_dIm);call dsavg(dudy_dIm); call dsavg(dudz_dIm)
      call dsavg(dvdx_dIm);call dsavg(dvdy_dIm); call dsavg(dvdz_dIm)
      call dsavg(dwdx_dIm);call dsavg(dwdy_dIm); call dsavg(dwdz_dIm)

!     real part of the adjoint mode
      call gradm1(dudx_aRe, dudy_aRe, dudz_aRe, vx_aRe, nelv)
      call gradm1(dvdx_aRe, dvdy_aRe, dvdz_aRe, vy_aRe, nelv)
      call gradm1(dwdx_aRe, dwdy_aRe, dwdz_aRe, vz_aRe, nelv)
      call dsavg(dudx_aRe);call dsavg(dudy_aRe); call dsavg(dudz_aRe)
      call dsavg(dvdx_aRe);call dsavg(dvdy_aRe); call dsavg(dvdz_aRe)
      call dsavg(dwdx_aRe);call dsavg(dwdy_aRe); call dsavg(dwdz_aRe)

!     imaginary part of the adjoint mode
      call gradm1(dudx_aIm, dudy_aIm, dudz_aIm, vx_aIm, nelv)
      call gradm1(dvdx_aIm, dvdy_aIm, dvdz_aIm, vy_aIm, nelv)
      call gradm1(dwdx_aIm, dwdy_aIm, dwdz_aIm, vz_aIm, nelv)
      call dsavg(dudx_aIm);call dsavg(dudy_aIm); call dsavg(dudz_aIm)
      call dsavg(dvdx_aIm);call dsavg(dvdy_aIm); call dsavg(dvdz_aIm)
      call dsavg(dwdx_aIm);call dsavg(dwdy_aIm); call dsavg(dwdz_aIm)

!     computation of real part of base flow sensitivity term related to downstream transport of perturbations
      call oprzero  (vx_tr,vy_tr,vz_tr)
      call opaddcol3(vx_tr,vy_tr,vz_tr,-vx_aRe,-vx_aRe,-vx_aRe,dudx_dRe,dudy_dRe,dudz_dRe)
      call opaddcol3(vx_tr,vy_tr,vz_tr,-vy_aRe,-vy_aRe,-vy_aRe,dvdx_dRe,dvdy_dRe,dwdz_dRe)
      call opaddcol3(vx_tr,vy_tr,vz_tr,-vz_aRe,-vz_aRe,-vz_aRe,dwdx_dRe,dwdy_dRe,dwdz_dRe)
      call opaddcol3(vx_tr,vy_tr,vz_tr,-vx_aIm,-vx_aIm,-vx_aIm,dudx_dIm,dudy_dIm,dudz_dIm)
      call opaddcol3(vx_tr,vy_tr,vz_tr,-vy_aIm,-vy_aIm,-vy_aIm,dvdx_dIm,dvdy_dIm,dwdz_dIm)
      call opaddcol3(vx_tr,vy_tr,vz_tr,-vz_aIm,-vz_aIm,-vz_aIm,dwdx_dIm,dwdy_dIm,dwdz_dIm)

!     computation of imaginary part of base flow sensitivity term related to downstream transport of perturbations
      call oprzero  (vx_ti,vy_ti,vz_ti)
      call opaddcol3(vx_ti,vy_ti,vz_ti,vx_aRe,vx_aRe,vx_aRe,dudx_dIm,dudy_dIm,dudz_dIm)
      call opaddcol3(vx_ti,vy_ti,vz_ti,vy_aRe,vy_aRe,vy_aRe,dvdx_dIm,dvdy_dIm,dwdz_dIm)
      call opaddcol3(vx_ti,vy_ti,vz_ti,vz_aRe,vz_aRe,vz_aRe,dwdx_dIm,dwdy_dIm,dwdz_dIm)
      call opaddcol3(vx_ti,vy_ti,vz_ti,-vx_aIm,-vx_aIm,-vx_aIm,dudx_dRe,dudy_dRe,dudz_dRe)
      call opaddcol3(vx_ti,vy_ti,vz_ti,-vy_aIm,-vy_aIm,-vy_aIm,dvdx_dRe,dvdy_dRe,dwdz_dRe)
      call opaddcol3(vx_ti,vy_ti,vz_ti,-vz_aIm,-vz_aIm,-vz_aIm,dwdx_dRe,dwdy_dRe,dwdz_dRe)

!     computation of real part of base flow sensitivity term related to perturbations production
      call oprzero  (vx_pr,vy_pr,vz_pr)
      call opaddcol3(vx_pr,vy_pr,vz_pr,vx_dRe,vx_dRe,vx_dRe,dudx_aRe,dvdx_aRe,dwdx_aRe)
      call opaddcol3(vx_pr,vy_pr,vz_pr,vy_dRe,vy_dRe,vy_dRe,dudy_aRe,dvdy_aRe,dwdy_aRe)
      call opaddcol3(vx_pr,vy_pr,vz_pr,vz_dRe,vz_dRe,vz_dRe,dudz_aRe,dvdz_aRe,dwdz_aRe)
      call opaddcol3(vx_pr,vy_pr,vz_pr,vx_dIm,vx_dIm,vx_dIm,dudx_aIm,dvdx_aIm,dwdx_aIm)
      call opaddcol3(vx_pr,vy_pr,vz_pr,vy_dIm,vy_dIm,vy_dIm,dudy_aIm,dvdy_aIm,dwdy_aIm)
      call opaddcol3(vx_pr,vy_pr,vz_pr,vz_dIm,vz_dIm,vz_dIm,dudz_aIm,dvdz_aIm,dwdz_aIm)

!     computation of imaginary part of base flow sensitivity term related to perturbations production
      call oprzero  (vx_pi,vy_pi,vz_pi)
      call opaddcol3(vx_pi,vy_pi,vz_pi,vx_dRe,vx_dRe,vx_dRe,dudx_aIm,dvdx_aIm,dwdx_aIm)
      call opaddcol3(vx_pi,vy_pi,vz_pi,vy_dRe,vy_dRe,vy_dRe,dudy_aIm,dvdy_aIm,dwdy_aIm)
      call opaddcol3(vx_pi,vy_pi,vz_pi,vz_dRe,vz_dRe,vz_dRe,dudz_aIm,dvdz_aIm,dwdz_aIm)
      call opaddcol3(vx_pi,vy_pi,vz_pi,-vx_dIm,-vx_dIm,-vx_dIm,dudx_aRe,dvdx_aRe,dwdx_aRe)
      call opaddcol3(vx_pi,vy_pi,vz_pi,-vy_dIm,-vy_dIm,-vy_dIm,dudy_aRe,dvdy_aRe,dwdy_aRe)
      call opaddcol3(vx_pi,vy_pi,vz_pi,-vz_dIm,-vz_dIm,-vz_dIm,dudz_aRe,dvdz_aRe,dwdz_aRe)

      ifvo=.true.; ifpo=.false.; ifto=.false.
!     call filter_s0(vx_tr,0.5,1,'vortx')
!     call filter_s0(vy_tr,0.5,1,'vortx')
!     call filter_s0(vz_tr,0.5,1,'vortx')
      call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'tr_')

!     call filter_s0(vx_ti,0.5,1,'vortx')
!     call filter_s0(vy_ti,0.5,1,'vortx')
!     call filter_s0(vz_ti,0.5,1,'vortx')
      call outpost(vx_ti, vy_ti, vz_ti, pr, t, 'ti_')

!     call filter_s0(vx_pr,0.5,1,'vortx')
!     call filter_s0(vy_pr,0.5,1,'vortx')
!     call filter_s0(vz_pr,0.5,1,'vortx')
      call outpost(vx_pr, vy_pr, vz_pr, pr, t, 'pr_')

!     call filter_s0(vx_pi,0.5,1,'vortx')
!     call filter_s0(vy_pi,0.5,1,'vortx')
!     call filter_s0(vz_pi,0.5,1,'vortx')
      call outpost(vx_pi, vy_pi, vz_pi, pr, t, 'pi_')

      call opadd2(vx_tr, vy_tr, vz_tr, vx_pr, vy_pr, vz_pr)
      call opadd2(vx_ti, vy_ti, vz_ti, vx_pi, vy_pi, vz_pi)

      call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'sr_')
!     call opcmult(vx_ti, vy_ti, vz_ti, kappa)

!     ! check
!     call opchsgn(vx_ti, vy_ti, vz_ti)
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





      subroutine biorthogonalize(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, 
     $     vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, 
     $     vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, 
     $     vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm)

      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

!     ----- Real part of the direct mode.
      real, dimension(lv) :: vx_dRe, vy_dRe, vz_dRe, t_dRe
      real, dimension(lp) :: pr_dRe

!     ----- Imaginary part of the direct mode.
      real, dimension(lv) :: vx_dIm, vy_dIm, vz_dIm, t_dIm
      real, dimension(lp) :: pr_dIm

!     ----- Real part of the adjoint mode.
      real, dimension(lv) :: vx_aRe, vy_aRe, vz_aRe, t_aRe
      real, dimension(lp) :: pr_aRe

!     ----- Imaginary part of the adjoint mode.
      real, dimension(lv) :: vx_aIm, vy_aIm, vz_aIm, t_aIm
      real, dimension(lp) :: pr_aIm

!     ----- Temporary arrays.
      real, dimension(lv) :: work1_vx, work1_vy, work1_vz, work1_t
      real, dimension(lp) :: work1_pr

      real, dimension(lv) :: work2_vx, work2_vy, work2_vz, work2_t
      real, dimension(lp) :: work2_pr

      real :: alpha, beta, gamma, delta

!     --> Ensure that the direct mode is normalize to || u || = 1
      call norm(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, alpha)
      alpha = alpha**2

      call norm(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, beta)
      beta = beta**2

      gamma = 1.0D+00 / sqrt(alpha + beta)

      call opcmult(vx_dRe, vy_dRe, vz_dRe, gamma)
      call opcmult(vx_dIm, vy_dIm, vz_dIm, gamma)

!     --> Compute the scalar product between the direct and adjoint mode.
      call inner_product(alpha, vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe)
      call inner_product(beta, vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm, vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm)
      gamma = alpha + beta      ! Real part of the inner product

      call inner_product(alpha, vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm)
      call inner_product(beta, vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm, vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe)
      delta = alpha - beta      ! Complex part of the inner product

!     --> Bi-orthonormalize the adjoint mode.
      work1_vx = (gamma * vx_aRe - delta * vx_aIm) / (gamma**2 + delta**2)
      work2_vx = (gamma * vx_aIm + delta * vx_aRe) / (gamma**2 + delta**2)

      work1_vy = (gamma * vy_aRe - delta * vy_aIm) / (gamma**2 + delta**2)
      work2_vy = (gamma * vy_aIm + delta * vy_aRe) / (gamma**2 + delta**2)

      work1_vz = (gamma * vz_aRe - delta * vz_aIm) / (gamma**2 + delta**2)
      work2_vz = (gamma * vz_aIm + delta * vz_aRe) / (gamma**2 + delta**2)

      work1_pr = (gamma * pr_aRe - delta * pr_aIm) / (gamma**2 + delta**2)
      work2_pr = (gamma * pr_aIm + delta * pr_aRe) / (gamma**2 + delta**2)

      work1_t = (gamma * t_aRe - delta * t_aIm) / (gamma**2 + delta**2)
      work2_t = (gamma * t_aIm + delta * t_aRe) / (gamma**2 + delta**2)

      call nopcopy(vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, work1_vx, work1_vy, work1_vz, work1_pr, work1_t)
      call nopcopy(vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm, work2_vx, work2_vy, work2_vz, work2_pr, work2_t)

      return
      end subroutine biorthogonalize
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
