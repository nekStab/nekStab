      !-----------------------------------------------------------------------
      ! sensitivity.f90 — Structural sensitivity and mode animation
      !
      ! Purpose:
      !   Computes wavemaker, baseflow sensitivity, steady force sensitivity,
      !   and eigenvalue variations. Also provides mode animation utilities
      !   for direct/adjoint eigenmodes (standard and Floquet).
      !
      ! Public interface:
      !   wave_maker, bf_sensitivity, ts_steady_force_sensitivity,
      !   initialize_rhs_ts_steady_force_sensitivity, biorthogonalize,
      !   delta_forcing, animate_mode_only, compute_omegaR,
      !   animate_mode
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL, ADJOINT
      !-----------------------------------------------------------------------

   module nekstab_sensitivity
   use krylov_subspace
   use nekstab_nek_bridge
   use nekstab_vectors
   use nekstab_io
   use nekstab_eigensolvers
   use nekstab_matvec
   use nekstab_diagnostics
   use nekstab_newton
   use nekstab_energy_budget
   use nekstab_torque_mod
   implicit none
   private
   public :: wave_maker, bf_sensitivity,&
      ts_steady_force_sensitivity,&
      initialize_rhs_ts_steady_force_sensitivity,&
      biorthogonalize, delta_forcing,&
      animate_mode_only, compute_omegaR,&
      animate_mode
   contains

      !-----------------------------------------------------------------------
      ! wave_maker — Compute the wavemaker from direct and adjoint modes
      !
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
      !-----------------------------------------------------------------------
    subroutine wave_maker

    use krylov_subspace

    !  All local arrays are allocatable (heap) to avoid stack overflow.
    !  With fixed-size `real, dimension(lv)`, a 3D case can put several GB
    !  on the stack per subroutine frame — exceeding OS stack limits.
    real, allocatable :: vx_dRe(:), vy_dRe(:), vz_dRe(:)
    real, allocatable :: vx_dIm(:), vy_dIm(:), vz_dIm(:)

    real, allocatable :: vx_aRe(:), vy_aRe(:), vz_aRe(:)
    real, allocatable :: vx_aIm(:), vy_aIm(:), vz_aIm(:)

    real, allocatable :: wavemaker(:), work1(:), work2(:)

   character(len=80) :: filename

   allocate(vx_dRe(lv), vy_dRe(lv), vz_dRe(lv))
   allocate(vx_dIm(lv), vy_dIm(lv), vz_dIm(lv))
   allocate(vx_aRe(lv), vy_aRe(lv), vz_aRe(lv))
   allocate(vx_aIm(lv), vy_aIm(lv), vz_aIm(lv))
   allocate(wavemaker(lv), work1(lv), work2(lv))

!     Defensive zero-init: opcopy/load_fld only fill [1:nv]; trailing
!     [nv+1:lv] (where lv = lelv max, nv = nelv*lx1**ndim runtime) was
!     left uninitialized -> whole-array `sqrt(a**2 + b**2 + ...)` ops at
!     line 115-117 squared garbage -> Inf/NaN in wm_ output.
   vx_dRe = 0.0d0; vy_dRe = 0.0d0; vz_dRe = 0.0d0
   vx_dIm = 0.0d0; vy_dIm = 0.0d0; vz_dIm = 0.0d0
   vx_aRe = 0.0d0; vy_aRe = 0.0d0; vz_aRe = 0.0d0
   vx_aIm = 0.0d0; vy_aIm = 0.0d0; vz_aIm = 0.0d0
   wavemaker = 0.0d0; work1 = 0.0d0; work2 = 0.0d0

   ifto = .false.; ifpo = .false.

      !     --> Load real part of the direct mode
   write (filename, '(a,a,a)') 'dRe', trim(SESSION), '0.f00001'
   call load_fld(filename)
   call opcopy(vx_dRe, vy_dRe, vz_dRe, vx, vy, vz)

      !     --> Load imaginary part of the direct mode
   write (filename, '(a,a,a)') 'dIm', trim(SESSION), '0.f00001'
   call load_fld(filename)
   call opcopy(vx_dIm, vy_dIm, vz_dIm, vx, vy, vz)

      !     --> Load real part of the adjoint mode
   write (filename, '(a,a,a)') 'aRe', trim(SESSION), '0.f00002'
   call load_fld(filename)
   call opcopy(vx_aRe, vy_aRe, vz_aRe, vx, vy, vz)

      !     --> Load imaginary part of the adjoint mode
   write (filename, '(a,a,a)') 'aIm', trim(SESSION), '0.f00002'
   call load_fld(filename)
   call opcopy(vx_aIm, vy_aIm, vz_aIm, vx, vy, vz)

      !     --> Normalize the adjoint mode.
   call biorthogonalize(vx_dRe, vy_dRe, vz_dRe, pr, t,&
      vx_dIm, vy_dIm, vz_dIm, pr, t,&
      vx_aRe, vy_aRe, vz_aRe, pr, t,&
      vx_aIm, vy_aIm, vz_aIm, pr, t)

      !     --> Compute the wavemaker.
   work1 = sqrt(vx_dRe**2 + vx_dIm**2 + vy_dRe**2 + vy_dIm**2 + vz_dRe**2 + vz_dIm**2)
   work2 = sqrt(vx_aRe**2 + vx_aIm**2 + vy_aRe**2 + vy_aIm**2 + vz_aRe**2 + vz_aIm**2)
   wavemaker = work1*work2

   ifto = .true.; ifvo = .false.
   call outpost(vx, vy, vz, pr, wavemaker, "wm_")

   deallocate(vx_dRe, vy_dRe, vz_dRe)
   deallocate(vx_dIm, vy_dIm, vz_dIm)
   deallocate(vx_aRe, vy_aRe, vz_aRe)
   deallocate(vx_aIm, vy_aIm, vz_aIm)
   deallocate(wavemaker, work1, work2)

   end subroutine wave_maker

      !-----------------------------------------------------------------------
      ! bf_sensitivity — Compute baseflow sensitivity from direct/adjoint modes
      !
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
      !-----------------------------------------------------------------------
   subroutine bf_sensitivity

   use krylov_subspace

   !  Allocatable (heap) — see wave_maker comment on stack overflow risk.
   real, allocatable :: vx_dRe(:), vy_dRe(:), vz_dRe(:)
   real, allocatable :: vx_dIm(:), vy_dIm(:), vz_dIm(:)
   real, allocatable :: vx_aRe(:), vy_aRe(:), vz_aRe(:)
   real, allocatable :: vx_aIm(:), vy_aIm(:), vz_aIm(:)
   real, allocatable :: dudx_dRe(:), dudy_dRe(:), dudz_dRe(:)
   real, allocatable :: dvdx_dRe(:), dvdy_dRe(:), dvdz_dRe(:)
   real, allocatable :: dwdx_dRe(:), dwdy_dRe(:), dwdz_dRe(:)
   real, allocatable :: dudx_dIm(:), dudy_dIm(:), dudz_dIm(:)
   real, allocatable :: dvdx_dIm(:), dvdy_dIm(:), dvdz_dIm(:)
   real, allocatable :: dwdx_dIm(:), dwdy_dIm(:), dwdz_dIm(:)
   real, allocatable :: dudx_aRe(:), dudy_aRe(:), dudz_aRe(:)
   real, allocatable :: dvdx_aRe(:), dvdy_aRe(:), dvdz_aRe(:)
   real, allocatable :: dwdx_aRe(:), dwdy_aRe(:), dwdz_aRe(:)
   real, allocatable :: dudx_aIm(:), dudy_aIm(:), dudz_aIm(:)
   real, allocatable :: dvdx_aIm(:), dvdy_aIm(:), dvdz_aIm(:)
   real, allocatable :: dwdx_aIm(:), dwdy_aIm(:), dwdz_aIm(:)
   real, allocatable :: vx_tr(:), vy_tr(:), vz_tr(:)
   real, allocatable :: vx_ti(:), vy_ti(:), vz_ti(:)
   real, allocatable :: vx_pr(:), vy_pr(:), vz_pr(:)
   real, allocatable :: vx_pi(:), vy_pi(:), vz_pi(:)

   character(len=80) :: filename

   allocate(vx_dRe(lv), vy_dRe(lv), vz_dRe(lv))
   allocate(vx_dIm(lv), vy_dIm(lv), vz_dIm(lv))
   allocate(vx_aRe(lv), vy_aRe(lv), vz_aRe(lv))
   allocate(vx_aIm(lv), vy_aIm(lv), vz_aIm(lv))
   allocate(dudx_dRe(lv), dudy_dRe(lv), dudz_dRe(lv))
   allocate(dvdx_dRe(lv), dvdy_dRe(lv), dvdz_dRe(lv))
   allocate(dwdx_dRe(lv), dwdy_dRe(lv), dwdz_dRe(lv))
   allocate(dudx_dIm(lv), dudy_dIm(lv), dudz_dIm(lv))
   allocate(dvdx_dIm(lv), dvdy_dIm(lv), dvdz_dIm(lv))
   allocate(dwdx_dIm(lv), dwdy_dIm(lv), dwdz_dIm(lv))
   allocate(dudx_aRe(lv), dudy_aRe(lv), dudz_aRe(lv))
   allocate(dvdx_aRe(lv), dvdy_aRe(lv), dvdz_aRe(lv))
   allocate(dwdx_aRe(lv), dwdy_aRe(lv), dwdz_aRe(lv))
   allocate(dudx_aIm(lv), dudy_aIm(lv), dudz_aIm(lv))
   allocate(dvdx_aIm(lv), dvdy_aIm(lv), dvdz_aIm(lv))
   allocate(dwdx_aIm(lv), dwdy_aIm(lv), dwdz_aIm(lv))
   allocate(vx_tr(lv), vy_tr(lv), vz_tr(lv))
   allocate(vx_ti(lv), vy_ti(lv), vz_ti(lv))
   allocate(vx_pr(lv), vy_pr(lv), vz_pr(lv))
   allocate(vx_pi(lv), vy_pi(lv), vz_pi(lv))

!     Defensive zero-init: same fix as in wave_maker (line ~85). Trailing
!     [nv+1:lv] elements are uninitialized after allocate; opcopy/opaddcol3
!     touch only [1:nv]; outpost writes element-bounded so the file is OK,
!     but any whole-array op (incl. unary minus expressions like `-vx_dIm`
!     creating Fortran temporaries) hits the garbage tail.
   vx_dRe = 0.0d0; vy_dRe = 0.0d0; vz_dRe = 0.0d0
   vx_dIm = 0.0d0; vy_dIm = 0.0d0; vz_dIm = 0.0d0
   vx_aRe = 0.0d0; vy_aRe = 0.0d0; vz_aRe = 0.0d0
   vx_aIm = 0.0d0; vy_aIm = 0.0d0; vz_aIm = 0.0d0
   dudx_dRe = 0.0d0; dudy_dRe = 0.0d0; dudz_dRe = 0.0d0
   dvdx_dRe = 0.0d0; dvdy_dRe = 0.0d0; dvdz_dRe = 0.0d0
   dwdx_dRe = 0.0d0; dwdy_dRe = 0.0d0; dwdz_dRe = 0.0d0
   dudx_dIm = 0.0d0; dudy_dIm = 0.0d0; dudz_dIm = 0.0d0
   dvdx_dIm = 0.0d0; dvdy_dIm = 0.0d0; dvdz_dIm = 0.0d0
   dwdx_dIm = 0.0d0; dwdy_dIm = 0.0d0; dwdz_dIm = 0.0d0
   dudx_aRe = 0.0d0; dudy_aRe = 0.0d0; dudz_aRe = 0.0d0
   dvdx_aRe = 0.0d0; dvdy_aRe = 0.0d0; dvdz_aRe = 0.0d0
   dwdx_aRe = 0.0d0; dwdy_aRe = 0.0d0; dwdz_aRe = 0.0d0
   dudx_aIm = 0.0d0; dudy_aIm = 0.0d0; dudz_aIm = 0.0d0
   dvdx_aIm = 0.0d0; dvdy_aIm = 0.0d0; dvdz_aIm = 0.0d0
   dwdx_aIm = 0.0d0; dwdy_aIm = 0.0d0; dwdz_aIm = 0.0d0
   vx_tr = 0.0d0; vy_tr = 0.0d0; vz_tr = 0.0d0
   vx_ti = 0.0d0; vy_ti = 0.0d0; vz_ti = 0.0d0
   vx_pr = 0.0d0; vy_pr = 0.0d0; vz_pr = 0.0d0
   vx_pi = 0.0d0; vy_pi = 0.0d0; vz_pi = 0.0d0

   ifto = .false.; ifpo = .false.

      !     load real part of the direct mode
   write (filename, '(a,a,a)') 'dRe', trim(SESSION), '0.f00001'
   call load_fld(filename)
   call opcopy(vx_dRe, vy_dRe, vz_dRe, vx, vy, vz)

      !     load imaginary part of the direct mode
   write (filename, '(a,a,a)') 'dIm', trim(SESSION), '0.f00001'
   call load_fld(filename)
   call opcopy(vx_dIm, vy_dIm, vz_dIm, vx, vy, vz)

      !     load real part of the adjoint mode
   write (filename, '(a,a,a)') 'aRe', trim(SESSION), '0.f00002'
   call load_fld(filename)
   call opcopy(vx_aRe, vy_aRe, vz_aRe, vx, vy, vz)

      !     load imaginary part of the adjoint mode
   write (filename, '(a,a,a)') 'aIm', trim(SESSION), '0.f00002'
   call load_fld(filename)
   call opcopy(vx_aIm, vy_aIm, vz_aIm, vx, vy, vz)

      !     --> Normalize the adjoint mode.
   call biorthogonalize(vx_dRe, vy_dRe, vz_dRe, pr, t,&
      vx_dIm, vy_dIm, vz_dIm, pr, t,&
      vx_aRe, vy_aRe, vz_aRe, pr, t,&
      vx_aIm, vy_aIm, vz_aIm, pr, t)

      !     gradient computation
      !     real part of the direct mode
   call compute_velocity_gradient_tensor(vx_dRe, vy_dRe, vz_dRe,&
      dudx_dRe, dudy_dRe, dudz_dRe,&
      dvdx_dRe, dvdy_dRe, dvdz_dRe,&
      dwdx_dRe, dwdy_dRe, dwdz_dRe)

      !     imaginary part of the direct mode
   call compute_velocity_gradient_tensor(vx_dIm, vy_dIm, vz_dIm,&
      dudx_dIm, dudy_dIm, dudz_dIm,&
      dvdx_dIm, dvdy_dIm, dvdz_dIm,&
      dwdx_dIm, dwdy_dIm, dwdz_dIm)

      !     real part of the adjoint mode
   call compute_velocity_gradient_tensor(vx_aRe, vy_aRe, vz_aRe,&
      dudx_aRe, dudy_aRe, dudz_aRe,&
      dvdx_aRe, dvdy_aRe, dvdz_aRe,&
      dwdx_aRe, dwdy_aRe, dwdz_aRe)

      !     imaginary part of the adjoint mode
   call compute_velocity_gradient_tensor(vx_aIm, vy_aIm, vz_aIm,&
      dudx_aIm, dudy_aIm, dudz_aIm,&
      dvdx_aIm, dvdy_aIm, dvdz_aIm,&
      dwdx_aIm, dwdy_aIm, dwdz_aIm)

      !     computation of real part of base flow sensitivity term related to downstream transport of perturbations
   call oprzero(vx_tr, vy_tr, vz_tr)
   call opaddcol3(vx_tr, vy_tr, vz_tr, -vx_aRe, -vx_aRe, -vx_aRe, dudx_dRe, dudy_dRe, dudz_dRe)
   call opaddcol3(vx_tr, vy_tr, vz_tr, -vy_aRe, -vy_aRe, -vy_aRe, dvdx_dRe, dvdy_dRe, dvdz_dRe)
   call opaddcol3(vx_tr, vy_tr, vz_tr, -vz_aRe, -vz_aRe, -vz_aRe, dwdx_dRe, dwdy_dRe, dwdz_dRe)
   call opaddcol3(vx_tr, vy_tr, vz_tr, -vx_aIm, -vx_aIm, -vx_aIm, dudx_dIm, dudy_dIm, dudz_dIm)
   call opaddcol3(vx_tr, vy_tr, vz_tr, -vy_aIm, -vy_aIm, -vy_aIm, dvdx_dIm, dvdy_dIm, dvdz_dIm)
   call opaddcol3(vx_tr, vy_tr, vz_tr, -vz_aIm, -vz_aIm, -vz_aIm, dwdx_dIm, dwdy_dIm, dwdz_dIm)

      !     computation of imaginary part of base flow sensitivity term related to downstream transport of perturbations
   call oprzero(vx_ti, vy_ti, vz_ti)
   call opaddcol3(vx_ti, vy_ti, vz_ti, vx_aRe, vx_aRe, vx_aRe, dudx_dIm, dudy_dIm, dudz_dIm)
   call opaddcol3(vx_ti, vy_ti, vz_ti, vy_aRe, vy_aRe, vy_aRe, dvdx_dIm, dvdy_dIm, dvdz_dIm)
   call opaddcol3(vx_ti, vy_ti, vz_ti, vz_aRe, vz_aRe, vz_aRe, dwdx_dIm, dwdy_dIm, dwdz_dIm)
   call opaddcol3(vx_ti, vy_ti, vz_ti, -vx_aIm, -vx_aIm, -vx_aIm, dudx_dRe, dudy_dRe, dudz_dRe)
   call opaddcol3(vx_ti, vy_ti, vz_ti, -vy_aIm, -vy_aIm, -vy_aIm, dvdx_dRe, dvdy_dRe, dvdz_dRe)
   call opaddcol3(vx_ti, vy_ti, vz_ti, -vz_aIm, -vz_aIm, -vz_aIm, dwdx_dRe, dwdy_dRe, dwdz_dRe)

      !     computation of real part of base flow sensitivity term related to perturbations production
   call oprzero(vx_pr, vy_pr, vz_pr)
   call opaddcol3(vx_pr, vy_pr, vz_pr, vx_dRe, vx_dRe, vx_dRe, dudx_aRe, dvdx_aRe, dwdx_aRe)
   call opaddcol3(vx_pr, vy_pr, vz_pr, vy_dRe, vy_dRe, vy_dRe, dudy_aRe, dvdy_aRe, dwdy_aRe)
   call opaddcol3(vx_pr, vy_pr, vz_pr, vz_dRe, vz_dRe, vz_dRe, dudz_aRe, dvdz_aRe, dwdz_aRe)
   call opaddcol3(vx_pr, vy_pr, vz_pr, vx_dIm, vx_dIm, vx_dIm, dudx_aIm, dvdx_aIm, dwdx_aIm)
   call opaddcol3(vx_pr, vy_pr, vz_pr, vy_dIm, vy_dIm, vy_dIm, dudy_aIm, dvdy_aIm, dwdy_aIm)
   call opaddcol3(vx_pr, vy_pr, vz_pr, vz_dIm, vz_dIm, vz_dIm, dudz_aIm, dvdz_aIm, dwdz_aIm)

      !     computation of imaginary part of base flow sensitivity term related to perturbations production
   call oprzero(vx_pi, vy_pi, vz_pi)
   call opaddcol3(vx_pi, vy_pi, vz_pi, vx_dRe, vx_dRe, vx_dRe, dudx_aIm, dvdx_aIm, dwdx_aIm)
   call opaddcol3(vx_pi, vy_pi, vz_pi, vy_dRe, vy_dRe, vy_dRe, dudy_aIm, dvdy_aIm, dwdy_aIm)
   call opaddcol3(vx_pi, vy_pi, vz_pi, vz_dRe, vz_dRe, vz_dRe, dudz_aIm, dvdz_aIm, dwdz_aIm)
   call opaddcol3(vx_pi, vy_pi, vz_pi, -vx_dIm, -vx_dIm, -vx_dIm, dudx_aRe, dvdx_aRe, dwdx_aRe)
   call opaddcol3(vx_pi, vy_pi, vz_pi, -vy_dIm, -vy_dIm, -vy_dIm, dudy_aRe, dvdy_aRe, dwdy_aRe)
   call opaddcol3(vx_pi, vy_pi, vz_pi, -vz_dIm, -vz_dIm, -vz_dIm, dudz_aRe, dvdz_aRe, dwdz_aRe)

   ifvo = .true.; ifpo = .false.; ifto = .false.
   call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'tr_')

   call outpost(vx_ti, vy_ti, vz_ti, pr, t, 'ti_')

   call outpost(vx_pr, vy_pr, vz_pr, pr, t, 'pr_')

   call outpost(vx_pi, vy_pi, vz_pi, pr, t, 'pi_')

   call opadd2(vx_tr, vy_tr, vz_tr, vx_pr, vy_pr, vz_pr)
   call opadd2(vx_ti, vy_ti, vz_ti, vx_pi, vy_pi, vz_pi)

   call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'sr_')
      !     ! check
      !     call opchsgn(vx_ti, vy_ti, vz_ti)
   call outpost(vx_ti, vy_ti, vz_ti, pr, t, 'si_')

   deallocate(vx_dRe, vy_dRe, vz_dRe, vx_dIm, vy_dIm, vz_dIm)
   deallocate(vx_aRe, vy_aRe, vz_aRe, vx_aIm, vy_aIm, vz_aIm)
   deallocate(dudx_dRe, dudy_dRe, dudz_dRe)
   deallocate(dvdx_dRe, dvdy_dRe, dvdz_dRe)
   deallocate(dwdx_dRe, dwdy_dRe, dwdz_dRe)
   deallocate(dudx_dIm, dudy_dIm, dudz_dIm)
   deallocate(dvdx_dIm, dvdy_dIm, dvdz_dIm)
   deallocate(dwdx_dIm, dwdy_dIm, dwdz_dIm)
   deallocate(dudx_aRe, dudy_aRe, dudz_aRe)
   deallocate(dvdx_aRe, dvdy_aRe, dvdz_aRe)
   deallocate(dwdx_aRe, dwdy_aRe, dwdz_aRe)
   deallocate(dudx_aIm, dudy_aIm, dudz_aIm)
   deallocate(dvdx_aIm, dvdy_aIm, dvdz_aIm)
   deallocate(dwdx_aIm, dwdy_aIm, dwdz_aIm)
   deallocate(vx_tr, vy_tr, vz_tr, vx_ti, vy_ti, vz_ti)
   deallocate(vx_pr, vy_pr, vz_pr, vx_pi, vy_pi, vz_pi)

   end subroutine bf_sensitivity

      !-----------------------------------------------------------------------
      ! ts_steady_force_sensitivity — Sensitivity to a steady force via GMRES
      !
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
      !-----------------------------------------------------------------------
   subroutine ts_steady_force_sensitivity

   use krylov_subspace

      !     ----- Right-hand side : baseflow sensitivity
   type(krylov_vector) :: rhs

      !     ----- Solution of the linear system.
   type(krylov_vector) :: sol

      !     ----- Misc.
   character(len=80) :: filename
   character(len=3) :: prefix
   real :: alpha
   integer :: calls, k_out, newton_iter
   real :: dtol

      !     --> Load base flow.
   write (filename, '(a, a, a)') 'BF_', trim(session), '0.f00001'
   call load_fld(filename)
   call opcopy(ubase, vbase, wbase, vx, vy, vz)

      !     --> Load the forcing term.
   if (uparam(01) == 4.41) then
   write (filename, '(a, a, a)') 'sr_', trim(session), '0.f00001'
   prefix = 'fsr'
   elseif (uparam(01) == 4.42) then
   write (filename, '(a, a, a)') 'si_', trim(session), '0.f00001'
   prefix = 'fsi'
   end if
   call load_fld(filename)
   call opcopy(rhs%vx, rhs%vy, rhs%vz, vx, vy, vz)

      !     --> Zero-out initial guess.
   call k_zero(sol)

      !     --> Recast rhs into time-stepper/discrete-time framework.
   call initialize_rhs_ts_steady_force_sensitivity(rhs)

      !     --> Normalize right-hand side for simplicity in gmres.
   call k_normalize(rhs, alpha)

      !     --> Solve the linear system.
   dtol = 1.0d-6
   newton_iter = 0
   call ts_gmres(rhs, sol, 10, k_dim, 1.0d-6, calls, k_out, newton_iter, dtol)

      !     -->
   call k_cmult(sol, alpha)

      !     --> Outpost solution.
   call outpost(sol%vx, sol%vy, sol%vz, sol%pr, sol%t, prefix)

   end subroutine ts_steady_force_sensitivity

      !-----------------------------------------------------------------------
      ! initialize_rhs_ts_steady_force_sensitivity — Prepare RHS for force sensitivity
      !-----------------------------------------------------------------------
   subroutine initialize_rhs_ts_steady_force_sensitivity(rhs)

   use krylov_subspace

   type(krylov_vector), intent(inout) :: rhs

      !     --> Setup the parameters for the linearized solver.
   ifpert = .true.; ifadj = .true.
   call bcast(ifpert, lsize); call bcast(ifadj, lsize)

      !     -->
   call nopcopy(vx, vy, vz, pr, t, ubase, vbase, wbase, pbase, tbase)

      !     --> General initialization of the linear solver.
   call prepare_linearized_solver()

      !     -->
   ifbase = .false.

      !     --> Zero-out the initial perturbation.
   call oprzero(vxp, vyp, vzp)

   time = 0.0d+00
   do istep = 1, nsteps
      !     --> Pass the forcing to nek.
   call opcopy(fcx, fcy, fcz, rhs%vx, rhs%vy, rhs%vz)

      !     --> Integrate forward in time.
    call nekstab_usrchk()
    call nek_advance()
    end do

    !     --> Copy the final solution as the new rhs for the time-stepper formulation.
    call opcopy(rhs%vx, rhs%vy, rhs%vz, vxp(1, 1), vyp(1, 1), vzp(1, 1))
    call zero_forcing()

   end subroutine initialize_rhs_ts_steady_force_sensitivity

      !-----------------------------------------------------------------------
      ! biorthogonalize — Bi-orthonormalize adjoint mode w.r.t. direct mode
      !-----------------------------------------------------------------------
   subroutine biorthogonalize( &
      vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, &
      vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, &
      vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, &
      vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm)

   use krylov_subspace

      !     ----- Real part of the direct mode.
   !     TODO(CHT): t_dRe etc. are dimension(lv) but norm/inner_product
   !     expect (lt, ldimt). For CHT (lelt > lelv, lt > lv) this is a
   !     stride mismatch that would corrupt inner products. Safe only
   !     when lelt == lelv (no conjugate heat transfer).
   real, dimension(lv), intent(inout) :: vx_dRe, vy_dRe, vz_dRe
   real, dimension(lv), intent(inout) :: t_dRe
   real, dimension(lp), intent(inout) :: pr_dRe

      !     ----- Imaginary part of the direct mode.
   real, dimension(lv), intent(inout) :: vx_dIm, vy_dIm, vz_dIm
   real, dimension(lv), intent(inout) :: t_dIm
   real, dimension(lp), intent(inout) :: pr_dIm

      !     ----- Real part of the adjoint mode.
   real, dimension(lv), intent(inout) :: vx_aRe, vy_aRe, vz_aRe
   real, dimension(lv), intent(inout) :: t_aRe
   real, dimension(lp), intent(inout) :: pr_aRe

      !     ----- Imaginary part of the adjoint mode.
   real, dimension(lv), intent(inout) :: vx_aIm, vy_aIm, vz_aIm
   real, dimension(lv), intent(inout) :: t_aIm
   real, dimension(lp), intent(inout) :: pr_aIm

      !     ----- Temporary arrays (allocatable to avoid stack overflow).
   real, allocatable, dimension(:) :: wk1_vx, wk1_vy, wk1_vz, wk1_t
   real, allocatable, dimension(:) :: wk1_pr

   real, allocatable, dimension(:) :: wk2_vx, wk2_vy, wk2_vz, wk2_t
   real, allocatable, dimension(:) :: wk2_pr

   real :: alpha, beta, gamma, delta

   allocate(wk1_vx(lv), wk1_vy(lv), wk1_vz(lv), wk1_t(lv), wk1_pr(lp))
   allocate(wk2_vx(lv), wk2_vy(lv), wk2_vz(lv), wk2_t(lv), wk2_pr(lp))

      !     --> Ensure that the direct mode is normalize to || u || = 1
   call norm(vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe, alpha)
   alpha = alpha**2

   call norm(vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm, beta)
   beta = beta**2

   if (alpha + beta < 1e-60) then
   if (nid == 0) write (6, *) 'WARNING: zero-norm direct mode in biorthogonalize'
   return
   end if
   gamma = 1.0d+00/sqrt(alpha + beta)

   call opcmult(vx_dRe, vy_dRe, vz_dRe, gamma)
   call opcmult(vx_dIm, vy_dIm, vz_dIm, gamma)

      !     --> Compute the scalar product between the direct and adjoint mode.
   call inner_product(alpha, vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe)
   call inner_product(beta, vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm, vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm)
   gamma = alpha + beta ! Real part of the inner product

   call inner_product(alpha, vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, vx_dIm, vy_dIm, vz_dIm, pr_dIm, t_dIm)
   call inner_product(beta, vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm, vx_dRe, vy_dRe, vz_dRe, pr_dRe, t_dRe)
   delta = alpha - beta ! Complex part of the inner product

      !     --> Bi-orthonormalize the adjoint mode.
   if (gamma**2 + delta**2 < 1e-60) then
   if (nid == 0) write (6, *) 'WARNING: degenerate adjoint-direct pairing in biorthogonalize'
   return
   end if
   wk1_vx = (gamma*vx_aRe - delta*vx_aIm)/(gamma**2 + delta**2)
   wk2_vx = (gamma*vx_aIm + delta*vx_aRe)/(gamma**2 + delta**2)

   wk1_vy = (gamma*vy_aRe - delta*vy_aIm)/(gamma**2 + delta**2)
   wk2_vy = (gamma*vy_aIm + delta*vy_aRe)/(gamma**2 + delta**2)

   wk1_vz = (gamma*vz_aRe - delta*vz_aIm)/(gamma**2 + delta**2)
   wk2_vz = (gamma*vz_aIm + delta*vz_aRe)/(gamma**2 + delta**2)

   wk1_pr = (gamma*pr_aRe - delta*pr_aIm)/(gamma**2 + delta**2)
   wk2_pr = (gamma*pr_aIm + delta*pr_aRe)/(gamma**2 + delta**2)

   wk1_t = (gamma*t_aRe - delta*t_aIm)/(gamma**2 + delta**2)
   wk2_t = (gamma*t_aIm + delta*t_aRe)/(gamma**2 + delta**2)

   call nopcopy(vx_aRe, vy_aRe, vz_aRe, pr_aRe, t_aRe, wk1_vx, wk1_vy, wk1_vz, wk1_pr, wk1_t)
   call nopcopy(vx_aIm, vy_aIm, vz_aIm, pr_aIm, t_aIm, wk2_vx, wk2_vy, wk2_vz, wk2_pr, wk2_t)

   end subroutine biorthogonalize

      !-----------------------------------------------------------------------
      ! delta_forcing — Eigenvalue variations from steady pointwise force
      !
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
      !-----------------------------------------------------------------------
   subroutine delta_forcing

   use krylov_subspace

   !  Allocatable (heap) — see wave_maker comment on stack overflow risk.
   real, allocatable, dimension(:) :: vx_bf, vy_bf, vz_bf
   real, allocatable, dimension(:) :: fsrx, fsry, fsrz
   real, allocatable, dimension(:) :: fsix, fsiy, fsiz
   real, allocatable, dimension(:) :: work, workr, worki
   real, allocatable, dimension(:) :: delta_lambda, delta_omega

      !     ----- Misc.
   character(len=80) :: filename
   real :: alpha

   allocate(vx_bf(lv), vy_bf(lv), vz_bf(lv))
   allocate(fsrx(lv), fsry(lv), fsrz(lv))
   allocate(fsix(lv), fsiy(lv), fsiz(lv))
   allocate(work(lv), workr(lv), worki(lv))
   allocate(delta_lambda(lv), delta_omega(lv))

   alpha = 1.0d0
      !     load base flow
   write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
   call load_fld(filename)
   call opcopy(vx_bf, vy_bf, vz_bf, vx, vy, vz)

      !     load growth rate sensitivity
   write (filename, '(a, a, a)') 'fsr', trim(session), '0.f00001'
   call load_fld(filename)
   call opcopy(fsrx, fsry, fsrz, vx, vy, vz)

      !     load frequency sensitivity
   write (filename, '(a, a, a)') 'fsi', trim(session), '0.f00001'
   call load_fld(filename)
   call opcopy(fsix, fsiy, fsiz, vx, vy, vz)

   work = sqrt(vx_bf**2 + vy_bf**2 + vz_bf**2)
   workr = fsrx*vx_bf + fsry*vy_bf + fsrz*vz_bf
   worki = fsix*vx_bf + fsiy*vy_bf + fsiz*vz_bf

   delta_lambda = -alpha*work*workr
   delta_omega = alpha*work*worki

   ifvo = .true.; ifpo = .false.; ifto = .false.
   call outpost(delta_lambda, delta_omega, t, pr, t, 'dfr')

   end subroutine delta_forcing

      !-----------------------------------------------------------------------
      ! animate_mode_only — Animate eigenmode snapshots at different phases
      !-----------------------------------------------------------------------
   subroutine animate_mode_only(num_steps, mode)

      !  Animate the eigenmode by computing and outputting snapshots
      !  of the perturbation field at different phases of the period.

   use krylov_subspace

   integer, intent(in) :: num_steps
   character(len=*), intent(in) :: mode ! 'd' or 'a'

   type(krylov_vector) :: BF, Re, Im, Re_sin
   character(len=80) :: filename
   character(len=32) :: mode_local ! ifx: trim() on assumed-length args can segfault
   real :: frequency, omega, sigma, u_max, A0
   integer :: i

   mode_local = mode ! copy to local before trim() for ifx compatibility

   if (nid == 0) then
   write (6, *) 'Animating mode function in mode:', mode
   write (6, *) 'Number of steps:', num_steps
   end if

   call read_eigenvalue(sigma, omega)

! frequency = 1.0 / param(10)  ! f = 1 / T (T = period)
! frequency = omega / (2.0d0 * NEKSTAB_PI)  ! f = omega / (2 * pi)
! omega = 2.0d0 * NEKSTAB_PI * frequency  ! omega = 2 * pi * f

   call k_load(BF, 'BF_'//trim(SESSION)//'0.f00001')
   call compute_omegaR(BF%vx, BF%vy, BF%vz, BF%t(:, 1))
   ifto = .true.
   call outpost(BF%vx, BF%vy, BF%vz, BF%pr, BF%t, 'BF_')

   A0 = 1.0e-3
   sigma = 1.0e-1 ! force a value of sigma

   call load_mode_pair(mode, Re, Im)

! Loop over num_steps to create snapshots
   do i = 1, num_steps

   time = i*(param(10)/num_steps)

   if (nid == 0) then
   write (6, '(A, I0, A, F8.4, A, I0)') 'i: ', i, ', time: ', time, ', n: ', num_steps
   write (6, '(A, F8.4, A, F8.4, A, F8.4, A)')&
      'Time/param(10): ', time/param(10),&
      ' (', time,&
      ' / ', param(10),&
      ') [time/period]'
   end if

   call k_copy(Re_sin, Re)
   call k_axpby(Re_sin, cos(omega*time), Im, -sin(omega*time))
   call outpost2(Re_sin%vx, Re_sin%vy, Re_sin%vz, Re_sin%pr, Re_sin%t, nof, trim(mode_local)//'Q_')

   end do

   end subroutine animate_mode_only

      !-----------------------------------------------------------------------
      ! compute_omegaR — Compute rotation-rate indicator field
      !-----------------------------------------------------------------------
   subroutine compute_omegaR(vx_in, vy_in, vz_in, omegaR)

   use krylov_subspace

   integer, parameter :: nxyz = lx1*ly1*lz1
   real, intent(in), dimension(lx1, ly1, lz1, lelv) :: vx_in, vy_in, vz_in
   real, intent(out), dimension(lx1, ly1, lz1, lelv) :: omegaR
   real, dimension(lx1, ly1, lz1, lelv) :: norm_a, norm_b, eps_field
   real, dimension(lx1*ly1*lz1, ldim, ldim) :: gije
   real, dimension(ldim, ldim) :: ss, oo
   real omegaR_max, omegaR_min, glmax, glmin
   real, save :: omega, optimal_eps = 2.0d0
   real :: eps = 0.0d0
   integer ie, l, i, j, iter

   do ie = 1, nelv
   call comp_gije(gije, vx_in(1, 1, 1, ie), vy_in(1, 1, 1, ie), vz_in(1, 1, 1, ie), ie)
   do l = 1, nxyz; do j = 1, ldim; do i = 1, ldim
   ss(i, j) = 0.50d0*(gije(l, i, j) + gije(l, j, i))
   oo(i, j) = 0.50d0*(gije(l, i, j) - gije(l, j, i))
   end do; end do
   norm_a(l, 1, 1, ie) = (norm2(ss)**2)**2
   norm_b(l, 1, 1, ie) = (norm2(oo)**2)**2
   eps_field(l, 1, 1, ie) = norm_b(l, 1, 1, ie) - norm_a(l, 1, 1, ie)
   end do
   end do

   if (optimal_eps >= 1.0d0) then

   eps = glmax(eps_field, nv)*0.001d0
   if (nid == 0) write (6, *) 'Initial eps = ', eps

   do iter = 1, 20
   do ie = 1, nelv
   do l = 1, nxyz
   omegaR(l, 1, 1, ie) = norm_b(l, 1, 1, ie)/(norm_a(l, 1, 1, ie) + norm_b(l, 1, 1, ie) + eps)
   end do
   end do

   omegaR_max = glmax(omegaR, nv)
   if (abs(omegaR_max - 1.0d0) < 1.0e-8) exit

   omegaR_min = glmin(omegaR, nv)
   if (nid == 0) write (6, *) 'Iteration ', iter, ': omegaR min,max:', omegaR_min, omegaR_max

   if (omegaR_max > 1.0d0) then
   if (nid == 0) write (6, *) 'Increasing eps: ', eps
   eps = eps*2.0d0
   else if (omegaR_max < 1.0d0) then
   if (nid == 0) write (6, *) 'Decreasing eps: ', eps
   eps = eps*0.50d0
   end if

   optimal_eps = eps
   if (nid == 0) write (6, *) 'Final eps value: ', optimal_eps
   call bcast(optimal_eps, wdsize)
   end do ! iter
   end if ! optimal_eps

   do ie = 1, nelv; do l = 1, nxyz
   omegaR(l, 1, 1, ie) = norm_b(l, 1, 1, ie)/(norm_a(l, 1, 1, ie) + norm_b(l, 1, 1, ie) + optimal_eps)
   end do; end do

   omegaR_min = glmin(omegaR, nv)
   omegaR_max = glmax(omegaR, nv)

   if (nid == 0) write (6, '(A,3ES12.4)') 'eps, omegaR min,max:', optimal_eps, omegaR_min, omegaR_max
   if (omegaR_max > 1.05d0) then
   optimal_eps = 2.0d0
   else if (omegaR_max < -0.05d0) then
   optimal_eps = 2.0d0
   end if

   call smooth_field(omegaR) ! smoothing causes out of bound values
   do ie = 1, nelv; do l = 1, nxyz ! clip [0,1] (tested on cylinder, produce better results)
   omegaR(l, 1, 1, ie) = min(1.0d0, max(0.0d0, merge(0.0d0, omegaR(l, 1, 1, ie),&
      abs(omegaR(l, 1, 1, ie)) < epsilon(1.0d0))))
   end do; end do

   end subroutine compute_omegaR

      !-----------------------------------------------------------------------
      ! animate_mode — Animate eigenmode with base flow superposition
      !-----------------------------------------------------------------------
   subroutine animate_mode(num_of_files, mode, floquet_flag)

      !  Animate the eigenmode by computing and outputting snapshots
      !  of the perturbation field at different phases of the period.

   use krylov_subspace

   integer, intent(in) :: num_of_files
   character(len=*), intent(in) :: mode ! 'd' or 'a'
   logical, intent(in) :: floquet_flag

   type(krylov_vector), save :: BF, Re, Im
   type(krylov_vector) :: Re_cos
   real, save :: frequency, omega, sigma, u_max, A0
   character(len=80) :: filename
   character(len=32) :: mode_local ! ifx: trim() on assumed-length args can segfault
   real :: amplitude, period
   integer :: i, iosteps, nfiles
   integer, save :: counter = 1

   mode_local = mode ! copy to local before trim() for ifx compatibility

   if (nid == 0) then
   write (6, *) 'Animating mode function in mode:', mode
   write (6, *) 'Number of steps:', num_of_files
   end if

   call read_eigenvalue(sigma, omega)

! frequency = 1.0 / param(10)  ! f = 1 / T (T = period)
! frequency = omega / (2.0d0 * NEKSTAB_PI)  ! f = omega / (2 * pi)
! omega = 2.0d0 * NEKSTAB_PI * frequency  ! omega = 2 * pi * f

   if (floquet_flag) then

   A0 = 1.0e-3
   sigma = 1.0e-1 ! force a value of sigma

   call load_mode_pair(mode, Re, Im)

   period = 2.0d0*NEKSTAB_PI/omega ! period of leading mode
   fintim = 1.0d0*period ! add more periods here !
! compute for fintim and not period, as we might have n periods
   timeio = fintim/dble(num_of_files)

   call load_fld('BF_'//trim(SESSION)//'0.f00001') ! load velocity field to compute ctarg
   call compute_cfl(ctarg, vx, vy, vz, 1.0d0)
   if (nid == 0) write (6, *) 'Maximum spatial restriction:', ctarg

   dt = param(26)/ctarg ! param(26) is the max CFL specified by the user
   if (.not. (dt > 0.0d0 .and. dt < 1.0d3)) then
      if (nid == 0) write (6, *) &
          'WARN: compute_cfl gave unusable result (ctarg=', ctarg, &
          '); falling back to dt = 1e-3. Set explicit dt in .par ', &
          'to silence this.'
      dt = 1.0d-3
   end if
   nsteps = ceiling(fintim/dt)
   dt = fintim/dble(nsteps)

   iosteps = nsteps/num_of_files
   if (iosteps*num_of_files /= nsteps) then
   nsteps = num_of_files*ceiling(dble(nsteps)/dble(num_of_files))
   dt = period/dble(nsteps)
   iosteps = nsteps/num_of_files
   end if

   timeio = dt*dble(iosteps)
   nfiles = num_of_files

   if (nid == 0) then
   write (6, '(A,E15.7)') 'Time step:', dt
   write (6, '(A,I8)') 'Total steps:', nsteps
   write (6, '(A,I8)') 'Steps between outputs:', iosteps
   write (6, '(A,I8)') 'Number of files to output:', nfiles
   write (6, '(A,E15.7)') 'Output interval:', timeio
   end if

   call compute_cfl(ctarg, vx, vy, vz, dt)
   if (nid == 0) write (6, '(A,E15.7)') 'Final CFL:', ctarg

   lastep = 0
   param(10) = period
   param(12) = -abs(dt)
   param(14) = timeio
   param(15) = 0.0d0

   call bcast(param, 200*wdsize)
   time = 0.0d0
   do istep = 1, nsteps
   call nek_advance
   if (istep >= nsteps) lastep = 1
   call check_ioinfo
   call set_outfld

   call hpts
   call nekStab_torque('lift_drag.dat')
   if (ifoutfld) then

   if (nid == 0) then
   write (6, *) 'mywrite counter, time, dt, fintim:', counter, time, dt, fintim
   write (6, *) 'mywrite counter*timeio:', counter*timeio, counter*timeio - time
   write (6, *) 'mywrite istep*nsteps:', istep*dt, istep*dt - time
   end if
   counter = counter + 1
   call compute_omegaR(vx, vy, vz, t(1, 1, 1, 1, 1))
   ifto = .true.

   call k_copy(Re_cos, Re)
   call k_axpby(Re_cos, cos(omega*time), Im, -sin(omega*time))
   call compute_omegaR(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%t(:, 1))
   call outpost(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%pr, Re_cos%t, trim(mode_local)//'Qm')
   amplitude = A0*exp(sigma*time)
   if (nid == 0) write (6, *) 'Amplitude: ', A0, amplitude
   call k_cmult(Re_cos, amplitude)
   ! deform rides on the *instantaneous* (evolving) base flow, not a static
   ! snapshot: copy the current DNS state into BF before superposing the mode
   call nopcopy(BF%vx, BF%vy, BF%vz, BF%pr, BF%t, vx, vy, vz, pr, t)
   call k_add2(Re_cos, BF)
   call compute_omegaR(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%t(:, 1))
   call outpost(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%pr, Re_cos%t, trim(mode_local)//'Qb')
   end if

   call prepost(ifoutfld, 'his')
   call in_situ_check()
   if (lastep == 1) exit
   end do

   else

   call k_load(BF, 'BF_'//trim(SESSION)//'0.f00001')
   call compute_omegaR(BF%vx, BF%vy, BF%vz, BF%t(:, 1))
   ifto = .true.
   call outpost(BF%vx, BF%vy, BF%vz, BF%pr, BF%t, 'BF_')

   A0 = 1.0e-3
   sigma = 1.0e-1 ! force a value of sigma

   call load_mode_pair(mode, Re, Im)

! Loop over num_of_files to create snapshots
   do i = 1, num_of_files

   time = i*(param(10)/num_of_files)

   if (nid == 0) then
   write (6, '(A, I0, A, F8.4, A, I0)') 'i: ', i, ', time: ', time, ', n: ', num_of_files
   write (6, '(A, F8.4, A, F8.4, A, F8.4, A)')&
      'Time/param(10): ', time/param(10),&
      ' (', time,&
      ' / ', param(10),&
      ') [time/period]'
   end if

   ifto = .true.

   call k_copy(Re_cos, Re)
   call k_axpby(Re_cos, cos(omega*time), Im, -sin(omega*time))

   call compute_omegaR(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%t(:, 1))
   call outpost(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%pr, Re_cos%t, trim(mode_local)//'Qm')

   amplitude = A0*exp(sigma*time)
   if (nid == 0) write (6, *) 'Amplitude: ', A0, amplitude

   call k_cmult(Re_cos, amplitude)
   call k_add2(Re_cos, BF)
   call compute_omegaR(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%t(:, 1))
   call outpost(Re_cos%vx, Re_cos%vy, Re_cos%vz, Re_cos%pr, Re_cos%t, trim(mode_local)//'Qb')

   end do

   end if

   end subroutine animate_mode

   end module nekstab_sensitivity
