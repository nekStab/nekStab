!-----------------------------------------------------------------------
      subroutine wave_maker

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt  = lx1*ly1*lz1*lelt

      real, dimension(lt) :: vx_dRe, vy_dRe, vz_dRe
      real, dimension(lt) :: vx_dIm, vy_dIm, vz_dIm
      real, dimension(lt) :: vx_aRe, vy_aRe, vz_aRe
      real, dimension(lt) :: vx_aIm, vy_aIm, vz_aIm
      real, dimension(lt) :: vx_aRet, vy_aRet, vz_aRet
      real, dimension(lt) :: vx_aImt, vy_aImt, vz_aImt

      real, dimension(lt) :: vx_dw,vy_dw,vz_dw
      real, dimension(lt) :: vx_aw,vy_aw,vz_aw
      real, dimension(lt) :: vx_wm,vy_wm,vz_wm

      character(len=80)   :: filename

      real :: alpha, beta, gamma, delta, epsilon, zeta, eta, theta, iota
      alpha   = 0.0d0
      beta    = 0.0d0
      gamma   = 0.0d0
      delta   = 0.0d0
      epsilon = 0.0d0
      zeta    = 0.0d0
      eta     = 0.0d0
      theta   = 0.0d0
      iota    = 0.0d0

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
      write(filename,'(a,a,a)')'aRe',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_aRe,vy_aRe,vz_aRe,vx,vy,vz)

!     load imaginary part of the adjoint mode
      write(filename,'(a,a,a)')'aIm',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_aIm,vy_aIm,vz_aIm,vx,vy,vz)

!     normalisation of direct mode
      call inner_product(alpha  ,vx_dRe, vy_dRe, vz_dRe, pr, t, vx_dRe, vy_dRe, vz_dRe, pr, t)
      call inner_product(beta   ,vx_dIm, vy_dIm, vz_dIm, pr, t, vx_dIm, vy_dIm, vz_dIm, pr, t)
      alpha   = 1.0d0/dsqrt(alpha+beta)
      call opcmult(vx_dRe,vy_dRe,vz_dRe, alpha)
      call opcmult(vx_dIm,vy_dIm,vz_dIm, alpha)

!     normalisation of adjoint mode
      call inner_product(theta  ,vx_aRe, vy_aRe, vz_aRe, pr, t, vx_aRe, vy_aRe, vz_aRe, pr, t)
      call inner_product(iota   ,vx_aIm, vy_aIm, vz_aIm, pr, t, vx_aIm, vy_aIm, vz_aIm, pr, t)
      theta   = 1.0d0/dsqrt(theta+iota)
      call opcmult(vx_aRe,vy_aRe,vz_aRe, theta)
      call opcmult(vx_aIm,vy_aIm,vz_aIm, theta)

!     biorthonormality condition
      call inner_product(gamma  ,vx_aRe, vy_aRe, vz_aRe, pr, t, vx_dRe, vy_dRe, vz_dRe, pr, t)
      call inner_product(delta  ,vx_aIm, vy_aIm, vz_aIm, pr, t, vx_dIm, vy_dIm, vz_dIm, pr, t)
      gamma   = gamma+delta
      call inner_product(epsilon,vx_aRe, vy_aRe, vz_aRe, pr, t, vx_dIm, vy_dIm, vz_dIm, pr, t)
      call inner_product(zeta   ,vx_aIm, vy_aIm, vz_aIm, pr, t, vx_dRe, vy_dRe, vz_dRe, pr, t)
      epsilon = epsilon-zeta
      eta     = gamma*gamma + epsilon*epsilon
      gamma   = gamma/eta
      epsilon = epsilon/eta
      call opcopy(vx_aRet,vy_aRet,vz_aRet,vx_aRe,vy_aRe,vz_aRe)
      call opcopy(vx_aImt,vy_aImt,vz_aImt,vx_aIm,vy_aIm,vz_aIm)
      call opcmult(vx_aRe,vy_aRe,vz_aRe, gamma)
      call opcmult(vx_aImt,vy_aImt,vz_aImt, epsilon)
      call opcmult(vx_aRet,vy_aRet,vz_aRet,-epsilon)
      call opcmult(vx_aIm,vy_aIm,vz_aIm, gamma)
      call opadd2(vx_aRe,vy_aRe,vz_aRe,vx_aImt,vy_aImt,vz_aImt)
      call opadd2(vx_aIm,vy_aIm,vz_aIm,vx_aRet,vy_aRet,vz_aRet)

!     wave maker computation
      call oprzero(vx_wm,vy_wm,vz_wm)
      call oprzero(vx_dw,vy_dw,vz_dw)
      call opaddcol3(vx_dw,vy_dw,vz_dw, vx_dRe,vy_dRe,vz_dRe, vx_dRe,vy_dRe,vz_dRe)
      call opaddcol3(vx_dw,vy_dw,vz_dw, vx_dIm,vy_dIm,vz_dIm, vx_dIm,vy_dIm,vz_dIm)
      call add4(vx_wm,vx_dw,vy_dw,vz_dw,lt)
      call vsqrt(vx_wm,lt)
      call oprzero(vx_aw,vy_aw,vz_aw)
      call opaddcol3(vx_aw,vy_aw,vz_aw, vx_aRe,vy_aRe,vz_aRe, vx_aRe,vy_aRe,vz_aRe)
      call opaddcol3(vx_aw,vy_aw,vz_aw, vx_aIm,vy_aIm,vz_aIm, vx_aIm,vy_aIm,vz_aIm)
      call add4(vy_wm,vx_aw,vy_aw,vz_aw,lt)
      call vsqrt(vy_wm,lt)
      call addcol3(vz_wm,vx_wm,vy_wm,lt)
      call filter_s0(vz_wm,0.5,1,'vortx')
      ifvo=.false.; ifpo=.false.; ifto=.true.
      call outpost(vz_wm, vz_dw, vz_aw, pr, vz_wm, 'wm_')
      ifvo=.true.; ifpo=.true.; ifto=.true.

      return
      end subroutine wave_maker
!-----------------------------------------------------------------------
      subroutine bf_sensitivity

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt  = lx1*ly1*lz1*lelt

      real, dimension(lt) :: vx_dRe, vy_dRe, vz_dRe
      real, dimension(lt) :: vx_dIm, vy_dIm, vz_dIm
      real, dimension(lt) :: vx_aRe, vy_aRe, vz_aRe
      real, dimension(lt) :: vx_aIm, vy_aIm, vz_aIm
      real, dimension(lt) :: vx_aRet, vy_aRet, vz_aRet
      real, dimension(lt) :: vx_aImt, vy_aImt, vz_aImt

      real, dimension(lt) :: dudx_dRe, dudy_dRe, dudz_dRe
      real, dimension(lt) :: dvdx_dRe, dvdy_dRe, dvdz_dRe
      real, dimension(lt) :: dwdx_dRe, dwdy_dRe, dwdz_dRe

      real, dimension(lt) :: dudx_dIm, dudy_dIm, dudz_dIm
      real, dimension(lt) :: dvdx_dIm, dvdy_dIm, dvdz_dIm
      real, dimension(lt) :: dwdx_dIm, dwdy_dIm, dwdz_dIm

      real, dimension(lt) :: dudx_aRe, dudy_aRe, dudz_aRe
      real, dimension(lt) :: dvdx_aRe, dvdy_aRe, dvdz_aRe
      real, dimension(lt) :: dwdx_aRe, dwdy_aRe, dwdz_aRe

      real, dimension(lt) :: dudx_aIm, dudy_aIm, dudz_aIm
      real, dimension(lt) :: dvdx_aIm, dvdy_aIm, dvdz_aIm
      real, dimension(lt) :: dwdx_aIm, dwdy_aIm, dwdz_aIm

      real, dimension(lt) :: vx_tr,vy_tr,vz_tr
      real, dimension(lt) :: vx_ti,vy_ti,vz_ti
      real, dimension(lt) :: vx_pr,vy_pr,vz_pr
      real, dimension(lt) :: vx_pi,vy_pi,vz_pi
      real, dimension(lt) :: vx_sr,vy_sr,vz_sr

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
      write(filename,'(a,a,a)')'aRe',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_aRe,vy_aRe,vz_aRe,vx,vy,vz)

!     load imaginary part of the adjoint mode
      write(filename,'(a,a,a)')'aIm',trim(SESSION),'0.f00001'
      call load_fld(filename)
      call opcopy(vx_aIm,vy_aIm,vz_aIm,vx,vy,vz)

!     normalisation of direct mode
      call inner_product(alpha  ,vx_dRe, vy_dRe, vz_dRe, pr, t, vx_dRe, vy_dRe, vz_dRe, pr, t)
      call inner_product(beta   ,vx_dIm, vy_dIm, vz_dIm, pr, t, vx_dIm, vy_dIm, vz_dIm, pr, t)
      alpha   = 1.0d0/dsqrt(alpha+beta)
      call opcmult(vx_dRe,vy_dRe,vz_dRe, alpha)
      call opcmult(vx_dIm,vy_dIm,vz_dIm, alpha)
!     call outpost(vx_dRe,vy_dRe,vz_dRe, pr, t, 'dr_')
!     call outpost(vx_dIm,vy_dIm,vz_dIm, pr, t, 'di_')

!     normalisation of adjoint mode
      call inner_product(theta  ,vx_aRe, vy_aRe, vz_aRe, pr, t, vx_aRe, vy_aRe, vz_aRe, pr, t)
      call inner_product(iota   ,vx_aIm, vy_aIm, vz_aIm, pr, t, vx_aIm, vy_aIm, vz_aIm, pr, t)
      theta   = 1.0d0/dsqrt(theta+iota)
      call opcmult(vx_aRe,vy_aRe,vz_aRe, theta)
      call opcmult(vx_aIm,vy_aIm,vz_aIm, theta)

!     biorthonormality condition
      call inner_product(gamma  ,vx_aRe, vy_aRe, vz_aRe, pr, t, vx_dRe, vy_dRe, vz_dRe, pr, t)
      call inner_product(delta  ,vx_aIm, vy_aIm, vz_aIm, pr, t, vx_dIm, vy_dIm, vz_dIm, pr, t)
      gamma   = gamma+delta
      call inner_product(epsilon,vx_aRe, vy_aRe, vz_aRe, pr, t, vx_dIm, vy_dIm, vz_dIm, pr, t)
      call inner_product(zeta   ,vx_aIm, vy_aIm, vz_aIm, pr, t, vx_dRe, vy_dRe, vz_dRe, pr, t)
      epsilon = epsilon-zeta
      eta     = gamma*gamma + epsilon*epsilon
      gamma   = gamma/eta
      epsilon = epsilon/eta
      call opcopy(vx_aRet,vy_aRet,vz_aRet,vx_aRe,vy_aRe,vz_aRe)
      call opcopy(vx_aImt,vy_aImt,vz_aImt,vx_aIm,vy_aIm,vz_aIm)
      call opcmult(vx_aRe,vy_aRe,vz_aRe, gamma)
      call opcmult(vx_aImt,vy_aImt,vz_aImt, epsilon)
      call opcmult(vx_aRet,vy_aRet,vz_aRet,-epsilon)
      call opcmult(vx_aIm,vy_aIm,vz_aIm, gamma)
      call opadd2(vx_aRe,vy_aRe,vz_aRe,vx_aImt,vy_aImt,vz_aImt)
      call opadd2(vx_aIm,vy_aIm,vz_aIm,vx_aRet,vy_aRet,vz_aRet)

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
      call filter_s0(vx_tr,0.5,1,'vortx')
      call filter_s0(vy_tr,0.5,1,'vortx')
      call filter_s0(vz_tr,0.5,1,'vortx')
      call outpost(vx_tr, vy_tr, vz_tr, pr, t, 'tr_')

      call filter_s0(vx_ti,0.5,1,'vortx')
      call filter_s0(vy_ti,0.5,1,'vortx')
      call filter_s0(vz_ti,0.5,1,'vortx')
      call outpost(vx_ti, vy_ti, vz_ti, pr, t, 'ti_')

      call filter_s0(vx_pr,0.5,1,'vortx')
      call filter_s0(vy_pr,0.5,1,'vortx')
      call filter_s0(vz_pr,0.5,1,'vortx')
      call outpost(vx_pr, vy_pr, vz_pr, pr, t, 'pr_')

      call filter_s0(vx_pi,0.5,1,'vortx')
      call filter_s0(vy_pi,0.5,1,'vortx')
      call filter_s0(vz_pi,0.5,1,'vortx')
      call outpost(vx_pi, vy_pi, vz_pi, pr, t, 'pi_')

      call opadd2(vx_tr,vy_tr,vz_tr, vx_pr,vy_pr,vz_pr)
      call opadd2(vx_ti,vy_ti,vz_ti, vx_pi,vy_pi,vz_pi)

      call outpost(vx_tr,vy_tr,vz_tr, pr, t, 'sr_')
      call opcmult(vx_ti,vy_ti,vz_ti, kappa)
      call outpost(vx_ti,vy_ti,vz_ti, pr, t, 'si_')

      return
      end subroutine bf_sensitivity
