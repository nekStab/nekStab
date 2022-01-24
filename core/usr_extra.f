c-----------------------------------------------------------------------
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

      bst_skp = 10              ! boostconv skip iterations
      bst_snp = 10              ! bootsconv residual subspace matrix size

      ifres  = .false.          ! outpost restart files (KRY*, HES*)
      ifvor  = .false.          ! outpost vorticity (vor* omega_x,omega_y,omega_z components)
      ifvox  = .false.          ! outpost vortex (vox*: q,lambda2,omega criterions)
      ifldbf = .true.           ! load base flow for stability computations
      ifbf2D = .false.          ! force 2D base flow solution
      ifstorebase = .true.      ! store base flow for Floquet analysis (dynamic allocated)
      ifdyntol = .false.        ! dynamical tolerances for SFD and Newton (potential speed-up)

      ifseed_nois = .true.      ! noise as initial seed
      ifseed_symm = .false.     ! symmetry initial seed
      ifseed_load = .false.     ! loading initial seed (e.g. Re_ )
!     if all false the 'useric' subroutine prescribes the initial seed

      xLspg   = 0.0d0; call bcast(xLspg, wdsize) ! x left
      xRspg   = 0.0d0; call bcast(xRspg, wdsize) ! x right
      yLspg   = 0.0d0; call bcast(yLspg, wdsize)
      yRspg   = 0.0d0; call bcast(yRspg, wdsize)
      zLspg   = 0.0d0; call bcast(zLspg, wdsize)
      zRspg   = 0.0d0; call bcast(zRspg, wdsize)
!     percentage for the acceleration phase in the sponge (e.g. 1/3)
      acc_spg = 0.333d0; call bcast(acc_spg, wdsize)

!     !Broadcast all defaults !
      call bcast(schur_tgt, isize) ! isize for integer
      call bcast(eigen_tol, wdsize) ! wdsize for real
      call bcast(schur_del, wdsize)
      call bcast(maxmodes, isize)
      call bcast(k_dim, isize)
      call bcast(bst_skp, isize)
      call bcast(bst_snp, isize)
      call bcast(glob_skip, isize)

      call bcast(ifres   , lsize) !lsize for boolean
      call bcast(ifvor   , lsize)
      call bcast(ifvox   , lsize)
      call bcast(ifseed_nois  , lsize)
      call bcast(ifseed_symm  , lsize)
      call bcast(ifseed_load  , lsize)
      call bcast(ifldbf  , lsize)
      call bcast(ifbf2D  , lsize)
      call bcast(ifstorebase  , lsize)
      call bcast(ifdyntol  , lsize)

      return
      end subroutine nekStab_setDefault
c-----------------------------------------------------------------------
      subroutine nekStab_init
      use krylov_subspace
!     initialize arrays and variables defaults
      implicit none
      include 'SIZE'
      include 'TOTAL'
      logical ifto_sav, ifpo_sav
      real glmin,glmax
      integer i
      n = nx1*ny1*nz1*nelv

      call nekStab_setDefault
      call nekStab_usrchk       ! where user change defaults
      call nekStab_printNEKParams

      xmn = glmin(xm1,n); xmx = glmax(xm1,n)
      ymn = glmin(ym1,n); ymx = glmax(ym1,n)
      zmn = glmin(zm1,n); zmx = glmax(zm1,n)

      if (nid==0) then
         print *,'                 __   _____  __          __  '
         print *,'   ____   ___   / /__/ ___/ / /_ ____ _ / /_ '
         print *,'  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \'
         print *,' / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /'
         print *,'/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/ '  
         print *,'COPYRIGHT (c) 2020-2022 DynFluid Laboratoire Paris ',NSVERSION
         print *,'Nek5000 ', NVERSION
         print *,''
!     print *,'==================================================='
!     print *,'== nekStab ========================================'
!     print *,'== Copyright (c) 2021 DynFluid Laboratoire ========'
!     print *,'==================================================='
!     write(6,'(A,A)')' Nek5000 version        : ', NVERSION
!     write(6,'(A,A)')' NekStab version        : ', NSVERSION
      endif

      call copy(bm1s, bm1, n)   ! never comment this !

      if(uparam(10).gt.0)then   !sponge on

         if(nid.eq.0)write(6,*)
         if(nid.eq.0)write(6,*)' Initializing sponge...'
         if(nid.eq.0)write(6,*)

         spng_str = uparam(10)
         call spng_init

!     applying sponge to the BM1 matrix to remove the sponge zone from eigensolver
         do i=1,n
            if( spng_fun( i ) .gt. 0 ) bm1s( i,1,1,1 )=0.0d0
         enddo

!     outposting BM1s to disk for check
!     ifto_sav = ifto; ifpo_sav = ifpo
!     ifvo=.false.; ifpo = .false.; ifto = .true.
!     call outpost(vx,vy,vz,pr,bm1s,'BMS')
!     ifvo=.true.; ifpo = ifpo_sav; ifto = ifto_sav
      endif
      ifbfcv = .false.

      return
      end subroutine nekStab_init
c-----------------------------------------------------------------------
      subroutine nekStab
!     nekStab main driver
      implicit none
      include 'SIZE'
      include 'TOTAL'

      if(istep.eq.0)call nekStab_init

      call oprzero(fcx,fcy,fcz) ! never comment this!
      call rzero(fct,nx1*ny1*nz1*nelv)

      select case (floor(uparam(1)))

      case(0)                   ! DNS

         call nekStab_outpost   ! outpost vorticity
         call nekStab_comment   ! print comments
         call nekStab_energy   (vx,vy,vz,t,'global_energy.dat'   ,glob_skip)
         call nekStab_enstrophy(vx,vy,vz,t,'global_enstrophy.dat',glob_skip)

      case(1)                   ! fixed points computation

         call nekStab_outpost   ! outpost vorticity
         call nekStab_comment   ! print comments

         if( uparam(1) .eq. 1.1)then
            if(nid.eq.0)write(6,*)'SFD'
            call SFD
         elseif( uparam(1) .eq. 1.2)then
            if(nid.eq.0)write(6,*)'BOOSTCONV'
            call BoostConv
         elseif( uparam(1) .eq. 1.3)then
            if(nid.eq.0)write(6,*)'DMT'
            if(nid.eq.0)write(6,*)'stopping ! not yet ported to this version'; call nek_end
            ! call DMT
         elseif( uparam(1) .eq. 1.4)then
            if(nid.eq.0)write(6,*)'TDF'
            call TDF
         endif

         if(ifbfcv)call nek_end

      case(2)                   ! Newton-Krylov solver

!     sanity check
         if(uparam(1).eq.2.0)then
            if(nid.eq.0)write(6,*)'Newton-Krylov for fixed points...'
         elseif(uparam(1).eq.2.1)then
            if(nid.eq.0)write(6,*)'Newton-Krylov for UPOs...'
         elseif(uparam(1).eq.2.2)then
            if(nid.eq.0)write(6,*)'Newton-Krylov for forced UPOs...'
         else
            if(nid.eq.0)write(6,*)'NEWTON MODE NOT CORRECTLY SPECIFIED!'
            call nek_end
         endif
         call newton_krylov
         call nek_end

      case(3)                   ! eigenvalue problem

!     sanity check
         if    (uparam(1).eq.3.1)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for direct LNSE...'
         elseif(uparam(1).eq.3.11)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for direct LNSE in Floquet...'
         elseif(uparam(1).eq.3.2)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for adjoint LNSE...'
         elseif(uparam(1).eq.3.21)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for adjoint LNSE in Floquet...'           
         elseif(uparam(1).eq.3.3)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for transient growth...'
         elseif(uparam(1).eq.3.31)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for transient growth in Floquet...'
         else
            if(nid.eq.0)write(6,*)'Krylov-Schur MODE NOT CORRECTLY SPECIFIED!'
            call nek_end
         endif
         call krylov_schur
         call nek_end

      case(4)                   ! in postprocessing.f

         if(uparam(01) .eq. 4.0) then ! all
            call stability_energy_budget
            call wave_maker
            call bf_sensitivity
         endif
!     -----> Direct mode kinetic energy budget.
         if(uparam(01) .eq. 4.1) call stability_energy_budget

!     -----> Wavemaker computation.
         if(uparam(01) .eq. 4.2) call wave_maker

!     -----> Baseflow sensitivity.
         if(uparam(01) .eq. 4.3) call bf_sensitivity

!     -----> Sensitivity to steady force.
         if(uparam(01) .eq. 4.4) call ts_steady_force_sensitivity

         call nek_end

      end select

      return
      end subroutine nekStab
c-----------------------------------------------------------------------
      subroutine nekStab_outpost
!     nekStab custom outpost routine
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real vort(lv,3),wo1(lv),wo2(lv)
      common /ugrad/ vort,wo1,wo2

      logical ifto_sav, ifpo_sav

      if((istep.eq.0).OR.ifoutfld.AND.(ifvor.or.ifvox))then

         ifto_sav = ifto; ifpo_sav = ifpo
         ifto = .false.; ifpo = .false.

!---  > Compute and oupost vorticity.
         if(ifvor)then
            call oprzero(vort(1,1),vort(1,2),vort(1,3))
            call comp_vort3(vort,wo1,wo2,vx,vy,vz)
            call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t, 'vor')
         endif

!---  > Compute and outpost vortex fields.
         if(ifvox.and.(ifaxis.eqv..false.))then
            call vortex_core(vort(:,1),'lambda2')
            call vortex_core(vort(:,2),'q')
            call vortex_core(vort(:,3),'omega')
            call outpost(vort(:,1),vort(:,2),vort(:,3),pr,t,'vox')

            if(.not.if3d)then
               ifvo=.false.; ifto = .true. ! just outposting one field ... v's and p ignored
               call outpost(vx,vy,vz,pr,vort(:,3),'omg')
               ifvo=.true.; ifto = ifto_sav
            endif

         endif

         ifto = ifto_sav ; ifpo = ifpo_sav

      endif
!---  > Outpost initial condition.
      if(istep.eq.0)call outpost(vx,vy,vz,pr,t,'   ')

      return
      end subroutine nekStab_outpost
c-----------------------------------------------------------------------
      subroutine nekStab_comment
!     Comment the evoltuion of the simulation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real*8,save :: eetime0,eetime1,eetime2
      data           eetime0,eetime1,eetime2 /0.0d0, 0.0d0, 0.0d0/
      real, save :: deltatime
      real telapsed,tpernondt,tmiss,dnekclock,ttime
      integer ttime_stp

!     if extrapolation is not OIFS -> ifchar = false
!     if OIFS active -> ifchar = .true. and CFL 2-5
!     some cases can have CFL>1 in initial time steps
      if (courno.gt.10) then
         if (nio.eq.0)then
            write(6,*)
            write(6,*)'    CFL > 10 stopping code'
            write(6,*)
         endif
         call nek_end
      endif

      if (nio.ne.0) return

      if (eetime0.eq.0.0 .and. istep.eq.1)then
         eetime0=dnekclock()
         deltatime=time
      endif
      eetime1=eetime2
      eetime2=dnekclock()

      if (istep.gt.0 .and. lastep.eq.0 .and. iftran) then

         ttime_stp = eetime2-eetime1 ! time per timestep
         ttime     = eetime2-eetime0 ! sum of all timesteps

         if(istep.eq.1)then
            ttime_stp = 0.0d0; ttime = 0.0d0
         endif

         if (mod(istep,5).eq.0) then

            telapsed = ttime/3600.0d0
            tpernondt = (ttime/(time-deltatime))
            tmiss = (param(10)-time)*tpernondt/3600.0d0

            print *,''
            write(6,"('      Mean time per timestep: ',F8.4,'  dev:',I8,'ms')") ttime/istep,ceiling(((ttime/istep)-ttime_stp)*1000) !to ms
            write(6,"('      Remaining time: ',I8,' h ',I2,' min')") int(tmiss),ceiling((tmiss-int(tmiss))*60.)
            if(tpernondt.gt.60.)then
               write(6,"('      Time per nondimensional time: ',F8.2,' sec')") tpernondt
            else
               write(6,"('      Time per nondimensional time: ',F8.2,' min ')") tpernondt/60.0d0
            endif
            write(6,"('      Local time: ',F8.4,'  File:',I8)") time-deltatime, int((time-deltatime)/param(14))+1
            print *,''

         endif
      endif

      return
      end subroutine nekStab_comment
c-----------------------------------------------------------------------
      subroutine nekStab_printNEKParams
!     print parameters at initialization for sanity check
      implicit none
      include 'SIZE'
      include 'TOTAL'
      if(nid.eq.0)then
         write(6,*)'P01=',param(1),'density'
         write(6,*)'P02=',param(2),'viscosity (1/Re)'
         write(6,*)'P07=',param(7),'rhoCp'
         write(6,*)'P08=',param(8),'conductivity (1/(Re*Pr))'
         write(6,*)'P10=',param(10),'stop at endTime'
         write(6,*)'P10=',param(11),'stop at numSteps'
         write(6,*)'P14=',param(14),'io step'
         write(6,*)'P15=',param(15),'io time'
         write(6,*)'P21=',param(21),'pressure sol tol'
         write(6,*)'P22=',param(22),'velocity sol tol'
         write(6,*)'P26=',param(26),'target CFL number'
         write(6,*)'P27=',param(27),'order in time'
         write(6,*)'P28=',param(28),'use same torder for mesh solver'
         write(6,*)'P31=',param(31),'numberOfPerturbations'
         write(6,*)'P41=',param(41),'1 for multiplicative SEMG'
         write(6,*)'P42=',param(42),'lin solv for the pres equation 0:GMRES,1:CG'
         write(6,*)'P43=',param(43),'0:additive multilevel scheme 1:orig 2lvl sch'
         write(6,*)'P44=',param(44),'0=E-based addit Schwarz PnPn-2;1=A-based'
         write(6,*)'P93=',param(93),'num vectors for projection'
         write(6,*)'P94 =',param(94),'num projection for helmholz solves'
         write(6,*)'P95=',param(95),'projection for pressure solver on/off'
         write(6,*)'P101=',param(101),'no additional modes'
         write(6,*)'P103=',param(103),'filter weight'
         write(6,*)
         write(6,*)'uparam01=',uparam(1)
         write(6,*)'uparam02=',uparam(02)
         write(6,*)'uparam03=',uparam(03)
         write(6,*)'uparam04=',uparam(04)
         write(6,*)'uparam05=',uparam(05)
         write(6,*)'uparam06=',uparam(06)
         write(6,*)'uparam07=',uparam(07)
         write(6,*)'uparam08=',uparam(08)
         write(6,*)'uparam09=',uparam(09)
         write(6,*)'uparam10=',uparam(10)
         write(6,*)
         write(6,*)'x min,max,tot=',xmn,xmx,xmx-xmn
         write(6,*)'y min,max,tot=',ymn,ymx,ymx-xmn
         write(6,*)'z min,max,tot=',zmn,zmx,zmx-zmn
         write(6,*)
      endif
      end subroutine nekStab_printNEKParams
c-----------------------------------------------------------------------
      subroutine nekStab_energy(px, py, pz, pt, fname, skip)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real, dimension(lv), intent(in) :: px, py, pz, pt
      integer, intent(in) :: skip
      character(len=20), intent(in) :: fname
      real glsc3,uek,vek,wek,eek,pot
      save eek
      logical, save :: initialized
      data             initialized /.FALSE./
      n = nx1*ny1*nz1*nelv
      eek = 0.50d0/volvm1
      uek = 0.0d0; vek = 0.0d0; wek = 0.0d0;  pot = 0.0d0

      if (.not. initialized) then
         if(nid.eq.0)open (730,file=fname,action='write',status='replace')
         initialized = .true.
      endif

      if (mod(istep,skip).eq.0) then
         uek = glsc3(px,bm1,px,n)
         vek = glsc3(py,bm1,py,n)
         if(if3d)wek = glsc3(pz,bm1,pz,n)
         if(ifheat)pot = glsc3(pt,bm1,ym1,n)
         if(nid.eq.0)write(730,"(6E15.7)")time,uek*eek,vek*eek,wek*eek,(uek+vek+wek)*eek,pot*eek
      endif

      return
      end subroutine nekStab_energy
c-----------------------------------------------------------------------
      subroutine nekStab_enstrophy(px, py, pz, pt, fname, skip)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real, dimension(lv), intent(in) :: px, py, pz, pt
      integer, intent(in) :: skip
      character(len=20), intent(in) :: fname
      real vort(lv,3),wo1(lv),wo2(lv)
      common /ugrad/ vort,wo1,wo2
      real glsc3,uek,vek,wek,eek
      save eek
      logical, save :: initialized
      data             initialized /.FALSE./
      n = nx1*ny1*nz1*nelv
      eek = 0.50d0/volvm1
      uek = 0.0d0; vek = 0.0d0; wek = 0.0d0

      if (.not. initialized) then
         if(nid.eq.0)open (736,file=fname,action='write',status='replace')
         initialized = .true.
      endif

      if (mod(istep,skip).eq.0) then

         call comp_vort3(vort,wo1,wo2,px,py,pz)

         uek = glsc3(vort(1,1),bm1,vort(1,1),n)
         vek = glsc3(vort(1,2),bm1,vort(1,2),n)
         if(if3d)wek = glsc3(vort(1,3),bm1,vort(1,3),n)
         if(nid.eq.0)write(736,"(5E15.7)")time,uek*eek,vek*eek,wek*eek,(uek+vek+wek)*eek

      endif

      return
      end subroutine nekStab_enstrophy
c-----------------------------------------------------------------------