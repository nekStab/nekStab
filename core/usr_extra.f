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

!     position for zero-crossing vertical velocity check ! 
      xck = 2.0D0  ; call bcast(xck, wdsize)
      yck = 0.0D0  ; call bcast(yck, wdsize)
      zck = 0.0D0  ; call bcast(zck, wdsize)

      xLspg   = 0.0d0; call bcast(xLspg, wdsize) ! x left
      xRspg   = 0.0d0; call bcast(xRspg, wdsize) ! x right
      yLspg   = 0.0d0; call bcast(yLspg, wdsize)
      yRspg   = 0.0d0; call bcast(yRspg, wdsize)
      zLspg   = 0.0d0; call bcast(zLspg, wdsize)
      zRspg   = 0.0d0; call bcast(zRspg, wdsize)
      acc_spg = 0.333d0; call bcast(acc_spg, wdsize) !percentage for the acceleration phase in the sponge (e.g. 1/3)
      spng_str = 0.0d0;  call bcast(spng_str, wdsize)

      evop = '_'

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
         print *,'COPYRIGHT (c) 2020-2023 DynFluid Laboratoire Paris ',NSVERSION
         print *,'Nek5000 ', NVERSION
         print *,''
      endif

      call copy(bm1s, bm1, n)   ! never comment this !

      if(spng_str.ne.0)then     !sponge on

         if(nid.eq.0)write(6,*)
         if(nid.eq.0)write(6,*)' Initializing sponge...'
         if(nid.eq.0)write(6,*)' Sponge strenght:',spng_str
         if(spng_str.lt.0)then
            spng_str=abs(spng_str) 
            if(nid.eq.0)write(6,*)' Ensure positive sponge strenght:',spng_str
         endif
         call spng_init

!     applying sponge to the BM1 matrix to remove the sponge zone from eigensolver
         do i=1,n
            if( spng_fun( i ) .ne. 0 ) bm1s( i,1,1,1 )=0.0d0
         enddo

!     outposting BM1s to disk for check
!     ifto_sav = ifto; ifpo_sav = ifpo
!     ifvo=.false.; ifpo = .false.; ifto = .true.
!     call outpost(vx,vy,vz,pr,bm1s,'BMS')
!     ifvo=.true.; ifpo = ifpo_sav; ifto = ifto_sav

         if(nid.eq.0)write(6,*)'Sponge activated.'
         if(nid.eq.0)write(6,*)
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
         call nekStab_energy   (vx,vy,vz,t,'total_energy.dat',glob_skip)
         call nekStab_enstrophy(vx,vy,vz,t,'total_enstrophy.dat',glob_skip)

      case(1)                   ! fixed points computation

         call nekStab_outpost   ! outpost vorticity
         call nekStab_comment   ! print comments

         if( uparam(1) .eq. 1.1)then
            if(nid.eq.0)write(6,*)'SFD'
            call SFD
            if(uparam(5).eq.0)call nekStab_energy(vx,vy,vz,t,'total_energy.dat',glob_skip)
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
         if(uparam(01) .eq. 4.41 .or. 
     $        uparam(01) .eq. 4.42) call ts_steady_force_sensitivity
         if(uparam(01) .eq. 4.43) call delta_forcing

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
            write(6,"('      Local time: ',F10.4,'  File:',I8)") time-deltatime, int((time-deltatime)/param(14))+1
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
      real, dimension(lv), intent(in) :: px, py, pz
      real, dimension(lv,ldimt), intent(in) :: pt
!     real, dimension(lx1,ly1,lz1,lelt,ldimt), intent(in) :: pt
      integer, intent(in) :: skip
      character(len=16), intent(in) :: fname
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
         uek = glsc3(px,bm1s,px,n)
         vek = glsc3(py,bm1s,py,n)
         if(if3d)wek = glsc3(pz,bm1s,pz,n)
         if(ifheat)pot = glsc3(pt(:,1),bm1s,ym1,n)
!     if(ifheat)pot = glsc3(pt(:,:,:,:,1),bm1s,ym1,n)
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
      real, dimension(lv), intent(in) :: px, py, pz
      real, dimension(lv,ldimt), intent(in) :: pt
!     real, dimension(lx1,ly1,lz1,lelt,ldimt), intent(in) :: pt
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

         uek = glsc3(vort(1,1),bm1s,vort(1,1),n)
         vek = glsc3(vort(1,2),bm1s,vort(1,2),n)
         if(if3d)wek = glsc3(vort(1,3),bm1s,vort(1,3),n)
         if(nid.eq.0)write(736,"(5E15.7)")time,uek*eek,vek*eek,wek*eek,(uek+vek+wek)*eek

      endif

      return
      end subroutine nekStab_enstrophy
c-----------------------------------------------------------------------
      subroutine nekStab_torque(fname, skip)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, intent(in) :: skip
      character(len=20), intent(in) :: fname

      logical, save :: initialized
      data             initialized /.FALSE./

      integer, save :: bIDs(1), iobj_wall(1)

      real, save :: x0(3), scale
      data x0 /3*0/
      
      integer :: nij,i,iobj,memtot,mem,ieg,ifc,ie
      integer, parameter :: lr=lx1*ly1*lz1

      real glmin,glmax,x1min,x2min,x3min,x1max,x2max,x3max,w1(0:maxobj)
      real flow_rate,base_flow,domain_length,xsec,scale_vf(3),sij,pm1,xm0,ym0,zm0

      real ur,us,ut,vr,vs,vt,wr,ws,wt
      common /scruz/ ur(lr),us(lr),ut(lr),vr(lr),vs(lr),vt(lr),wr(lr),ws(lr),wt(lr)
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec,scale_vf
      common /scrns/ sij(lx1*ly1*lz1*6*lelv)
      common /scrcg/ pm1(lx1,ly1,lz1,lelv)
      common /scrsf/ xm0(lx1,ly1,lz1,lelt),ym0(lx1,ly1,lz1,lelt),zm0(lx1,ly1,lz1,lelt)
      n = lx1*ly1*lz1*nelv

      if (.not. initialized) then
         if(nid.eq.0)write(6,*)'Initializing torque routine...'
         if(nid.eq.0)open (737,file=fname,action='write',status='replace')
         bIDs(1) = 1
         call create_obj(iobj_wall(1), bIDs, 1)
         scale = 2
         initialized = .true.
      endif

      if (mod(istep,skip).eq.0) then

         call mappr(pm1,pr,xm0,ym0) ! map pressure onto Mesh 1
         if (param(55).ne.0) then
            dpdx_mean = -scale_vf(1)
            dpdy_mean = -scale_vf(2)
            dpdz_mean = -scale_vf(3)
         endif
         call add2s2(pm1,xm1,dpdx_mean,n) ! Doesn't work if object is cut by 
         call add2s2(pm1,ym1,dpdy_mean,n) ! periodicboundary.  In this case,
         call add2s2(pm1,zm1,dpdz_mean,n) ! set ._mean=0 and compensate in
         nij = 3
         if (if3d.or.ifaxis) nij=6
         call comp_sij(sij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)
         if (istep.lt.1) call cfill(vdiff,param(2),n)
         call cadd2(xm0,xm1,-x0(1),n)
         call cadd2(ym0,ym1,-x0(2),n)
         call cadd2(zm0,zm1,-x0(3),n)
         x1min=glmin(xm0(1,1,1,1),n)
         x2min=glmin(ym0(1,1,1,1),n)
         x3min=glmin(zm0(1,1,1,1),n)
         x1max=glmax(xm0(1,1,1,1),n)
         x2max=glmax(ym0(1,1,1,1),n)
         x3max=glmax(zm0(1,1,1,1),n)
         do i=0,maxobj
            dragpx(i) = 0       ! BIG CODE  :}
            dragvx(i) = 0
            dragx (i) = 0
            dragpy(i) = 0
            dragvy(i) = 0
            dragy (i) = 0
            dragpz(i) = 0
            dragvz(i) = 0
            dragz (i) = 0
            torqpx(i) = 0
            torqvx(i) = 0
            torqx (i) = 0
            torqpy(i) = 0
            torqvy(i) = 0
            torqy (i) = 0
            torqpz(i) = 0
            torqvz(i) = 0
            torqz (i) = 0
         enddo
         ifield = 1
         do iobj = 1,nobj
            memtot = nmember(iobj)
            do mem  = 1,memtot
               ieg   = object(iobj,mem,1)
               ifc   = object(iobj,mem,2)
               if (gllnid(ieg).eq.nid) then ! this processor has a contribution
                  ie = gllel(ieg)
                  call drgtrq(dgtq,xm0,ym0,zm0,sij,pm1,vdiff,ifc,ie)
                  call cmult(dgtq,scale,12)
                  dragpx(iobj) = dragpx(iobj) + dgtq(1,1) ! pressure 
                  dragpy(iobj) = dragpy(iobj) + dgtq(2,1)
                  dragpz(iobj) = dragpz(iobj) + dgtq(3,1)
                  dragvx(iobj) = dragvx(iobj) + dgtq(1,2) ! viscous
                  dragvy(iobj) = dragvy(iobj) + dgtq(2,2)
                  dragvz(iobj) = dragvz(iobj) + dgtq(3,2)
                  torqpx(iobj) = torqpx(iobj) + dgtq(1,3) ! pressure 
                  torqpy(iobj) = torqpy(iobj) + dgtq(2,3)
                  torqpz(iobj) = torqpz(iobj) + dgtq(3,3)
                  torqvx(iobj) = torqvx(iobj) + dgtq(1,4) ! viscous
                  torqvy(iobj) = torqvy(iobj) + dgtq(2,4)
                  torqvz(iobj) = torqvz(iobj) + dgtq(3,4)
               endif
            enddo
         enddo
         call gop(dragpx,w1,'+  ',maxobj+1)
         call gop(dragpy,w1,'+  ',maxobj+1)
         call gop(dragpz,w1,'+  ',maxobj+1)
         call gop(dragvx,w1,'+  ',maxobj+1)
         call gop(dragvy,w1,'+  ',maxobj+1)
         call gop(dragvz,w1,'+  ',maxobj+1)
         call gop(torqpx,w1,'+  ',maxobj+1)
         call gop(torqpy,w1,'+  ',maxobj+1)
         call gop(torqpz,w1,'+  ',maxobj+1)
         call gop(torqvx,w1,'+  ',maxobj+1)
         call gop(torqvy,w1,'+  ',maxobj+1)
         call gop(torqvz,w1,'+  ',maxobj+1)
         do i=1,nobj
            dragx(i) = dragpx(i) + dragvx(i)
            dragy(i) = dragpy(i) + dragvy(i)
            dragz(i) = dragpz(i) + dragvz(i)
            torqx(i) = torqpx(i) + torqvx(i)
            torqy(i) = torqpy(i) + torqvy(i)
            torqz(i) = torqpz(i) + torqvz(i)
            dragpx(0) = dragpx (0) + dragpx (i)
            dragvx(0) = dragvx (0) + dragvx (i)
            dragx (0) = dragx  (0) + dragx  (i)
            dragpy(0) = dragpy (0) + dragpy (i)
            dragvy(0) = dragvy (0) + dragvy (i)
            dragy (0) = dragy  (0) + dragy  (i)
            dragpz(0) = dragpz (0) + dragpz (i)
            dragvz(0) = dragvz (0) + dragvz (i)
            dragz (0) = dragz  (0) + dragz  (i)
            torqpx(0) = torqpx (0) + torqpx (i)
            torqvx(0) = torqvx (0) + torqvx (i)
            torqx (0) = torqx  (0) + torqx  (i)
            torqpy(0) = torqpy (0) + torqpy (i)
            torqvy(0) = torqvy (0) + torqvy (i)
            torqy (0) = torqy  (0) + torqy  (i)
            torqpz(0) = torqpz (0) + torqpz (i)
            torqvz(0) = torqvz (0) + torqvz (i)
            torqz (0) = torqz  (0) + torqz  (i)
         enddo
         do i=1,nobj
            if (nio.eq.0) then
               if (if3d.or.ifaxis) then
                  write(737,"(i8,19E15.7)") istep,time,
     $                 dragx(i),dragpx(i),dragvx(i),dragy(i),dragpy(i),dragvy(i),dragz(i),dragpz(i),dragvz(i),
     $                 torqx(i),torqpx(i),torqvx(i),torqy(i),torqpy(i),torqvy(i),torqz(i),torqpz(i),torqvz(i)
               else
                  write(737,"(i8,10E15.7)") istep,time,
     $                 dragx(i),dragpx(i),dragvx(i),dragy(i),dragpy(i),dragvy(i),torqz(i),torqpz(i),torqvz(i)
               endif
            endif
         enddo
      endif
      return
      end subroutine nekStab_torque
c-----------------------------------------------------------------------
      subroutine nekStab_define_obj
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer iel,ifc

      do iel=1,nelt
         do ifc=1,2*ndim
            if (cbc(ifc,iel,1) .eq. 'W  ') boundaryID(ifc,iel) = 1
         enddo
      enddo

      return
      end subroutine nekStab_define_obj
c-----------------------------------------------------------------------
      subroutine zero_crossing(v_mean_init)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real, intent(in) :: v_mean_init
      real, save :: T_delayed(lv,3),do1(lv),do2(lv),do3(lv),v_mean
      integer, parameter :: plor = 2
      real h1,l2,semi,linf
      save l2
      real               :: glsum, dtime, vdot, vddot
      real, save         :: velp(plor), v_sum, time0
      real, save         :: t_cross, t_cross_old
      real, save         :: v_cross, v_cross_old
      real, save         :: p_now, p_sum, p_mean, p_old
      integer, save      :: probe_nel, probe_nid, t_cross_count
      integer            :: i

      if(istep.eq.0) then
         if(nid.eq.0)write(6,*)'Initializing zero-crossing routine...'
         probe_nel = 0;         probe_nid = 0; vdot=0.0d0; vddot=0.0d0
         velp(:) = 0.0d0;       v_sum = 0.0d0;          v_mean = v_mean_init
         p_now = 0.0d0;         p_sum = 0.0d0;          p_mean = 0.0d0
         t_cross = 0.0d0;       t_cross_old = 0.0d0; t_cross_count = 0
         time0 = time;          p_old = 0.0d0;                l2=0.0d0
         call pointcheck(probe_nel,probe_nid) !alter xck, yck, zck in usrchck
         if(nid.eq.0)open(unit=17,file='zc_period.dat')
         if(nid.eq.0)open(unit=19,file='zc_poincare.dat')
         call opcopy(T_delayed(:,1),T_delayed(:,2),T_delayed(:,3),vx,vy,vz)
      endif

      velp(plor)= 0.0d0
      if(probe_nid.eq.1)velp(plor) = vy(probe_nel,1,1,1)
      velp(plor) = glsum(velp(plor),1) !broadcast
      dtime = time - time0

      if(istep.gt.1)then
         v_sum = v_sum + velp(plor)*dt !(v_mean*(dtime-dt)+v_now*dt)/dtime
         v_mean = v_mean_init + v_sum/dtime
!     if(nid.eq.0)write(6,*)' probe: v_mean, v_now= ',v_mean,velp(plor)
!     if(nid.eq.0)write(6,*)'         v_sum, dtime= ',v_sum,dtime
      endif
      if( velp(plor-1) .le. v_mean .AND. velp(plor) .ge. v_mean ) then !period found

         p_old = p_now          !save old value
         t_cross_old = t_cross  !save old value
         t_cross = dtime        !update new value
         p_now = t_cross - t_cross_old !compute period

         call opsub3(do1,do2,do3,vx,vy,vz,T_delayed(:,1),T_delayed(:,2),T_delayed(:,3)) !ub=v-vold
         call normvc(h1,semi,l2,linf,do1,do2,do3)
         call opcopy(T_delayed(:,1),T_delayed(:,2),T_delayed(:,3),vx,vy,vz)
         if(nid.eq.0)write(6,*)' Zero-crossing T=',p_now,abs(p_now-p_old),l2
         v_cross_old = v_cross; v_cross = velp(plor)
         if(nid.eq.0)write(17,"(5E15.7)")time, p_now, abs(p_now-p_old), v_mean, l2

      endif
!     https://en.wikipedia.org/wiki/Finite_difference_coefficient
      if(istep.gt.3)then
         vdot  =( (11./6.)*velp(plor)-3*velp(plor-1) +1.5*velp(plor-2) -(1./3.)*velp(plor-3) )*dt**-1
         vddot =(        2*velp(plor)-5*velp(plor-1) +4.0*velp(plor-2)     -1.0*velp(plor-3) )*dt**-2
      else
         vdot=0.0d0; vddot=0.0d0
      endif
      if(nid.eq.0)write(19,"(4E15.7)")time, velp(plor), vdot, vddot

      do i = 1,plor-1
         velp(i)  = velp(i+1)
      enddo

      return
      end subroutine zero_crossing
c-----------------------------------------------------------------------
