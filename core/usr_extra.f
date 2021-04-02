c-----------------------------------------------------------------------
      subroutine nekStab_setDefault
!     use nek_stab
      implicit none
      include 'SIZE'
      include 'TOTAL'

      k_dim = 1
      schur_tgt = 2             ! schur target
      eigen_tol = 1.0e-6        !
      schur_del = 0.10d0        !
      maxmodes = k_dim          ! max modes to outpost

      bst_skp = 1               ! boostconv skip
      bst_snp = 1               ! bootsconv matrix size

      ifres  = .false.          ! outpost restart files for eig
      ifvor  = .false.          ! outpost vorticity
      ifvox  = .false.          ! outpost vortex
      ifldbf = .true.           ! load base flow for stability
      ifbf2D = .false.          ! force 2D solution

!     ifnwt = .false. ! newton-iteration flag
!     ifsfd = .false.
!     ifbcv = .false.
!     iftdf = .false.

      ifseed_nois = .true.      ! noise as initial seed
      ifseed_symm = .false.     ! symmetry initial seed
      ifseed_load = .false.     ! loading initial seed (e.g. Re_ )
!     else all fase -> prescribed by usric

!     arnD_ = .false.
!     arnA_ = .false.
!     arnDA = .false.
!     arnAD = .false.

      xLspg   = 0.0d0; call bcast(xLspg, wdsize)
      xRspg   = 0.0d0; call bcast(xRspg, wdsize)
      yLspg   = 0.0d0; call bcast(yLspg, wdsize)
      yRspg   = 0.0d0; call bcast(yRspg, wdsize)
      zLspg   = 0.0d0; call bcast(zLspg, wdsize)
      zRspg   = 0.0d0; call bcast(zRspg, wdsize)
      acc_spg = 0.0d0; call bcast(acc_spg, wdsize)

!     !Broadcast all defaults !
      call bcast(schur_tgt, isize) ! isize for integer
      call bcast(eigen_tol, wdsize) ! wdsize for real
      call bcast(schur_del, wdsize)
      call bcast(maxmodes, isize)
      call bcast(k_dim, isize)
      call bcast(bst_skp, isize)
      call bcast(bst_snp, isize)

      call bcast(ifres   , lsize) !lsize for boolean
      call bcast(ifvor   , lsize)
      call bcast(ifvox   , lsize)
      call bcast(ifseed_nois  , lsize)
      call bcast(ifseed_symm  , lsize)
      call bcast(ifseed_load  , lsize)
      call bcast(ifldbf  , lsize)
      call bcast(ifbf2D  , lsize)
!     call bcast(ifnewton  , lsize)

      return
      end subroutine nekStab_setDefault
c-----------------------------------------------------------------------
      subroutine nekStab_init
!     initialize arrays and variables defaults
      implicit none
      include 'SIZE'
      include 'TOTAL'
      logical ifto_sav, ifpo_sav
      real glmin,glmax
      integer n,i
      n = nx1*ny1*nz1*nelv

      if (nid==0) then
         print *,'==================================================='
         print *,'== nekStab ========================================'
         print *,'== Copyright (c) 2021 DynFluid Laboratoire ========'
         print *,'==================================================='
         write(6,'(A,A)')' Nek5000 version        : ', NVERSION
         write(6,'(A,A)')' NekStab version        : ', NSVERSION
      endif

      call nekStab_setDefault
      call nekStab_usrchk       ! where user change defaults
      call nekStab_printNEKParams

      xmn = glmin(xm1,n); xmx = glmax(xm1,n)
      ymn = glmin(ym1,n); ymx = glmax(ym1,n)
      zmn = glmin(zm1,n); zmx = glmax(zm1,n)

      if(nid.eq.0)then
         write(6,*)' x min max =',xmn,xmx
         write(6,*)' x   total =',xmx-xmn
         write(6,*)' y min max =',ymn,ymx
         write(6,*)' y   total =',ymx-xmn
         write(6,*)' z min max =',zmn,zmx
         write(6,*)' z   total =',zmx-zmn
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

      ifbfcv = .false.          !

      return
      end subroutine nekStab_init
c-----------------------------------------------------------------------
      subroutine nekStab
!     nekStab main driver
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt) :: qx, qy, qz, qt
      real, dimension(lt2) :: qp

      integer n

      n = nx1*ny1*nz1*nelv

      if(istep.eq.0)call nekStab_init

      call oprzero(fcx,fcy,fcz) ! never comment this!
      call rzero(fct,n)

!     think on a better place for this part!
      if(if3d .AND. ifbf2d)then
         if(nid.eq.0)write(6,*)' Forcing vz = 0'
         call rzero(vz,nx1*ny1*nz1*nelv)
      endif

      select case (floor(uparam(1)))

      case(0)                   ! DNS

         call nekStab_outpost   ! outpost vorticity
         call nekStab_comment   ! print comments
         call nekStab_energy(vx,vy,vz,t(1,1,1,1,1),'global_energy.dat',20)

      case(1)                   ! fixed points computation

         call nekStab_outpost   ! outpost vorticity
         call nekStab_comment   ! print comments

         if(uparam(1).ge.1)then !compose forcings to fcx,fcy,fcz

            if(uparam(3).eq.1)then
               call SFD         !ifSFD=.true.
            elseif(uparam(3).eq.2)then
               call BoostConv   !ifbst=.true.
            elseif(uparam(3).eq.3)then
!     --> Copy initial guess into newton-krylov array.
               call nopcopy(qx, qy, qz, qp, qt, vx, vy, vz, pr, t)
!     --> Newton-Krylov solver.
               call newton_krylov(qx, qy, qz, qp, qt)
!     --> Outpost solution.
               call outpost(qx, qy, qz, qp, qt, "BF_")
               call nek_end
            endif

         endif
         if(ifbfcv)call nek_end

      case(2)                   ! limit cycle computation

         write(6,*) 'NOT IMPLEMENTED'; call nek_end

      case(3)                   ! eigenvalue problem

         if(uparam(3).eq.3)then
            if(nid.eq.0)write(6,*) 'forcing uparam(3) to ZERO!!!!! OTHERWISE CRASH'
            uparam(3)=0; call bcast(uparam, 1*wdsize)

         endif
         call Krylov_Schur
         call nek_end

      case(4)                   ! in postprocessing.f

         if(uparam(01) .eq. 4.1) call wave_maker

         if(uparam(01) .eq. 4.2) call BF_sensitivity

         if( (uparam(01) .eq. 4.31) .or. (uparam(01) .eq. 4.32) ) then
            call ts_steady_force_sensitivity()
         end if

         call nek_end

      end select

      return
      end subroutine nekStab
c-----------------------------------------------------------------------
      subroutine nekStab_outpost
!     nekStab custom outpost routine
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lt=lx1*ly1*lz1*lelt

      real wo1(lt),wo2(lt),wo3(lt),vort(lt,3)
      common /ugrad/ wo1,wo2,wo3,vort

      logical ifto_sav, ifpo_sav

      if((istep.eq.0).OR.ifoutfld.AND.(ifvor.or.ifvox))then

         ifto_sav = ifto; ifpo_sav = ifpo
         ifto = .false.; ifpo = .false.

!---  > Compute and oupost vorticity.
         if(ifvor)then
            call oprzero(wo1,wo2,wo3)
            call oprzero(vort(:,1),vort(:,2),vort(:,3))
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
         write(6,*)'P94=',param(94),'num projection for helmholz solves (controled by ifprojfld)'
         write(6,*)'P95=',param(95),'projection for pressure solver on/off'
         write(6,*)'P101=',param(101),'no additional modes'
         write(6,*)'P103=',param(103),'filter weight'
         write(6,*)'uparam1=',uparam(1)
         write(6,*)'uparam01=',uparam(01)
         write(6,*)'uparam02=',uparam(02)
         write(6,*)'uparam03=',uparam(03)
         write(6,*)'uparam04=',uparam(04)
         write(6,*)'uparam05=',uparam(05)
         write(6,*)'uparam06=',uparam(06)
         write(6,*)'uparam07=',uparam(07)
         write(6,*)'uparam08=',uparam(08)
         write(6,*)'uparam09=',uparam(09)
         write(6,*)'uparam10=',uparam(10)
      endif
      end subroutine nekStab_printNEKParams
c-----------------------------------------------------------------------
      subroutine nekStab_energy(px, py, pz, pt, fname, skip)
!     
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter :: lt=lx1*ly1*lz1*lelt
      real, dimension(lt), intent(in) :: px, py, pz, pt
      integer, intent(in) :: skip
      character(len=20), intent(in) :: fname
      real glsc3,uek,vek,wek,eek,pot
      integer n
      save eek,n
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
         if(nid.eq.0)write(730,"(6E15.7)")time,uek,vek,wek,(uek+vek+wek)*eek,pot*eek
      endif

      return
      end subroutine nekStab_energy
c-----------------------------------------------------------------------
