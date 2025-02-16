!
!     OTD implementation by Simon Kern
!     Based on the method proposed in
!     Babaee,H. and Sapsis,T.P., "A minimization principle for the
!      description of modes associated with finite-time instabilities",
!      Proc. R. Soc. A 474: 20150779,
!      https://dx.doi.org/10.1098/rspa.2015.0779  
!      
!     Email:     skern@mech.kth.se      
!
!     List of subroutines:
!       FRAMEWORK:
!         otd_register
!         otd_init
!       MAIN INTERFACE:
!         run_otd
!       SETUP:
!         white_noise_IC 
!         read_OTDIC
!         OTD_ON
!       LINEARIZED OPERATOR:
!         gen_LU
!         proj_OTDmodes
!         gen_OTD_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)       -> userbc
!       OUTPOSTING
!         outpost_projmodes
!         outpost_OTDbasis
!       FTLEs
!         get_FTLE
!         reset_FTLE
!         compute_FTLE_Phi
!         compute_FTLE_Blanchard
!
!======================================================================

!======================================================================
!> @brief Register OTD module
!! @note This routine should be called in frame_usr_register

      subroutine otd_register()     
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'OTD'

      ! local variables
      integer lpmid
      real ltim

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,otd_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(otd_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'Parent module ['//'FRAME'//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(otd_id,lpmid,otd_name,
     $          'Optimally time-dependent (OTD) basis.')

      ! register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      ! total time
      call mntr_tmr_reg(otd_ttot_id,lpmid,otd_id,
     $     'OTD_TOT','OTD module total time',.false.)
      ! Initialisation time
      call mntr_tmr_reg(otd_tini_id,otd_ttot_id,otd_id,
     $     'OTD_INI','OTD initialisation time',.true.)
      ! Orthonormalisation time
      call mntr_tmr_reg(otd_tGS_id,otd_ttot_id,otd_id,
     $     'OTD_GS','Orthonormalisation time',.true.)
      ! LU & forcing building time
      call mntr_tmr_reg(otd_tgen_id,otd_ttot_id,otd_id,
     $     'OTD_gen','LU and foring computation time',.true.)
      ! OTD mode computation time
      call mntr_tmr_reg(otd_tget_id,otd_ttot_id,otd_id,
     $     'OTD_get','OTD mode computation time',.true.)
      ! FTLE computation time
      call mntr_tmr_reg(otd_tFTLE_id,otd_ttot_id,otd_id,
     $     'OTD_FLTE','FTLE computation time',.true.)
      ! IO
      call mntr_tmr_reg(otd_tIO_id,otd_ttot_id,otd_id,
     $     'OTD_IO','OTD module IO time',.true.)
       
      ! register and set active section
      call rprm_sec_reg(otd_sec_id,otd_id,
     $     '_'//adjustl(otd_name),
     $     'Runtime parameter section for FST module')
      call rprm_sec_set_act(.true.,otd_sec_id)

      ! register parameters
      ! otd_nusrIC
      call rprm_rp_reg(otd_nusric_id,otd_sec_id,'OTD_NUSRIC',
     $     'Number of IC fields to be read ',
     $     rpar_int,0,0.0,.false.,'  ')
      ! otd_IOstep
      call rprm_rp_reg(otd_iostep_id,otd_sec_id,'OTD_IOSTEP',
     $     'IO frequency for the OTD module ',
     $     rpar_int,5,0.0,.false.,'  ')
      ! otd_rststp
      call rprm_rp_reg(otd_rststp_id,otd_sec_id,'OTD_RSTSTP',
     $     'IO frequency for OTD basis (restart) ',
     $     rpar_int,1000,0.0,.false.,'  ')
      ! otd_prntsp
      call rprm_rp_reg(otd_prntsp_id,otd_sec_id,'OTD_PRNTSP',
     $     'Print frequency for the OTD module ',
     $     rpar_int,100,0.0,.false.,'  ')
      ! otd_GSstep
      call rprm_rp_reg(otd_gsstep_id,otd_sec_id,'OTD_GSSTEP',
     $     'Frequency of GS orthonormalization of OTD basis ',
     $     rpar_int,5,0.0,.false.,'  ')
      ! otd_FTLEpd
      call rprm_rp_reg(otd_FTLEpd_id,otd_sec_id,'OTD_FTLEpd',
     $     'Time horizon for FTLE computation ',
     $     rpar_real,0,0.0,.false.,'  ')
      ! otd_ifFTLE
      call rprm_rp_reg(otd_ifFTLE_id,otd_sec_id,'OTD_IFFTLE',
     $     'Compute FTLEs? ',
     $     rpar_log,0,0.0,.false.,'  ')
      ! otd_ifvlog
      call rprm_rp_reg(otd_ifvlog_id,otd_sec_id,'OTD_IFVLOG',
     $     'Verbose logging for OTD module ',
     $     rpar_log,0,0.0,.false.,'  ')
      ! otd_debug
      call rprm_rp_reg(otd_debug_id,otd_sec_id,'OTD_DEBUG',
     $     'Debugging mode for OTD module ',
     $     rpar_log,0,0.0,.false.,'  ')

      ! set initialisation flag
      otd_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_ttot_id,1,ltim)

      return
      end subroutine otd_register

!======================================================================
!> @brief Initilise OTD module
!! @note This routine should be called in frame_usr_init
!! @remark This routine uses global scratch space \a SCRUZ      

      subroutine otd_init()
      implicit none

      include 'SIZE'
      include 'INPUT'               ! param(59), initc
      include 'TSTEP'               ! time
      include 'GEOM'                ! [xyz]m1
      include 'FRAMELP'
      include 'OTD'

      ! local variables
      integer       itmp
      real          rtmp, ltim
      logical       ltmp
      character*20  ctmp
      character*2   str1, str2
      character*200 lstring
      real          xtmp(lx1,ly1,lz1,lelv)
      real          ytmp(lx1,ly1,lz1,lelv)
      real          ztmp(lx1,ly1,lz1,lelv)

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the module was already initialised
      if (otd_ifinit) then
         call mntr_warn(otd_id,
     $        'module ['//trim(otd_name)//'] already initiaised.')
         return
      endif

      ! get runtime parameters
! otd_nusric
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_nusric_id,rpar_int)
      otd_nusric = itmp
! otd_IOstep
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_iostep_id,rpar_int)
      otd_iostep = itmp
! otd_rststp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_rststp_id,rpar_int)
      otd_rststp = itmp
! otd_prntsp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_prntsp_id,rpar_int)
      otd_prntsp = itmp
! otd_GSstep
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_gsstep_id,rpar_int)
      otd_gsstep = itmp
! otd_FTLEpd
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_FTLEpd_id,rpar_real)
      otd_FTLEpd = rtmp
! otd_ifFTLE
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_ifFTLE_id,rpar_log)
      otd_ifFTLE = ltmp
! otd_ifvlog
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_ifvlog_id,rpar_log)
      otd_ifvlog = ltmp
! otd_debug
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_debug_id,rpar_log)
      otd_debug = ltmp

!---------------------------------------------------------
!    Initialise all IC fields to white noise
!
      call white_noise_IC()

!---------------------------------------------------------
!    Read IC fields 
!
      if (otd_nusric.lt.0) then
        call mntr_abort(otd_id, 'Choose a valid number of IC files.')
      else if (otd_nusric.eq.0) then
        call mntr_log(otd_id,lp_inf,
     $ 'No IC files read. ICs will be white noise.')
      else
        if (otd_nusric.gt.npert) then
          call mntr_abort(otd_id,
     $ 'otd_nusric > npert! Increase the number of perturbations.')
        endif
        ! save mesh
        call opcopy(xtmp,ytmp,ztmp,xm1,ym1,zm1)
        ! read ICs
        call read_OTDIC()
        ! restore mesh
        call opcopy(xm1,ym1,zm1,xtmp,ytmp,ztmp)
      endif

!---------------------------------------------------------
!    Sanity check
!
      if (otd_ifFTLE) then
        if (otd_FTLEpd.eq.0.0) then
          call mntr_warn(otd_id,
     $ 'The FLTEs will be computed continuously. Set otd_FTLEpd in the
     $ par-file to define a finite time horizon.')
        endif
        pcount = 0
      endif

!---------------------------------------------------------
!    Ensure right variable size
!
      if (npert.gt.lpert) then
        call mntr_abort(otd_id,
     $ 'npert > lpert! Set lpert=npert for the tool to work.')
      endif

!---------------------------------------------------------
!    Set ICs 
!
      call blank(initc,132)
      call setics
      call mntr_log(otd_id,lp_inf,'INIT - set ICs :: done')
      
!---------------------------------------------------------
!    Set default timestep at which to start OTD computation
!
!    we use the first timestep to make sure that the perturbations
!     1. are divergence free
!     2. satisfy the orthonormality constraint
!     3. satisfy the boundary conditions of the problem
!
      startOTD = 1

      ! everything is initialised
      otd_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tini_id,1,ltim)

      return
      end subroutine otd_init

!=======================================================================
!> @brief Check if module was initialised
!! @return otd_is_initialised

      logical function otd_is_initialised()
      implicit none

      include 'SIZE'
      include 'OTD'
!-----------------------------------------------------------------------
      otd_is_initialised = otd_ifinit

      return
      end function otd_is_initialised

!=======================================================================
!> @brief Main interface for the OTD module
!  
      subroutine run_OTD()
      implicit none

      include 'SIZE'
      include 'TSTEP'           ! istep,iostep
      include 'FRAMELP'
      include 'OTD'

      include 'MASS'

      real op_glsc2_wt
      integer n9,ntot,i
      real iz(lpert)

!----------------------------------------------------------------------

      if (istep.eq.0) then              ! OTD initialisation
!-----------------------------------------------------------------------
!       Perform orthonormalisation of the basis vectors
!       Initialise fundamental solution matrix for FTLE computation
!       Compute one time step 
!
        gsstep_override = .true.        ! initialisation
        call OTD_ON()
!        call outpost_OTDbasis()
        if (otd_ifFTLE) then
          call reset_FTLE
        endif
!
      elseif (istep.ge.startOTD) then   ! standard OTD

!-----------------------------------------------------------------------
!       Orthonormalize perturbations (if necessary)     
!
        gsstep_override = .false.
        call OTD_ON()

!-----------------------------------------------------------------------
!       Build the action of the linearized NS-operator and compute the
!       reduced operator Lr_{ij} 
!
        call gen_Lu()

!-----------------------------------------------------------------------
!       Compute additional forcing for perturbation equation 
!
        call gen_OTD_forcing()
        
!-----------------------------------------------------------------------
!       Compute the OTD modes
!
        if (mod(istep,otd_prntsp).eq.0) then
          call proj_OTDmodes()
        endif

!-----------------------------------------------------------------------
!       Output OTD modes
!
        if (mod(istep,otd_iostep).eq.0) then
          call outpost_projmodes()
        endif

!-----------------------------------------------------------------------
!       Output perturbation fields (OTD basis) for restart
!
        if (mod(istep,otd_rststp).eq.0
     $      .or.istep.eq.nsteps.or.lastep.eq.1) then
          call outpost_OTDbasis()
        endif

!-----------------------------------------------------------------------
!       FTLEs
!
        if (otd_ifFTLE) then
          call get_FTLE()
        endif
!
      endif     ! istep 0 --> init


      return
      end subroutine run_OTD

!=======================================================================
!> @brief Construct action of linearized NS operator on perturbation
!  field
!
      subroutine gen_Lu()
      implicit none
              
      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'FRAMELP'
      include 'MASS'            ! BM1
      include 'SOLN'            ! V[XYZ]P
      include 'TSTEP'           ! istep
      include 'OTD'             ! conv[xyz],diff[xyz],gradp[xyz]
                                ! Lu[xyz], Lr

      ! local variables
      integer ipert, jpert, ntot, i

      ! timing
      real ltim

      ! function
      real op_ip
      real dnekclock

!#define DEBUG
#ifdef DEBUG
      character*1 str
      character*3 oname
#endif

      ! timing
      ltim = dnekclock()

!----------------------------------------------------------------------
!     Build the elements of the linearized NS-operator L_{NS} (u_j)
!
!     L_{NS} (u_j) = 1/Re (grad^2 u)_j - (grad p)_j - (Ub.grad) u_j - (u_j.grad) Ub
!
      do ipert=1,npert

        ! Convective terms 
        call Lu_op_conv(vxp(1,ipert),vyp(1,ipert),vzp(1,ipert),ipert)

        ! Perturbation pressure gradient
        call Lu_op_gradp(prp(1,ipert),ipert)

        ! Diffusive term
        call Lu_op_diff(vxp(1,ipert),vyp(1,ipert),vzp(1,ipert),ipert)
      
      enddo

#ifdef DEBUG
! DIAGNOSTICS of the elements of LU_{NS}
      do i=1,npert
        if (mod(istep,otd_iostep).eq.0) then
          write(str,'(I1)') i
          oname = 'p'//trim(str)//'d'
          call outpost(diffx(1,i),diffy(1,i),diffz(1,i),pr,t,oname)
          oname = 'p'//trim(str)//'p'
          call outpost(gradpx(1,i),gradpy(1,i),gradpz(1,i),pr,t,oname)
          oname = 'p'//trim(str)//'c'
          call outpost(convx(1,i),convy(1,i),convz(1,i),pr,t,oname)
! To compute the individual convective terms, activate in Lu_op_conv
!          oname = 'c'//trim(str)//'1'
!          call outpost(convx1(1,i),convy1(1,i),convz1(1,i),pr,t,oname)
!          oname = 'c'//trim(str)//'2'
!          call outpost(convx2(1,i),convy2(1,i),convz2(1,i),pr,t,oname)
        endif
      enddo
#endif
#undef DEBUG

!----------------------------------------------------------------------
!     Assemble the action of the operator
!
!     L_{NS} (u_j) = 1/Re grad^2 u_j - grad p - (Ub.grad) u_j - (u_j.grad) Ub
!
      ntot = lx1*ly1*lz1*nelv
!      
      do jpert=1,npert
        do i=1,ntot
          Lux(i,jpert) = diffx(i,jpert)- gradpx(i,jpert)- convx(i,jpert)
          Luy(i,jpert) = diffy(i,jpert)- gradpy(i,jpert)- convy(i,jpert)
          if (if3d) then
            Luz(i,jpert)=diffz(i,jpert)- gradpz(i,jpert)- convz(i,jpert)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     Compute the innner product < L_{NS}(u_i),u_j > with i,j = 1,...,r
!
      call rzero(Lr,lpert*lpert)
      do ipert=1,npert
        do jpert=1,npert
          Lr(ipert,jpert) = op_ip(ipert,jpert,2)
        enddo
      enddo
!
! Debug output: Print reduced operator Lr_{ij}
!
      if (nid.eq.0.and.otd_debug) then
        if (mod(istep,otd_prntsp).eq.0) then
          ! write out Lr
          call outmat(Lr,npert,npert,'Lrmat  ',istep)
        endif    ! otd_prntsp
      endif      ! otd_debug.and.nid.eq.0

!----------------------------------------------------------------------
!     Here you could add internal rotations phi_rot into the method.
!     This does NOT change the subspace.
!
!       dU / dt = L_{NS}U - U (Lr - phi_rot)
!
      call rzero(phi_rot,lpert*lpert)
!
!     The rotation matrix Phi_ij must be skew-symmetric (but is otherwise
!     arbitrary
!
!       phi_rot(i,j) = -phi_rot(j,i)
!      
!     e.g. to obtain an evolution that corresponds to continuously
!     performing Gram-Schmidt on the basis (i.e. turning Lr into a lower
!     triangular matrix), set the rotation matrix to
! 
!                  / -<Lu_j,u_i>     j < i
!       phi_rot = {   0              j = i
!                  \  <Lu_j,u_i>     j > i
!
      if (npert.gt.1) then
        do jpert=1,npert
          do ipert=jpert+1,npert
            phi_rot(ipert,jpert) = op_ip(ipert,jpert,2)
            phi_rot(jpert,ipert) = -phi_rot(ipert,jpert)
          enddo
        enddo
      endif
      
      if (nid.eq.0.and.otd_debug) then
        if (mod(istep,otd_prntsp).eq.0) then
          call outmat(Lr,npert,npert,'Lr-mat',istep)
          call outmat(phi_rot,npert,npert,'phimat',istep)
        endif
      endif
      ! add internal rotation if defined
      call sub2(Lr,phi_rot,lpert*lpert)
      if (nid.eq.0.and.otd_debug) then
        if (mod(istep,otd_prntsp).eq.0) then
          call outmat(Lr,npert,npert,'Lrpmat',istep)
        endif
      endif

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tgen_id,1,ltim)

      return
      end subroutine gen_Lu

!=======================================================================
!> @brief Compute eigenspectrum of the reduced operator Lr_{ij} and
!         project the velocity perturbations onto the eigendirections to
!         obtain the most unstable modes
! 
      subroutine proj_OTDmodes()
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'FRAMELP'
      include 'TSTEP'           ! istep
      include 'SOLN'            ! v[xyz]p
      include 'WLAPACK'
      include 'OTD'             ! EIG[RI],EV[RL],Lr,OTDmr[xyz]

      ! local variables
      real    tmp(lpert,lpert)
      integer ntot,i,j
      character*1 str
      character*3 oname
      character*20 fmtr, fmti
      ! timing
      real ltim

      ! function
      real dnekclock
!----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

!-----------------------------------------------------------------------
!     Compute the eigenspectrum of the reduced operator, sort the
!     eigenvalues in decreasing order
!
      ! LAPACK
      ! compute eigenvalues of Lsym = (Lr+Lr^T)/2
      call copy(tmp,Lr,lpert*lpert)                   ! Save Lr
      do i=1,npert
        do j=1,npert
          Lr(i,j) = 0.5*(tmp(i,j)+tmp(j,i))           ! Compute Lsym
        enddo
      enddo
      call eig_wrapper(npert,'r')                     ! Compute lambdas
      call sorteigs('Ls ')
! Output
      if (mod(istep,otd_prntsp).eq.0) then
        if (nid.eq.0) then
          write(fmtr,'("(", I0, "(E15.7,1X))")') npert
          write(6,100,ADVANCE='NO') istep, time, 'Ls | Re'
          write(6,fmtr) (EIGR(i), i=1,npert)
        endif
      endif
      
      ! compute eigenvalues of Lr
      call copy(Lr,tmp,lpert*lpert)                   ! Restore Lr
      call eig_wrapper(npert,'r')                     ! Compute lambdas
      call sorteigs('Lr ')
      call copy(Lr,tmp,lpert*lpert)                   ! Restore Lr

! Output
      if (mod(istep,otd_prntsp).eq.0) then
        if (nid.eq.0) then
          write(6,100,ADVANCE='NO') istep, time, 'Lr | Re '
          write(6,fmtr) (EIGR(i), i=1,npert)
          write(6,100,ADVANCE='NO') istep, time, 'Lr | Im '
          write(6,fmtr) (EIGI(i), i=1,npert)
          ! print order for reference
          write(fmti,'("(", I0, "(I4,1X))")') npert
          write(6,100,ADVANCE='NO') istep, time, 's-idx   '
          write(6,fmti) (idx(i), i=1,npert)
          ! print out non-zero elements of rotated Lr
          write(fmtr,'("(", I0, "(E15.7,1X))")') npert*(npert+1)/2
          write(6,100,ADVANCE='NO') istep, time, 'Lrmat   '
          write(6,fmtr) ( ( Lr(i,j) , j=i,npert ), i=1,npert )
        endif
      endif
  100 format('  [OTD] ',I7,1X,E14.7,1X,A8)

!-----------------------------------------------------------------------
!     Project the perturbation velocity field (OTD basis) onto the 
!     eigendirections of the reduced operator to obtain the most
!     unstable directions
!
      ntot = lx1*ly1*lz1*lelv
      call mxm(vxp,ntot,EVRr,lpert,OTDmrx,npert)
      call mxm(vxp,ntot,EVRi,lpert,OTDmix,npert)
      call mxm(vyp,ntot,EVRr,lpert,OTDmry,npert)
      call mxm(vyp,ntot,EVRi,lpert,OTDmiy,npert)
      if (if3d) then
        call mxm(vzp,ntot,EVRr,lpert,OTDmrz,npert)
        call mxm(vzp,ntot,EVRi,lpert,OTDmiz,npert)
      endif

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tget_id,1,ltim)

      return
      end subroutine proj_OTDmodes

!=======================================================================
!> @brief Output the projection of the OTD basis on the eigendirection
!     as fields
!
      subroutine outpost_projmodes()
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'OTD'

      ! local variables
      integer ipert
      character*2 str
      character*3 oname
      ! timing
      real ltim

      ! function
      real dnekclock
!-----------------------------------------------------------------------
      if (mod(istep,otd_prntsp).ne.0) then
        call proj_OTDmodes()
      endif

      ! timing
      ltim = dnekclock()

      do ipert=1,npert
        if (lpert .ge. 10) then
          write(str,'(I2.2)') ipert
          oname = 'o'//trim(str)
        else
          write(str,'(I1)') ipert
          oname = 'ot'//trim(str)
        endif
        call outpost(OTDmrx(1,ipert),OTDmry(1,ipert),OTDmrz(1,ipert)
     $ ,             prp(1,ipert),t,oname)
      enddo

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tIO_id,1,ltim)

      return
      end subroutine outpost_projmodes

!=======================================================================
!> @brief Output the OTD basis directly to restart.
!     We could alternatively reconstruct the OTD basis from the modes
!     but for this we would need both real and imaginary part. Since we
!     currently only outpost the real part, it's cheaper to just outpost
!     the OTD basis directly when we also outpost the baseflow.     
!
      subroutine outpost_OTDbasis()
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'OTD'

      ! local variables
      integer ipert
      character*2 str
      character*3 oname
      ! timing
      real ltim

      ! function
      real dnekclock
!-----------------------------------------------------------------------
      ! orthonormalize
      gsstep_override = .true.
      call OTD_ON()
      ! timing
      ltim = dnekclock()

      do ipert=1,npert
        write(str,'(I2.2)') ipert
        oname = 'r'//trim(str)
        call outpost(vxp(1,ipert),vyp(1,ipert),vzp(1,ipert)
     $ ,             prp(1,ipert),t,oname)
      enddo

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tIO_id,1,ltim)

      return
      end subroutine outpost_OTDbasis

!=======================================================================
!> @brief Initialise all perturbation fields to random noise
!
      subroutine white_noise_IC()
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'FRAMELP'
      include 'NEKUSE'
      include 'GEOM'            ! [xyz]m1
      include 'PARALLEL'        ! lglel
      include 'OTD'

      ! local variables
      integer ix,iy,iz,ie,ieg,i,ijke
      real    mth_rand,xl(LDIM),fcoeff(3)
!-----------------------------------------------------------------------
      call mntr_log(otd_id,lp_inf,'INIT - creating white noise IC')
      do i=1,npert
        do ix=1,lx1
        do iy=1,ly1
        do iz=1,lz1
          do ie=1,lelv
            xl(1) = XM1(ix,iy,iz,ie)
            xl(2) = YM1(ix,iy,iz,ie)
            if (if3d) then
              xl(NDIM) = ZM1(ix,iy,iz,ie)
            endif
            ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(ie-1)))
            ieg=lglel(ie)
            fcoeff(1)= sin(real(i))**2* 3.0e4
            fcoeff(2)= sin(real(i))**2*(-1.5e3)
            fcoeff(3)= sin(real(i))**2* 0.5e5
            vxpic(ijke,i)=mth_rand(ix,iy,iz,ieg,xl,fcoeff)
            fcoeff(1)= sin(real(i))**2* 2.3e4
            fcoeff(2)= sin(real(i))**2* 2.3e3
            fcoeff(3)= sin(real(i))**2*(-2.0e5)
            vypic(ijke,i)=mth_rand(ix,iy,iz,ieg,xl,fcoeff)
            if (if3d) then
              fcoeff(1)= sin(real(i))**2*2.e4
              fcoeff(2)= sin(real(i))**2*1.e3
              fcoeff(3)= sin(real(i))**2*1.e5
              vzpic(ijke,i)=mth_rand(ix,iy,iz,ieg,xl,fcoeff)
            endif
          enddo
        enddo
        enddo
        enddo
      enddo
      call mntr_log(otd_id,lp_inf,'INIT - white noise :: done')

      return
      end subroutine white_noise_IC

!=======================================================================
!> @brief Read OTD IC fields and run setics 
!
      subroutine read_OTDIC()
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'
      include 'SOLN'
      include 'OTD'
      include 'TSTEP'

      ! local variables
      logical       exist_IC
      integer       i
      character*2   istr
      character*132 ifile
!-----------------------------------------------------------------------

      call mntr_log(otd_id,lp_inf,'INIT - read ICs')
      do i=1,otd_nusric
        write(istr,'(I0.2)') i
        ifile='OTDIC_'//trim(istr)//'.fld'
        inquire (file=ifile,exist=exist_IC)
        if (exist_IC) then
          call mntr_log(OTD_id,lp_inf,
     $             'INIT - Reading IC file '//trim(ifile))
          call load_fld(ifile)
          ! Copy the initial conditions into the fields v[xyz]pic
          call opcopy(vxpic(1,i),vypic(1,i),vzpic(1,i),vx,vy,vz)
          OTDrsttime = time
          call mntr_logr(otd_id,lp_inf,
     $                   'IC file '//trim(ifile)//' time:',time)
        else
          call mntr_abort(otd_id,'Cannot open '//trim(ifile)//' !')
        endif
      enddo
      call mntr_log(otd_id,lp_inf,'INIT - read ICs :: done')

      return
      end subroutine read_OTDIC

!=======================================================================
!> @brief Create forcing for OTD evolution equation 
!
      subroutine gen_OTD_forcing()
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'OTD'
      include 'INPUT'   ! if3d
      include 'TSTEP'

      ! local variable list
      integer ntot

      ! timing
      real ltim

      ! function 
      real dnekclock

!#define DEBUG
#ifdef DEBUG
      character*1 str
      character*3 oname
      integer i
      real glmax 
      real gl(3)
      real dudtx (lx1*ly1*lz1*lelv,lpert)
     $ ,   dudty (lx1*ly1*lz1*lelv,lpert)
     $ ,   dudtz (lx1*ly1*lz1*lelv,lpert)
#endif
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ntot = lx1*ly1*lz1*lelv

      call mxm(VXP,ntot,Lr,lpert,OTDfx,npert)
      call mxm(VYP,ntot,Lr,lpert,OTDfy,npert)
      call mxm(VZP,ntot,Lr,lpert,OTDfz,npert)

#ifdef DEBUG
      if (mod(istep+1,otd_iostep).eq.0) then
        do i=1,npert
          write(str,'(I1)') i
          oname = 'p'//trim(str)//'f'
          call outpost(OTDfx(1,i),OTDfy(1,i),OTDfz(1,i),pr,t,oname)
          call opcopy(dudtx(1,i),dudty(1,i),dudtz(1,i)
     $ ,              Lux(1,i),Luy(1,i),Luz(1,i))
          call sub2(dudtx(1,i),OTDfx(1,i),ntot)
          call sub2(dudty(1,i),OTDfy(1,i),ntot)
          call sub2(dudtz(1,i),OTDfz(1,i),ntot)
          oname = 'p'//trim(str)//'u'
          call outpost(dudtx(1,i),dudty(1,i),dudtz(1,i),prp,t,oname)
          oname = 'p'//trim(str)//'l'
          call outpost(Lux(1,i),Luy(1,i),Luz(1,i),prp,t,oname)
        enddo
      endif
#endif
#undef DEBUG
      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tgen_id,0,ltim)

      return 
      end subroutine gen_OTD_forcing

!=======================================================================
!> @brief Set forcing for OTD evolution equation 
!     This function is called in userf for each GLL point. Therefore we
!     need a manual switch for when to start setting the forcing.
!
      subroutine set_OTD_forcing(FFX,FFY,FFZ,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'
      include 'SOLN'           ! jp
      include 'TSTEP'          ! istep
      include 'PARALLEL'       ! gllel
      include 'OTD'

      ! argument list
      integer ix,iy,iz,ieg
     
      ! output 
      real    ffx,ffy,ffz

      ! local variable list
      integer ijke,e
!-----------------------------------------------------------------------

      e = gllel(ieg)

      ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(e-1)))
      if (jp.ne.0) then
!       only for the perturbations
        FFX = FFX - OTDfx(ijke,jp)
        FFY = FFY - OTDfy(ijke,jp)
        FFZ = FFZ - OTDfz(ijke,jp)
      endif

      return 
      end subroutine set_OTD_forcing
 
!=======================================================================
!> @brief Orthonormalize OTD basis 
!
      subroutine OTD_ON()
      implicit none

      include 'SIZE'
      include 'OTD'
     
      include 'SOLN'    ! v[xyz]p 
      include 'TSTEP'   ! istep

      ! local variables
      real N,O
      logical runON

      ! timing
      real ltim

      ! function
      real dnekclock
!----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()
!
      call compute_NO(N,O,'pre   ',.true.)
!
      runON = .false.
      if (otd_gsstep.ne.0) then                 ! gsstep=0 => no GS
        if (mod(istep,otd_gsstep).eq.0) then
          runON = .true.
        endif
      endif
!     override
      if (gsstep_override) runON = .true.       ! override when needed

!     runON
      if (runON) then
!        call CGS()        ! Classical Gram-Schmidt
        call MGS()        ! Modified Gram-Schmidt
      endif
!
      call compute_NO(N,O,'post  ',.false.)
!
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tGS_id,1,ltim)
 
      return     
      end subroutine OTD_ON

!=======================================================================
!> @brief Reset FTLE computation 
!
      subroutine reset_FTLE()
      implicit none

      include 'SIZE'
      include 'OTD'
      
      call rzero(FTLEv,lpert)
      ! Phi
      !call ident(Phi,lpert)
      ! Blanchard
      call rzero(trapz,lpert)

      end subroutine reset_FTLE

!=======================================================================
!> @brief Compute FTLEs (Based on Babaee et al., 2017) 
!
!     1. Advect fundamental solution matrix Phi
!
!       dPhi/dt = Lr    , with Phi(t=t0) = I(rxr)
!
!     2. Compute FTLEs
!
!       FTLE(i) = 1/T * log(svd(Phi(t)))
!
      subroutine compute_FTLE_Phi(Lrp,Lrc,deltat,ftledt)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'

      ! argument list
      real deltat                       ! time interval for FTLE computation 
      real ftledt                       ! dt for Phi advection
      real Lrp   (lpert,lpert)          ! Lr matrix from previous step
      real Lrc   (lpert,lpert)          ! Lr matrix for current step

      ! local variables
      real lhs   (lpert,lpert)          ! lhs of Phi advection equation
     $ ,   invlhs(lpert,lpert)          ! inverse of lhs
     $ ,   tmp   (lpert,lpert)           
     $ ,   rhs   (lpert,lpert)          ! rhs of Phi advection equation
      real fact
      integer i

      ! Advect the fundamental solution matrix (using the implicit CN scheme)
      fact = 0.5*ftledt
      ! build LHS
      call ident(lhs,lpert)
      call add2s2(lhs,Lrc,-fact,lpert*lpert)
      ! build RHS
      call ident(tmp,npert)
      call add2s2(tmp,Lrp, fact,lpert*lpert)
      call mxm(tmp,lpert,Phi,lpert,rhs,lpert)
      ! invert LHS and solve system
      call invmt(lhs,invlhs,tmp,npert)
      call mxm(invlhs,lpert,rhs,lpert,Phi,lpert)

      ! compute SVD
      call copy(VMATX,Phi,lpert*lpert)
      call svd_wrapper(lpert,lpert,'A')

      ! compute FTLEs
      do i=1,npert
        FTLEv(i) = log(OSIGMA(i))/deltat
      enddo

      end subroutine compute_FTLE_Phi

!=======================================================================
!> @brief Compute FTLEs (Based on Blanchard & Sapsis, 2019) 
!
!       FTLE(i) = 1/T * ( int_(t_0)^t <Lu_i,u_i> d tau )
!
      subroutine compute_FTLE_Blanchard(Lrp,Lrc,deltat,ftledt)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'

      ! argument list
      real deltat                       ! time interval for FTLE computation 
      real ftledt                       ! dt for Phi advection
      real Lrp  (lpert,lpert)           ! Lr matrix from previous step
      real Lrc  (lpert,lpert)           ! Lr matrix for current step

      ! local variables
      integer i

      ! compute FTLEs
      do i=1,npert
        ! Trapezoid rule, this could be done better using the AB scheme
        trapz(i) = trapz(i) + 0.5*ftledt*(Lrp(i,i) + Lrc(i,i))
        FTLEv(i) = trapz(i)/deltat
      enddo

      end subroutine compute_FTLE_Blanchard

!=======================================================================
!> @brief Get FTLEs
!      
      subroutine get_FTLE()

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'             ! Lr, FTLEv

      ! local variables
      real pfrac                ! current fraction of the FTLE comp. period
      real ftledt               ! dt for FTLE computation
      real Lrp(lpert,lpert)     ! Lr from previous step
      real Lrc(lpert,lpert)     ! (linear approx.) of Lr at end of int. interval
      real fact
      real t0
      integer icalld
      ! timing
      real ltim
      ! function
      real dnekclock
      ! save variables
      save Lrp
      data t0 /0.0/
      save t0
      data icalld /0/
      save icalld

!----------------------------------------------------------------------

      ! timing
      ltim = dnekclock()

      ! determine FTLE horizon
      if (otd_FTLEpd.eq.0.0) then
        if (icalld.eq.0) then
          t0 = time
          icalld = 1
        endif
        period = time-t0
        pfrac = time-t0
      else
        period = otd_FTLEpd
        pfrac = mod(time,period)
      endif
      ! initialisation
      if (istep.eq.startOTD) then
        call copy(Lrp,Lr,lpert*lpert)
      endif

      if (pfrac.lt.dt) then
        ftledt = dt-pfrac
        ! compute approximation of Lrc at end of period (linear interp.)
        call copy(Lrc,Lrp,lpert*lpert)
        fact = ftledt/dt
        call add2s2(Lrc,Lrp,-fact,lpert*lpert)
        call add2s2(Lrc,Lr , fact,lpert*lpert)
        !call compute_FTLE_Phi(Lrp,Lrc,period,ftledt)
        call compute_FTLE_Blanchard(Lrp,Lrc,period,ftledt)
        pcount = pcount + 1
        if (nid.eq.0) then
          write(6,100,ADVANCE='NO') pcount,istep,time
          do i=1,npert
            write(6,102,ADVANCE='NO') FTLEv(i)
          enddo
          write(6,*)
        endif
        call reset_FTLE
        call copy(Lrp,Lrc,lpert*lpert)
        call copy(Lrc,Lr ,lpert*lpert)
        !call compute_FTLE_Phi(Lrp,Lrc,period,pfrac)
        call compute_FTLE_Blanchard(Lrp,Lrc,period,pfrac)
      else
        call copy(Lrc,Lr,lpert*lpert)
        !call compute_FTLE_Phi(Lrp,Lrc,pfrac,dt)
        call compute_FTLE_Blanchard(Lrp,Lrc,pfrac,dt)
      endif
      ! update Lrp
      call copy(Lrp,Lr,lpert*lpert)
 
      if (nid.eq.0.and.mod(istep,otd_prntsp).eq.0) then
        write(6,101,ADVANCE='NO') istep,time,pfrac
        do i=1,npert
          write(6,102,ADVANCE='NO') FTLEv(i)
        enddo
        write(6,*)
      endif
  100 format(' [OTD] FTLE PRD',1X,I5,1X,I7,' t=',1X,E15.7)
  101 format(' [OTD] FTLE (t)',1X,I7,' t=',2(1X,E15.7))
  102 format(1X,E15.7)

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tFTLE_id,1,ltim)

      end subroutine get_FTLE


!
!     Tools used in the computation of the OTD modes
!     Based on the method proposed in
!     Babaee,H. and Sapsis,T.P., "A minimization principle for the
!      descripertion of modes associated with finite-time instabilities",
!      Proc. R. Soc. A 474: 20150779,
!      https://dx.doi.org/10.1098/rspa.2015.0779  
!      
!     Simon Kern
!     Email:     skern@mech.kth.se 
!
!     List of subroutines:
!       Lu_op_conv(uxp,uyp,uzp,ipert)
!                  IN:  uxp,uyp,uzp         perturbation velocity fields   
!                  IN:  ipert               index
!                  OUT: conv[xyz] (OTD)     convective part of Lu
!       Lu_op_gradp(prpert,ipert)
!                  IN:  prpert              perturbation pressure field
!                  IN:  ipert               index
!                  OUT: gradp[xyz] (OTD)    pressure gradient part of Lu
!       Lu_op_diff(uxp,uyp,uzp,ipert)    
!                  IN:  uxp,uyp,uzp         perturbation velocity fields   
!                  IN:  ipert               index
!                  OUT: diff[xyz] (OTD)     diffusive part of Lu
!       laplacian(lapu,up)
!                  IN:  up                  perturbation velocity field
!                  OUT: lapu                laplacian of up
!       op_ip_vp(uxp,uyp,uzp,jpert)
!                  IN:  uxp,uyp,uzp         input velocity fields   
!                  IN:  jpert               index
!       op_ip(ipert,jpert,iflag)
!                  IN:  ipert               index of first pert. field   
!                  IN:  jpert               index of second pert. field
!                  IN:  iflag               flag to switch u_i / Lu_i
!       op_norm(uxp,uyp,uzp)
!                  IN:  uxp,uyp,uzp         input velocity fields   
!       CGS (classical Gram-Schmidt orthonormalisation)  
!                  IN:  ---
!                  OUT: ---
!       MGS (modified Gram-Schmidt orthonormalisation)
!                  IN:  ---
!                  OUT: ---
!       computeNO(N,O,info,flag) (compute measures for normality and orthogonality)
!                  IN:  info (6 chars)      info
!                  IN:  flag (logical)      print N,L
!                  OUT: N,O
!       sorteigs (sort eigenvalues in decreasing order)
!                  IN:  ---
!                  OUT: ---
!     
!======================================================================

!======================================================================
!> @brief Construct the convective terms for Lu
!
!     Lu_conv = (u.grad) Ub + (Ub.grad) u
!
!     Using the convop routine takes care of the dealiasing.
!
      subroutine Lu_op_conv(uxp,uyp,uzp,ipert)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'SOLN'            ! v[xyz]
      include 'OTD'             ! conv[xyz]

      ! argument list
      real uxp (lx1*ly1*lz1*lelv)   ! perturbation velocity components
     $ ,   uyp (lx1*ly1*lz1*lelv)
     $ ,   uzp (lx1*ly1*lz1*lelv)
      integer ipert                 ! index of the considered perturbation

      ! local variables
      integer ntot,i
      real TA1 (LX1,LY1,LZ1,LELV)
     $ ,   TA2 (LX1,LY1,LZ1,LELV)
     $ ,   TA3 (LX1,LY1,LZ1,LELV)
     $ ,   TB1 (LX1,LY1,LZ1,LELV)
     $ ,   TB2 (LX1,LY1,LZ1,LELV)
     $ ,   TB3 (LX1,LY1,LZ1,LELV)
!-----------------------------------------------------------------------

      ntot = lx1*ly1*lz1*lelv
!
      if (if3d) then
        call opcopy  (tb1,tb2,tb3,vx,vy,vz)         ! Save velocity
        call opcopy  (vx,vy,vz,uxp,uyp,uzp)         ! U <-- u
! convop(conv,fld): builds the convective term for the scalar field fld
!       conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
        call convop  (ta1,tb1)                      ! (u.grad) Ub
        call convop  (ta2,tb2)                                
        call convop  (ta3,tb3)
        ! Copy fields into the correct variables
        call opcopy  (convx(1,ipert),convy(1,ipert),convz(1,ipert)
     $ ,              ta1,ta2,ta3)
        call opcopy  (vx,vy,vz,tb1,tb2,tb3)         ! Restore velocity
c
!       conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
        call convop  (tb1,uxp)                      ! (Ub.grad) u
        call convop  (tb2,uyp)
        call convop  (tb3,uzp)
        ! Add fields to the convective term
        call opadd2  (convx(1,ipert),convy(1,ipert),convz(1,ipert)
     $ ,              tb1,tb2,tb3)

!   DIAGNOSTICS      
!        call opcopy  (convx1(1,ipert),convy1(1,ipert),convz1(1,ipert)
!     $ ,              ta1,ta2,ta3)
!        call opcopy  (convx2(1,ipert),convy2(1,ipert),convz2(1,ipert)
!     $ ,              tb1,tb2,tb3)
      else ! 2D
        call opcopy  (tb1,tb2,tb3,vx,vy,vz)         ! Save velocity
        call opcopy  (vx,vy,vz,uxp,uyp,uzp)         ! U <-- u
! convop(conv,fld): builds the convective term for the scalar field fld
!       conv_i = (v_j.grad_j)*fld_i                 => (vp_j.grad_j)*v_i
        call convop  (ta1,tb1)                      ! (u.grad) Ub
        call convop  (ta2,tb2)                                
        ! Copy fields into the correct variables
        call opcopy  (convx(1,ipert),convy(1,ipert),ta3,ta1,ta2,ta3)
        call opcopy  (vx,vy,vz,tb1,tb2,tb3)         ! Restore velocity
c
!       conv_i = (v_j.grad_j)*fld_i                 => (v_j.grad_j)*vp_i
        call convop  (tb1,uxp)                      ! (Ub.grad) u
        call convop  (tb2,uyp)
        ! Add fields to the convective term
        call opadd2  (convx(1,ipert),convy(1,ipert),tb3,tb1,tb2,tb3)
!   DIAGNOSTICS      
!        call opcopy  (convx1(1,ipert),convy1(1,ipert),ta3,ta1,ta2,ta3)
!        call opcopy  (convx2(1,ipert),convy2(1,ipert),tb1,tb1,tb2,tb3)
      endif ! if3d

      return
      end subroutine Lu_op_conv

!======================================================================
!> @brief Construct the pressure gradient term for Lu
!
!     Lu_gradp = grad p
!
!     Note: The pressure gradient is computed directly on the v-mesh!
!
      subroutine Lu_op_gradp(prpert,ipert)
      implicit none

      include 'SIZE'
      include 'OTD'             ! gradp[xyz]

      ! argument list
      real prpert (lx2*ly2*lz2*lelv,1) ! perturbation pressure field
      integer ipert                    ! number of the considered pert.

      ! local variables
      integer ntot
      real ta1 (lx1,ly1,lz1,lelv)  
     $ ,   ta2 (lx1,ly1,lz1,lelv)  
     $ ,   wrk (lx1,ly1,lz1,lelv)  
!-----------------------------------------------------------------------

      ntot = lx2*ly2*lz2*lelv
!     Map the perturbation pressure to the velocity mesh
      call mappr(wrk,prpert,ta1,ta2)
!     compute the gradient on the velocity mesh directly
      call gradm1(gradpx(1,ipert),gradpy(1,ipert),gradpz(1,ipert),wrk)

      return
      end subroutine Lu_op_gradp

!======================================================================
!> @brief Construct the diffusive term for Lu
!
!     Lu_op_diff = 1/Re*grad^2 u
!
      subroutine Lu_op_diff(uxp,uyp,uzp,ipert)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'SOLN'            ! vdiff
      include 'OTD'             ! diff[xyz]

      ! argument list
      real uxp (lx1*ly1*lz1*lelv,1) ! perturbation velocity components
     $ ,   uyp (lx1*ly1*lz1*lelv,1)
     $ ,   uzp (lx1*ly1*lz1*lelv,1)
      integer ipert                 ! number of the considered pert.
      ! local variables
      integer ntot
!-----------------------------------------------------------------------

      ntot = lx1*ly1*lz1*lelv
      ! compute laplacian
      call laplacian(diffx(1,ipert),uxp)
      call laplacian(diffy(1,ipert),uyp)
      if (if3d) call laplacian(diffz(1,ipert),uzp) 
      ! multiply by 1/Re                > remove for operator diagnostics
      call col2(diffx(1,ipert),vdiff,ntot) 
      call col2(diffy(1,ipert),vdiff,ntot) 
      if (if3d) call col2(diffz(1,ipert),vdiff,ntot) 

      return
      end subroutine Lu_op_diff

!======================================================================
!> @brief Construct the diffusion term (laplacian of u) for direction i
!
      subroutine laplacian(lapu,up)
      implicit none
!
      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'DXYZ'            ! dxm1,d[xy]tm1
      include 'GEOM'            ! r[xy]m1,s[xy]m1,t[xy]m1,jacmi
!
      ! argument list
      real up (lx1*ly1*lz1*lelv,1)       ! perturbation velocity component
!
      ! output
      real lapu (lx1*ly1*lz1,lelv)
!
      ! local variables
      real ux  (lx1*ly1*lz1,lelv)
     $ ,   uy  (lx1*ly1*lz1,lelv)
     $ ,   uz  (lx1*ly1*lz1,lelv)
     $ ,   ur  (lx1*ly1*lz1)
     $ ,   us  (lx1*ly1*lz1)
     $ ,   ut  (lx1*ly1*lz1)
!
      common /ctmp1/ ur,us,ut
!
      integer e,i,lxyz,nel,N
!-----------------------------------------------------------------------

      lxyz = lx1*ly1*lz1
      nel = nx1-1
      call gradm1(ux,uy,uz,up)
      do e=1,lelt
        if (if3d) then
          call local_grad3(ur,us,ut,ux,nel,e,dxm1,dxtm1)
          do i=1,lxyz
            lapu(i,e) = jacmi(i,e)*(  ur(i)*rxm1(i,1,1,e)
     $                              + us(i)*sxm1(i,1,1,e)
     $                              + ut(i)*txm1(i,1,1,e) )
          enddo
          call local_grad3(ur,us,ut,uy,nel,e,dxm1,dxtm1)
          do i=1,lxyz
            lapu(i,e) = lapu(i,e) + jacmi(i,e)*(  ur(i)*rym1(i,1,1,e)
     $                                          + us(i)*sym1(i,1,1,e)
     $                                          + ut(i)*tym1(i,1,1,e) )
          enddo
          call local_grad3(ur,us,ut,uz,nel,e,dxm1,dxtm1)
          do i=1,lxyz   
            lapu(i,e) = lapu(i,e) + jacmi(i,e)*(  ur(i)*rzm1(i,1,1,e)
     $                                          + us(i)*szm1(i,1,1,e)
     $                                          + ut(i)*tzm1(i,1,1,e) )
          enddo
        else ! 2D
          call local_grad2(ur,us,ux,nel,e,dxm1,dytm1)
          do i=1,lxyz
            lapu(i,e) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                            + us(i)*sxm1(i,1,1,e) )
          enddo
          call local_grad2(ur,us,uy,nel,e,dxm1,dytm1)
          do i=1,lxyz
            lapu(i,e) = lapu(i,e)
     $                  + jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                              + us(i)*sym1(i,1,1,e) )
          enddo
        endif ! if3d
      enddo
!      
      return
      end subroutine laplacian
!
!======================================================================
!> @brief Compute the global inner product with the pert. velocities 
!
!     < vc_i, v[xyz]p_j >
!      
      real function op_ip_vp(vcx,vcy,vcz,jpert)
      implicit none

      include 'SIZE'
      include 'SOLN'            ! v[xyz]p, jp
      include 'TSTEP'           ! ifield
      include 'FRAMELP'         ! mntr_abort
      include 'MASS'            ! bm1

      ! argument list
      real vcx(lx1*ly1*lz1*lelv)
      real vcy(lx1*ly1*lz1*lelv)
      real vcz(lx1*ly1*lz1*lelv)
      integer jpert
      ! functions and local variables
      real op_glsc2_wt

      ifield = 1
      op_ip_vp = 0.5*op_glsc2_wt(VCX,         VCY,         VCZ,
     $                           VXP(1,jpert),VYP(1,jpert),VZP(1,jpert),
     $                           bm1)

      return
      end function op_ip_vp
!
!======================================================================
!> @brief Compute the global inner product for the perturbations 
!
!     < v[xyz]p_i , v[xyz]p_j >  or  < LU[xyz]_i , v[xyz]p_j >
!      
      real function op_ip(ipert,jpert,iflag)
      implicit none

      include 'SIZE'
      include 'SOLN'            ! v[xyz]p, jp
      include 'FRAMELP'         ! mntr_abort
      include 'OTD'             ! LU[xyz], otd_id

      ! argument list
      integer ipert,jpert       ! [IN]  perturbation indices
      integer iflag             ! [IN]  iflag = 1: < U_i, U_j >
      !                                 iflag = 2: < LU_i, U_j >
      ! functions and local variables
      real op_ip_vp

      if (iflag.eq.1) then
        op_ip = op_ip_vp(VXP(1,ipert),VYP(1,ipert),VZP(1,ipert),jpert)
      elseif (iflag.eq.2) then
        op_ip = op_ip_vp(LUx(1,ipert),LUy(1,ipert),LUz(1,ipert),jpert)
      else
        call mntr_abort(otd_id,'Error in op_ip!')
      endif

      return
      end function op_ip
!
!======================================================================
!> @brief Normalize vector field 
!
!     u_i = v_i/||v_i||
!      
      subroutine op_norm(uxp,uyp,uzp)
      implicit none

      include 'SIZE'
      include 'TSTEP'           ! ifield
      include 'MASS'            ! bm1
      include 'FRAMELP'         ! mntr_abort
      include 'INPUT'           ! if3d
      include 'OTD'

      ! argument list
      real uxp (lx1*ly1*lz1*lelv,1) ! perturbation velocity components
     $ ,   uyp (lx1*ly1*lz1*lelv,1)
     $ ,   uzp (lx1*ly1*lz1*lelv,1)

      ! functions and local variables
      real invnorm,n2,ntot
      real op_glsc2_wt

      ifield = 1
      ntot = lx1*ly1*lz1*lelv
      n2 = 0.5*op_glsc2_wt(uxp,uyp,uzp,uxp,uyp,uzp,bm1)
      if (n2.le.0.0) then
        call mntr_abort(otd_id,'Error in op_norm!')
      endif
      invnorm = 1./sqrt(n2)
      call cmult(uxp,invnorm,ntot)
      call cmult(uyp,invnorm,ntot)
      if (if3d) call cmult(uzp,invnorm,ntot)

      return
      end subroutine op_norm
!
!======================================================================
!> @brief Perform Classical Gram-Schmidt (CGS) orthonormalization on the 
!         perturbation velocity field 
!
!     u_k = v_k - sum_{j = 1 -> k-1} proj_{u_j} (u_k)
!
!        with proj_{u_j} (u_k) = < u_k , u_j >/||u_j|| * u_j
!
      subroutine CGS()
      implicit none
!
      include 'SIZE'
      include 'INPUT'   ! if3d
      include 'TSTEP'   ! istep
      include 'SOLN'    ! V[XYZ]P
      include 'FRAMELP' ! lp_inf
      include 'OTD'     ! otd_id
!
      ! local variables
      integer i,j,ntot
      real invnorm, proj

      ! functions
      real op_ip
!-----------------------------------------------------------------------
      ntot = lx1*ly1*lz1*lelv
      ! orthonormalize
      do i=1,npert
        do j=1,i-1
          proj = -op_ip(i,j,1)/op_ip(j,j,1)
          call add2s2(vxp(1,i),vxp(1,j),proj,ntot)
          call add2s2(vyp(1,i),vyp(1,j),proj,ntot)
          if (if3d) then
            call add2s2(vzp(1,i),vzp(1,j),proj,ntot)
          endif
        enddo
        invnorm = 1/sqrt(op_ip(i,i,1))
        call cmult(vxp(1,i),invnorm,ntot)
        call cmult(vyp(1,i),invnorm,ntot)
        if (if3d) then
          call cmult(vzp(1,i),invnorm,ntot)
        endif
      enddo
! 
!     stamp logs
      call mntr_logi(otd_id,lp_inf,
     $               'V[XYZ]P orthonormalized (CGS)',istep)

      return 
      end subroutine CGS
!
!======================================================================
!> @brief Perform Modified Gram-Schmidt (MGS) orthonormalization on the 
!         perturbation velocity field for improved numerical stability 
!
!     do i=1,npert
!       u_i = v_i/||v_i||
!       do j=i+1,npert
!         u_j = v_j - proj_{u_i} (v_j)
!       enddo
!     enddo
!
!        with proj_{u_i} (v_j) = < v_j , u_i >/||u_i|| * v_j
!                              = < v_j , u_i > * v_j   since ||u_i|| = 1
!
      subroutine MGS()
      implicit none
!
      include 'SIZE'
      include 'INPUT'   ! if3d
      include 'TSTEP'   ! istep
      include 'SOLN'    ! V[XYZ]P
      include 'FRAMELP' ! lp_inf
      include 'OTD'     ! otd_id
!
      ! local variables
      integer i,j,ntot
      real invnorm, proj

      ! functions
      real op_ip
!-----------------------------------------------------------------------
      ntot = lx1*ly1*lz1*lelv
      ! orthonormalize
      do i=1,npert
        invnorm = 1/sqrt(op_ip(i,i,1))
        call cmult(vxp(1,i),invnorm,ntot)
        call cmult(vyp(1,i),invnorm,ntot)
        if (if3d) then
          call cmult(vzp(1,i),invnorm,ntot)
        endif
        do j=i+1,npert
          proj = -op_ip(i,j,1)
          call add2s2(vxp(1,j),vxp(1,i),proj,ntot)
          call add2s2(vyp(1,j),vyp(1,i),proj,ntot)
          if (if3d) then
            call add2s2(vzp(1,j),vzp(1,i),proj,ntot)
          endif
        enddo
      enddo
! 
!     stamp logs
      call mntr_logi(otd_id,lp_inf,
     $               'V[XYZ]P orthonormalized (MGS)',istep)

      return 
      end subroutine MGS
!
!======================================================================
!> @brief Compute measures of orthogonality and normality of the
!         perturbations 
!
      subroutine compute_NO(N,O,info,flag)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'

      ! argument list
      real        N,O
      logical     flag
      character*6 info

      ! local variables
      integer     i,j,n9
      real        ip(npert,npert)

      ! function
      real        op_ip

      do i=1,npert
        do j=1,npert
          ip(i,j) = op_ip(i,j,1)
        enddo
      enddo

      if (nid.eq.0) then
        ! compute measure for basis vector normality
        N = 0
        do i=1,npert
          N = N + ip(i,i)**2
        enddo
        N = sqrt(N/npert)

        ! compute measure for basis vector orthogonality
        if (npert.gt.1) then
         O = 0
          do i=1,npert
            do j=i+1,npert
              O = O + ip(i,j)**2
            enddo
          enddo
          O = sqrt(2*O)/(npert*(npert-1))
        endif

        ! output
        if (otd_debug) then            ! debug output pre MGS
          n9 = min(npert,9)
          write(6,*) ' [OTD] ipout   istp ipert
     $                <u_i,u_j> ON '//trim(info)
          do i=1,npert
            write(6,100) istep,i,(ip(i,j),j=1,n9)
          enddo
        endif
        if (flag) write(6,101) istep,N-1.0,O
      endif
 100  format('  [OTD] ipout',I7,1X,I5,2X,9(1X,F9.6))
 101  format('  [OTD] NOout',I7,1X,'N-1.0',1X,E15.7,1X,'O',1X,E15.7)
    
      return
      end subroutine compute_NO

!======================================================================
!> @brief 1. Sort the eigenvalues l_i such that their real parts are 
!            ranked in decreasing order
!
!            Re(l_1) .ge. Re(l_i) .ge. Re(l_r), i = 1,...,r
!
!         2. Apply the same sorting to the columns of the right 
!            eigenvector matrix and separate real and imaginary parts
!
!            EVR => EVRR + i*EVRI
! 
      subroutine sorteigs(str)
      implicit none
!
      include 'SIZE'
      include 'OTD'     ! EIG[IR], EVR, EVR[RI]
      include 'TSTEP'

      ! argument list
      character*3 str
!
      ! local variables
      integer i,j,id
      logical mk(lpert)
      real    wrk1(lpert)
      real    wrk2(lpert)
!-----------------------------------------------------------------------
      ! zero out indices, output and mask
      call izero(idx,lpert)
      call rzero(EVRR,lpert*lpert)
      call rzero(EVRI,lpert*lpert)
      do i=1,npert
        mk(i) = .true.
      enddo
      call copy(wrk1,EIGR,lpert)
      call copy(wrk2,EIGI,lpert)

      ! we need to exclude the trailing zeros for the sorting to work
      if (lpert.gt.npert) then
        do i=npert+1,lpert
          mk(i) = .false.
        enddo
      endif
      if (otd_debug) then
        if (nid.eq.0) then
          call print_eigenvalues(' [OTD debug] e-vals: '//trim(str)
     $ ,                         npert,EIGR,EIGI)
        call print_eigenvectors(' [OTD debug] right e-vecs: '//trim(str)
     $ ,                         npert,EIGI,EVR,npert)
        endif
      endif
      ! sorting
      j = 1
      do while (j.le.npert)
        EIGR(j)    = maxval(wrk1,mask=mk)       ! find largest real eigenvalue in remaining list
        id         = maxloc(wrk1,1,mk)          ! find its index
        EIGI(j)    = wrk2(id)                   ! extract corresponding imaginary part
        mk(id)     = .false.                    ! update mask
        idx(j)     = id
        if (EIGI(j).eq.0.0) then
          do i=1,npert
            EVRR(i,j) = EVR(i,id)
            ! EVRI(i,j) = 0.0
          enddo
          j=j+1
        else  ! complex conjugate eigenvectors!
          EIGI(j+1)    = -EIGI(j)
          EIGR(j+1)    =  EIGR(j)
          mk(id+1)     = .false.
          idx(j+1)     = id+1
          do i=1,npert
            EVRR(i,j)   =  EVR(i,id)
            EVRR(i,j+1) =  EVR(i,id)
            EVRI(i,j)   =  EVR(i,id+1)
            EVRI(i,j+1) = -EVR(i,id+1)
          enddo
          j=j+2
        endif
      enddo

! debugging output
      if (otd_debug) then
        if (nid.eq.0) then
        call print_eigenvalues(' [OTD debug] sorted e-vals: '//trim(str)
     $ ,                         npert,EIGR,EIGI)
          write(6,*)
          write(6,*) ' [OTD debug] sorted right e-vecs: '//trim(str)
          do i=1,npert
            do j=1,npert
              write(6,200,ADVANCE='NO') EVRR(i,j), EVRI(i,j)
            enddo
            write(6,*)
          enddo
          write(6,*)
        endif
      endif
  200 format( 9(:,3X,F6.2,' + i*',F6.2) )

      return
      end subroutine sorteigs


!=======================================================================
! Name        : linalg
! Author      : Prabal S. Negi, Mattias Brynjell-Rahkola
! Version     : last modification 2018.02.22
! Copyright   : GPL
! Description : set of linear algebra routines and wrappers for LAPACK
!======================================================================
!---------------------------------------------------------------------- 
!     
!     LAPACK interface for SVD.
!     Upon finishing, SIGMA contains the singular values, VMATX n left
!     singular vectors and VMATXT n transposed right singular vectors.
!
!     Dongarra et al. (1999)
!
      subroutine svd_wrapper(nrow,ncol,dtype)

      implicit none

      include 'SIZE'
      include 'OTD'
      include 'WLAPACK'      
      
      character*1 jobz    !  Specifies options for computing all or part of the matrix U:
                          ! 'A':  all M columns of U and all N rows of V**T are
                          !    returned in the arrays U and VMATXT;
                          ! 'S':  the first min(M,N) columns of U and the first
                          !   min(M,N) rows of V**T are returned in the arrays U and VMATXT;
                          ! 'O':  If M >= N, the first N columns of U are overwritten
                          !   on the array A and all rows of V**T are returned in the array VMATXT;
                          !   otherwise, all columns of U are returned in the array U 
                          !   and the first M rows of V**T are overwritten in the array A;
                          ! 'N':  no columns of U or rows of V**T are computed.

      integer info        ! = 0:  successful exit.
                          ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                          ! > 0:  DBDSDC did not converge, updating process failed.                          

      character*1 job     !  Specifies for which problem the reciprocal condition numbers
                          !  should be computed:
                          ! = 'E':  the eigenvectors of a symmetric/Hermitian matrix;
                          ! = 'L':  the left singular vectors of a general matrix;
                          ! = 'R':  the right singular vectors of a general matrix.
      
      integer m, n, lda, ldu, ldvt,
     $     nrow, ncol
      real u,
     $     epsmch, serrbnd
      character*3 dtype
      character*23 fname

      integer i
      real dlamch       ! function

!     Define LAPACK-variables
      jobz = 'O'                ! return the n first singular vectors 
      m    = nrow               ! no rows in matrix
      n    = ncol               ! no columns in matrix
      lda  = LPERT              ! leading dimension of the matrix
      u    = 0.0                ! left singular vectors (not referenced)
      ldu  = 1                  ! leading dimension of u
      ldvt = LPERT              ! leading dimension of VMATXT

!     VMATX                      ! compute SVD for this matrix.
!     SIGMA                     ! Singular values of VMATX      

!     Compute SVD in double precision with divide-and-conquer
      call dgesdd(jobz,m,n,VMATX,lda,OSIGMA,u,ldu,VMATXT,ldvt,
     $     RWORK,LWORKR,IWORK,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: DBDSDC did not converge, updating process failed.'
     $        , info
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGESDD: successful exit!'
         if (nid.eq.0) write(6,*) '         Optimal LWORKR=',
     $        int(RWORK(1)), LWORKR
      endif

      epsmch = dlamch('E')      ! Machine epsilon in double precision
      if (nid.eq.0) write(6,*) 'Relative machine precision:', epsmch

!     Compute reciprocal condition numbers for singular vectors in double precision
      job = 'l'
      call ddisna(job,m,n,OSIGMA,RCL,info)
!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DDISNA: successful exit!'
      endif

      job = 'r'
      call ddisna(job,m,n,OSIGMA,RCR,info)
!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DDISNA: successful exit!'
      endif
            
!     Output singular values and estimated error bounds
      write(fname,'(A3,A20)') dtype,'_singular_values.txt'
      if (nid.eq.0) then
         open(61,file=fname,action='write',status='unknown')
         write(61,'(TR2,A2,TR18,A6,TR16,A8,TR16,A8,TR16,A8)')
     $        'I;','sigma;','serrbnd;','uerrbnd;','verrbnd'

!     see ch.4.9 Dongarra et al. (1999)
!         serrbnd = epsmch*SIGMA(1)
         do i=1,min(m,n)
            write(61,'(I4,4G24.16)') i, OSIGMA(i), serrbnd,
     $           serrbnd/RCL(i), serrbnd/RCR(i)            
         enddo
         close(61)
      endif
      
      return
      end
!----------------------------------------------------------------------
!      
!     LAPACK interface for the non-symmetric eigenvalue solver.
!     Upon finishing, RITZR and RITZI contain the real and imaginary part,
!     of the computed eigenvalues. Complex conjugate pairs of the
!     eigenvalues appear with the eigenvalue having the positive
!     imaginary part first.
!     The corresponding eigenvectors are stored in EVEC. If the j:th and
!     (j+1):th eigenvalue form a complex conjugate pair, then:
!     v(j) = EVEC(:,j)+i*EVEC(:,j+1), v(j+1) = EVEC(:,j)-i*EVEC(:,j+1)
!
!     Dongarra et al. (1999)
!
      subroutine eig_wrapper(n,kind)

      implicit none

      include 'SIZE'
      include 'OTD'
      include 'WLAPACK'
      
      character*1 jobvl, jobvr, kind
      integer lda, ldvl, ldvr, info,
     $     n,
     $     i, i0

!     Input parameter 'kind' determines whether left and/or right
!     eigenvectors should be computed
!      
!     Define LAPACK-variables
      if (kind.eq.'r') then
         jobvl = 'N'            ! don't compute left eigenvectors
         jobvr = 'V'            ! compute right eigenvectors
      elseif (kind.eq.'l') then
         jobvl = 'V'            ! compute left eigenvectors
         jobvr = 'N'            ! don't compute right eigenvectors
      elseif (kind.eq.'b') then
         jobvl = 'V'            ! compute left eigenvectors
         jobvr = 'V'            ! compute right eigenvectors
      else
         if (nid.eq.0)
     $        write(6,*) 'ERROR: choose left/right/both eigenvectors',
     $        kind
      endif
      lda   = npert              ! leading dimension of LR
      ldvl  = npert              ! leading dimension of EVECL
      ldvr  = npert              ! leading dimension of EVECR
!     Compute the eigenvalues/-vectors in double precision
      call dgeev(jobvl,jobvr,n,LR,lda,EIGR,EIGI,EVL,ldvl,
     $        EVR,ldvr,RWORK,LWORKR,info)
         
!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) then
            write(6,*) 'ERROR: the QR algorithm failed.', info
            write(6,*) '         Converged eigenvalues:'
            i0 = info+1
            do i=i0,n
               write(6,*) EIGR(i), EIGI(i)
            enddo
         endif
         
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGEEV: successful exit!'
         if (nid.eq.0) write(6,*) '        Optimal LWORKR=',
     $        int(RWORK(1)), LWORKR
      endif
      
      return
      end
!----------------------------------------------------------------------      
!                                              
!     LAPACK interface for determining the reciprocal condition number
!     for a square (nflds x nflds) input matrix. The choice of
!     condition number is determined by norm.
!
!     Note: nflds <= LPERT
!      
      subroutine cond_wrapper(rcond,anorm,nflds,norm)

      implicit none

      include 'SIZE'
      include 'OTD'
      include 'WLAPACK'

      character*1 norm    ! Specifies whether the 1-norm condition number
                          ! or the infinity-norm condition number is requested:
                          ! = '1' or 'O': 1-norm
                          ! = 'I'       : Infinity-norm

      integer info        ! = 0:  successful exit.
                          ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                          ! > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                          !       has been completed, but the factor U is exactly singular,
                          !       and division by zero will occur if it is used to solve a
                          !       system of equations. (For DGETRF only)
      
      integer m, n, lda,
     $     nflds
      real anorm, rcond
      real dlange       ! function
      
      if (nflds.ge.LPERT) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: nflds>=LPERT. Increase LPERT.', nflds, LPERT
         call exitt
      endif
      
!     Define LAPACK-variables
      m     = nflds             ! no rows in matrix
      n     = nflds             ! no columns in matrix/the order of the matrix
      lda   = LPERT              ! leading dimension of the matrix

!     Compute the selected norm of the input matrix in double precision
      anorm = dlange(norm,m,n,LR,lda,RWORK)
         
      call copy(ALU,LR,LPERT*LPERT)
!     LU factorize the matrix A in double precision
!     NOTE: Elements 1:min(m,n) of IWORK contain the pivot indices
      call dgetrf(m,n,ALU,lda,IWORK,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) write(6,*) 'WARNING: U(i,i) is exactly zero.',
     $        info
      else
         if (nid.eq.0) write(6,*) 'DGETRF: successful exit!'
      endif
      
!     Compute the resiprocal of the condition number in double precision
      call dgecon(norm,n,ALU,lda,anorm,rcond,RWORK,IWORK,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argument had an illegal value.', abs(info)
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGECON: successful exit!'
      endif

      return
      end
!----------------------------------------------------------------------

!=======================================================================
*
*     Auxiliary routine: printing eigenvalues.
*
      SUBROUTINE PRINT_EIGENVALUES( DESC, N, WR, WI )
      CHARACTER*(*)    DESC
      INTEGER          N
      REAL             WR( * ), WI( * )
*
      REAL             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO J = 1, N
         IF( WI( J ).EQ.ZERO ) THEN
            WRITE(*,9998,ADVANCE='NO') WR( J )
         ELSE
            WRITE(*,9999,ADVANCE='NO') WR( J ), WI( J )
         END IF
      END DO
      WRITE(*,*)
*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
*
*     Auxiliary routine: printing eigenvectors.
*
      SUBROUTINE PRINT_EIGENVECTORS( DESC, N, WI, V, LDV )
      CHARACTER*(*)    DESC
      INTEGER          N, LDV
      REAL             WI( * ), V( LDV, * )
*
      REAL             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, N
         J = 1
         DO WHILE( J.LE.N )
            IF( WI( J ).EQ.ZERO ) THEN
               WRITE(*,9998,ADVANCE='NO') V( I, J )
               J = J + 1
            ELSE
               WRITE(*,9999,ADVANCE='NO') V( I, J ), V( I, J+1 )
               WRITE(*,9999,ADVANCE='NO') V( I, J ), -V( I, J+1 )
               J = J + 2
            END IF
         END DO
         WRITE(*,*)
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END 
