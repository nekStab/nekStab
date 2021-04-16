!-----------------------------------------------------------------------





      subroutine inner_product(alpha, px,py,pz,pp,pt, qx,qy,qz,qp,qt)

!     This function provides the user-defined inner product to be used throughout
!     the computation.
!
!     INPUTS
!     ------
!
!     px, py, pz : nek arrays of size lt = lx1*ly1*lz1*lelt.
!     Arrays containing the velocity fields of the first vector.
!
!     pp : nek array of size lt2 = lx2*ly2*lz2*lelt
!     Array containing the pressure field of the first vector.
!
!     qx, qy, qz : nek arrays of size lt = lx1*ly1*lz1*lelt.
!     Arrays containing the velocity fields of the second vector.
!
!     qp : nek array of size lt2 = lx2*ly2*lz2*lelt
!     Array containing the pressure field of the second vector.
!
!     RETURN
!     ------
!
!     alpha : real
!     Value of the inner-product alpha = <p, q>.
!
!     Last edit : April 2nd 2020 by JC Loiseau.

      implicit none
      include "SIZE"
      include "TOTAL"

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt), intent(in) :: px, py, pz, pt
      real, dimension(lt), intent(in) :: qx, qy, qz, qt
      real, dimension(lt2), intent(in) :: pp, qp !not used

      real, intent(out) :: alpha
      real :: glsc3
      integer :: n

      n = nx1 * ny1 * nz1 * nelt

      alpha = glsc3(px, bm1s, qx, n) + glsc3(py, bm1s, qy, n)
      if (if3d) alpha = alpha + glsc3(pz, bm1s, qz, n)
      if (ifheat) alpha = alpha + glsc3(pt, bm1s, qt, n)

      return
      end subroutine inner_product





!--------------------------------------------------------------------------





      subroutine norm(qx, qy, qz, qp, qt, alpha) ! Compute vector norm
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter               :: lt  = lx1*ly1*lz1*lelt
      integer, parameter               :: lt2 = lx2*ly2*lz2*lelt
      real, intent(in), dimension(lt)  :: qx, qy, qz, qt
      real, intent(in), dimension(lt2) :: qp
      real, intent(out)                :: alpha

      call inner_product(alpha, qx,qy,qz,qp,qt, qx,qy,qz,qp,qt)
      alpha = dsqrt(alpha)

      return
      end subroutine norm





!----------------------------------------------------------------------





      subroutine normalize(qx, qy, qz, qp, qt, alpha)

!     This function normalizes the state vector [qx, qy, qz, qp]^T where
!     qx, qy and qz are the streamwise, cross-stream and spanwise velocity
!     components while qp is the corresponding pressure field.
!
!     INPUTS / OUTPUTS
!     ----------------
!
!     qx, qy, qz : nek arrays of size lt = lx1*ly1*lz1*lelt.
!     Arrays storing the velocity components.
!
!     qp : nek array of size lt2 = lx2*ly2*lz2*lelt
!     Array storing the corresponding pressure field.
!
!     alpha : real
!     Norm of the vector.
!
!     Last edit : April 2nd 2020 by JC Loiseau.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                  :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                  :: lt2 = lx2*ly2*lz2*lelt
      real, dimension(lt), intent(inout)  :: qx, qy, qz, qt
      real, dimension(lt2), intent(inout) :: qp
      real, intent(out)                   :: alpha
      real                                :: beta

!     --> Compute the user-defined norm.
      call norm(qx, qy, qz, qp, qt, alpha)
      beta = 1.0d0/alpha

!     --> Normalize the vector.
      call nopcmult(qx, qy, qz, qp, qt, beta)

      return
      end subroutine normalize

!-----------------------------------------------------------------------


      subroutine krylov_schur
!     use KS
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

!     -----Krylov basis V for the projection M*V = V*H -----
      real, allocatable,dimension(:,:) :: qx, qy, qz, qt
      real, allocatable,dimension(:,:) :: qp

!     ----- Upper Hessenberg matrix -----
      real, allocatable,dimension(:,:) :: H
      real, allocatable,dimension(:,:) :: b_vec

!     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----
      complex*16, allocatable,dimension(:) :: vals
      complex*16, allocatable,dimension(:,:) :: vecs

      real, allocatable, dimension(:) :: residual

!     ----- Miscellaneous -----
      real wo1(lt),wo2(lt),wo3(lt)
      common /ugrad/ wo1,wo2,wo3

      integer :: mstart, cnt
!     real, dimension(k_dim)             :: residual
      real                               :: alpha, beta, glsc3
      logical                            :: converged
      integer                            :: n, i, j
      character(len=30)                  :: filename

!     ----- Allocate arrays -----
      allocate(qx(lt,k_dim+1), qy(lt,k_dim+1), qz(lt,k_dim+1), qt(lt,k_dim+1),qp(lt2,k_dim+1))
      allocate(H(k_dim+1,k_dim),b_vec(1,k_dim),vals(k_dim),vecs(k_dim,k_dim),residual(k_dim))

      n      = nx1*ny1*nz1*nelt
      time   = 0.0d0
      H(:,:)  = 0.0d0
      b_vec  = 0.0d0
      residual = 0.0d0

      call oprzero(wo1, wo2, wo3)
      do i = 1, k_dim+1
         call noprzero(qx(:, i), qy(:, i), qz(:, i), qp(:, i), qt(:, i))
      enddo

!     ----- Loading baseflow from disk (optional) -----

      if(ifldbf)then            !skip loading if single run
         write(filename,'(a,a,a)')'BF_',trim(SESSION),'0.f00001'
         if(nid.eq.0)write(*,*)'Loading base flow: ',filename
         call load_fld(filename)
      endif

!     ----- Save baseflow to disk (recommended) -----
      call opcopy(ubase,vbase,wbase,vx,vy,vz)
      if(ifheat) call copy(tbase,t(1,1,1,1,1),n)
!     call outpost(vx,vy,vz,pr,t,'BF_') !outpost for sanity check

!     ----- Prepare stability parameters -----

      call krylov_schur_prepare

!     ----- First vector (new from noise or restart) -----

      if (uparam(2).eq.0) then

         if(nid.eq.0)write(6,*)'Starting first Arnoldi decomposition...'

!     ----- Creates seed vector for the Krylov subspace -----

         if(ifseed_nois)then    ! noise as initial seed

            if(nid.eq.0)write(6,*)'Filling fields with noise...'
            call add_noise(vxp(:,1),vyp(:,1),vzp(:,1),tp(:,1,1))
            call normalize(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1),alpha)
!     call outpost(vxp,vyp,vzp,pr,tp,'SS_') !outpost for sanity check
            call matvec(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1),
     $           vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1))

         elseif(ifseed_symm)then ! symmetry initial seed

            if(nid.eq.0)write(6,*)'Enforcing symmetric seed perturb...'
            call add_symmetric_seed(vxp(:,1),vyp(:,1),vzp(:,1),tp(:,1,1))
!     call outpost(vxp,vyp,vzp,pr,t,'SS_') !outpost for sanity check

         elseif(ifseed_load)then ! loading initial seed (e.g. Re_ )

            write(filename,'(a,a,a)')'Re_',trim(SESSION),'0.f00001'
            if(nid.eq.0)write(*,*)'Load real part of leading mode as seed: ',filename
            call load_fld(filename)
            call nopcopy(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1), vx,vy,vz,pr,t(1,1,1,1,1))

         else

            call opcopy(vxp(:,1),vyp(:,1),vzp(:,1),ubase,vbase,wbase)
            if(ifheat) call copy(tp(:,1,1),tbase,n)

         endif

!     ----- Normalized to unit-norm -----

         call normalize(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1),alpha)

         mstart = 1; istep = 1; time = 0.0d0

         call nopcopy(qx(:,1),qy(:,1),qz(:,1),qp(:,1),qt(:,1), vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1))

         call whereyouwant('KRY',1)
         time = 0.0d0
         call outpost(qx(:,1), qy(:,1), qz(:,1), qp(:,1), qt(:,1), 'KRY')

      elseif(uparam(2).gt.0)then

         mstart = int(uparam(2))

         if(nid.eq.0)then
            write(6,*)'Restarting from:',mstart
            write(6,'(a,a,i4.4)')' Loading Hessenberg matrix: HES',trim(SESSION),mstart
            write(filename,'(a,a,i4.4)')'HES',trim(SESSION),mstart

            open(67, file=trim(filename), status='unknown', form='formatted')
            if (k_dim.lt.mstart) then !subsampling
               do i = 1,k_dim+1
                  do j = 1,mstart
                     if (j.le.k_dim)read(67,"(1E15.7)") H(i,j)
                  enddo
               enddo
            else
               read(67, *)  (( H(i, j), j=1, mstart ), i = 1, mstart+1 )
            endif
            close(67)
            write(6,*)'Broadcast H matrix to all procs...'
         endif                  !nid.eq.0
         call bcast(H,(k_dim+1)*k_dim*wdsize) !broadcast H matrix to all procs

!     if(nid.eq.2)then !-> debug only
!     write(filename,'(a,a,i4.4)')'HESloaded',trim(SESSION),mstart
!     write(6,*) ''
!     open(67,file=trim(filename),status='unknown',form='formatted')
!     write(67,*) beta
!     do i = 1,mstart+1
!     do j = 1,mstart
!     write(67,*) H(i,j)
!     enddo
!     enddo
!     close(67)
!     endif

         mstart=mstart+1        !careful here!
         call load_files(qx, qy, qz, qp, qt, mstart, k_dim+1, 'KRY')
!     k_dim+1 is the dim of V_x
         if(nid.eq.0) write(6,*)'Restart fields loaded to memory!'

      endif

!     ======================================
!     =====                            =====
!     ===== Krylov-Schur decomposition =====
!     =====                            =====
!     ======================================

      schur_cnt = 0
      converged = .false.
      do while ( .not. converged )
!     --> Arnoldi factorization.
         call arnoldi_factorization(qx, qy, qz, qp, qt, H, mstart, k_dim, k_dim)

!     --> Compute the eigenspectrum of the Hessenberg matrix.
         call eig(H(1:k_dim, 1:k_dim), vecs, vals, k_dim)

!     --> Check the residual of the eigenvalues.
         residual = abs(H(k_dim+1, k_dim) * vecs(k_dim, :))
         cnt = count(residual .lt. eigen_tol)
         if (nid .eq. 0) write(6, *) 'converged eigenvalues:',cnt

!     --> Select whether to stop or apply Schur condensation depending on schur_tgt.
         select case (schur_tgt)

!     --> k-step Arnoldi factorization completed.
         case (:0)
            converged = .true.
            if (nid .eq. 0) write(6, *) 'Arnoldi factorization completed.'

!     --> Krylov-Schur factorization.
         case (1:)
            if (cnt .ge. schur_tgt) then ! Krylov-Schur factorization completed.
               converged = .true.
            else                ! Apply Schur condensation before restarting the factorization.
               schur_cnt = schur_cnt + 1
               if (nid .eq. 0) write(6, *) 'Starting Schur condensation phase.',schur_cnt
               call schur_condensation(mstart, H, qx, qy, qz, qp, qt, k_dim)
            endif

         end select

      enddo

!     ----- Output all the spectrums and converged eigenmodes -----
      if(nid.eq.0)write(6,*)'Exporting modes...'
      if(nid.eq.0)print *,''
      call outpost_ks(vals, vecs, qx, qy, qz, qp, qt, residual)

      if(nid.eq.0)write(6,*)'converged eigenmodes:',cnt
      if(nid.eq.0)write(6,*)'Eigenproblem solver finished.'

!     --> Deallocation.
      deallocate(qx, qy, qz, qt, qp)
      deallocate(H, b_vec, vals, vecs)

      return
      end subroutine krylov_schur














c-----------------------------------------------------------------------





      subroutine schur_condensation(mstart, H, qx, qy, qz, qp, qt, ksize)

      implicit none
      include 'SIZE'
      include 'TOTAL'

c     =================================================
c     =====                                       =====
c     ===== Declaration of the required variables =====
c     =====                                       =====
c     =================================================

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

      integer :: ksize

c     ----- Krylov basis V for the projection M*V = V*H -----

      real, dimension(lt,ksize+1)        :: qx, qy, qz, qt
      real, dimension(lt2,ksize+1)       :: qp

c     ----- Upper Hessenberg matrix -----

      real, dimension(ksize+1, ksize)     :: H
      real, dimension(ksize)           :: b_vec

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(ksize)       :: vals

c     ----- Miscellaneous -----

      integer :: mstart,n
      logical, dimension(ksize)          :: selected

c     ----- Schur and Hessenberg decomposition -----

      real, dimension(ksize, ksize)       :: vecs

      n  = nx1*ny1*nz1*nelt

!     --> Initialize arrays.
      b_vec = 0.0D0 ; b_vec(ksize) = H(ksize+1, ksize)
      vals = ( 0.0D0 , 0.0D0 ) ; vecs = 0.0D0

!     --> Perform the Schur decomposition.
      call schur(H(1:ksize, 1:ksize), vecs, vals, ksize)

!     --> Partition the eigenvalues in wanted / unwanted.
      call select_eigenvalues(selected, mstart, vals, schur_del, schur_tgt, ksize)
      if ( nid.eq.0 ) write(6,*) mstart, 'Ritz eigenpairs have been selected.'

!     --> Re-order the Schur decomposition based on the partition.
      call ordschur(H(1:ksize, 1:ksize), vecs, selected, ksize)

!     --> Zero-out the unwanted blocks of the Schur matrix.
      H(1:mstart,mstart+1:ksize) = 0.0D0
      H(mstart+1:ksize+1, :)     = 0.0D0

!     --> Re-order the Krylov basis accordingly.
      qx(:, 1:ksize) = matmul(qx(:, 1:ksize), vecs)
      qy(:, 1:ksize) = matmul(qy(:, 1:ksize), vecs)
      if (if3d) qz(:, 1:ksize) = matmul(qz(:, 1:ksize), vecs)
      if (ifpo) qp(:, 1:ksize) = matmul(qp(:, 1:ksize), vecs)
      if (ifheat) qt(:, 1:ksize) = matmul(qt(:, 1:ksize), vecs)

!     --> Update the Schur matrix with b.T @ Q corresponding to
!     the residual beta in the new basis.
      b_vec = matmul(b_vec, vecs)
      H(mstart+1, :) = b_vec

!     --> Add the last generated Krylov vector as the new starting one.
      mstart = mstart + 1

      call nopcopy(qx(:,mstart),  qy(:,mstart),  qz(:,mstart),  qp(:,mstart),  qt(:,mstart),
     $             qx(:,ksize+1), qy(:,ksize+1), qz(:,ksize+1), qp(:,ksize+1), qt(:,ksize+1))

      return
      end subroutine schur_condensation





c----------------------------------------------------------------------





      subroutine prepare_baseflow
!     preparing for Floquet upgrade

      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,i
      character(len=30) filename
      logical ifto_sav, ifpo_sav
      n  = nx1*ny1*nz1*nelt

      if(    uparam(3).eq.0.1)then !steady prescribe

!     if(nid.eq.0)write(6,*)'Copying base flow!'
         call opcopy(vx,vy,vz,ubase,vbase,wbase)
         if(ifheat) call copy(t(1,1,1,1,1), tbase, n)

!     ----- Save baseflow to disk (recommended) -----
         call opcopy(ubase,vbase,wbase,vx,vy,vz)
         if(ifheat) call copy(tbase,t(1,1,1,1,1),n)
         ifto_sav=ifto;ifpo_sav=ifpo
         ifto=.true.;ifpo=.true.
         call outpost(vx,vy,vz,pr,t,'BFL') !outpost for sanity check
         ifto=ifto_sav;ifpo=ifpo_sav

      elseif(uparam(3).eq.0.2)then !steady load


!     if(.not.iffloq .and. (.not. ifbfcv .or. ifldbf))then
         write(filename,'(a,a,a)')'BF_',trim(SESSION),'0.f00001'
         if(nid.eq.0)write(*,*)'Loading base flow: ',filename
         call load_fld(filename)

!     ----- Save baseflow to disk (recommended) -----
         call opcopy(ubase,vbase,wbase,vx,vy,vz)
         if(ifheat) call copy(tbase,t(1,1,1,1,1),n)
         ifto_sav=ifto;ifpo_sav=ifpo
         ifto=.true.;ifpo=.true.
         call outpost(vx,vy,vz,pr,t,'BFL') !outpost for sanity check
         ifto=ifto_sav;ifpo=ifpo_sav


      elseif(uparam(3).eq.0.3)then !steady compute

         ifbase = .true.
         if(nid.eq.0)write(6,*)'Running DNS alongside stability!'

      elseif(uparam(3).eq.1.1)then !periodic prescribe

         if(nid.eq.0)write(6,*)'Prescribed periodic baseflow!'

      elseif(uparam(3).eq.1.2)then !periodic reconstruct

      elseif(uparam(3).eq.1.21)then !periodic load all

      elseif(uparam(3).eq.1.3)then !periodic compute

         ifbase = .true.
         if(nid.eq.0)write(6,*)'Running DNS alongside stability!'

      endif


      return
      end subroutine prepare_baseflow



c----------------------------------------------------------------------

      subroutine krylov_schur_prepare

!     This
!
!     INPUT
!     -----
!
!     RETURNS
!     -------
!
!     Last edit : March 26th 2021 by RAS Frantz.

      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'

!     forcing npert to unity
      if(param(31).gt.1)then
         write(6,*)'nekStab not ready for npert>1 -- jp loops NOT implemented! STOPPING!'
         call nek_end
      endif
      param(31) = 1 ; npert = param(31)

!     force OIFS deactivation
      if(ifchar.and.nid.eq.0)write(6,*)'OIFS not working with linearized and adjoint -> turning OFF!'
      ifchar = .false.
      call bcast(ifchar, lsize)

!     enforce CFL target for EXTk
      if( param(26).gt.0.5 )then
         if(nid.eq.0)write(6,*)'reducing target CFL to 0.5!'
         param(26)=0.50d0 ; ctarg = param(26)
      endif

      if(param(10).gt.0)then
!     if param(10)=endTime=0 -> param(11) = numSteps
         call compute_cfl(dt,vx,vy,vz,1.0) ! dt=1 ! vx at this point is base flow
         dt = ctarg/dt
         nsteps = ceiling(param(10)/dt)
         dt = param(10)/nsteps
         if(nid.eq.0)write(6,*)'endTime specified! computing CFL from BASE FLOW!'
         if(nid.eq.0)write(6,*)' computing timeStep dt=',dt
         if(nid.eq.0)write(6,*)' computing numSteps=',nsteps
         if(nid.eq.0)write(6,*)' sampling period =',nsteps*dt
         param(12) = dt
      endif

!     deactivate variable time step! !freeze dt
      param(12) = -abs(param(12))

!     broadcast all parameters to processors
      call bcast(param,200*wdsize)

      return
      end subroutine krylov_schur_prepare














c----------------------------------------------------------------------




      subroutine outpost_ks(vals, vecs, qx, qy, qz, qp, qt, residual)
!     outposting vectors
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

c     ----- Krylov basis V for the projection M*V = V*H -----

      real wo1(lt),wo2(lt),wo3(lt),vort(lt,3)
      real, dimension(lt,k_dim+1)        :: qx, qy, qz, qt
      real, dimension(lt2,k_dim+1)       :: qp

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: vals
      complex*16, dimension(k_dim,k_dim) :: vecs

c     ----- Arrays for the storage/output of a given eigenmode of the NS operator -----

      complex*16, dimension(lt)          :: fp_cx, fp_cy, fp_cz, fp_ct
      complex*16, dimension(lt2)         :: fp_cp

c     ----- Miscellaneous -----
      integer :: n, i

      real                               :: sampling_period
      real, dimension(k_dim)             :: residual
      real                               :: alpha, alpha_r, alpha_i, beta
      complex :: log_transform
      logical ifto_sav, ifpo_sav

      character(len=80) :: filename
      character(len=20) :: fich1,fich2,fich3,fmt2,fmt3,fmt4,fmt5,fmt6
      character(len=3)  :: nre,nim,nrv
      integer :: outposted

      sampling_period = dt*nsteps
      n = nx1*ny1*nz1*nelt

c     ----- Output all the spectrums and converged eigenmodes -----

!     evop defined in matrix_vector_product

      write(nre,'(A,A)')trim(evop),'Re'
      write(nim,'(A,A)')trim(evop),'Im'
      write(nrv,'(A,A)')trim(evop),'Rv'

      write(fich1,'(A,A,A)')'Spectre_H',trim(evop),'.dat' !,trim(SESSION)
      write(fich2,'(A,A,A)')'Spectre_NS',trim(evop),'.dat' !,trim(SESSION)
      write(fich3,'(A,A,A)')'Spectre_NS',trim(evop),'_conv.dat' !,trim(SESSION)

      if(nid.eq.0) then
         open(unit=10, file=fich1, form='formatted')
         open(unit=20, file=fich2, form='formatted')
         open(unit=30, file=fich3, form='formatted')
      endif

      outposted = 0
      do i = 1, k_dim

         time = real(i)         !here for the outposted files have different times

         if(nid.eq.0) then

!     --> Outpost the eigenspectrum of the Hessenberg matrix.
            write(10,"(3E15.7)") real(vals(i)), aimag(vals(i)), residual(i)

!     --> Outpost the log-transform spectrum (i.e. eigenspectrum of the linearized Navier-Stokes operator).
            write(20, "(3E15.7)") real(log_transform(vals(i))) / sampling_period,
     $           aimag(log_transform(vals(i))) / sampling_period,
     $           residual(i)

         endif

         if (residual(i) .lt. eigen_tol .and. outposted.lt.maxmodes ) then !just the converged ones
         outposted = outposted + 1

!     --> Outpost only the converged part of the log-transformed spectrum.
         if(nid.eq.0)write(30, "(2E15.7)") real(log_transform(vals(i))) / sampling_period,
     $        aimag(log_transform(vals(i))) / sampling_period

c     ----- Computation of the corresponding eigenmode -----
         fp_cx = matmul(qx(:, 1:k_dim), vecs(:, i))
         fp_cy = matmul(qy(:, 1:k_dim), vecs(:, i))
         if (if3d) fp_cz = matmul(qz(:, 1:k_dim), vecs(:, i))
         if (ifpo) fp_cp = matmul(qp(:, 1:k_dim), vecs(:, i))
         if (ifheat) fp_ct = matmul(qt(:, 1:k_dim), vecs(:, i))

c     ----- Normalization to be unit-norm -----
c     Note: volume integral of FP*conj(FP) = 1.
         call norm(real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), real(fp_ct), alpha_r)
         call norm(aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), aimag(fp_ct), alpha_i)
         alpha = alpha_r**2 + alpha_i**2
         beta = 1.0d0/sqrt(alpha)

c     ----- Output the real part -----
         call nopcopy(vx,vy,vz,pr,t(1,1,1,1,1), real(fp_cx),real(fp_cy),real(fp_cz),real(fp_cp),real(fp_ct))
         call nopcmult(vx,vy,vz,pr,t(1,1,1,1,1), beta)
         call outpost(vx, vy, vz, pr, t(1,1,1,1,1), nRe)

c     ----- Output the imaginary part -----
         call nopcopy(vx,vy,vz,pr,t(1,1,1,1,1), aimag(fp_cx),aimag(fp_cy),aimag(fp_cz),aimag(fp_cp),aimag(fp_ct))
         call nopcmult(vx,vy,vz,pr,t(1,1,1,1,1), beta)
         call outpost(vx, vy, vz, pr, t(1,1,1,1,1), nIm)

         if(ifvor)then
c     ----- Output vorticity from real part -----
            call oprzero(wo1, wo2, wo3)
            call comp_vort3(vort, wo1, wo2, vx, vy, vz)

            ifto_sav = ifto; ifpo_sav = ifpo; ifto = .false.; ifpo = .false.
            call outpost(vort(1,1), vort(1,2), vort(1,3), pr, t, nRv)
            ifto = ifto_sav ; ifpo = ifpo_sav
         endif

      endif

      enddo

      if (nid .eq. 0) then

         close(10) ; close(20) ;  close(30)
!
         write(fmt2,'("(A,I16)")')
         write(fmt3,'("(A,F16.4)")')
         write(fmt4,'("(A,F16.12)")')
         write(fmt5,'("(A,E15.7)")') ! max precision
         write(fmt6,'("(A,E13.4)")') ! same as hmhlz
!
         write(filename,'(A,A,A)')'Spectre_',trim(evop),'.info'
!     write(filename,"(',I7.7,'.info')") itime/ioutput
         open (844,file=filename,action='write',status='replace')

         write(844,'(A,A)')'Nek5000 version:',NVERSION
         write(844,'(A,A)')'nekStab version:',NSVERSION
         write(844,'(A)')  '[mesh]'
         write(844,fmt2)   'lx1=             ',lx1
         write(844,fmt2)   'polyOrder N=     ',lx1-1
         write(844,fmt2)   'tot elemts=      ',nelgv
         write(844,fmt2)   'tot points=      ',nelgv*(lx1)**ldim
         write(844,fmt2)   'MPI ranks=       ',np
         write(844,fmt2)   'e/rank=          ',nelgv/np
         write(844,'(A)')  '[userParams]'
         write(844,fmt4)   'uparam01=        ',uparam(01)
         write(844,fmt4)   'uparam02=        ',uparam(02)
         write(844,fmt4)   'uparam03=        ',uparam(03)
         write(844,fmt4)   'uparam04=        ',uparam(04)
         write(844,fmt4)   'uparam05=        ',uparam(05)
         write(844,fmt4)   'uparam06=        ',uparam(06)
         write(844,fmt4)   'uparam07=        ',uparam(07)
         write(844,fmt4)   'uparam08=        ',uparam(08)
         write(844,fmt4)   'uparam09=        ',uparam(09)
         write(844,fmt4)   'uparam10=        ',uparam(10)
         write(844,'(A)')  '[solver]'
         write(844,fmt3)   'ctarg=           ',ctarg
         write(844,fmt2)   'nsteps=          ',nsteps
         write(844,fmt5)   'dt=              ',dt
         write(844,fmt3)   'Re=              ',1.0/param(2)
         write(844,fmt6)   'residualTol PRE= ',param(21)
         write(844,fmt6)   'residualTol VEL= ',param(22)
         if(ifheat)then
            write(844,fmt6)   'residualTol TEM= ',param(22)
            write(844,fmt3)   'Pe=              ',1.0/param(8)
         endif
         write(844,'(A)')  '[eigensolver]'
         write(844,fmt4)   'sampling period =',sampling_period
         write(844,fmt2)   'k_dim=           ',k_dim
         write(844,fmt6)   'eigentol=        ',eigen_tol
         write(844,fmt2)   'schur_target=    ',schur_tgt
         write(844,fmt3)   'schur_del=       ',schur_del
         write(844,fmt2)   'schur iterations=',schur_cnt
         write(844,fmt2)   'outposted=       ',outposted
         close(844)
      endif

      return
      end subroutine outpost_ks





c-----------------------------------------------------------------------





      subroutine select_eigenvalues(selected, cnt, vals, delta, nev, n)

!     This function selects the eigenvalues to be placed in the upper left corner
!     during the Schur condensation phase.
!
!     INPUTS
!     ------
!
!     vals : n-dimensional complex array.
!     Array containing the eigenvalues.

!     delta : real
!     All eigenvalues outside the circle of radius 1-delta will be selected.
!
!     nev : integer
!     Number of desired eigenvalues. At least nev+4 eigenvalues will be selected
!     to ensure "smooth" convergence of the Krylov-Schur iterations.
!
!     n : integer
!     Total number of eigenvalues.
!
!     RETURNS
!     -------
!
!     selected : n-dimensional logical array.
!     Array indicating which eigenvalue has been selected (.true.).
!
!     cnt : integer
!     Number of selected eigenvalues. cnt >= nev + 4.
!
!     Last edit : April 2nd 2020 by JC Loiseau.

      implicit none

!     ----- Input arguments -----
      integer :: nev
      integer :: n
      integer, dimension(n) :: idx
      complex*16, dimension(n) :: vals
      real :: delta

!     ----- Output argument -----
      logical, dimension(n) :: selected
      integer :: cnt

!     ----- Miscellaneous -----
      integer :: i

!     --> Sort eigenvalues based on their magnitudes (Note : the array vals itself is not sorted).
      do i = 1, n
         idx(i) = i
      enddo
      call quicksort2(n, abs(vals), idx)

!     --> Select eigenvalues closer to the unit circle.
      selected = abs(vals) .ge. (1.0d0-delta)

!     --> Select at least the nev+4 largest eigenvalues.
      selected(idx(n-(nev+3):n)) = .true.
      if (aimag(vals(idx(n-(nev+3)))) .eq. -aimag(vals(idx(n-(nev+4))))) then
      selected(idx(n-(nev+4))) = .true.
      endif

      cnt = count(selected)

      return
      end subroutine select_eigenvalues










!     ------------------------------------------------------------------------------------





      subroutine arnoldi_checkpoint(f_xr, f_yr, f_zr, f_pr, f_tr, H, k)

!     This function implements a fairly simple checkpointing procedure in case one
!     would need to restart the computation (e.g. in case of cluster shutdown).
!
!     INPUTS
!     ------
!
!     f_xr, f_yr, f_zr : nek arrays of size lt = lx1*ly1*lz1*lelt
!     Velocity components of the latest Krylov vector.
!
!     f_pr : nek array of size lt2 = lx2*ly2*lz2*lelt
!     Pressure field of the latest Krylov vector.
!
!     H : k+1 x k real matrix.
!     Current upper Hessenberg matrix resulting from the k-step Arnoldi factorization.
!
!     k : int
!     Current iteration of the Arnoldi factorization.
!
!     Last edit : April 3rd 2020 by JC Loiseau.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt), intent(in) :: f_xr, f_yr, f_zr, f_tr
      real, dimension(lt2), intent(in) :: f_pr

      integer, intent(in) :: k
      real, dimension(k+1,k), intent(in) :: H

      integer :: i, j, cnt
      real :: beta
      complex*16, dimension(k) :: vals
      complex*16, dimension(k,k) :: vecs
      real, dimension(k) :: residual
      character(len=80) filename
      complex :: log_transform

!     -->
      if(nid .eq. 0) then
         write(6, *)
         write(6, *) 'Outposting Krylov vector to'
      endif

!     --> Outpost the latest Krylov vector.
      call whereyouwant("KRY", k+1)
      time = time * k           !order in ParaView
      call outpost(f_xr, f_yr, f_zr, f_pr, f_tr, "KRY")

!     --> Compute the eigenvalues and eigenvectors of the current Hessenberg matrix.
      call eig(H(1:k, 1:k), vecs, vals, k)

!     --> Compute the residual (assume H results from classical Arnoldi factorization).
      residual = abs(H(k+1, k) * vecs(k, :))
      cnt = count(residual .lt. eigen_tol)

      if (nid.eq.0) then
!     --> Outpost the eigenspectrum and residuals of the current Hessenberg matrix.
!     write(filename, '(a, i4.4)') "H", k

         write(filename, '(A,A,A)') 'Spectre_H',evop,'.dat'
         write(6, *) 'Outposting eigenspectrum of current Hessenberg matrix to : ', filename

         open(67, file=trim(filename), status='unknown', form='formatted')
         write(67, '(3E15.7)') (real(vals(i)), aimag(vals(i)), residual(i), i=1, k)
         close(67)

!     --> Outpost the log-transform spectrum (i.e. eigenspectrum of the linearized Navier-Stokes operator).
!     write(filename, '(a, i4.4)') "S", k
         write(filename, '(A,A,A)') 'Spectre_NS',evop,'.dat'
         write(6, *) 'Outposting eigenspectrum of current log-transform spectrum matrix to : ', filename

         open(67, file=trim(filename), status='unknown', form='formatted')
         write(67, '(3E15.7)') (real(log_transform(vals(i))) / (dt*nsteps),
     $        aimag(log_transform(vals(i))) / (dt*nsteps),
     $        residual(i), i=1, k)
         close(67)

!     --> Outpost the Hessenberg matrix for restarting purposes (if needed).
         write(filename, '(a, a, i4.4)') 'HES', trim(SESSION), k
         write(6, *) 'Outposting current Hessenberg matrix to :', filename

         open(67, file=trim(filename), status='unknown', form='formatted')
         write(67, *) ((H(i, j), j=1, k), i=1, k+1)
         close(67)

!     --> Write to logfile the current number of converged eigenvalues.
         write(6, *) 'iteration converged and target :',cnt,'/',schur_tgt !keep small caps to ease grep
      endif


      !if(cnt.ge.schur_tgt)then
      !  if(nid.eq.0)write(6,*) 'Target reached! exporting and stopping'
      !  !call outpost_ks(vals, vecs, qx, qy, qz, qp, qt, residual)
      !  call nek_end
      !endif


      return
      end subroutine arnoldi_checkpoint

!     ------------------------------------------------------------------------------------
      function log_transform(x)
      implicit none

      complex :: log_transform,x
      log_transform = log(x)
      if (aimag(x) .eq. 0) log_transform = real(log_transform)

      return
      end function log_transform


      subroutine power_iteration

          implicit none
          include 'SIZE'
          include 'TOTAL'

          integer, parameter :: lt = lx1*ly1*lz1*lelt
          integer, parameter :: lt2 = lx2*ly2*lz2*lelt

          ! -->
          real, dimension(lt) :: qx, qy, qz, qt
          real, dimension(lt2) :: qp

          ! -->
          real, dimension(lt) :: fx, fy, fz, ft
          real, dimension(lt2) :: fp

          real :: alpha, beta
          integer :: n, i, j
          character(len=80) :: filename

          n = nx1 * ny1 * nz1 * nelt

          write(filename,'(a,a,a)')'BF_',trim(SESSION),'0.f00001'
          if(nid.eq.0)write(*,*)'Loading base flow: ',filename
          call load_fld(filename)
          call opcopy(ubase,vbase,wbase,vx,vy,vz)
          if(ifheat) call copy(tbase,t(1,1,1,1,1),n)

          call krylov_schur_prepare

          call add_noise(vxp(:,1),vyp(:,1),vzp(:,1),tp(:,1,1))
          call normalize(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,1,1),alpha)
          call nopcopy(qx, qy, qz, qp, qt, vxp(:,1), vyp(:,1), vzp(:,1), prp(:,1), tp(:,1,1))

          call outpost(qx, qy, qz, qp, qt, "PRT")

          do i = 1, 10
              call matvec(fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
              call inner_product(alpha, fx, fy, fz, fp, ft, qx, qy, qz, qp, qt)
              call norm(qx, qy, qz, qp, qt, beta)
              if (nid.EQ.0) write(*, *) "Rayleigh Quotient : ", alpha/beta, alpha, beta
              call nopcopy(qx, qy, qz, qp, qt, fx, fy, fz, fp, ft)
              call normalize(qx, qy, qz, qp, qt, alpha)
              call outpost(qx, qy, qz, qp, qt, "PRT")
          enddo

          call nek_end

          return
      end subroutine power_iteration
