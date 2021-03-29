!-----------------------------------------------------------------------





      subroutine arnoldi_factorization(qx, qy, qz, qp, qt, H, mstart, mend, ksize)

!     This function implements the k-step Arnoldi factorization of the linearized
!     Navier-Stokes operator. The rank k of the Arnoldi factorization is set as a user
!     parameter in x_SIZE.f (see parameter k_dim).
!     
!     INPUT
!     -----
!     
!     mstart : integer
!     Index at which to start the Arnoldi factorization. By default, it should be set to 1.
!     Note that it changes when the Arnoldi factorization is used as part of the Krylov-Schur
!     algorithm.
!     
!     mend : integer
!     Index at which to stop the Arnoldi factorization. By default, it should be set to kdim.
!     Note that it changes when the Arnoldi factorization is used as part of the GMRES solver
!     
!     ksize : integer
!     Size of the Krylov subspace (same as k_dim).
!     
!     RETURNS
!     -------
!     
!     qx, qy, qz : nek arrays of size (lx1*ly1*lz1*lelt, ksize).
!     Arrays containing the various Krylov vectors associated to each velocity component.
!     
!     qp : nek arrays of size (lx2*ly2*lz2*lelt, ksize)
!     Arrays containing the various Krylov vectors associated to the pressure field.
!     
!     H : k x k real matrix.
!     Upper Hessenberg matrix resulting from the Arnoldi factorization of the linearized
!     Navier-Stokes operator.
!     
!     Last edit : April 3rd 2020 by JC Loiseau.
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

!     ----- Miscellaneous -----
      real                               :: alpha
      integer                            :: mstart, mend, ksize !, mstep moved to NEKSTAB.inc
      integer                            :: n

!     ----- Timer -----
      real*8 :: eetime0,eetime1
      real   :: telapsed,tmiss,dnekclock

!     ----- Orthogonal residual f = w - (Q,w)*Q -----
      real, dimension(lt)                :: f_xr, f_yr, f_zr, f_tr
      real, dimension(lt2)               :: f_pr
c     
!     ----- Krylov basis V for the projection MQ = QH -----
      real, dimension(lt,ksize+1)        :: qx, qy, qz, qt
      real, dimension(lt2,ksize+1)       :: qp

!     ----- Upper Hessenberg matrix -----
      real, dimension(ksize+1, ksize)    :: H

      n  = nx1*ny1*nz1*nelt

!     --> Initialize arrays.
      call noprzero(f_xr, f_yr, f_zr, f_pr, f_tr)
      alpha = 0.0d0

!     --> Arnoldi factorization.
      do mstep = mstart, mend

         if (nid.eq.0) write(6,*) 'iteration current and total:', mstep , '/' , mend

         eetime0=dnekclock()

!     --> Matrix-vector product f = M * v (e.g. calling the linearized Navier-Stokes solver).
         call matrix_vector_product(f_xr, f_yr, f_zr, f_pr, f_tr, qx(:,mstep), qy(:,mstep), qz(:,mstep), qp(:,mstep), qt(:,mstep))

!     --> Update Hessenberg matrix and compute the orthogonal residual f.
         call update_hessenberg_matrix(H(1:mstep, 1:mstep), f_xr, f_yr, f_zr, f_pr, f_tr, qx(:, 1:mstep), qy(:, 1:mstep), qz(:, 1:mstep), qp(:, 1:mstep), qt(:, 1:mstep), mstep)

!     --> Normalise the residual vector.
         call normalize(f_xr, f_yr, f_zr, f_pr, f_tr, alpha)

!     --> Update the Hessenberg matrix.
         H(mstep+1, mstep) = alpha

!     --> Add the residual vector as the new Krylov vector.
         call nopcopy(qx(:,mstep+1),qy(:,mstep+1),qz(:,mstep+1),qp(:,mstep+1),qt(:,mstep+1),  f_xr,f_yr,f_zr,f_pr,f_tr)

!     --> Save checkpoint for restarting/run-time analysis.
         if(ifres) call arnoldi_checkpoint(f_xr, f_yr, f_zr, f_pr, f_tr, H(1:mstep+1, 1:mstep), mstep)

!     --> Output timing statistics
         eetime1=dnekclock() ; telapsed = (eetime1-eetime0)/3600.0d0 ; tmiss = telapsed*(ksize-mstep)

         if(nid.eq.0) then
            write(6,"(' Time per iteration/remaining:',I3,'h ',I2,'min /',I3,'h ',I2,'min')"),
     $int(telapsed),ceiling((telapsed-int(telapsed))*60.),
     $int(tmiss),ceiling((tmiss-int(tmiss))*60.)
            print *, ''
         endif

      enddo

      return
      end subroutine arnoldi_factorization





!     -------------------------------------------------------------------------





      subroutine update_hessenberg_matrix(H, f_xr, f_yr, f_zr, f_pr, f_tr, qx, qy, qz, qp, qt, k)

!     This function orthonormalizes the latest Krylov vector f w.r.t. all of the
!     previous ones and updates the entries of the Hessenberg matrix accordingly.
!     
!     INPUTS
!     ------
!     
!     k : int
!     Current step of the Arnoldi factorization.
!     
!     f_xr, f_yr, f_zr : nek arrays of size lt = lx1*ly1*lz1*lelt.
!     Velocity components of the latest Krylov vector.
!     When returned, it has been orthonormalized w.r.t. to all previous
!     Krylov vectors.
!     
!     f_pr : nek array of size lt2 = lx2*ly2*lz2*lelt.
!     Pressure component of the latest Krylov vector.
!     
!     qx, qy, qz : nek arrays of size (lt, k)
!     Velocity components of the Krylov basis.
!     
!     qp : nek array of size (lt2, k).
!     Pressure component of the Krylov basis.
!     
!     H : k x k real matrix.
!     Upper Hessenberg matrix.
!     
!     Last edit : April 3rd 2020 by JC Loiseau.

      implicit none
      include "SIZE"
      include "TOTAL"

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      integer, intent(in) :: k
      real, dimension(k, k), intent(inout) :: H

      real, dimension(lt, k), intent(in) :: qx, qy, qz, qt
      real, dimension(lt2, k), intent(in) :: qp

      real, dimension(lt), intent(inout) :: f_xr, f_yr, f_zr, f_tr
      real, dimension(lt2), intent(inout) :: f_pr

      integer i, n
      real alpha

      real, dimension(lt) :: wk1, wk2, wk3, wkt
      real, dimension(lt2) :: wkp
      real, dimension(k) :: h_vec

      n = nx1*ny1*nz1*nelt

!     --> Initialize array.
      call rzero(h_vec, k)

!     --> Orthonormalize f w.r.t the Krylov basis.
      do i = 1, k

!     --> Copy the i-th Krylov vector to the working arrays.
         call nopcopy(wk1,wk2,wk3,wkp,wkt, qx(:,i),qy(:,i),qz(:,i),qp(:,i),qt(:,i))

!     --> Orthogonalize f w.r.t. to q_i.
         call inner_product(alpha, f_xr, f_yr, f_zr, f_pr, f_tr, wk1, wk2, wk3, wkp, wkt)
         call nopcmult(wk1, wk2, wk3, wkp, wkt, alpha)
         call nopsub2(f_xr,f_yr,f_zr,f_pr,f_tr, wk1,wk2,wk3,wkp,wkt)

!     --> Update the corresponding entry in the Hessenberg matrix.
         H(i, k) = alpha

      enddo

      return
      end subroutine update_hessenberg_matrix
