      
      subroutine arnoldi_factorization(Q, H, mstart, mend, ksize)
      
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
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      !     ----- Miscellaneous -----
         real :: alpha
         integer, intent(in) :: mstart, mend, ksize
      
      !     ----- Timer -----
         real :: eetime0, eetime1, telapsed, tmiss, avg_time
         real, external :: dnekclock
      
      !     ----- Orthogonal residual f = w - (Q,w)*Q -----
         type(krylov_vector) :: f
      !     ----- Krylov basis V for the projection MQ = QH -----
         type(krylov_vector), dimension(ksize + 1), intent(inout) :: Q
      
      !     ----- Upper Hessenberg matrix -----
         real, dimension(ksize + 1, ksize), intent(out) :: H
      
      !     ----- Check k_dim -----
         if (ksize == 0) then
            if (nid == 0) write (6, *) 'Krylov base dimension == 0! Increase it.. STOP'
            call nek_end
         end if
      
      !     --> Initialize arrays.
         call k_zero(f); alpha = 0.0d0
      
      !     --> Arnoldi factorization.
         eetime0 = dnekclock() ! Start time for entire process
         do mstep = mstart, mend
      
            if (mstart < mend) then
               if (nid == 0) write (6, "('        ARNOLDI - Starting iteration ',I3,'/',I3)") mstep, mend
            end if
      
      !     --> Matrix-vector product f = M * v (e.g. calling the linearized Navier-Stokes solver).
            call matvec(f, Q(mstep))
      
      !     --> Update Hessenberg matrix and compute the orthogonal residual f.
            call update_hessenberg_matrix(H(1:mstep + 1, 1:mstep), f, Q(1:mstep), mstep)
      
      !     --> Add the residual vector as the new Krylov vector.
            call k_copy(Q(mstep + 1), f)
      
      !     --> Save checkpoint for restarting/run-time analysis. #  --> not in newton !
            if (ifres .and. (mstart < mend)) call arnoldi_checkpoint(f%vx, f%vy, f%vz, f%pr, f%t, H(1:mstep + 1, 1:mstep), mstep)
      
            if (nid == 0) then
               if (mstart == mend) then
                  write (6, "('        ARNOLDI - Finished iteration ',I3,'/',I3)") mstep, mend
               else
                  eetime1 = dnekclock()
                  telapsed = (eetime1 - eetime0)/3600.0d0
                  avg_time = telapsed/mstep ! Average time per iteration
                  tmiss = avg_time*(mend - mstep) ! More accurate ETA based on average
                  write (6, "('        ARNOLDI - Finished iteration:',I3,'/',I3,' elapsed:',I3,'h',I2,'m / ETA:',I3,'h',I2,'m')")
     $   mstep, mend, int(telapsed), ceiling((telapsed - int(telapsed))*60.0d0),
     $   int(tmiss), ceiling((tmiss - int(tmiss))*60.0d0)
               end if
            end if
         end do
      
         return
      end subroutine arnoldi_factorization
      
      !     -------------------------------------------------------------------------
      
      subroutine update_hessenberg_matrix(H, f, q, k)
      
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
      !     f_pr : nek array of size lp = lx2*ly2*lz2*lelt.
      !     Pressure component of the latest Krylov vector.
      !
      !     qx, qy, qz : nek arrays of size (lv, k)
      !     Velocity components of the Krylov basis.
      !
      !     qp : nek array of size (lp, k).
      !     Pressure component of the Krylov basis.
      !
      !     H : k x k real matrix.
      !     Upper Hessenberg matrix.
      !
      !     Last edit : April 3rd 2020 by JC Loiseau.
      
         use krylov_subspace
         implicit none
         include "SIZE"
         include "TOTAL"
      
         integer, intent(in) :: k
         real, dimension(k + 1, k), intent(inout) :: H
      
         type(krylov_vector), dimension(k), intent(in) :: q
         type(krylov_vector), intent(inout) :: f
         type(krylov_vector) :: wrk
      
         integer :: i
         real :: alpha, beta
      
         call k_norm(beta, f)
      
      ! --> Orthonormalize f w.r.t the Krylov basis.
         do i = 1, k
      
      ! --> Copy the i-th Krylov vector to the working arrays.
            call k_copy(wrk, q(i))
      
      ! --> Orthogonalize f w.r.t. to q_i.
            call k_dot(alpha, f, wrk)
            call k_cmult(wrk, alpha)
            call k_sub2(f, wrk)
      
      ! --> Update the corresponding entry in the Hessenberg matrix.
            H(i, k) = alpha
      
         end do
      
      ! --> Perform full re-orthogonalization (see instability of MGS process).
         do i = 1, k
            call k_copy(wrk, q(i))
            call k_dot(alpha, f, wrk)
            call k_cmult(wrk, alpha)
            call k_sub2(f, wrk)
            H(i, k) = H(i, k) + alpha
      
      !if (nid == 0) then
      !   write (*, *) "ALPHA REORTH :", alpha
      !end if
      
         end do
      
      ! --> Normalise the residual vector.
         call k_normalize(f, alpha)
      
      ! --> Update the Hessenberg matrix.
         H(k + 1, k) = alpha
      
         return
      end subroutine update_hessenberg_matrix
