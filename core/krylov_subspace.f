      module krylov_subspace
      implicit none
      include 'SIZE'

      private

      integer, public, parameter :: lv = lx1*ly1*lz1*lelv
      integer, public, parameter :: lp = lx2*ly2*lz2*lelv
      integer, save, public :: n,n2
      type, public :: krylov_vector
      real, dimension(lv) :: vx, vy, vz
      real, dimension(lp) :: pr
      real, dimension(lv,ldimt) :: theta
      real :: time
      end type krylov_vector

      type(krylov_vector), save, public :: ic_nwt, fc_nwt
      real,save,allocatable,dimension(:, :), public :: uor,vor,wor,tor

      contains
      end module krylov_subspace

      subroutine krylov_outpost(q, prefix)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      character(len=3), intent(in) :: prefix
      type(krylov_vector), intent(in) :: q

      call outpost(q%vx, q%vy, q%vz, q%pr, q%theta, prefix)

      return
      end subroutine krylov_outpost

      subroutine krylov_inner_product(alpha, p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector), intent(in) :: p, q
      real, intent(out) :: alpha
      real :: glsc3
      integer m

      n = nx1*ny1*nz1*nelv

!     --> Kinetic energy.
      alpha = glsc3(p%vx, bm1s, q%vx, n) + glsc3(p%vy, bm1s, q%vy, n)
      if (if3d) alpha = alpha + glsc3(p%vz, bm1s, q%vz, n)

!     --> Potential energy.
      if (ldimt.gt.0) then
       do m = 1,ldimt
        alpha = alpha + glsc3(p%theta(:,m), bm1s, q%theta(:,m), n)
       enddo
      endif

!     --> Time component.
      if ( uparam(1) .eq. 2.1 ) then
         alpha = alpha + p%time * q%time
      end if

!     --> Check integrity.
      if ( isnan(alpha) ) call nek_end

      return
      end subroutine krylov_inner_product

      subroutine krylov_norm(alpha, p)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector), intent(in) :: p
      real, intent(out) :: alpha

      call krylov_inner_product(alpha, p, p)
      alpha = dsqrt(alpha)
      return
      end subroutine krylov_norm

      subroutine krylov_normalize(p, alpha)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector),intent(inout) :: p
      real, intent(out) :: alpha
      real :: inv_alpha

!     --> Compute the user-defined norm.
      call krylov_norm(alpha, p)
      inv_alpha = 1.0D+00 / alpha

!     --> Normalize the vector.
      call krylov_cmult(p, inv_alpha)

      return
      end subroutine krylov_normalize

      subroutine krylov_cmult(p, alpha)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p
      real alpha
      integer i

      n = nx1*ny1*nz1*nelv
      n2 = nx2*ny2*nz2*nelv

                call cmult(p%vx(:),alpha,n)
                call cmult(p%vy(:),alpha,n)
                call cmult(p%pr(:),alpha,n2)
      if (if3d) call cmult(p%vz(:),alpha,n)
      if (ldimt.gt.0) then
        do i = 1,ldimt
          call cmult(p%theta(:,i),alpha,n)
        enddo
      endif
      p%time = p%time * alpha

      return
      end subroutine krylov_cmult

      subroutine krylov_add2(p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q
      integer i

      n = nx1*ny1*nz1*nelv
      n2 = nx2*ny2*nz2*nelv

                call add2(p%vx(:),q%vx(:),n)
                call add2(p%vy(:),q%vy(:),n)
                call add2(p%pr(:),q%pr(:),n2)
      if (if3d) call add2(p%vz(:),q%vz(:),n)
      if (ldimt.gt.0) then
        do i = 1,ldimt
                call add2(p%theta(:,i),q%theta(:,i),n)
        enddo
      endif
      p%time = p%time + q%time

      return
      end subroutine krylov_add2


      subroutine krylov_sub2(p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q
      integer i
      n = nx1*ny1*nz1*nelv
      n2 = nx2*ny2*nz2*nelv

                call sub2(p%vx(:),q%vx(:),n)
                call sub2(p%vy(:),q%vy(:),n)
                call sub2(p%pr(:),q%pr(:),n2)
      if (if3d) call sub2(p%vz(:),q%vz(:),n)
      if (ldimt.gt.0) then
       do i = 1,ldimt
                call sub2(p%theta(:,i),q%theta(:,i),n)
       enddo
      endif
      p%time = p%time - q%time

      return
      end subroutine krylov_sub2

      subroutine krylov_zero(p)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p
      integer i
      n = nx1*ny1*nz1*nelv
      n2 = nx2*ny2*nz2*nelv

                call rzero(p%vx(:),n)
                call rzero(p%vy(:),n)
                call rzero(p%pr(:),n2)
      if (if3d) call rzero(p%vz(:),n)
      if (ldimt.gt.0) then
         do i = 1,ldimt
                call rzero(p%theta(:,i),n)
         enddo
      endif
      p%time = 0.0D+00

      return
      end subroutine krylov_zero

      subroutine krylov_copy(p, q)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'
      type(krylov_vector) :: p, q
      integer i
      n = nx1*ny1*nz1*nelv
      n2 = nx2*ny2*nz2*nelv

                call copy(p%vx(:),q%vx(:),n)
                call copy(p%vy(:),q%vy(:),n)
                call copy(p%pr(:),q%pr(:),n2)
      if (if3d) call copy(p%vz(:),q%vz(:),n)
      if (ldimt.gt.0) then
        do i = 1,ldimt
                call copy(p%theta(:,i),q%theta(:,i),n)
        enddo
      endif
      p%time = q%time

      return
      end subroutine krylov_copy

      subroutine krylov_matmul(dq, Q, yvec, k)
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer :: i, j, k
      type(krylov_vector) :: dq
      type(krylov_vector), dimension(k) :: Q
      real, dimension(k) :: yvec

      real, dimension(lv, k) :: qx, qy, qz
      real, dimension(lp, k) :: qp
      real, dimension(lv, k, ldimt) :: qt

      real, dimension(k) :: time_comp

      do i = 1, k
                        qx(:, i) = Q(i)%vx(:)
                        qy(:, i) = Q(i)%vy(:)
                        qp(:, i) = Q(i)%pr(:)
         if (if3d)      qz(:, i) = Q(i)%vz(:)
         if (ldimt.gt.0) then
            do j = 1,ldimt
                        qt(:, i, j) = Q(i)%theta(:,j)
            enddo
         endif
         time_comp(i) = Q(i)%time
      enddo

      call krylov_zero(dq)

                        dq%vx(:) = matmul(qx(:,:), yvec(:))
                        dq%vy(:) = matmul(qy(:,:), yvec(:))
                        dq%pr(:) = matmul(qp(:,:), yvec(:))
      if(if3d)          dq%vz(:) = matmul(qz(:,:), yvec(:))
      if(ldimt.gt.0)then
       do j = 1,ldimt
                dq%theta(:,j) = matmul(qt(:,:,j), yvec(:))
       enddo
      endif
      dq%time = dot_product(time_comp(:), yvec(:))

      return
      end subroutine krylov_matmul

      subroutine krylov_biorthogonalize(real_p, imag_p, real_q, imag_q)
        use krylov_subspace
        implicit none
        include "SIZE"
        include "TOTAL"

        type(krylov_vector), intent(inout) :: real_p, imag_p
        type(krylov_vector), intent(inout) :: real_q, imag_q
        type(krylov_vector) :: wrk1, wrk2, wrk3, wrk4
        real :: alpha, beta, gamma, delta

        ! --> Normalize the direct mode to unit-norm.
        call krylov_norm(alpha, real_p) ; alpha = alpha**2
        call krylov_norm(beta, imag_p)  ; beta  = beta**2

        delta = 1.0D+00 / sqrt(alpha + beta)
        call krylov_cmult(real_p, delta)
        call krylov_cmult(imag_p, delta)

        ! --> Inner product between direct and adjoint modes.
        call krylov_inner_product(alpha, real_p, real_q)
        call krylov_inner_product(beta, imag_p, imag_q)
        gamma = alpha + beta

        call krylov_inner_product(alpha, real_q, imag_p)
        call krylov_inner_product(beta, imag_q, real_p)
        delta = alpha - beta

        ! --> Bi-orthogonalize the adjoint mode.
        call krylov_copy(wrk1, real_q) ; call krylov_copy(wrk2, imag_q)
        call krylov_cmult(wrk1, gamma) ; call krylov_cmult(wrk2, delta)
        call krylov_sub2(wrk1, wrk2)   ; call krylov_copy(wrk3, wrk1)
        call krylov_cmult(wrk3, 1.D+00 / (gamma**2 + delta**2))

        call krylov_copy(wrk1, real_q) ; call krylov_copy(wrk2, imag_q)
        call krylov_cmult(wrk1, delta) ; call krylov_cmult(wrk2, gamma)
        call krylov_add2(wrk1, wrk2)   ; call krylov_copy(wrk4, wrk1)
        call krylov_cmult(wrk4, 1.D+00 / (gamma**2 + delta**2))

        call krylov_copy(real_q, wrk3) ; call krylov_copy(imag_q, wrk4)

        return
      end subroutine krylov_biorthogonalize

      subroutine krylov_gradient(dxp, dyp, dzp, p)
        use krylov_subspace
        implicit none
        include "SIZE"
        include "TOTAL"

        type(krylov_vector), intent(in) :: p
        type(krylov_vector), intent(out) :: dxp, dyp, dzp

        call gradm1(dxp%vx, dyp%vx, dzp%vx, p%vx, nelv)
        call dsavg(dxp%vx) ; call dsavg(dyp%vx) ; call dsavg(dzp%vx)

        call gradm1(dxp%vy, dyp%vy, dzp%vy, p%vy, nelv)
        call dsavg(dxp%vy) ; call dsavg(dyp%vy) ; call dsavg(dzp%vy)

        call gradm1(dxp%vz, dyp%vz, dzp%vz, p%vz, nelv)
        call dsavg(dxp%vz) ; call dsavg(dyp%vz) ; call dsavg(dzp%vz)

        return
      end subroutine krylov_gradient
