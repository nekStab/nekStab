!-----------------------------------------------------------------------
! base_mut.f90 — Base eddy viscosity for quasilaminar RANS stability
!
! Purpose:
!   Implements the frozen-eddy-viscosity ("quasi-laminar") linearisation of
!   a RANS base flow. The eddy viscosity mu_t is captured ONCE on the base
!   flow and then held fixed through every Frechet evaluation, so the
!   perturbation operator acts on velocity alone and the turbulence model
!   does not respond to the perturbation (Mettot & Sipp 2014; Meliga et al.
!   2012; Pickering et al. 2021). This is the complement of the coupled
!   operator (the full RANS Jacobian) selected by ifquasilaminar in the
!   .usr.
!
! Why the capture is fed from the .usr:
!   The eddy-viscosity getter (rans_mut) is part of experimental/rans_komg.f,
!   which is include'd into the case .usr — it lives in the case object, not
!   in libnek5000.a. A src/ routine therefore cannot reference it directly
!   (that would leave an unresolved symbol for every non-RANS example). The
!   general capture loop lives here and receives the getter as a procedure
!   argument, so the reusable machinery is in src/ while the only
!   model-specific reference stays in the .usr.
!
! Faithfulness:
!   The snapshot is the literature definition of a frozen eddy viscosity: a
!   fixed spatial field mu_t,base(x). Indexing into it makes the freeze
!   model-agnostic (correct even for a strain-dependent SST model, where
!   solver=none alone would still let mu_t track the perturbed velocity) and
!   immune to in-place limit_ktau clipping of the (now unused) k,tau fields.
!   For the standard k-tau model (m_id=4, mu_t = rho*alp_str*k*tau) mu_t is
!   velocity-independent, so capturing on the first finite-difference substep
!   — where k,tau already equal the base values — recovers the base mu_t
!   exactly; for SST the capture is off by O(epsilon)=O(1e-6), negligible.
!
! Public interface:
!   base_mut_capture(mut_fun)    — snapshot the base eddy viscosity once,
!                                  calling the passed getter mut_fun(ix,iy,iz,iel)
!   base_mut_get(ix,iy,iz,iel)   — fetch the frozen value at a point
!   base_mut_ready()             — .true. once the snapshot has been taken
!   base_mut_reset()             — invalidate the snapshot (force re-capture)
!
! Dependencies:
!   nekstab_nek_bridge (sizes, nelv, nid)
!-----------------------------------------------------------------------
module nekstab_base_mut

   use nekstab_nek_bridge, only: lx1, ly1, lz1, lelv, nelv, nid
   implicit none

   private

   ! Compile-time storage extent — matches the (lx1,ly1,lz1,lelv) layout of
   ! the RANS mut field so the linear index below is stride-consistent.
   integer, parameter :: lmut = lx1*ly1*lz1*lelv

   real, save :: mut0(lmut) = 0.0d0  ! captured base eddy-viscosity field
   logical, save :: mut0_ready = .false.

   public :: base_mut_capture, base_mut_get, base_mut_ready, &
             base_mut_reset

contains

!-----------------------------------------------------------------------
! Snapshot the base eddy viscosity. mut_fun is the case-local getter
! (rans_mut), passed in so this routine carries no RANS dependency.
! The loop walks the full compile-time (lx1,ly1,lz1) extent in
! ix-fastest order so the running index matches base_mut_get below;
! mut_fun's recompute-on-(1,1,1,1) refreshes the field from the current
! (base) k,tau before the first read.
!-----------------------------------------------------------------------
   subroutine base_mut_capture(mut_fun)
      real, external :: mut_fun
      integer :: ix, iy, iz, iel, idx

      idx = 0
      do iel = 1, nelv
         do iz = 1, lz1
            do iy = 1, ly1
               do ix = 1, lx1
                  idx = idx + 1
                  mut0(idx) = mut_fun(ix, iy, iz, iel)
               end do
            end do
         end do
      end do
      mut0_ready = .true.
      if (nid == 0) write (6, *) &
         'base_mut: captured base eddy viscosity (quasilaminar)'
   end subroutine base_mut_capture

!-----------------------------------------------------------------------
! Fetch the frozen eddy viscosity at a point. Index uses compile-time
! lx1,ly1,lz1 strides to match the captured (lx1,ly1,lz1,lelv) layout.
!-----------------------------------------------------------------------
   real function base_mut_get(ix, iy, iz, iel)
      integer, intent(in) :: ix, iy, iz, iel
      integer :: idx

      idx = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(iel - 1)))
      base_mut_get = mut0(idx)
   end function base_mut_get

   logical function base_mut_ready()
      base_mut_ready = mut0_ready
   end function base_mut_ready

   subroutine base_mut_reset()
      mut0_ready = .false.
   end subroutine base_mut_reset

end module nekstab_base_mut
