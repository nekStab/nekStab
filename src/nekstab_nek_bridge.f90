      !-----------------------------------------------------------------------
      !     nekstab_nek_bridge.f90 — Fixed-format bridge to Nek5000 globals
      !
      ! Purpose:
      !   This module is the ONLY nekStab file that includes Nek5000 headers.
      !   All other nekStab files use this module instead of include statements.
      !   Compiled with -ffixed-form because Nek5000 headers are fixed-format.
      !   Any change to shared NEKSTAB common blocks (for example adding
      !   thermal_norm_weight) must rebuild this bridge so the used-module
      !   interface seen by free-form Fortran sources stays in sync.
      !
      ! Public interface:
      !   (bridge module — provides Nek globals to free-form sources via use)
      !
      ! Dependencies:
      !   Nek5000 headers: SIZE, TOTAL, ADJOINT
      !-----------------------------------------------------------------------
      module nekstab_nek_bridge
         implicit none
         include 'SIZE'
         include 'TOTAL'
         include 'ADJOINT'
         integer, parameter :: nekStab_dp = kind(0.0d0)

!  nekStab-owned Fortran I/O units. Keep these in one high
!  range to avoid collisions with Nek5000 and legacy units.

!  Fixed-point / residual logs
         integer, parameter :: NEKSTAB_UNIT_RESIDU      = 700
         integer, parameter :: NEKSTAB_UNIT_DYNTOL      = 701

!  Energy / diagnostics
         integer, parameter :: NEKSTAB_UNIT_ENERGY      = 710
         integer, parameter :: NEKSTAB_UNIT_ENSTRO      = 711
         integer, parameter :: NEKSTAB_UNIT_PKE         = 712
         integer, parameter :: NEKSTAB_UNIT_TORQUE      = 713
         integer, parameter :: NEKSTAB_UNIT_ZC1         = 714
         integer, parameter :: NEKSTAB_UNIT_ZC2         = 715

!  Eigenvalue / Krylov output
         integer, parameter :: NEKSTAB_UNIT_HESS        = 730
         integer, parameter :: NEKSTAB_UNIT_FICH1       = 731
         integer, parameter :: NEKSTAB_UNIT_FICH2       = 732
         integer, parameter :: NEKSTAB_UNIT_FICH3       = 733
         integer, parameter :: NEKSTAB_UNIT_FICH4       = 734
         integer, parameter :: NEKSTAB_UNIT_EXPORT      = 735

!  Modal output
         integer, parameter :: NEKSTAB_UNIT_POD         = 750
         integer, parameter :: NEKSTAB_UNIT_DMD         = 751
         integer, parameter :: NEKSTAB_UNIT_SPOD        = 752

!  Input / auxiliary data
         integer, parameter :: NEKSTAB_UNIT_FST         = 760

!  Solver logs
         integer, parameter :: NEKSTAB_UNIT_NEWTON_LOG  = 770
         integer, parameter :: NEKSTAB_UNIT_GMRES_LOG   = 771
         integer, parameter :: NEKSTAB_UNIT_ARNOLDI_LOG = 772

         public :: nekStab_error, nekStab_log

      contains

!-----------------------------------------------------------------------
! nekStab_error — Uniform loud failure
!  Fixed-form layout for this bridge file.
!-----------------------------------------------------------------------
      subroutine nekStab_error(msg)
         character(len=*), intent(in) :: msg
         if (nid == 0) write (6, *) 'ERROR: ', trim(msg)
         call nek_end
      end subroutine nekStab_error

!-----------------------------------------------------------------------
! nekStab_log — Centralized rank-0 status with short prefix
!-----------------------------------------------------------------------
      subroutine nekStab_log(msg)
         character(len=*), intent(in) :: msg
         if (nid == 0) write (6, *) 'nekStab: ', trim(msg)
      end subroutine nekStab_log

      end module nekstab_nek_bridge
