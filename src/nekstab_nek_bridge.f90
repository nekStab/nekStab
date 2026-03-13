c-----------------------------------------------------------------------
c     nekstab_nek_bridge.f90 -- Fixed-format bridge to Nek5000 globals
c
c     This module is the ONLY nekStab file that includes Nek5000 headers.
c     All other nekStab files use this module instead of include statements.
c     Compiled with -ffixed-form because Nek5000 headers are fixed-format.
c-----------------------------------------------------------------------
       module nekstab_nek_bridge
          implicit none
           include 'SIZE'
           include 'TOTAL'
           include 'ADJOINT'
          integer, parameter :: nekStab_dp = kind(0.0d0)
       end module
