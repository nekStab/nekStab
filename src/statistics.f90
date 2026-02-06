      !-----------------------------------------------------------------------
      ! statistics.f90 — Time-averaging and RMS statistics
      !
      ! Purpose:
      !   Computes running mean and RMS statistics of velocity, pressure,
      !   and temperature fields during DNS. Periodically outputs averaged
      !   fields and fluctuation fields to disk.
      !
      ! Public interface:
      !   nekStab_avg — compute and output running statistics
      !
      ! Dependencies:
      !   krylov_subspace, SIZE, TOTAL, AVG
      !-----------------------------------------------------------------------

      module nekstab_statistics
         use krylov_subspace
         private
         public :: nekStab_avg
      contains

      !-----------------------------------------------------------------------
      ! nekStab_avg — Running mean and RMS statistics
      !
      ! Purpose:
      !   Incrementally computes time-averaged fields and (optionally)
      !   second-order statistics (RMS). Outputs at user-defined intervals.
      !
      ! Arguments:
      !   ifstatis [in] — if .true., also compute fluctuation statistics
      !-----------------------------------------------------------------------
      subroutine nekStab_avg(ifstatis)
         use krylov_subspace
         include 'SIZE'
         include 'TOTAL'
         include 'AVG'

         logical, intent(in) :: ifstatis
         real do1(lv), do2(lv), do3(lv)

         real, save :: x0(3)
         data x0/0.0, 0.0, 0.0/

         integer, save :: icalld
         data icalld/0/

         logical ifverbose

         if (ax1 /= lx1 .or. ay1 /= ly1 .or. az1 /= lz1) then
            if (nid == 0) write (6, *) 'ABORT: wrong size of ax1,ay1,az1 in avg_all(), check SIZE!'
            call exitt
         end if
         if (ax2 /= lx2 .or. ay2 /= ay2 .or. az2 /= lz2) then
            if (nid == 0) write (6, *) 'ABORT: wrong size of ax2,ay2,az2 in avg_all(), check SIZE!'
            call exitt
         end if

         ntot = lx1*ly1*lz1*nelv
         nto2 = lx2*ly2*lz2*nelv

      !     initialization
         if (icalld == 0) then
            icalld = icalld + 1
            atime = 0.; timel = time
            call oprzero(uavg, vavg, wavg)
            call rzero(pavg, nto2)
            call oprzero(urms, vrms, wrms)
            call rzero(prms, nto2)
            call oprzero(vwms, wums, uvms)
         end if

         dtime = time - timel
         atime = atime + dtime

      !     dump freq
         iastep = param(68)
         if (iastep == 0) iastep = param(15) ! same as iostep
         if (iastep == 0) iastep = 500

         ifverbose = .false.
         if (istep <= 10) ifverbose = .true.
         if (mod(istep, iastep) == 0) ifverbose = .true.

         if (atime /= 0. .and. dtime /= 0.) then
            if (nio == 0) write (6, *) 'Computing statistics ...'
            beta = dtime/atime
            alpha = 1.-beta

      !     compute averages E(X) !avg(k) = alpha*avg(k) + beta*f(k)
            call avg1(uavg, vx, alpha, beta, ntot, 'um  ', ifverbose)
            call avg1(vavg, vy, alpha, beta, ntot, 'vm  ', ifverbose)
            call avg1(wavg, vz, alpha, beta, ntot, 'wm  ', ifverbose)
            call avg1(pavg, pr, alpha, beta, nto2, 'prm ', ifverbose)
            call avg1(tavg, t(1, 1, 1, 1, 1), alpha, beta, ntot, 'tm ', ifverbose)

            if (ifstatis) then !compute fluctuations
      !     compute averages E(X^2)
               call avg2(urms, vx, alpha, beta, ntot, 'ums ', ifverbose)
               call avg2(vrms, vy, alpha, beta, ntot, 'vms ', ifverbose)
               call avg2(wrms, vz, alpha, beta, ntot, 'wms ', ifverbose)
               call avg2(prms, pr, alpha, beta, nto2, 'prms', ifverbose)
               call avg2(trms, t(1, 1, 1, 1, 1), alpha, beta, ntot, 'trms', ifverbose)

      !     compute averages E(X*Y)

            end if


         end if

         if ((mod(istep, iastep) == 0 .and. istep > 1) .or. lastep == 1) then

            time_temp = time
            time = atime ! Output the duration of this avg
            dtmp = param(63)
            param(63) = 1 ! Enforce 64-bit output

      !     mean fluctuation fields
            ifto = .false.; ifpo = .false.
            call opsub3(do1, do2, do3, uavg, vavg, wavg, ubase, vbase, wbase)
            call outpost(do1, do2, do3, pr, t, 'avt')

            call outpost2(uavg, vavg, wavg, pavg, tavg, ldimt, 'avg')
            call outpost2(urms, vrms, wrms, prms, trms, ldimt, 'rms')
            call outpost(uvms, vwms, wums, prms, trms, 'rm2')


            param(63) = dtmp
            atime = 0.
            time = time_temp ! Restore clock

         end if
         timel = time

         return
      end subroutine nekStab_avg
      !-----------------------------------------------------------------------

      end module nekstab_statistics
