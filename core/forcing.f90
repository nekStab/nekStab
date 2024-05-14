      !--------------------------------------------------------------------
      subroutine nekStab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real ffx, ffy, ffz
         integer ix, iy, iz, ieg
         integer iel, ip
      
      ! local element number
         iel = gllel(ieg)
      
      ! here we add other forcings (SFD, TDF, if any)
         ffx = ffx + fcx(ix, iy, iz, iel)
         ffy = ffy + fcy(ix, iy, iz, iel)
         if (if3D) ffz = ffz + fcz(ix, iy, iz, iel)
      
         if (spng_st /= 0) then
      
      ! compute the corresponding index in the pertubation arrays
            ip = ix + nx1*(iy - 1 + ny1*(iz - 1 + nz1*(iel - 1)))
      
            if (jp == 0) then ! dns
               ffx = ffx + spng_fn(ip)*(spng_vr(ip, 1) - vx(ix, iy, iz, iel))*spng_st
               ffy = ffy + spng_fn(ip)*(spng_vr(ip, 2) - vy(ix, iy, iz, iel))*spng_st
               if (if3D) ffz = ffz + spng_fn(ip)*(spng_vr(ip, ndim) - vz(ix, iy, iz, iel))*spng_st
            else ! perturbation
               ffx = ffx - spng_fn(ip)*vxp(ip, jp) ! spng_st alaways = 1
               ffy = ffy - spng_fn(ip)*vyp(ip, jp) ! spng_st alaways = 1
               if (if3D) ffz = ffz - spng_fn(ip)*vzp(ip, jp) ! spng_st alaways = 1
            end if
      
         end if
         return
      end subroutine nekStab_forcing
      !--------------------------------------------------------------------
      subroutine nekStab_forcing_temp(temp, ix, iy, iz, ieg, m)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real temp
         integer ix, iy, iz, ieg
         integer iel, m, ip
      
      ! local element number
         iel = gllel(ieg)
         if (jp == 0) temp = temp + fct(ix, iy, iz, iel, m)
      
         if (spng_st /= 0) then   !!!HERE SPONGE STRENGHT ALWAYS UNITY!
      
      ! compute the corresponding index in the pertubation arrays
            ip = ix + nx1*(iy - 1 + ny1*(iz - 1 + nz1*(iel - 1)))
      
            if (jp == 0) then      ! dns ! t(1,1,1,1,ifield-1)
               temp = temp + spng_fn(ip)*(spng_vt(ip, m) - t(ix, iy, iz, iel, m))
            else                   ! perturbation   ! tp(lpx1*lpy1*lpz1*lpelt,ldimt,lpert)
               temp = temp - spng_fn(ip)*tp(ip, m, jp)
            end if
      
         end if
         return
      end subroutine nekStab_forcing_temp
      
      !--------------------------------------------------------------------
      subroutine activate_sponge
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer i
         nv = nx1*ny1*nz1*nelv
         if (spng_st /= 0) then !sponge on
      
            if (nid == 0) write (6, *)
            if (nid == 0) write (6, *) ' Initializing sponge...'
            if (nid == 0) write (6, *) ' Sponge strenght:', spng_st
            if (spng_st < 0) then
               spng_st = abs(spng_st)
               if (nid == 0) write (6, *) ' Ensure positive sponge strenght:', spng_st
            end if
            call spng_init
      
      !     applying sponge function to BM1 matrix to remove the sponge zone from eigensolver
            do i = 1, nv
               if (spng_fn(i) /= 0) bm1s(i, 1, 1, 1) = 0.0d0
            end do
      
      !  outposting BM1s to disk for check
      !  ifto_sav = ifto; ifpo_sav = ifpo
      !  ifvo=.false.; ifpo = .false.; ifto = .true.
      !  call outpost(vx,vy,vz,pr,bm1s,'BMS')
      !  ifvo=.true.; ifpo = ifpo_sav; ifto = ifto_sav
      
            if (nid == 0) write (6, *) 'Sponge activated.'
            if (nid == 0) write (6, *)
         end if
      end subroutine activate_sponge
      !--------------------------------------------------------------------
      subroutine spng_init
      ! credits to KTH Toolbox https://github.com/KTH-Nek5000/KTH_Toolbox/blob/b7dc43a92bb6759132a1baae9d290727de29c257/utility/forcing/sponge_box/spongebx.f
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer n
         n = nx1*ny1*nz1*nelv
         acc_spg = abs(acc_spg)
      
         spng_wl(1) = (1.0d0 - acc_spg)*xLspg ! Sponge left section width; dimension X
         spng_wl(2) = (1.0d0 - acc_spg)*yLspg
         if (if3D) spng_wl(3) = (1.0d0 - acc_spg)*zLspg
      
         spng_wr(1) = (1.0d0 - acc_spg)*xRspg ! Sponge right section width; dimension X
         spng_wr(2) = (1.0d0 - acc_spg)*yRspg
         if (if3D) spng_wr(3) = (1.0d0 - acc_spg)*zRspg
      
         spng_dl(1) = (acc_spg)*xLspg ! Sponge left drop/rise section width; dimension X
         spng_dl(2) = (acc_spg)*yLspg
         if (if3D) spng_dl(3) = (acc_spg)*zLspg
      
         spng_dr(1) = (acc_spg)*xRspg !Sponge right drop/rise section width; dimension X
         spng_dr(2) = (acc_spg)*yRspg
         if (if3D) spng_dr(3) = (acc_spg)*zRspg
      
         if (nid == 0) then
            write (6, *) '  Left  section width x y z:', spng_wl
            write (6, *) '  Right section width x y z:', spng_wr
            write (6, *) '  Left  drop/rise width x y z:', spng_dl
            write (6, *) '  Right drop/rise width x y z:', spng_dr
         end if
      
      !     save reference field -> sponge value reference
         call opcopy(spng_vr(1, 1), spng_vr(1, 2), spng_vr(1, NDIM), vx, vy, vz) !only DNS
         if (ifto) call copy(spng_vt(:, 1), t(1, 1, 1, 1, 1), n) !only DNS - temperature
         call spng_set             ! -> compute spng_fn
      
         return
      end subroutine spng_init
      !--------------------------------------------------------------------
      subroutine spng_set
      ! credits to KTH Toolbox https://github.com/KTH-Nek5000/KTH_Toolbox/blob/b7dc43a92bb6759132a1baae9d290727de29c257/utility/forcing/sponge_box/spongebx.f
      !     set sponge function and refernece fields
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         real lcoord(LX1*LY1*LZ1*LELV)
         common/SCRUZ/lcoord
      
         integer ntot, il, jl
         real rtmp, bmin(LDIM), bmax(LDIM)
         real xxmax, xxmax_c, xxmin, xxmin_c, arg, mth_stepf
         logical ltmp, ltmp2
      
         ntot = NX1*NY1*NZ1*NELV
         bmin(1) = xmn
         bmax(1) = xmx
         bmin(2) = ymn
         bmax(2) = ymx
         if (if3D) then
            bmin(NDIM) = zmn
            bmax(NDIM) = zmx
         end if
      
         call rzero(spng_fn, ntot)
      !     for every dimension
         do il = 1, NDIM
      
            if (spng_wl(il) > 0.0 .or. spng_wr(il) > 0.0) then
      
               if (spng_wl(il) < spng_dl(il) .or. spng_wr(il) < spng_dr(il)) then
                  write (6, *) 'Wrong sponge parameters!'
               end if
      
               xxmax = bmax(il) - spng_wr(il) ! sponge beginning (rise at xmax; right)
               xxmin = bmin(il) + spng_wl(il) ! end (drop at xmin; left)
               xxmax_c = xxmax + spng_dr(il) ! beginnign of constant part (right)
               xxmin_c = xxmin - spng_dl(il) ! beginnign of constant part (left)
      
      !     get SPNG_FUN
               if (xxmax <= xxmin) then
                  write (6, *) 'Sponge too wide'
               else
      !     this should be done by pointers, but for now I avoid it
                  if (il == 1) then
                     call copy(lcoord, XM1, ntot)
                  elseif (il == 2) then
                     call copy(lcoord, YM1, ntot)
                  elseif (il == 3) then
                     call copy(lcoord, ZM1, ntot)
                  end if
      
                  do jl = 1, ntot
                     rtmp = lcoord(jl)
                     if (rtmp <= xxmin_c) then ! constant; xmin
                        rtmp = 1.0d0
                     elseif (rtmp < xxmin) then ! fall; xmin
                        arg = (xxmin - rtmp)/spng_wl(il)
                        rtmp = mth_stepf(arg)
                     elseif (rtmp <= xxmax) then ! zero
                        rtmp = 0.0
                     elseif (rtmp < xxmax_c) then ! rise
                        arg = (rtmp - xxmax)/spng_wr(il)
                        rtmp = mth_stepf(arg)
                     else             ! constant
                        rtmp = 1.0d0
                     end if
                     spng_fn(jl) = max(spng_fn(jl), rtmp)
                  end do
               end if                  ! xxmax.le.xxmin
            end if                     ! spng_w(il).gt.0.0
         end do
      
         ltmp = ifto; ltmp2 = ifpo
         ifto = .true.; ifpo = .false.
         call outpost2(spng_vr(1, 1), spng_vr(1, 2), spng_vr(1, NDIM), spng_fn, spng_fn, 1, 'SPG')
         ifto = ltmp; ifpo = ltmp2
      
         return
      end subroutine spng_set
      !--------------------------------------------------------------------
      real function mth_stepf(x) ! credits to KTH Toolbox
      !     compute sponge function
         implicit none
         real x, xdmin, xdmax
         parameter(xdmin=0.0010d0, xdmax=0.9990d0)
         if (x <= xdmin) then
            mth_stepf = 0.0d0
         elseif (x <= xdmax) then
            mth_stepf = 1.0d0/(1.0d0 + exp(1.0d0/(x - 1.0d0) + 1.0d0/x))
         else
            mth_stepf = 1.0d0
         end if
      end function mth_stepf
