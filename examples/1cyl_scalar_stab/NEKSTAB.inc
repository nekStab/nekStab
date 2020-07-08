      ! LINEAR STABILITY FRAMEWORK
      ! Shared memory bounded arrays and parameters

      logical, parameter :: ifnois = .true.
      logical, parameter :: ifldbf = .true.

      character(len=80), parameter :: bf_handle = 'BF_1cyl0.f00001'

      integer, parameter :: k_dim  = 85
      integer, parameter :: schur_tgt  = 1 !schur on if >1

      real, parameter    :: eigen_tol  = 1.0E-6 !standard value
      real, parameter    :: schur_del  = 0.10D0
      logical              ifbfcv
      common /lst_logical/ ifbfcv

      real ubase(lx1,ly1,lz1,lelt)
     $    ,vbase(lx1,ly1,lz1,lelt)
     $    ,wbase(lx1,ly1,lz1,lelt)
     $    ,tbase(lx1,ly1,lz1,lelt)
      common /baseflow/ ubase,vbase,wbase,tbase

      real fcx(lx1,ly1,lz1,lelt)
     $    ,fcy(lx1,ly1,lz1,lelt)
     $    ,fcz(lx1,ly1,lz1,lelt)
      common /force/ fcx, fcy, fcz !general forcing array

      real spng_fun(LX1*LY1*LZ1*LELV)     ! sponge function
     $    ,spng_vr(LX1*LY1*LZ1*LELV,LDIM) ! reference velocity field
      common /SPONGEV/ spng_fun, spng_vr

      real spng_str      !var sponge strength
      real spng_wl(LDIM) !var sponge width (left section; every dmension separately)
      real spng_wr(LDIM) !var sponge width (right section)
      real spng_dl(LDIM) !var sponge drop/rise width (left section)
      real spng_dr(LDIM) !var sponge drop/rise width (right section)
      common /SPONGER/ spng_str, spng_wl, spng_wr, spng_dl, spng_dr

      integer, parameter :: stat_lvar = 44 ! total number of statistical fields
      real stat(lx1*ly1*lz1*lelt,stat_lvar)
     $    ,stat_temp(lx1*ly1*lz1*lelt)
     $    ,stat_uu(lx1*ly1*lz1*lelt)
     $    ,stat_vv(lx1*ly1*lz1*lelt)
     $    ,stat_ww(lx1*ly1*lz1*lelt)
     $    ,stat_pp(lx1*ly1*lz1*lelt)
     $    ,stat_uuu(lx1*ly1*lz1*lelt)
     $    ,stat_vvv(lx1*ly1*lz1*lelt)
     $    ,stat_www(lx1*ly1*lz1*lelt)
     $    ,stat_ppp(lx1*ly1*lz1*lelt)
     $    ,stat_atime,stat_tstart,stat_comp

      common /statsr/ stat,stat_temp,stat_uu,stat_vv,stat_ww,stat_pp
     $  ,stat_uuu,stat_vvv,stat_www,stat_ppp,stat_atime,stat_tstart,stat_comp
