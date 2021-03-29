c-----------------------------------------------------------------------
      subroutine whereyouwant(resnam,posfil) !file numbering suffix counter
      character*3 resnam
      integer posfil
      common /RES_WANT/ nopen(99,2)

!     change prepost.f line 1094 from "save nopen" to "common /RES_WANT/ nopen"

      iprefix          = i_find_prefix(resnam,99)
      nopen(iprefix,1) = posfil -1

      return
      end
c-----------------------------------------------------------------------
      subroutine load_files(V_x,V_y,V_z,P_r,V_t,mstart,kd,fname)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelv
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelv
      integer                            :: kd,mstart,n,n2,i,j

c     ----- Krylov basis V for the projection M*V = V*H -----

      real, dimension(lt,kd+1)        :: V_x, V_y, V_z, V_t
      real, dimension(lt2,kd+1)       :: P_r
      character*3  fname
      character*60  filename
      character(len=7)   :: tl
      character(len=20)  :: fmt

      n     = nx1*ny1*nz1*nelv
      n2    = nx2*ny2*nz2*nelv

c     ----- Upload the snapshots -----

      do i = 1, mstart

         if(i.lt.10)                        fmt = '(i1.1)' ! 1 to 9
         if(i.ge.10.and.i.lt.100)           fmt = '(i2.2)' ! 10 to 99
         if(i.ge.100.and.i.lt.1000)         fmt = '(i3.3)' ! 100 to 999
         if(i.ge.1000.and.i.lt.10000)       fmt = '(i4.4)' ! 1000 to 9999
         if(i.ge.10000.and.i.lt.100000)     fmt = '(i5.5)' ! 10000 to 99999

         j=i
         write(tl,fmt) j
         if(i.lt.10)
     $        filename = trim(fname)//trim(SESSION)//'0.f0000'//trim(tl)
         if(i.ge.10.and.i.lt.100)
     $        filename = trim(fname)//trim(SESSION)//'0.f000'//trim(tl)
         if(i.ge.100.and.i.lt.1000)
     $        filename = trim(fname)//trim(SESSION)//'0.f00'//trim(tl)
         if(i.ge.1000.and.i.lt.10000)
     $        filename = trim(fname)//trim(SESSION)//'0.f0'//trim(tl)
         if(i.ge.10000.and.i.lt.100000)
     $        filename = trim(fname)//trim(SESSION)//'0.f'//trim(tl)

         call load_fld(filename)
         call opcopy(V_x(:,i),V_y(:,i),V_z(:,i),vx,vy,vz)
         if(ifpo)call copy(P_r(:,i),pr,n2)
         if(ifto)call copy(V_t(:,i),t(1,1,1,1,1),n)

      enddo
      return
      end
c-----------------------------------------------------------------------
