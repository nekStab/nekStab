c-----------------------------------------------------------------------
      subroutine fourier_reconstruction
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer, parameter  :: lt = lx1*ly1*lz1*lelt
      real, dimension(lt) :: d_x, d_y, d_z
      integer i

      call oprzero(vx,vy,vz)
      call oprzero(d_x,d_y,d_z)
      do i = 1,f_modes
       call add3s2(d_x,ubr(:,i),ubi(:,i),amplr(i)*cos(om(i)*time),ampli(i)*sin(OM(i)*time),lt)
       call add3s2(d_y,vbr(:,i),vbi(:,i),amplr(i)*cos(om(i)*time),ampli(i)*sin(OM(i)*time),lt)
       if(if3D)
     $ call add3s2(d_z,wbr(:,i),wbi(:,i),amplr(i)*cos(om(i)*time),ampli(i)*sin(OM(i)*time),lt)
       call opadd2(vx, vy, vz, d_x, d_y, d_z)
      enddo
      !call outpost(vx,vy,vz,pr,t,'MFR')!debug only
      return
      end
c-----------------------------------------------------------------------
      subroutine fourier_decomposition (vx_in, vy_in, vz_in, time_in)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt=lx1*ly1*lz1*lelt
      integer n,i,jj

      real,    dimension(lt,norbit)  :: vx_in, vy_in, vz_in
      real,    dimension(   norbit)  :: time_in

      real,    dimension(lt,norbit/2)  :: v_xr, v_yr, v_zr
      real,    dimension(lt,norbit/2)  :: v_xi, v_yi, v_zi
      real,    dimension(   norbit/2)  :: ff, ampl, amp_real, amp_img
      integer, dimension(   norbit/2)  :: ord
      real,    dimension(lt,mxfourier)      :: rex, rey, rez
      real,    dimension(lt,mxfourier)      :: imx, imy, imz
      real                             :: energy, glsum, glsc3

      if(nid.eq.0)write(6,*)' starting fourier decomposition...'
      call opfft3(norbit,vx_in, vy_in, vz_in,v_xr,v_xi,v_yr,v_yi,v_zr,v_zi,time_in,ff)

      energy = 0.0d0
      do i = 1,norbit/2
        amp_real(i) = + glsc3(v_xr(:,i),bm1,v_xr(:,i),lt) + glsc3(v_yr(:,i),bm1,v_yr(:,i),lt)
        if (if3d) amp_real(i) = amp_real(i) + glsc3(v_zi(:,i),bm1,v_zi(:,i),lt)
        amp_img(i) = + glsc3(v_xi(:,i),bm1,v_xi(:,i),lt) + glsc3(v_yi(:,i),bm1,v_yi(:,i),lt)
        if (if3d) amp_img(i) = amp_img(i) + glsc3(v_zi(:,i),bm1,v_zi(:,i),lt)
        amp_real(i) = sqrt(amp_real(i))
        amp_img(i) = sqrt(amp_img(i))
        ampl(i) = sqrt(amp_real(i)**2+amp_img(i)**2)
        energy  = energy + ampl(i)
      enddo

      call sortvec(ampl,ord,norbit/2)!sortvec[a(n),ord(n),n
      if(nid.eq.0)write(6,*)' sorted modes...',energy
      do i = 1,norbit/2
        if(nid.eq.0)write(6,*)'mode: ',i,ampl(i),ord(i)
      enddo

      jj = 2
      do while (sum(ampl(2:jj))/(energy-ampl(1))*100.lt.99.0)
        jj = jj+1
      enddo
      open(unit=25,file='fft_ampl.dat')
      write(25,*) porbit, sum(ampl(2:min(mxfourier,jj))) /(energy-ampl(1))*100, min(mxfourier,jj)

      do i = 1,min(mxfourier,jj)
        call opcopy(rex(:,i),rey(:,i),rez(:,i), v_xr(:,ord(i)),v_yr(:,ord(i)),v_zr(:,ord(i)))
        call opcopy(imx(:,i),imy(:,i),imz(:,i), v_xi(:,ord(i)),v_yi(:,ord(i)),v_zi(:,ord(i)))
        write(25,*) amp_real(i), amp_img(i), ff(ord(i))

      enddo
      close(25)

      do i = 1,min(mxfourier,jj)
        time=i*dt!just to open in paraview!
        ifpo=.false.;ifto=.false.
        call opcmult(rex(:,i),rey(:,i),rez(:,i),1.0d0/amp_real(i))!normalize
        call outpost(rex(:,i),rey(:,i),rez(:,i),pr,t,'fRe')!outpost
        call opcmult(imx(:,i),imy(:,i),imz(:,i),1.0d0/amp_img(i))!normalize
        call outpost(imx(:,i),imy(:,i),imz(:,i),pr,t,'fIm')!outpost
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine sortvec (a,ord,n)
      implicit none
      integer n, i, j, ord(n)
      real a(n), a1(n), temp(2)
      a1(1) = a(1)
      ord(1) = 1
      do i = 2,n
        temp = 0.0d0
        do j = 2,n
          if (a(j).gt.temp(1)) then
            temp(1) = a(j)
            temp(2) = j
          endif
        enddo
        a1(i)  = temp(1)
        ord(i) = int(temp(2))
        a(int(temp(2))) = -1
      enddo
      a = a1
      return
      end
c----------------------------------------------------------------------
c     FFT discretization in time
c     INPUT: a1,a2,a3 =  velocity field
c     sz     = number of snapshots
c     Time     = temporal discretization
c     OUTPUT : b1,b2  = real and imag part of fft(a1)
c     c1,c2  = real and imag part of fft(a2)
c     d1,d2  = real and imag part of fft(a3)
c     OM     = pulsation discretization
      subroutine opfft3 (sz,a1,a2,a3,b1,b2,c1,c2,d1,d2,time,ff)
      implicit none
      include 'SIZE'
      integer, parameter :: lt=lx1*ly1*lz1*lelt
      integer            :: sz
      real, dimension(lt,sz)   :: a1,a2,a3
      real, dimension(lt,sz/2) :: b1,c1,d1,b2,c2,d2
      real, dimension(   sz)      :: time,ff
      call fftd_3d (sz,a1,time,b1,b2,ff,lt)
      call fftd_3d (sz,a2,time,c1,c2,ff,lt)
      if(ndim.eq.3)call fftd_3d (sz,a3,time,d1,d2,ff,lt)
      return
      end
c----------------------------------------------------------------------
      subroutine fftd_3d (snap,a1,time,b1,b2,om,lt)
      implicit none
      integer snap, i, lt
      real, dimension(lt,snap)   :: a1
      real, dimension(lt,snap/2) :: b1,b2,om
      real, dimension(snap)      :: time
      do i = 1,lt
        call fftd(snap,a1(i,:),time,b1(i,:),b2(i,:),om)
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine fftd (n,a1,time,b1,b2,om)
      implicit none
      integer i,j,l,n
      real, dimension(n)   :: a1,time,temp
      real, dimension(n/2) :: b1,b2,om
      temp = a1
      call f_f_t_r_2(n,a1,a1,+1)
      a1 = a1*2; j = 1; b1 = 0.0d0; b2 = 0.0d0; om = 0.0d0
      do i = 1, n-1
        if (mod(i,2).eq.0) then
          b2(j) = a1(i)
        else
          b1(j) = a1(i)
          om(j) = (j-1)*8.0d0*atan(1.0d0)/(time(n)-time(1))
          j = j+1
        endif
      enddo
      a1 = temp
      b1(1) = b1(1)/2
      return
      end
c=======================================================================
c     SUBROUTINE FFTI(n,a1,Time,b1,b2,OM)
c     IMPLICIT NONE
c     INTEGER n,i,j
c     REAL*8 a1(n), Time(n)
c     REAL*8 b1(n/2), b2(n/2), OM(n/2)
c     a1 = b1(1)
c     DO i = 2,n/2
c     do j = 1,n
c     a1(j) = a1(j) + b1(i)*dcos(OM(i)*Time(j))+
c     $                      b2(i)*dsin(OM(i)*Time(j))
c     enddo
c     enddo
c     RETURN
c     END
c=======================================================================

!     ------------------------------------------------------------
!     Contenuto del file:
!     f_F_t_row_2 -------> routine di interfaccia per la
!     |        la libreria di Temperton
!     |--> Derive_f_c ---> Effettua la derivazione di
!     funzione nello spazio delle
!     trasformate
!     ------------------------------------------------------------

!     **********************************************************************
!     *                                                                    *
!     *     ***** ***** *****       *****        ***                       *
!     *     *     *       *         *   *       *   *                      *
!     *     ***   ***     *         ****           *                       *
!     *     *     *       *         *  *          *                        *
!     *     *     *       *   ***** *   * *****  ****                      *
!     *                                                                    *
!     **********************************************************************

      SUBROUTINE f_F_t_r_2(n,a1,a2,idir)
!     ------------------------------------------------------------
!     ATTENZIONE: FUNZIONA SOLO PER UN NUMERO PARI DI PUNTI  o COEFF
!     Questa routine costituisce un front-end per la routine che
!     calcola la trasformata di Fourier complessa.
!     E' realizzato per avere in ingresso due segnali reali o
!     i rispettivi coefficienti.
!     I segnali sono immagazzinati dal punto x= 0 in poi
!     e' ovviamente escluso l'ultimo punto dell'intervallo.
!     Dati in ingresso:
!     n: numero di coefficienti di Fourier--numero di punti
!     a1(0:n-1): vettore dei coefficienti
!     ordinati in questo modo: elementi pari corrispondono alla parte
!     reale elementi dispari alla parte immaginaria
!     o della funzione 1
!     a2(0:n-1): vettore dei coefficienti (vedi sopra) o della funzione 2
!     idir: se uguale a 1 effettua la trasformata se -1
!     antitrasforma
!
!     Dati in uscita:
!     a1(0:n-1): vettore della funzione o dei coefficienti 1
!     a2(0:n-1): vettore della funzione o dei coefficienti 2
C
C     Memoria di lavoro:
C     TRIGS(2*n), xr(0:n-1), xi(0:n-1)
!     ------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  idir, n, i
      REAL*8 scale, a1(0:n-1), a2(0:n-1)
      REAL*8  TRIGS(2*n),xr(0:n-1), xi(0:n-1)

      CALL SETGPFA(TRIGS,N)


!     Prepara il vettore da sottoporre alla trasformazione
      IF (idir.EQ.1) THEN
         DO i = 0, n-1
            xr(i) = a1(i)
            xi(i) = a2(i)
         END DO

      ELSE IF (idir.EQ.-1) THEN

         xr(0)  = a1(0)
         xr(n/2)  = a1(n-1)
         xi(0)  = a2(0)
         xi(n/2)  = a2(n-1)

         DO i=1,n/2-1

            xr (i) = a1(2*i) - a2(2*i-1)
            xr(n-i) = a1(2*i) + a2(2*i-1)
            xi (i) = a2(2*i) + a1(2*i-1)
            xi(n-i) = a2(2*i) - a1(2*i-1)

         END DO

      ELSE
         STOP 'Errore: il valore di idir e'' errato'
      END IF

!     Chiama la routine di libreria che calcola la trasformata
!     di Fourier
      scale = 1.0d0/(n**((max(0,idir))))
      DO i = 0, n-1
         xr(i) = scale * xr(i)
         xi(i) = scale * xi(i)
      END DO
      CALL GPFA(xr,xi,TRIGS,1,N,N,1,idir)

!     Separa i  risultati della FFT
!     Prepara i risultati da mandare in uscita
      IF (idir.EQ.1) THEN

         a1(0) = xr(0)
         a2(0) = xi(0)

         DO i=1,n/2-1

            a1(2*i)   = 0.5d0 * (   xr(i) + xr(n-i) )
            a2(2*i)   = 0.5d0 * (   xi(i) + xi(n-i) )

            a1(2*i-1) = 0.5d0 * (   xi(i) - xi(n-i) )
            a2(2*i-1) = 0.5d0 * ( - xr(i) + xr(n-i) )

         END DO

         a1(n-1) = xr(n/2)
         a2(n-1) = xi(n/2)

      ELSE IF (idir.EQ.-1) THEN
         DO i = 0, n-1
            a1(i) = xr(i)
            a2(i) = xi(i)
         END DO

      ELSE
         STOP 'Errore: il valore di idir e'' errato'
      END IF

      RETURN
      END

!     **********************************************************************
!     *                                                                    *
!     *     ****  ***** ****  * *   * *****       *****       *****        *
!     *     *   * *     *   * * *   * *           *           *            *
!     *     *   * ***   ****  * *   * ***         ***         *            *
!     *     *   * *     *  *  *  * *  *           *           *            *
!     *     ****  ***** *   * *   *   ***** ***** *     ***** *****        *
!     *                                                                    *
!     **********************************************************************

      SUBROUTINE Derive_f_c(n,deltax,f)
!     ----------------------------------------------------------------------
!     Versione 1.0. Release 9/12/1998.
!     Effettua la derivazione nello spazio delle frequenze della
!     funzione reale periodica f fornita mediante i coefficienti di Fourier.
!     In place!!
!     Dati in ingresso:
!     n: numero di coefficienti di Fourier
!     deltax: ampiezza dell'intervallo
!     f(0:n-1): coefficienti di Fourier cosi' ordinati:
!     0, 2, ..., 2i, ..., n-2, n-1 => parti reali dei coeff. di
!     frequenza i/2
!     1, 3, ..., 2i-1, ..., n-3    => parti immaginarie dei coeff.
!     di frequenza i/2
!     Risultati:
!     f(0:n-1): coefficienti di Fourier della derivata ordinati come sopra
!     ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  n, i
      REAL*8  temp, deltax, scale, pi
      PARAMETER (pi =  3.14159265358979323846d0)
      REAL*8  f(0:n-1)

      scale = 2.d0 * pi/deltax

      f(0) = 0.d0

      DO i = 2, n-2, 2
         temp   =  f(i)
         f(i)   =  f(i-1) * dfloat(i/2) * scale
         f(i-1) = -temp   * dfloat(i/2) * scale
      END DO

      f(n-1) = 0.d0

      RETURN
      END

!     **********************************************************************
!     *                                                                    *
!     *     ****  ***** ****  * *   * *****       *****       *****        *
!     *     *   * *     *   * * *   * *           *           *            *
!     *     *   * ***   ****  * *   * ***         ***         *            *
!     *     *   * *     *  *  *  * *  *           *           *            *
!     *     ****  ***** *   * *   *   ***** ***** *     ***** ***** *****  *
!     *                                                                    *
!     **********************************************************************

      SUBROUTINE Derive_f_c_n_i(n,deltax,f,df)
!     ----------------------------------------------------------------------
!     Versione 1.0. Release 26/1/1999.
!     Effettua la derivazione nello spazio delle frequenze della
!     funzione reale periodica f fornita mediante i coefficienti di Fourier.
!     In place!!
!     Dati in ingresso:
!     n: numero di coefficienti di Fourier
!     deltax: ampiezza dell'intervallo
!     f(0:n-1): coefficienti di Fourier cosi' ordinati:
!     0, 2, ..., 2i, ..., n-2, n-1 => parti reali dei coeff. di
!     frequenza i/2
!     1, 3, ..., 2i-1, ..., n-3    => parti immaginarie dei coeff.
!     di frequenza i/2
!     Risultati:
!     f(0:n-1): coefficienti di Fourier della derivata ordinati come sopra
!     ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  n, i
      REAL*8 deltax, scale, pi
      PARAMETER (pi =  3.14159265358979323846d0)
      REAL*8 f(0:*), df(0:*)

      scale = 2.d0 * pi/deltax

      df(0) = 0.d0

      DO i = 2, n-2, 2
         df(i)   =  f(i-1) * dfloat(i/2) * scale
         df(i-1) = -f(i)   * dfloat(i/2) * scale
      END DO

      df(n-1) = 0.d0

      RETURN
      END

!=======================================================================
C     SUBROUTINE 'SETGPFA'
C     SETUP ROUTINE FOR SELF-SORTING IN-PLACE
C     GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
C
C     CALL SETGPFA(TRIGS,N)
C
C     INPUT :
C     -----
C     N IS THE LENGTH OF THE TRANSFORMS. N MUST BE OF THE FORM:
C     -----------------------------------
C     N = (2**IP) * (3**IQ) * (5**IR)
C     -----------------------------------
C
C     OUTPUT:
C     ------
C     TRIGS IS A TABLE OF TWIDDLE FACTORS,
C     OF LENGTH 2*IPQR (REAL) WORDS, WHERE:
C     --------------------------------------
C     IPQR = (2**IP) + (3**IQ) + (5**IR)
C     --------------------------------------
C
C     WRITTEN BY CLIVE TEMPERTON 1990
C
C----------------------------------------------------------------------
C
      SUBROUTINE SETGPFA(TRIGS,N)
C
      IMPLICIT NONE
C     Incoming data
      INTEGER N

C     Outcoming data
      REAL*8 TRIGS(*)

C     Internal data
      INTEGER NJ(3), KK, LL, NN, NI, IFAC, I, K, IROT, IP, IQ, IR, KINK
      REAL*8  ANGLE, TWOPI, DEL
C
C     DECOMPOSE N INTO FACTORS 2,3,5
C     ------------------------------
      NN = N
      IFAC = 2
C
      DO 30 LL = 1 , 3
         KK = 0
 10      CONTINUE
         IF (MOD(NN,IFAC).NE.0) GO TO 20
         KK = KK + 1
         NN = NN / IFAC
         GO TO 10
 20      CONTINUE
         NJ(LL) = KK
         IFAC = IFAC + LL
 30   CONTINUE
C
      IF (NN.NE.1) THEN
         WRITE(6,40) N
 40      FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
         RETURN
      ENDIF
C
      IP = NJ(1)
      IQ = NJ(2)
      IR = NJ(3)
C
C     COMPUTE LIST OF ROTATED TWIDDLE FACTORS
C     ---------------------------------------
      NJ(1) = 2**IP
      NJ(2) = 3**IQ
      NJ(3) = 5**IR
C
      TWOPI = 4.0 * ASIN(1.d0)
      I = 1
C
      DO 60 LL = 1 , 3
         NI = NJ(LL)
         IF (NI.EQ.1) GO TO 60
C
         DEL = TWOPI / DFLOAT(NI)
         IROT = N / NI
         KINK = MOD(IROT,NI)
         KK = 0
C
         DO 50 K = 1 , NI
            ANGLE = DFLOAT(KK) * DEL
            TRIGS(I) = COS(ANGLE)
            TRIGS(I+1) = SIN(ANGLE)
            I = I + 2
            KK = KK + KINK
            IF (KK.GT.NI) KK = KK - NI
 50      CONTINUE
 60   CONTINUE
C
      RETURN
      END
!=======================================================================
C     **********************************************************************
C     **********************************************************************
C     SUBROUTINE 'GPFA'
C     SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT
C
C     *** THIS IS THE ALL-FORTRAN VERSION ***
C     -------------------------------
C
C     CALL GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)
C
C     A IS FIRST REAL INPUT/OUTPUT VECTOR
C     B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR
C     TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
C     BY CALLING SUBROUTINE 'SETGPFA'
C     INC IS THE INCREMENT WITHIN EACH DATA VECTOR
C     JUMP IS THE INCREMENT BETWEEN DATA VECTORS
C     N IS THE LENGTH OF THE TRANSFORMS:
C     -----------------------------------
C     N = (2**IP) * (3**IQ) * (5**IR)
C     -----------------------------------
C     LOT IS THE NUMBER OF TRANSFORMS
C     ISIGN = +1 FOR FORWARD TRANSFORM
C     = -1 FOR INVERSE TRANSFORM
C
C     WRITTEN BY CLIVE TEMPERTON
C     RECHERCHE EN PREVISION NUMERIQUE
C     ATMOSPHERIC ENVIRONMENT SERVICE, CANADA
C
C----------------------------------------------------------------------
C
C     DEFINITION OF TRANSFORM
C     -----------------------
C
C     X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
C
C---------------------------------------------------------------------
C
C     FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
C     SEE:
C
C     C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
C     FOR ANY N = (2**P)(3**Q)(5**R)",
C     SIAM J. SCI. STAT. COMP., MAY 1992.
C
C----------------------------------------------------------------------
C
      SUBROUTINE GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)
C
      IMPLICIT NONE

C     Incoming data
      INTEGER INC, JUMP, N, LOT, ISIGN
      REAL*8 TRIGS(*)

C     In-Outcoming data
      REAL*8 A(*), B(*)

C     Internal data
      INTEGER NJ(3), NN, IFAC, LL, KK, IP, IQ, IR, I
C
C     DECOMPOSE N INTO FACTORS 2,3,5
C     ------------------------------
      NN = N
      IFAC = 2
C
      DO 30 LL = 1 , 3
         KK = 0
 10      CONTINUE
         IF (MOD(NN,IFAC).NE.0) GO TO 20
         KK = KK + 1
         NN = NN / IFAC
         GO TO 10
 20      CONTINUE
         NJ(LL) = KK
         IFAC = IFAC + LL
 30   CONTINUE
C
      IF (NN.NE.1) THEN
         WRITE(6,40) N
 40      FORMAT(' *** WARNING!!',I10,' IS NOT A LEGAL VALUE OF N ***')
         RETURN
      ENDIF
C
      IP = NJ(1)
      IQ = NJ(2)
      IR = NJ(3)
C
C     COMPUTE THE TRANSFORM
C     ---------------------
      I = 1
      IF (IP.GT.0) THEN
         CALL GPFA2F(A,B,TRIGS,INC,JUMP,N,IP,LOT,ISIGN)
         I = I + 2 * ( 2**IP)
      ENDIF
      IF (IQ.GT.0) THEN
         CALL GPFA3F(A,B,TRIGS(I),INC,JUMP,N,IQ,LOT,ISIGN)
         I = I + 2 * (3**IQ)
      ENDIF
      IF (IR.GT.0) THEN
         CALL GPFA5F(A,B,TRIGS(I),INC,JUMP,N,IR,LOT,ISIGN)
      ENDIF
C
      RETURN
      END
!=======================================================================
C-------------------------------------------------------------------
C
      subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)
      IMPLICIT NONE

C     Incoming data
      INTEGER inc, jump, n, mm, lot, isign
      REAL*8 trigs(*)

C     In-Outcoming data
      REAL*8 a(*), b(*)

C     Internal data
      INTEGER lvr/1024/, n2, inq, jstepx, ninc, ink, m, m2, m8, mh,
     &     nblox, left, istart, nb, nvex, la, mu, nu, ipass, jstep,
     &     jstepl, jjj, ja, jb, jc, jd, je, jf, jg, jh, ji, jj, jk,
     &     jl, jm, jn, jo, jp,
     &     k, j, l, kk, laincl, ll
      REAL*8  s, ss, aja, ajb, ajc, ajd, aje, ajf, ajg, ajh, aji, ajj,
     &     ajk, ajl, ajm, ajn, ajo, ajp,
     &     t0, t1, t2, t3, bja, bjb, bjc, bjd, bje, bjf, bjg, bjh,
     &     bji, bjj, bjk, bjl, bjm, bjn, bjo, bjp,
     &     u0, u1, u2, u3, co1, co2, co3, co4, co5, co6, co7,
     &     si1, si2, si3, si4, si5, si6, si7,
     &     c1, c2, c3

C
C     ***************************************************************
C     *                                                             *
C     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
C     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
C     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
C     *                                                             *
C     ***************************************************************
C
      n2 = 2**mm
      inq = n/n2
      jstepx = (n2-n) * inc
      ninc = n * inc
      ink = inc * inq
C
      m2 = 0
      m8 = 0
      if (mod(mm,2).eq.0) then
         m = mm/2
      else if (mod(mm,4).eq.1) then
         m = (mm-1)/2
         m2 = 1
      else if (mod(mm,4).eq.3) then
         m = (mm-3)/2
         m8 = 1
      endif
      mh = (m+1)/2
C
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = dfloat(isign)
      istart = 1
C
C     loop on blocks of lvr transforms
C     --------------------------------
      do 500 nb = 1 , nblox
C
         if (left.le.lvr) then
            nvex = left
         else if (left.lt.(2*lvr)) then
            nvex = left/2
            nvex = nvex + mod(nvex,2)
         else
            nvex = lvr
         endif
         left = left - nvex
C
         la = 1
C
C     loop on type I radix-4 passes
C     -----------------------------
         mu = mod(inq,4)
         if (isign.eq.-1) mu = 4 - mu
         ss = 1.0d0
         if (mu.eq.3) ss = -1.0d0
C
         if (mh.eq.0) go to 200
C
         do 160 ipass = 1 , mh
            jstep = (n*inc) / (4*la)
            jstepl = jstep - ninc
C
C     k = 0 loop (no twiddle factors)
C     -------------------------------
            do 120 jjj = 0 , (n-1)*inc , 4*jstep
               ja = istart + jjj
C
C     "transverse" loop
C     -----------------
               do 115 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc.lt.istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd.lt.istart) jd = jd + ninc
                  j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
                  do 110 l = 1 , nvex
                     aja = a(ja+j)
                     ajc = a(jc+j)
                     t0 = aja + ajc
                     t2 = aja - ajc
                     ajb = a(jb+j)
                     ajd = a(jd+j)
                     t1 = ajb + ajd
                     t3 = ss * ( ajb - ajd )
                     bja = b(ja+j)
                     bjc = b(jc+j)
                     u0 = bja + bjc
                     u2 = bja - bjc
                     bjb = b(jb+j)
                     bjd = b(jd+j)
                     u1 = bjb + bjd
                     u3 = ss * ( bjb - bjd )
                     a(ja+j) = t0 + t1
                     a(jc+j) = t0 - t1
                     b(ja+j) = u0 + u1
                     b(jc+j) = u0 - u1
                     a(jb+j) = t2 - u3
                     a(jd+j) = t2 + u3
                     b(jb+j) = u2 + t3
                     b(jd+j) = u2 - t3
                     j = j + jump
 110              continue
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
 115           continue
 120        continue
C
C     finished if n2 = 4
C     ------------------
            if (n2.eq.4) go to 490
            kk = 2 * la
C
C     loop on nonzero k
C     -----------------
            do 150 k = ink , jstep-ink , ink
               co1 = trigs(kk+1)
               si1 = s*trigs(kk+2)
               co2 = trigs(2*kk+1)
               si2 = s*trigs(2*kk+2)
               co3 = trigs(3*kk+1)
               si3 = s*trigs(3*kk+2)
C
C     loop along transform
C     --------------------
               do 140 jjj = k , (n-1)*inc , 4*jstep
                  ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                  do 135 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd.lt.istart) jd = jd + ninc
                     j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep,shortloop
                     do 130 l = 1 , nvex
                        aja = a(ja+j)
                        ajc = a(jc+j)
                        t0 = aja + ajc
                        t2 = aja - ajc
                        ajb = a(jb+j)
                        ajd = a(jd+j)
                        t1 = ajb + ajd
                        t3 = ss * ( ajb - ajd )
                        bja = b(ja+j)
                        bjc = b(jc+j)
                        u0 = bja + bjc
                        u2 = bja - bjc
                        bjb = b(jb+j)
                        bjd = b(jd+j)
                        u1 = bjb + bjd
                        u3 = ss * ( bjb - bjd )
                        a(ja+j) = t0 + t1
                        b(ja+j) = u0 + u1
                        a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                        b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                        a(jc+j) = co2*(t0-t1) - si2*(u0-u1)
                        b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
                        a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
                        b(jd+j) = si3*(t2+u3) + co3*(u2-t3)
                        j = j + jump
 130                 continue
C-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
 135              continue
 140           continue
C-----( end of loop along transforms )
               kk = kk + 2*la
 150        continue
C-----( end of loop on nonzero k )
            la = 4*la
 160     continue
C-----( end of loop on type I radix-4 passes)
C
C     central radix-2 pass
C     --------------------
 200     continue
         if (m2.eq.0) go to 300
C
         jstep = (n*inc) / (2*la)
         jstepl = jstep - ninc
C
C     k=0 loop (no twiddle factors)
C     -----------------------------
         do 220 jjj = 0 , (n-1)*inc , 2*jstep
            ja = istart + jjj
C
C     "transverse" loop
C     -----------------
            do 215 nu = 1 , inq
               jb = ja + jstepl
               if (jb.lt.istart) jb = jb + ninc
               j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
               do 210 l = 1 , nvex
                  aja = a(ja+j)
                  ajb = a(jb+j)
                  t0 = aja - ajb
                  a(ja+j) = aja + ajb
                  a(jb+j) = t0
                  bja = b(ja+j)
                  bjb = b(jb+j)
                  u0 = bja - bjb
                  b(ja+j) = bja + bjb
                  b(jb+j) = u0
                  j = j + jump
 210           continue
C-----(end of loop across transforms)
               ja = ja + jstepx
               if (ja.lt.istart) ja = ja + ninc
 215        continue
 220     continue
C
C     finished if n2=2
C     ----------------
         if (n2.eq.2) go to 490
C
         kk = 2 * la
C
C     loop on nonzero k
C     -----------------
         do 260 k = ink , jstep - ink , ink
            co1 = trigs(kk+1)
            si1 = s*trigs(kk+2)
C
C     loop along transforms
C     ---------------------
            do 250 jjj = k , (n-1)*inc , 2*jstep
               ja = istart + jjj
C
C     "transverse" loop
C     -----------------
               do 245 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  j = 0
C
C     loop across transforms
C     ----------------------
                  if (kk.eq.n2/2) then
c     dir$ ivdep, shortloop
                     do 230 l = 1 , nvex
                        aja = a(ja+j)
                        ajb = a(jb+j)
                        t0 = ss * ( aja - ajb )
                        a(ja+j) = aja + ajb
                        bjb = b(jb+j)
                        bja = b(ja+j)
                        a(jb+j) = ss * ( bjb - bja )
                        b(ja+j) = bja + bjb
                        b(jb+j) = t0
                        j = j + jump
 230                 continue
C
                  else
C
c     dir$ ivdep, shortloop
                     do 240 l = 1 , nvex
                        aja = a(ja+j)
                        ajb = a(jb+j)
                        t0 = aja - ajb
                        a(ja+j) = aja + ajb
                        bja = b(ja+j)
                        bjb = b(jb+j)
                        u0 = bja - bjb
                        b(ja+j) = bja + bjb
                        a(jb+j) = co1*t0 - si1*u0
                        b(jb+j) = si1*t0 + co1*u0
                        j = j + jump
 240                 continue
C
                  endif
C
C-----(end of loop across transforms)
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
 245           continue
 250        continue
C-----(end of loop along transforms)
            kk = kk + 2 * la
 260     continue
C-----(end of loop on nonzero k)
C-----(end of radix-2 pass)
C
         la = 2 * la
         go to 400
C
C     central radix-8 pass
C     --------------------
 300     continue
         if (m8.eq.0) go to 400
         jstep = (n*inc) / (8*la)
         jstepl = jstep - ninc
         mu = mod(inq,8)
         if (isign.eq.-1) mu = 8 - mu
         c1 = 1.0d0
         if (mu.eq.3.or.mu.eq.7) c1 = -1.0d0
         c2 = sqrt(0.5d0)
         if (mu.eq.3.or.mu.eq.5) c2 = -c2
         c3 = c1 * c2
C
C     stage 1
C     -------
         do 320 k = 0 , jstep - ink , ink
            do 315 jjj = k , (n-1)*inc , 8*jstep
               ja = istart + jjj
C
C     "transverse" loop
C     -----------------
               do 312 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc.lt.istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd.lt.istart) jd = jd + ninc
                  je = jd + jstepl
                  if (je.lt.istart) je = je + ninc
                  jf = je + jstepl
                  if (jf.lt.istart) jf = jf + ninc
                  jg = jf + jstepl
                  if (jg.lt.istart) jg = jg + ninc
                  jh = jg + jstepl
                  if (jh.lt.istart) jh = jh + ninc
                  j = 0
c     dir$ ivdep, shortloop
                  do 310 l = 1 , nvex
                     aja = a(ja+j)
                     aje = a(je+j)
                     t0 = aja - aje
                     a(ja+j) = aja + aje
                     ajc = a(jc+j)
                     ajg = a(jg+j)
                     t1 = c1 * ( ajc - ajg )
                     a(je+j) = ajc + ajg
                     ajb = a(jb+j)
                     ajf = a(jf+j)
                     t2 = ajb - ajf
                     a(jc+j) = ajb + ajf
                     ajd = a(jd+j)
                     ajh = a(jh+j)
                     t3 = ajd - ajh
                     a(jg+j) = ajd + ajh
                     a(jb+j) = t0
                     a(jf+j) = t1
                     a(jd+j) = c2 * ( t2 - t3 )
                     a(jh+j) = c3 * ( t2 + t3 )
                     bja = b(ja+j)
                     bje = b(je+j)
                     u0 = bja - bje
                     b(ja+j) = bja + bje
                     bjc = b(jc+j)
                     bjg = b(jg+j)
                     u1 = c1 * ( bjc - bjg )
                     b(je+j) = bjc + bjg
                     bjb = b(jb+j)
                     bjf = b(jf+j)
                     u2 = bjb - bjf
                     b(jc+j) = bjb + bjf
                     bjd = b(jd+j)
                     bjh = b(jh+j)
                     u3 = bjd - bjh
                     b(jg+j) = bjd + bjh
                     b(jb+j) = u0
                     b(jf+j) = u1
                     b(jd+j) = c2 * ( u2 - u3 )
                     b(jh+j) = c3 * ( u2 + u3 )
                     j = j + jump
 310              continue
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
 312           continue
 315        continue
 320     continue
C
C     stage 2
C     -------
C
C     k=0 (no twiddle factors)
C     ------------------------
         do 330 jjj = 0 , (n-1)*inc , 8*jstep
            ja = istart + jjj
C
C     "transverse" loop
C     -----------------
            do 328 nu = 1 , inq
               jb = ja + jstepl
               if (jb.lt.istart) jb = jb + ninc
               jc = jb + jstepl
               if (jc.lt.istart) jc = jc + ninc
               jd = jc + jstepl
               if (jd.lt.istart) jd = jd + ninc
               je = jd + jstepl
               if (je.lt.istart) je = je + ninc
               jf = je + jstepl
               if (jf.lt.istart) jf = jf + ninc
               jg = jf + jstepl
               if (jg.lt.istart) jg = jg + ninc
               jh = jg + jstepl
               if (jh.lt.istart) jh = jh + ninc
               j = 0
c     dir$ ivdep, shortloop
               do 325 l = 1 , nvex
                  aja = a(ja+j)
                  aje = a(je+j)
                  t0 = aja + aje
                  t2 = aja - aje
                  ajc = a(jc+j)
                  ajg = a(jg+j)
                  t1 = ajc + ajg
                  t3 = c1 * ( ajc - ajg )
                  bja = b(ja+j)
                  bje = b(je+j)
                  u0 = bja + bje
                  u2 = bja - bje
                  bjc = b(jc+j)
                  bjg = b(jg+j)
                  u1 = bjc + bjg
                  u3 = c1 * ( bjc - bjg )
                  a(ja+j) = t0 + t1
                  a(je+j) = t0 - t1
                  b(ja+j) = u0 + u1
                  b(je+j) = u0 - u1
                  a(jc+j) = t2 - u3
                  a(jg+j) = t2 + u3
                  b(jc+j) = u2 + t3
                  b(jg+j) = u2 - t3
                  ajb = a(jb+j)
                  ajd = a(jd+j)
                  t0 = ajb + ajd
                  t2 = ajb - ajd
                  ajf = a(jf+j)
                  ajh = a(jh+j)
                  t1 = ajf - ajh
                  t3 = ajf + ajh
                  bjb = b(jb+j)
                  bjd = b(jd+j)
                  u0 = bjb + bjd
                  u2 = bjb - bjd
                  bjf = b(jf+j)
                  bjh = b(jh+j)
                  u1 = bjf - bjh
                  u3 = bjf + bjh
                  a(jb+j) = t0 - u3
                  a(jh+j) = t0 + u3
                  b(jb+j) = u0 + t3
                  b(jh+j) = u0 - t3
                  a(jd+j) = t2 + u1
                  a(jf+j) = t2 - u1
                  b(jd+j) = u2 - t1
                  b(jf+j) = u2 + t1
                  j = j + jump
 325           continue
               ja = ja + jstepx
               if (ja.lt.istart) ja = ja + ninc
 328        continue
 330     continue
C
         if (n2.eq.8) go to 490
C
C     loop on nonzero k
C     -----------------
         kk = 2 * la
C
         do 350 k = ink , jstep - ink , ink
C
            co1 = trigs(kk+1)
            si1 = s * trigs(kk+2)
            co2 = trigs(2*kk+1)
            si2 = s * trigs(2*kk+2)
            co3 = trigs(3*kk+1)
            si3 = s * trigs(3*kk+2)
            co4 = trigs(4*kk+1)
            si4 = s * trigs(4*kk+2)
            co5 = trigs(5*kk+1)
            si5 = s * trigs(5*kk+2)
            co6 = trigs(6*kk+1)
            si6 = s * trigs(6*kk+2)
            co7 = trigs(7*kk+1)
            si7 = s * trigs(7*kk+2)
C
            do 345 jjj = k , (n-1)*inc , 8*jstep
               ja = istart + jjj
C
C     "transverse" loop
C     -----------------
               do 342 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc.lt.istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd.lt.istart) jd = jd + ninc
                  je = jd + jstepl
                  if (je.lt.istart) je = je + ninc
                  jf = je + jstepl
                  if (jf.lt.istart) jf = jf + ninc
                  jg = jf + jstepl
                  if (jg.lt.istart) jg = jg + ninc
                  jh = jg + jstepl
                  if (jh.lt.istart) jh = jh + ninc
                  j = 0
c     dir$ ivdep, shortloop
                  do 340 l = 1 , nvex
                     aja = a(ja+j)
                     aje = a(je+j)
                     t0 = aja + aje
                     t2 = aja - aje
                     ajc = a(jc+j)
                     ajg = a(jg+j)
                     t1 = ajc + ajg
                     t3 = c1 * ( ajc - ajg )
                     bja = b(ja+j)
                     bje = b(je+j)
                     u0 = bja + bje
                     u2 = bja - bje
                     bjc = b(jc+j)
                     bjg = b(jg+j)
                     u1 = bjc + bjg
                     u3 = c1 * ( bjc - bjg )
                     a(ja+j) = t0 + t1
                     b(ja+j) = u0 + u1
                     a(je+j) = co4*(t0-t1) - si4*(u0-u1)
                     b(je+j) = si4*(t0-t1) + co4*(u0-u1)
                     a(jc+j) = co2*(t2-u3) - si2*(u2+t3)
                     b(jc+j) = si2*(t2-u3) + co2*(u2+t3)
                     a(jg+j) = co6*(t2+u3) - si6*(u2-t3)
                     b(jg+j) = si6*(t2+u3) + co6*(u2-t3)
                     ajb = a(jb+j)
                     ajd = a(jd+j)
                     t0 = ajb + ajd
                     t2 = ajb - ajd
                     ajf = a(jf+j)
                     ajh = a(jh+j)
                     t1 = ajf - ajh
                     t3 = ajf + ajh
                     bjb = b(jb+j)
                     bjd = b(jd+j)
                     u0 = bjb + bjd
                     u2 = bjb - bjd
                     bjf = b(jf+j)
                     bjh = b(jh+j)
                     u1 = bjf - bjh
                     u3 = bjf + bjh
                     a(jb+j) = co1*(t0-u3) - si1*(u0+t3)
                     b(jb+j) = si1*(t0-u3) + co1*(u0+t3)
                     a(jh+j) = co7*(t0+u3) - si7*(u0-t3)
                     b(jh+j) = si7*(t0+u3) + co7*(u0-t3)
                     a(jd+j) = co3*(t2+u1) - si3*(u2-t1)
                     b(jd+j) = si3*(t2+u1) + co3*(u2-t1)
                     a(jf+j) = co5*(t2-u1) - si5*(u2+t1)
                     b(jf+j) = si5*(t2-u1) + co5*(u2+t1)
                     j = j + jump
 340              continue
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
 342           continue
 345        continue
            kk = kk + 2 * la
 350     continue
C
         la = 8 * la
C
C     loop on type II radix-4 passes
C     ------------------------------
 400     continue
         mu = mod(inq,4)
         if (isign.eq.-1) mu = 4 - mu
         ss = 1.0d0
         if (mu.eq.3) ss = -1.0d0
C
         do 480 ipass = mh+1 , m
            jstep = (n*inc) / (4*la)
            jstepl = jstep - ninc
            laincl = la * ink - ninc
C
C     k=0 loop (no twiddle factors)
C     -----------------------------
            do 430 ll = 0 , (la-1)*ink , 4*jstep
C
               do 420 jjj = ll , (n-1)*inc , 4*la*ink
                  ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                  do 415 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd.lt.istart) jd = jd + ninc
                     je = ja + laincl
                     if (je.lt.istart) je = je + ninc
                     jf = je + jstepl
                     if (jf.lt.istart) jf = jf + ninc
                     jg = jf + jstepl
                     if (jg.lt.istart) jg = jg + ninc
                     jh = jg + jstepl
                     if (jh.lt.istart) jh = jh + ninc
                     ji = je + laincl
                     if (ji.lt.istart) ji = ji + ninc
                     jj = ji + jstepl
                     if (jj.lt.istart) jj = jj + ninc
                     jk = jj + jstepl
                     if (jk.lt.istart) jk = jk + ninc
                     jl = jk + jstepl
                     if (jl.lt.istart) jl = jl + ninc
                     jm = ji + laincl
                     if (jm.lt.istart) jm = jm + ninc
                     jn = jm + jstepl
                     if (jn.lt.istart) jn = jn + ninc
                     jo = jn + jstepl
                     if (jo.lt.istart) jo = jo + ninc
                     jp = jo + jstepl
                     if (jp.lt.istart) jp = jp + ninc
                     j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
                     do 410 l = 1 , nvex
                        aja = a(ja+j)
                        ajc = a(jc+j)
                        t0 = aja + ajc
                        t2 = aja - ajc
                        ajb = a(jb+j)
                        ajd = a(jd+j)
                        t1 = ajb + ajd
                        t3 = ss * ( ajb - ajd )
                        aji = a(ji+j)
                        ajc =  aji
                        bja = b(ja+j)
                        bjc = b(jc+j)
                        u0 = bja + bjc
                        u2 = bja - bjc
                        bjb = b(jb+j)
                        bjd = b(jd+j)
                        u1 = bjb + bjd
                        u3 = ss * ( bjb - bjd )
                        aje = a(je+j)
                        ajb =  aje
                        a(ja+j) = t0 + t1
                        a(ji+j) = t0 - t1
                        b(ja+j) = u0 + u1
                        bjc =  u0 - u1
                        bjm = b(jm+j)
                        bjd =  bjm
                        a(je+j) = t2 - u3
                        ajd =  t2 + u3
                        bjb =  u2 + t3
                        b(jm+j) = u2 - t3
C----------------------
                        ajg = a(jg+j)
                        t0 = ajb + ajg
                        t2 = ajb - ajg
                        ajf = a(jf+j)
                        ajh = a(jh+j)
                        t1 = ajf + ajh
                        t3 = ss * ( ajf - ajh )
                        ajj = a(jj+j)
                        ajg =  ajj
                        bje = b(je+j)
                        bjg = b(jg+j)
                        u0 = bje + bjg
                        u2 = bje - bjg
                        bjf = b(jf+j)
                        bjh = b(jh+j)
                        u1 = bjf + bjh
                        u3 = ss * ( bjf - bjh )
                        b(je+j) = bjb
                        a(jb+j) = t0 + t1
                        a(jj+j) = t0 - t1
                        bjj = b(jj+j)
                        bjg =  bjj
                        b(jb+j) = u0 + u1
                        b(jj+j) = u0 - u1
                        a(jf+j) = t2 - u3
                        ajh =  t2 + u3
                        b(jf+j) = u2 + t3
                        bjh =  u2 - t3
C----------------------
                        ajk = a(jk+j)
                        t0 = ajc + ajk
                        t2 = ajc - ajk
                        ajl = a(jl+j)
                        t1 = ajg + ajl
                        t3 = ss * ( ajg - ajl )
                        bji = b(ji+j)
                        bjk = b(jk+j)
                        u0 = bji + bjk
                        u2 = bji - bjk
                        ajo = a(jo+j)
                        ajl =  ajo
                        bjl = b(jl+j)
                        u1 = bjg + bjl
                        u3 = ss * ( bjg - bjl )
                        b(ji+j) = bjc
                        a(jc+j) = t0 + t1
                        a(jk+j) = t0 - t1
                        bjo = b(jo+j)
                        bjl =  bjo
                        b(jc+j) = u0 + u1
                        b(jk+j) = u0 - u1
                        a(jg+j) = t2 - u3
                        a(jo+j) = t2 + u3
                        b(jg+j) = u2 + t3
                        b(jo+j) = u2 - t3
C----------------------
                        ajm = a(jm+j)
                        t0 = ajm + ajl
                        t2 = ajm - ajl
                        ajn = a(jn+j)
                        ajp = a(jp+j)
                        t1 = ajn + ajp
                        t3 = ss * ( ajn - ajp )
                        a(jm+j) = ajd
                        u0 = bjd + bjl
                        u2 = bjd - bjl
                        bjn = b(jn+j)
                        bjp = b(jp+j)
                        u1 = bjn + bjp
                        u3 = ss * ( bjn - bjp )
                        a(jn+j) = ajh
                        a(jd+j) = t0 + t1
                        a(jl+j) = t0 - t1
                        b(jd+j) = u0 + u1
                        b(jl+j) = u0 - u1
                        b(jn+j) = bjh
                        a(jh+j) = t2 - u3
                        a(jp+j) = t2 + u3
                        b(jh+j) = u2 + t3
                        b(jp+j) = u2 - t3
                        j = j + jump
 410                 continue
C-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
 415              continue
 420           continue
 430        continue
C-----( end of double loop for k=0 )
C
C     finished if last pass
C     ---------------------
            if (ipass.eq.m) go to 490
C
            kk = 2*la
C
C     loop on nonzero k
C     -----------------
            do 470 k = ink , jstep-ink , ink
               co1 = trigs(kk+1)
               si1 = s*trigs(kk+2)
               co2 = trigs(2*kk+1)
               si2 = s*trigs(2*kk+2)
               co3 = trigs(3*kk+1)
               si3 = s*trigs(3*kk+2)
C
C     double loop along first transform in block
C     ------------------------------------------
               do 460 ll = k , (la-1)*ink , 4*jstep
C
                  do 450 jjj = ll , (n-1)*inc , 4*la*ink
                     ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                     do 445 nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        jd = jc + jstepl
                        if (jd.lt.istart) jd = jd + ninc
                        je = ja + laincl
                        if (je.lt.istart) je = je + ninc
                        jf = je + jstepl
                        if (jf.lt.istart) jf = jf + ninc
                        jg = jf + jstepl
                        if (jg.lt.istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh.lt.istart) jh = jh + ninc
                        ji = je + laincl
                        if (ji.lt.istart) ji = ji + ninc
                        jj = ji + jstepl
                        if (jj.lt.istart) jj = jj + ninc
                        jk = jj + jstepl
                        if (jk.lt.istart) jk = jk + ninc
                        jl = jk + jstepl
                        if (jl.lt.istart) jl = jl + ninc
                        jm = ji + laincl
                        if (jm.lt.istart) jm = jm + ninc
                        jn = jm + jstepl
                        if (jn.lt.istart) jn = jn + ninc
                        jo = jn + jstepl
                        if (jo.lt.istart) jo = jo + ninc
                        jp = jo + jstepl
                        if (jp.lt.istart) jp = jp + ninc
                        j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
                        do 440 l = 1 , nvex
                           aja = a(ja+j)
                           ajc = a(jc+j)
                           t0 = aja + ajc
                           t2 = aja - ajc
                           ajb = a(jb+j)
                           ajd = a(jd+j)
                           t1 = ajb + ajd
                           t3 = ss * ( ajb - ajd )
                           aji = a(ji+j)
                           ajc =  aji
                           bja = b(ja+j)
                           bjc = b(jc+j)
                           u0 = bja + bjc
                           u2 = bja - bjc
                           bjb = b(jb+j)
                           bjd = b(jd+j)
                           u1 = bjb + bjd
                           u3 = ss * ( bjb - bjd )
                           aje = a(je+j)
                           ajb =  aje
                           a(ja+j) = t0 + t1
                           b(ja+j) = u0 + u1
                           a(je+j) = co1*(t2-u3) - si1*(u2+t3)
                           bjb =  si1*(t2-u3) + co1*(u2+t3)
                           bjm = b(jm+j)
                           bjd =  bjm
                           a(ji+j) = co2*(t0-t1) - si2*(u0-u1)
                           bjc =  si2*(t0-t1) + co2*(u0-u1)
                           ajd =  co3*(t2+u3) - si3*(u2-t3)
                           b(jm+j) = si3*(t2+u3) + co3*(u2-t3)
C----------------------------------------
                           ajg = a(jg+j)
                           t0 = ajb + ajg
                           t2 = ajb - ajg
                           ajf = a(jf+j)
                           ajh = a(jh+j)
                           t1 = ajf + ajh
                           t3 = ss * ( ajf - ajh )
                           ajj = a(jj+j)
                           ajg =  ajj
                           bje = b(je+j)
                           bjg = b(jg+j)
                           u0 = bje + bjg
                           u2 = bje - bjg
                           bjf = b(jf+j)
                           bjh = b(jh+j)
                           u1 = bjf + bjh
                           u3 = ss * ( bjf - bjh )
                           b(je+j) = bjb
                           a(jb+j) = t0 + t1
                           b(jb+j) = u0 + u1
                           bjj = b(jj+j)
                           bjg =  bjj
                           a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                           b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                           a(jj+j) = co2*(t0-t1) - si2*(u0-u1)
                           b(jj+j) = si2*(t0-t1) + co2*(u0-u1)
                           ajh =  co3*(t2+u3) - si3*(u2-t3)
                           bjh =  si3*(t2+u3) + co3*(u2-t3)
C----------------------------------------
                           ajk = a(jk+j)
                           t0 = ajc + ajk
                           t2 = ajc - ajk
                           ajl = a(jl+j)
                           t1 = ajg + ajl
                           t3 = ss * ( ajg - ajl )
                           bji = b(ji+j)
                           bjk = b(jk+j)
                           u0 = bji + bjk
                           u2 = bji - bjk
                           ajo = a(jo+j)
                           ajl =  ajo
                           bjl = b(jl+j)
                           u1 = bjg + bjl
                           u3 = ss * ( bjg - bjl )
                           b(ji+j) = bjc
                           a(jc+j) = t0 + t1
                           b(jc+j) = u0 + u1
                           bjo = b(jo+j)
                           bjl =  bjo
                           a(jg+j) = co1*(t2-u3) - si1*(u2+t3)
                           b(jg+j) = si1*(t2-u3) + co1*(u2+t3)
                           a(jk+j) = co2*(t0-t1) - si2*(u0-u1)
                           b(jk+j) = si2*(t0-t1) + co2*(u0-u1)
                           a(jo+j) = co3*(t2+u3) - si3*(u2-t3)
                           b(jo+j) = si3*(t2+u3) + co3*(u2-t3)
C----------------------------------------
                           ajm = a(jm+j)
                           t0 = ajm + ajl
                           t2 = ajm - ajl
                           ajn = a(jn+j)
                           ajp = a(jp+j)
                           t1 = ajn + ajp
                           t3 = ss * ( ajn - ajp )
                           a(jm+j) = ajd
                           u0 = bjd + bjl
                           u2 = bjd - bjl
                           a(jn+j) = ajh
                           bjn = b(jn+j)
                           bjp = b(jp+j)
                           u1 = bjn + bjp
                           u3 = ss * ( bjn - bjp )
                           b(jn+j) = bjh
                           a(jd+j) = t0 + t1
                           b(jd+j) = u0 + u1
                           a(jh+j) = co1*(t2-u3) - si1*(u2+t3)
                           b(jh+j) = si1*(t2-u3) + co1*(u2+t3)
                           a(jl+j) = co2*(t0-t1) - si2*(u0-u1)
                           b(jl+j) = si2*(t0-t1) + co2*(u0-u1)
                           a(jp+j) = co3*(t2+u3) - si3*(u2-t3)
                           b(jp+j) = si3*(t2+u3) + co3*(u2-t3)
                           j = j + jump
 440                    continue
C-----(end of loop across transforms)
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
 445                 continue
 450              continue
 460           continue
C-----( end of double loop for this k )
               kk = kk + 2*la
 470        continue
C-----( end of loop over values of k )
            la = 4*la
 480     continue
C-----( end of loop on type II radix-4 passes )
C-----( nvex transforms completed)
 490     continue
         istart = istart + nvex * jump
 500  continue
C-----( end of loop on blocks of transforms )
C
      return
      end
!=======================================================================
C     fortran version of *gpfa3* -
C     radix-3 section of self-sorting, in-place
C     generalized PFA
C
C-------------------------------------------------------------------
C
      subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)

      IMPLICIT NONE
C     Incoming data
      INTEGER inc, jump, n, mm, lot, isign
      REAL*8 trigs(*)

C     In-Outcoming data
      REAL*8 a(*), b(*)

C     Internal data
      INTEGER lvr/128/, n3, inq, jstepx, ninc, ink, mu, m, mh, nblox,
     &     left, istart, nb, nvex, la, ipass, jstep, jstepl, jjj,
     &     ja, nu, jb, jc, jd, je, jf, jg, jh, ji, j, l, ll, k, kk,
     &     laincl
      REAL*8  sin60, s, c1,
     &     aja, ajb, ajc, ajd, aje, ajf, ajg, ajh, aji,
     &     t1, t2, t3,
     &     bja, bjb, bjc, bjd, bje, bjf, bjg, bjh, bji,
     &     u1, u2, u3, co1, co2, si1, si2

      sin60 = SQRT(3.0d0)/2.0d0

C     ***************************************************************
C     *                                                             *
C     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
C     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
C     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
C     *                                                             *
C     ***************************************************************
C
      n3 = 3**mm
      inq = n/n3
      jstepx = (n3-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,3)
      if (isign.eq.-1) mu = 3-mu
      m = mm
      mh = (m+1)/2
      s = dfloat(isign)
      c1 = sin60
      if (mu.eq.2) c1 = -c1
C
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = dfloat(isign)
      istart = 1
C
C     loop on blocks of lvr transforms
C     --------------------------------
      do 500 nb = 1 , nblox
C
         if (left.le.lvr) then
            nvex = left
         else if (left.lt.(2*lvr)) then
            nvex = left/2
            nvex = nvex + mod(nvex,2)
         else
            nvex = lvr
         endif
         left = left - nvex
C
         la = 1
C
C     loop on type I radix-3 passes
C     -----------------------------
         do 160 ipass = 1 , mh
            jstep = (n*inc) / (3*la)
            jstepl = jstep - ninc
C
C     k = 0 loop (no twiddle factors)
C     -------------------------------
            do 120 jjj = 0 , (n-1)*inc , 3*jstep
               ja = istart + jjj
C
C     "transverse" loop
C     -----------------
               do 115 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc.lt.istart) jc = jc + ninc
                  j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
                  do 110 l = 1 , nvex
                     ajb = a(jb+j)
                     ajc = a(jc+j)
                     t1 = ajb + ajc
                     aja = a(ja+j)
                     t2 = aja - 0.5d0 * t1
                     t3 = c1 * ( ajb - ajc )
                     bjb = b(jb+j)
                     bjc = b(jc+j)
                     u1 = bjb + bjc
                     bja = b(ja+j)
                     u2 = bja - 0.5d0 * u1
                     u3 = c1 * ( bjb - bjc )
                     a(ja+j) = aja + t1
                     b(ja+j) = bja + u1
                     a(jb+j) = t2 - u3
                     b(jb+j) = u2 + t3
                     a(jc+j) = t2 + u3
                     b(jc+j) = u2 - t3
                     j = j + jump
 110              continue
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
 115           continue
 120        continue
C
C     finished if n3 = 3
C     ------------------
            if (n3.eq.3) go to 490
            kk = 2 * la
C
C     loop on nonzero k
C     -----------------
            do 150 k = ink , jstep-ink , ink
               co1 = trigs(kk+1)
               si1 = s*trigs(kk+2)
               co2 = trigs(2*kk+1)
               si2 = s*trigs(2*kk+2)
C
C     loop along transform
C     --------------------
               do 140 jjj = k , (n-1)*inc , 3*jstep
                  ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                  do 135 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep,shortloop
                     do 130 l = 1 , nvex
                        ajb = a(jb+j)
                        ajc = a(jc+j)
                        t1 = ajb + ajc
                        aja = a(ja+j)
                        t2 = aja - 0.5d0 * t1
                        t3 = c1 * ( ajb - ajc )
                        bjb = b(jb+j)
                        bjc = b(jc+j)
                        u1 = bjb + bjc
                        bja = b(ja+j)
                        u2 = bja - 0.5d0 * u1
                        u3 = c1 * ( bjb - bjc )
                        a(ja+j) = aja + t1
                        b(ja+j) = bja + u1
                        a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                        b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                        a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
                        b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
                        j = j + jump
 130                 continue
C-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
 135              continue
 140           continue
C-----( end of loop along transforms )
               kk = kk + 2*la
 150        continue
C-----( end of loop on nonzero k )
            la = 3*la
 160     continue
C-----( end of loop on type I radix-3 passes)
C
C     loop on type II radix-3 passes
C     ------------------------------
 400     continue
C
         do 480 ipass = mh+1 , m
            jstep = (n*inc) / (3*la)
            jstepl = jstep - ninc
            laincl = la*ink - ninc
C
C     k=0 loop (no twiddle factors)
C     -----------------------------
            do 430 ll = 0 , (la-1)*ink , 3*jstep
C
               do 420 jjj = ll , (n-1)*inc , 3*la*ink
                  ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                  do 415 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     jd = ja + laincl
                     if (jd.lt.istart) jd = jd + ninc
                     je = jd + jstepl
                     if (je.lt.istart) je = je + ninc
                     jf = je + jstepl
                     if (jf.lt.istart) jf = jf + ninc
                     jg = jd + laincl
                     if (jg.lt.istart) jg = jg + ninc
                     jh = jg + jstepl
                     if (jh.lt.istart) jh = jh + ninc
                     ji = jh + jstepl
                     if (ji.lt.istart) ji = ji + ninc
                     j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
                     do 410 l = 1 , nvex
                        ajb = a(jb+j)
                        ajc = a(jc+j)
                        t1 = ajb + ajc
                        aja = a(ja+j)
                        t2 = aja - 0.5d0 * t1
                        t3 = c1 * ( ajb - ajc )
                        ajd = a(jd+j)
                        ajb =  ajd
                        bjb = b(jb+j)
                        bjc = b(jc+j)
                        u1 = bjb + bjc
                        bja = b(ja+j)
                        u2 = bja - 0.5d0 * u1
                        u3 = c1 * ( bjb - bjc )
                        bjd = b(jd+j)
                        bjb =  bjd
                        a(ja+j) = aja + t1
                        b(ja+j) = bja + u1
                        a(jd+j) = t2 - u3
                        b(jd+j) = u2 + t3
                        ajc =  t2 + u3
                        bjc =  u2 - t3
C----------------------
                        aje = a(je+j)
                        ajf = a(jf+j)
                        t1 = aje + ajf
                        t2 = ajb - 0.5d0 * t1
                        t3 = c1 * ( aje - ajf )
                        ajh = a(jh+j)
                        ajf =  ajh
                        bje = b(je+j)
                        bjf = b(jf+j)
                        u1 = bje + bjf
                        u2 = bjb - 0.5d0 * u1
                        u3 = c1 * ( bje - bjf )
                        bjh = b(jh+j)
                        bjf =  bjh
                        a(jb+j) = ajb + t1
                        b(jb+j) = bjb + u1
                        a(je+j) = t2 - u3
                        b(je+j) = u2 + t3
                        a(jh+j) = t2 + u3
                        b(jh+j) = u2 - t3
C----------------------
                        aji = a(ji+j)
                        t1 = ajf + aji
                        ajg = a(jg+j)
                        t2 = ajg - 0.5d0 * t1
                        t3 = c1 * ( ajf - aji )
                        t1 = ajg + t1
                        a(jg+j) = ajc
                        bji = b(ji+j)
                        u1 = bjf + bji
                        bjg = b(jg+j)
                        u2 = bjg - 0.5d0 * u1
                        u3 = c1 * ( bjf - bji )
                        u1 = bjg + u1
                        b(jg+j) = bjc
                        a(jc+j) = t1
                        b(jc+j) = u1
                        a(jf+j) = t2 - u3
                        b(jf+j) = u2 + t3
                        a(ji+j) = t2 + u3
                        b(ji+j) = u2 - t3
                        j = j + jump
 410                 continue
C-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
 415              continue
 420           continue
 430        continue
C-----( end of double loop for k=0 )
C
C     finished if last pass
C     ---------------------
            if (ipass.eq.m) go to 490
C
            kk = 2*la
C
C     loop on nonzero k
C     -----------------
            do 470 k = ink , jstep-ink , ink
               co1 = trigs(kk+1)
               si1 = s*trigs(kk+2)
               co2 = trigs(2*kk+1)
               si2 = s*trigs(2*kk+2)
C
C     double loop along first transform in block
C     ------------------------------------------
               do 460 ll = k , (la-1)*ink , 3*jstep
C
                  do 450 jjj = ll , (n-1)*inc , 3*la*ink
                     ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                     do 445 nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        jd = ja + laincl
                        if (jd.lt.istart) jd = jd + ninc
                        je = jd + jstepl
                        if (je.lt.istart) je = je + ninc
                        jf = je + jstepl
                        if (jf.lt.istart) jf = jf + ninc
                        jg = jd + laincl
                        if (jg.lt.istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh.lt.istart) jh = jh + ninc
                        ji = jh + jstepl
                        if (ji.lt.istart) ji = ji + ninc
                        j = 0
C
C     loop across transforms
C     ----------------------
c     dir$ ivdep, shortloop
                        do 440 l = 1 , nvex
                           ajb = a(jb+j)
                           ajc = a(jc+j)
                           t1 = ajb + ajc
                           aja = a(ja+j)
                           t2 = aja - 0.5d0 * t1
                           t3 = c1 * ( ajb - ajc )
                           ajd = a(jd+j)
                           ajb =  ajd
                           bjb = b(jb+j)
                           bjc = b(jc+j)
                           u1 = bjb + bjc
                           bja = b(ja+j)
                           u2 = bja - 0.5d0 * u1
                           u3 = c1 * ( bjb - bjc )
                           bjd = b(jd+j)
                           bjb =  bjd
                           a(ja+j) = aja + t1
                           b(ja+j) = bja + u1
                           a(jd+j) = co1*(t2-u3) - si1*(u2+t3)
                           b(jd+j) = si1*(t2-u3) + co1*(u2+t3)
                           ajc =  co2*(t2+u3) - si2*(u2-t3)
                           bjc =  si2*(t2+u3) + co2*(u2-t3)
C----------------------
                           aje = a(je+j)
                           ajf = a(jf+j)
                           t1 = aje + ajf
                           t2 = ajb - 0.5d0 * t1
                           t3 = c1 * ( aje - ajf )
                           ajh = a(jh+j)
                           ajf =  ajh
                           bje = b(je+j)
                           bjf = b(jf+j)
                           u1 = bje + bjf
                           u2 = bjb - 0.5d0 * u1
                           u3 = c1 * ( bje - bjf )
                           bjh = b(jh+j)
                           bjf =  bjh
                           a(jb+j) = ajb + t1
                           b(jb+j) = bjb + u1
                           a(je+j) = co1*(t2-u3) - si1*(u2+t3)
                           b(je+j) = si1*(t2-u3) + co1*(u2+t3)
                           a(jh+j) = co2*(t2+u3) - si2*(u2-t3)
                           b(jh+j) = si2*(t2+u3) + co2*(u2-t3)
C----------------------
                           aji = a(ji+j)
                           t1 = ajf + aji
                           ajg = a(jg+j)
                           t2 = ajg - 0.5d0 * t1
                           t3 = c1 * ( ajf - aji )
                           t1 = ajg + t1
                           a(jg+j) = ajc
                           bji = b(ji+j)
                           u1 = bjf + bji
                           bjg = b(jg+j)
                           u2 = bjg - 0.5d0 * u1
                           u3 = c1 * ( bjf - bji )
                           u1 = bjg + u1
                           b(jg+j) = bjc
                           a(jc+j) = t1
                           b(jc+j) = u1
                           a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                           b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                           a(ji+j) = co2*(t2+u3) - si2*(u2-t3)
                           b(ji+j) = si2*(t2+u3) + co2*(u2-t3)
                           j = j + jump
 440                    continue
C-----(end of loop across transforms)
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
 445                 continue
 450              continue
 460           continue
C-----( end of double loop for this k )
               kk = kk + 2*la
 470        continue
C-----( end of loop over values of k )
            la = 3*la
 480     continue
C-----( end of loop on type II radix-3 passes )
C-----( nvex transforms completed)
 490     continue
         istart = istart + nvex * jump
 500  continue
C-----( end of loop on blocks of transforms )
C
      return
      end
!=======================================================================
C     fortran version of *gpfa5* -
C     radix-5 section of self-sorting, in-place,
C     generalized pfa
C
C     **********************************************************************
C     **********************************************************************
C     fortran version of *gpfa5* -
C     radix-5 section of self-sorting, in-place,
C     generalized pfa
C
C-------------------------------------------------------------------
C
      subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)

      IMPLICIT NONE

C     Incoming data
      INTEGER inc, jump, n, mm, lot, isign
      REAL*8 trigs(*)

C     In-Outcoming data
      REAL*8 a(*), b(*)

C     Internal data
      INTEGER lvr/128/, n5, inq, jstepx, ink, mu, m, mh, nblox,
     &     left, istart, nb, nvex, la, ipass, jstep, jstepl,
     &     kk, k, jjj, j, ja, jb, jc, jd, je, jf, jg, jh, ji,
     &     jj, jk, jl, jm, jn, jo, jp, jq, jr, js, jt, ju, jv,
     &     jw, jx, jy, nu, l, laincl,
     &     ll, ninc


      REAL*8  sin36/0.587785252292473/, sin72/0.951056516295154/,
     &     qrt5/0.559016994374947/, s, c1, c2, c3,
     &     co1, co2, co3, co4, si1, si2, si3, si4, ax, bx,
     &     aja, ajb, ajc, ajd, aje, ajf, ajg, ajh, aji, ajj, ajk,
     &     ajl, ajm, ajn, ajo, ajp, ajq, ajr, ajs, ajt, aju, ajv,
     &     ajw, ajx, ajy,
     &     bja, bjb, bjc, bjd, bje, bjf, bjg, bjh, bji, bjj, bjk,
     &     bjl, bjm, bjn, bjo, bjp, bjq, bjr, bjs, bjt, bju, bjv,
     &     bjw, bjx, bjy,
     &     t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11,
     &     u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11
C     ***************************************************************
C     *                                                             *
C     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
C     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
C     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
C     *                                                             *
C     ***************************************************************
C
      n5 = 5 ** mm
      inq = n / n5
      jstepx = (n5-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,5)
      if (isign.eq.-1) mu = 5 - mu
C
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = qrt5
      c2 = sin72
      c3 = sin36
      if (mu.eq.2.or.mu.eq.3) then
         c1 = -c1
         c2 = sin36
         c3 = sin72
      endif
      if (mu.eq.3.or.mu.eq.4) c2 = -c2
      if (mu.eq.2.or.mu.eq.4) c3 = -c3
C
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = float(isign)
      istart = 1
C
C     loop on blocks of lvr transforms
C     --------------------------------
      do 500 nb = 1 , nblox
C
         if (left.le.lvr) then
            nvex = left
         else if (left.lt.(2*lvr)) then
            nvex = left/2
            nvex = nvex + mod(nvex,2)
         else
            nvex = lvr
         endif
         left = left - nvex
C
         la = 1
C
C     loop on type I radix-5 passes
C     -----------------------------
         do 160 ipass = 1 , mh
            jstep = (n*inc) / (5*la)
            jstepl = jstep - ninc
            kk = 0
C
C     loop on k
C     ---------
            do 150 k = 0 , jstep-ink , ink
C
               if (k.gt.0) then
                  co1 = trigs(kk+1)
                  si1 = s*trigs(kk+2)
                  co2 = trigs(2*kk+1)
                  si2 = s*trigs(2*kk+2)
                  co3 = trigs(3*kk+1)
                  si3 = s*trigs(3*kk+2)
                  co4 = trigs(4*kk+1)
                  si4 = s*trigs(4*kk+2)
               endif
C
C     loop along transform
C     --------------------
               do 140 jjj = k , (n-1)*inc , 5*jstep
                  ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                  do 135 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd.lt.istart) jd = jd + ninc
                     je = jd + jstepl
                     if (je.lt.istart) je = je + ninc
                     j = 0
C
C     loop across transforms
C     ----------------------
                     if (k.eq.0) then
C
c     dir$ ivdep, shortloop
                        do 110 l = 1 , nvex
                           ajb = a(jb+j)
                           aje = a(je+j)
                           t1 = ajb + aje
                           ajc = a(jc+j)
                           ajd = a(jd+j)
                           t2 = ajc + ajd
                           t3 = ajb - aje
                           t4 = ajc - ajd
                           t5 = t1 + t2
                           t6 = c1 * ( t1 - t2 )
                           aja = a(ja+j)
                           t7 = aja - 0.25d0 * t5
                           a(ja+j) = aja + t5
                           t8 = t7 + t6
                           t9 = t7 - t6
                           t10 = c3 * t3 - c2 * t4
                           t11 = c2 * t3 + c3 * t4
                           bjb = b(jb+j)
                           bje = b(je+j)
                           u1 = bjb + bje
                           bjc = b(jc+j)
                           bjd = b(jd+j)
                           u2 = bjc + bjd
                           u3 = bjb - bje
                           u4 = bjc - bjd
                           u5 = u1 + u2
                           u6 = c1 * ( u1 - u2 )
                           bja = b(ja+j)
                           u7 = bja - 0.25d0 * u5
                           b(ja+j) = bja + u5
                           u8 = u7 + u6
                           u9 = u7 - u6
                           u10 = c3 * u3 - c2 * u4
                           u11 = c2 * u3 + c3 * u4
                           a(jb+j) = t8 - u11
                           b(jb+j) = u8 + t11
                           a(je+j) = t8 + u11
                           b(je+j) = u8 - t11
                           a(jc+j) = t9 - u10
                           b(jc+j) = u9 + t10
                           a(jd+j) = t9 + u10
                           b(jd+j) = u9 - t10
                           j = j + jump
 110                    continue
C
                     else
C
c     dir$ ivdep,shortloop
                        do 130 l = 1 , nvex
                           ajb = a(jb+j)
                           aje = a(je+j)
                           t1 = ajb + aje
                           ajc = a(jc+j)
                           ajd = a(jd+j)
                           t2 = ajc + ajd
                           t3 = ajb - aje
                           t4 = ajc - ajd
                           t5 = t1 + t2
                           t6 = c1 * ( t1 - t2 )
                           aja = a(ja+j)
                           t7 = aja - 0.25d0 * t5
                           a(ja+j) = aja + t5
                           t8 = t7 + t6
                           t9 = t7 - t6
                           t10 = c3 * t3 - c2 * t4
                           t11 = c2 * t3 + c3 * t4
                           bjb = b(jb+j)
                           bje = b(je+j)
                           u1 = bjb + bje
                           bjc = b(jc+j)
                           bjd = b(jd+j)
                           u2 = bjc + bjd
                           u3 = bjb - bje
                           u4 = bjc - bjd
                           u5 = u1 + u2
                           u6 = c1 * ( u1 - u2 )
                           bja = b(ja+j)
                           u7 = bja - 0.25d0 * u5
                           b(ja+j) = bja + u5
                           u8 = u7 + u6
                           u9 = u7 - u6
                           u10 = c3 * u3 - c2 * u4
                           u11 = c2 * u3 + c3 * u4
                           a(jb+j) = co1*(t8-u11) - si1*(u8+t11)
                           b(jb+j) = si1*(t8-u11) + co1*(u8+t11)
                           a(je+j) = co4*(t8+u11) - si4*(u8-t11)
                           b(je+j) = si4*(t8+u11) + co4*(u8-t11)
                           a(jc+j) = co2*(t9-u10) - si2*(u9+t10)
                           b(jc+j) = si2*(t9-u10) + co2*(u9+t10)
                           a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
                           b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
                           j = j + jump
 130                    continue
C
                     endif
C
C-----( end of loop across transforms )
C
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
 135              continue
 140           continue
C-----( end of loop along transforms )
               kk = kk + 2*la
 150        continue
C-----( end of loop on nonzero k )
            la = 5*la
 160     continue
C-----( end of loop on type I radix-5 passes)
C
         if (n.eq.5) go to 490
C
C     loop on type II radix-5 passes
C     ------------------------------
 400     continue
C
         do 480 ipass = mh+1 , m
            jstep = (n*inc) / (5*la)
            jstepl = jstep - ninc
            laincl = la * ink - ninc
            kk = 0
C
C     loop on k
C     ---------
            do 470 k = 0 , jstep-ink , ink
C
               if (k.gt.0) then
                  co1 = trigs(kk+1)
                  si1 = s*trigs(kk+2)
                  co2 = trigs(2*kk+1)
                  si2 = s*trigs(2*kk+2)
                  co3 = trigs(3*kk+1)
                  si3 = s*trigs(3*kk+2)
                  co4 = trigs(4*kk+1)
                  si4 = s*trigs(4*kk+2)
               endif
C
C     double loop along first transform in block
C     ------------------------------------------
               do 460 ll = k , (la-1)*ink , 5*jstep
C
                  do 450 jjj = ll , (n-1)*inc , 5*la*ink
                     ja = istart + jjj
C
C     "transverse" loop
C     -----------------
                     do 445 nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        jd = jc + jstepl
                        if (jd.lt.istart) jd = jd + ninc
                        je = jd + jstepl
                        if (je.lt.istart) je = je + ninc
                        jf = ja + laincl
                        if (jf.lt.istart) jf = jf + ninc
                        jg = jf + jstepl
                        if (jg.lt.istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh.lt.istart) jh = jh + ninc
                        ji = jh + jstepl
                        if (ji.lt.istart) ji = ji + ninc
                        jj = ji + jstepl
                        if (jj.lt.istart) jj = jj + ninc
                        jk = jf + laincl
                        if (jk.lt.istart) jk = jk + ninc
                        jl = jk + jstepl
                        if (jl.lt.istart) jl = jl + ninc
                        jm = jl + jstepl
                        if (jm.lt.istart) jm = jm + ninc
                        jn = jm + jstepl
                        if (jn.lt.istart) jn = jn + ninc
                        jo = jn + jstepl
                        if (jo.lt.istart) jo = jo + ninc
                        jp = jk + laincl
                        if (jp.lt.istart) jp = jp + ninc
                        jq = jp + jstepl
                        if (jq.lt.istart) jq = jq + ninc
                        jr = jq + jstepl
                        if (jr.lt.istart) jr = jr + ninc
                        js = jr + jstepl
                        if (js.lt.istart) js = js + ninc
                        jt = js + jstepl
                        if (jt.lt.istart) jt = jt + ninc
                        ju = jp + laincl
                        if (ju.lt.istart) ju = ju + ninc
                        jv = ju + jstepl
                        if (jv.lt.istart) jv = jv + ninc
                        jw = jv + jstepl
                        if (jw.lt.istart) jw = jw + ninc
                        jx = jw + jstepl
                        if (jx.lt.istart) jx = jx + ninc
                        jy = jx + jstepl
                        if (jy.lt.istart) jy = jy + ninc
                        j = 0
C
C     loop across transforms
C     ----------------------
                        if (k.eq.0) then
C
c     dir$ ivdep, shortloop
                           do 410 l = 1 , nvex
                              ajb = a(jb+j)
                              aje = a(je+j)
                              t1 = ajb + aje
                              ajc = a(jc+j)
                              ajd = a(jd+j)
                              t2 = ajc + ajd
                              t3 = ajb - aje
                              t4 = ajc - ajd
                              ajf = a(jf+j)
                              ajb =  ajf
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              aja = a(ja+j)
                              t7 = aja - 0.25d0 * t5
                              a(ja+j) = aja + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              ajk = a(jk+j)
                              ajc =  ajk
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              bjb = b(jb+j)
                              bje = b(je+j)
                              u1 = bjb + bje
                              bjc = b(jc+j)
                              bjd = b(jd+j)
                              u2 = bjc + bjd
                              u3 = bjb - bje
                              u4 = bjc - bjd
                              bjf = b(jf+j)
                              bjb =  bjf
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              bja = b(ja+j)
                              u7 = bja - 0.25d0 * u5
                              b(ja+j) = bja + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              bjk = b(jk+j)
                              bjc =  bjk
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jf+j) = t8 - u11
                              b(jf+j) = u8 + t11
                              aje =  t8 + u11
                              bje =  u8 - t11
                              a(jk+j) = t9 - u10
                              b(jk+j) = u9 + t10
                              ajd =  t9 + u10
                              bjd =  u9 - t10
C----------------------
                              ajg = a(jg+j)
                              ajj = a(jj+j)
                              t1 = ajg + ajj
                              ajh = a(jh+j)
                              aji = a(ji+j)
                              t2 = ajh + aji
                              t3 = ajg - ajj
                              t4 = ajh - aji
                              ajl = a(jl+j)
                              ajh =  ajl
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              t7 = ajb - 0.25d0 * t5
                              a(jb+j) = ajb + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              ajq = a(jq+j)
                              aji =  ajq
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              bjg = b(jg+j)
                              bjj = b(jj+j)
                              u1 = bjg + bjj
                              bjh = b(jh+j)
                              bji = b(ji+j)
                              u2 = bjh + bji
                              u3 = bjg - bjj
                              u4 = bjh - bji
                              bjl = b(jl+j)
                              bjh =  bjl
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              u7 = bjb - 0.25d0 * u5
                              b(jb+j) = bjb + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              bjq = b(jq+j)
                              bji =  bjq
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jg+j) = t8 - u11
                              b(jg+j) = u8 + t11
                              ajj =  t8 + u11
                              bjj =  u8 - t11
                              a(jl+j) = t9 - u10
                              b(jl+j) = u9 + t10
                              a(jq+j) = t9 + u10
                              b(jq+j) = u9 - t10
C----------------------
                              ajo = a(jo+j)
                              t1 = ajh + ajo
                              ajm = a(jm+j)
                              ajn = a(jn+j)
                              t2 = ajm + ajn
                              t3 = ajh - ajo
                              t4 = ajm - ajn
                              ajr = a(jr+j)
                              ajn =  ajr
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              t7 = ajc - 0.25d0 * t5
                              a(jc+j) = ajc + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              ajw = a(jw+j)
                              ajo =  ajw
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              bjo = b(jo+j)
                              u1 = bjh + bjo
                              bjm = b(jm+j)
                              bjn = b(jn+j)
                              u2 = bjm + bjn
                              u3 = bjh - bjo
                              u4 = bjm - bjn
                              bjr = b(jr+j)
                              bjn =  bjr
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              u7 = bjc - 0.25d0 * u5
                              b(jc+j) = bjc + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              bjw = b(jw+j)
                              bjo =  bjw
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jh+j) = t8 - u11
                              b(jh+j) = u8 + t11
                              a(jw+j) = t8 + u11
                              b(jw+j) = u8 - t11
                              a(jm+j) = t9 - u10
                              b(jm+j) = u9 + t10
                              a(jr+j) = t9 + u10
                              b(jr+j) = u9 - t10
C----------------------
                              ajt = a(jt+j)
                              t1 = aji + ajt
                              ajs = a(js+j)
                              t2 = ajn + ajs
                              t3 = aji - ajt
                              t4 = ajn - ajs
                              ajx = a(jx+j)
                              ajt =  ajx
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              ajp = a(jp+j)
                              t7 = ajp - 0.25d0 * t5
                              ax = ajp + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jp+j) = ajd
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(jd+j) = ax
                              bjt = b(jt+j)
                              u1 = bji + bjt
                              bjs = b(js+j)
                              u2 = bjn + bjs
                              u3 = bji - bjt
                              u4 = bjn - bjs
                              bjx = b(jx+j)
                              bjt =  bjx
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              bjp = b(jp+j)
                              u7 = bjp - 0.25d0 * u5
                              bx = bjp + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jp+j) = bjd
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(jd+j) = bx
                              a(ji+j) = t8 - u11
                              b(ji+j) = u8 + t11
                              a(jx+j) = t8 + u11
                              b(jx+j) = u8 - t11
                              a(jn+j) = t9 - u10
                              b(jn+j) = u9 + t10
                              a(js+j) = t9 + u10
                              b(js+j) = u9 - t10
C----------------------
                              ajv = a(jv+j)
                              ajy = a(jy+j)
                              t1 = ajv + ajy
                              t2 = ajo + ajt
                              t3 = ajv - ajy
                              t4 = ajo - ajt
                              a(jv+j) = ajj
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              aju = a(ju+j)
                              t7 = aju - 0.25d0 * t5
                              ax = aju + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(ju+j) = aje
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(je+j) = ax
                              bjv = b(jv+j)
                              bjy = b(jy+j)
                              u1 = bjv + bjy
                              u2 = bjo + bjt
                              u3 = bjv - bjy
                              u4 = bjo - bjt
                              b(jv+j) = bjj
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              bju = b(ju+j)
                              u7 = bju - 0.25d0 * u5
                              bx = bju + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(ju+j) = bje
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(je+j) = bx
                              a(jj+j) = t8 - u11
                              b(jj+j) = u8 + t11
                              a(jy+j) = t8 + u11
                              b(jy+j) = u8 - t11
                              a(jo+j) = t9 - u10
                              b(jo+j) = u9 + t10
                              a(jt+j) = t9 + u10
                              b(jt+j) = u9 - t10
                              j = j + jump
 410                       continue
C
                        else
C
c     dir$ ivdep, shortloop
                           do 440 l = 1 , nvex
                              ajb = a(jb+j)
                              aje = a(je+j)
                              t1 = ajb + aje
                              ajc = a(jc+j)
                              ajd = a(jd+j)
                              t2 = ajc + ajd
                              t3 = ajb - aje
                              t4 = ajc - ajd
                              ajf = a(jf+j)
                              ajb =  ajf
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              aja = a(ja+j)
                              t7 = aja - 0.25d0 * t5
                              a(ja+j) = aja + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              ajk = a(jk+j)
                              ajc =  ajk
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              bjb = b(jb+j)
                              bje = b(je+j)
                              u1 = bjb + bje
                              bjc = b(jc+j)
                              bjd = b(jd+j)
                              u2 = bjc + bjd
                              u3 = bjb - bje
                              u4 = bjc - bjd
                              bjf = b(jf+j)
                              bjb =  bjf
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              bja = b(ja+j)
                              u7 = bja - 0.25d0 * u5
                              b(ja+j) = bja + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              bjk = b(jk+j)
                              bjc =  bjk
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jf+j) = co1*(t8-u11) - si1*(u8+t11)
                              b(jf+j) = si1*(t8-u11) + co1*(u8+t11)
                              aje =  co4*(t8+u11) - si4*(u8-t11)
                              bje =  si4*(t8+u11) + co4*(u8-t11)
                              a(jk+j) = co2*(t9-u10) - si2*(u9+t10)
                              b(jk+j) = si2*(t9-u10) + co2*(u9+t10)
                              ajd =  co3*(t9+u10) - si3*(u9-t10)
                              bjd =  si3*(t9+u10) + co3*(u9-t10)
C----------------------
                              ajg = a(jg+j)
                              ajj = a(jj+j)
                              t1 = ajg + ajj
                              ajh = a(jh+j)
                              aji = a(ji+j)
                              t2 = ajh + aji
                              t3 = ajg - ajj
                              t4 = ajh - aji
                              ajl = a(jl+j)
                              ajh =  ajl
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              t7 = ajb - 0.25d0 * t5
                              a(jb+j) = ajb + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              ajq = a(jq+j)
                              aji =  ajq
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              bjg = b(jg+j)
                              bjj = b(jj+j)
                              u1 = bjg + bjj
                              bjh = b(jh+j)
                              bji = b(ji+j)
                              u2 = bjh + bji
                              u3 = bjg - bjj
                              u4 = bjh - bji
                              bjl = b(jl+j)
                              bjh =  bjl
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              u7 = bjb - 0.25d0 * u5
                              b(jb+j) = bjb + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              bjq = b(jq+j)
                              bji =  bjq
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jg+j) = co1*(t8-u11) - si1*(u8+t11)
                              b(jg+j) = si1*(t8-u11) + co1*(u8+t11)
                              ajj =  co4*(t8+u11) - si4*(u8-t11)
                              bjj =  si4*(t8+u11) + co4*(u8-t11)
                              a(jl+j) = co2*(t9-u10) - si2*(u9+t10)
                              b(jl+j) = si2*(t9-u10) + co2*(u9+t10)
                              a(jq+j) = co3*(t9+u10) - si3*(u9-t10)
                              b(jq+j) = si3*(t9+u10) + co3*(u9-t10)
C----------------------
                              ajo = a(jo+j)
                              t1 = ajh + ajo
                              ajm = a(jm+j)
                              ajn = a(jn+j)
                              t2 = ajm + ajn
                              t3 = ajh - ajo
                              t4 = ajm - ajn
                              ajr = a(jr+j)
                              ajn =  ajr
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              t7 = ajc - 0.25d0 * t5
                              a(jc+j) = ajc + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              ajw = a(jw+j)
                              ajo =  ajw
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              bjo = b(jo+j)
                              u1 = bjh + bjo
                              bjm = b(jm+j)
                              bjn = b(jn+j)
                              u2 = bjm + bjn
                              u3 = bjh - bjo
                              u4 = bjm - bjn
                              bjr = b(jr+j)
                              bjn =  bjr
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              u7 = bjc - 0.25d0 * u5
                              b(jc+j) = bjc + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              bjw = b(jw+j)
                              bjo =  bjw
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jh+j) = co1*(t8-u11) - si1*(u8+t11)
                              b(jh+j) = si1*(t8-u11) + co1*(u8+t11)
                              a(jw+j) = co4*(t8+u11) - si4*(u8-t11)
                              b(jw+j) = si4*(t8+u11) + co4*(u8-t11)
                              a(jm+j) = co2*(t9-u10) - si2*(u9+t10)
                              b(jm+j) = si2*(t9-u10) + co2*(u9+t10)
                              a(jr+j) = co3*(t9+u10) - si3*(u9-t10)
                              b(jr+j) = si3*(t9+u10) + co3*(u9-t10)
C----------------------
                              ajt = a(jt+j)
                              t1 = aji + ajt
                              ajs = a(js+j)
                              t2 = ajn + ajs
                              t3 = aji - ajt
                              t4 = ajn - ajs
                              ajx = a(jx+j)
                              ajt =  ajx
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              ajp = a(jp+j)
                              t7 = ajp - 0.25d0 * t5
                              ax = ajp + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jp+j) = ajd
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(jd+j) = ax
                              bjt = b(jt+j)
                              u1 = bji + bjt
                              bjs = b(js+j)
                              u2 = bjn + bjs
                              u3 = bji - bjt
                              u4 = bjn - bjs
                              bjx = b(jx+j)
                              bjt =  bjx
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              bjp = b(jp+j)
                              u7 = bjp - 0.25d0 * u5
                              bx = bjp + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jp+j) = bjd
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(jd+j) = bx
                              a(ji+j) = co1*(t8-u11) - si1*(u8+t11)
                              b(ji+j) = si1*(t8-u11) + co1*(u8+t11)
                              a(jx+j) = co4*(t8+u11) - si4*(u8-t11)
                              b(jx+j) = si4*(t8+u11) + co4*(u8-t11)
                              a(jn+j) = co2*(t9-u10) - si2*(u9+t10)
                              b(jn+j) = si2*(t9-u10) + co2*(u9+t10)
                              a(js+j) = co3*(t9+u10) - si3*(u9-t10)
                              b(js+j) = si3*(t9+u10) + co3*(u9-t10)
C----------------------
                              ajv = a(jv+j)
                              ajy = a(jy+j)
                              t1 = ajv + ajy
                              t2 = ajo + ajt
                              t3 = ajv - ajy
                              t4 = ajo - ajt
                              a(jv+j) = ajj
                              t5 = t1 + t2
                              t6 = c1 * ( t1 - t2 )
                              aju = a(ju+j)
                              t7 = aju - 0.25d0 * t5
                              ax = aju + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(ju+j) = aje
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(je+j) = ax
                              bjv = b(jv+j)
                              bjy = b(jy+j)
                              u1 = bjv + bjy
                              u2 = bjo + bjt
                              u3 = bjv - bjy
                              u4 = bjo - bjt
                              b(jv+j) = bjj
                              u5 = u1 + u2
                              u6 = c1 * ( u1 - u2 )
                              bju = b(ju+j)
                              u7 = bju - 0.25d0 * u5
                              bx = bju + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(ju+j) = bje
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(je+j) = bx
                              a(jj+j) = co1*(t8-u11) - si1*(u8+t11)
                              b(jj+j) = si1*(t8-u11) + co1*(u8+t11)
                              a(jy+j) = co4*(t8+u11) - si4*(u8-t11)
                              b(jy+j) = si4*(t8+u11) + co4*(u8-t11)
                              a(jo+j) = co2*(t9-u10) - si2*(u9+t10)
                              b(jo+j) = si2*(t9-u10) + co2*(u9+t10)
                              a(jt+j) = co3*(t9+u10) - si3*(u9-t10)
                              b(jt+j) = si3*(t9+u10) + co3*(u9-t10)
                              j = j + jump
 440                       continue
C
                        endif
C
C-----(end of loop across transforms)
C
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
 445                 continue
 450              continue
 460           continue
C-----( end of double loop for this k )
               kk = kk + 2*la
 470        continue
C-----( end of loop over values of k )
            la = 5*la
 480     continue
C-----( end of loop on type II radix-5 passes )
C-----( nvex transforms completed)
 490     continue
         istart = istart + nvex * jump
 500  continue
C-----( end of loop on blocks of transforms )
C
      return
      end
