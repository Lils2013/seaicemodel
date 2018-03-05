      SUBROUTINE
     * dzcal4(dz,f,diag,km2,adz,il,jl,il1,jl1,ilp,jlp, maxite,
     *        len,krylov,krp1,b,sl,r0,s0,c0,giv,w0,z0,v0)

      REAL*8   B(LEN), SL(LEN), R0(LEN)
      REAL*8   S0(LEN), C0(KRP1,KRP1), GIV(2,KRP1),W0(LEN)
      REAL*8   Z0(LEN,KRP1),V0(LEN)
      dimension dz(-1:ilp,-1:jlp), f(il,jl), diag(19,il,jl),
     *          km2(0:il1,0:jl1), adz(-1:ilp,-1:jlp)
      real*8 dz, f, diag, adz

*     Переписываем f -> adz
      do j=-1,jlp
      do i=-1,ilp
      adz(i,j)=0.d0
      end do
      end do
      do j=1,jl
      do i=1,il
      adz(i,j)=f(i,j)
      end do
      end do

*     Переписываем вектор dz в sl(len) и f(i,j) в b(len).
*     Используется диагональный предобуславливатель.

      !write(*,*) "adz", adz
      call prdiag(diag,km2,adz,il,jl,ilp,jlp,il1,jl1)
      

      m= 0
      do j=1,jl
      do i=1,il
      if(km2(i,j) .GT. 0) then
      m= m+1
      sl(m)= dz (i,j)
cc      sl(m)= 0.d0
      b (m)= adz(i,j)
      end if
      end do
      end do
      if( m .NE. len) write(*,*)'Problems with LENа...'
      
      !write(*,*) "V0", V0
      !write(*,*) "DZBEFORE", dz
      !write(*,*) "diag", diag
      !write(*,*) "b", b
      !write(*,*) "adz", adz
      !write(*,*) "R0", R0
      !write(*,*) "S0", S0
      !write(*,*) "C0", C0
      !write(*,*) "GIV", GIV
      !write(*,*) "W0", W0
      !write(*,*) "Z0", Z0
      
      call gmres(b,sl,dz,km2,diag,adz,il,jl,il1,jl1,ilp,jlp,maxite
     &           ,LEN,Krylov,KRP1,R0,S0,C0,GIV,W0,V0,Z0)
      !write(*,*) "DZAFTER", dz
      !write(*,*) "V0", V0
      
*     Переписываем результат в двумерный массив dz.
      m= 0
      
      do j=1,jl
      do i=1,il
      if (km2(i,j) .GT. 0) then
      m= m+1
      dz(i,j)= sl(m)
      end if
      end do
      end do

      return
      end

      SUBROUTINE GMRES(B,X,dz,km2,diag,adz,il,jl,il1,jl1,ilp,jlp,
     &                 MAXITE,LEN,K,KRP1,R,S,C,GIV,W,V,Z)
*     LEN is a length of the seeked vector with 1D indexing.
*     K is a Krylov space dimension. One of the most critical parms.
*     B is a right hand side, X is the seeked vector.
*     MAXITE is a max number of restarts of GMRES.

      dimension km2(0:il1,0:jl1), diag(19,il,jl)

*     Some service arrays.
      dimension dz(-1:ilp,-1:jlp), adz(-1:ilp,-1:jlp)

      REAL*8   B(LEN), X(LEN), R(LEN)
      REAL*8   RESI, RESIN, SUM
      REAL*8   S(LEN), C(KRP1,KRP1), GIV(2,KRP1),W(LEN)
      REAL*8   Z(LEN,KRP1),  V(LEN)
      real*8 dz, adz, diag, err

      err=1.d-19

      C = 0.0D0
      Z = 0.d0

      DO 100 ITER = 1,MAXITE

*     Rewriting vector X to vector R.
      R= X

*     Multiplication R by Matrix. The result is in the R vector.
*     DZ, ADZ are used as a service array.
      CALL PROD( diag,km2,R,dz,adz,il,jl,il1,jl1,ilp,jlp,len)

         RESI    = 0.0D0

         DO 101 I = 1,LEN
            R(I) = B(I) - R(I)
            RESI = RESI + R(I)**2
101      CONTINUE

         SUM     = DSQRT(RESI)

*     Controlling the accuracy and the formation of the criterion to
*     stop the calculation.
cc         WRITE(*,*) ' ABSERR :',SUM

        IF (sum .LT. err) GOTO 200

         W(1)    = SUM
         RESI    = 2.0D0 * SUM *(SUM -  R(1))
         RESI    = DSQRT(RESI)
         RESIN=1.D0/RESI
         Z(1,1)  = (R(1) - SUM)*RESIN

         DO 105 I = 2,LEN
            Z(I,1) = R(I) / RESI
            W(I)   = 0.0D0
105      CONTINUE

         KK = 0

         DO 107  I  = 1,K
            DO 108 J = 1,LEN
               R(J) = 0.d0
108         CONTINUE
            R(I)= 1.d0
            DO 109 J = I,1,-1
               SUM  = 0.0D0
               DO 110 M = 1,LEN
                  SUM  = SUM + Z(M,J)*R(M)
110            CONTINUE
               DO 111 M = 1,LEN
                  R(M) = R(M) - 2.0D0 * SUM* Z(M,J)
111            CONTINUE
109         CONTINUE

*     Vector R is rewritten in V.
      V= R

*     Vector V is multiplied by the Matrix. Result is in the V.
*     DZ, ADZ are used as a service array.
            CALL PROD(diag,km2,V,dz,adz,il,jl,il1,jl1,ilp,jlp,len)

            DO 115 J= 1,I
               SUM= 0.0D0
               DO 116 M = 1,LEN
                  SUM   = SUM + Z(M,J) * V(M)
116            CONTINUE
               DO 117 M = 1,LEN
                  V(M)  = V(M) - 2.0D0 * SUM * Z(M,J)
117            CONTINUE
115         CONTINUE

            SUM      = 0.0D0
            DO 118 J = I+1,LEN
               SUM   = SUM + V(J) ** 2
118         CONTINUE

            IF (SUM .NE. 0.0D0) THEN
               SUM         = DSQRT(SUM)
               RESI        = 2.0D0 * SUM ** 2 - 2.0D0 * SUM * V(I+1)
               RESI        = DSQRT(RESI)
               RESIN=1.D0/RESI
               Z(I+1,I+1)  = (V(I+1) - SUM) / RESI
               DO  119  J  = I+2,LEN
                  Z(J,I+1) = V(J) / RESI
                  V(J)     = 0.d0
119            CONTINUE
               V(I+1)      = SUM
            END IF


            IF (I .GT. 1) THEN
               DO 121  J = 1,I-1
                  RESI   = V(J)
                  SUM    = V(J+1)
                  V(J)   = RESI * GIV(1,J) - SUM * GIV(2,J)
                  V(J+1) = RESI * GIV(2,J) + SUM * GIV(1,J)
121            CONTINUE
            END IF

            IF (V(I+1) .EQ. 0.d0) THEN
               GIV(1,I)  = 0.0D0
               GIV(2,I)  = 0.0D0
            ELSE
               RESI      =   DSQRT(V(I) ** 2 + V(I+1) ** 2)
               RESIN=1.D0/RESI
               GIV(1,I)  =   V(I)  * RESIN
               GIV(2,I)  = - V(I+1)* RESIN
               V(I)      = V(I) * GIV(1,I) - V(I+1) * GIV(2,I)
               V(I+1)    = 0.d0
               RESI      = W(I)
               SUM       = W(I+1)
               W(I)      = RESI * GIV(1,I) - SUM * GIV(2,I)
               W(I+1)    = RESI * GIV(2,I) + SUM * GIV(1,I)
            END IF

            DO 122  J    = 1,I
               C(J,I)    = V(J)
122         CONTINUE
            KK           = I

C     Checking the error in the inner loop.
            RESI         = W(I+1)**2


cccc          WRITE(*,*)'OutIt=',iter,' InIt=', I,' INN. RESI:',W(I+1)


*     Check the accuracy and finish the loop.
         IF(DABS(W(I+1)) .LT. err) GO TO 123
107      CONTINUE

123      S(KK)       = W(KK) / C(KK,KK)
         DO 124  I   = 1,KK-1
            SUM      = 0.0D0
            DO 125 J = KK-I+1,KK
               SUM   = SUM + C(KK-I,J) * S(J)
125         CONTINUE
            S(KK-I)  = (W(KK-I) - SUM) / C(KK-I,KK-I)
124      CONTINUE

         DO 126 I    = 1,KK
            DO 127 J = 1,LEN
               V(J)  = 0.d0
127         CONTINUE
            V(I)     = 1.d0
            DO 128 J = I,1,-1
               SUM   = 0.0D0
               DO 129 M = 1,LEN
                  SUM   = SUM + Z(M,J) * V(M)
129            CONTINUE
               DO 130 M = 1,LEN
                  V(M)  = V(M) - 2.0D0 * SUM* Z(M,J)
130            CONTINUE
128         CONTINUE
            DO 131 M = 1,LEN
               X(M)  = X(M) + S(I) * V(M)
131         CONTINUE
126      CONTINUE


100   CONTINUE

c     Here is the end of outer loop iteration after accur. criterion.
200   CONTINUE


cc      WRITE(*,*)'Sea Level Residual:',DSQRT(resi)
      RETURN
      END

      SUBROUTINE PROD(diag,km2,V,dz,adz,il,jl,il1,jl1,ilp,jlp,len)
      dimension diag(19,il,jl), km2(0:il1,0:jl1),
     *          v(len), dz(-1:ilp,-1:jlp)
      dimension adz(-1:ilp,-1:jlp)
      real*8 diag, v, dz, adz

c     Разворачиваем одномерный массив в двумерный v -> dz.
      m= 0
      
      do j=1,jl
      do i=1,il
      if( km2(i,j) .GT. 0) then
      if (i .eq. 11 .and. j .eq. 1) then
      end if
      m= m+1
      dz(i,j)= v(m)
      end if
      end do
      end do
      
      
c     Умножаем на матрицу. Результат в ADZ.
      do j=1,jl
      do i=1,il
      if( km2(i,j) .GT. 0) then
      adz(i,j)=
     *  diag(1 ,i,j)*dz(i  ,j-2) + diag(2 ,i,j)*dz(i+1,j-2)+
     *  diag(3 ,i,j)*dz(i+2,j-2) + diag(4 ,i,j)*dz(i-1,j-1)+
     *  diag(5 ,i,j)*dz(i  ,j-1) + diag(6 ,i,j)*dz(i+1,j-1)+
     *  diag(7 ,i,j)*dz(i+2,j-1) + diag(8 ,i,j)*dz(i-2,j  )+
     *  diag(9 ,i,j)*dz(i-1,j  ) + diag(10,i,j)*dz(i  ,j  )+
     *  diag(11,i,j)*dz(i+1,j  ) + diag(12,i,j)*dz(i+2,j  )+
     *  diag(13,i,j)*dz(i-2,j+1) + diag(14,i,j)*dz(i-1,j+1)+
     *  diag(15,i,j)*dz(i  ,j+1) + diag(16,i,j)*dz(i+1,j+1)+
     *  diag(17,i,j)*dz(i-2,j+2) + diag(18,i,j)*dz(i-1,j+2)+
     *  diag(19,i,j)*dz(i  ,j+2)
      if (adz(i,j) .ne. adz(i,j)) then
      !  write(*,*) "ADZZZZZZZZZZZZ", i, j
        end if
      else
      adz(i,j)= 0.d0
      end if
      end do
      end do

c     Делаем предобуславливание.

      call prdiag(diag,km2,adz,il,jl,ilp,jlp,il1,jl1)

c     Переписываем в одномерный массив v.
      m= 0
      do j=1,jl
      do i=1,il
      if( km2(i,j) .GT. 0) then
      m= m+1
      v(m)= adz(i,j)
      end if
      end do
      end do

      return
      end

      subroutine prdiag(diag,km2,adz,il,jl,ilp,jlp,il1,jl1)

      dimension diag(19,il,jl), km2(0:il1,0:jl1), adz(-1:ilp,-1:jlp)
      real*8 diag,adz

      do j=1,jl
      do i=1,il
      if(km2(i,j) .GT. 0) then
      adz(i,j)= adz(i,j)/ diag(10,i,j)
      else
      adz(i,j)= 0.d0
      end if
      end do
      end do

      return
      end

