      subroutine dzvel
*     Ocean Velocity Caused by Sea Level.
*     Version 05.02.2010 for time-splitting scheme
      include 'Slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

      Ax= G/Hx
      Ay= G/Hy
      asr=hx/hy
	asr2=asr*asr

      DO 41 j=1,jL
      S0=Si(j)
      SM=Si(j-1)
      SP=Si(j+1)
      S1=(2.*S0+SM)/3.
      S2=(2.*SM+S0)/3.
      S3=S1
      S4=(2.*S0+SP)/3.
      S5=(2.*SP+S0)/3.
      S6=S4

      rs0= r*s0

      DO 41 i=1,iL

      cor= -2.*om*co(i)*S0
      CB= 0.5*cor*dt
cc      CB= cor*dt

      kk2= km2(i,j)
      IF( kk2 .GT. 0) THEN

      DO 42 k=1,kk2

      n= ABS(nt3(i,j,k))
      np=n
      if(k.lt.kk2) np=ABS(nt3(i,j,k+1))

      IF( k .EQ. 1) THEN
      hzk1= 0.
      ELSE
      hzk1= hz(k-1)
      END IF

      IF( k .EQ. kk2) THEN
      hzk= 0.
      ELSE
      hzk= hz(k)
      END IF

      CA= 12.*DT/(RS0*(CG(n)*HZK1+CG(np)*HZK))
      CD= CA/(1.+CB**2)

      if( k .eq. 1) then
      derx= 0.
      dery= 0.
      else
      call graddz(i,j,N ,derx ,dery ,dz,ilp,jlp,KT)
      end if

      if( k.lt.kk2) then
      call graddz(i,j,NP,derxp,deryp,dz,ilp,jlp,KT)
      else
      derxp=0.
      deryp=0.
      end if

      fu= 0.25*Ax*( hzk1*derx +hzk*derxp)
      fv= 0.25*Ay*( hzk1*dery +hzk*deryp)

c      fu= 0.5*Ax*( hzk1*derx +hzk*derxp)
c      fv= 0.5*Ay*( hzk1*dery +hzk*deryp)
      if (i .eq. 11 .and. j .eq. 1) then
      !  write(*,*) "U(i,j,K)", U(i,j,K), CD, FU, CB, FV
      !  write(*,*) "Ax", Ax, hzk1, derx, hzk, derxp
      !  write(*,*) "KT", KT
      !  write(*,*) "dz", dz
      end if
      U(i,j,K)=U(i,j,K) +CD*(FU+CB*FV)
      
      V(i,j,K)=V(i,j,K) +CD*(FV-CB*FU)

42    CONTINUE
      END IF
41    CONTINUE
      RETURN
      END

      subroutine rside(ncond,rs)
*     -----------------------------------
*     Sea Level Equation Right Hand Side
*     -----------------------------------

*     Version 03.07.2015 for free slip condition
*     and Flather b.c. at open boundary modified for ice,
*     snow and Pa.

      INCLUDE 'Slo2.fi'
      dimension rs(il,jl)
	double precision rs,asr,rhx,rrhx,Hmean,dzmean,PME_Mean,
     *                 s0,sm,sp,s1,s2,s3,s4,s5,s6,F,divk,divkp
	Include 'Tparm.fi'

      Hmean= 1.d5
     
      AsR = DBLE(Hx/Hy)
      RHx = DBLE(R*Hx)
	RRHX= 1.d0/RHx

      DO 1 j=1,jL
      S0=DBLE(Si(j))
      SM=DBLE(Si(j-1))
      SP=DBLE(Si(J+1))
      S1=(2.d0*S0+SM)/3.d0
      S2=(2.d0*SM+S0)/3.d0
      S3=S1
      S4=(2.d0*S0+SP)/3.d0
      S5=(2.d0*SP+S0)/3.d0
      S6=S4

      DO 1 i=1,iL
      kk2= km2(i,j)
      n=abs(nt3(i,j,1))
      IF( kk2 .GT. 0) THEN

*     -----------------------------------
      IF( Ncond .EQ. 1) THEN

	DzMean=( DBLE(KT(1,n))*S1*(2.d0*dz(i,j)+dz(i+1,j)+dz(i+1,j-1))
     +        +DBLE(KT(2,n))*S2*(2.d0*dz(i,j)+dz(i+1,j-1)+dz(i,j-1))
     +        +DBLE(KT(3,n))*S3*(2.d0*dz(i,j)+dz(i,j-1)+dz(i-1,j))
     +        +DBLE(KT(4,n))*S4*(2.d0*dz(i,j)+dz(i-1,j)+dz(i-1,j+1))
     +        +DBLE(KT(5,n))*S5*(2.d0*dz(i,j)+dz(i,j+1)+dz(i-1,j+1))
     +        +DBLE(KT(6,n))*S6*(2.d0*dz(i,j)+dz(i+1,j)+dz(i,j+1))
     *                                                    )/4.d0
     
      PME_mean= DBLE(PME(i,j))
     *         *(DBLE(KT(1,n))*S1+DBLE(KT(2,n))*S2+DBLE(KT(3,n))*S3+
     *           DBLE(KT(4,n))*S4+DBLE(KT(5,n))*S5+DBLE(KT(6,n))*S6) 


      F=(DzMean+PME_mean)/(6.d0*DBLE(DT))
      !write(*,*) "DzMean", DzMean, "PME_mean", PME_mean, i, j
      !write(*,*) "DBLE(PME(i,j))", DBLE(PME(i,j)) !! 16, 41
      else
      F= 0.d0
      end if      

*     Pressure by snow and ice   
	hi= 0.
      do mg=1,mgrad
	hi=hi+rosdry*hsnow(mg,i,j)/row ! Normalized snow mass
	hi=hi+roi   *hice (mg,i,j)/row ! Normalized ice mass
	end do
	
*     Atmospheric pressure	
*	hi= hi +10.*PA(i,j)/(g*row)
	
      nb= 1
      N= NT3(i,j,1)
      IF( n .LT. 0) nb=-1

      DO 2 K=1,KK2-1

      N= ABS(NT3(i,j,k+1))
      call divergD(n ,nb,i,j,k  ,divk ,u,v,ub,dzext,hi,km2,z,
     *             il1,jl1,il,jl,kl,klp,KT,g,s0,s1,s2,s3,s4,s5,s6,asr)
      call divergD(n ,nb,i,j,k+1,divkp,u,v,ub,dzext,hi,km2,z,
     *             il1,jl1,il,jl,kl,klp,KT,g,s0,s1,s2,s3,s4,s5,s6,asr)

2     F= F+ RRHx*DBLE(hz(k))*(divkp+divk)/12.d0
      rs(i,j)= F/Hmean

      END IF
1     CONTINUE

c     Check up the Fredholm alternative for the case of
c     only inflow waters forcing: INT(rs)dxdy = 0.

c      do i=1,il
c      do j=1,jl
c      q= q+rs(i,j)
c      end do
c      end do
c      write(*,*)' Mean RS: ', q

      return
      end

      subroutine gradDz(i,j, N,derx,dery,dz,ilp,jlp,KT)
      dimension dz(-1:ilp,-1:jlp), KT(6,13)
      real KT
      double precision dz
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

      derx=-(   KT(1,n)*REAL(dz(i+1,j  )-dz(i  ,j  ))+
     +          KT(2,n)*REAL(dz(i+1,j-1)-dz(i  ,j-1))+
     +          KT(3,n)*REAL(dz(i  ,j  )-dz(i-1,j  ))+
     +          KT(4,n)*REAL(dz(i  ,j  )-dz(i-1,j  ))+
     +          KT(6,n)*REAL(dz(i+1,j  )-dz(i  ,j  ))+
     +          KT(5,n)*REAL(dz(i  ,j+1)-dz(i-1,j+1)))/6.

      dery=-(s1*KT(1,n)*REAL(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*KT(2,n)*REAL(dz(i  ,j  )-dz(i  ,j-1))+
     +       s3*KT(3,n)*REAL(dz(i  ,j  )-dz(i  ,j-1))+
     +       s4*KT(4,n)*REAL(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*KT(5,n)*REAL(dz(i  ,j+1)-dz(i  ,j  ))+
     +       s6*KT(6,n)*REAL(dz(i  ,j+1)-dz(i  ,j  )))/6.
      return
      end

      subroutine DivergD 
     *(m,nb,i,j,k,divk,u,v,ub,dzext,hi,km2,z,il1,jl1,il,jl,kl,klp,KT,g,
     *s0,s1,s2,s3,s4,s5,s6,asr)

*     Version 27/05/15 for Flather b.c. at open boundary modified for ice.

      dimension u(0:il1,0:jl1,kl), v(0:il1,0:jl1,kl)
      dimension KT(6,13), ub(0:il1,0:jl1,kl)
	dimension dzext(il,jl), km2(0:il1,0:jl1), z(klp)
      real KT
      Double precision divk, s0,s1,s2,s3,s4,s5,s6,asr,add,coef

	coef=1.0d0

ccc	DZ_Trelax= 6.4e-8
***	DZ_Trelax= SQRT(g*z(km2(i,j))) * 1.e-7 !! Horizontal resolution

      DivK=    -DBLE(KT(1,m)*(U(I,J,K)+U(I+1,J,K)+U(I+1,J-1,K)))
     +     + S2*DBLE(KT(2,m)*(V(I,J,K)+V(I,J-1,K)+V(I+1,J-1,K)))*AsR
     +         +DBLE(KT(3,m)*(U(I,J,K)+U(I-1,J,K)+U(I  ,J-1,K)))
     +     + S3*DBLE(KT(3,m)*(V(I,J,K)+V(I-1,J,K)+V(I  ,J-1,K)))*AsR
     +         +DBLE(KT(4,m)*(U(I,J,K)+U(I-1,J,K)+U(I-1,J+1,K)))
     -     - S5*DBLE(KT(5,m)*(V(I,J,K)+V(I-1,J+1,K)+V(I,J+1,K)))*AsR
     -         -DBLE(KT(6,m)*(U(I,J,K)+U(I+1,J,K)+U(I  ,J+1,K)))
     -     - S6*DBLE(KT(6,m)*(V(I,J,K)+V(I+1,J,K)+V(I  ,J+1,K)))*AsR

           IF( nb .EQ. -1) THEN

***      add= - DBLE( (dzext(i,j)-hi) * DZ_Trelax)
      add= - coef*DBLE( (dzext(i,j)-hi) * SQRT(g/z(km2(i,j))))
ccc      add= - coef*DBLE( -hi * SQRT(g/z(km2(i,j))))
      IF(m.GE.10)
     *DivK=DivK-3.d0*(1.d0+S0*AsR)*(DBLE(Ub(i,j,k))+add)
      IF(m.EQ.6) DivK=DivK-6.d0*(DBLE(Ub(i,j,K)) +add)
      IF(m.EQ.7) DivK=DivK-6.d0*S0*AsR*(DBLE(Ub(i,j,k))+add)
      IF(m.EQ.8) DivK=DivK-6.d0*(DBLE(Ub(i,j,k))+add)
      IF(m.EQ.9) DivK=DivK-6.d0*S0*AsR*(DBLE(Ub(i,j,k))+add)

           END IF

      return
      end

