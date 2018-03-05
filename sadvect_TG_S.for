      SUBROUTINE Sadvect
******************************************************
*     Salinity 3D trasport and Horizontal diffusion
*     FEM Streamline Upwind Scheme Taylor-Galerkin
*     No Lumping.
*     Explicit time stepping.
*     Special treatment of open boundaries.
*     Version 15.08.2012
******************************************************

      parameter (omega=1.0, itermax=20)  ! iteration parameters

      INCLUDE 'slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      common /diffus/ al,alt
	common /consts/ c3,c6,c9,c12,c18,Rhx,R2hx2
cc	double precision dzmean, dzmmean
      INCLUDE 'tparm.fi'

      Tdamp_In =       2.* dt   ! Time to restore at inflow
      Tdamp_Out=180.*24.*3600.  ! and outflow points.

	c6= 1./6.
	c3= 1./3.
	c9 =c3*c3
	c2=0.5
	c4=0.25
	c12=.5*c6
	c18=c3*c6
	c24=0.5*c12
	R2=R*R
      asr=hx/hy
	asr2=asr*asr
	Rhx=R*hx
	Rhy=R*hy
	R2hx2=c18*Rhx*Rhy 
      HXX= 1./HX
      COEFX=c18*HXX/R
      HYY= 1./HY
      COEFY=c18*HYY/R


      RsT(:,:,:)=0.
      S  (:,:,:)=0.
      Sm1(:,:,:)=0.


      DO J=1,JL
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
      A=.5*HXX*HXX/S0/R2
      r2s12= c12*S0
                     
      DO I=1,IL

      KB = KM2(I,J)
      IF(nt3(i,j,1).GT.0) THEN

*-------------------------------------------------
*         Surface and Nonbottom points
*-------------------------------------------------
      n= ABS(nt3(i,j,1))

      DO K=1, KB -1
      n= abs(nt3(i,j,k))
      np=abs(nt3(i,j,k+1))
	hzk = hz(k)

      IF( k .EQ. 1) THEN
        hzk1= 0.
      ELSE
        hzk1= hz(k-1)
      END IF

      CA= dt


      call DIFFUW_TG(Sm2,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,dt,KT,ALT)

      if( k .GT. 1) then 
      call transp(FXTK,FYTK,i,j,k,n,u,v,Sm2,KT,il1,jl1,kl)
      else
      FXTK =0.
      FYTK =0.
      end if

      call transp(FXTK1,FYTK1,i,j,k,np,u,v,Sm2,KT,il1,jl1,kl)

      FXTK =COEFX*FXTK
      FYTK =COEFY*FYTK
      FXTK1=COEFX*FXTK1
      FYTK1=COEFY*FYTK1
      DIFT =A*DIFT
      DIFTP=A*DIFTP
      FXT=0.5*(FXTK*HZK1 +FXTK1*HZK)
      FYT=0.5*(FYTK*HZK1 +FYTK1*HZK)

      IF( k .EQ. 1) THEN
*     -------------- Surface -----------------
      FZT=-r2s12*cg(n)*((SM2(I,J,2)+SM2(I,J,1))*W(I,J,2)
     *    -2.*Sm2(I,J,1)*(W(i,j,1)-PME(i,j)/dt)
     *  )
      ELSE
*     --------------  Deep water -------------
      FZT=-r2s12*(cg(np)*(SM2(I,J,K+1)+SM2(I,J,K)) *W(I,J,K+1)
     -           -cg( n)*(SM2(I,J,K  )+SM2(I,J,K-1))*W(I,J,K)
     * )
      END IF

      FT=-FXT-FYT+FZT  +DIFT+DIFTP

      RsT(i,j,k)= CA*FT
	Sm1(i,j,k)= RsT(i,j,k)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

	end do   ! k

*     --------------- Bottom ------------------
      N = NP
      CA= dt
      call DIFFUW_TG(Sm2,u,v,w,hz,i,j,kb,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,dt,KT,ALT)
      DIFT=A*DIFT
      call transp(FXTK,FYTK,i,j,kb,n,u,v,Sm2,KT,il1,jl1,kl)

      FXTK=COEFX*FXTK
      FYTK=COEFY*FYTK
      FXT=.5*FXTK*HZ(KB-1)
      FYT=.5*FYTK*HZ(KB-1)

      FZT=r2s12*cg(n)*(SM2(I,J,KB-1)+SM2(I,J,KB))*W(I,J,KB)

      FT=-FXT-FYT+FZT  +DIFT

      RsT(i,j,kb)= CA*FT
	Sm1(i,j,kb)= RsT(i,j,kb)/( r2s12*CG(n)*HZ(KB-1))

c	if(i.eq.33.and.j.eq.30.) write(*,*)'init ',Sm1(i,j,kb)

      END IF   ! nt3 >0
	end do
	end do


*     Mass matrix invertion

      DO iter=1,itermax

	do j=1,jl
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
      r2s12= c12*S0
	
	do i=1,il
	IF( nt3(i,j,1) .GT. 0) then

	kb=km2(i,j)
      do k=1,kb
	n=abs(nt3(i,j,k))

	sumk=c24*(KT(1,n)*S1*(2.*Sm1(i,j,k)+Sm1(i+1,j,k)+Sm1(i+1,j-1,k))
     +         +KT(2,n)*S2*(2.*Sm1(i,j,k)+Sm1(i,j-1,k)+Sm1(i+1,j-1,k))
     +         +KT(3,n)*S3*(2.*Sm1(i,j,k)+Sm1(i,j-1,k)+Sm1(i-1,j  ,k))
     +         +KT(4,n)*S4*(2.*Sm1(i,j,k)+Sm1(i-1,j,k)+Sm1(i-1,j+1,k))
     +         +KT(5,n)*S5*(2.*Sm1(i,j,k)+Sm1(i-1,j+1,k)+Sm1(i,j+1,k))
     +         +KT(6,n)*S6*(2.*Sm1(i,j,k)+Sm1(i,j+1,k)+Sm1(i+1,j  ,k)))
c	aiik=c12*( KT(1,n)*S1+KT(2,n)*S2+KT(3,n)*S3
c     +          +KT(4,n)*S4+KT(5,n)*S5+KT(6,n)*S6)

	if(k.LT.kb)then
	kp=k+1
	np=abs(nt3(i,j,kp))
	hzk=hz(k)
	sumkp=
     *c24*(KT(1,np)*S1*(2.*Sm1(i,j,k)+Sm1(i+1,j,k)+Sm1(i+1,j-1,k))
     +    +KT(2,np)*S2*(2.*Sm1(i,j,k)+Sm1(i,j-1,k)+Sm1(i+1,j-1,k))
     +    +KT(3,np)*S3*(2.*Sm1(i,j,k)+Sm1(i,j-1,k)+Sm1(i-1,j  ,k))
     +    +KT(4,np)*S4*(2.*Sm1(i,j,k)+Sm1(i-1,j,k)+Sm1(i-1,j+1,k))
     +    +KT(5,np)*S5*(2.*Sm1(i,j,k)+Sm1(i-1,j+1,k)+Sm1(i,j+1,k))
     +    +KT(6,np)*S6*(2.*Sm1(i,j,k)+Sm1(i,j+1,k)+Sm1(i+1,j  ,k)))
	sumkpp=
     *c24*(KT(1,np)*S1*(2.*Sm1(i,j,kp)+Sm1(i+1,j,kp)+Sm1(i+1,j-1,kp))
     +    +KT(2,np)*S2*(2.*Sm1(i,j,kp)+Sm1(i,j-1,kp)+Sm1(i+1,j-1,kp))
     +    +KT(3,np)*S3*(2.*Sm1(i,j,kp)+Sm1(i,j-1,kp)+Sm1(i-1,j  ,kp))
     +    +KT(4,np)*S4*(2.*Sm1(i,j,kp)+Sm1(i-1,j,kp)+Sm1(i-1,j+1,kp))
     +    +KT(5,np)*S5*(2.*Sm1(i,j,kp)+Sm1(i-1,j+1,kp)+Sm1(i,j+1,kp))
     +    +KT(6,np)*S6*(2.*Sm1(i,j,kp)+Sm1(i,j+1,kp)+Sm1(i+1,j  ,kp)))
c	aiikp=c12*( KT(1,np)*S1+KT(2,np)*S2+KT(3,np)*S3
c     +           +KT(4,np)*S4+KT(5,np)*S5+KT(6,np)*S6)
	else
	hzk= 0.
	np=n
	sumkp=0.
	sumkpp=0.
c	aiikp=0.
	end if

	if(k.GT.1)then
	km=k-1
	hzk1=hz(km)
	sumkm=
     *c24*(KT(1,n)*S1*(2.*Sm1(i,j,km)+Sm1(i+1,j,km)+Sm1(i+1,j-1,km))
     +    +KT(2,n)*S2*(2.*Sm1(i,j,km)+Sm1(i,j-1,km)+Sm1(i+1,j-1,km))
     +    +KT(3,n)*S3*(2.*Sm1(i,j,km)+Sm1(i,j-1,km)+Sm1(i-1,j  ,km))
     +    +KT(4,n)*S4*(2.*Sm1(i,j,km)+Sm1(i-1,j,km)+Sm1(i-1,j+1,km))
     +    +KT(5,n)*S5*(2.*Sm1(i,j,km)+Sm1(i-1,j+1,km)+Sm1(i,j+1,km))
     +    +KT(6,n)*S6*(2.*Sm1(i,j,km)+Sm1(i,j+1,km)+Sm1(i+1,j  ,km)))
      
	else
	hzk1=0.
	sumkm=0.
	end if

	sum= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp))!/dt
c	aii= (HZK1*aiik +HZK*aiikp)!/dt

cc      S(i,j,k)= Sm1(i,j,k) +omega*(RsT(i,j,k)-sum)/aii
      S(i,j,k)= Sm1(i,j,k) 
     +        + omega*(RsT(i,j,k)-sum)/(r2s12*(CG(n)*HZK1+CG(np)*HZK))


c	if(i.eq.33.and.j.eq.30.and.k.eq.5) 
c     * write(*,*) iter, S(i,j,k), RsT(i,j,k)-sum


	end do   ! k

      END IF   ! nt3 >0
	end do
	end do

	Sm1(:,:,:)= S(:,:,:)

	END DO   ! iter

	do k=1,kl
	do j=1,jl
	do i=1,il
	if(nt3(i,j,k).GT.0) S(i,j,k)=Sm2(i,j,k)+S(i,j,k)
	end do
	end do
	end do

*     Liquid points
      call TSbc(Tdamp_In,Tdamp_Out,
     *          dt,S,Sm2,Sobs,u,v,w,nt3,Si,il1,jl1,kl,hx,hy,hz,R)

      do k=1,kl
	do j=1,jl
	do i=1,il
	s(i,j,k)=max(0.,s(i,j,k))
	end do
	end do
	end do
	 

      RETURN
      END
