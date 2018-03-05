      SUBROUTINE Sadvect
******************************************************
*     Salinity 3D trasport and Horizontal diffusion
*     FEM Streamline Upwind Scheme Taylor-Galerkin
*     No Lumping.
*     Explicit time stepping.
*     Special treatment of open boundaries.
*     Version 10.05.2012
******************************************************

      parameter (itermax=10, omega=1.0)  ! iteration parameters

      INCLUDE 'slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      common /diffus/ al,alt
	common /consts/ c3,c6,c9,c12,c18,Rhx,R2hx2
cc	double precision dzmean, dzmmean
      INCLUDE 'tparm.fi'

      Tdamp_In = 24.*3600.      ! Time to restore at inflow
      Tdamp_Out=180.*24.*3600.  ! and outflow points.

	c6= 1./6.
	c3= 1./3.
	c9 =c3*c3
	c2=0.5
	c4=0.25
	c12=.5*c6
	c18=c3*c6
	c36=c6*c6
	c24=0.5*c12
	R2=R*R
      asr=hx/hy
	asr2=asr*asr
	Rhx=R*hx
	Rhy=R*hy
	R2hx2=c18*Rhx*Rhy 
	R2hx2=c18*Rhx*Rhy 

      COEFX=c18/Rhx
      COEFY=c18/Rhy


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
      A=.5/(S0*Rhx*Rhy)
      r2s12= c12*S0
                     
      DO I=1,IL

      KB = KM2(I,J)
      IF(nt3(i,j,1).GT.0) THEN

*-------------------------------------------------
*         Surface and Nonbottom points
*-------------------------------------------------

      DO K=1, KB-1


	kp= k+1
      n= nt3(i,j,k)
      np=nt3(i,j,k+1)
	hzk = hz(k)

      IF( k .EQ. 1) THEN
        hzk1= 0.
      ELSE
        hzk1= hz(k-1)
      END IF

      CA= dt

      call DIFFUW_TG(Sm2,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,il,jl,dt,KT,ALT,hx,hy,r)

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
	sumkp=
     *c4*(KT(1,np)*(2.*Sm2(i,j,k)+Sm2(i+1,j,k)+Sm2(i+1,j-1,k))
     +   +KT(2,np)*(2.*Sm2(i,j,k)+Sm2(i,j-1,k)+Sm2(i+1,j-1,k))
     +   +KT(3,np)*(2.*Sm2(i,j,k)+Sm2(i,j-1,k)+Sm2(i-1,j  ,k))
     +   +KT(4,np)*(2.*Sm2(i,j,k)+Sm2(i-1,j,k)+Sm2(i-1,j+1,k))
     +   +KT(5,np)*(2.*Sm2(i,j,k)+Sm2(i-1,j+1,k)+Sm2(i,j+1,k))
     +   +KT(6,np)*(2.*Sm2(i,j,k)+Sm2(i,j+1,k)+Sm2(i+1,j  ,k)))
	sumkpp=
     *c4*(KT(1,np)*(2.*Sm2(i,j,kp)+Sm2(i+1,j,kp)+Sm2(i+1,j-1,kp))
     +   +KT(2,np)*(2.*Sm2(i,j,kp)+Sm2(i,j-1,kp)+Sm2(i+1,j-1,kp))
     +   +KT(3,np)*(2.*Sm2(i,j,kp)+Sm2(i,j-1,kp)+Sm2(i-1,j  ,kp))
     +   +KT(4,np)*(2.*Sm2(i,j,kp)+Sm2(i-1,j,kp)+Sm2(i-1,j+1,kp))
     +   +KT(5,np)*(2.*Sm2(i,j,kp)+Sm2(i-1,j+1,kp)+Sm2(i,j+1,kp))
     +   +KT(6,np)*(2.*Sm2(i,j,kp)+Sm2(i,j+1,kp)+Sm2(i+1,j  ,kp)))


      FZT=-r2s12*((sumkp+sumkpp)*W(I,J,2)
     *                -2.*sumkp*(W(i,j,1)-PME(i,j)/dt) 
     *                                                     )
      ELSE
*     --------------  Deep water -------------
	sumkp=
     *c4*(KT(1,np)*(2.*Sm2(i,j,k)+Sm2(i+1,j,k)+Sm2(i+1,j-1,k))
     +   +KT(2,np)*(2.*Sm2(i,j,k)+Sm2(i,j-1,k)+Sm2(i+1,j-1,k))
     +   +KT(3,np)*(2.*Sm2(i,j,k)+Sm2(i,j-1,k)+Sm2(i-1,j  ,k))
     +   +KT(4,np)*(2.*Sm2(i,j,k)+Sm2(i-1,j,k)+Sm2(i-1,j+1,k))
     +   +KT(5,np)*(2.*Sm2(i,j,k)+Sm2(i-1,j+1,k)+Sm2(i,j+1,k))
     +   +KT(6,np)*(2.*Sm2(i,j,k)+Sm2(i,j+1,k)+Sm2(i+1,j  ,k)))
	sumkpp=
     *c4*(KT(1,np)*(2.*Sm2(i,j,kp)+Sm2(i+1,j,kp)+Sm2(i+1,j-1,kp))
     +   +KT(2,np)*(2.*Sm2(i,j,kp)+Sm2(i,j-1,kp)+Sm2(i+1,j-1,kp))
     +   +KT(3,np)*(2.*Sm2(i,j,kp)+Sm2(i,j-1,kp)+Sm2(i-1,j  ,kp))
     +   +KT(4,np)*(2.*Sm2(i,j,kp)+Sm2(i-1,j,kp)+Sm2(i-1,j+1,kp))
     +   +KT(5,np)*(2.*Sm2(i,j,kp)+Sm2(i-1,j+1,kp)+Sm2(i,j+1,kp))
     +   +KT(6,np)*(2.*Sm2(i,j,kp)+Sm2(i,j+1,kp)+Sm2(i+1,j  ,kp)))

	km=k-1

	sumk=
     *c4*(KT(1,n)*(2.*Sm2(i,j,k)+Sm2(i+1,j,k)+Sm2(i+1,j-1,k))
     +   +KT(2,n)*(2.*Sm2(i,j,k)+Sm2(i,j-1,k)+Sm2(i+1,j-1,k))
     +   +KT(3,n)*(2.*Sm2(i,j,k)+Sm2(i,j-1,k)+Sm2(i-1,j  ,k))
     +   +KT(4,n)*(2.*Sm2(i,j,k)+Sm2(i-1,j,k)+Sm2(i-1,j+1,k))
     +   +KT(5,n)*(2.*Sm2(i,j,k)+Sm2(i-1,j+1,k)+Sm2(i,j+1,k))
     +   +KT(6,n)*(2.*Sm2(i,j,k)+Sm2(i,j+1,k)+Sm2(i+1,j  ,k)))

	sumkm=
     *c4*(KT(1,n)*(2.*Sm2(i,j,km)+Sm2(i+1,j,km)+Sm2(i+1,j-1,km))
     +   +KT(2,n)*(2.*Sm2(i,j,km)+Sm2(i,j-1,km)+Sm2(i+1,j-1,km))
     +   +KT(3,n)*(2.*Sm2(i,j,km)+Sm2(i,j-1,km)+Sm2(i-1,j  ,km))
     +   +KT(4,n)*(2.*Sm2(i,j,km)+Sm2(i-1,j,km)+Sm2(i-1,j+1,km))
     +   +KT(5,n)*(2.*Sm2(i,j,km)+Sm2(i-1,j+1,km)+Sm2(i,j+1,km))
     +   +KT(6,n)*(2.*Sm2(i,j,km)+Sm2(i,j+1,km)+Sm2(i+1,j  ,km)))


      FZT=-r2s12*((sumkpp+sumkp)*W(I,J,K+1)-(sumk+sumkm)*W(I,J,K))      
      END IF

      FT=-FXT-FYT+FZT  +DIFT+DIFTP

      RsT(i,j,k)= CA*FT
c	Sm1(i,j,k)= omega*RsT(i,j,k)/sclamp(i,j,k)
	Sm1(i,j,k)= RsT(i,j,k)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

	end do   ! k

*     --------------- Bottom ------------------
      N = NP
      CA= dt
      call DIFFUW_TG(Sm2,u,v,w,hz,i,j,kb,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,il,jl,dt,KT,ALT,hx,hy,r)
      DIFT=A*DIFT
      call transp(FXTK,FYTK,i,j,kb,n,u,v,Sm2,KT,il1,jl1,kl)

      FXTK=COEFX*FXTK
      FYTK=COEFY*FYTK
      FXT=.5*FXTK*HZ(KB-1)
      FYT=.5*FYTK*HZ(KB-1)

	km=kb-1

	sumk=
     *c4*(KT(1,n)*(2.*Sm2(i,j,kb)+Sm2(i+1,j,kb)+Sm2(i+1,j-1,kb))
     +   +KT(2,n)*(2.*Sm2(i,j,kb)+Sm2(i,j-1,kb)+Sm2(i+1,j-1,kb))
     +   +KT(3,n)*(2.*Sm2(i,j,kb)+Sm2(i,j-1,kb)+Sm2(i-1,j  ,kb))
     +   +KT(4,n)*(2.*Sm2(i,j,kb)+Sm2(i-1,j,kb)+Sm2(i-1,j+1,kb))
     +   +KT(5,n)*(2.*Sm2(i,j,kb)+Sm2(i-1,j+1,kb)+Sm2(i,j+1,kb))
     +   +KT(6,n)*(2.*Sm2(i,j,kb)+Sm2(i,j+1,kb)+Sm2(i+1,j  ,kb)))

	sumkm=
     *c4*(KT(1,n)*(2.*Sm2(i,j,km)+Sm2(i+1,j,km)+Sm2(i+1,j-1,km))
     +   +KT(2,n)*(2.*Sm2(i,j,km)+Sm2(i,j-1,km)+Sm2(i+1,j-1,km))
     +   +KT(3,n)*(2.*Sm2(i,j,km)+Sm2(i,j-1,km)+Sm2(i-1,j  ,km))
     +   +KT(4,n)*(2.*Sm2(i,j,km)+Sm2(i-1,j,km)+Sm2(i-1,j+1,km))
     +   +KT(5,n)*(2.*Sm2(i,j,km)+Sm2(i-1,j+1,km)+Sm2(i,j+1,km))
     +   +KT(6,n)*(2.*Sm2(i,j,km)+Sm2(i,j+1,km)+Sm2(i+1,j  ,km)))


      FZT=r2s12*(sumk+sumkm)*W(I,J,KB)

      FT=-FXT-FYT+FZT  +DIFT

      RsT(i,j,kb)= CA*FT
c	Sm1(i,j,kb)= RsT(i,j,kb)/sclamp(i,j,kb)
	Sm1(i,j,kb)= RsT(i,j,kb)/( r2s12*CG(n)*HZ(KB-1))

      END IF   ! nt3 >0
	end do
	end do

cc      GOTO 1000

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
	n=nt3(i,j,k)

	if(k.LT.kb)then
	kp=k+1
	np=nt3(i,j,kp)
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


	else
	hzk= 0.
	np= n
	sumkp  =0.
	sumkpp =0.
	end if

	if(k.GT.1)then
	km=k-1
	hzk1=hz(km)

	sumk=
     *c24*(KT(1,n)*S1*(2.*Sm1(i,j,k)+Sm1(i+1,j,k)+Sm1(i+1,j-1,k))
     +    +KT(2,n)*S2*(2.*Sm1(i,j,k)+Sm1(i,j-1,k)+Sm1(i+1,j-1,k))
     +    +KT(3,n)*S3*(2.*Sm1(i,j,k)+Sm1(i,j-1,k)+Sm1(i-1,j  ,k))
     +    +KT(4,n)*S4*(2.*Sm1(i,j,k)+Sm1(i-1,j,k)+Sm1(i-1,j+1,k))
     +    +KT(5,n)*S5*(2.*Sm1(i,j,k)+Sm1(i-1,j+1,k)+Sm1(i,j+1,k))
     +    +KT(6,n)*S6*(2.*Sm1(i,j,k)+Sm1(i,j+1,k)+Sm1(i+1,j  ,k)))


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
	sumk =0.
	end if

	sum= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp))

      S(i,j,k)= Sm1(i,j,k) 
     *        + omega*(RsT(i,j,k)-sum)/sclamp(i,j,k)

c	if(i.eq.22.and.j.eq.44.and.k.eq.32) 
c     *write(*,*)iter,S(i,j,k),sum,RsT(i,j,k),Sm2(i,j,k)


	end do   ! k

      END IF   ! nt3 >0
	end do
	end do

	Sm1(:,:,:)= S(:,:,:)

	END DO   ! Mass matrix iteration


1000  continue


	do k=1,kl
	do j=1,jl
	do i=1,il
	if(nt3(i,j,k).GT.0) S(i,j,k)=Sm2(i,j,k)+Sm1(i,j,k)
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
