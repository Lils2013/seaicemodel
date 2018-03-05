      SUBROUTINE Sadvect
******************************************************
*     Salinity 3D trasport and Horizontal diffusion
*     FEM Streamline Upwind Scheme Taylor-Galerkin
*     Lumping.
*     Explicit time stepping.
*     Special treatment of open boundaries.
*     Version 19.06.2012
******************************************************

      INCLUDE 'slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      common /diffus/ al,alt
	common /consts/ c3,c6,c9,c12,c18,Rhx,R2hx2
cc	double precision dzmean, dzmmean
      INCLUDE 'tparm.fi'

      Tdamp_In = dt             ! Time to restore at inflow
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
      FZT=-r2s12*cg(n)*((SM2(I,J,2)+SM2(I,J,1))*W(I,J,2)
     *    -2.*Sm2(I,J,1)*(W(i,j,1)-PME(i,j)/dt)
     *  )
      ELSE
*     --------------  Deep water -------------
      FZT=-r2s12*(cg(np)*(SM2(I,J,K+1)+SM2(I,J,K))  *W(I,J,K+1)
     -           -cg( n)*(SM2(I,J,K  )+SM2(I,J,K-1))*W(I,J,K)
     * )
      END IF

      FT=-FXT-FYT+FZT  +DIFT+DIFTP

	S(i,j,k)= Sm2(i,j,k)+CA*FT/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

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

      FZT= r2s12*cg(n)*(SM2(I,J,KB-1)+SM2(I,J,KB))*W(I,J,KB)

      FT=-FXT-FYT+FZT  +DIFT

	S(i,j,kb)= Sm2(i,j,kb)+CA*FT/( r2s12*CG(n)*HZ(KB-1))

      END IF   ! nt3 >0
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
