      SUBROUTINE Tadvect
******************************************************
*     Temperature 3D trasport and Horizontal diffusion
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
      INCLUDE 'tparm.fi'

      Tdamp_In = dt             ! Time to restore at inflow
      Tdamp_Out=180.*24.*3600.  ! and outflow points.

	c6= 1./6.
	c3= 1./3.
	c9 =c3*c3
	c2=0.5
	c4=0.25
	c12=0.5*c6
	c18=c3*c6
	c24=0.5*c12
	c36=c6*c6
      asr=hx/hy
	asr2=asr*asr
	R2=R*R
	Rhx=R*hx
	Rhy=R*hy
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

      DO K=1, KB -1
      n= nt3(i,j,k)
      np=nt3(i,j,k+1)
	hzk = hz(k)

      IF( k .EQ. 1) THEN
        hzk1= 0.
      ELSE
        hzk1= hz(k-1)
      END IF
	 
      CA= dt

      call DIFFUW_TG(Tm2,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,il,jl,dt,KT,ALT,hx,hy,r)

      if( k .GT. 1) then 
      call transp(FXTK,FYTK,i,j,k,n,u,v,tm2,KT,il1,jl1,kl)
      else
      FXTK =0.
      FYTK =0.
      end if
      call transp(FXTK1,FYTK1,i,j,k,np,u,v,tm2,KT,il1,jl1,kl)

      FXTK =COEFX*FXTK
      FYTK =COEFY*FYTK
      FXTK1=COEFX*FXTK1
      FYTK1=COEFY*FYTK1
      DIFT =A*DIFT
      DIFTP=A*DIFTP
      FXT=.5*(FXTK*HZK1 +FXTK1*HZK)
      FYT=.5*(FYTK*HZK1 +FYTK1*HZK)

      IF( k .EQ. 1) THEN
*     -------------- Surface -----------------
      FZT=-r2s12*cg(n)*((TM2(I,J,2)+TM2(I,J,1))*W(I,J,2)
     *                 -2.*Tm2(I,J,1)*(W(i,j,1)-PME(i,j)/dt) 
     *                                         )
      ELSE
*     --------------  Deep water -------------
      FZT=-r2s12*(cg(np)*(TM2(I,J,K+1)+TM2(I,J,K))*W(I,J,K+1)
     -               -cg( n)*(TM2(I,J,K  )+TM2(I,J,K-1))*W(I,J,K) 
     * )
      END IF

      FT=-FXT-FYT+FZT  +DIFT+DIFTP

	T(i,j,k)=Tm2(i,j,k) +CA*FT/( r2s12*(CG(n)*HZK1+CG(np)*HZK))
      end do

*     --------------- Bottom ------------------
      N = NP
      CA= dt   
      call DIFFUW_TG(Tm2,u,v,w,hz,i,j,kb,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,il,jl,dt,KT,ALT,hx,hy,r)
      DIFT=A*DIFT
      call transp(FXTK ,FYTK ,i,j,kb  ,n,u,v,tm2,KT,il1,jl1,kl)

      FXTK=COEFX*FXTK
      FYTK=COEFY*FYTK
      FXT=.5*FXTK*HZ(KB-1)
      FYT=.5*FYTK*HZ(KB-1)

      FZT=r2s12*cg(n)*(TM2(I,J,KB-1)+TM2(I,J,KB))*W(I,J,KB)

      FT=-FXT-FYT+FZT  +DIFT

	T(i,j,kb)= Tm2(i,j,kb)+CA*FT/( r2s12*CG(n)*HZ(KB-1))

      END IF
     
      end do
	end do

*     Liquid points
      call TSbc(Tdamp_In,Tdamp_Out,
     *          dt,T,Tm2,Tobs,u,v,w,nt3,Si,il1,jl1,kl,hx,hy,hz,R)

	do j=1,jl
	do i=1,il
	do k=1,km2(i,j)
	T(i,j,k)= max( T(i,j,k), Tfr(S(i,j,k),0.))
	end do
	end do
	end do


      RETURN
      END

      subroutine transp(FXTK,FYTK,i,j,k,n,u,v,t,KT,il,jl,kl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl),
     *          t(0:il,0:jl,kl), KT(6,13)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

*     Version 16.10.2003.

      FXTK= 0.
      FYTK= 0.
*---------------     Triangle # 1.  ---------------------------------
      FXTK=FXTK+(U(I,J,K)+U(I+1,J,K)+U(I+1,J-1,K))*
     *          (T(I  ,J  ,K)+T(I+1,J  ,K)+T(I+1,J-1,K))
     * *KT(1,n)
*---------------     Triangle # 2.  -----------------------------------
      FYTK=FYTK-S2*(V(I,J,K)+V(I,J-1,K)+V(I+1,J-1,K))*
     *             (T(I  ,J  ,K)+T(I  ,J-1,K)+T(I+1,J-1,K))
     *  *KT(2,n)
*---------------     Triangle # 3.  -------------------------------
      FXTK=FXTK-(U(I,J,K)+U(I-1,J,K)+U(I,J-1,K))*
     *          (T(I  ,J  ,K)+T(I-1,J  ,K)+T(I  ,J-1,K))
     *  *KT(3,n)
      FYTK=FYTK-S3*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K))*
     *             (T(I  ,J  ,K)+T(I  ,J-1,K)+T(I-1,J  ,K))
     *  *KT(3,n)
*---------------     Triangle # 4.  ---------------------------------
      FXTK=FXTK-(U(I,J,K)+U(I-1,J,K)+U(I-1,J+1,K))*
     *          (T(I  ,J  ,K)+T(I-1,J  ,K)+T(I-1,J+1,K))
     *  *KT(4,n)
*---------------     Triangle # 5.  ----------------------------------
      FYTK=FYTK+S5*(V(I,J,K)+V(I,J+1,K)+V(I-1,J+1,K))*
     *             (T(I  ,J  ,K)+T(I  ,J+1,K)+T(I-1,J+1,K))
     *  *KT(5,n)
*---------------     Triangle # 6.  -----------------------------------
      FXTK=FXTK+(U(I,J,K)+U(I+1,J,K)+U(I,J+1,K))*
     *          (T(I  ,J  ,K)+T(I+1,J  ,K)+T(I  ,J+1,K))
     *  *KT(6,n)
      FYTK=FYTK+S6*(V(I,J,K)+V(I,J+1,K)+V(I+1,J,K))*
     *             (T(I  ,J  ,K)+T(I  ,J+1,K)+T(I+1,J  ,K))
     *  *KT(6,n)



      return
      end

	subroutine TSbc(Tdamp_In,Tdamp_Out,
     *           dt,T,Tm2,Tobs,u,v,w, nt3,Si,il,jl,kl,hx,hy,hz,R)

*     Version 14/11/2006 with vertical advection

      dimension T(0:il,0:jl,kl), Tm2(0:il,0:jl,kl), hz(kl),
     &          Tobs(0:il,0:jl,kl), nt3(0:il,0:jl,kl),Si(0:jl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)

	do j=1,jl-1
	S0=Si(j)
      do i=1,il-1
	do k=1,kl

	np=nt3(i,j,k)

*     ------------ Vertical advection - explicit upwinding ----------

      if(np .NE. 0) then

      FZ = 0.
	if( k.LT. kl) then
      FZ= FZ + min(0., w(i,j,k+1))*(Tm2(i,j,k+1)-Tm2(i,j,k))/hz(k)
	end if
	if (k.GT.1) then
      FZ= FZ + max(0., w(i,j,k  ))*(Tm2(i,j,k)-Tm2(i,j,k-1))/hz(k-1)
	end if

	end if  ! np .ne. 0

*     ------------------- Vertical advection ------------------------

	if(np.EQ.-6) then

      delta_t=(T(i+1,j,k)-Tm2(i+1,j,k))
	delta_x=(T(i+1,j,k)-T(i+2,j,k))/(R*hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))/(R*hy)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))/(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/(R*hy)

c	umean= (2.*u(i,j,k)+2.*u(i+1,j,k)+u(i+1,j-1,k)+u(i,j+1,k))/6.
c	vmean= (2.*v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)
c     &+v(i+1,j-1,k)+v(i,j+1,k))/6.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

	Rx= Rx -dt*u(i+1,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i+1,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-6

	if(np.EQ.-9) then

      delta_t=(T(i,j-1,k)-Tm2(i,j-1,k))
	delta_x=(T(i,j-1,k)-T(i,j-2,k))/(R*hy)
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))/(R*hx*S0)
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

	Rx= Rx +dt*v(i,j-1,k)/(R*hy)
	Ry= Ry +dt*u(i,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-9

	if(np.EQ.-8) then

      delta_t=(T(i-1,j,k)-Tm2(i-1,j,k))
	delta_x=(T(i-1,j,k)-T(i-2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k)) /(R*hy)
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/ (R*hy) 

	Rx= Rx +dt*u(i-1,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i-1,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-8

	if(np.EQ.-7) then

      delta_t=(T(i,j+1,k)-Tm2(i,j+1,k))
	delta_x=(T(i,j+1,k)-T(i,j+2,k)) /(R*hy)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(R*hx*S0)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

	Rx= Rx -dt*v(i,j+1,k)/(R*hy)
	Ry= Ry -dt*u(i,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-7

	if(np.EQ.-10) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to W

      delta_t=(T(i+1,j,k)-Tm2(i+1,j,k))
	delta_x=(T(i+1,j,k)-T(i+2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k)) /(R*hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hx*S0)
	Ry=-delta_t*delta_y*D /(R*hy)      

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/3.
c	vmean= (v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/3.

c      Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Rx +dt*vmean/(R*hy)

	Rx= Rx -dt*u(i+1,j-1,k)/(R*hx*S0)
	Ry= Ry +dt*v(i+1,j-1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=(T(i,j-1,k)-Tm2(i,j-1,k))
	delta_x=(T(i,j-1,k)-T(i,j-2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

	Rx= Rx +dt*v(i+1,j-1,k)/(R*hy)
	Ry= Ry +dt*u(i+1,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To S

	end if  ! n=-10


	if(np.EQ.-11) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to E
      delta_t=(T(i-1,j,k)-Tm2(i-1,j,k))
	delta_x=(T(i-1,j,k)-T(i-2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))/(R*hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/(R*hy)

	Rx= Rx +dt*u(i-1,j-1,k)/(R*hx*S0)
	Ry= Ry -dt*v(i-1,j-1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=(T(i,j-1,k)-Tm2(i,j-1,k))
	delta_x=(T(i,j-1,k)-T(i,j-2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k)) /(R*hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

	Rx= Rx +dt*v(i-1,j-1,k)/(R*hy)
	Ry= Ry +dt*u(i-1,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
    
	end if  ! To S

	end if  ! n=-11

	if(np.EQ.-12) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to E
      delta_t=(T(i-1,j,k)-Tm2(i-1,j,k))
	delta_x=(T(i-1,j,k)-T(i-2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/(R*hy)

	Rx= Rx +dt*u(i-1,j+1,k)/(R*hx*S0)
	Ry= Ry -dt*v(i-1,j+1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=(T(i,j+1,k)-Tm2(i,j+1,k))
	delta_x=(T(i,j+1,k)-T(i,j+2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(R*hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hy)
	Ry=-delta_t*delta_y*D /(R*hx*S0)

	Rx= Rx -dt*v(i-1,j+1,k)/(R*hy)
	Ry= Ry -dt*u(i-1,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To N

	end if  ! n=-12

	if(np.EQ.-13) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to W

      delta_t=(T(i+1,j,k)-Tm2(i+1,j,k))
	delta_x=(T(i+1,j,k)-T(i+2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k)) /(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hx*S0)
	Ry=-delta_t*delta_y*D /(R*hy) 

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))/3.
c	vmean= (v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))/3.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

	Rx= Rx -dt*u(i+1,j+1,k)/(R*hx*S0)
	Ry= Ry +dt*v(i+1,j+1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=(T(i,j+1,k)-Tm2(i,j+1,k))
	delta_x=(T(i,j+1,k)-T(i,j+2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i,j+1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hy)
	Ry=-delta_t*delta_y*D /(R*hx*S0)

	Rx= Rx -dt*v(i+1,j+1,k)/(R*hy)
	Ry= Ry -dt*u(i+1,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if

      end if  ! To N
	end if  ! n=-13

	end do
	end do
	end do

	Return
	End

      SUBROUTINE DIFFUW_TG(T,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &                     il,jl,kl,ilm,jlm,dt,KT,A,hx,hy,R)

*     Version 01.06.2012.

*     Galerkin artificial diffusion approximation by Taylor-Galerkin
*     A is the "physical" or rather "background" diffusivity.
*     Boundary conditions are taken into account in the Main Program.
*     Additional crosswind and induced by vertical convection diffusion.

      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)
      dimension T(0:il,0:jl,kl), KT(6,13), 
     *          hz(kl)
	real a11(6),a12(6),a22(6)
	real a13(6),a23(6),a33(6)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
	common /consts/ c3,c6,c9,c12,c18,Rhx,R2hx2
      real KT

	TAU = 0.5*dt

      C1  = 0.25                 ! crosswind diffusion a=C1*C*u

c     Vertical diffusion - influence on horizontal diffusion
      Depth=1.e4         ! Vertical Length Scale
      C2= 0. !!!(Rhx/Depth)**2


	Rhy=Rhx/asr
	km=k-1
	kp=k+1	
      DIFT=0.
      DIFTP=0.
      Q1K=T(i,j,k)+T(i+1,j,k)+T(i+1,j-1,k)
      Q2K=T(i,j,k)+T(i+1,j-1,k)+T(i,j-1,k)
      Q3K=T(i,j,k)+T(i-1,j,k)+T(i,j-1,k)
      Q4K=T(i,j,k)+T(i-1,j,k)+T(i-1,j+1,k)
      Q5K=T(i,j,k)+T(i,j+1,k)+T(i-1,j+1,k)
      Q6K=T(i,j,k)+T(i+1,j,k)+T(i,j+1,k)

      IF( K .GT. 1) THEN
C     UPPER HALF

      Q1Km=T(i,j,km)+T(i+1,j,km)+T(i+1,j-1,km)
      Q2Km=T(i,j,km)+T(i+1,j-1,km)+T(i,j-1,km)
      Q3Km=T(i,j,km)+T(i-1,j,km)+T(i,j-1,km)
      Q4Km=T(i,j,km)+T(i-1,j,km)+T(i-1,j+1,km)
      Q5Km=T(i,j,km)+T(i,j+1,km)+T(i-1,j+1,km)
      Q6Km=T(i,j,km)+T(i+1,j,km)+T(i,j+1,km)

c     Triangle # 1.

      u1=c6*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)+
     +    u(i,j,km)+u(i+1,j,km)+u(i+1,j-1,km))
      v1=c6*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k)+
     +    v(i,j,km)+v(i+1,j,km)+v(i+1,j-1,km))
      w1=c3*(w(i,j,k)+w(i+1,j,k)+w(i+1,j-1,k))

      a_11= TAU*u1*u1
      a_12= TAU*u1*v1 
      a_22= TAU*v1*v1
      a13(1)= TAU*u1*w1
      a23(1)= TAU*w1*v1 
      a33(1)= TAU*w1*w1


c     Crosswind and background diffusion
	a11(1)= a_11 +C1*a_22 +C2*a33(1) +a
	a12(1)= a_12 -C1*a_12
	a22(1)= a_22 +C1*a_11 +C2*a33(1) +a

c     -------------------------------

c     Triangle # 2.

      u2=c6*(u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k)+
     +    u(i,j,km)+u(i,j-1,km)+u(i+1,j-1,km))
      v2=c6*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k)+
     +    v(i,j,km)+v(i,j-1,km)+v(i+1,j-1,km))
      w2=c3*(w(i,j,k)+w(i,j-1,k)+w(i+1,j-1,k))

      a_11= TAU*u2*u2
      a_12= TAU*u2*v2 
      a_22= TAU*v2*v2
      a13(2)= TAU*u2*w2
      a23(2)= TAU*w2*v2 
      a33(2)= TAU*w2*w2

c     Crosswind and background diffusion
	a11(2)= a_11 +C1*a_22 +C2*a33(2) +a
	a12(2)= a_12 -C1*a_12
	a22(2)= a_22 +C1*a_11 +C2*a33(2) +a

c     -------------------------------

c     Triangle # 3.

      u3=c6*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k)+
     +    u(i,j,km)+u(i-1,j,km)+u(i,j-1,km))
      v3=c6*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k)+
     +    v(i,j,km)+v(i-1,j,km)+v(i,j-1,km))
      w3=c3*(w(i,j,k)+w(i-1,j,k)+w(i,j-1,k))

      a_11= TAU*u3*u3
      a_12= TAU*u3*v3 
      a_22= TAU*v3*v3
      a13(3)= TAU*u3*w3
      a23(3)= TAU*w3*v3 
      a33(3)= TAU*w3*w3


c     Crosswind and background diffusion
	a11(3)= a_11 +C1*a_22 +C2*a33(3) +a
	a12(3)= a_12 -C1*a_12
	a22(3)= a_22 +C1*a_11 +C2*a33(3) +a

c     -------------------------------

c     Triangle # 4.

      u4=c6*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k)+
     +    u(i,j,km)+u(i-1,j,km)+u(i-1,j+1,km))
      v4=c6*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k)+
     +    v(i,j,km)+v(i-1,j,km)+v(i-1,j+1,km))
      w4=c3*(w(i,j,k)+w(i-1,j,k)+w(i-1,j+1,k))

      a_11= TAU*u4*u4
      a_12= TAU*u4*v4 
      a_22= TAU*v4*v4
      a13(4)= TAU*u4*w4
      a23(4)= TAU*w4*v4 
      a33(4)= TAU*w4*w4

c     Crosswind and background diffusion
	a11(4)= a_11 +C1*a_22 +C2*a33(4) +a
	a12(4)= a_12 -C1*a_12
	a22(4)= a_22 +C1*a_11 +C2*a33(4) +a

c     -------------------------------

c     Triangle # 5.

      u5=c6*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k)+
     +    u(i,j,km)+u(i-1,j+1,km)+u(i,j+1,km))
      v5=c6*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)+
     +    v(i,j,km)+v(i-1,j+1,km)+v(i,j+1,km))
      w5=c3*(w(i,j,k)+w(i-1,j+1,k)+w(i,j+1,k))

      a_11= TAU*u5*u5
      a_12= TAU*u5*v5 
      a_22= TAU*v5*v5
      a13(5)= TAU*u5*w5
      a23(5)= TAU*w5*v5 
      a33(5)= TAU*w5*w5

c     Crosswind and background diffusion
	a11(5)= a_11 +C1*a_22 +C2*a33(5) +a
	a12(5)= a_12 -C1*a_12
	a22(5)= a_22 +C1*a_11 +C2*a33(5) +a

c     -------------------------------

c     Triangle # 6.

      u6=c6*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k)+
     +    u(i,j,km)+u(i+1,j,km)+u(i,j+1,km))
      v6=c6*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k)+
     +    v(i,j,km)+v(i+1,j,km)+v(i,j+1,km))
      w6=c3*(w(i,j,k)+w(i+1,j,k)+w(i,j+1,k))

      a_11= TAU*u6*u6
      a_12= TAU*u6*v6 
      a_22= TAU*v6*v6
      a13(6)= TAU*u6*w6
      a23(6)= TAU*w6*v6 
      a33(6)= TAU*w6*w6

c     Crosswind and background diffusion
	a11(6)= a_11 +C1*a_22 +C2*a33(6) +a
	a12(6)= a_12 -C1*a_12
	a22(6)= a_22 +C1*a_11 +C2*a33(6) +a

c     -------------------------------

c     Triangle # 1.

c     Horizontal part of the operator
      DIFT= DIFT +KT(1,n)*0.5*(
     &     a11(1)*(T(I+1,J,K)-T(I  ,J,  K))/S0
     &+asr*a12(1)*(T(i+1,j,k)-T(i+1,j-1,k))   )*hz(km)  

c     dT/dz*dfi/dlambda
      DIFT= DIFT +KT(1,n)*Rhx*c6*a13(1)*(Q1K-Q1Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(1,n)*Rhx*c12*a13(1)*(T(i+1,j,k )-T(i,j,k)+
     &                                   T(i+1,j,km)-T(i,j,km))

c     dfi/dz*dT/dteta
      DIFT=DIFT-KT(1,n)*Rhy*c12*S1*a23(1)*(T(i+1,j,k )-T(i+1,j-1,k)+
     &                                     T(i+1,j,km)-T(i+1,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(1,n)*R2hx2*S1*a33(1)*(Q1K-Q1Km)/hz(km)

c     Triangle # 2.

c     Horizontal part of the operator
      DIFT= DIFT -KT(2,n)*0.5*(
     &asr2*a22(2)*S2*(T(I,J,K)-T(I,J-1,K))
     &+asr*a12(2)*   (T(i+1,j-1,k)-T(i,j-1,k)))*hz(km)


c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(2,n)*Rhx*c12*a13(2)*(T(i+1,j-1,k )-T(i,j-1,k)+
     &                                    T(i+1,j-1,km)-T(i,j-1,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT -KT(2,n)*Rhy*S2*c6*a23(2)*asr*(Q2K-Q2Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(2,n)*Rhy*c12*S2*a23(2)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,km)-T(i,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(2,n)*R2hx2*S2*a33(2)*(Q2K-Q2Km)/hz(km)

c     Triangle # 3.

c     Horizontal part of the operator
      DIFT= DIFT +KT(3,n)*0.5*(
     &         a11(3)*(T(I-1,J,K)-T(I,J,K))/S0
     &    +asr*a12(3)*(T(i,j-1,k)-T(i  ,j,k))
     &    -asr*a12(3)*(T(i,j  ,k)-T(i-1,j,k))
     &-S3*asr2*a22(3)*(T(I,J  ,K)-T(I,J-1,K)))*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT -KT(3,n)*Rhx*c6*a13(3)*(Q3K-Q3Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(3,n)*Rhx*c12*a13(3)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,km)-T(i-1,j,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT -KT(3,n)*Rhy*S3*c6*a23(3)*asr*(Q3K-Q3Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(3,n)*Rhy*c12*S3*a23(3)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,km)-T(i,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(3,n)*R2hx2*S3*a33(3)*(Q3K-Q3Km)/hz(km)
c      write(*,*) 'point 4'

c     Triangle # 4.

c     Horizontal part of the operator
      DIFT= DIFT +KT(4,n)*0.5*(
     &     a11(4)*(T(I-1,J,K)-T(I,J,K))/S0
     &+asr*a12(4)*(T(i-1,j,k)-T(i-1,j+1,k)))*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT -KT(4,n)*Rhx*c6*a13(4)*(Q4K-Q4Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(4,n)*Rhx*c12*a13(4)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,km)-T(i-1,j,km))

c     dfi/dz*dT/dteta
      DIFT=DIFT-KT(4,n)*Rhy*c12*S4*a23(4)*(T(i-1,j+1,k )-T(i-1,j,k)+
     &                                     T(i-1,j+1,km)-T(i-1,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(4,n)*R2hx2*S4*a33(4)*(Q4K-Q4Km)/hz(km)

c     Triangle # 5.

c     Horizontal part of the operator
      DIFT= DIFT -KT(5,n)*0.5*(
     &asr2*S5*a22(5)*(T(I,J,K)-T(I,J+1,K))
     &   +asr*a12(5)*(T(i-1,j+1,k)-T(i,j+1,k)))*hz(km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(5,n)*Rhx*c12*a13(5)*(T(i,j+1,k )-T(i-1,j+1,k)+
     &                                    T(i,j+1,km)-T(i-1,j+1,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT +KT(5,n)*Rhy*S5*c6*a23(5)*asr*(Q5K-Q5Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(5,n)*Rhy*c12*S5*a23(5)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,km)-T(i,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(5,n)*R2hx2*S5*a33(5)*(Q5K-Q5Km)/hz(km)

c     Triangle # 6.

c     Horizontal part of the operator
      DIFT= DIFT +KT(6,n)*0.5*(
     &         a11(6)*(T(I+1,J,K)-T(I,J,K))/S0
     &    +asr*a12(6)*(T(i,j+1,k)-T(i,j,k))   
     &    -asr*a12(6)*(T(i,j,k)-T(i+1,j,k))
     *-asr2*S6*a22(6)*(T(I,J,K)-T(I,J+1,K)))*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT +KT(6,n)*Rhx*c6*a13(6)*(Q6K-Q6Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(6,n)*Rhx*c12*a13(6)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,km)-T(i,j,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT +KT(6,n)*Rhy*S6*c6*a23(6)*asr*(Q6K-Q6Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(6,n)*Rhy*c12*S6*a23(6)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,km)-T(i,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(6,n)*R2hx2*S6*a33(6)*(Q6K-Q6Km)/hz(km)

      END IF

      IF( K .LT. KB) THEN 
C     LOWER HALF

      Q1Kp=T(i,j,kp)+T(i+1,j,kp)+T(i+1,j-1,kp)
      Q2Kp=T(i,j,kp)+T(i+1,j-1,kp)+T(i,j-1,kp)
      Q3Kp=T(i,j,kp)+T(i-1,j,kp)+T(i,j-1,kp)
      Q4Kp=T(i,j,kp)+T(i-1,j,kp)+T(i-1,j+1,kp)
      Q5Kp=T(i,j,kp)+T(i,j+1,kp)+T(i-1,j+1,kp)
      Q6Kp=T(i,j,kp)+T(i+1,j,kp)+T(i,j+1,kp)


c     Triangle # 1.

      u1=c6*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)+
     +    u(i,j,kp)+u(i+1,j,kp)+u(i+1,j-1,kp))
      v1=c6*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k)+
     +    v(i,j,kp)+v(i+1,j,kp)+v(i+1,j-1,kp))
      w1=c3*(w(i,j,kp)+w(i+1,j,kp)+w(i+1,j-1,kp))

      a_11= TAU*u1*u1
      a_12= TAU*u1*v1 
      a_22= TAU*v1*v1
      a13(1)= TAU*u1*w1
      a23(1)= TAU*w1*v1 
      a33(1)= TAU*w1*w1

c     Crosswind and background diffusion
	a11(1)= a_11 +C1*a_22 +C2*a33(1) +a
	a12(1)= a_12 -C1*a_12
	a22(1)= a_22 +C1*a_11 +C2*a33(1) +a

c     -------------------------------

c     Triangle # 2.

      u2=c6*(u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k)+
     +    u(i,j,kp)+u(i,j-1,kp)+u(i+1,j-1,kp))
      v2=c6*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k)+
     +    v(i,j,kp)+v(i,j-1,kp)+v(i+1,j-1,kp))
      w2=c3*(w(i,j,kp)+w(i,j-1,kp)+w(i+1,j-1,kp))

      a_11= TAU*u2*u2
      a_12= TAU*u2*v2 
      a_22= TAU*v2*v2
      a13(2)= TAU*u2*w2
      a23(2)= TAU*w2*v2 
      a33(2)= TAU*w2*w2

c     Crosswind and background diffusion
	a11(2)= a_11 +C1*a_22 +C2*a33(2) +a
	a12(2)= a_12 -C1*a_12
	a22(2)= a_22 +C1*a_11 +C2*a33(2) +a

c     -------------------------------

c     Triangle # 3.

      u3=c6*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k)+
     +    u(i,j,kp)+u(i-1,j,kp)+u(i,j-1,kp))
      v3=c6*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k)+
     +    v(i,j,kp)+v(i-1,j,kp)+v(i,j-1,kp))
      w3=c3*(w(i,j,kp)+w(i-1,j,kp)+w(i,j-1,kp))

      a_11= TAU*u3*u3
      a_12= TAU*u3*v3 
      a_22= TAU*v3*v3
      a13(3)= TAU*u3*w3
      a23(3)= TAU*w3*v3 
      a33(3)= TAU*w3*w3

c     Crosswind and background diffusion
	a11(3)= a_11 +C1*a_22 +C2*a33(3) +a
	a12(3)= a_12 -C1*a_12
	a22(3)= a_22 +C1*a_11 +C2*a33(3) +a

c     -------------------------------

c     Triangle # 4.

      u4=c6*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k)+
     +    u(i,j,kp)+u(i-1,j,kp)+u(i-1,j+1,kp))
      v4=c6*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k)+
     +    v(i,j,kp)+v(i-1,j,kp)+v(i-1,j+1,kp))
      w4=c3*(w(i,j,kp)+w(i-1,j,kp)+w(i-1,j+1,kp))

      a_11= TAU*u4*u4
      a_12= TAU*u4*v4 
      a_22= TAU*v4*v4
      a13(4)= TAU*u4*w4
      a23(4)= TAU*w4*v4 
      a33(4)= TAU*w4*w4

c     Crosswind and background diffusion
	a11(4)= a_11 +C1*a_22 +C2*a33(4) +a
	a12(4)= a_12 -C1*a_12
	a22(4)= a_22 +C1*a_11 +C2*a33(4) +a

c     -------------------------------

c     Triangle # 5.

      u5=c6*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k)+
     +    u(i,j,kp)+u(i-1,j+1,kp)+u(i,j+1,kp))
      v5=c6*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)+
     +    v(i,j,kp)+v(i-1,j+1,kp)+v(i,j+1,kp))
      w5=c3*(w(i,j,kp)+w(i-1,j+1,kp)+w(i,j+1,kp))

      a_11= TAU*u5*u5
      a_12= TAU*u5*v5 
      a_22= TAU*v5*v5
      a13(5)= TAU*u5*w5
      a23(5)= TAU*w5*v5 
      a33(5)= TAU*w5*w5

c     Crosswind and background diffusion
	a11(5)= a_11 +C1*a_22 +C2*a33(5) +a
	a12(5)= a_12 -C1*a_12
	a22(5)= a_22 +C1*a_11 +C2*a33(5) +a

c     -------------------------------

c     Triangle # 6.

      u6=c6*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k)+
     +    u(i,j,kp)+u(i+1,j,kp)+u(i,j+1,kp))
      v6=c6*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k)+
     +    v(i,j,kp)+v(i+1,j,kp)+v(i,j+1,kp))
      w6=c3*(w(i,j,kp)+w(i+1,j,kp)+w(i,j+1,kp))

      a_11= TAU*u6*u6
      a_12= TAU*u6*v6 
      a_22= TAU*v6*v6
      a13(6)= TAU*u6*w6
      a23(6)= TAU*w6*v6 
      a33(6)= TAU*w6*w6

c     Crosswind and background diffusion
	a11(6)= a_11 +C1*a_22 +C2*a33(6) +a
	a12(6)= a_12 -C1*a_12
	a22(6)= a_22 +C1*a_11 +C2*a33(6) +a

c     -------------------------------

c     Triangle # 1.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(1,np)*(a11(1)*(T(I+1,J,K)-T(I,J,K))/S0
     &            +asr*a12(1)*(T(i+1,j,k)-T(i+1,j-1,k)))*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP +KT(1,np)*Rhx*c6*a13(1)*(Q1Kp-Q1K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(1,np)*Rhx*c12*a13(1)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,kp)-T(i,j,kp))

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(1,np)*Rhy*c12*S1*a23(1)*(T(i+1,j,k )-T(i+1,j-1,k)+
     &                                       T(i+1,j,kp)-T(i+1,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(1,np)*R2hx2*S1*a33(1)*(Q1K-Q1Kp)/hz(k)

c     Triangle # 2.

c     Horizontal part of the operator
      DIFTP= DIFTP-0.5*KT(2,np)*(a22(2)*S2*(T(I,J,K)-T(I,J-1,K))*ASR2
     &           +        asr*a12(2)*(T(i+1,j-1,k)-T(i,j-1,k)))*hz(k)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(2,np)*Rhx*c12*a13(2)*
     &                                   (T(i+1,j-1,k )-T(i,j-1,k)+
     &                                    T(i+1,j-1,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP -KT(2,np)*Rhy*S2*c6*a23(2)*asr*(Q2Kp-Q2K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(2,np)*Rhy*c12*S2*a23(2)*
     &                                      (T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(2,np)*R2hx2*S2*a33(2)*(Q2K-Q2Kp)/hz(k)

c     Triangle # 3.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(3,np)*(a11(3)*(T(I-1,J,K)-T(I,J,K))/S0
     &            +asr*a12(3)*(T(i,j-1,k)-T(i,j,k))
     &            -asr*a12(3)*(T(i,j,k)-T(i-1,j,k))
     &             -S3*a22(3)*(T(I,J,K)-T(I,J-1,K))*ASR2)*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP -KT(3,np)*Rhx*c6*a13(3)*(Q3Kp-Q3K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(3,np)*Rhx*c12*a13(3)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,kp)-T(i-1,j,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP -KT(3,np)*Rhy*S3*c6*a23(3)*asr*(Q3Kp-Q3K)

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(3,np)*Rhy*c12*S3*a23(3)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(3,np)*R2hx2*S3*a33(3)*(Q3K-Q3Kp)/hz(k)

c     Triangle # 4.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(4,np)*(a11(4)*(T(I-1,J,K)-T(I,J,K))/S0
     &            +asr*a12(4)*(T(i-1,j,k)-T(i-1,j+1,k)))*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP -KT(4,np)*Rhx*c6*a13(4)*(Q4Kp-Q4K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(4,np)*Rhx*c12*a13(4)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,kp)-T(i-1,j,kp))

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(4,np)*Rhy*c12*S4*a23(4)*
     &                                  (T(i-1,j+1,k )-T(i-1,j,k)+
     &                                   T(i-1,j+1,kp)-T(i-1,j,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(4,np)*R2hx2*S4*a33(4)*(Q4K-Q4Kp)/hz(k)


c     Triangle # 5.

c     Horizontal part of the operator
      DIFTP= DIFTP-0.5*KT(5,np)*(a22(5)*S5*(T(I,J,K)-T(I,J+1,K))*ASR2
     &           +        asr*a12(5)*(T(i-1,j+1,k)-T(i,j+1,k)))*hz(k)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(5,np)*Rhx*c12*a13(5)*
     &                                   (T(i,j+1,k )-T(i-1,j+1,k)+
     &                                    T(i,j+1,kp)-T(i-1,j+1,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP +KT(5,np)*Rhy*S5*c6*a23(5)*asr*(Q5Kp-Q5K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(5,np)*Rhy*c12*S5*a23(5)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,kp)-T(i,j,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(5,np)*R2hx2*S5*a33(5)*(Q5K-Q5Kp)/hz(k)


c     Triangle # 6.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(6,np)*(a11(6)*(T(I+1,J,K)-T(I,J,K))/S0
     &               +a12(6)*asr*(T(i,j+1,k)-T(i,j,k))   
     &               -a12(6)*asr*(T(i,j,k)-T(i+1,j,k))
     *               -S6*a22(6)*(T(I,J,K)-T(I,J+1,K))*ASR2)*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP +KT(6,np)*Rhx*c6*a13(6)*(Q6Kp-Q6K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(6,np)*Rhx*c12*a13(6)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,kp)-T(i,j,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP +KT(6,np)*Rhy*S6*c6*a23(6)*asr*(Q6Kp-Q6K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(6,np)*Rhy*c12*S6*a23(6)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,kp)-T(i,j,kp))

c     dT/dz*dfi/dz 
      DIFTP= DIFTP -KT(6,np)*R2hx2*S6*a33(6)*(Q6K-Q6Kp)/hz(k)

      END IF

      RETURN
      END
