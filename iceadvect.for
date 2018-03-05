      SUBROUTINE iceadvect(m)
*     *****************************************
*     Ice advection by Upwinding Scheme + Crosswind.
*     Matsuno time scheme.
*     Special treatment of open boundaries.
*     Version 26.03.2009.
*     *****************************************

      INCLUDE 'slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
	Include 'tparm.fi'

	Adiff= 5.0e6 ! Background diffusivity to suppress oscillations.
                   ! This oscillations occur primarility at solid bndrs.

      c3 =1./3.
	c9 =c3*c3
	c6= 1./6.
	c12=0.5*c6
	c18= c6/3.
	c24=0.5*c12
      asr= hx/hy
	asr2=asr*asr
      HXX= 1./HX/R
      COEFX=HXX/18.
	HYY= 1./HY/R
      COEFY=HYY/18.
      A= 0.5*HXX*HXX

*     Divergence.
	do j=1,jl
      do i=1,il
	if(km2(i,j).GT.0) then
      div_ice_tr(i,j)= 0.0
      do mg=1,mgrad
      div_ice_tr(i,j)= div_ice_tr(i,j)+aice2(mg,i,j)
      end do 
	end if
	end do
	end do

*     Internal points and solid boundaries.

      DO 2 J=1,JL
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4

      DO 2 I=1,IL
      ni= nt3(i,j,1)

      IF( ni.GT.0 .OR. (ni.LT.0.AND.ice_mask(i,j).EQ.0) ) then

	n=abs(ni)
      CA=6.*dt/(CG(n)*S0)

	do mg=1,mgrad

      IF(m.EQ.1) then
	call transp_ice(FXa,FYa,i,j,mg,n,uice,vice,aice2,KT,
     *                il1,jl1,mgrad)
      FXa =COEFX*FXa
      FYa =COEFY*FYa

      call transp_ice(FXi,FYi,i,j,mg,n,uice,vice,Hice2,KT,
     *                il1,jl1,mgrad)
      FXi =COEFX*FXi
      FYi =COEFY*FYi

      call transp_ice(FXs,FYs,i,j,mg,n,uice,vice,Hsnow2,KT,
     *                il1,jl1,mgrad)
      FXs =COEFX*FXs
      FYs =COEFY*FYs

	else

	call transp_ice(FXa,FYa,i,j,mg,n,uice,vice,aice1,KT,
     *                il1,jl1,mgrad)
      FXa =COEFX*FXa
      FYa =COEFY*FYa

      call transp_ice(FXi,FYi,i,j,mg,n,uice,vice,Hice1,KT,
     *                il1,jl1,mgrad)
      FXi =COEFX*FXi
      FYi =COEFY*FYi

      call transp_ice(FXs,FYs,i,j,mg,n,uice,vice,Hsnow1,KT,
     *                il1,jl1,mgrad)
      FXs =COEFX*FXs
      FYs =COEFY*FYs

      end if

      call DIFFUW_ICE(Aice2, uice,vice,i,j,mg,n,DIFFa,il1,jl1,
     &                mgrad,KT,Adiff,hx,hy,r)
      call DIFFUW_ICE(Hice2, uice,vice,i,j,mg,n,DIFFi,il1,jl1,
     &                mgrad,KT,Adiff,hx,hy,r)
      call DIFFUW_ICE(Hsnow2,uice,vice,i,j,mg,n,DIFFs,il1,jl1,
     &                mgrad,KT,Adiff,hx,hy,r)

      Fa=-FXa-FYa+a*DIFFa
      Fi=-FXi-FYi+a*DIFFi
      Fs=-FXs-FYs+a*DIFFs

      Aice (MG,i,j)=Aice2(MG,i,j) +CA*Fa
      Hice (MG,i,j)=Hice2(MG,i,j) +CA*Fi
      Hsnow(MG,i,j)=Hsnow2(MG,i,j)+CA*Fs

	end do  ! mgrad

      end if    ! n>0
2     continue


*     Open boundaries, n<0.
	call
     $ icebc(Aice,Aice2,uice,vice,nt3,Si,ice_mask,
     &       mgrad,il1,jl1,kl,R,hx,hy,dt)
	call
     $ icebc(hsnow,hsnow2,uice,vice,nt3,Si,ice_mask,
     &       mgrad,il1,jl1,kl,R,hx,hy,dt)
	call 
     & icebc(hice,hice2,uice,vice,nt3,Si,ice_mask,
     &       mgrad,il1,jl1,kl,R,hx,hy,dt)

c     Angle Points Special treatment.
c     Cutting off negative concentrations.

      do j=1,jl
      do i=1,il
      ni= abs(nt3(i,j,1))

c      if (ni.eq.13) then
c      do mg=1,mgrad
c      Hice (mg,i,j)= Hice (mg,i+1,j+1)
c      Aice (mg,i,j)= Aice (mg,i+1,j+1)
c      Hsnow(mg,i,j)= Hsnow(mg,i+1,j+1)
c      end do
c      end if

c      if (ni.eq.11) then
c      do mg=1,mgrad
c      Hice (mg,i,j)= Hice (mg,i-1,j-1)
c      Aice (mg,i,j)= Aice (mg,i-1,j-1)
c      Hsnow(mg,i,j)= Hsnow(mg,i-1,j-1)
c      end do
c      end if

      aice(0,i,j)=1.
	do mg=1,mgrad

	if(aice(mg,i,j).LT.aimin .OR. hice(mg,i,j).LT.himin)then
ccc	write(*,*) 'Negative ice:', mg,i,j,aice(mg,i,j),hice(mg,i,j)
	aice(mg,i,j) =0.
	hice(mg,i,j) =0.
	hsnow(mg,i,j)=0.
	end if

	aice(0,i,j)=aice(0,i,j)-aice(mg,i,j)

	if(hsnow(mg,i,j).LT.hsmin) hsnow(mg,i,j)=0.	
	end do
     
      aice(0,i,j)= MAX(0.,aice(0,i,j))

      end do
      end do

*     Divergence.
	do j=1,jl
      do i=1,il
	if(km2(i,j).GT.0) then
      do mg=1,mgrad
      div_ice_tr(i,j)= div_ice_tr(i,j)-aice(mg,i,j)
      end do 
	div_ice_tr(i,j)=div_ice_tr(i,j)/dt
	end if
	end do
	end do

      return
      end

      subroutine transp_ice(FXTK,FYTK,i,j,mg,n,u,v,t,KT,il,jl,mgrad)
      dimension u(0:il,0:jl),v(0:il,0:jl),t(0:mgrad,0:il,0:jl),KT(6,13)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

*     Version 10.10.2008.

      FXTK= 0.
      FYTK= 0.

*---------------     Triangle # 1.  ---------------------------------
      FXTK=FXTK+(U(I,J)+U(I+1,J)+U(I+1,J-1))*
     *          (T(mg,i  ,J  )+T(mg,i+1,J  )+T(mg,i+1,J-1))
     *  *KT(1,n)
*---------------     Triangle # 2.  -----------------------------------
      FYTK=FYTK-S2*(V(I,J)+V(I,J-1)+V(I+1,J-1))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J-1)+T(mg,i+1,J-1))
     *  *KT(2,n)
*---------------     Triangle # 3.  -------------------------------
      FXTK=FXTK-(U(I,J)+U(I-1,J)+U(I,J-1))*
     *          (T(mg,i  ,J  )+T(mg,i-1,J  )+T(mg,i  ,J-1))
     *  *KT(3,n)
      FYTK=FYTK-S3*(V(I,J)+V(I,J-1)+V(I-1,J))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J-1)+T(mg,i-1,J  ))
     *  *KT(3,n)
*---------------     Triangle # 4.  ---------------------------------
      FXTK=FXTK-(U(I,J)+U(I-1,J)+U(I-1,J+1))*
     *          (T(mg,i  ,J  )+T(mg,i-1,J  )+T(mg,i-1,J+1))
     *  *KT(4,n)
*---------------     Triangle # 5.  ----------------------------------
      FYTK=FYTK+S5*(V(I,J)+V(I,J+1)+V(I-1,J+1))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J+1)+T(mg,i-1,J+1))
     *  *KT(5,n)
*---------------     Triangle # 6.  -----------------------------------
      FXTK=FXTK+(U(I,J)+U(I+1,J)+U(I,J+1))*
     *          (T(mg,i  ,J   )+T(mg,i+1,J  )+T(mg,i  ,J+1))
     *  *KT(6,n)
      FYTK=FYTK+S6*(V(I,J)+V(I,J+1)+V(I+1,J))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J+1)+T(mg,i+1,J))
     *  *KT(6,n)

      return
      end

	subroutine icebc
     &(hice,hice2,u,v,nt3,Si,ice_mask,mgrad,il1,jl1,kl,R,hx,hy,dt)

*     Only transport as a phase velocity
*     19.04.2006.

	dimension hice(0:mgrad,0:il1,0:jl1),hice2(0:mgrad,0:il1,0:jl1),
     &          nt3(0:il1,0:jl1,kl), Si(0:jl1),
     &          ice_mask(0:il1,0:jl1),
     &          u(0:il1,0:jl1),v(0:il1,0:jl1)

	do j=1,jl1-1
	S0=Si(j)
      do i=1,il1-1

	IF(ice_mask(i,j).EQ.1) then
	ni=nt3(i,j,1)

	do mg=1,mgrad
	if(ni.EQ.-6) then

	Rx=-dt*u(i+1,j)/(R*hx*S0)
	Ry=+dt*v(i+1,j)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i+1,j)
     $                             +Ry*(hice2(mg,i,j-1)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i+1,j)
     $                             -Ry*(hice2(mg,i,j+1)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
	end if  ! n=-6

	if(ni.EQ.-7) then

	Rx=-dt*v(i,j+1)/(R*hy)
	Ry=-dt*u(i,j+1)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j+1)
     $                             +Ry*(hice2(mg,i-1,j)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j+1)
     $                             -Ry*(hice2(mg,i+1,j)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if

	end if  ! n=-7

	if(ni.EQ.-8) then

	Rx=+dt*u(i-1,j)/(R*S0*hx)
	Ry=-dt*v(i-1,j)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i-1,j)
     $                             +Ry*(hice2(mg,i,j-1)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i-1,j)
     $                             -Ry*(hice2(mg,i,j+1)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
	end if  ! n=-8

	if(ni.EQ.-9) then

	Rx=+dt*v(i,j-1)/(R*hy)
	Ry=+dt*u(i,j-1)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j-1)
     &                             +Ry*(hice2(mg,i-1,j)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j-1)
     &                             -Ry*(hice2(mg,i+1,j)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
	end if  ! n=-9

	if(ni.EQ.-10) then

	if(nt3(i+1,j,1).GT.0.AND.nt3(i,j-1,1).LT.0)then !transport to W

	Rx=-dt*u(i+1,j)/(R*hx*S0)
	Ry=+dt*v(i+1,j)/(R*hy)    

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i+1,j)
     $                             +Ry*(hice2(mg,i,j-1)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i+1,j))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if

      end if  ! To W

	if(nt3(i+1,j,1).LT.0.AND.nt3(i,j-1,1).GT.0)then !transport to S

	Rx=+dt*v(i,j-1)/(R*hy)
	Ry=+dt*u(i,j-1)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j-1))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j-1)
     &                             -Ry*(hice2(mg,i+1,j)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
      end if  ! To S

	end if  ! n=-10

	if(ni.EQ.-11) then

	if(nt3(i-1,j,1).GT.0.AND.nt3(i,j-1,1).LT.0)then !transport to E

	Rx=+dt*u(i-1,j)/(R*S0*hx)
	Ry=-dt*v(i-1,j)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i-1,j)
     $                             +Ry*(hice2(mg,i,j-1)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i-1,j))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
      
	end if ! To E

	if(nt3(i-1,j,1).LT.0.AND.nt3(i,j-1,1).GT.0)then !transport to S

	Rx=+dt*v(i,j-1)/(R*hy)
	Ry=+dt*u(i,j-1)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j-1)
     &                             +Ry*(hice2(mg,i-1,j)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j-1))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
      
	end if  ! To S

	end if  ! n=-11

	if(ni.EQ.-12) then

	if(nt3(i-1,j,1).GT.0.AND.nt3(i,j+1,1).LT.0)then !transport to E

	Rx=+dt*u(i-1,j)/(R*S0*hx)
	Ry=-dt*v(i-1,j)/(R*hy) 

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i-1,j))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i-1,j)
     $                             -Ry*(hice2(mg,i,j+1)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
      
	end if   ! To E

	if(nt3(i-1,j,1).LT.0.AND.nt3(i,j+1,1).GT.0)then !transport to N

	Rx=-dt*v(i,j+1)/(R*hy)
	Ry=-dt*u(i,j+1)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j+1)
     $                             +Ry*(hice2(mg,i-1,j)-hice2(mg,i,j)))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j+1))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
      
	end if  ! To N

	end if  ! n=-12

	if(ni.EQ.-13) then

	if(nt3(i+1,j,1).GT.0.AND.nt3(i,j+1,1).LT.0)then !transport to W

	Rx=-dt*u(i+1,j)/(R*hx*S0)
	Ry=+dt*v(i+1,j)/(R*hy)     

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i+1,j))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i+1,j)
     $                             -Ry*(hice2(mg,i,j+1)-hice2(mg,i,j)))
	end if
	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if

      end if

	if(nt3(i+1,j,1).LT.0.AND.nt3(i,j+1,1).GT.0)then !transport to N

	Rx=-dt*v(i,j+1)/(R*hy)
	Ry=-dt*u(i,j+1)/(R*hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j+1))
	else
	hice(mg,i,j)=D*(hice2(mg,i,j)+Rx*hice(mg,i,j+1)
     $                             -Ry*(hice2(mg,i+1,j)-hice2(mg,i,j)))
	end if

	else
      hice(mg,i,j)=hice2(mg,i,j)
	end if
      
	end if  ! To N

	end if  ! n=-13

      end do  ! mgrad
      end if  ! mask=1
	end do  ! j
	end do  ! i

	return
	end


      SUBROUTINE DIFFUW_ICE(T,u,v,i,j,mg,n,DIFT,il,jl,mgrad,
     &                      KT,A,hx,hy,R)
*     Version 10.10.08
*     Shock capture 
*     Boundary conditions are taken into account in the Main Program.

      dimension u(0:il,0:jl), v(0:il,0:jl)
      dimension T(0:mgrad,0:il,0:jl), KT(6,13)
	dimension a11(6),a12(6),a22(6)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      real KT


      Upwind= 0.25

	c3=1./3.
      DIFT=0.

      CL= SQRT(0.5*hx*s0 *hy)  ! Length scale in radians.
      C=  R*CL                 ! Length scale.
	wmin=1.e-12

      u1=c3*(u(i,j)+u(i+1,j)+u(i+1,j-1))
      v1=c3*(v(i,j)+v(i+1,j)+v(i+1,j-1))

      gradu= (T(mg,i+1,j)-T(mg,i,j))/s0
	gradv= (T(mg,i+1,j)-T(mg,i+1,j-1))*asr

	g_norm= gradu**2 +gradv**2
	IF(g_norm .GT. wmin) then
	scalprod=(gradu*u1 +gradv*v1)/g_norm
	ushock= scalprod*gradu
	vshock= scalprod*gradv
	else
	ushock=u1
	vshock=v1
	end if
      w1=MAX(wmin,SQRT(ushock**2 +vshock**2))
	PeH= 0.5*w1*C/a
	alpha= MIN( 0.333*PeH, 1.0)
	C2= Upwind*alpha*C

	a11(1)= C2*ushock*ushock/w1 +a
	a12(1)= C2*ushock*vshock/w1 
	a22(1)= C2*vshock*vshock/w1 +a
c     -------------------------------

      u2=c3*(u(i,j)+u(i,j-1)+u(i+1,j-1))
      v2=c3*(v(i,j)+v(i,j-1)+v(i+1,j-1))

      gradu= (T(mg,i+1,j-1)-T(mg,i,j-1))/s0
	gradv= (T(mg,i,j)-T(mg,i,j-1))*asr

	g_norm= gradu**2 +gradv**2
	IF(g_norm .GT. wmin) then
	scalprod=(gradu*u2 +gradv*v2)/g_norm
	ushock= scalprod*gradu
	vshock= scalprod*gradv
	else
	ushock=u2
	vshock=v2
	end if
      w2=MAX(wmin,SQRT(ushock**2 +vshock**2))
	PeH= 0.5*w2*C/a
	alpha= MIN( 0.333*PeH, 1.0)
	C2= Upwind*alpha*C

	a11(2)= C2*ushock*ushock/w2 +a
	a12(2)= C2*ushock*vshock/w2 
	a22(2)= C2*vshock*vshock/w2 +a
c     -------------------------------

      u3=c3*(u(i,j)+u(i-1,j)+u(i,j-1))
      v3=c3*(v(i,j)+v(i-1,j)+v(i,j-1))

      gradu= (T(mg,i,j)-T(mg,i-1,j))/s0
	gradv= (T(mg,i,j)-T(mg,i,j-1))*asr

	g_norm= gradu**2 +gradv**2
	IF(g_norm .GT. wmin) then
	scalprod=(gradu*u3 +gradv*v3)/g_norm
	ushock= scalprod*gradu
	vshock= scalprod*gradv
	else
	ushock=u3
	vshock=v3
	end if
      w3=MAX(wmin,SQRT(ushock**2 +vshock**2))
	PeH= 0.5*w3*C/a
	alpha= MIN( 0.333*PeH, 1.0)
	C2= Upwind*alpha*C

	a11(3)= C2*ushock*ushock/w3 +a
	a12(3)= C2*ushock*vshock/w3 
	a22(3)= C2*vshock*vshock/w3 +a
c     -------------------------------

      u4=c3*(u(i,j)+u(i-1,j)+u(i-1,j+1))
      v4=c3*(v(i,j)+v(i-1,j)+v(i-1,j+1))

c     Shock capturing
      gradu= (T(mg,i,j)-T(mg,i-1,j))/s0
	gradv= (T(mg,i-1,j+1)-T(mg,i-1,j))*asr

	g_norm= gradu**2 +gradv**2
	IF(g_norm .GT. wmin) then
	scalprod=(gradu*u4 +gradv*v4)/g_norm
	ushock= scalprod*gradu
	vshock= scalprod*gradv
	else
	ushock=u4
	vshock=v4
	end if
      w4=MAX(wmin,SQRT(ushock**2 +vshock**2))
	PeH= 0.5*w4*C/a
	alpha= MIN( 0.333*PeH, 1.0)
	C2= Upwind*alpha*C

	a11(4)= C2*ushock*ushock/w4 +a
	a12(4)= C2*ushock*vshock/w4 
	a22(4)= C2*vshock*vshock/w4 +a
c     -------------------------------

      u5=c3*(u(i,j)+u(i-1,j+1)+u(i,j+1))
      v5=c3*(v(i,j)+v(i-1,j+1)+v(i,j+1))

      gradu= (T(mg,i,j+1)-T(mg,i-1,j+1))/s0
	gradv= (T(mg,i,j+1)-T(mg,i,j))*asr

	g_norm= gradu**2 +gradv**2
	IF(g_norm .GT. wmin) then
	scalprod=(gradu*u5 +gradv*v5)/g_norm
	ushock= scalprod*gradu
	vshock= scalprod*gradv
	else
	ushock=u5
	vshock=v5
	end if
      w5=MAX(wmin,SQRT(ushock**2 +vshock**2))
	PeH= 0.5*w5*C/a
	alpha= MIN( 0.333*PeH, 1.0)
	C2= Upwind*alpha*C

	a11(5)= C2*ushock*ushock/w5 +a
	a12(5)= C2*ushock*vshock/w5 
	a22(5)= C2*vshock*vshock/w5 +a
c     -------------------------------

      u6=c3*(u(i,j)+u(i+1,j)+u(i,j+1))
      v6=c3*(v(i,j)+v(i+1,j)+v(i,j+1))

      gradu= (T(mg,i+1,j)-T(mg,i,j))/s0
	gradv= (T(mg,i,j+1)-T(mg,i,j))*asr

	g_norm= gradu**2 +gradv**2
	IF(g_norm .GT. wmin) then
	scalprod=(gradu*u6 +gradv*v6)/g_norm
	ushock= scalprod*gradu
	vshock= scalprod*gradv
	else
	ushock=u6
	vshock=v6
	end if
      w6=MAX(wmin,SQRT(ushock**2 +vshock**2))
	PeH= 0.5*w6*C/a
	alpha= MIN( 0.333*PeH, 1.0)
	C2= Upwind*alpha*C

	a11(6)= C2*ushock*ushock/w6 +a
	a12(6)= C2*ushock*vshock/w6 
	a22(6)= C2*vshock*vshock/w6 +a
c     -------------------------------

c     Triangle # 1.
      DIFT= DIFT+KT(1,n)*(a11(1)*(T(mg,I+1,J)-T(mg,I,J))/S0
     &               +asr*a12(1)*(T(mg,i+1,j)-T(mg,i+1,j-1)))
c     Triangle # 2.
      DIFT= DIFT-KT(2,n)*(a22(2)*S2*(T(mg,I,J)-T(mg,I,J-1))*ASR2
     &               +asr*a12(2)*(T(mg,i+1,j-1)-T(mg,i,j-1)))
c     Triangle # 3.
      DIFT= DIFT+KT(3,n)*(a11(3)*(T(mg,I-1,J)-T(mg,I,J))/S0
     &               +asr*a12(3)*(T(mg,i,j-1)-T(mg,i,j))
     &               -asr*a12(3)*(T(mg,i,j)-T(mg,i-1,j))
     &                -S3*a22(3)*(T(mg,I,J)-T(mg,I,J-1))*ASR2)
c     Triangle # 4.
      DIFT= DIFT+KT(4,n)*(a11(4)*(T(mg,I-1,J)-T(mg,I,J))/S0
     &               +asr*a12(4)*(T(mg,i-1,j)-T(mg,i-1,j+1)))
c     Triangle # 5.
      DIFT= DIFT-KT(5,n)*(a22(5)*S5*(T(mg,I,J)-T(mg,I,J+1))*ASR2
     &               +asr*a12(5)*(T(mg,i-1,j+1)-T(mg,i,j+1)))
c     Triangle # 6.
      DIFT= DIFT+KT(6,n)*(a11(6)*(T(mg,I+1,J)-T(mg,I,J))/S0
     &               +asr*a12(6)*(T(mg,i,j+1)-T(mg,i,j))   
     &               -asr*a12(6)*(T(mg,i,j)-T(mg,i+1,j))
     &                -S6*a22(6)*(T(mg,I,J)-T(mg,I,J+1))*ASR2)

      RETURN
      END

