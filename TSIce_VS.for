      SUBROUTINE TSICE_VIS

c     Version 09.11.2012 for ice thickness gradations.
c     Potential temperature and Freezing point as function on depth.

*-----------------------------------------------------
c     Variable ice salinity
*-----------------------------------------------------

c     Nonuniform heat flux ice-ocean
c     Spatially Varied Cloudiness
c     AOMIP forcing - Net LW radiation.
c     Penetration of the Short Wave Radiation into the Ocean.
c     Friction velocity as in Los Alamos model?
c     Evaporation and snow/ice sublimation. 

**    Temperature flux due to Rain and Rivers

      INCLUDE 'slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
	REAL dHiceT, dHsnowT, SiceT
      INCLUDE 'tparm.fi'

	c3=1./3.
	
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

      do i=1,il

      dHiceT = 0.
      dHsnowT= 0.
	SiceT  = 0.

      KB = KM2(I,J)
	n=ABS(nt3(i,j,1))
*-------------------------------------------------------
*     Calculation in all points.

      IF( KB .GT. 0) THEN

*     -------------- Ocean Surface -----------------

*     Precipitation: Rain to the open water. Weatherly climate model.
	IF(TA(i,j).GT. 5.) Rain= Pr(i,j)
	IF(TA(i,j).LT.-5.) Rain= Pr(i,j)*Aice(0,i,j)
	IF(TA(i,j).GE.-5..AND.TA(i,j).LE.5.) then
	Rain= Pr(i,j)*(Aice(0,i,j)+(1.-Aice(0,i,j))*
     &                0.1*(TA(i,j)+5.) )
	end if

 
      CA=2.*DT/hz(1)

      AM(1)= 0.
      CM(1)= -CA*AZT(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) 
      FM(1)= T(i,j,1)  +CA*( Rain*TA(i,j) +River(i,j)*Tobs(i,j,1))

*     Heat fluxes through ocean and ice.
*     Wind is in m/s.
      WIND=100.*(SQRT(wx(i,j)**2+wy(i,j)**2) +E0)


*     1. Sesible heat at open water.
      BM(1)=BM(1)+CA*roa*Cpa*CDH(T(i,j,1),TA(i,j))*
     &WIND        *Aice(0,i,j)/(row*cpw)
      FM(1)=FM(1)+CA*roa*Cpa*CDH(T(i,j,1),TA(i,j))*
     &WIND*TA(i,j)*Aice(0,i,j)/(row*cpw)


*     3. Incoming radiation .
*     3.1. Short wave at open water.

	Qrad= Aice(0,i,j)*(1.-aopw)*SW(i,j)

      aI0= (0.18 +0.17*cloud(i,j))    ! Grenfell & Maykut, 1977. 
cc	aI0= 0.4*aI0


      do mg=1,mgrad
	Ai=MAX(Aimin,Aice(mg,i,j))
	Hi=Hice(mg,i,j)/Ai

	ppp= 0. !!!!!1.e-5*g*row*Hi                   
      TFC= TFr(S(i,j,1), ppp) ! Freezing point at the depth of ice bottom


*     2. Sensible heat at ice-water interface. 

*     McPhee, M.G. 1992. Turbulent heat fluxes in the upper ocean
*     under ice. JGR, 97, 5365-5379.
*     Q=rw*Cp*CDw*CoefW*dT, CoefW=SQRT(CDw*((uice-u)**2) -friction velocity.
*     Skin velocity W=1-3 cm/s, Mellor&Kantha 0.5.
*     Los Alamos sea ice model - cut off 0.075 -> dU=1cm/s
*     Drag coefficient


*     Cdw=5.5e-3 or hi dependent - Mellor&Kantha 1989
c      Hscale= 100.
c      CDw1=( 0.4/alog(Hscale*300./(max(5.,0.5*hi))) )**2

	Cdw1=Cdw

      CoefW=SQRT(CDw1*( (uice(i,j)-um2(i,j,1))**2
     &                 +(vice(i,j)-vm2(i,j,1))**2) )
      CoefW= MIN( max(CoefW,0.075),50.0)
	CTb=6.0e-3*CoefW

      FM(1)= FM(1) +CA*TFC*Aice(mg,i,j)*CTb
      BM(1)= BM(1) +CA*    Aice(mg,i,j)*CTb

*     5. Radiation through semi-clear ice    

      Hs= Hsnow(mg,i,j)/ai
c     Snow albedo. 
      asnow=F_as(Tsnow(mg,i,j))

c     If snow is thin - one can see ice
      albedoi=F_ai(Tice(mg,i,j),Hi)
	sf=hs/(hs+2.)

*     SW Dumping scale is 1.5m 
      Qrad= Qrad+ 
     +   aI0*Aice(mg,i,j)*(1.-sf)*(1.-albedoi)*SW(i,j)*EXP(-Hi/150.)


	end do   ! mg - ice thickness gradations

      FM(1)= FM(1)+ CA*Qrad
     #*(1.-(Ra*dzi1*(1.-exps1(2))+(1.-Ra)*dzi2*(1.-exps2(2)))/hz(1))
     #       /(row*cpw)

c     3.2. Net longwave radiation (Rosati & Miyakoda,1988)
c     Atmosphere vapor pressure epsa:
      TC= T(i,j,1)
      EPSs=0.98* (10.**( (0.7859+0.03477*TC)/(1.+0.00412*TC)+2.))
c      EPSs=10.**( (0.7859+0.03477*TC)/(1.+0.00412*TC)+2.)
      epsa=Q2m(i,j)*Pa(i,j)/0.622
	TK= TC+273.16
	TAK=Ta(i,j)+273.16

c     Rosati&Miyakoda, 1988.

      AA= Ew*sigma*(0.39-0.005*sqrt(epsa))*(1.-0.8*cloud(i,j))
	BB= 4.*Ew*sigma
c     Semimplicit time scheme 1.
      CC= BB*TAK
      RLWb= (4.*(AA+BB)*TK-3.*CC)*TK**2
      RLWf= (3.*(AA+BB)*TK-2.*CC)*TK**3 -273.16*RLWb
      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

c     Semimplicit time scheme 2.
c      RLWf= (-AA*TK +BB*Ta(i,j) )*TK**3
c      RLWb=  BB*TK**3     
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
c      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

c     Explicit time scheme
c      RLWf= (-AA*TK +BB*(Ta(i,j)-TC))*TK**3
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)

c     Marshunova, 1966 Explicit scheme.
c      RLWf= -Ew*sigma*(TK**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)

c     Marshunova, 1966 Semi-implicit scheme.
c      RLWf= Ew*sigma*(3.*TC-273.16)*TK**3 + 
c     &sigma*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))*TAK**4
c	RLWb= 4.*Ew*Sigma*TK**3
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
c      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

c     Koenig-Langlo & Augstein, 1994 (AWI) Semiimplicit scheme.
c      RLWf= Ew*sigma*(3.*TC-273.16)*(TK**3) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c	RLWb= 4.*Ew*Sigma*TK**3
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
c      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

*     4. Latent heat for open water.
c     Saturation vapor pressure EPSs:
	Qs=0.622*EPSs/Pa(i,j)
      Evapor= -roa*CDL(T(i,j,1),TA(i,j))*wind*(Qs-Q2m(i,j))  

      LH= QLw*Evapor
	Evapor= Evapor/row

      FM(1)= FM(1)+CA*Aice(0,i,j)*LH/(row*cpw)

*     ------------------- Deep water ---------------
      DO K=2, KB -1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

      AM(K)= -CA*CG( n)*AZT(i,j,k-1)/HZ(K-1)
      CM(K)= -CA*CG(np)*AZT(i,j,k  )/HZ(k)
      BM(K)= 1. -AM(K) -CM(K)
      FM(K)= T(i,j,k) 

*     Short wave.
      FM(K)=FM(K)+CA*Qrad*(
     *     Ra*dzi1*
     *( CG(n )*( EXPS1(k-1) -EXPS1(k  ))/hz(k-1)-
     *  CG(np)*( EXPS1(k  ) -EXPS1(k+1))/hz(k  )  )+	
     *(1.-Ra)*dzi2*
     *( CG(n )*( EXPS2(k-1) -EXPS2(k  ))/hz(k-1)-
     *  CG(np)*( EXPS2(k  ) -EXPS2(k+1))/hz(k  )  )
     *                                          )/(row*cpw)

      END DO
*     ----------------- Ocean bottom ----------------

      CA=2.*DT
      Q=1./HZ(KB-1)

      AM(KB)= -CA*Q*Q*AZT(i,j,kb-1)
      BM(KB)= 1. -AM(KB)
      CM(KB)= 0.
      FM(KB)= T(i,j,kb)

*     Short wave. Some flux to the bottom!
      FM(KB)=FM(KB)+CA*Q*Qrad*(
     *     Q*Ra*dzi1*(EXPS1(kb-1)-EXPS1(kb))+
     *(1.-Ra)*Q*dzi2*(EXPS2(kb-1)-EXPS2(kb))
     *     - Ra*EXPS1(kb) -(1.-Ra)*EXPS2(kb)
     *                                    )/(row*cpw)

*     3 Diagonal matrix factorization.
*     Vertical temperature distribution.
      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)
      Do k=1,kb
      T(i,j,k)= RKSI(K)
      End do


*     River + Precipitation - Evaporation.
      PME(i,j)= dt*( River(i,j) + Rain + Evapor*Aice(0,i,j))

*     Sea ice thermodynamics.
*     dHsnowT, dHiceT are total melted snow and ice mass,
*     used for vertical velocity calculation.
*     Ocean temperature is calculated already
*     (Except the case of complete ice melting and side melting - see IceTH).

      CALL ICETH(i,j,Wind,dHiceT,dHsnowT,SiceT)

*     Salinity changes if only snow melts dHsnowT<0.
      DHS= MIN(dHsnowT,0.)

*     River + Precipitation - Evaporation + Ice + Snow

      PME_DZ(i,j)= -roi*dHiceT/row ! Ice doesn't change the dyn. level
      PME(i,j)   = PME(i,j)-(roi*dHiceT+rosdry*DHS)/row

*     Vertical distribution of salinity.

*     -------------- Ocean surface -----------------
      CA=2.*DT/hz(1)
      AM(1)= 0.
      CM(1)= -CA*AZS(i,j,1)/hz(1)
      BM(1)= 1. -CM(1)
      FM(1)= S(i,j,1) -2.*Sice*dHiceT*roi/(row*hz(1)) ! salt from ice
ccc      FM(1)= S(i,j,1) -2.*SiceT*roi/(row*hz(1)) ! salt from ice

*     -------------- Deep water --------------------
      DO K=2, KB -1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

      AM(K)= -CA*CG( n)*AZS(i,j,k-1)/HZ(k-1)
      CM(K)= -CA*CG(np)*AZS(i,j,k  )/HZ(k)
      BM(K)= 1. -AM(K) -CM(K)
      FM(K)= S(i,j,k)

      END DO
*     ------------- Ocean bottom -----------------
      CA=2.*DT
      Q=1./HZ(KB-1)

      AM(KB)= -CA*Q*Q*AZS(i,j,kb-1)
      BM(KB)= 1. -AM(KB)
      CM(KB)= 0.
      FM(KB)= S(i,j,kb)

*     Matrix factorization
      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)

c	if(i.eq.13 .and. j.eq.4) then
c	write(*,*) '------'
c	write(*,*) (RKSI(k), k=1,kb)
c	end if


      Do k=1,kb

	if(RKSI(k) .LE. 0.) then
	write(*,*) 
     *i,j,k,RKSI(k),am(k),bm(k),cm(k),fm(k),S(i,j,k),SiceT,dHiceT
	end if

      S(i,j,k)= MAX(0.,RKSI(K))
      End do


c	IF(i.eq.8 .and. j.eq.13)  then
c	IF(i.eq.10 .and. j.eq.48) then
c	write(*,*) 'am:',am
c	write(*,*) 'bm:',bm
c	write(*,*) 'fm:',fm
c	write(*,*) 'cm:',cm
c	write(*,*) 'S:',rksi(1),rksi(2),rksi(3),rksi(4)
c	write(*,*) Ta(i,j),Rain, Evapor, dHiceT, DHS
c	end if



*     -----------  Check the stability  ----------------------------

c      do k=1,kb-1

c	kp=k+1

c	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
c      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
c      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)

c      if( Ro1 .LT. Ro2 ) then
c	write(*,*) 'Instability:',i,j,k, t(i,j,kp),t(i,j,k)

c      Az (i,j,k)= 1.e4
c      AzT(i,j,k)= 1.e3
c      AzS(i,j,k)= 1.e3

c	end if   ! isoneutral Ro is unstable

c      end do   ! k

*     --------------------------------- end of stability check
  
      END IF
      END DO
      END DO

*     Open Water "Compactness".
	do j=1,jl
      do i=1,il
      aice(0,i,j)= 1.0
      do m=1,mgrad
      aice(0,i,j)= MAX(0.,aice(0,i,j)-aice(m,i,j))
      end do 
	end do
	end do

      RETURN
      END

c----------------------------------------------------------------
      SUBROUTINE IceTH(i,j,Wind,dHiceT,dHsnowT,SiceT)
c----------------------------------------------------------------
c     Version #7.14 - 10.11.2012.
c     Potential Temperature and freezing point as a function of depth.
c     Variable conductivity coefficient ckice = ckice0 +beta*S/T
c     Method of Sections to Solve the Nonlinear Equation.
c     AOMIP Net LW radiation.
c     Flooding and aging with correspondent freshwater flux to ocean.
c     Semi-Transparent snow and ice.

**    Modification of side melt.

      include 'slo2.fi'
	REAL dHiceT, dHsnowT, SiceT
ccc	double precision f0,f1,f2,t0,t1,t2
      include 'tparm.fi'

	TAK=Ta(i,j)+273.16

c     Accuracy - 1 erg/cm2 = 0.001 W/m2 
      accur= 1.
      miter= 20

c     Atmosphere vapor pressure epsa:
      epsa=Q2m(i,j)*Pa(i,j)/0.622

c     Initial salt content od sea ice in the cell

	SI_init= 0.
      do m=1,mgrad
	if(aice(m,i,j).GT.aimin) then
	hi=Hice(m,i,j)/Aice(m,i,j)
	SI_init= SI_init +SaIce(hi,S(i,j,1))*Hice(m,i,j)
	end if
	end do

c     1. Old snow and ice.


      do m=1,mgrad

      dHsnow =0.
      dHice  =0.

      if( aice(m,i,j) .GT. Aimin) then

      Hi= Hice(m,i,j)/aice(m,i,j)
	ppp= 0. !!!1.e-5*g*row*Hi                   
      TFC= TFr(S(i,j,1), ppp) ! Freezing point at the depth of ice bottom

*     T,S dependent conductivity
c      c2= (ckice0 +2.*0.117e5*SaIce(hi,S(i,j,1))/(Tice(m,i,j)+TFC))/Hi

*     Constant conductivity 
      c2= ckice0/Hi 

      aI0= (0.18 +0.17*cloud(i,j))    ! Grenfell & Maykut, 1977. 
cc	aI0= 0.4*aI0                    ! fraction of penetration.     

c     SnowFall. Weatherly climate model.
      IF(TA(i,j).LT.-5.) then
	Hsnow(m,i,j)=Hsnow(m,i,j)+dt*row*Pr(i,j)*Aice(m,i,j)/rosdry
	end if
	IF(TA(i,j).GE.-5..AND.TA(i,j).LE.5.) then
	Hsnow(m,i,j)=Hsnow(m,i,j)+dt*row*Pr(i,j)*Aice(m,i,j)*
     &                               (1.-0.1*(TA(i,j)+5.))/rosdry
	end if

c     SnowFall. Simplified model.
c      IF(TA(i,j).LT. -1.) then
c	Hsnow(m,i,j)=Hsnow(m,i,j)+dt*row*Pr(i,j)*Aice(m,i,j)/rosdry
c	end if

c     1.1. Upper surface.

c     1.1.1. The case with snow cover.
c            Volumetric fusion heat.
      if( HSnow(m,i,j) .GT. Himin) then  
      Hs= Hsnow(m,i,j)/aice(m,i,j)
      c1= cksnow/Hs

c     Snow albedo. 
      asnow=F_as(Tsnow(m,i,j))

c     If snow is thin - one can see ice
      albedoi=F_ai(Tice(m,i,j),Hi)
	sf=hs/(hs+2.)

c     The case of linear temp. profile -> Semtner ice.
c     In a case of multilevel ice should be some
c     solution of the thermoconductivity equation.

c     Shortwave radiation.
      Qrad= (sf*(1.-asnow)+
     &       (1.-sf)*(1.-ai0)*(1.-albedoi)*(1.-EXP(-Hi/150.)))*SW(i,j)

c     Method of Sections.
c     T0 and T1 - initial estimates, f(T0)*f(T1)<0.
      T0=Ta(i,j)-20.

c     Saturation vapor pressure EPSs:
cc      EPSs=611.*10.**(9.5*T0/(273.16-7.66+T0))
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Es*sigma*((273.16+T0)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Es*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Es*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))
c     Total heat flux:
      F0=roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0)
     *   +c2*(TFC-T0)/(1.+c2/c1)
     *   +RLW + Qrad +LH

c     Initial estimate of T1.

      do L=1,5
	T1= 10.*REAL(1-L)
c     Saturation vapor pressure EPSs:
c      EPSs=611.*10.**(9.5*T1/(273.16-7.66+T1))
      EPSs=10.**( (0.7859+0.03477*T1)/(1+0.00412*T1)+0.00422*T1+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T1,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Es*sigma*((273.16+T1)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Es*sigma*(T1-Ta(i,j))*((273.16+T1)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Es*sigma*((273.16+T1)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T1)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F1=roa*Cpa*CDH(T1,TA(i,j))*wind*(TA(i,j)-T1)
     *   +c2*(TFC-T1)/(1.+c2/c1)
     *   +RLW + Qrad +LH
      if(f0*f1 .LT. 0.) goto 123
	end do 

123   iter=0
14    iter= iter +1

c     Regularization in a case of method stagnation.
	if( abs(F1-F0) .LT. 1.e-6) then
	F1=F1+accur
	T1=T1+0.05
	end if

	T2= T1 -F1*(T1-T0)/(F1-F0)

c     Saturation vapor pressure EPSs:
c      EPSs=611.*10.**(9.5*T2/(273.16-7.66+T2))
      EPSs=10.**( (0.7859+0.03477*T2)/(1+0.00412*T2)+0.00422*T2+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T2,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Es*sigma*((273.16+T2)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Es*sigma*(T2-Ta(i,j))*((273.16+T2)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Es*sigma*((273.16+T2)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T2)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F2=roa*Cpa*CDH(T2,TA(i,j))*wind*(TA(i,j)-T2)
     *   +c2*(TFC-T2)/(1.+c2/c1)
     *   +RLW + Qrad +LH 

	T0=T1
	T1=T2
	F0=F1
	F1=F2

      IF( ABS(F2).GT.accur.AND.iter.LT.miter) GOTO 14

      Tsnow(m,i,j)= T2

c     End of iterations to determine surface temperature.

      if (Tsnow(m,i,j) .LT. Tmelt_s) then
c     No snow melting.
      Tice(m,i,j)= (Tsnow(m,i,j) +c2*TFC/c1)/(1.+c2/c1)
      else

c     Snow melting.

      Tsnow(m,i,j)= Tmelt_s

c     New ice temperature.
      Tice(m,i,j)= (Tmelt_s +c2*TFC/c1)/(1.+c2/c1)

c     Wet snow albedo: Thin snow is transparent
      asnow=F_as(Tmelt_s)
      albedoi=F_ai(Tice(m,i,j),Hi)
	sf=hs/(hs+2.)
c     Shortwave radiation.
      Qrad= (sf*(1.-asnow)+(1.-sf)*(1.-ai0)*(1.-albedoi))*SW(i,j)

c     Heat flux in ice.
      Qsi= c1*(Tice(m,i,j)-Tmelt_s)
c     Total heat flux at snow surface.
c     Saturation vapor pressure EPSs:
      T0= Tmelt_s
c      EPSs=611.*10.**(9.5*T0/(273.16-7.66+T0))
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	Qs=0.622*EPSs/Pa(i,j)
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Es*sigma*((273.16+T0)**4) *(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Es*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Es*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

      Qas= roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0)+Qsi
     #     +RLW +LH+Qrad

c     Snow mass change.
      dHsnow= aice(m,i,j)*MIN(-dt*Qas/Qsnow, 0.)

c     If the mass of melted snow is large, the remaining heat
c     Qdelta melts underlying ice.
c     Volumetric specific heat of fusion.
c     Ice temperature is set to Tmelt_i.
         IF( Hsnow(m,i,j)+dHsnow .LT. 0.) then
           dHice= QSnow*(Hsnow(m,i,j)+dHsnow)/Qice
           dHsnow = -Hsnow(m,i,j)
           Tice(m,i,j)= Tmelt_i
         end if

      Hsnow(m,i,j)= Hsnow(m,i,j) + dHsnow
      Hice(m,i,j)= Hice(m,i,j) + dHice !!! Added 21.06.2008
      dHsnowT    = dHsnowT +dHsnow
      dHiceT     = dHiceT +dHice


	if(Hsnow(m,i,j).LT.0.1) Hsnow(m,i,j)=0. !Reg. for thin melting snow.
      end if ! T<Tmelt_s


c     Snow Sublimation
        if(Hsnow(m,i,j).GT.0.) then
c     Saturation vapor pressure EPSs:
          T0= Tsnow(m,i,j)
          EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	    Qs=0.622*EPSs/Pa(i,j)
	    Sublim= -roa*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))/rosdry
	    Hsnow(m,i,j)=Hsnow(m,i,j)+Aice(m,i,j)*Sublim*dt
	  end if

      else

c     1.1.2. No snow.
      
      albedoi=F_ai(Tice(m,i,j),Hi)
c     Shortwave radiation.
      Qrad= (1.-aI0)*(1.-albedoi)*SW(i,j)*(1.-EXP(-hi/150.)) 

c     Method of Sections.
c     T0 and T1 - initial estimates, f(T0)*f(T1)<0.
      T0=Ta(i,j)-20.

c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T0)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Ei*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F0= roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0) +c2*(TFC-T0)
     *   +RLW + Qrad + LH

c     Initial estimate of T1.
      do L=1,5
	T1= 10.*REAL(1-L)
c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T1)/(1+0.00412*T1)+0.00422*T1+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T1,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T1)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T1-Ta(i,j))*((273.16+T1)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Ei*sigma*((273.16+T1)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T1)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F1= roa*Cpa*CDH(T1,TA(i,j))*wind*(TA(i,j)-T1) +c2*(TFC-T1)
     *   +RLW + Qrad + LH
      if(f0*f1 .LT. 0.) goto 124
	end do 

c      if(f0*f1 .GT. 0.) then
c	write(*,*) 'Ice. Initial estimates are bad!!!'
c	write(*,*) f0,f1,t0,t1
c	end if

124   iter= 0
140   iter= iter+1

c     Regularization in a case of method stagnation.
	if( abs(F1-F0) .LT. 1.e-6) then
	F1=F1+accur
	T1=T1+0.05
	end if

	T2= T1 -F1*(T1-T0)/(F1-F0)

c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T2)/(1+0.00412*T2)+0.00422*T2+2.)
c     Specific humidity at the surface Qs:
	Qs=0.622*EPSs/Pa(i,j)
c     Latent heat:
	LH=-roa*QLi*CDL(T2,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T2)**4)*(0.39-0.005*sqrt(epsa))*
     &    (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T2-Ta(i,j))*((273.16+T2)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Ei*sigma*((273.16+T2)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T2)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

      F2= roa*Cpa*CDH(T2,TA(i,j))*wind*(TA(i,j)-T2) +c2*(TFC-T2)
     *   +RLW + Qrad + LH

	T0=T1
	T1=T2
	F0=F1
	F1=F2

c      if(kprint.eq.1) write(*,*)iter, T0,T1,T2,F2

      IF(ABS(F2).GT.accur.AND.iter.LT.miter) GOTO 140

c	if(iter.eq.miter.and.abs(f2).GT.accur)then
c	write(*,*) '============================='
c	write(*,*)'Ice temp.!', T1,T2,F2
c	write(*,*) '============================='
c	end if

      Tice(m,i,j)= T2

      if( Tice(m,i,j) .GE. Tmelt_i) then
      Tice(m,i,j)= Tmelt_i
      T0=Tmelt_i
c     New energy balance at ice surface.
c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	Qs=0.622*EPSs/Pa(i,j)
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T0)**4)*(0.39-0.005*sqrt(epsa))*
     &    (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo & Augstein, 1994 (AWI)
c      RLW= -Ei*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &                          (0.765+0.22*cloud(i,j)**3)
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

      Qrad=LH+(1.-aI0)*(1.-F_ai(Tmelt_i,Hi))*SW(i,j)*(1.-EXP(-hi/150.))
      Qai= roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0) + Qrad
     #     +RLW +c2*(TFC-T0)

c       Ice mass change on the ice top.
        dHice= aice(m,i,j)*MIN(-dt*Qai/Qice, 0.)
      end if

      if( (Hice(m,i,j)+dHice) .LT. 0.)then
        Qdelta= -Qice*(dHice+Hice(m,i,j))
        T(i,j,1)= T(i,j,1) + 2.*Qdelta/(hz(1)*row*Cpw)
        dHice= -Hice(m,i,j)
        Aice(m,i,j)=0.
      end if

      Hice(m,i,j)= Hice(m,i,j) + dHice
      dHiceT= dHiceT + dHice

c     Ice Sublimation
        IF(Hice(m,i,j) .GT. 0.) then
c     Saturation vapor pressure EPSs:
          T0= Tice(m,i,j)
          EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	    Qs=0.622*EPSs/Pa(i,j)
	    Sublim= -roa*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))/roi
	    Hice(m,i,j)=Hice(m,i,j)+Aice(m,i,j)*Sublim*dt
	  end if

      end if

c     1.2. Heat fluxes at the basement of the ice.
c          Ice temperature is not lower TFC.
*     McPhee, M.G. 1992. Turbulent heat fluxes ib the upper ocean
*     under ice. JGR, 97, 5365-5379.
*     5.5e-3 or hi dependent - Mellor&Kantha 1989
c      Hscale= 100.
c      CDw1= ( 0.4/alog(Hscale*300./(max(5.,0.5*hi))) )**2

      CDw1= CDw

      CoefW=SQRT(CDw1*((uice(i,j)-um2(i,j,1))**2
     &                +(vice(i,j)-vm2(i,j,1))**2))
      CoefW= MIN( max(CoefW,0.075),50.0)
	CTb=6.0e-3*CoefW

	ppp= 0. !!!1.e-5*g*row*Hi 
      TFC= TFr(S(i,j,1), ppp) ! Freezing point at the depth of ice bottom

      TW= MAX(TFC,T(i,j,1))     ! If T<TFC - ice formation
      Qiw= row*cpw*CTb*(TW-TFC)
      Qi = c2*(Tice(m,i,j)-TFC)

c     Ice bottom mass change. 
cc    dHice= -aice(m,i,j)*dt*(Qiw+Qi)/Qice

c     Heat of fusion is less factor 0.92 because
c     of brine pockets at the ice bottom - 
c     see Los Alamos Ice Model, Bitz & Limpscomb 1999.
      dHice= -1.087*aice(m,i,j)*dt*(Qiw+Qi)/Qice

c     If too much ice was melted, the remaining heat is going to
c     the very upper ocean layer. Snow gets to water and chills it.
c     WARNING: This is essentially explicit time scheme.

      if( (Hice(m,i,j)+dHice) .LT. 0.)then
        Qdelta= Qice*(-dHice-Hice(m,i,j)) -Qsnow*Hsnow(m,i,j)
        T(i,j,1)= T(i,j,1) + 2.*Qdelta/(hz(1)*row*Cpw)
        dHice= -Hice(m,i,j)
        Aice(m,i,j)=0.
	  dHsnow= -Hsnow(m,i,j) ! Complete Snow Melting
	  Hsnow(m,i,j)=0.
        dHsnowT    = dHsnowT +dHsnow
      end if

      Hice(m,i,j)= Hice(m,i,j) + dHice
      dHiceT= dHiceT + dHice

      end if
      end do  ! Loop over ice thickness gradations

	do m=1,mgrad

cc     Change due to side melting.     

cc      if(Hice(m,i,j).GT.Himin) then

cc      Hi= Hice(m,i,j)/max(Aimin,aice(m,i,j))

cc	ppp= 0.5e-5*g*row*Hi                  
cc      TFC= TFr(S(i,j,1), ppp) ! Freezing point at the half-depth of ice bottom
 
cc      TW= MAX(TFC,T(i,j,1))     ! If T<TFC - ice formation
cc	Rside= 3.14*1.6e-6*SQRT((TW-TFC)**3)*DT/(0.66*300.) ! See NCAR manual
cc      dAice=  -Rside*Aice(m,i,j)
cc      dHice=  -Rside*Hice(m,i,j)
cc      dHsnow= -Rside*Hsnow(m,i,j)
cc      Hsnow(m,i,j)= Hsnow(m,i,j)+dHsnow

cc      Hice(m,i,j)= Hice(m,i,j)+dHice
cc      Aice(m,i,j)= Aice(m,i,j)+dAice

cc      dHsnowT    = dHsnowT +dHsnow
cc      dHiceT= dHiceT +dHice

cc      Q_melt= Qsnow*dHsnow +Qice*dHice
cc      T(i,j,1)= T(i,j,1) +2.*Q_melt/(hz(1)*row*Cpw)

cc      else  
   
cc      Hice(m,i,j) = 0.
cc      Aice(m,i,j) = 0.
cc      Hsnow(m,i,j)= 0.
cc      end if ! Ice Side Melting.

      if(Hsnow(m,i,j).LT.Hsmin) Hsnow(m,i,j)=0.
      if(Hice(m,i,j).LT.Himin) then
	Hsnow(m,i,j)=0.
	Hice(m,i,j) =0.
	Aice(m,i,j) =0.
	end if

c     Snow to Ice conversion.

      if(Aice(m,i,j). GT. Aimin)then
	Hs=Hsnow(m,i,j)/Aice(m,i,j)
      Hi= Hice(m,i,j)/aice(m,i,j)

      if(Hs.GT.Hsmin) then
      
cc	if(Hs.GT.10.) then   ! Snow aging
	StoIce=GAMMA*dt*Hs
	Hs= Hs -StoIce
	Hi= Hi +StoIce*rosdry/roi
      dHsnow= -StoIce*aice(m,i,j)
	dHice = +StoIce*rosdry*aice(m,i,j)/roi
cc	end if

      Hsnow(m,i,j)= Hsnow(m,i,j)+dHsnow
      Hice(m,i,j) = Hice(m,i,j) +dHice
      dHsnowT    = dHsnowT +dHsnow
      dHiceT= dHiceT +dHice


c     Snow to ice conversion due to flooding.
      hdraft= (rosdry*Hs+roi*Hi)/row
	hflood= MAX(0.,hdraft-Hi)
	Hs= Hs -hflood*roi/rosdry
	Hi= Hi +hflood

	Hsnow(m,i,j)= Hs*aice(m,i,j)
	Hice (m,i,j)= Hi*aice(m,i,j)

      dHsnow= -hflood*aice(m,i,j)*roi/rosdry
	dHice = +hflood*aice(m,i,j)

      Hsnow(m,i,j)= Hsnow(m,i,j)+dHsnow
      Hice(m,i,j) = Hice(m,i,j) +dHice
      dHsnowT     = dHsnowT +dHsnow
      dHiceT      = dHiceT +dHice

	end if       
	end if

      end do  ! Loop over ice thickness gradations

c  2. New ice formation.

      dHiceN= 0.

	kb=km2(i,j)
      do k=1,kb

      if (k.lt.kb) then
         if(k.EQ.1)then
         vol=0.5*hz(1)*cg(abs(nt3(i,j,1)))
         else
         vol=0.5*(hz(k-1) *cg(abs(nt3(i,j,k  )))
     *           +hz(k)   *cg(abs(nt3(i,j,k+1))) )
         end if
      else
         vol=0.5*hz(kb-1)*cg(abs(nt3(i,j,kb)))
      end if

	ppp= 1.e-5* g*row*z(k) 
      TFC= TFr(S(i,j,k),ppp) ! potential temperature

      if(T(i,j,k).LT.TFC ) then 
cc        dHiceN=dHiceN + Cpw*row*vol*(TFC-T(i,j,k))/
cc     #                  Qice/cg(abs(nt3(i,j,1))) 
c     Heat of fusion is less factor 0.92 because
c     of brine pockets at the ice bottom - 
c     see Los Alamos Ice Model, Bitz & Limpscomb 1999.
      dHiceN=dHiceN + 
     *1.087*Cpw*row*vol*(TFC-T(i,j,k))/Qice/cg(abs(nt3(i,j,1)))
      T(i,j,k)= TFC                 
      end if  ! T<TFC
      end do  ! k=1,kb loop

	IF(dHiceN .GT. 0.) then

c     New ice temperature - for T in situ only!
c      if(hice(1,i,j).GT.himin) then
c      Q= (Tice(1,i,j)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
c      Tice(1,i,j)=Q/Hice(1,i,j)
c      else
c      Tice(1,i,j)=TFC
c      end if

c     In a case of sufficient open water: frazil ice to the open water.
	hhh0= dHiceN/MAX(aimin,aice(0,i,j))

c     Thickness of the new ice

      IF( hhh0 .LE. Hmax(2)) then
	delta=dHiceN/Href

      Hice(1,i,j)=Hice(1,i,j)+dHiceN
      aice(1,i,j)= aice(1,i,j)+ MIN(delta, aice(0,i,j))
      aice(0,i,j)= MAX(0.,aice(0,i,j)- delta)

      else

c     New Ice mass spreaded unifirmly over entire cell
c     if no sufficient open water
      Hice(1,i,j)=Hice(1,i,j)+dHiceN*(Aice(0,i,j)+Aice(1,i,j))
      do m=2,mgrad
	hi= Hice(m,i,j)/max(aimin,Aice(m,i,j))
	Hice(m,i,j)=Hice(m,i,j)+dHiceN*Aice(m,i,j)
	end do

      aice(1,i,j)= aice(1,i,j)+ aice(0,i,j)
      aice(0,i,j)= 0.

	end if

      dHiceT= dHiceT + dHiceN

	end if ! New ice > Himin

C     Change of salt content in ice

	SI_final= 0.
      do m=1,mgrad
	if(aice(m,i,j).GT.aimin) then
	hi=Hice(m,i,j)/Aice(m,i,j)
	SI_final= SI_final +SaIce(hi,S(i,j,1))*Hice(m,i,j)
	end if
	end do

	SiceT= SI_final -SI_init

	if(SiceT .LT. -1.e4) then
      write(*,*) 'in TSice', i,j, SI_init, SI_final
      end if


      return
      end

C---------------------------------------------------------------------
      SUBROUTINE FACTOR(ND,A,B,C,F,RKSI,II,JJ)
C---------------------------------------------------------------------
C UP-DOWN : A(N)*Fi(N-1) +B(N)*Fi(N) +C(N)*Fi(N+1)=F(N) ; N=II,JJ
C A(II) AND C(JJ) ARE NOT USED (I HOPE SO), SO THEY ARE SET = 0.
      Parameter (MDIM=100)
      REAL A(ND),B(ND),C(ND),F(ND),RKSI(ND)
      REAL X(MDIM),Y(MDIM)

      A(II) = 0.
      C(JJ) = 0.
      X(II) = - C(II) /B(II)
      Y(II) =   F(II) /B(II)

      JJ1 = II+1
      DO J1=JJ1,JJ
      denom = A(J1) * X(J1-1) + B(J1)
ccc      IF(ABS(Denom) .LT. 1.e-29) denom= C(j1)
      X(J1) = -C(J1) / Denom
      Y(J1) = (F(J1) - A(J1)*Y(J1-1)) / Denom
      END DO

      RKSI(JJ)=Y(JJ)
      
	JJ1 = JJ - 1
      
	DO J1=JJ1,II,-1
      RKSI(J1) = X(J1) * RKSI(J1+1) + Y(J1)      
	END DO

      RETURN
      END

      function TFr(S,p)

*     Version 06.06.2012.
*     For global models it's recommended to use simplified formula 
*     to estimate the Tfr

	real S,p
	real*8 s8,p8,pn,pd,sqs

*     Algorithms for in situ temperature/
*
*     Makshtas A.P. The heat budget of Arctic ice in the winter.
*     Publ. by Int. Glaciological Soc. Cambridge CB2 1ER. UK. -
*     1991. 77 p.
c      TF= -0.054*S

*     Mellor&Kantha 1989 -0.0543    

c     TF= -0.0545*S

*     Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
*     AOMIP - they also refers to Millero.  But Linear.
*     Millero, F. J., 1978: Annex 6, freezing point of seawater. 
*     Unesco Tech. Papers in the Marine Sciences, 28, 29-35.

c      TFr=-0.0575*S+1.710523e-3*SQRT(S**3)-2.155e-4*S*S-7.53e-9*p

*     ------------------------------------------------------------------
*     Algorithm for potential temperature
!     Jackett, D. R., McDougall, T. J., Feistel, R., Wright, 
!     D. G., and Griffies, S. M.: 
!     Algorithms for density, potential temperature, conservative
!     temperature, and freezing temperature of seawater, 
!     Journal of Atmospheric and Oceanic Technology, 23, 1709–1728, 2006.
!
!    s                : salinity                           (psu)
!    p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)

      s8=DBLE(s)
	p8=DBLE(p)
      sqs=DSQRT(S8)


      Pn =2.5180516744541290d-03                  +
     *    s8*(-5.8545863698926184d-02             +
     *        sqs*( 2.2979985780124325d-03        -     
     *              sqs*3.0086338218235500d-04))  +     
     *        p8*(-7.0023530029351803d-04         +     
     *            p8*( 8.4149607219833806d-09     +     
     *            s8*1.1845857563107403d-11))

      Pd =1.0d0                            +
     *    p8*(-3.8493266309172074d-05      +
     *        p8*9.1686537446749641d-10)   +
     *    s8*s8*sqs*1.3632481944285909d-06 


      TFr= REAL(Pn/Pd) ! Pure water

	TFr=TFr - 2.518052e-3 + 1.428571e-5*S ! Air saturated water

      return
      end

      Function F_ai(T,h)
*     Mellor and others
*	parameter (aow=0.10, ami=0.55, c11=0.44, c12=0.075)

	parameter (aow=0.10) ! AOMIP
cc	parameter (aow=0.06) ! NCAR, LANL

	if(T .LT. -1.0)then  ! -1C - see BoulderIce
cc  	F_ai= MAX(aow, c11*SQRT(SQRT(0.01*h))+aow)
cc      ai= MIN(1.,F_ai)
cc      ai= 0.62
cc      ai= 0.60  ! AOMIP, approx. 3m.
cc      ai= 0.70  ! New AOMIP, 15.05.2003 approx. 4m.
cc      ai= 0.65  ! Boulder Ice CCSM3 (Ice version 4, 2002) 
      ai= 0.55  ! Boulder Ice CCSM3 (Ice version 5, 2004) 

	if(h.LE.50.)then  ! see BoulderIce - this is a good approx.
	F_ai= aow +(ai-aow)*h/50.
	else
	F_ai= ai
	end if

	else

cc      ai= 0.65  ! Boulder Ice CCSM3 (Ice version 4, 2002)
cc      ai= 0.60  ! AOMIP, approx. 3m.
      ai= 0.55  ! Boulder Ice CCSM3 (Ice version 5, 2004) 

	if(h.LE.50.)then  ! see BoulderIce - this is a good approx.
	F_ai= aow +(ai-0.075*(T+1.)-aow)*h/50. 
	else
	F_ai= ai -0.075*(T+1.)
	end if

cc	F_ai= MIN(ami,c12*1.e-4*h**2 +aow)
cc      ai= 0.50  ! Weatherly and AOMIP, SHEBA, ~2.5m
cc      ai= 0.52  ! Boulder Ice
cc      ai= 0.55  ! Pinto, 1999, SHEBA  
cc      ai= 0.60
cc      ai= 0.68  ! New AOMIP, 15.05.2003
	end if

	return
	end
	
	Function F_as(T)
*     AOMIP lbedoes 21.05.03. or SHEBA=Weatherly
	if(T .LT. -1.0)then   ! -1C - see BoulderIce
cc	F_as=0.75  ! Flato&Brown 1986
	F_as=0.80  ! Old Aomip, close to CCSM3
cc	F_as=0.82  ! SHEBA
cc	F_as=0.85  ! Pinto, 1999, SHEBA, Boulder Ice
cc	F_as=0.81  ! New AOMIP, 15.05.2003
	else
cc	F_as=0.65  ! F&B 1986
cc	F_as=0.70  ! Old AOMIP SHEBA 2001 -SHCE
cc	F_as=0.72  ! Boulder
	F_as=0.80 -0.1*(T+1.0)  ! Pinto, 1999, SHEBA + Boulder Ice, CCSM3
cc	F_as=0.77  ! New AOMIP, 15.05.2003
	end if
	return
	end 

	Function CDH(T,TA)
*     version 17.12.2009. Regularisation for sections method.
	if(T.GT.TA+1.0) then
	CDH=1.75e-3       ! Unstable boundary layer, Makshtas, SP-23.
	else
	  if(T.LT.TA-1.0) then
        CDH=1.2e-3      ! Stable boundary layer, AOMIP
        else
        CDH=1.2e-3 +0.5*(T-TA+1.0)*0.55e-3 ! Transition case
	  end if
	end if
      return
	end
	
	Function CDL(T,TA)
*     version 17.12.2009. Regularisation for sections method.
	if(T.GT.TA+1.0) then
	CDL=1.75e-3       ! Unstable boundary layer, Makshtas, SP-23.
	else      
	  if(T.LT.TA-1.0) then
        CDL=1.5e-3      ! Stable boundary layer, AOMIP
        else
        CDL=1.5e-3 +0.5*(T-TA+1.0)*0.25e-3 ! Transition case
	  end if
      end if
      return
	end

	Function SaIce(hi,S)
*     Version 19.10.2010
*     sea ice thickness at hi cm.
*     S - ocean salinity

*     Compilated according data by A. Kovacs, 1996. Sea Ice. Part I.
*     Balk salinity versus floe thickness. US Army Corps of Engeneers,
*     CRREL Rep. 96-7., 16 p.

*     Data for multiyear ice were correted to make the function continious.

      hi_in= max(hi,5.)

      if(hi_in .LE. 200.) then
      a= 4.606+91.603/hi_in
	else
	a= 1.85+128560.6/(hi_in**2)
	end if

	a=4.   ! In a case of constant ice salinity

      SaIce= MIN(a, max(0.,S)) ! ice salinity is <= the ocean one

      return
	end