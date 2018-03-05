      SUBROUTINE TSICE

c     Version 07.07.2015 for ice thickness gradations.
c     Potential temperature and Freezing point as function on depth.

c     Nonuniform heat flux ice-ocean
c     Spatially Varied Cloudiness
c     AOMIP forcing - Net LW radiation.
c     Penetration of the Short Wave Radiation into the Ocean.
c     Friction velocity as in Los Alamos model.
c     Evaporation and snow/ice sublimation. 

**    Temperature flux due to Rain and Rivers
c     Heat flux by snow fusion in precipitation.

      INCLUDE 'Slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
	REAL dHiceT, dHsnowT
      INCLUDE 'Tparm.fi'

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

      KB = KM2(I,J)
	n=ABS(nt3(i,j,1))
*-------------------------------------------------------
*     Calculation in all points.

      IF( KB .GT. 0) THEN
      
*     -------------- Ocean Surface -----------------

*     Precipitation: Rain (or snow) to the open water.
*     A.Dai. GRL, v.35, 2008.
	IF(TA(i,j).GT. 5.) then
	Rain = Pr(i,j)
	Psnow= 0.
	end if
	IF(TA(i,j).LT.-1.) then
	Psnow= Pr(i,j)*Aice(0,i,j)
	Rain = 0.
	end if
	IF(TA(i,j).GE.-1..AND.TA(i,j).LE.5.) then
	Rain= Pr(i,j)*(Aice(0,i,j)+(1.-Aice(0,i,j))*
     &                           0.166666*(TA(i,j)+1.) )
      Psnow= Pr(i,j)*Aice(0,i,j)*(1.- 0.166666*(TA(i,j)+1.) )
	end if	
	
*     Wind is in m/s.
      WIND=100.*(SQRT(wx(i,j)**2+wy(i,j)**2) +E0)
	
      CA=2.*DT/hz(1)

*     Rain and snow fluxes. Temperature of snow < 0C
      AM(1)= 0.
      CM(1)= -CA*AZT(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) +CA*(Rain+Psnow)
      FM(1)= T(i,j,1)+CA*(Rain*TA(i,j)+Psnow*min(0.,TA(i,j)))

*     Snow melting on the open water    
      FM(1)= FM(1) -CA*Psnow*rosdry*Qsnow/(row*cpw)

*     Heat fluxes through ocean and ice.

*     1. Sensible heat at open water.
      BM(1)=BM(1)+CA*roa*Cpa*CDH(T(i,j,1),TA(i,j))*
     &            WIND        *Aice(0,i,j)/(row*cpw)
      FM(1)=FM(1)+CA*roa*Cpa*CDH(T(i,j,1),TA(i,j))*
     &            WIND*TA(i,j)*Aice(0,i,j)/(row*cpw)
     
     
*     Antidiffusive fluxes from previous time step
      BM(1)=BM(1) +CA*PME(i,j)/dt     


*     3. Incoming radiation.
*     3.1. Short wave at open water.

	Qrad= Aice(0,i,j)*(1.-aopw)*SW(i,j)

      aI0= 0.18 +0.17*cloud(i,j)    ! Grenfell & Maykut, 1977. 
cc	aI0= 0.4*aI0                  ! To account volume heating for 0-layer


      do mg=1,mgrad
	Ai=MAX(Aimin,Aice(mg,i,j))
	Hi=Hice(mg,i,j)/Ai

	ppp= 1.e-5*g*row*Hi
      TFC= TFr(S(i,j,1), ppp) ! Freez T at the depth of ice bottom


*     2. Sensible heat at ice-water interface. 

* --------------------------------------------------------------------------
*     McPhee, M.G. 1992. Turbulent heat fluxes in the upper ocean
*     under ice. JGR, 97, 5365-5379.
*     Q=rw*Cp*CDw*CoefW*dT, CoefW=SQRT(CDw*((uice-u)**2) -friction velocity.
*     Skin velocity W=1-3 cm/s, Mellor&Kantha 0.5.
*     Los Alamos sea ice model - cut off 0.075 -> dU=1cm/s
*     Drag coefficient
*     Cdw= const or hi dependent - Mellor&Kantha 1989
c      Hscale= 100.
c      CDw1=( 0.4/alog(Hscale*300./(max(5.,0.5*hi))) )**2

	u2= (uice(i,j)-um2(i,j,1))**2
	v2= (vice(i,j)-vm2(i,j,1))**2

*     Friction velocity
	u_star= SQRT(CDw*(u2 + v2))
	u_star=MAX(5.e-2, u_star)

****	u_star=min(100.,SQRT(0.21*Q2turb(i,j,1)))

	CTb=6.0e-3*u_star           ! CICE
	
**	write(*,*) 'Ctb:', CTb, u_star
**	CTb=7.27e-3                 ! Ebert & Curry, 1993

* ----------------------------------------------------------------
* Omstedt, A., and J. S. Wettlaufer, 
* Ice growth and oceanic heat-flux: Models and measurements,
* J. Geophys. Res., 97, 9383–9390, 1992.

c      CTb= 2.e-4*
c     * SQRT((uice(i,j)-um2(i,j,1))**2 +(vice(i,j)-vm2(i,j,1))**2)
* ----------------------------------------------------------------    


      FM(1)= FM(1) +CA*TFC*Aice(mg,i,j)*CTb
      BM(1)= BM(1) +CA*    Aice(mg,i,j)*CTb

*     5. Radiation through semi-clear ice. We need snow\ice temperatures    

*      Hs= Hsnow(mg,i,j)/ai
c     Snow albedo. 
*      asnow=F_as(Tsnow(mg,i,j))

c     If snow is thin - one can see ice
*      albedoi=F_ai(Tice(mg,i,j),Hi)
*	sf=hs/(hs+2.)

*     SW Dumping scale is 1.5m 
c      Qrad= Qrad+ 
c     +   aI0*Aice(mg,i,j)*(1.-sf)*(1.-albedoi)*SW(i,j)*EXP(-Hi/150.)


	end do   ! mg - ice thickness gradations

      FM(1)= FM(1)+ CA*Qrad
     #*(1.-(Ra*dzi1*(1.-exps1(2))+(1.-Ra)*dzi2*(1.-exps2(2)))/hz(1))
     #       /(row*cpw)

c     3.2. Net longwave radiation
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

c     Semiimplicit time scheme 3.
      RLWb= (4.*AA*TK-3.*BB*(TAK-TK))*TK**2
      RLWf= (3.*AA*TK-2.*BB*(TAK-TK))*TK**3-273.16*RLWb
      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

c     Semiimplicit time scheme 1.
c      CC= BB*TAK
c      RLWb= (4.*(AA+BB)*TK-3.*CC)*TK**2
c      RLWf= (3.*(AA+BB)*TK-2.*CC)*TK**3 -273.16*RLWb
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
c      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

c     Explicit time scheme
c      RLWf= (-AA*TK +BB*(TaK-TK))*TK**3
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

c     Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973, Marsunova, 1966
c     Semiimplicit
c	RLWb= 4.*Ew*Sigma*TK**3
c     RLWf= 3.*Ew*sigma*TK**4 + sigma*(TAK**4)*
c     &     (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     &     -273.16*RLWb
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
c      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

*     4. Latent heat and evaporation for open water.
c     Saturation vapor pressure EPSs:

	Qs=0.622*EPSs/Pa(i,j)
      Evapor= -roa*CDL(T(i,j,1),TA(i,j))*wind*(Qs-Q2m(i,j))  

      LH= QLw*Evapor    ! Latent Heat
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
      
*     Short wave. Some flux is going to the bottom!
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
      PME(i,j)= dt*(River(i,j) +Rain +Psnow +Evapor*Aice(0,i,j))

*     Sea ice thermodynamics.
*     dHsnowT, dHiceT are total melted snow and ice mass,
*     used for vertical velocity calculation.
*     Ocean temperature is calculated already
*     (Except the case of complete ice melting and 
*     side melting - see IceTH).

      CALL IceTH(i,j,Wind,dHiceT,dHsnowT)

*     Salinity and level changes if only snow melts dHsnowT<0.
      DHS= MIN(dHsnowT,0.)

*     River + Precipitation - Evaporation + Ice + Snow

c      PME_dz(i,j)= -roi*dHiceT/row ! ice trans doesn't change the dyn.level.
      PME(i,j)= PME(i,j)-(roi*dHiceT+rosdry*DHS)/row

*     Vertical distribution of salinity.

*     -------------- Ocean surface -----------------
      CA=2.*DT/hz(1)
      AM(1)= 0.
      CM(1)= -CA*AZS(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) -River(i,j)*dt/hz(1) ! River s=0.5*Ocean?
      FM(1)= S(i,j,1)  -2.*Sice*dHiceT*roi/(row*hz(1))
      
*     Antidiffusive fluxes
      BM(1)=BM(1) +CA*PME(i,j)/dt     


*     -------------- Deep water --------------------
      DO K=2, KB-1
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

*     ------------- Matrix factorization ---------

      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)

      Do K=1,KB
      S(i,j,k)= MAX(0.,RKSI(K))
      End do


      END IF
      END DO
      END DO


      RETURN
      END

c----------------------------------------------------------------
      SUBROUTINE IceTH(i,j,Wind,dHiceT,dHsnowT)
c----------------------------------------------------------------
c     Version #10.4 - 07.07.2015.
c     Potential Temperature and freezing point as a function of depth.
c     Variable conductivity coefficient ckice = ckice0 +beta*S/T
c     Method by Ridders to Solve the Nonlinear Equation.
c     C.J.F. Ridders, IEEE Transactions on Circuits and Systems, 
c     v. CAS-26, p.979, 1979.
c     AOMIP Net LW radiation.
c     Flooding and aging with correspondent freshwater flux to ocean.
c     Semi-Transparent snow and ice.

**    Modification of side melt.

      include 'Slo2.fi'
	REAL dHiceT, dHsnowT
ccc	double precision f0,f1,f2,t0,t1,t2
      include 'Tparm.fi'


      dHiceT = 0.
      dHsnowT= 0.

	TAK=Ta(i,j)+273.16

c     Accuracy - 1000 erg/cm2 = 1.0 W/m2 
      accur= 1.e2
      miter= 100

c     Atmosphere vapor pressure epsa:
      epsa=Q2m(i,j)*Pa(i,j)/0.622

c     1. Old snow and ice.

      do m=1,mgrad

      dHsnow =0.
      dHice  =0.

      if( aice(m,i,j) .GT. Aimin) then

      Hi= Hice(m,i,j)/aice(m,i,j)
	ppp= 1.e-5*g*row*Hi
      TFC= TFr(S(i,j,1), ppp) ! Freezing T at the depth of ice bottom

*     T,S dependent conductivity
c      c2= (ckice0 +2.*0.117e5*SaIce(hi,S(i,j,1))/(Tice(m,i,j)+TFC))/Hi

*     Constant conductivity 
      c2= ckice0/max(0.1,Hi)

      aI0= 0.18 +0.17*cloud(i,j)    ! Grenfell & Maykut, 1977. 
cc	aI0= 0.4*aI0                  ! To account volume heating for 0-layer


c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, before snowfall', m, Hsnow(m,i,j), TA(i,j)
c	end if


c     SnowFall. A.Dai. GRL, v.35, 2008.
      
	deltaHS=0.

      IF(TA(i,j).LT.-1.) then
	deltaHS = dt*row*Pr(i,j)*Aice(m,i,j)/rosdry
	end if

	IF(TA(i,j).GE.-1..AND.TA(i,j).LE.5.) then
	deltaHS = dt*row*Pr(i,j)*Aice(m,i,j)*
     &                     (1.-0.16666*(TA(i,j)+1.))/rosdry
	end if


c      if(Hsnow(m,i,j) .LT. hsmin .and. deltaHS .gt. hsmin) then
*     Temperature of thin new snow
c      Tsnow(m,i,j) = TA(i,j)-5.
c      end if

	Hsnow(m,i,j) = Hsnow(m,i,j)+deltaHS


c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, after snowfall', m, Hsnow(m,i,j), Tsnow(m,i,j)
c	end if



c     1.1. Upper surface.

c     1.1.1. The case with snow cover.
c            Volumetric fusion heat.
      if( HSnow(m,i,j) .GT. Hsmin) then 
	 
      Hs= Hsnow(m,i,j)/aice(m,i,j)
      c1= cksnow/max(0.1,Hs)

	sf= hs/(hs+2.) !! May be problems with too fast melting


c     The case of linear temp. profile -> Semtner ice.
c     In a case of multilevel ice should be some
c     solution of the thermoconductivity equation.


c     Method by Ridders.
c     T0 and T1 - initial estimates, f(T0)*f(T1)<0.

      T0= -100. !!!TA(i,j)-30.

c     Snow albedo. 
      asnow=F_as(T0)
c     If snow is thin - one can see ice
c     New ice temperature.
      Ti= (c1*T0 +c2*TFC)/(c1+c2)
      albedoi=F_ai(Ti,Hi)

c     Shortwave radiation.
*      Qrad= (sf*(1.-asnow)+
*     &      (1.-sf)*(1.-ai0*EXP(-Hi/150.))*(1.-albedoi))*SW(i,j)
      Qrad= ( sf*(1.-asnow)+
     &       (1.-sf)*(1.-ai0)*(1.-albedoi) )*SW(i,j)
**      Qrad= (1.-asnow)*SW(i,j)


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
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Es*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))
c     Total heat flux:
      F0=roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0)
     *   +c2*(TFC-T0)/(1.+c2/c1)
     *   +RLW + Qrad +LH

c     Initial estimate of T1.

***      L=0
***	F1=F0

***      do while (F1*F0 .GT. 0. .and. L.LT.10)
***	L= L+1
***	T1= T0 +20.*REAL(L)

      T1= 20.

c     Snow albedo. 
      asnow=F_as(T1)
c     If snow is thin - one can see ice
c     New ice temperature.
      Ti= (c1*T1 +c2*TFC)/(c1+c2)
      albedoi=F_ai(Ti,Hi)

c     Shortwave radiation.
c      Qrad= (sf*(1.-asnow)+
c     &      (1.-sf)*(1.-ai0*EXP(-Hi/150.))*(1.-albedoi))*SW(i,j)
      Qrad= ( sf*(1.-asnow)+
     &       (1.-sf)*(1.-ai0)*(1.-albedoi) )*SW(i,j)
**      Qrad= (1.-asnow)*SW(i,j)


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
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Es*sigma*((273.16+T1)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T1)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F1=roa*Cpa*CDH(T1,TA(i,j))*wind*(TA(i,j)-T1)
     *   +c2*(TFC-T1)/(1.+c2/c1)
     *   +RLW + Qrad +LH


***	end do 

      if(f0*f1 .GT. 0.) then
cc	write(*,*) 'Snow. Initial estimates are bad!!!'
cc	write(*,*) i,j,m,t0,t1,hsnow(m,i,j),tsnow(m,i,j)

	T2= 20. !  No solution in physically reasonable interval

	else

c     Iterations to determine surface temperature

      iter=0
	F2=F0

	do while (ABS(F2).GT.accur.AND.iter.LT.miter)

      iter= iter +1

	T05= 0.5*(T0+T1)

c     Snow albedo. 
      asnow=F_as(T05)
c     If snow is thin - one can see ice
c     New ice temperature.
      Ti= (c1*T05 +c2*TFC)/(c1+c2)
      albedoi=F_ai(Ti,Hi)

c     Shortwave radiation.
cc      Qrad= (sf*(1.-asnow)+
cc     &      (1.-sf)*(1.-ai0*EXP(-Hi/150.))*(1.-albedoi))*SW(i,j)
      Qrad= ( sf*(1.-asnow)+
     &       (1.-sf)*(1.-ai0)*(1.-albedoi) )*SW(i,j)
**      Qrad= (1.-asnow)*SW(i,j)
     

c     Saturation vapor pressure EPSs:
c      EPSs=611.*10.**(9.5*T1/(273.16-7.66+T1))
      EPSs=10.**( (0.7859+0.03477*T05)/(1+0.00412*T05)+0.00422*T05+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T05,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Es*sigma*((273.16+T05)**4)*(0.39-0.005*sqrt(epsa))*
     &(1.-0.8*cloud(i,j))-4.*Es*sigma*(T05-Ta(i,j))*((273.16+T05)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Es*sigma*((273.16+T1)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T1)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F05=roa*Cpa*CDH(T05,TA(i,j))*wind*(TA(i,j)-T05)
     *   +c2*(TFC-T05)/(1.+c2/c1)
     *   +RLW + Qrad +LH


	T2=T05 
     *  +(T05-T0)*SIGN(1.,F0-F1)*F05/SQRT(max(1.e-12,F05**2-F1*F0))

c     Snow albedo. 
      asnow=F_as(T2)
c     If snow is thin - one can see ice
c     New ice temperature.
      Ti= (c1*T2 +c2*TFC)/(c1+c2)
      albedoi=F_ai(Ti,Hi)

c     Shortwave radiation.
c      Qrad= (sf*(1.-asnow)+
c     &      (1.-sf)*(1.-ai0*EXP(-Hi/150.))*(1.-albedoi))*SW(i,j)
      Qrad= ( sf*(1.-asnow)+
     &       (1.-sf)*(1.-ai0)*(1.-albedoi) )*SW(i,j)
**      Qrad= (1.-asnow)*SW(i,j)
     


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
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Es*sigma*((273.16+T2)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T2)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F2=roa*Cpa*CDH(T2,TA(i,j))*wind*(TA(i,j)-T2)
     *   +c2*(TFC-T2)/(1.+c2/c1)
     *   +RLW + Qrad +LH 

	if(iter.eq.miter.and.abs(f2).GT.accur)then
	write(*,*) '============================='
	write(*,*)'Snow temp.!', T0,T1,T05,T2,F0,F1,F05,F2,i,j,m
	write(*,*) '============================='
	end if


      IF(F05*F2 .LT. 0.) then

	T0= T05
	T1= T2
	F0= F05
	F1= F2

	else

	IF(F0*F2 .LT. 0.) then
	T1=T2
	F1=F2
	end if

	IF(F1*F2 .LT. 0.) then
	T0=T2
	F0=F2
	end if

	end if

	end do

	end if

      Tsnow(m,i,j)= T2

!==========     End of iterations to determine surface temperature ========


      if (T2 .LT. Tmelt_s) then   ! No snow melting.

      Tice(m,i,j)= (c1*T2 +c2*TFC)/(c1+c2)

      else                                  ! Snow melting.

      T0= Tmelt_s
      Tsnow(m,i,j)= Tmelt_s

c     New ice temperature.
      Ti= (c1*Tmelt_s +c2*TFC)/(c1+c2)
	Tice(m,i,j) = Ti

c     Wet snow albedo: Thin snow is transparent
      asnow=F_as(Tmelt_s)
      albedoi=F_ai(Ti,Hi)

c     Shortwave radiation.
      Qrad= (sf*(1.-asnow)+(1.-sf)*(1.-ai0)*(1.-albedoi))*SW(i,j)
**      Qrad= (1.-asnow)*SW(i,j)
      

*      Qrad= (sf*(1.-asnow)+
*     &      (1.-sf)*(1.-ai0*EXP(-Hi/150.))*(1.-albedoi))*SW(i,j)


c     Heat flux in ice.
*      Qsi= c1*(Tice(m,i,j)-Tmelt_s)

	Qsi= c1*c2*(TFC-Tmelt_s)/(c1+c2)   ! another form - no arrays

c     Total heat flux at snow surface.
c     Saturation vapor pressure EPSs:
c      EPSs=611.*10.**(9.5*T0/(273.16-7.66+T0))
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	Qs=0.622*EPSs/Pa(i,j)
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Es*sigma*((273.16+T0)**4) *(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Es*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Es*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Es*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

      Qas= roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0)+Qsi
     #     +RLW +LH+Qrad

c     Snow mass change.
      dHsnow= aice(m,i,j)*MIN(-dt*Qas/Qsnow, 0.)
	dHice=0.

c     If the mass of melted snow is large, the remaining heat
c     Qdelta melts underlying ice.
c     Volumetric specific heat of fusion.
c     Ice temperature is set to Tmelt_i.
         IF( Hsnow(m,i,j)+dHsnow .LT. 0.) then
           dHice  = QSnow*(Hsnow(m,i,j)+dHsnow)/Qice
           dHsnow = -Hsnow(m,i,j)
         end if

      Hsnow(m,i,j)= Hsnow(m,i,j) + dHsnow
      Hice (m,i,j)= Hice(m,i,j) + dHice
      dHsnowT     = dHsnowT +dHsnow
      dHiceT      = dHiceT  +dHice



c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow melting', Hsnow(m,i,j), dHsnow, dHice, 
c     &   Qrad,Qsi,RLW,LH,roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0)
c	end if


      end if ! T<Tmelt_s


c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, after', m, Hsnow(m,i,j)
c	end if


c     Snow Sublimation - no flux to ocean!
        if(Hsnow(m,i,j).GT.0.) then
c     Saturation vapor pressure EPSs:
          T0= Tsnow(m,i,j)
          EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	    Qs=0.622*EPSs/Pa(i,j)
	    Sublim= -roa*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))/rosdry
	    Hsnow(m,i,j)=Hsnow(m,i,j)+Aice(m,i,j)*Sublim*dt
	  end if

c	if(i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, evaporation', m, Hsnow(m,i,j)
c	end if

      else   ! Hsnow < Hsmin

c     1.1.2. No snow.


c     Method by Ridders.
c     T0 and T1 - initial estimates, f(T0)*f(T1)<0.
      T0= -100. !!!TA(i,j)-50.

      albedoi=F_ai(T0,Hi)
c     Shortwave radiation.
*      Qrad= (1.-aI0*EXP(-hi/150.))*(1.-albedoi)*SW(i,j) 
      Qrad= (1.-aI0)*(1.-albedoi)*SW(i,j)       


c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T0)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Ei*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F0= roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0) +c2*(TFC-T0)
     *   +RLW + Qrad + LH

c     Initial estimate of T1.

***      L=0
***	F1=F0

***      do while (F1*F0 .GT. 0. .and. L.LT.10)
***	L= L+1
***	T1= T0 +20.*REAL(L)

      T1= 20.

      albedoi=F_ai(T1,Hi)
c     Shortwave radiation.
*      Qrad= (1.-aI0*EXP(-hi/150.))*(1.-albedoi)*SW(i,j) 
      Qrad= (1.-aI0)*(1.-albedoi)*SW(i,j)       


c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T1)/(1+0.00412*T1)+0.00422*T1+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T1,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T1)**4)*(0.39-0.005*sqrt(epsa))*
     &     (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T1-Ta(i,j))*((273.16+T1)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Ei*sigma*((273.16+T1)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T1)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F1= roa*Cpa*CDH(T1,TA(i,j))*wind*(TA(i,j)-T1) +c2*(TFC-T1)
     *   +RLW + Qrad + LH


***	end do 

      if(f0*f1 .GT. 0.) then
cc	write(*,*) 'Ice. Initial estimates are bad!!!'
cc	write(*,*) i,j,m,f0,f1,t0,t1

	T2= 20. !  No solution in physically reasonable interval

	else


c     Iterations to determine surface temperature

      iter=0
	F2=F0

	do while (ABS(F2).GT.accur.AND.iter.LT.miter)

	iter= iter+1

	T05= 0.5*(T0+T1)

      albedoi=F_ai(T05,Hi)
c     Shortwave radiation.
*      Qrad= (1.-aI0*EXP(-hi/150.))*(1.-albedoi)*SW(i,j) 
      Qrad= (1.-aI0)*(1.-albedoi)*SW(i,j)       


c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T05)/(1+0.00412*T05)+0.00422*T05+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(T05,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T05)**4)*(0.39-0.005*sqrt(epsa))*
     &(1.-0.8*cloud(i,j))-4.*Ei*sigma*(T05-Ta(i,j))*((273.16+T05)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Ei*sigma*((273.16+T1)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T1)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

c     Total heat flux:
      F05= roa*Cpa*CDH(T05,TA(i,j))*wind*(TA(i,j)-T05) +c2*(TFC-T05)
     *   +RLW + Qrad + LH


	T2=T05 
     *  +(T05-T0)*SIGN(1.,F0-F1)*F05/SQRT(max(1.e-12,F05**2-F1*F0))

      albedoi=F_ai(T2,Hi)
c     Shortwave radiation.
*      Qrad= (1.-aI0*EXP(-hi/150.))*(1.-albedoi)*SW(i,j) 
      Qrad= (1.-aI0)*(1.-albedoi)*SW(i,j)       


c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T2)/(1+0.00412*T2)+0.00422*T2+2.)
c     Specific humidity at the surface Qs:
	Qs=0.622*EPSs/Pa(i,j)
c     Latent heat:
	LH=-roa*QLi*CDL(T2,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T2)**4)*(0.39-0.005*sqrt(epsa))*
     &    (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T2-Ta(i,j))*((273.16+T2)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Ei*sigma*((273.16+T2)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T2)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))

      F2= roa*Cpa*CDH(T2,TA(i,j))*wind*(TA(i,j)-T2) +c2*(TFC-T2)
     *   +RLW + Qrad + LH


	if(iter.eq.miter.and.abs(f2).GT.accur)then
	write(*,*) '============================='
	write(*,*)'Ice temp.!', T0,T1,T05,T2,F0,F1,F05,F2,i,j,m
	write(*,*) '============================='
c      if(kprint.eq.1) write(*,*)iter, T0,T1,T2,F2
	end if


      IF(F05*F2 .LT. 0.) then

	T0= T05
	T1= T2
	F0= F05
	F1= F2

	else

	IF(F0*F2 .LT. 0.) then
	T1=T2
	F1=F2
	end if

	IF(F1*F2 .LT. 0.) then
	T0=T2
	F0=F2
	end if

	end if

	end do

	end if

      Tice(m,i,j)= T2

c     End of iterations to determine surface temperature.


      if( T2 .GE. Tmelt_i) then
      T0=Tmelt_i
      Tice(m,i,j)= T0

c     New energy balance at ice surface.
c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*T0)/(1+0.00412*T0)+0.00422*T0+2.)
c     Specific humidity at the surface Qs:
	Qs=0.622*EPSs/Pa(i,j)
	LH=-roa*QLi*CDL(T0,TA(i,j))*wind*(Qs-Q2m(i,j))
c     Net longwave radiation (Rosati & Miyakoda,1988)
      RLW= -Ei*sigma*((273.16+T0)**4)*(0.39-0.005*sqrt(epsa))*
     &    (1.-0.8*cloud(i,j))-4.*Ei*sigma*(T0-Ta(i,j))*((273.16+T0)**3)
c     Net longwave radiation: Koenig-Langlo&Augstein, 1994, Maykut&Church, 1973
c      RLW= -Ei*sigma*((273.16+T0)**4) + sigma*(TAK**4)*
c     &      (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     Marshunova, 1966. 
c      RLW= -Ei*sigma*((273.16+T0)**4) + 
c     &sigma*(TAK**4)*(0.67+0.005*sqrt(epsa))*(1.+0.29*cloud(i,j))


*      Qrad=(1.-aI0*EXP(-hi/150.))*(1.-F_ai(Tmelt_i,Hi))*SW(i,j)
      Qrad=(1.-aI0)*(1.-F_ai(T0,Hi))*SW(i,j)

      Qai= roa*Cpa*CDH(T0,TA(i,j))*wind*(TA(i,j)-T0) + LH +Qrad
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


c     Ice Sublimation - no flux to ocean!
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

* ----------------------------------------------------------------
*     McPhee, M.G. 1992. Turbulent heat fluxes ib the upper ocean
*     under ice. JGR, 97, 5365-5379.
*     Cdw = Const or hi dependent - Mellor&Kantha 1989
c      Hscale= 100.
c      CDw1= ( 0.4/alog(Hscale*300./(max(5.,0.5*hi))) )**2

	u2= (uice(i,j)-um2(i,j,1))**2
	v2= (vice(i,j)-vm2(i,j,1))**2

*     Friction velocity
	u_star= SQRT(Cdw*(u2 + v2))
	u_star=MAX(5.e-2, u_star)
	
****	u_star=min(100.,SQRT(0.21*Q2turb(i,j,1)))

	CTb=6.0e-3*u_star
**	CTb=7.27e-3                 ! Ebert & Curry, 1993
	
* ----------------------------------------------------------------
* Omstedt, A., and J. S. Wettlaufer, 
* Ice growth and oceanic heat-flux: Models and measurements,
* J. Geophys. Res., 97, 9383–9390, 1992.

c      CTb= 2.e-4*
c     * SQRT((uice(i,j)-um2(i,j,1))**2 +(vice(i,j)-vm2(i,j,1))**2)
* ----------------------------------------------------------------    
	
	ppp= 1.e-5*g*row*Hi
      TFC= TFr(S(i,j,1), ppp) ! Freez T at the depth of ice bottom

      TW= MAX(TFC,T(i,j,1))     ! If T<TFC - ice formation

      Qiw= row*cpw*CTb*(TW-TFC)
      Qi = c2*(Tice(m,i,j)-TFC)

c     Ice bottom mass change. 
cc    dHice= -aice(m,i,j)*dt*(Qiw+Qi)/Qice

c     Heat of fusion is less factor 0.92 because
c     of brine pockets at the ice bottom - 
c     see Los Alamos Ice Model, Bitz & Limpscomb 1999.
c     or Ebert & Curry, 1993
      dHice= -1.087*aice(m,i,j)*dt*(Qiw+Qi)/Qice

c     If too much ice was melted, the remaining heat is going to
c     the very upper ocean layer. Snow gets to water and chills it.
c     WARNING: This is essentially explicit time scheme.

      if( (Hice(m,i,j)+dHice) .LT. 0.)then
        Qdelta= Qice*(-dHice-Hice(m,i,j)) -Qsnow*Hsnow(m,i,j)
        T(i,j,1)= T(i,j,1) + 2.*Qdelta/(hz(1)*row*Cpw)
        dHice= -Hice(m,i,j)   ! Complete Ice Melting
        Aice(m,i,j)=0.
	  dHsnow= -Hsnow(m,i,j) ! Complete Snow Melting
	  Hsnow(m,i,j)=0.
        dHsnowT    = dHsnowT +dHsnow
      end if

      Hice(m,i,j)= Hice(m,i,j) + dHice
      dHiceT= dHiceT + dHice

	if(dHiceT .LT. -1.e1) then
	write(*,*)'base', i,j,m, dHiceT, dHice, aice(m,i,j), 
     *          Qiw,Qi,CTb,TFC,S(i,j,1),TW,Tice(m,i,j),c2
	end if

	else

	Tice(m,i,j)= 0.
	Tsnow(m,i,j)=0.

	end if  ! ai(m) > 0
      end do  ! m, Loop over ice thickness gradations

	do m=1,mgrad

c     Change due to side melting.     

      if(Aice(m,i,j).GT.Aimin) then

      Hi= Hice(m,i,j)/aice(m,i,j)

	ppp= 1.e-5*g*row*0.5*Hi              
      TFC= TFr(S(i,j,1), ppp) ! Freezing T at half-depth of ice bottom
 
      TW= MAX(TFC,T(i,j,1))     ! If T<TFC - ice formation
	Rside= 3.14*1.6e-6*SQRT((TW-TFC)**3)*DT/(0.66*300.) !See NCAR manual
      dAice=  -Rside*Aice (m,i,j)
      dHice=  -Rside*Hice (m,i,j)
      dHsnow= -Rside*Hsnow(m,i,j)

      Hsnow(m,i,j)= Hsnow(m,i,j)+dHsnow
      Hice (m,i,j)= Hice (m,i,j)+dHice
      Aice (m,i,j)= Aice (m,i,j)+dAice

      dHsnowT = dHsnowT +dHsnow
      dHiceT  = dHiceT +dHice

      Q_melt= Qsnow*dHsnow +Qice*dHice
      T(i,j,1)= T(i,j,1) +2.*Q_melt/(hz(1)*row*Cpw)

      else  
   
      Hice (m,i,j) = 0.
      Aice (m,i,j) = 0.
      Tice (m,i,j) = 0.
      Hsnow(m,i,j) = 0.
      Tsnow(m,i,j) = 0.
      
      end if ! Ice Side Melting.

c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, side melting', Hsnow(m,i,j)
c	end if


      if(Hsnow(m,i,j).LT.Hsmin) then
	Hsnow(m,i,j)= 0.
	Tsnow(m,i,j)=-10.
	end if
      if(Hice(m,i,j).LT.Himin) then
	Hsnow(m,i,j)=0.
	Tsnow(m,i,j)=-10.
	Hice(m,i,j) =0.
	Aice(m,i,j) =0.
	Tice(m,i,j) =-10.
	end if

c     Snow to Ice conversion.

      if(Aice(m,i,j). GT. Aimin)then
	Hs=Hsnow(m,i,j)/Aice(m,i,j)
      Hi= Hice(m,i,j)/Aice(m,i,j)

      if(Hs.GT.Hsmin) then
      
c	if(Hs.GT.10.) then   ! Snow aging -check it! Salinity?
c	StoIce=GAMMA*dt*Hs
c	Hs= Hs -StoIce
c	Hi= Hi +StoIce*rosdry/roi
c      dHsnow= -StoIce*aice(m,i,j)
c	dHice = +StoIce*rosdry*aice(m,i,j)/roi
c	end if

c      Hsnow(m,i,j)= Hsnow(m,i,j)+dHsnow
c      Hice(m,i,j) = Hice(m,i,j) +dHice
**      dHsnowT     = dHsnowT +dHsnow ! No salinity change!
**      dHiceT      = dHiceT  +dHice


c     Snow to ice conversion due to flooding.
      hdraft= (rosdry*Hs+roi*Hi)/row
	hflood= MAX(0.,hdraft-Hi)
	Hs= Hs -hflood*roi/rosdry
	Hi= Hi +hflood

      dHsnow= -hflood*aice(m,i,j)*roi/rosdry
	dHice = +hflood*aice(m,i,j)

      Hsnow(m,i,j)= Hsnow(m,i,j)+dHsnow
      Hice(m,i,j) = Hice(m,i,j) +dHice
      dHsnowT     = dHsnowT +dHsnow
      dHiceT      = dHiceT +dHice

	end if         ! Hs>Hsmin
	end if         ! Aice>Aimin

c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*)'Snow, after flooding',Hsnow(m,i,j),hflood,hi,hdraft
c	end if


      end do  ! Loop over ice thickness gradations

c  2. New ice formation.

      dHiceN= 0.

	kb=km2(i,j)
      do k=1,kb

	ppp= 1.e-5* g*row*z(k)
      TFC= TFr(S(i,j,k),ppp) ! Freezing temperature

      if(T(i,j,k).LT.TFC ) then 
      
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

c     In a case of sufficient open water: frazil ice to the open water.
	hhh0= dHiceN/MAX(aimin,aice(0,i,j))

c     Thickness of the new ice

      IF( hhh0 .LE. Hmax(2)) then
	delta=dHiceN/Href

      Hice(1,i,j)=Hice(1,i,j)+dHiceN
      aice(1,i,j)= aice(1,i,j)+ MIN(delta, aice(0,i,j))
      aice(0,i,j)= MAX(0.,aice(0,i,j)- delta)

c     New ice temperature - for T in situ only!
      Q= (Tice(1,i,j)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
      Tice(1,i,j)=Q/Hice(1,i,j)

      Tsnow(1,i,j)=-10.

      else

c     New Ice mass spreaded unifirmly over entire cell
c     if no sufficient open water
      aice(1,i,j)= aice(1,i,j)+ aice(0,i,j)
      aice(0,i,j)= 0.

      do m=1,mgrad
	Hice(m,i,j)=Hice(m,i,j)+dHiceN*Aice(m,i,j)
	end do

c     New ice temperature - for T in situ only! Should be corrected
      Q= (Tice(1,i,j)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
      Tice(1,i,j)=Q/Hice(1,i,j)
      Tsnow(1,i,j)=-10.

	end if

      dHiceT= dHiceT + dHiceN

	end if ! New ice > Himin

c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, return', Hsnow(m,i,j)
c	end if


      return
      end

C---------------------------------------------------------------------
      SUBROUTINE FACTOR(ND,A,B,C,F,RKSI,II,JJ)
C---------------------------------------------------------------------
C     UP-DOWN : A(N)*Fi(N-1) +B(N)*Fi(N) +C(N)*Fi(N+1)=F(N) ; N=II,JJ
C     A(II) AND C(JJ) ARE NOT USED (I HOPE SO), SO THEY ARE SET = 0.
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

      real function TFr(S,p)

*     Version 26.11.2012.
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

      s8=DBLE(max(0.,s))
	p8=DBLE(p) !!!- 10.1325d0
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

      Real Function F_ai(Ti,h)
*     Mellor and others
*	parameter (aow=0.10, ami=0.55, c11=0.44, c12=0.075)

	parameter (aow=0.10) ! AOMIP
cc	parameter (aow=0.06) ! NCAR, LANL

      T=min(0.,Ti)

	if(T .LT. -1.0)then  ! -1C - see BoulderIce
cc  	F_ai= MAX(aow, c11*SQRT(SQRT(0.01*h))+aow)
cc      ai= MIN(1.,F_ai)
cc      ai= 0.62
cc      ai= 0.60  ! AOMIP, approx. 3m.
cc      ai= 0.70  ! New AOMIP, 15.05.2003 approx. 4m.
cc      ai= 0.65  ! Boulder Ice CCSM3 (Ice version 4, 2002) 
      ai= 0.575  ! LANL CICE 4.1 

	if(h.LE.50.)then  ! see BoulderIce or CICE 4.1.
	F_ai= aow +(ai-aow)*h/50.
	else
	F_ai= ai
	end if

	else

cc      ai= 0.70  ! New AOMIP, 15.05.2003 approx. 4m.
cc      ai= 0.65  ! Boulder Ice CCSM3 (Ice version 4, 2002)
cc      ai= 0.60  ! AOMIP, approx. 3m.
      ai= 0.575  ! LANL CICE 4.1

	if(h.LE.50.)then  ! see BoulderIce - this is a good approx.
	F_ai= aow +(ai-0.075*(T+1.)-aow)*h/50. 
cc	F_ai= aow +(ai-0.1  *(T+1.)-aow)*h/50. 
	else
	F_ai= ai -0.075*(T+1.)
cc	F_ai= ai -0.1  *(T+1.)
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
	
	Real Function F_as(Ts)

      T= min(Ts,0.)

	if(T .LT. -1.0)then   ! -1C - see BoulderIce
cc	F_as=0.75  ! Flato&Brown 1986
	F_as=0.80  ! Old Aomip, close to CCSM3
cc	F_as=0.82  ! SHEBA
cc	F_as=0.85  ! Pinto, 1999, SHEBA, LANL CICE 4.1
cc	F_as=0.81  ! New AOMIP, 15.05.2003
	else
cc	F_as=0.65  ! F&B 1986
cc	F_as=0.70  ! Old AOMIP SHEBA 2001 -SHCE
cc	F_as=0.72  ! Boulder
cc	F_as=0.80 -0.125*(T+1.0)  ! AOMIP modified
cc	F_as=0.85 -0.125*(T+1.0)  ! SHEBA + LANL CICE 4.1
	F_as=0.80 -0.1*(T+1.0)  ! SHEBA + LANL CICE 4.1
cc	F_as=0.85 -0.1*(T+1.0)  ! Pinto,1999, SHEBA+Boulder Ice, CCSM3
cc	F_as=0.77  ! New AOMIP, 15.05.2003
	end if

	return
	end 

	Real Function CDH(T,TA)
*     version 17.12.2009. Regularisation for continuity.
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
	
	Real Function CDL(T,TA)
*     version 17.12.2009. Regularisation for continuity.
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
