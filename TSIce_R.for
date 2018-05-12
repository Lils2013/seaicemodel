      SUBROUTINE TSICE

c     Version 15.01.2016 for ice thickness gradations.
c     Potential temperature and Freezing point as function on depth.

c     Nonuniform heat flux ice-ocean
c     Spatially Varied Cloudiness
c     AOMIP forcing - Net LW radiation.
c     Penetration of the Short Wave Radiation into the Ocean.
c     Friction velocity as in Los Alamos model?
c     Evaporation and snow/ice sublimation. 

**    Temperature flux due to Rain and Rivers

      INCLUDE 'Slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
	REAL dHiceT, dHsnowT, DHI
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
      
      if(aice(0,i,j).LT.0.) aice(0,i,j)=0.      

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
      FM(1)= T(i,j,1)  +CA*(Rain*max(0.,TA(i,j))+Psnow*min(0.,TA(i,j)))

*     Atm Snow melting on the open water    
      FM(1)= FM(1) -CA*Psnow*Qsnow/(row*cpw)
      
c     Heat to melt the snow from ridging
      Qdelta= -Qsnow*STWR(i,j)/rosdry 
      FM(1)= FM(1) +2.*Qdelta/(hz(1)*row*Cpw)

*     Heat fluxes through ocean and ice.

*     1. Sensible heat at open water.
      BM(1)=BM(1)+CA*roa*Cpa*CDH(T(i,j,1),TA(i,j))*
     &WIND        *Aice(0,i,j)/(row*cpw)
      FM(1)=FM(1)+CA*roa*Cpa*CDH(T(i,j,1),TA(i,j))*
     &WIND*TA(i,j)*Aice(0,i,j)/(row*cpw)


*     3. Incoming radiation .
*     3.1. Short wave at open water.

	Qrad= Aice(0,i,j)*(1.-aopw)*SW(i,j)

      aI0= 0.18 +0.17*cloud(i,j)    ! Grenfell & Maykut, 1977. 
cc	aI0= 0.4*aI0                  ! To account volume heating for 0-layer

      FM(1)= FM(1)+ CA*Qrad
     #*(1.-(Ra*dzi1*(1.-exps1(2))+(1.-Ra)*dzi2*(1.-exps2(2)))/hz(1))
     #       /(row*cpw)
     
*     2. Sensible heat at ice-water interface. 

*     McPhee, M.G. 1992. Turbulent heat fluxes in the upper ocean
*     under ice. JGR, 97, 5365-5379.
*     Q=rw*Cp*CDw*CoefW*dT, CoefW=SQRT(CDw*((uice-u)**2) -friction velocity.
*     Skin velocity W=1-3 cm/s, Mellor&Kantha 0.5.
*     Los Alamos sea ice model - cut off 0.075 -> dU=1cm/s
*     Drag coefficient


*     Cdw= const or hi dependent - Mellor&Kantha 1989
c      Hscale= 100.
c      CDw1=( 0.4/alog(Hscale*300./(max(5.,0.5*hi))) )**2

*     Friction velocity
      u_star=SQRT(CDw)*delta_u(i,j)
	u_star=MAX(1.e-2, u_star)

	CTb=6.0e-3*u_star           ! CICE
**	CTb=7.27e-3                 ! Ebert & Curry, 1993
**    CTb=1.5 e-3                 ! Mellor&Kantha, 1989
* ----------------------------------------------------------------
* Omstedt, A., and J. S. Wettlaufer, 
* Ice growth and oceanic heat-flux: Models and measurements,
* J. Geophys. Res., 97, 9383–9390, 1992.

c      CTb= 2.e-4*
c     * SQRT((uice(i,j)-um2(i,j,1))**2 +(vice(i,j)-vm2(i,j,1))**2)
* ----------------------------------------------------------------    
     

      do mg=1,mgrad
	Ai=MAX(Aimin,Aice(mg,i,j))
	Hi=Hice(mg,i,j)/Ai

	ppp= 1.e-5*g*row*Hi                   
      TFC= TFr(S(i,j,1), ppp) ! Freezing point at the depth of ice bottom

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
c      RLWf= 3.*Ew*sigma*TK**4 + sigma*(TAK**4)*
c     &     (0.5893+5.351e-3*sqrt(epsa)+0.2176*SQRT(cloud(i,j)**3))
c     &     -273.16*RLWb
c      FM(1)= FM(1)+ CA*Aice(0,i,j)*RLWf/(row*cpw)
c      BM(1)= BM(1)+ CA*Aice(0,i,j)*RLWb/(row*cpw)

*     4. Latent heat for open water.
c     Saturation vapor pressure EPSs:
	Qs=0.622*EPSs/Pa(i,j)
      Evapor= -roa*CDL(T(i,j,1),TA(i,j))*wind*(Qs-Q2m(i,j))  

      LH= QLw*Evapor
      if (i .eq. 16 .and. j .eq. 41) then
       ! write(*,*) "Pa(i,j)", Pa(i,j), EPSs, TC
       ! write(*,*) "wind", wind, Qs, Q2m(i,j), CDL(T(i,j,1),TA(i,j))
       ! write(*,*) "Evapor", Evapor, T(i,j,1), TA(i,j)
      end if
	Evapor= Evapor/row

      
      FM(1)= FM(1)+CA*Aice(0,i,j)*LH/(row*cpw)

*     ------------------- Deep water ---------------
      DO K=2, KB -1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

      AM(K)= -CA*CG( n)*AZT(i,j,k-1)/HZ(K-1)
      if (i .eq. 11 .and. j .eq. 1 .and. K .eq. 2) then
      !  write(*,*) "AM(1)", CA, CG( n), AZT(i,j,k-1), HZ(K-1)
        !write(*,*) "AZTT", AZT
      end if
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

      if (i .eq. 16 .and. j .eq. 41) then
        Do k=1,kb
            !write(*,*) "T", T(i,j,k)
        end do
        !write(*,*) "AM", AM
        !write(*,*) "BM", BM
        !write(*,*) "CM", CM
        !write(*,*) "FM", FM  
      end if
      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)
      Do k=1,kb
      
      T(i,j,k)= RKSI(K)
      if (T(i,j,k) .ne. T(i,j,k)) then
        !write(*,*) "QWEQWE", AM, BM, CM, FM
        !write(*,*) T(i,j,k), "TTT", i,j,k
        end if
      End do
      if (i .eq. 16 .and. j .eq. 41) then
        Do k=1,kb
            !write(*,*) "T", T(i,j,k)
        end do
      end if

c      if(i.eq.17 .and. j.eq.24) then
c      write(*,*)'T thermo:', TA(i,j), T(i,j,1), LH, RLWf
c      end if

*     River + Precipitation - Evaporation.
      PME(i,j)= dt*(River(i,j) +Rain +Psnow +Evapor*Aice(0,i,j))
      if (i .eq. 16 .and. j .eq. 41) then
      !write(*,*) "IMPO", Aice(0,i,j), River(i,j), Rain, Psnow, Evapor
      !write(*,*) "rtant", T(i,j,1),TA(i,j)
      end if
*     Sea ice thermodynamics.
*     dHsnowT, dHiceT are total melted snow and ice mass,
*     used for vertical velocity calculation.
*     Ocean temperature is calculated already
*     (Except the case of complete ice melting and 
*     lateral melting - see IceTH).

      CALL IceTH(i,j,Wind,dHiceT,dHsnowT)
      if (dHiceT .ne. dHicet) then
        write(*,*) "ij", i, j, dHicet
      end if
      if (dHiceT > 100.) then
        write(*,*) dHiceT, "dHiceT"
      endif
      if (dHiceT < -100.) then
        write(*,*) dHiceT, "dHiceT"
      endif
      !write(*,*) dHiceT, "dHiceT asasdasdasd"
      !write(*,*) dHsnowT, "dHsnowT"

*     Salinity changes if only snow melts dHsnowT<0.
      DHS= MIN(dHsnowT,0.)
      DHI= MIN(dHiceT,0.)

*     River + Precipitation - Evaporation + Ice + Snow

ccc   PME_DZ(i,j)= -roi*dHiceT/row ! Ice doesn't change the dyn. level
      if (i .eq. 16 .and. j .eq. 41) then
       ! write(*,*) "IMPORTANT!!", PME(i,j), roi, DHI, rosdry, DHS, row
      end if
      PME(i,j)   = PME(i,j)-(roi*DHI+rosdry*DHS)/row
      
      
c     Ridging
      PME(i,j)= PME(i,j)+STWR(i,j)

      

*     Vertical distribution of salinity.

*     -------------- Ocean surface -----------------
      CA=2.*DT/hz(1)
      AM(1)= 0.
      CM(1)= -CA*AZS(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) -River(i,j)*dt/hz(1) ! River_S = 0.5*Ocean_S!
      FM(1)= S(i,j,1)  -2.*Sice*dHiceT*roi/(row*hz(1))
      
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
      if (RKSI(K) > 100000.) then
      write(*,*) RKSI(K), "RKSI(", K, ")", i, "i", j, "j"
      write(*,*) dHiceT, "dHiceT", S(i,j,k)
      endif
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
c     Version #11.3 - 05.02.2016.
c     Potential Temperature and freezing point as a function of depth.
c     Method by Ridders to Solve the Nonlinear Equation.
c     C.J.F. Ridders, IEEE Transactions on Circuits and Systems, 
c     v. CAS-26, p.979, 1979.
c     AOMIP Net LW radiation.
c     Flooding and aging with correspondent freshwater flux to ocean.

**    Modifications of lateral melt rate.

      include 'Slo2.fi'
	REAL dHiceT, dHsnowT
	real(8) :: tzero  = 273.15d0
      real(8), dimension(0:maxlay) :: rhs, a, b, diagon, wmesh, wold, 
     & a1,b1,c1,d1
      real(8) tseafrz, q0, fsens, CTb, flat
      real(8)  sstz, bshf, Fnet, dzf, fwf, 
     & es, zrchu1,zrchu2,zref,zssdqw
      real(8), dimension(0:maxlay) :: tiold, tinew, timid, ziold
      integer ni,ns
      real(8) :: hslim=0.0005d0
      integer    lice_top, layer
      real(8), dimension(maxlay) :: kice,
     &        he, heold, henew,
     &        re, tme, dhe, 
     &        em, cm, qm, km, rad, sali
      real(8), parameter :: rhosno    = 330.0_8,
     &        cp_ice    = 2.062e+03_8,
     &        cp_wat    = 3.99e+03_8,
     &        mlfus     = 3.335e+05_8,
     &        cond_sno  = 0.31_8,
     &        cond_ice  = 2.034_8, 
     &        emi       = 0.99_8,  
     &        stefa     = 5.6697e-08_8,
     &        rhowat    = 1025.0_8 ,
     &        latvap    = 2.5d6, 
     &        rhoice    = 917.0_8,
     &        ai2       = 21.8746,
     &        bi        =-265.5,
     &        rad_kappa_i = 0.8
       real(8) ftot2, sume2, flux2, pi, energy_snow_melt, 
     &       sublimation_speed, bottom_speed, rom
      real(8) dzi(maxlay), dzi_sw(maxlay)
      real(8) dzs(maxlay), dzs_sw(maxlay)
      real(8), dimension  (maxlay) :: dziold
      real(8) blhf, 
     &     submer_is, hdraft, sice, dhi, dhs, dq, 
     &   m1, m2, k1, theta_ther, eps, topflux, sum0, alb, estef, hmelt, 
     &   ftot,Cbrine,zn,salin,
     &        Tside,deltaT,wlat,rside,
     &  dhi_underice, dhi_openwater, dho, hi_accr, capa, efus, tt, sim,
     &        qsr, snow_precip, rain_precip, w1, w2, w3, w4, xi, slope,
     &        hiold, hsold, fac, maxvel, slop2,
     &        hsnew, hinew, tm, sn, sume, heat_melt_to_ocean, 
     &        energy_ice_growth, tt1, tt2, xis, xib       
      integer :: iteration, maxiter = 50
      real(8) tsu, tbo, ttt, ssss, aavg
      real(8) netlw, dwnlw, pres
      real(8) fac_transmi, swrad
      logical fixed_tb,melt_ts,melt_tb,panic
      real(8) Fprec, ssnow ! heat flux due to precipitation [W/m2]
      real(8) tair, qair, uair
      logical :: energy_check=.false., debug=.false.
      logical please_stop, melted, thin_snow_active
      real(8) :: cheat=4.2174d3*1.025d3
      real(8) dh_sni, deltaz_i_phy(0:maxlay)
      real(8) zi(0:maxlay), radtr_s(0:maxlay), radab_s(0:maxlay)
      real(8) radtr_i(0:maxlay), radab_phy_i(0:maxlay)
      real(8) ftrice, ab, isnow, zzs, zraext_i, zdummy
      include 'Tparm.fi'

      please_stop = .false.
      melted = .false.
      ni=nlice
      ns=ni+nlsno
           
      do k=0,ni-1
      zi(k) = dble(k)/dble(nlice)
      enddo
      zi(ni) = 1.d0
      do k=ni+1,ns-1
      zi(k) = 1.d0 + dble(k-ni)/dble(nlsno)
      enddo
      zi(ns) = 2.d0
      
      do k=1,ns
          dzi(k)=zi(k)-zi(k-1)
      enddo
      ziold =  zi
      dziold = dzi
      
      
      
      dHiceT = 0.
      dHsnowT= 0.

	TAK=Ta(i,j)+273.16

c     Accuracy - 1000 erg/cm2 = 1.0 W/m2 
      accur= 10.
      miter= 50

c     Atmosphere vapor pressure epsa:
      epsa=Q2m(i,j)*Pa(i,j)/0.622

c     1. Old snow and ice.
      do m=1,mgrad
      
      bottom_speed       = 0.d0
      energy_snow_melt   = 0.d0
      heat_melt_to_ocean = 0.d0
      melt_ts=.false. 
      thin_snow_active =.false. 
      fwf=0.d0
      theta_ther=1.d0 
      snosub = 0 !!!!!!!!!!!!!!!!!!!
      fac_transmi = 1

      dHsnow =0.
      dHice  =0.
      bshf=0.d0

      if (i .eq. 14 .and. j .eq. 1 .and. m .eq. 13) then
        ! write(*,*) aice(m,i,j)
         endif
         if (aice(m,i,j) .ne. aice(m,i,j)) then
         !write(*,*) m, i, j
         !write(*,*) aice(m,i,j)
         end if
      if( aice(m,i,j) .GT. Aimin) then

      Hi= (Hice(m,i,j)/aice(m,i,j))/100.
	!ppp= 1.e-5*g*row*Hi                   
      !TFC= TFr(S(i,j,1), ppp) ! Freezing point at the depth of ice bottom

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
!!!
!!!!!      IF(TA(i,j).LT.-1.) then
!!!!!      snow_precip = row*Pr(i,j)*Aice(m,i,j)/rosdry
!!!!!	deltaHS = dt*row*Pr(i,j)*Aice(m,i,j)/rosdry
!!!!!	end if
!!!!!
!!!!!	IF(TA(i,j).GE.-1..AND.TA(i,j).LE.5.) then
!!!!!	deltaHS = dt*row*Pr(i,j)*Aice(m,i,j)*
!!!!!     &                     (1.-0.16666*(TA(i,j)+1.))/rosdry
!!!!!      snow_precip = row*Pr(i,j)*Aice(m,i,j)*
!!!!!     &                     (1.-0.16666*(TA(i,j)+1.))/rosdry
!!!!!	end if


c      if(Hsnow(m,i,j) .LT. hsmin .and. deltaHS .gt. hsmin) then
*     Temperature of thin new snow
c      Tsnow(m,i,j) = TA(i,j)-5.
c      end if

!!!!      Hsnow(m,i,j) = Hsnow(m,i,j)+deltaHS


c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, after snowfall', m, Hsnow(m,i,j), Tsnow(m,i,j)
c	end if
      
      
      do k=1,ni
          tiold(k) = TiceFE(m,i,j,k)
      enddo
      do k=ni+1,ns
        tiold(k) = TsnowFE(m,i,j,k-ni)
        if (Hsnow(m,i,j) .gt. 0. .and. tiold(k) .eq. 0. .and. 
     &k .ne. ns) then
            tiold(k) = -10 !!!!!!1Tsnow(m,i,j)
        end if
      enddo
      kice(1:ni)   =cond_ice
      kice(ni+1:ns)=cond_sno
      
      seasal = 30 !S(i,j,1)
      tseafrz=-0.054*seasal
      
      tiold(0) = tseafrz
      tinew = tiold
      timid = tiold
      
      re(1:ni   )=rhoice
      re(ni+1:ns)=rhosno
      
      do k=1,ni
          sali(k) = 30 !S(i,j,1)
      enddo
      DO k=ni+1,ns
          sali(k) = 0.d0
      ENDDO
            
      do k=1,ns
          ssss=sali(k)
          tme(k)=Tfreeze1(ssss) !-0.054*sali(k)
      enddo
      
      do k=1,ni
        tt = 0.5d0*(tiold(k-1)+tiold(k)) 
        ssss=sali(k)
        kice(k) = func_ki(ssss,tt)
        if (i .eq. 13 .and. j .eq. 6) then
          !write(*,*) kice(k), "----", sali(k), tt, tseafrz, seasal
        end if 
      enddo
      
      !---------------------------------------------------------------------
      ! do computation if only enough ice
      !---------------------------------------------------------------------

      !if (hi.gt.hminice) then
      hsold=(Hsnow(m,i,j)/(aice(m,i,j)))/100.
      hiold=(Hice(m,i,j)/(aice(m,i,j)))/100.
      
      heold(   1:ni)=hiold*dzi(   1:ni)
      heold(ni+1:ns)=hsold*dzi(ni+1:ns)
      he   =heold
      henew=he
      
      wmesh(0:ns)=0.d0
      dhe = 0.d0
      dhs = 0.d0
      
      do k=1,ni
       dzi_sw(ni-k+1) = dziold(k) * hiold
      enddo
      do k=ni+1,ns
        dzs(ns-k+1) = dziold(k) * hsold
        if (dzs(ns-k+1) .eq. 0.) then
            !write(*,*) "dzs(ns-k+1)", dziold(k), hs
            end if
      enddo
      
      lice_top=ni
      if (hsold > hslim ) then
          lice_top=ns
          isnow = 1.
          else 
          isnow = 0.
      endif
            
      !---------------------------------------------------------------------
      ! latent heat calculation, sublimation processes and precipitation
      !---------------------------------------------------------------------
c     If snow is thin - one can see ice
c     New ice temperature.
      ti = tiold(lice_top)
      tsu    = tiold(lice_top) + tzero
      albedoi=F_ai(ti,Hi*100)

      sf= hsold/(hsold+0.02)
      
      Qrad= ab*( sf*(1.-0.85)+
     &       (1.-sf)*(1.-albedoi) )*SW(i,j)/10000.
     
      ab = 1.0 - ( 1.0 - isnow ) * 0.30 - isnow * 0.15
      !-----------------------------------------------------
      ! Solar radiation transmitted below the surface layer
      !-----------------------------------------------------
      ftrice      =  (SW(i,j)) * ( 1.0 - ab)*( sf*(1.-0.85)+
     &       (1.-sf)*(1.-albedoi) ) / 10000.
      radtr_s(0) =  ftrice
      zzs = 0.0
      DO layer = 1, nlsno
      if (ftrice .gt. 0) then
        !write(*,*) zzs, "zzs"
        end if
         zzs = zzs + dzs(layer)
!         radtr_s(layer) = radtr_s(0) * exp( - 5.8*( MAX( 0.0 ,
         radtr_s(layer) = radtr_s(0) * exp( - 10.8*( MAX( 0.0 ,
     &                    zzs ) ) )
         radab_s(layer) = radtr_s(layer-1) - radtr_s(layer)
      END DO
      
      radtr_i(0) = isnow * radtr_s(nlsno)  + ( 1. - isnow ) * ftrice
!!      zzs = 0.0
!!      DO layer = 1, nlice
!!            if (ftrice .gt. 0) then
!!        !write(*,*) zzs, "zzs"
!!        end if
!!         zzs = zzs + dzi_sw(layer)/100
!!!         radtr_i(layer) = radtr_i(0) * exp( - 1.4*( MAX( 0.0 ,
!!!     &                    zzs ) ) )
!!           radtr_i(layer) = radtr_i(0) * exp( - 0.80*( MAX( 0.0 ,
!!     &                    zzs ) ) )
!!         radab_phy_i(layer) = radtr_i(layer-1) - radtr_i(layer)
!!      END DO
      
       DO layer = 1, nlice
         ! extinction coefficient
         ! layer thickness
         zdummy    = radtr_i(layer-1) * 
     &               EXP ( - rad_kappa_i * dzi_sw(layer)) *
     &               dzi_sw(layer)

         ! physicallay absorbed
         radab_phy_i(layer) = ( rad_kappa_i )
     &                      * zdummy 
         radtr_i(layer)     = radtr_i(layer-1) - radab_phy_i(layer)
      END DO
      
      
      sumrad=0.d0
      DO k=lice_top,ni+1,-1
         Rad(k) = radab_s(ns-k+1) !k/20. !swradab_s(ns-k+1)
        if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10
     &  .and. Rad(k) .ne. 0.) then
        write(*,*) "Rad(k) snow", Rad(k), k
      end if
      ENDDO
      DO k=ni,1,-1
         Rad(k) = radab_phy_i(ni-k+1) !k/20. !swradab_i(ni-k+1)
         if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10
     &  .and. Rad(k) .ne. 0.) then
        write(*,*) "Rad(k) ice", Rad(k), k
      end if
      ENDDO
      do k=1,lice_top
         sumrad = sumrad + rad(k)
      enddo
            
      do k=1,ns
          tt1 = 0.5d0*(tiold(k-1)+tiold(k))
          ttt = tme(k)
          qm(k) = func_qmelt(ttt,tt1)
          em(k) = func_el(ttt,tt1)
      enddo

      
c     Saturation vapor pressure EPSs:
      EPSs=10.**( (0.7859+0.03477*ti)/(1+0.00412*ti)+0.00422*ti+2.)
c     Specific humidity at the surface Qs.
	Qs=0.622*EPSs/Pa(i,j)
c     Latent Heat:
	LH=-roa*QLi*CDL(tsu,TA(i,j))*wind*(Qs-Q2m(i,j))
      ! net longwave radiative flux
      
      epsa=Q2m(i,j)*Pa(i,j)/0.622
      
      RLW= (-emi*stefa*((273.16+ti)**4)*(0.39-0.005*sqrt(epsa))*
     &   (1.-0.8*cloud(i,j))-4.*emi*stefa*(ti-Ta(i,j))*((273.16+ti)**3))
      
!      dwnlw=200.
!        netlw = emi*(dwnlw - stefa*tsu*tsu*tsu*tsu)
        netlw = RLW
      
      ! pressure of water vapor saturation (Pa)
      es         =  611.d0*10.d0**(9.5d0*(tsu-273.16d0)/(tsu-7.66d0))
            ! intermediate variable
      q0  = 0.622*6.11/1013.*exp(min(ai2*ti/(ti-bi),10.))
      zssdqw     =  q0*q0*pres/ 
     &              (0.622d0*es)*log(10.d0)*9.5d0*
     &              ((273.16d0-7.66d0)/(tsu-7.66d0)**2.d0)
     
      ! derivative of the surface atmospheric net flux
      dzf    =  4.d0*emi*stefa*((tsu)**3) +
     &     + rho2*cp*CDH(tsu,TA(i,j))*(wind/100.)
     &     +rho2*lv*CDL(tsu,TA(i,j))*(wind/100.)*zssdqw
     
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
          !write(*,*) "tsu",  4.d0*emi*stefa*((tsu)**3), 
      !&     rho2,cp,CDH(tsu,TA(i,j)),(wind/100.),
      !&       rho2*lv*CDL(tsu,TA(i,j))*(wind/100.)*zssdqw
        end if
      ! surface atmospheric net flux
!      if (tsu-(TA(i,j)+tzero) .lt. -5.) then
!        fsens=-rho2*cp*CDH(tsu,TA(i,j))*(wind/100.)*(-5)
!      else if (tsu-(TA(i,j)+tzero) .gt. 5.) then
!        fsens=-rho2*cp*CDH(tsu,TA(i,j))*(wind/100.)*(5)
!        else
!       fsens=-rho2*cp*CDH(tsu,TA(i,j))*(wind/100.)*(tsu-(TA(i,j)+tzero))
!      end if
      fsens=-rho2*cp*CDH(tsu,TA(i,j))*(wind/100.)*(tsu-(TA(i,j)+tzero))
      flat= rho2*lv*CDL(tsu,TA(i,j))*(wind/100.)*(q0-Q2m(i,j))
      flat   =  MIN( -flat , 0.d0 )
      Fnet=Qrad + (netlw + fsens + flat)/10
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        write(*,*) "Ta(i,j)", Ta(i,j)
        write(*,*) "tiold(lice_top)", tiold(lice_top)
        write(*,*) "SW", SW(i,j) * ab, SW(i,j), ab
        write(*,*) "Qrad", Qrad, albedoi, sf
        write(*,*) "fsens", fsens,
     &roa*Cpa*CDH(tsu,TA(i,j))*wind*(tzero+TA(i,j)-tsu)
        write(*,*) "netlw", netlw, RLW
        write(*,*) "flat", flat, q0, Q2m(i,j)
        write(*,*) "hiold, hsold", hiold, hsold
        write(*,*) "TSNOW", Tsnow(m,i,j)
      end if
      !fac_transmi * swrad + netlw + fsens + flat
      !roa*Cpa*CDH(tsu,TA(i,j))*wind*(TA(i,j)-tsu)+RLW + Qrad +LH
      
      !write(*,*) "fsens", fsens, TA(i,j)
      
      !---------------------------------------------------------------------
      ! accumulation at the surface
      !---------------------------------------------------------------------
      if (TA(i,j) < 0.0 ) then
        snow_precip = row*Pr(i,j)*Aice(m,i,j)/(rosdry*10.)
        rain_precip = 0.d0
      else
        snow_precip = 0.d0
        rain_precip = Pr(i,j)*Aice(m,i,j)/10.
      endif
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        write(*,*) "precip",snow_precip,rain_precip,Pr(i,j),Aice(m,i,j)
      end if

   !   fwf = fwf + rain_precip ! add liquid precipitation to the freshwater flux to ocean

      !---------------------------------------------------------------------
      ! need to get rid of snow if too thin during melting/sublimating
      ! thin snow case
      !---------------------------------------------------------------------
       if ( hsold <= hslim ) then
       
       if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        write(*,*) "hsold <= hslim", hsold, hslim 
      end if
      ! then 1- compute total energy in snow
      !      2- if Fnet>0, we assume that we can use some of the heat flux to melt the thin snow
      !      3- if Fnet<=0,we assume that some of the sublimation can be used to sublimate the snow

        sume2 = 0.d0
        do k=ni+1,ns
          sume2 = sume2 + em(k) * re(k) * heold(k)
        enddo

        if ( Fnet > 0.d0 ) then ! use Fnet if positive

          energy_snow_melt = min( - sume2, Fnet * dt )! energy (negative) required to melt the thin snow (in entirity if Fnet large enough)
          dhs = 0.d0
          if ( hsold > 0.d0 ) then
            dhs = energy_snow_melt / sume2 * hsold ! melted snow thickness (negative)
          endif
          snow_precip = 0.d0

        else !sublimation can be used too

          sume2 = sume2 - latvap * re(k) * hsold ! the snow needs to be sublimated, so the total required energy of melting is higher!
          dq=max(sume2,flat*dt*snosub) ! negative energy required to sublimate the thin snow (in entirity if flat large enough, partial otherwise)
          if ( hsold > 0.d0 ) then
            dhs = snow_precip * dt - dq / sume2 * hsold
          else
            dhs = snow_precip * dt
          endif
      ! recalculate the top flux
          Flat =  Flat - dq / dt ! (-dq>0 so this is in effect reducing the latent heat)
          Fnet =  Fnet - dq / dt ! (-dq>0 so this is in effect reducing the latent heat)

        endif

        thin_snow_active = .true.

! new velocity at the snow-air interface
        wmesh(ns) = - dhs / dt
      
!---------------------------------------------------------------------
      endif ! end thin snow
!---------------------------------------------------------------------

      ! sublimation at the surface
      sublimation_speed = snosub * flat/(re(lice_top)
     &*(-em(lice_top)+latvap))

	ppp= 1.e-5*g*row*Hi*100.
       ! Freezing point at the depth of ice bottom

          ! If T<TFC - ice formation
!      TW = -1.5
!      TFC = tseafrz
      
!     Friction velocity
      u_star=SQRT(CDw)*delta_u(i,j)
	u_star=MAX(1.e-2, u_star)

	CTb=6.0e-3*u_star 
!	CTb=7.27e-3 

        TFC= TFr(S(i,j,1), ppp)
        !TFC = Tfreeze1(ssss)
        TW= MAX(TFC,T(i,j,1)) 
      Qiw= row*cpw*CTb*(TW-TFC)
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        write(*,*) "QIW", Qiw, TW, TFC, CTb
!        TW= MAX(TFC,T(i,j,1)) 
!        TFC= TFr(S(i,j,1), ppp)
!        write(*,*) "QIWORIG", row*cpw*CTb*(TW-TFC)
!        write(*,*) "CTb", 6.0e-3*u_star 
      end if
      bshf = -Qiw/2000
      
      !-----------------------------------------------------------------------
      !     Update energy flux due to snow precipitation
      !-----------------------------------------------------------------------
      
      if ( .not.thin_snow_active ) then
          tm    = tme(lice_top)
          tt    = TA(i,j)
          efus  = cp_ice * ( tt - tm ) - mlfus + cp_wat * tm
          Fprec = snow_precip * rhosno * efus
      else ! FD no snow
          Fprec = snow_precip * rhosno * em(ns)
      endif
      
      ! surface flux
      topflux = Fnet
     
      tiold(0) = tseafrz
      
      if (lice_top==ns) then
          wmesh(ns) =  - snow_precip - sublimation_speed
          wmesh(ni) =0.d0
          wmesh(0 ) = bottom_speed
      else
          wmesh(ni) = - sublimation_speed
          wmesh(0 ) = bottom_speed ! ice growth or melt
      endif
      
      fac=0.d0
      do k=1,ni-1
          fac = fac + dzi(k)
          wmesh(k)=wmesh(0 ) * (1.d0-fac) + wmesh(ni) * fac
      enddo
      fac=0.d0
      do k=ni+1,ns-1
          fac = fac + dzi(k)
          wmesh(k)=wmesh(ni) * (1.d0-fac) + wmesh(ns) * fac
      enddo

      wold=wmesh
    !  write(*,*) wmesh
      
      
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
          write(*,*) "theta_ther * dt * dzf", dzf
          write(*,*) "bshf", bshf
       write(*,*) "Fnet", Fnet ,"energy_snow_melt", energy_snow_melt
      end if
      
      do iteration = 1, maxiter ! convergence on temperature and wmesh

      rhs=0.d0
      a=0.d0
      b=0.0d0
      diagon=0.d0
      cm=0.d0
      
      do k=1,lice_top
       !   tt = 0.25d0*(tiold(k-1)+tinew(k-1)+tiold(k)+tinew(k))
          tt1 = 0.5d0*(tiold(k-1)+tiold(k))
          tt2 = 0.5d0*(tinew(k-1)+tinew(k))

          ttt = tme(k)
          capa = func_cp(ttt,tt1,tt2)
          m1=heold(k)*re(k)/3.d0*capa
          m2=heold(k)*re(k)/6.d0*capa
          k1=dt/heold(k)*kice(k)
          rhs(k-1) =  rhs(k-1) + k1 *( tiold(k)-tiold(k-1)) +
     & 0.5d0 * rad(k) * dt
          rhs(k)   =  rhs(k)   - k1 *( tiold(k)-tiold(k-1)) +
     & 0.5d0 * rad(k) * dt
          
          w1=wmesh(k-1)*re(k)*dt
          w2=wmesh(k  )*re(k)*dt
          rhs(k-1) =  rhs(k-1)  - w1 * em(k)
          rhs(k  ) =  rhs(k  ) + w2 * em(k)
          
          w1=w1 * capa * 0.5d0
          w2=w2 * capa * 0.5d0
      !    rhs(k-1) =  rhs(k-1) - w1 * ( tinew(k) + tinew(k-1) )
      !    rhs(k  ) =  rhs(k  ) + w2 * ( tinew(k) + tinew(k-1) )

          a  (k-1) =  a  (k-1) - k1 + m2 + w1
          diagon(k-1)= diagon(k-1) + k1 + m1 + w1
          diagon(k  )= diagon(k  ) + k1 + m1 - w2
          b  (k)   =  b  (k  ) - k1 + m2 - w2
      enddo !1,lice_top
     
      
      !---------------------------------------------------------------------
      ! implicit terms
      
      ! implicit terms  at surface
      diagon(lice_top) = diagon(lice_top) + theta_ther * dt * dzf
      !write(*,*) "theta_ther * dt * dzf", dzf

      !---------------------------------------------------------------------
      
      !---------------------------------------------------------------------
      ! add BCs:
      rhs(0       )= rhs(0       ) - dt * bshf
      !write(*,*) "bshf", bshf
      rhs(lice_top)= rhs(lice_top) + dt * Fnet - energy_snow_melt
      !write(*,*) "Fnet", Fnet ,"energy_snow_melt", energy_snow_melt
      flux2= 0.d0

      k=1
        tinew(k-1) = tseafrz ! keep basal temperature at sea temperature
        sum0 = em(k) * re(k) * dt
        w1=wmesh(k-1)*re(k)*dt
        rhs(k-1) =  rhs(k-1)  + w1 * em(k) ! add back
        rhs(k)   = rhs(k) - b(k) * ( tinew(k-1) - tiold(k-1) )
        b(k)= 0.d0
      k=0
        rhs(k) =  rhs(k) - diagon(k) * ( tinew(k) - tiold(k) )
        diagon(k)= sum0
      
      if ( melt_ts ) then
        k=lice_top
        w2=wmesh(k  )*re(k)*dtice
        rhs(k  ) =  rhs(k  ) - w2 * em(k) ! add back

        tinew(k) = tme(k) ! keep top surface temperature at melt point
        sum0 = em(k) * re(k) * dt
        rhs(k) =  rhs(k) - diagon(k) * ( tinew(k) - tiold(k) )
        diagon(k)= - sum0

        rhs(k-1)   = rhs(k-1) - a(k-1) * ( tinew(k) - tiold(k) )
        a(k-1)= 0.d0
      endif
      !---------------------------------------------------------------------
      ! add BCs for wmesh:
      
      !  efus = cp_ice * ( tme(lice_top) - tiold(lice_top) ) + mlfus * ( 1.d0 - tme(lice_top)/tiold(lice_top) ) - cp_wat * tme(lice_top)
      ! precipitation energy transfer
      if (lice_top == ns ) then
         flux2= flux2 + Fprec
      endif
      ! sublimation energy transfer
      if (.not.melt_ts) then
         efus  = em(lice_top)
         flux2 = flux2 + sublimation_speed * re(lice_top) * efus
      endif

      ! final combination
      rhs (lice_top) = rhs (lice_top) + dt * flux2
      
      !---------------------------------------------------------------------
      ! tri-diagonal solver call
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        do k=1,lice_top
            !write(*,*) rhs(k), "rhs(", k, ")"
            !write(*,*) a(k), "a(", k, ") "
            !write(*,*) b(k), "b(",k,")"
            !write(*,*) diagon(k),"d(", k, ")" 
        enddo !1,lice_top
      end if 

      call trisolverD(lice_top+1,a(0),b(0),diagon(0),rhs(0))
      
      ! new temperature profile
      tinew(0:lice_top)=rhs(0:lice_top)+tiold(0:lice_top)

      
      if (thin_snow_active) then
          tinew(lice_top+1:ns)=tiold(lice_top+1:ns)
      endif
        if (melt_ts ) then
          k=lice_top
          wmesh(k) = rhs(k)
          tinew(k) = tme(k)
          rhs  (k) = tinew(k) - tiold(k) ! needed for posterio diagnostic
        endif

      ! fixed bottom temperature at sea water temperature
      k=0
      wmesh(k) = rhs(k)
      tinew(k) = tseafrz
      rhs  (k) = tinew(k) - tiold(k) ! needed for posterio diagnostic
      if (i .eq. 10 .and. j .eq. 46) then
        !write(*,*) "wmesh(k)", wmesh(k)
      end if
      
      !---------------------------------------------------------------------
      ! Calculation of ice and snow mass changes due to melting
      !---------------------------------------------------------------------

      ! top
      k=lice_top
      slope = tinew(k)-tme(k)
      ! FD debug
      !write(*,*) 'ice',wmesh(ni),tinew(ni)

      if (slope >0.d0) then    ! condifion on top melting  (first time)
          melt_ts=.true.
          tinew(k) = tme(k)
      endif
      
      !---------------------------------------------------------------------
        ! reset top melting to 'off' if velocity has the wrong sign
        !---------------------------------------------------------------------

        if ( melt_ts ) then
            !write(*,*) "melt_ts", m,i,j
          if ( lice_top==ns ) then
             if ( wmesh(lice_top) < -snow_precip ) then
!               melt_ts = .false.
!               wmesh(lice_top)=-sublimation_speed-snow_precip
!               write(*,*) "wmesh(lice_top) < -snow_precip", 
!     &wmesh(lice_top), snow_precip,m,i,j
             endif
          else
             if ( wmesh(lice_top) < 0.d0 ) then
!               melt_ts = .false.
!               wmesh(lice_top)=-sublimation_speed
!               write(*,*) "wmesh(lice_top) < 0", wmesh(lice_top),m,i,j
             endif
          endif
        endif
        
        
      do k=0,lice_top
        if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
            !write(*,*) tinew(k),"tiNEW", k,i,j,m 
        end if
        if (tinew(k) .gt. 0.) then
        !write(*,*) tinew(k),"!!!!!!!!!!!!!!!!!!!", k,i,j,m 
            please_stop = .true.
         end if
      enddo !1,lice_top
      do k=0,lice_top
        if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
            !write(*,*) tiold(k),"tiOLD", k,i,j,m 
        end if
      enddo !1,lice_top
      !---------------------------------------------------------------------
      ! need to get rid of snow if wmesh too large
      !---------------------------------------------------------------------

      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        !write(*,*) "melt_ts", melt_ts, thin_snow_active
        end if
        if ( wmesh(ns)*dt > hsold + 1d-18 ) then
          sume2 = 0.d0
          do k=ni+1,ns
            sume2 = sume2 + em(k) * re(k) * heold(k)
          enddo
          energy_snow_melt = -sume2
          fwf = fwf + hsold / dt / rhowat * rhosno
          dhs = - hsold
          lice_top=ni ! do not treat the snow layer
          melt_ts = .false.
          thin_snow_active = .true.
          wmesh(ns)=hsold/dt
          sumrad=0.d0
          do k=1,lice_top
         !    sumrad = sumrad + rad(k) !!!!!!!!!!!!!!!1
          enddo
       endif
   
       !---------------------------------------------------------------------
      ! compute vertical velocity in sigma coordinate
      ! and new thickness for ice and snow
      !---------------------------------------------------------------------
      ! vertical regridding due to the sublimation/condensation at the top of the ice/snow
      ! vertical regridding due to the ice formation at the bottom of the ice

      ! apply a linear scaling between bottom and top ice growth/decay
        fac=0.d0
        do k=1,ni-1
          fac = fac + dzi(k)
          wmesh(k)=wmesh(0 ) * (1.d0-fac) + wmesh(ni) * fac
        enddo
        fac=0.d0
        do k=ni+1,ns-1
          fac = fac + dzi(k)
          wmesh(k)=wmesh(ni) * (1.d0-fac) + wmesh(ns) * fac
        enddo

      ! find final thickness after regridding

        dhi = ( wmesh(0 ) - wmesh(ni ) ) * dt
        dhs = ( 0.d0      - wmesh(ns ) ) * dt

        if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
          !write(*,*) "wmesh", wmesh(ns ), wmesh(ni ), wmesh(0 )
          !write(*,*) "dhi", dhi
          !write(*,*) "dhs", dhs
        end if
        if (-dhs > hsold) then
            !write(*,*) "-dhs > hsold", dhs, hsold
            hsnew=0
        end if
        if (-dhi > hiold) then
        !write(*,*) "-dhi > hiold", dhi, hiold, wmesh(0 )
            melted = .true.
            hsnew=0
            hinew=0
            exit
        end if
        hsnew = max(hsold+dhs,0.d0)
        hinew = hiold+dhi
        if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
            !write(*,*) "hinew", hinew 
            !write(*,*) "hiold", hiold 
        end if
      ! FD debug
      !write(*,*) 'snow',hsnew,dhs,wmesh(ns)*dtice

        he(   1:ni) = hinew*dzi(   1:ni)
        he(ni+1:ns) = hsnew*dzi(ni+1:ns)
        
        !---------------------------------------------------------------------
      ! reset ice/snow temperature if needed

        do k=lice_top,1,-1
         tinew(k  ) = min(tinew(k  ),tme(k))
         tinew(k-1) = min(tinew(k-1),tme(k))
       enddo

      do k=1,lice_top
        if (i .eq. 10 .and. j .eq. 46) then
            !write(*,*) tinew(k),"tttNEW", k,i,j,m 
        end if
      enddo
        maxvel=0.d0
        do k=0,lice_top
          maxvel=max(maxvel,abs(tinew(k)-timid(k)))
        enddo

        if (maxvel.lt.1d-12) then
        exit
        endif

       timid  = tinew
       wold   = wmesh
      
      enddo !1, maxiter
      
      if (i .eq. 14 .and. j .eq. 1 .and. m .eq. 13) then
        !write(*,*) "ITERATIONS FINISHED IN: ", iteration
        !write(*,*) "hiold: ", hiold
        
      end if
      
      !---------------------------------------------------------------------
      ! freshwater and heat flux on melting to the ocean
      !---------------------------------------------------------------------
        dq    = 0.0d0
        hmelt = 0.0d0

    !   if ( melt_ts ) then
    !  k=lice_top
    !     dhi   =   wmesh(k  ) * re(k) / rhowat
    !     hmelt = hmelt + dhi
    !     dq    = dq + dhi *cheat*(sstz-tme(k)) 
    !     fwf  = fwf  + hmelt   ! positive flux for input to ocean
    !     heat_melt_to_ocean = heat_melt_to_ocean  - dq
    !    endif

      ! surface flux
      topflux = Fnet - dzf * rhs(lice_top) - energy_snow_melt / dt
      
      !---------------------------------------------------------------------
      ! reset top temperature for thin snow for uniformity
      if ( thin_snow_active) then
        do k=lice_top+1,ns
          tinew(k)=tinew(k-1)
        enddo
      endif
      
      !---------------------------------------------------------------------
      ! preparation for ocean thermodynamics

        bshf = bshf +  heat_melt_to_ocean
        
        hi=hinew
        hs=hsnew
       ! write(*,*) hi, "hi"
        sim=0.D0
        do k=1,lice_top
          sim=sali(k)*he(k)
        enddo
        if (hi>0.d0) sim=sim/hi
        
        !
      ! When snow load excesses Archimede's limit, snow-ice interface goes down
      ! under sea-level, flooding of seawater transforms snow into ice
      ! dh_snowice is positive for the ice

      dh_sni = MAX( 0.d0 , ( rhosno * hs + (rhoice - rhowat ) * hi)
     &  / ( rhosno + rhowat - rhoice ) )/100.

      if (dh_sni > 0.d0 ) then
       hi  = hi + dh_sni
       hs  = hs - dh_sni
      endif
      
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
        write(*,*) "hi", hi, m, i, j
        write(*,*) "hs", hs, m, i, j
        write(*,*) "dh_sni", dh_sni
      end if 

      dHice = hi*Aice(m,i,j)*100 - Hice (m,i,j)
      dHsnow = hs*Aice(m,i,j)*100 - Hsnow(m,i,j)

      Hice (m,i,j) = hi*Aice(m,i,j)*100
      Hsnow(m,i,j) = hs*Aice(m,i,j)*100
      
      dHiceT = dHiceT + dHice
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then 
        write(*,*) "dHice", dHice
        write(*,*) "dHsnow", dHsnow
      end if
      if (dHice .ne. dHice .or. dHice .gt. 1000 
     & .or. dHice .lt. -1000) then
        write(*,*) "mmmm", m, i, j
!!!       	Hsnow(m,i,j)=0.
!!!	do k=ni+1,ns
!!!          TsnowFE(m,i,j,k-ni)=-10.  !!!!!!!!!!!!!!tme(k)
!!!      enddo
!!!	Hice(m,i,j) =0.
!!!	Aice(m,i,j) =0.
!!!	do k=1,ni
!!!          TiceFE(m,i,j,k)= tme(k)
!!!      enddo
        stop
      end if
      
      do k=0,lice_top
        if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
            write(*,*) tinew(k),"tiNEW", k,i,j,m 
        end if
      enddo !1,lice_top
      
        if (i .eq. 34 .and. j .eq. 43) then
            write(*,*) "Aice", Aice(m,i,j), i,j,m 
        end if
      dHsnowT = dHsnowT + dHsnow
      
      do k=1,ni
          TiceFE(m,i,j,k)= tinew(k)
        enddo
	if ( lice_top==ns ) then   
        do k=ni+1,ns
          TsnowFE(m,i,j,k-ni)= tinew(k)
        enddo
	end if
	
	Tice(m,i,j) = TiceFE(m,i,j,1)
	Tsnow(m,i,j) = TsnowFE(m,i,j,1)
	
	if(Aice(m,i,j).GT.Aimin) then

cc	ppp= 0.5e-5*g*row*Hi                  
      TFC= TFr(S(i,j,1), 0.) 
      TW =T(i,j,1)
      
      IF( TW .GT. TFC) then
 
***	Side= 3.14*1.6e-6*SQRT((TW-TFC)**3)*DT/(0.66*300.) !See NCAR manual
***	Side= 4.*1.6e-6*SQRT((TW-TFC)**3)*DT/300.          !My estimate
	Side= 1.44e-8*SQRT((TW-TFC)**3)*DT/100.              !Shmidt et.al., 2003
      dAice=  -side*Aice (m,i,j)
      dHice=  -side*Hice (m,i,j)
      dHsnow= -side*Hsnow(m,i,j)
      
      if (i .eq. 34 .and. j .eq. 43  .and. m .eq. 10) then
!            write(*,*) "dAice",dAice,i,j,m 
!            write(*,*) "dHice",dHice,i,j,m 
!            write(*,*) "dHsnow",dHsnow,i,j,m 
        end if

!      Hsnow(m,i,j)= Hsnow(m,i,j) +dHsnow
!      Hice (m,i,j)= Hice (m,i,j) +dHice
!      Aice (m,i,j)= Aice (m,i,j) +dAice
!
!      dHsnowT = dHsnowT +dHsnow
!      dHiceT  = dHiceT  +dHice

!      Q_melt= Qsnow*dHsnow +Qice*dHice
!      T(i,j,1)= T(i,j,1) +2.*Q_melt/(hz(1)*row*Cpw)
      
      end if  ! T > Tf

      else  
   
      Hice (m,i,j) = 0.
      Aice (m,i,j) = 0.
	Tice (m,i,j) = 0.
      Hsnow(m,i,j) = 0.
      Tsnow(m,i,j) = 0.
	do k=ni+1,ns
          TsnowFE(m,i,j,k-ni)=-10.  !!!!!!!!!!!!!!tme(k)
      enddo
	do k=1,ni
          TiceFE(m,i,j,k)= tme(k)
      enddo
      end if ! Ice Lateral Melting.
      
      !end if !hi.gt.hminice
	end if  ! ai(m) > 0
      enddo
      
	do m=1,mgrad

c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, lateral melting', Hsnow(m,i,j)
c	end if


      if(Hsnow(m,i,j)/100 .LT.Hsmin) then
      dHsnowT = dHsnowT +Hsnow(m,i,j)
	Hsnow(m,i,j)=0.
	do k=ni+1,ns
          TsnowFE(m,i,j,k-ni)=-10. !!!!!!!!!!!!tme(k)
      enddo
      Tsnow(m,i,j) = 0.
	end if
      if(Hice(m,i,j)/100. .LT.Himin) then
      dHsnowT = dHsnowT +Hsnow(m,i,j)
      dHiceT  = dHiceT  +Hice(m,i,j)
      !write(*,*) "Hice(m,i,j).LT.Himin", Hice(m,i,j), m, i, j
	Hsnow(m,i,j)=0.
	do k=ni+1,ns
          TsnowFE(m,i,j,k-ni)=-10.  !!!!!!!!!!!!!!tme(k)
      enddo
!      Qdelta= Qice*(-Hice(m,i,j)) -Qsnow*Hsnow(m,i,j)
!        T(i,j,1)= T(i,j,1) + 2.*Qdelta/(hz(1)*row*Cpw)
	Hice(m,i,j) =0.
	Aice(m,i,j) =0.
		Tice (m,i,j) = 0.
      Tsnow(m,i,j) = 0.
	do k=1,ni
          TiceFE(m,i,j,k)= tme(k)
      enddo
	end if

c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*)'Snow, after flooding',Hsnow(m,i,j),hflood,hi,hdraft
c	end if


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
     * 2*          !!!!!!!!!!!!
     &1.087*Cpw*row*vol*(TFC-T(i,j,k))/Qice/cg(abs(nt3(i,j,1)))
 
         
      if (dHiceN .gt. 1) then
       ! write(*,*) "TFC", TFC, T(i,j,k), i, j, k
      end if
*     New ocean temperature    
      T(i,j,k)= TFC
      
*     Correction of salinity

*      delta_S= S(i,j,1)-Sice
*      S(i,j,k) = S(i,j,k) +dHiceN*(S(i,j,k)-Sice)*roi/(row*vol)
*     Compensation      
*      S(i,j,1) = S(i,j,1)-2.*dHiceN*delta_S*roi/(row*hz(1))
      
      end if  ! T<TFC
      end do  ! k=1,kb loop

      if (aice(0,i,j) .gt. 0.9) then
        !write(*,*) "aice(0,i,j)", aice(0,i,j)
      end if
	IF(dHiceN .GT. 0.) then
!
!      if (dHiceN .GT. 10) then
!        dHiceN = 10
!      end if

c     In a case of sufficient open water: frazil ice to the open water.
	hhh0= dHiceN/MAX(aimin,aice(0,i,j))

c     Thickness of the new ice

      IF( hhh0 .LE. Hmax(2)) then
	delta=dHiceN/Href

      Hice(1,i,j)=Hice(1,i,j)+dHiceN
      aice(1,i,j)= aice(1,i,j)+ MIN(delta, aice(0,i,j))
      aice(0,i,j)= MAX(0.,aice(0,i,j)- delta)

c     New ice temperature - for T in situ only!
   !!!!   Q= (Tice(1,i,j)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
      Q= (TiceFE(1,i,j,1)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
      do k=1,ni
!        if (TiceFE(1,i,j,k) .eq. 0.) then
!            TiceFE(1,i,j,k)= (min(Q/Hice(1,i,j), -2.))
!            end if
          TiceFE(1,i,j,k)=  tinew(k)
      enddo
      Tice(1,i,j)=tinew(1)
  !!!    Tice(1,i,j)=Q/Hice(1,i,j)
      Tsnow(1,i,j)=-10.
      do k=ni+1,ns
            TsnowFE(1,i,j,k)=-10.
      enddo
      else

c     New Ice mass spreaded unifirmly over entire cell
c     if no sufficient open water
      Hice(1,i,j)=Hice(1,i,j)+dHiceN*(Aice(0,i,j)+Aice(1,i,j))
      do m=2,mgrad
	hi= Hice(m,i,j)/max(aimin,Aice(m,i,j))
	Hice(m,i,j)=Hice(m,i,j)+dHiceN*Aice(m,i,j)
	end do

c     New ice temperature - for T in situ only! Should be corrected
    !!!!   Q= (Tice(1,i,j)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
      Q= (TiceFE(1,i,j,1)+TFC)*Hice(1,i,j)/2. +TFC*dHiceN*Aice(1,i,j)
      do k=1,ni
!        if (TiceFE(1,i,j,k) .eq. 0.) then
!            TiceFE(1,i,j,k)= (min(Q/Hice(1,i,j), -2.))
!            end if
          TiceFE(1,i,j,k)=  tinew(k)
          if (i .eq. 26 .and. j .eq. 10 .and. m .eq. 1) then 
       ! write(*,*) "TiceFE(1,i,j,k)", TiceFE(1,i,j,k)
      end if
      enddo
      Tice(1,i,j)=tinew(1)
  !!!    Tice(1,i,j)=Q/Hice(1,i,j)
      Tsnow(1,i,j)=-10.
      do k=ni+1,ns
            TsnowFE(1,i,j,k)=-10.
      enddo


      aice(1,i,j)= aice(1,i,j)+ aice(0,i,j)
      aice(0,i,j)= 0.

	end if

      dHiceT= dHiceT + dHiceN
      !write(*,*) "!!!!!!!!!!!!!!!dHiceT", dHiceT, "dHiceN", dHiceN

	end if !dHiceN .GT. 0.

c	if(m.eq.10.and.i.eq.17 .and. j.eq.24) then
c	write(*,*) 'Snow, return', Hsnow(m,i,j)
c	end if
      if (please_stop .and. .not. melted) then
        !stop
      end if

      if (i .eq. 34 .and. j .eq. 43) then
      aavg=0.
            do m=1,mgrad
            aavg = Aice(m,i,j) + aavg
            end do
            write(*,*) "aavg",aavg,i,j,m 
        end if
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
ccc      IF(ABS(Denom) .LT. 1.e-19) denom= C(j1)
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
      ai= 0.60  ! AOMIP, approx. 3m.
cc      ai= 0.70  ! New AOMIP, 15.05.2003 approx. 4m.
cc      ai= 0.65  ! Boulder Ice CCSM3 (Ice version 4, 2002) 
cc      ai= 0.575  ! LANL CICE 4.1 

	if(h.LE.50.)then  ! see BoulderIce or CICE 4.1.
	F_ai= aow +(ai-aow)*h/50.
	else
	F_ai= ai
	end if

	else

cc      ai= 0.70  ! New AOMIP, 15.05.2003 approx. 4m.
cc      ai= 0.65  ! Boulder Ice CCSM3 (Ice version 4, 2002)
      ai= 0.60  ! AOMIP, approx. 3m.
cc      ai= 0.575  ! LANL CICE 4.1

	if(h.LE.50.)then  ! see BoulderIce - this is a good approx.
cc	F_ai= aow +(ai-0.075*(T+1.)-aow)*h/50. 
	F_ai= aow +(ai-0.1  *(T+1.)-aow)*h/50. 
	else
cc	F_ai= ai -0.075*(T+1.)
	F_ai= ai -0.1  *(T+1.)
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
cc	F_as=0.80  ! Old Aomip, close to CCSM3
cc	F_as=0.82  ! SHEBA
	F_as=0.85  ! Pinto, 1999, SHEBA, LANL CICE 4.1
cc	F_as=0.81  ! New AOMIP, 15.05.2003
	else
cc	F_as=0.65  ! F&B 1986
cc	F_as=0.70  ! Old AOMIP SHEBA 2001 -SHCE
cc	F_as=0.72  ! Boulder
cc	F_as=0.80 -0.125*(T+1.0)  ! AOMIP modified
cc	F_as=0.85 -0.125*(T+1.0)  ! SHEBA + LANL CICE 4.1
cc	F_as=0.80 -0.1*(T+1.0)  ! SHEBA + LANL CICE 4.1
	F_as=0.85 -0.1*(T+1.0)  ! Pinto,1999, SHEBA+Boulder Ice, CCSM3
cc	F_as=0.77  ! New AOMIP, 15.05.2003
	end if

	return
	end 

	Real Function CDH(tsu,TA)
*     version 17.12.2009. Regularisation for sections method.
	!if(tsu.GT.TA+1.0) then
	CDH=1.75e-3       ! Unstable boundary layer, Makshtas, SP-23.
	!else
	!  if(tsu.LT.TA-1.0) then
      !  CDH=1.2e-3      ! Stable boundary layer, AOMIP
      !  else
      !  CDH=1.2e-3 +0.5*(tsu-TA+1.0)*0.55e-3 ! Transition case
	!  end if
	!end if
      return
	end
	
	Real Function CDL(tsu,TA)
*     version 17.12.2009. Regularisation for sections method.
	!if(tsu.GT.TA+1.0) then
	CDL=1.75e-3       ! Unstable boundary layer, Makshtas, SP-23.
	!else      
	!  if(tsu.LT.TA-1.0) then
      !  CDL=1.5e-3      ! Stable boundary layer, AOMIP
      !  else
      !  CDL=1.5e-3 +0.5*(tsu-TA+1.0)*0.25e-3 ! Transition case
	!  end if
      !end if
      return
	end
	
	real function Tfreeze1(ssss)
	
	real(8) ssss
	real(8) :: a = -0.054d0
        Tfreeze1=ssss*a
      return
      end

      real FUNCTION func_cp(ttt,tt1,tt2)
      
      real(8) ttt, tt1, tt2
      real(8) :: cp_ice = 2.062e+03_8
      real(8) :: mlfus     = 3.335e+05_8
      real(8) :: TT

      TT = MIN (tt1 , -1d-10 ) * MIN ( tt2, -1d-10 )
      
      func_cp = MIN( cp_ice - mlfus * ttt / TT, 1d9 )
      return
      END
      
      FUNCTION func_qmelt(ttt,tt1)

      real(8) tt1, ttt
      real(8) :: cp_ice = 2.062e+03_8
      real(8) :: mlfus     = 3.335e+05_8

      func_qmelt =cp_ice*(ttt - tt1)+mlfus*(1.d0-ttt/MIN( tt1, -1d-20 ))

      return
      END
      
      
      FUNCTION func_El(ttt,tt1)

      real(8) ttt, tt1
      real(8) :: cp_ice = 2.062e+03_8
      real(8) :: mlfus     = 3.335e+05_8
      real(8) :: cp_wat    = 3.99e+03_8
      
      func_El =cp_ice*(tt1-ttt)-mlfus*(1.d0-ttt /MIN(tt1,-1d-20)) +
     &cp_wat*ttt

      return
      END
      
      FUNCTION func_bf(Tf,T)

      implicit none
      real(8)  func_bf, Tf, T

      func_bf = Tf / MIN ( T, -1d-20 )

      END FUNCTION func_bf
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      SUBROUTINE trisolverD(la,a,b,diag,rhs)
!--------------------------------------------------------------
! TRIAGONAL SYSTEM SOLVER
! for double float
!--------------------------------------------------------------
!
      IMPLICIT NONE
      integer la
      real(8) A(la), B(la),DIAG(la),RHS(la)
      real(8) tmp,tmpdiag
      INTEGER K
! first Gauss elimination
!
! At surface

          tmp      = 1.d0 / DIAG(1)
          a(1)     = A(1)   * tmp
          rhs(1)   = rhs(1) * tmp 
!
! at mid depth
!
          DO K = 2, LA-1
            tmpdiag  = DIAG(K) - A(K-1) * B(K)
            tmp      = 1.d0 / tmpdiag
            a(K)   = A(K) * tmp
            rhs(K) = (rhs(K) -rhs(K-1)*B(K) ) * tmp
          END DO
!
! at bottom, solution in rhs
!
          tmpdiag    = DIAG(LA) - A(LA-1) * B(LA)
          tmp        = 1.d0 / tmpdiag
          rhs(LA)    = (rhs(LA) -rhs(LA-1)*B(LA) ) * tmp

!--------------------------------------------------------------
! second and final Gauss elimination to surface
! solution in rhs

          DO K = LA-1,1,-1
            rhs(K)   = rhs(K) -rhs(K+1)*A(K)
          END DO

      RETURN
      END SUBROUTINE trisolverD

!---------------------------------------------------------------------


      FUNCTION func_ki(ssss,tt)

      real(8)  func_ki, ssss, tt, betanew
      real(8) :: ki0 = 2.034d+0        ! thermal conductivity of fresh ice	[W/m/C]
      real(8) :: mu = 0.054d0          ! empirical constant relating S and Tf	[C/psu]

      betanew = (ki0-0.02d0)*mu                 ! see HEADER
      func_ki = ki0 + betanew * ssss / MIN( tt, -1d-20 )

      END FUNCTION func_ki