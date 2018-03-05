      subroutine CH4Rivers
     *          (month,CH4,nt3,km2,si,cg,il1,jl1,kl,r,hx,hy,hz,dt)

*     version 14.10.2013.
*     AOMIP River runoff, Prange, cm3/s.
*     Monthly mean climatology.

      dimension CH4(0:il1,0:jl1,kl)
	dimension nt3(0:il1,0:jl1,kl), km2(0:il1,0:jl1), cg(13)
	dimension si(0:jl1), hz(kl)
	dimension TW(12)                 ! Month time weight
	data TW/4*0.02,0.1,0.31,0.17,0.12,0.1,0.07,0.03,0.02/

	c6= 1./6.
	RR= 1.e-6*R*R               ! Normalization
      C= 12.*TW(month)/(RR*hx*hy) ! year mean by planar area element

c     Ob and Yenisey, Puhr, Taymyra, Pyasina

**    CH4_0 = 50.0   ! Methane concentration in river Golubeva
      CH4_0 = 150.0  ! Book Lein Tabl. 3.5.1

      i=34
	Volume= 0.
	do j=27,29
	do k=1,km2(i,j)

	if(k.eq.1) then
	hzk= hz(k)
	hzk1=0.
	else

      	if(k.lt.km2(i,j)) then
	    hzk= hz(k  )
	    hzk1=hz(k-1)
	    end if

	end if

	if(k.eq.km2(i,j)) then
	hzk= 0.
	hzk1=hz(k-1)
	end if

	Volume= Volume+ 0.5*(hzk+hzk1)*c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	end do

	do j=27,29
	do k=1,km2(i,j)
      CH4(i,j,k)= CH4(i,j,k) + CH4_0*45918.9410*C*dt/Volume
	end do
	end do

c     Lena, Khatanga, Olenek, Yana

      CH4_0 = 300.0   ! Methane concentration in river

      i= 31
	Volume= 0.
	do j=11,14
      do k=1,km2(i,j)

	if(k.eq.1) then
	hzk= hz(k)
	hzk1=0.
	else

      	if(k.lt.km2(i,j)) then
	    hzk= hz(k  )
	    hzk1=hz(k-1)
	    end if

	end if

	if(k.eq.km2(i,j)) then
	hzk= 0.
	hzk1=hz(k-1)
	end if


	Volume= Volume+ 0.5*(hzk+hzk1)*c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	end do

	do j=11,14
      do k=1,km2(i,j)
      CH4(i,j,k)= CH4(i,j,k) +CH4_0*26964.1602*C*dt/Volume
	end do
	end do

c     Mackenzie

c     De Angelis MA, Lilley MD (1987) Methane in surface waters of Oregon estuaries
c     and rivers. Limnol Oceanogr 32:716-722

      CH4_0 = 5.0   ! Methane concentration in river

      i=2
	j=10

      do k=1,km2(i,j)
	if(k.eq.1) then
	hzk= hz(k)
	hzk1=0.
	else

      	if(k.lt.km2(i,j)) then
	    hzk= hz(k  )
	    hzk1=hz(k-1)
	    end if

	end if

	if(k.eq.km2(i,j)) then
	hzk= 0.
	hzk1=hz(k-1)
	end if

	Volume= 0.5*(hzk+hzk1)*c6*cg(abs(nt3(i,j,1)))*si(j)
	end do

      do k=1,km2(i,j)
      CH4(i,j,k)= CH4(i,j,k) +CH4_0*11838.9104*C/Volume
	end do

c     Dvina and Mezen

*     Russian book, Tabl. 3.5.1.

      CH4_0 = 201.3   ! Methane concentration in river 

      i=34
	j=43

      do k=1,km2(i,j)
	if(k.eq.1) then
	hzk= hz(k)
	hzk1=0.
	else

      	if(k.lt.km2(i,j)) then
	    hzk= hz(k  )
	    hzk1=hz(k-1)
	    end if

	end if

	if(k.eq.km2(i,j)) then
	hzk= 0.
	hzk1=hz(k-1)
	end if

	Volume= 0.5*(hzk+hzk1)*c6*cg(abs(nt3(i,j,1)))*si(j)
	end do

      do k=1,km2(i,j)
      CH4(i,j,k)= CH4(i,j,k) +CH4_0*4328.8674*C/Volume
	end do

c     Pechora

      CH4_0 = 0.0   ! No data

      i=35
	j=36

      do k=1,km2(i,j)
	if(k.eq.1) then
	hzk= hz(k)
	hzk1=0.
	else

      	if(k.lt.km2(i,j)) then
	    hzk= hz(k  )
	    hzk1=hz(k-1)
	    end if

	end if

	if(k.eq.km2(i,j)) then
	hzk= 0.
	hzk1=hz(k-1)
	end if

	Volume= 0.5*(hzk+hzk1)*c6*cg(abs(nt3(i,j,1)))*si(j)
	end do

      do k=1,km2(i,j)
      CH4(i,j,k)= CH4(i,j,k) +CH4_0*6072.7550*C/Volume
	end do

c     Kolyma+Indigirka

      CH4_0 = 200.00   ! Methane concentration in river

      i=24
	j=6

      do k=1,km2(i,j)
	if(k.eq.1) then
	hzk= hz(k)
	hzk1=0.
	else

      	if(k.lt.km2(i,j)) then
	    hzk= hz(k  )
	    hzk1=hz(k-1)
	    end if

	end if

	if(k.eq.km2(i,j)) then
	hzk= 0.
	hzk1=hz(k-1)
	end if

	Volume= 0.5*(hzk+hzk1)*c6*cg(abs(nt3(i,j,1)))*si(j)
	end do

      do k=1,km2(i,j)
      CH4(i,j,k)= CH4(i,j,k) +CH4_0*6110.2248*C/Volume
	end do

      return
	end

	subroutine VertCH4

*     25.03.2015

*     Parameters for Bunsen solubility coefficients,
*     Detaile description in Wanninkhof, 1992 
	parameter(A1=-68.8862,A2=101.4956,A3=28.7314)
	parameter(B1=-0.076146,B2=0.043970,B3=-0.0068672)


      INCLUDE 'Slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
      INCLUDE 'Tparm.fi'


*     Shelf/slope methane production.
*     Black Sea sediment production (100-1500 m)
*     WILLIAM S. REEBURGH, et. al., Black Sea methane geochemistry.
*     Deep-Sea Research.Vol. 38. Suppl. 2. pp. SI189-S1210. 1991

cc      Flux_From_Bottom= 4.8516  ! nM cm/s =1.53 mol/m2/y=4.85e-8mole/m2/s 

*     Global estimate 3.e-9 mole/m2/s for 100-1000 m
*     Judd, A. G. (2003), The global importance and context of methane escape
*     from the seabed, Geo Mar. Lett., 23, 147–154, doi:10.1007/s00367-003-
*     0136-z.

*     Used in:
*     Elliott, S., M. Maltrud, M. Reagan, G. Moridis, and P. Cameron-Smith (2011), 
*     Marine methane cycle simulations for the period of early global warming, 
*     J. Geophys. Res., 116, G01010, doi:10.1029/2010JG001300.
*     Sensitivity studies showed that it should be 3-30 times less !!!
*     Happens at 100-1000m

*     Savvichev 2.6-61mcMol/m2/day = 3.e-3 - 7.e-2 nM cm/s

      Flux_From_Bottom= 3.e-1  !!!! nM cm/s = 3.e-9 mole/m2/s

      Flux_From_Bottom= 0.1*Flux_From_Bottom

*     Possible clathrate-driven input at some specified locations
*     1.e-6 mole/m2/s = 100. nM cm/s
*     Reagan, M. T., and G. J. Moridis (2008), Dynamic response of oceanic
*     hydrate deposits to ocean temperature change,
*     J. Geophys. Res., 113, C12023, doi:10.1029/2008JC004938. 

*     Should be located only at some places like ESAS

      Flux_From_Bottom_MH= 0. !!!  100.  ! nM cm/s 


*     Fraction of bubbles - 0 - all methane is dissolved at bottom,
*     1 - all methane in bubbles immediately ascending to the sea surface.

      Bubble_fraction = 0.

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
*-----------------------------------------------------------------
*     Calculation in all points.

      IF( KB .GT. 0) THEN


*     Wind (wx,wy) is in m/s.
      WIND=wx(i,j)**2+wy(i,j)**2

*     Flux to atmosphere Wanninkhof, 1992 

      TTT = T(i,j,1)
	TK100  = (TTT +273.16)/100.

      Schmidt= 2039.2 +(-120.31 +(3.4209 -0.040437*TTT)*TTT )*TTT
	FTA= 0.31*WIND*SQRT(660./Schmidt)   !! [cm/h]
	FTA= FTA/3600.                      !! [cm/s]

*     Bunsen solubility

      beta= A1 +A2/TK100 +A3*ALOG(TK100) +  
     +      S(i,j,1) * (B1 +B2*TK100 + B3*TK100**2) 
	beta = EXP(beta)                    !! mL/mL/Atm    
	beta = 1.e9*beta/22.4               !! nmol/L/Atm    

	CH4_Saturated= beta*1.75e-6   ! Arctic CH4 concentration 1.7-1.8ppm


*     -------------- Ocean surface -----------------
      CA=2.*DT/hz(1)
      AM(1)= 0.
      CM(1)= -CA*AZS(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) + CA*PME(i,j)/dt

*     Flux through open water Aice(0,j,j)

      FM(1)= CH4(i,j,1)  +CA* FTA*CH4_Saturated *Aice(0,i,j) 
      BM(1)= BM(1)       +CA* FTA               *Aice(0,i,j)

*     Bubbles
      FM(1)= FM(1)  +CA*Bubble_fraction*Flux_From_Bottom_MH ! ESAS flux


*     -------------- Deep water -----------------------------------
      DO K=2, KB-1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

      AM(K)= -CA*CG( n)*AZS(i,j,k-1)/HZ(k-1)
      CM(K)= -CA*CG(np)*AZS(i,j,k  )/HZ(k)
      BM(K)= 1. -AM(K) -CM(K)
      FM(K)= CH4(i,j,k)

      END DO

*     ------------- Ocean bottom ----------------------------------
      CA=2.*DT
      Q=1./HZ(KB-1)

      AM(KB)= -CA*Q*Q*AZS(i,j,kb-1)
      BM(KB)= 1. -AM(KB)
      CM(KB)= 0.

cc	if(z(kb) .GE. 100.e2 .and. z(kb) .LE. 1000.e2) then ! Judd, A. G. (2003)
cc	if(z(kb) .LE. 1000.e2) then           ! Variations of Judd, A. G. (2003)
	if(z(kb) .LE. 100.e2) then                   ! Malakhova&Golubeva (2013)
      FM(KB)= CH4(i,j,kb)  +CA*Q*Flux_From_Bottom ! specified
	else
      FM(KB)= CH4(i,j,kb)  
	end if


	if( i.gt.17 .and. j.lt.24 .and. z(kb) .LE. 20.e2) then           
      FM(KB)= FM(KB)  +CA*Q*(1.-Bubble_fraction)*Flux_From_Bottom_MH ! ESAS flux
	end if

*     ------------- Matrix factorization --------------------------

      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)

      Do K=1,KB
      CH4(i,j,k)= MAX(0.,RKSI(K))
      End do

*     Flux to atmosphere nM cm/s and then in nM cm/s 

      FlCH4A= - FTA*(CH4_Saturated -CH4(i,j,1))*Aice(0,i,j)

      Flux_CH4_To_A(i,j)= FlCH4A* 13.824 ! mg[CH4]/m2/day

cc	write(*,*) i,j,FTA,WIND

      END IF
      END DO
      END DO


      RETURN
      END

	subroutine CH4Oxidation

*     Version 17.02.2014

      INCLUDE 'Slo2.fi'

*     Parameter CTime - is the characterist time of oxidation
*     by various mechanisms.

c     CTime = 3600.*24./0.015   ! Kara Sea, Lein and Ivanov, 2009
c	CTime = 1.e4 * 24.* 3600. ! 10 000 days in deep ocean
c	CTime = 1.e3 * 24.* 3600. ! 1  000 days in deep ocean

      do j=1,jl
	do i=1,il
	if(km2(i,j) .GT. 0) then

	do k=1,km2(i,j)    


      if(o2(i,j,k).GE.1.e3  .and. bio(i,j,k).GT.1.e-3) then ! O2 > 10 microM
c     Quadratic formula by Elliott, et.al., 2011
c     10 days for 1 microM [CH4]
      
	Ctime= 3.6*2.4e8/ max(1.e-3,bio(i,j,k)*CH4(i,j,k))    ! [sec]

	else                             ! No oxidation

      Ctime= 1.e19
	bio(i,j,k) =-1.e-3               ! No bacteria

	end if

	OldCH4 = CH4(i,j,k)

	CH4(i,j,k) = OldCH4 /(1.+dt/CTime)
	deltaCH4 = CH4(i,j,k)-OldCH4

*     Oxigen consuption for methane oxydation 
*     stoichiometry CH4:O2:carbonate -1:-2:+1

      o2(i,j,k) = max(0., o2(i,j,k) +2.*deltaCH4)

	end do

	end if

	end do
	end do

      return
	end

	subroutine CH4_Atlantic_Boundary

*     version 03.11.13

*     CH4 at saturated level

*     Parameters for Bunsen solubility coefficients,
*     Detaile description in Wanninkhof, 1992 
	parameter(A1=-68.8862,A2=101.4956,A3=28.7314)
	parameter(B1=-0.076146,B2=0.043970,B3=-0.0068672)


      INCLUDE 'Slo2.fi'

*     Norwegian Sea

      j = 49
      do i= 11,23
	do k=1,km2(i,j)    


      TTT = T(i,j,k)
	TK100  = (TTT +273.16)/100.

*     Bunsen solubility

      beta= A1 +A2/TK100 +A3*ALOG(TK100) +  
     +      S(i,j,1) * (B1 +B2*TK100 + B3*TK100**2) 
	beta = EXP(beta)                    !! mL/mL/Atm    
	beta = 1.e9*beta/22.4               !! nmol/L/Atm    

	CH4Obs(i,j,k)= beta*1.75e-6   ! Arctic CH4 concentration 1.7-1.8ppm

      end do
	end do

*     Denmark Strait
    
      i = 7
      do j=43,46
	do k=1,km2(i,j)  
	

      TTT = T(i,j,k)
	TK100  = (TTT +273.16)/100.

*     Bunsen solubility

      beta= A1 +A2/TK100 +A3*ALOG(TK100) +  
     +      S(i,j,1) * (B1 +B2*TK100 + B3*TK100**2) 
	beta = EXP(beta)                    !! mL/mL/Atm    
	beta = 1.e9*beta/22.4               !! nmol/L/Atm    

	CH4Obs(i,j,k)= beta*1.75e-6   ! Arctic CH4 concentration 1.7-1.8ppm

      end do
	end do
	
      return	  
      end

	subroutine O2_Boundary

*     version 15.11.13

*     O2 at saturated level

*     Parameters for Bunsen solubility coefficients,
*     Detaile description in Wanninkhof, 1992 
	parameter(A1=-58.3877,A2=85.8077,A3=23.8439)
	parameter(B1=-0.034892,B2=0.015568,B3=-0.0019387)


      INCLUDE 'Slo2.fi'

	do j=1,jl
	do i=1,il
	if(km2(i,j) .GT. 0) then
	do k=1,km2(i,j)

	if( nt3(i,j,k) .LT. 0) then

      TTT = T(i,j,k)
	TK100  = (TTT +273.16)/100.

*     Bunsen solubility

      beta= A1 +A2/TK100 +A3*ALOG(TK100) +  
     +      S(i,j,1) * (B1 +B2*TK100 + B3*TK100**2) 
	beta = EXP(beta)                    !! mL/mL/Atm    
	beta = 1.e9*beta/22.4               !! nmol/L/Atm    

	o2Obs(i,j,k)= beta*209.46e-3   ! O2 saturated concentration

	end if   ! nt <0

      end do   ! k

	end if   ! km2 >0

      end do   ! i
	end do   ! j

      return	  
      end

	subroutine VertO2

*     25.03.2015

*     Parameters for Bunsen solubility coefficients,
*     Detaile description in Wanninkhof, 1992 
	parameter(A1=-58.3877,A2=85.8077,A3=23.8439)
	parameter(B1=-0.034892,B2=0.015568,B3=-0.0019387)


      INCLUDE 'Slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
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
*-----------------------------------------------------------------
*     Calculation in all points.

      IF( KB .GT. 0) THEN


*     Wind (wx,wy) is in m/s.
      WIND=wx(i,j)**2+wy(i,j)**2

*     Flux to atmosphere Wanninkhof, 1992 

      TTT = T(i,j,1)
	TK100  = (TTT +273.16)/100.

      Schmidt= 2039.2 +(-120.31 +(3.4209 -0.040437*TTT)*TTT )*TTT
	FTA= 0.31*WIND*SQRT(660./Schmidt)   !! [cm/h]
	FTA= FTA/3600.                      !! [cm/s]

*     Bunsen solubility

      beta= A1 +A2/TK100 +A3*ALOG(TK100) +  
     +      S(i,j,1) * (B1 +B2*TK100 + B3*TK100**2) 
	beta = EXP(beta)                    !! mL/mL/Atm    
	beta = 1.e9*beta/22.4               !! nmol/L/Atm    

	O2_Saturated= beta*209.46e-3   ! O2 - 20,946% of atmosphere

cc	write(*,*) 'O2 saturated', i,j, O2_Saturated, o2(i,j,1)


*     -------------- Ocean surface -----------------
      CA=2.*DT/hz(1)
      AM(1)= 0.
      CM(1)= -CA*AZS(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) + CA*PME(i,j)/dt

*     Flux through open water Aice(0,j,j)

      FM(1)= o2(i,j,1)  +CA* FTA*O2_Saturated *Aice(0,i,j) 
      BM(1)= BM(1)      +CA* FTA              *Aice(0,i,j)

*     -------------- Deep water -----------------------------------
      DO K=2, KB-1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

      AM(K)= -CA*CG( n)*AZS(i,j,k-1)/HZ(k-1)
      CM(K)= -CA*CG(np)*AZS(i,j,k  )/HZ(k)
      BM(K)= 1. -AM(K) -CM(K)
      FM(K)= o2(i,j,k)

      END DO

*     ------------- Ocean bottom ----------------------------------
      CA=2.*DT
      Q=1./HZ(KB-1)

      AM(KB)= -CA*Q*Q*AZS(i,j,kb-1)
      BM(KB)= 1. -AM(KB)
      CM(KB)= 0.
      FM(KB)= o2(i,j,kb)  

*     ------------- Matrix factorization --------------------------

      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)

      Do K=1,KB
      o2(i,j,k)= MAX(0.,RKSI(K))
      End do

      END IF
      END DO
      END DO

      RETURN
      END

	subroutine VertBIO

*     25.03.2015

      INCLUDE 'Slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
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
*-----------------------------------------------------------------
*     Calculation in all points.

      IF( KB .GT. 0) THEN


*     -------------- Ocean surface -----------------
      CA=2.*DT/hz(1)
      AM(1)= 0.
      CM(1)= -CA*AZS(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) + CA*PME(i,j)/dt
      FM(1)= bio(i,j,1)  

*     -------------- Deep water -----------------------------------
      DO K=2, KB-1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

      AM(K)= -CA*CG( n)*AZS(i,j,k-1)/HZ(k-1)
      CM(K)= -CA*CG(np)*AZS(i,j,k  )/HZ(k)
      BM(K)= 1. -AM(K) -CM(K)
      FM(K)= bio(i,j,k)

      END DO

*     ------------- Ocean bottom ----------------------------------
      CA=2.*DT
      Q=1./HZ(KB-1)

      AM(KB)= -CA*Q*Q*AZS(i,j,kb-1)
      BM(KB)= 1. -AM(KB)
      CM(KB)= 0.
      FM(KB)= bio(i,j,kb)  

*     ------------- Matrix factorization --------------------------

      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)

      Do K=1,KB
      bio(i,j,k)= MAX(0.,RKSI(K))
      End do

      END IF
      END DO
      END DO

      RETURN
      END

      subroutine bio_recovery

*     Version 17.02.2014

      INCLUDE 'Slo2.fi'

*     Methanetporhs recovery

	Ctime= 3600.*24.*365.    ! Recovery time = 1 year, [sec]
	cc   = dt/Ctime

      do j=1,jl
	do i=1,il
	if(km2(i,j) .GT. 0) then

	do k=1,km2(i,j)    

      if(o2(i,j,k).GE.1.e3  .and. bio(i,j,k).GT.1.e-3) then ! O2 > 10 microM
	bio(i,j,k) = (cc +bio(i,j,k))/(1.+ cc)
	end if

	end do

	end if

	end do
	end do

      return
	end
