      Subroutine Forcing3(iyear,iday,L,dthr,TA,Q2m,WX,WY,Pa,il,jl,
     #                    il1,jl1,si,co,hx,hy,r,om,serv,serv1) 

*     Version 4.6 09.02.2014.

c     AOMIP forcing, independent on the solution:
c     Air temperature Ta [C], pressure Pa [Pa], 
c     wind (WX,WY) [m/s],
c     humidity at 2m Q2m [g/g], precipitation Pr [g/cm2/s].
c     Cloudness with spatial variations.
c     Rivers are treated by special procedures.

c ---------------------------------------------------------------------------
c
c Actually this procedure works for dthr < 24 and time is shifted half a day!
c
c ---------------------------------------------------------------------------

	real TA(il,jl),Pa(il,jl),
     #     WX(0:il1,0:jl1),WY(0:il1,0:jl1), Q2m(il,jl),
     #     si(0:jl1),co(0:il1), serv(il,jl), serv1(il,jl)
	character*80 name


      aaa= (24./dthr-L)/(24./dthr -1.)

	iyearS=iyear
	idayS =iday

      IF(iyear .EQ. 2014) then  ! The end of the Run!
	iyearS=2013
	idayS =365
	END IF

      iyear1S=iyear
	IF(iday .LT. 365) then
	iday1S=iday +1
	else
	iyear1S=iyear+1
	iday1S=1
	ENDIF

	IF(iyear1S.EQ.2014) then  ! The end of the Run!
	iyear1S=2013
	iday1S= 365
	END IF

c     Atmosphere pressure according to the time.
	name='C:\SLPMOD\pa1948.dat'
      write(name(13:16),'(i4.4)') iyearS
	open (51,file=name,access='direct',form='unformatted',
     *      recl=il*jl)
	read(51,rec=idayS) Serv
	close(51)

	name='C:\SLPMOD\pa1948.dat'
      write(name(13:16),'(i4.4)') iyear1S
	open (51,file=name,access='direct',form='unformatted',
     *      recl=il*jl)
	read(51,rec=iday1S) Serv1
	close(51)

	do j=1,jl
	do i=1,il
      Pa(i,j)=aaa*Serv(i,jl-j+1)+(1.-aaa)*Serv1(i,jl-j+1)
cc	call random_seed()
cc	call random_number(q)
cc	Pa(i,j)=Pa(i,j)*(1.0+0.1*q)
	end do
	end do

c     Atmosphere temperature according to the time.
	name='C:\TAMOD\ta1948.dat'
      write(name(12:15),'(i4.4)') iyearS
	open (51,file=name,access='direct',form='unformatted',
     *      recl=il*jl)
	read(51,rec=idayS) Serv
	close(51)

	name='C:\TAMOD\ta1948.dat'
      write(name(12:15),'(i4.4)') iyear1S
	open (51,file=name,access='direct',form='unformatted',
     *      recl=il*jl)
	read(51,rec=iday1S) Serv1
	close(51)

	do i=1,il
	do j=1,jl
      Ta(i,j)= aaa*Serv(i,jl-j+1)+(1.-aaa)*Serv1(i,jl-j+1)
	end do
	end do

      do i=1,il
      do j=1,jl
      TC= Ta(i,j)
c     The -12C limit is by E.V. Volodin. AMIP. See QMAX function.
	if (tc .GE. 0.) then ! Warm atmosphere.
c     In AOMIP - only warm atmosphere.
c      ea= 0.9* 611.*(10.**(7.5*TC/(273.16-35.86+TC)))   

	ea=10.**(2.+(0.7859+0.03477*TC)/(1.+0.00412*TC)) 
	ea= 0.9*ea                                         

	else ! Atmosphere with snow and ice. 
      ea=10.**(0.00422*Tc+2.+(0.7859+0.03477*Tc)/(1.+0.00412*Tc)) 
	ea= 0.9*ea                                                  
c      ea= 0.9* 611.*10.**(9.5*TC/(273.16- 7.66+TC))
	end if  ! TA > 0C

      Q2m(i,j)= 0.622*ea/Pa(i,j)
      end do
      end do
	
      call pwind(Pa,si,co,wx,wy,il,jl,il1,jl1,
     #           hx,hy,r,om)

      return
      end

      subroutine shortwave(SW,Pa,Q2m,m1day,mhour,dthr,cloud,il,jl)

c     version 01.02.2013 with correct Fi and Stime, real day

c     Short wave radiation calculation. Cloud is the cloud cover.
c     Zillman, 1972, Parkinson & Washington, 1979.

	parameter (Solar=1353.*1.e3)  ! Solar constant 
	parameter (hx=1.,hy=1.)       ! Grid size in degrees.
	parameter (pi= ACOS(-1.0))    ! pi=3.14159265
	parameter(Alpha=0., Beta=-90., Gamma=0.)
	dimension SW(il,jl),Pa(il,jl),Q2m(il,jl),cloud(il,jl)
	dimension R_Forward(3,3)
	parameter(FiSmodel=-25., aLoWmodel=-16.)


c     Rotation matrix into geographical coordinates
      cosa= COS(Alpha*pi/180.)
	sina= SIN(Alpha*pi/180.)
      cosb= COS(Beta*pi/180.)
	sinb= SIN(Beta*pi/180.)
      cosg= COS(Gamma*pi/180.)
	sing= SIN(Gamma*pi/180.)
 
      R_Forward(1,1)=  cosg*cosb*cosa - sing*sina
      R_Forward(1,2)= -sing*cosb*cosa - cosg*sina
      R_Forward(1,3)=  sinb*cosa
      R_Forward(2,1)=  cosg*cosb*sina + sing*cosa
      R_Forward(2,2)= -sing*cosb*sina + cosg*cosa
      R_Forward(2,3)=  sinb*sina
      R_Forward(3,1)= -cosg*sinb
      R_Forward(3,2)=  sing*sinb
      R_Forward(3,3)=  cosb


      do i=1,il
	do j=1,jl

	Xin= aLoWmodel+(i-1.)*hx
	Yin= FiSmodel +(jL-j)*hy
	call rotation(Xin,Yin,Xout,Fi,R_Forward)

c     Vapor pressure [Pa]
      evapor= Q2m(i,j)*Pa(i,j)/0.622

c     Solar time at the location
	stime=dthr*(REAL(mhour)-0.5) ! GMT
	stime=stime +Xout/15.

	day= REAL(m1day) + stime/24.
	if(stime .GT. 24.) stime= stime -24.

c     Declination angle
      delta= 23.44*cos((172.-day)*pi/180.)
	  sd=sind(delta)
	  cd=cosd(delta)

	CZ= MAX(0.,sind(fi)*sd+cosd(fi)*cd*cos((12.-stime)*pi/12.))
	Q0=Solar*CZ**2/((CZ+2.7)*evapor*1.e-5 +1.085*CZ+0.10)

c     Clouds
      SW(i,j)=Q0*(1.-0.6*(cloud(i,j)**3))

c	if(i.eq.17.and.j.eq.24) .and.SW(i,j).GT.0.) 
c     &write (*,*) day, stime, Q0, SW(i,j)

	end do
	end do
	return
	end

	subroutine Rotation(Xin,Yin,Xout,Yout,R)
	dimension R(3,3)
	parameter (Pi= ACOS(-1.0))

c     A coordinate transformation from coordinates
c     with rotated grid, Euler angles (X_IN_ANGLE,Y_IN_ANGLE,0),
c     to geographical coordinates (long, teta).
c     R - Euler matrix of rotation.
c     (Xin, Yin) - initial model coordinates,
c     (Xout,Yout) - geographical coordinates.

c     avoid troble of an exactly zero angle by adding the offset.      
c     and converting to radians.
      eps= 1.e-5
	  Xin= Pi*(Xin +eps)/180.
	  Yin= Pi*(Yin +eps)/180.

c     Spherical coordinates to cartesian coordinates
      xx= COS(Yin)*COS(Xin)
	yy= COS(Yin)*SIN(Xin)
	zz= SIN(Yin)

c     New cartesian coordinates
      Xnew= R(1,1) * xx
     &     +R(1,2) * yy
     &     +R(1,3) * zz

      Ynew= R(2,1) * xx
     &     +R(2,2) * yy
     &     +R(2,3) * zz

      Znew= R(3,1) * xx
     &     +R(3,2) * yy
     &     +R(3,3) * zz

c     New angles teta_new, costn, phi_new
      teta_new = ASIN(Znew)
	costn= MAX(eps,SQRT(1.0 - Znew**2))

	IF      ((Xnew.GT.0.) .AND. (Ynew.GT.0.)) THEN

	     IF(Xnew.LT.Ynew) THEN
		 phi_new= ACOS(Xnew/costn)
		 ELSE
		 phi_new= ASIN(Ynew/costn)
		 ENDIF

      ELSEIF ((Xnew.LT.0.) .AND. (Ynew.GT.0.)) THEN
      
	      IF(ABS(Xnew).LT.Ynew) THEN
		  phi_new= pi-ACOS(ABS(Xnew)/costn)
		  ELSE
		  phi_new= pi-ASIN(    Ynew /costn)
		  ENDIF

      ELSEIF ((Xnew.LT.0.) .AND. (Ynew.LT.0.)) THEN

	      IF(ABS(Xnew).LT.ABS(Ynew)) THEN
		  phi_new= -pi+ACOS(ABS(Xnew)/costn)
		  ELSE
		  phi_new= -pi+ASIN(ABS(Ynew)/costn)
		  ENDIF

      ELSEIF ((Xnew.GT.0.) .AND. (Ynew.LT.0.)) THEN

	      IF( Xnew.LT.ABS(Ynew)) THEN
		  phi_new= -ACOS(ABS(Xnew)/costn)
		  ELSE
		  phi_new= -ASIN(ABS(Ynew)/costn)
		  ENDIF

	ENDIF

c     New spherical coordinates.

      Xout= Phi_new  *180./Pi
	Yout= teta_new *180./Pi
ccc        IF( Xout.LT. 0.) Xout= Xout +360.

c     avoid troble of an exactly zero angle by subtracting the offset.      
	Xout= Xout -eps
	Yout= Yout -eps

	END

      subroutine pwind(Pa,si,co,windx,windy,il,jl,il1,jl1,
     #                 hx,hy,r,om)
*     AOMIP Wind velocity (in M/Sec) Calculation.
*     version 12.02.2016.


      dimension  windx(0:il1,0:jl1), windy(0:il1,0:jl1),
     #           si(0:jl1), co(0:il1)
      dimension  Pa(il,jl)
      roa=1.30e-2

      do j=1,jl
      s0=si(j)

      do i=1,il
      cor= -2.*om*co(i)*S0
      if(i.gt.1.and.i.lt.il) then
      derx= .5*(pa(i+1,j)-pa(i-1,j))
      else
          if(i.eq.1 ) derx=pa( 2,j)-pa(   1,j)
          if(i.eq.il) derx=pa(il,j)-pa(il-1,j)
      end if

      if(j.gt.1.and.j.lt.jl) then
      dery= .5*(pa(i,j+1)-pa(i,j-1))
      else
          if(j.eq.1 ) dery=pa(i,2 )-pa(i,   1)
          if(j.eq.jl) dery=pa(i,jl)-pa(i,jl-1)
      end if

c     Geostrophical wind
      wgy=  derx/(roa*cor*s0*r*hx)
      wgx= -dery/(roa*cor*r*hy)
	wgmod= SQRT(wgx**2+wgy**2)

	if(wgmod. GT. 60.) then  ! Wind cut-off
	wgy= wgy*60./wgmod
	wgx= wgx*60./wgmod
	end if

*-----  Scheme 1 - AOMIP  ------------------
* Hunke, E. C., and M. M. Holland (2007), 
* Global atmospheric forcing data for Arctic ice-ocean modeling, 
* J. Geophys. Res., 112, C04S14, doi:10.1029/2006JC003640.

      if( wgmod .LT. 15.)then
	a= -30.
      reduc=0.8
	else
	a= -20.
      reduc=0.7
	end if
		
c	IF(wgmod .LT. 20.)then
c	a=-30.+wgmod
c	else
c	    IF(wgmod .LT. 30.) then
c	    a=-20.+0.5*wgmod
c	    else
c	      IF(wgmod .LT. 40.) then
c		  a=-12.5+0.25*wgmod
c	      else
c	      a=-2.5
c	      END IF
c	    END IF
c	END IF

*-----  Scheme 2   -------------------------

c      G1=-0.0709*wgmod
c	G2=-0.4840*wgmod
c	a= -(41.264*EXP(G1)-11.268*EXP(G2))*pi/180.
c	reduc=0.8

c	if( wgmod .LT. 15.)then
c	reduc=0.7
c	else
c	reduc=0.8
c	end if

*-------------------------------------------
      sina= sind(a)
      cosa= cosd(a)

      windx(i,j)= reduc*(wgx*cosa -wgy*sina)
      windy(i,j)= reduc*(wgx*sina +wgy*cosa)

      end do
      end do

      return
      end

      subroutine TSInt(imonth)

c     Observed T,S - for liquid boundaries and flux-correction.
c     Version 07.06.2015

c     PHC 3.0 data are for in situ temperature - so we need to recalculate it to potential one.
c     Bryden, 1973, 

c     or 

c     McDougall, T. J.,D. R. Jackett, D. G. Wright, and R. Feistel, 2003: Accurate
c     and computationally efficient algorithms for potential temperature
c     and density of seawater. J. Atmos. Oceanic Technol., 20, 730–741.
c     Tpot (S, T, p, pref )= T +(p - pref)P(S, T, [p + pref])
c     P(S, T, [p + pref ]) = a1 + a2*S + a3*[p + pref] + a4*T + a5*ST + a6*TT +a7*T[p + pref],  
c     Coefficients corrected in
c     Jackett, D. R., McDougall, T. J., Feistel, R., Wright, 
c     D. G., and Griffies, S. M.: 
c     Algorithms for density, potential temperature, conservative
c     temperature, and freezing temperature of seawater, 
c     Journal of Atmospheric and Oceanic Technology, 23, 1709–1728, 2006.	

      include 'Slo2.fi'

      row = 1.020     ! Water density, g/cm3
	g=980.

c     Polynomial coefficients.

      a1 =  8.654839e-6
      a2 = -1.416363e-6
      a3 = -7.382865e-9
      a4 = -8.382414e-6
      a5 =  2.839334e-8
      a6 =  1.778040e-8
      a7 =  1.711556e-10

	open (51,file='SPHC.dat',access='direct',form='unformatted',
     *      recl=il*jl*kl)
	read(51,rec=imonth) Serv2
	close(51)
      do i=1,il
      do j=1,jl
      do k=1,kl
      sobs(i,j,k)= serv2(i,jl-j+1,k)
      end do
      end do
      end do
*     GIN Sea
c      do k=1,kl
c      sobs(10,jl,k)= sobs(11,jl,k)
c      end do


	open (51,file='TPHC.dat',access='direct',form='unformatted',
     *      recl=il*jl*kl)
	read(51,rec=imonth) Serv2
	close(51)
      do i=1,il
      do j=1,jl
      do k=1,kl
	Tin = serv2(i,jl-j+1,k)
	Tin2= Tin**2
	Tin3= Tin**3
	ppp = 1.e-5*g*row*z(k)
	sss = sobs(i,j,k)

c     Potential temperature is by Bryden H.L., 1973. See A. Gill textbook.
c     Too narrow salinity interval [30,40] and temperature [2,30].

*      tobs(i,j,k)= Tin 
*     &   - ppp*(3.6504e-4 +8.3198e-5*Tin-5.4065e-7*Tin2+4.0274e-9*Tin3)
*     &   - ppp*(sss-35.)*(1.7439e-5-2.9778e-7*Tin)
*     &   - ppp*ppp*(8.9309e-7-3.1628e-8*Tin+2.1987e-10*Tin2)
*     &   + 4.1057e-9*(sss-35.)*ppp*ppp
*     &   - (-1.6056e-10 +5.0484e-12*Tin)*ppp**3

c     Potential temperature by MJWF03 and Jackett, et. al. 2006. 
c     Pressure converted in dbar. Reference pressure is 0 bar.
c     Wide T,S,p interval.

      Poly=a1 +a2*sss +a3*ppp+
     &     a4*Tin +a5*Sss*Tin +a6*Tin2 +a7*Tin*ppp
      tobs(i,j,k)= Tin +Poly*ppp
     
      end do
      end do
      end do
*     GIN Sea
c	do k=1,kl
c      tobs(10,jl,k)= tobs(11,jl,k)
c      end do

      return
      end

	      
	subroutine river_as_rain
     *           (month,River,il,jl,nt3,si,cg,il1,jl1,kl,r,hx,hy)

*     version 29.03.2015.

*     AOMIP River runoff, Prange, cm3/s.
*     Gordeev, Peterson et.al. 2002, Gordeev, Rachold, 2003 - *0.5
*     Monthly mean climatology.
*     Rivers as rain, incerted in the specified grid cell, cm/s

	dimension nt3(0:il1,0:jl1,kl), cg(13)
	dimension si(0:jl1), River(il,jl)
	dimension TW(12)                 ! Month time weight
	data TW/4*0.02,0.1,0.31,0.17,0.12,0.1,0.07,0.03,0.02/

	c6= 1./6.
	RR= 1.e-6*R*R               ! Normalization
      C= 12.*TW(month)/(RR*hx*hy) ! year mean by planar area element

c     Ob and Yenisey, Puhr, Taymyra, Pyasina

	Sect= 0.
	do i=33,34
	do j=27,29
	Sect= Sect+ c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	end do

      do i=33,34
	do j=27,29
      River(i,j)= 45918.9410*C/Sect
	end do
	end do
	
	write(*,*) River(34,28), Sect

c     Lena, Khatanga, Olenek, Yana

	Sect= 0.
	do i=30,31
	do j=11,14
	Sect= Sect+ c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	end do

      do i=30,31
      do j=11,14
      River(i,j)= 26964.1602*C/Sect
	end do
	end do

c     McKenzie
	Sect= 0.
      i=2
	do j=9,10
	Sect= Sect +c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	
	do j=9,10
      River(i,j)= 11838.9104*C/Sect
      end do

c     Dvina and Mezen
      Sect=0.
	j=43
      do i=33,34
	Sect= Sect +c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	
	do i=33,34
      River(i,j)= 4328.8674*C/Sect
      end do

c     Pechora
      Sect= 0.
      i=35
      do j=35,36
	Sect= Sect +c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	
	do j=35,36
      River(i,j)= 6072.7550*C/Sect
      end do

c     Kolyma+Indigirka
      i=24
      Sect= 0.
	do j=5,6
	Sect= Sect +c6*cg(abs(nt3(i,j,1)))*si(j)
	end do
	
	do j=5,6
      River(i,j)= 6110.2248*C/Sect
      end do

      return
	end

      Subroutine UB_SLO(ub,unept,vnept,nt3,km2,dzext,Pext,PA,wx,wy,
     *                  tobs, sobs, az,il,jl,il1,jl1,kl,klp,
     *                  aice,u,v,uice,vice,CDgwd,mgrad,
     *                  hx,hy,r,si,co,hz,z,g,om,row,roa,CDw)

*     ------------ Version 5.8 for free slip conditions -----------
*           Geostrophic + Ekman velocity at long open boundaries
*               Either observed level or dynamical level

*                                 31.01.16
*
*                      Positive values are OUT the area.
*
*            Narrow passages velocity is specified by my own hands -
*                   it's the main feature of the program.
*
      Dimension ub(0:il1,0:jl1,kl),nt3(0:il1,0:jl1,kl),km2(0:il1,0:jl1)
      Dimension si(0:jl1),co(0:il1), hz(kl), z(klp), unept(0:il1,0:jl1),
     &          vnept(0:il1,0:jl1), dzext(il,jl), Pext(il,jl,kl), 
     &          PA(il,jl),
     * tobs(0:il1,0:jl1,kl),sobs(0:il1,0:jl1,kl),
     * wx(0:il1,0:jl1),wy(0:il1,0:jl1),
     * aice(0:mgrad,0:il1,0:jl1), u(0:il1,0:jl1,kl), v(0:il1,0:jl1,kl),
     * CDgwd(il,jl),az(il,jl,kl),
     * uice(0:il1,0:jl1), vice(0:il1,0:jl1)
     
    

*     Parameters for neptune parameterization
      RN2= (3.e5)**2        ! length scale squared - usually 5 km
	alpha_neptune = 1.0   ! To account for resolution - stream width 30km 
      
*     Parameter of pressure calculation - npres=1 - Level data, elese - dynamical meth.      
	npres =1

*     Vertical viscisity - assumed to be constant for Ekman layer calculation.
*     Small values - shallow wind layer, large values - deep layer.
*     May be estimated also by Az(i,j,k) in upper ocean.
*      Vert_visc= 10. - 100.

      asr=hx/hy
	ub= 0.


*     Norwegian sea
*     -------------------------------------------------

      s0=si(jl)

      do i=11,23

	if(npres. eq. 1) then
*     Pressure based on monthly mean density and mean sea level
      Pext(i,jl,1)= row*dzext(i,jl) +10.*PA(i,jl)/g

	else

*     Pressure based on dynamical method formulas - integration up to the bottom
      Pext(i,jl,1)= 0.

      do k=2,min(11,km2(i,jl))
      IF( nt3(i,jl,k) .LT. 0) then
	ppp= 1.e-5*(g*row*z(k) +10.*PA(i,jl))
      ro1= sigma_t(tobs(i,jl,k),sobs(i,jl,k),ppp)
	ppp= 1.e-5*(g*row*z(k-1)+10.*PA(i,jl))
      ro2= sigma_t(tobs(i,jl,k-1),sobs(i,jl,k-1),ppp)
	Ro1= 0.5*(Ro1+Ro2)    ! Density in the middle point
	Pext(i,jl,1)= Pext(i,jl,1) -1.e-3*REAL(Ro1)*hz(k-1)
	end if
	end do

	end if  ! npres

      do k=2,km2(i,jl)
      IF( nt3(i,jl,k) .LT. 0) then
	ppp= 1.e-5*(g*row*z(k)+10.*PA(i,jl))
      ro1= sigma_t(tobs(i,jl,k),sobs(i,jl,k),ppp)
	ppp= 1.e-5*(g*row*z(k-1)+10.*PA(i,jl))
      ro2= sigma_t(tobs(i,jl,k-1),sobs(i,jl,k-1),ppp)
	Ro1= 0.5*(Ro1+Ro2)    ! Density in the middle point
	Pext(i,jl,k)= Pext(i,jl,k-1) +1.e-3*REAL(Ro1)*hz(k-1)

	end if
	end do ! k


	end do ! i


      do i=11,23

	cor=  -2.*om*s0*co(i)   
	coef= g/(r*s0*hx*cor*row)

*     Wind Stress Components
      wmod= SQRT(wx(i,jl)**2 + wy(i,jl)**2)
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      tx= cdrag*wx(i,jl)
      ty= cdrag*wy(i,jl)

*     Ice drift stress components
	drag=CDgwd(i,jl) +row*Cdw*SQRT( (uice(i,jl)-u(i,jl,1))**2
     +                               +(vice(i,jl)-v(i,jl,1))**2 )
	tx = Aice(0,i,jl)*tx  
     +   + drag*(1.-Aice(0,i,jl))*(uice(i,jl)-u(i,jl,1))
	ty = Aice(0,i,jl)*ty 
     +   + drag*(1.-Aice(0,i,jl))*(vice(i,jl)-v(i,jl,1))


*     Vertical viscosity estimation upper 100m
      vert_visc=0.
      do k=1,10
	vert_visc= vert_visc +az(i,jl,k)
	end do
	vert_visc= vert_visc/10.

      alpha= SQRT(ABS(cor)/(2.*vert_visc))


      do k=1,km2(i,jl)
      IF( nt3(i,jl,1) .LT. 0) then


*     Pressure gradient
      derx=0.
      n=abs(nt3(i,jl,k))
	if(n.eq.10) derx= (Pext(i+1,jl,k)-Pext(i,jl,k))
      if(n.EQ.9)  derx= 0.5*(Pext(i+1,jl,k)-Pext(i-1,jl,k))
	if(n.EQ.11) derx= (Pext(i,jl,k)-Pext(i-1,jl,k))

*     Ekman component - the real wind.  
*     Formula is for Northern hemisphere.

      eaz  = EXP(-alpha*z(k)) /(2.*vert_visc*alpha)  ! Az=10
      sinaz= SIN(alpha*z(k))
      cosaz= COS(alpha*z(k))
      Vijk = eaz*((ty+tx)*cosaz - (ty-tx)*sinaz)
      Uijk = eaz*((tx-ty)*cosaz - (tx+ty)*sinaz)

*     Normal projection
	if(n.eq.11) Ekman= ( hy*Uijk + hx*s0*Vijk)/(hy +hx*s0)
      if(n.EQ.9)  Ekman= Vijk
	if(n.eq.10) Ekman= (-hy*Uijk + hx*s0*Vijk)/(hy +hx*s0)

	WW= -alpha_neptune*cor*RN2*vnept(i,jl)
      ub(i,jl,k) =  WW +coef*derx +Ekman

c      write(*,*) 'NS:', i,k,ww,coef*derx,Ekman
	END IF
      end do

      end do

cc      write(1,*)' Norw. sea:', (ub(i,jl,10),i=11,22)

*     Canadian Archipelago - ASOF - total ~ 2.4 Sv
*     Zhang 1998 - Nares 0.7 Now - 0.8.
*     --------------------------
      Sect=0.
      j=28
      s0=si(j)

      kk2=km2(10,j)
      IF( nt3(10,j,1) .LT. 0) then

      do k=1,kk2
      n= nt3(10,j,k)
      if( k.lt.kk2) np=nt3(10,j,k+1)

      IF( k.eq.1) then
      divk= 0.0
      hzk1= 0.0
      else
      hzk1= hz(k-1)
      DivK=0.
      IF(n.EQ.-6) DivK=-6.
      IF(n.EQ.-13) DivK=-3.*(1.+s0*asr)
      IF(n.EQ.-10) DivK=-3.*(1.+s0*asr)
      end if

      if( k.eq.kk2) then
      hzk= 0.0
      divkp= 0.0
      else
      hzk= hz(k)
      DivKP=0.
      IF(np.EQ.-6) DivKP=-6.
      IF(np.EQ.-13) DivKP=-3.*(1.+s0*asr)
      IF(np.EQ.-10) DivKP=-3.*(1.+s0*asr)
      end if

      Sect = Sect + (hzk1*divk+hzk*divkp)/12.

      end do
      END IF
                   b= -0.8e+12/(Sect*r*hy)

      j=28
      do k=1,km2(10,j)
      ub(10,j,k) = b
      end do
cc      write(1,*)' Nares:', b

*     Zhang 1998 - M'Clure 0.8
*     --------------------------
      Sect=0.
      do j=16,17
      s0=si(j)

      kk2=km2(3,j)
      IF( nt3(3,j,1) .LT. 0) then

      do k=1,kk2
      n= nt3(3,j,k)
      if( k.lt.kk2) np=nt3(3,j,k+1)

      IF( k.eq.1) then
      divk= 0.0
      hzk1= 0.0
      else
      hzk1= hz(k-1)
      DivK=0.
      IF(n.EQ.-6) DivK=-6.
      IF(n.EQ.-13) DivK=-3.*(1.+s0*asr)
      IF(n.EQ.-10) DivK=-3.*(1.+s0*asr)
      end if

      if( k.eq.kk2) then
      hzk= 0.0
      divkp= 0.0
      else
      hzk= hz(k)
      DivKP=0.
      IF(np.EQ.-6) DivKP=-6.
      IF(np.EQ.-13) DivKP=-3.*(1.+s0*asr)
      IF(np.EQ.-10) DivKP=-3.*(1.+s0*asr)
      end if

      Sect = Sect + (hzk1*divk+hzk*divkp)/12.

      end do
      END IF

	end do
                   b= -0.8e+12/(Sect*r*hy)

      do j=16,17
      do k=1,km2(3,j)
      ub(3,j,k) = b
      end do
	end do
cc      write(1,*)' MClure:', b

*     Middle part of the Canadian Archipelago - 0.8
*     ---------------------------------------------
      Sect=0.
      do j=18,21
      s0=si(j)

      kk2=km2(6,j)
      IF( nt3(6,j,1) .LT. 0) then

      do k=1,kk2
      n= nt3(6,j,k)
      if( k.lt.kk2) np=nt3(6,j,k+1)

      IF( k.eq.1) then
      divk= 0.0
      hzk1= 0.0
      else
      hzk1= hz(k-1)
      DivK=0.
      IF(n.EQ.-6) DivK=-6.
      IF(n.EQ.-13) DivK=-3.*(1.+s0*asr)
      IF(n.EQ.-10) DivK=-3.*(1.+s0*asr)
      end if

      if( k.eq.kk2) then
      hzk= 0.0
      divkp= 0.0
      else
      hzk= hz(k)
      DivKP=0.
      IF(np.EQ.-6) DivKP=-6.
      IF(np.EQ.-13) DivKP=-3.*(1.+s0*asr)
      IF(np.EQ.-10) DivKP=-3.*(1.+s0*asr)
      end if

      Sect = Sect + (hzk1*divk+hzk*divkp)/12.

      end do ! k
      END IF
	end do ! j
                   b= -0.8e+12/(Sect*r*hy)

      do j=18,21
      do k=1,km2(6,j)
      ub(6,j,k) = b
      end do
	end do

cc      write(1,*)' Middle part of Archipelago:', b



*     Bering passage 0.6 Sv
*     Zhang - 0.8
*     ---------------------

      s0=si(1)

      do i=11,15

	if(npres. eq. 1) then
*     Pressure based on monthly mean density and mean sea level
      Pext(i,1,1)= row*dzext(i,1) +10.*PA(i,1)/g

	else

*     Pressure based on dynamical method formulas - integration up to the bottom
      Pext(i,1,1)= 0.

      do k=2,min(11,km2(i,1))
      IF( nt3(i,1,k) .LT. 0) then
	ppp= 1.e-5*(g*row*z(k) +10.*PA(i,1))
      ro1= sigma_t(tobs(i,1,k),sobs(i,1,k),ppp)
	ppp= 1.e-5*(g*row*z(k-1) +10.*PA(i,1))
      ro2= sigma_t(tobs(i,1,k-1),sobs(i,1,k-1),ppp)
	Ro1= 0.5*(Ro1+Ro2)    ! Density in the middle point
	Pext(i,1,1)= Pext(i,1,1) -1.e-3*REAL(Ro1)*hz(k-1)
	end if
	end do

	end if

      do k=2,km2(i,1)
      IF( nt3(i,1,k) .LT. 0) then
	ppp= 1.e-5*(g*row*z(k)  +10.*PA(i,1))
      ro1= sigma_t(tobs(i,1,k),sobs(i,1,k),ppp)
	ppp= 1.e-5*(g*row*z(k-1)  +10.*PA(i,1))
      ro2= sigma_t(tobs(i,1,k-1),sobs(i,1,k-1),ppp)
	Ro1= 0.5*(Ro1+Ro2)    ! Density in the middle point
	Pext(i,1,k)= Pext(i,1,k-1) +1.e-3*REAL(Ro1)*hz(k-1)

	end if
	end do
	end do

      do i=11,15

	cor= -2.*om*s0*co(i)
	coef= g/(r*s0*hx*cor*row)

*     Wind Stress Components
      wmod= SQRT(wx(i,1)**2 + wy(i,1)**2)
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      tx= cdrag*wx(i,1)
      ty= cdrag*wy(i,1)

*     Ice drift stress components
	drag=CDgwd(i,1) +row*Cdw*SQRT( (uice(i,1)-u(i,1,1))**2
     +                              +(vice(i,1)-v(i,1,1))**2 )
	
	tx = Aice(0,i,1)*tx  
     +   + drag*(1.-Aice(0,i,1))*(uice(i,1)-u(i,1,1))
	ty = Aice(0,i,1)*ty 
     +   + drag*(1.-Aice(0,i,1))*(vice(i,1)-v(i,1,1))


*     Vertical viscosity estimation upper 100m
      vert_visc=0.
      do k=1,10
	vert_visc= vert_visc +az(i,1,k)
	end do
	vert_visc= vert_visc/10.

      alpha= SQRT(ABS(cor)/(2.*vert_visc))

      do k=1,km2(i,1)
      IF( nt3(i,1,1) .LT. 0) then

*     Pressure gradient
      derx=0.
      n=abs(nt3(i,1,k))
	if(n.eq.13) derx= (Pext(i+1,1,k)-Pext(i,1,k))
      if(n.EQ.7)  derx= 0.5*(Pext(i+1,1,k)-Pext(i-1,1,k))
	if(n.EQ.12) derx= (Pext(i,1,k)-Pext(i-1,1,k))

*     Ekman component - the real wind. Visc. coeff. assumed constant and 100. 
*     Formula is for Northern hemisphere.

      eaz  = EXP(-alpha*z(k)) /(2.*vert_visc*alpha)
      sinaz=SIN(alpha*z(k))
      cosaz=COS(alpha*z(k))
      Vijk =eaz*((ty+tx)*cosaz - (ty-tx)*sinaz)
      Uijk =eaz*((tx-ty)*cosaz - (tx+ty)*sinaz)

*     Normal projection
	if(n.eq.13) Ekman= (-hy*Uijk - hx*s0*Vijk)/(hy +hx*s0)
      if(n.EQ.7)  Ekman= -Vijk
	if(n.eq.12) Ekman= ( hy*Uijk - hx*s0*Vijk)/(hy +hx*s0)

	WW= alpha_neptune*cor*RN2*vnept(i,1)
      ub(i,1,k) =  WW -coef*derx +Ekman

c      write(*,*) 'BP:', i,k,ww,-coef*derx,Ekman

	END IF
      end do

      end do

*     Denmark Strait 
*     -------------------------------------------------------

      i=7

      do j=43,46

      if(npres .eq. 1) then
*     Pressure based on monthly mean density and mean sea level
      Pext(i,j,1)= row*dzext(i,j) +10.*PA(i,j)/g

	else

*     Dynamical level
      Pext(i,j,1)= 0.

      do k=2,min(4,km2(i,j))
      IF( nt3(i,j,k) .LT. 0) then
	ppp= 1.e-5*(g*row*z(k)  +10.*PA(i,j))
      ro1= sigma_t(tobs(i,j,k),sobs(i,j,k),ppp)
	ppp= 1.e-5*(g*row*z(k-1)  +10.*PA(i,j))
      ro2= sigma_t(tobs(i,j,k-1),sobs(i,j,k-1),ppp)
	Ro1= 0.5*(Ro1+Ro2)    ! Density in the middle point
	Pext(i,j,1)= Pext(i,j,1) -1.e-3*REAL(Ro1)*hz(k-1)

	end if
	end do

	end if


      do k=2,km2(i,j)
      IF( nt3(i,j,k) .LT. 0) then
	ppp= 1.e-5*(g*row*z(k)   +10.*PA(i,j))
      ro1= sigma_t(tobs(i,j,k),sobs(i,j,k),ppp)
	ppp= 1.e-5*(g*row*z(k-1)  +10.*PA(i,j))
      ro2= sigma_t(tobs(i,j,k-1),sobs(i,j,k-1),ppp)
	Ro1= 0.5*(Ro1+Ro2)    ! Density in the middle point
	Pext(i,j,k)= Pext(i,j,k-1) +1.e-3*REAL(Ro1)*hz(k-1)

	end if
	end do
	end do


      do j=43,46

	cor= -2.*om*si(j)*co(i)
	coef= g/(r*hy*cor*row)

*     Wind Stress Components
      wmod= SQRT(wx(i,j)**2 + wy(i,j)**2)
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      tx= cdrag*wx(i,j)
      ty= cdrag*wy(i,j)

*     Ice drift stress components
	drag=CDgwd(i,j) +row*Cdw*SQRT( (uice(i,j)-u(i,j,1))**2
     +                              +(vice(i,j)-v(i,j,1))**2 )
	tx = Aice(0,i,j)*tx  
     +   + drag*(1.-Aice(0,i,j))*(uice(i,j)-u(i,j,1))
	ty = Aice(0,i,j)*ty 
     +   + drag*(1.-Aice(0,i,j))*(vice(i,j)-v(i,j,1))

*     Vertical viscosity estimation upper 100m
      vert_visc=0.
      do k=1,10
	vert_visc= vert_visc +az(i,j,k)
	end do
	vert_visc= vert_visc/10.

      alpha= SQRT(ABS(cor)/(2.*vert_visc))

      do k=1,km2(i,j)
      IF( nt3(i,j,1) .LT. 0) then

*     Pressure gradient
      derx=0.
      n=abs(nt3(i,j,k))
	if(n.eq.13) derx= (Pext(i,j+1,k)-Pext(i,j,k))
      if(n.EQ.6)  derx= 0.5*(Pext(i,j+1,k)-Pext(i,j-1,k))
	if(n.EQ.10) derx= (Pext(i,j,k)-Pext(i,j-1,k))

*     Ekman component - the real wind. Visc. coeff. assumed constant and 100. 
*     Formula is for Northern hemisphere.

      eaz  = EXP(-alpha*z(k)) /(2.*vert_visc*alpha)
      sinaz= SIN(alpha*z(k))
      cosaz= COS(alpha*z(k))
      Vijk =eaz*((ty+tx)*cosaz - (ty-tx)*sinaz)
      Uijk =eaz*((tx-ty)*cosaz - (tx+ty)*sinaz)

*     Normal projection
	if(n.eq.13) Ekman= (-hy*Uijk - hx*si(j)*Vijk)/(hy +hx*si(j))
      if(n.EQ.6)  Ekman= -Uijk
	if(n.eq.10) Ekman= (-hy*Uijk + hx*si(j)*Vijk)/(hy +hx*si(j))

	WW= alpha_neptune*cor*RN2*unept(i,j)
      ub(i,j,k) =  WW +coef*derx + Ekman 

c      write(*,*) 'DStr:', j,k,ww,coef*derx,Ekman

	END IF
      end do

      end do


      RETURN
      END
