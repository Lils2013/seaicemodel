      program FEMAO1_Methane

*     Version 17.02.2014.

*     methane evolution included + o2 + bio

*     New scheme of reading input data.
*     PHC monthly mean data.
*     Monthly mean data are saved on disc.
*     Fast ice and stop at narrow passages.
*     GW drag for ice-ocean.
*     Suitable for Flather, 1976, radiation condition.

      parameter (pi=3.1415926536)
      parameter (hxgr=1., hygr=1., dthr=1.0, finord=23., alwest=-16.)
      parameter (ncond=1)
      parameter (NRestor= 1)  ! T-S Restoring [0/1]
      parameter (lprint= 24)

*     Parameters for GMRES
      PARAMETER (LEN= 1084)
      PARAMETER (Krylov=7,KRP1=Krylov+1)

*     Watcom and NDP Fortrans Lfort=4, Compaq and UNIX - Lfort=1
*     For Real*8 - 2 times more
      Parameter (Lfort=1)

      INCLUDE 'Slo2.fi'
      common /diffus/ al,alt
	dimension MEC(12)
      dimension rs(il,jl),adz(-1:ilp,-1:jlp),serv(il,jl),serv1(il,jl)

      real Agm(2,0:il1,0:jl1)   ! Variable GMG coefficient at the triangle

*     Arrays for 3Dtransport

	dimension Sclr (0:il1,0:jl1,kl),Sclr1(0:il1,0:jl1,kl),
     *          Sclr2(0:il1,0:jl1,kl),SclrObs(0:il1,0:jl1,kl),
     *          Qcascad(il,jl,kl)



      REAL*8   B(LEN), SL(LEN), R0(LEN)
      REAL*8   SK(LEN), CK(KRP1,KRP1), GIV(2,KRP1),W0(LEN)
      REAL*8   Z0(LEN,KRP1),V0(LEN)

c     Monthly means.
      dimension umonth(0:il1,0:jl1,kl), vmonth(0:il1,0:jl1,kl), 
     *          wmonth(0:il1,0:jl1,kl),
     *          Tmonth(0:il1,0:jl1,kl), Smonth(0:il1,0:jl1,kl),
     *          CH4month(0:il1,0:jl1,kl),
     *          O2month (0:il1,0:jl1,kl),biomonth (0:il1,0:jl1,kl),
     *          Flux_CH4_To_A_month(0:il1,0:jl1),
     *  uimonth(0:il1,0:jl1), vimonth(0:il1,0:jl1),
     *  himonth(0:mgrad,0:il1,0:jl1), aimonth(0:mgrad,0:il1,0:jl1),
     *  hsmonth(0:mgrad,0:il1,0:jl1), dzmonth(-1:ilp,-1:jlp)

cc	dimension cloud_PW(12)

c     mean Dyn. level AVISO data
	dimension dzext_mean(il,jl)

*     Tides, Normal velocity
	DIMENSION AMPL(14), PH(14)   ! Incoming Tidal wave - Norwegian Sea.
	DIMENSION AMPL2(6), PH2(6)   ! Incoming Tidal wave - Denmark Str.
	DIMENSION AMPL3(5), PH3(5)   ! Incoming Tidal wave - Canadian A.
	DIMENSION AMPL4(5), PH4(5)   ! Incoming Tidal wave - Bering Pas.
*     Tides, Sea level
	DIMENSION DAMPL(14),DPH(14)  ! Incoming Tidal wave - Norwegian Sea.
	DIMENSION DAMPL2(6),DPH2(6)  ! Incoming Tidal wave - Denmark Str.
	DIMENSION DAMPL3(5),DPH3(5)  ! Incoming Tidal wave - Canadian A.
	DIMENSION DAMPL4(5),DPH4(5)  ! Incoming Tidal wave - Bering Pas.

      double precision rs, adz, time, dzmonth   !!, Tide_M2
      character*72 udat,tsdat,dzdat,diagdat,nt3dat,ubdat,km2dat,
     *             neptunedat,imask, CH4dat, CH4obsdat,FluxCH4dat,
     *             O2dat, O2obsdat, aidat,hisdat,uidat,tmpisdat



*======================================================================
      data udat/'u.dat'/ ,tsdat/'t.dat'/,
     *     dzdat/'dz.dat'/, diagdat/'diag.dat'/, nt3dat/'nt.dat'/,
     *     ubdat/'ub.dat'/, km2dat/'km2.dat'/, 
     *     aidat/'ai.dat'/, hisdat/'his.dat'/,uidat/'ui.dat'/,
     *     tmpisdat/'tmpis.dat'/,neptunedat/'neptune.dat'/,
     *     imask/'icemask.dat'/,
     *     CH4dat/'CH4.dat'/, CH4obsdat/'CH4obs.dat'/,
     *     O2dat/'O2.dat'/, O2obsdat/'O2obs.dat'/,
     *     FluxCH4dat/'FluxCH4.dat'/


      data cg/6.,4.,5.,4.,5.,4*3.,2.,1.,2.,1./
	data MEC/31,28,31,30,31,30,31,31,30,31,30,31/


cc	data cloud_PW/.81,.78,.81,.71,.78,.83,.85,.87,.85,.89,.78,.78/

c     Tidal Velocity
*     Norwegian sea
c      data ampl/19.91,31.06,19.39,5.86,3.82,3.07,2.56,2.10,
c     *	      1.95,1.83,2.00,2.00,3.11,4.00/

*     Modul
      data ampl/30.00,37.52,24.26,8.32,5.23,3.85,3.24,2.72,
     *	      2.33,3.49,6.59,5.90,4.83,4.24/

	data ph/287.89,241.27,217.64,244.17,260.17,274.58,280.48,
     *        289.80,307.68,284.47,210.11,250.86,285.71,279.18/
*     Denmark Strait
c      data ampl2/3.70,12.86,13.17,12.29,10.95,3.93/

*     Modul
      data ampl2/7.51,13.87,15.87,16.71,14.10,4.66/

	data ph2/214.01,256.57,273.54,275.12,188.97,106.20/


*     Middle Part of the Canadian Archipelago
c      data ampl3/1.59,4.63,3.23,2.79,2.34/

*     Modul
      data ampl3/2.26,2.18,5.52,3.27,1.93/

	data ph3/105.57,112.46,115.19,129.49,143.78/
*     Bering Strait
c      data ampl4/6.14,1.76,2.12,1.23,4.35/
      data ampl4/2.00,1.83,2.87,2.02,4.34/
	data ph4/298.31,318.92,332.64,272.39,282.34/

c     Tidal Sea level

*     Norwegian sea
      data dampl/66.63,31.66,7.85,15.19,25.45,33.03,39.36,45.40,
     *           51.08,57.04,63.83,71.61,79.61,81.23/
	data dph/134.23,107.76,64.28,319.01,311.42,3*310.20,310.72,
     *           310.90,310.93,311.17,311.63,314.19/
*     Denmark Strait
      data dampl2/27.42,24.42,25.58,44.21,89.17,120.24/
	data dph2/123.00,153.90,207.24,239.29,246.38,208.84/
*     Middle Part of the Canadian Archipelago
      data dampl3/2.00,5.24,11.02,11.35,8.97/
	data dph3/210.26,248.60,242.16,270.97,271.73/
*     Bering Strait
      data dampl4/7.69,5.60,5.60,6.53,9.17/
	data dph4/264.58,220.89,176.71,129.91,170.80/

c	dampl (:)=2.*dampl (:)
c	dampl2(:)=2.*dampl2(:)
c	dampl3(:)=2.*dampl3(:)
c	dampl4(:)=2.*dampl4(:)

	ph (:)=dph(:)  -180.
	ph2(:)=dph2(:) -180.
	ph3(:)=dph3(:) -180.
	ph4(:)=dph4(:) -180.

      INCLUDE 'Tparm.fi'

	PME     =0.
	River   =0.
	Q2Turb  =1.
	Q2Turb1 =0.
	Q2Turb2 =1.
c	CDgwd   =5.5e-3
c	Q2L     =1.e3
c	cL      =1.e3

	unept= 0.
	vnept= 0.

	bio =1.0

*     Climate restoring time (sec).
*     -----------------------------------
      tscale(1)= 365.*24.*3600.
      do k=2,kl
      tscale(k)= 36500.*24.*3600.
      end do

*     Background Horizontal turbulent coefficients
*     Mainly for numerical purposes
*    ----------------------------------
                al = 5.0e7
                alt= 1.0e5
*    ----------------------------------

*     Background Vertical turbulent coefficients
*     Mainly for numerical purposes
*    ----------------------------------
                azbg = 1.e+0
                aztbg= 1.e-2
*    ----------------------------------

	Az= azbg 
	Azt=aztbg
	Azs=aztbg

*     Vertical grid
        z(1)= 0.
        z(2)= 2.5
        z(3)= 5.
        z(4)= 7.5
        z(5)= 10.
        z(6)= 15.
        z(7)= 20.
        z(8)= 27.5
        z(9)= 37.5
        z(10)=50.
        z(11)=65.
        z(12)=82.5
        z(13)=100.
        z(14)=125.
        z(15)=150.
        z(16)=175.
        z(17)=200.
        z(18)=225.
        z(19)=250.
        z(20)=275.
        z(21)=300.
        z(22)=325.
        z(23)=350.
        z(24)=375.
        z(25)=410.
        z(26)=450.
        z(27)=500.
        z(28)=600.
        z(29)=700.
        z(30)=850.
        z(31)=1050.
        z(32)=1250.
        z(33)=1500.
        z(34)=1750.
        z(35)=2100.
        z(36)=2500.
        z(37)=3000.
        z(38)=3500.
        z(39)=4000.
        z(40)=4500.
        z(41)= 9999999.

*     -----------------
      volume=  1.1047875E+22

      call ktform(KT)

1     format(10f4.0)
1119  format(27f3.0)

      do 33 k=1,klp
33    z(k)=100.*z(k)

      do 34 k=1,kl

	if( z(k)/dzi1 .LT. 20.) then
	exps1(k)=EXP(-z(k)/dzi1)
	else
	exps1(k)=0.
	end if

	if( z(k)/dzi2 .LT. 20.) then
	exps2(k)=EXP(-z(k)/dzi2)
	else
	exps2(k)=0.
	end if

34    hz(k)=(z(k+1)-z(k))

      open (31,file='start.txt')
*     Number of years per run and number of GMRES iterations.
      read (31,*) nyear, maxite
      close(31)

      nnn=(il+2)*(jl+2)*kl*Lfort
      nn=(il+2)*(jl+2)*Lfort
      nni=nn*mgrad
      nna=nn*(mgrad+1)
      nn2= 2*Lfort*(il+4)*(jl+4)
      nd=38*il*jl*Lfort

c     Ocean Data
      open (unit=11,file=nt3dat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=12,file=udat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=13,file=ubdat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=14,file=tsdat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=16,file=dzdat,status='old',access='direct',
     *      form='unformatted',recl=nn2)
      open (unit=17,file=km2dat,status='old',access='direct',
     *      form='unformatted',recl=nn)
      open (unit=18,file=diagdat,status='old',access='direct',
     *      form='unformatted',recl=nd)


c     Methane and oxygen data

      open (unit=44,file=CH4dat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=45,file=CH4obsdat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=46,file=O2dat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=47,file=O2obsdat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=48,file='bio.dat',status='old',access='direct',
     *      form='unformatted',recl=nnn)

c     Ice Data
      open (unit=26,file=hisdat,status='old',access='direct',
     *      form='unformatted',recl=nna)
      open (unit=27,file=aidat,status='old',access='direct',
     *      form='unformatted',recl=nna)
      open (unit=28,file=uidat,status='old',access='direct',
     *      form='unformatted',recl=nn)
      open (unit=29,file=tmpisdat,status='old',access='direct',
     *      form='unformatted',recl=nna)
c      open (unit=30,file=imask,status='old',access='direct',
c     *      form='unformatted',recl=nn)
	 

c     Neptune components      
      open (unit=33,file=neptunedat,status='old',access='direct',
     *      form='unformatted',recl=nn)

                         dt=3600.*dthr
                         hx=hxgr*pi/180.
                         hy=hygr*pi/180.
                         om=.729e-4
                         r=.637e9
                         g=980.
	                   Tide_M2= (2.*pi/(12.42*3600.))

      nst=0
      do 52 j=0,jl1
      teta=pi*.5 - finord*pi/180.+(j-1.)*hy
52    si(j)=sin(teta)
      do 521 i=0,il1
      teta= alwest*pi/180. +(i-1.)*hx
521   co(i)=cos(teta)

      read (12,rec=1) u
      read (12,rec=2) v
      read (12,rec=3) w
cc      read (13,rec=1) ub
      read (14,rec=1) t
      read (14,rec=2) s

      read (44,rec=1) CH4
      read (45,rec=1) CH4obs

      read (46,rec=1) o2
      read (47,rec=1) o2obs
      read (48,rec=1) bio

      read (16,rec=1) dz
      read (16,rec=1) dzm
      read (18,rec=1) diag
      read (11,rec=1) nt3
      read (17,rec=1) km2

      read (28,rec=1) uice
      read (28,rec=2) vice

      read (26,rec=1) Hice
      read (26,rec=2) Hsnow
      read (27,rec=1) Aice
      read (29,rec=1) Tice
      read (29,rec=2) Tsnow
cccc      read (30,rec=1) ice_mask_p
	read (33,rec=1) unept
	read (33,rec=2) vnept


      close(33)
      close(13)
      close(18)
      close(11)
      close(17)
ccc      close(30)


*     Lamped matrix for scalars
      call scalar_lamped(Sclr1)


* -------------     Landfast and boundary ice     ------------------

      ice_mask_p =1

      do j=1,jl
	do i=1,il
	if(nt3(i,j,1) .eq. 0 .or. nt3(i,j,1) .ge. 2) ice_mask_p(i,j) =0
	end do
	end do

c     Narrow passages of Canadian Archipelago

      ice_mask_p(10,28)= 0
	ice_mask_p(3,16) = 0

	do j=18,22
	ice_mask_p(6,j)  = 0
	end do

c     Ice mask for fast ice
	ice_mask= ice_mask_p

* ------------    End of landfast and boundary ice    --------------

      open (32, file='time.dat')
      read (32,*) time
      close(32)

c     Global number of the day.
      day= REAL(time/(3600.d0*24.d0))
c     Number of the year.
	iybeg= INT(day/365.) + 1948

	write(*,*) 'Int.year and ## years-', iybeg, nyear

*     Baroclinic pressure I(Ro,z)
      PBclin =0.

*     Outside sea level 

      open(32,file='dzext.dat',status='old',access='direct',
     *      form='unformatted',recl=il*jl)
	read(32, rec=1) dzext_mean
	close(32)

*     Canadian Passages

	dzext(10,28)= -50.              ! Nares
	dzext(3,16:17) = dzext_mean(3,16)  -20. ! McClure
	dzext(6,:)  = dzext_mean(6,:)   -20. ! Sverdrup Basin. See the ref.:
*     Melling, H., Sea ice of the northern Canadian Arctic Archipelago, 
*     J. Geophys. Res.,
*     107(C11), 3181, doi:10.1029/2001JC001102, 2002.


	GOTO 9999

*     ---------------- Initial Tides -------------
	IF( time.EQ.0.d0) then

*     Norwegian Sea

*     Velocity
      do i=11,22
	T_new= cos(               ph(i-9)*pi/180.)-
     *       cos(Tide_M2*(-dt) +ph(i-9)*pi/180.)
	do k=1,km2(i,jl)
	ub(i,jl,k)= ub(i,jl,k)
     *+   AMPL (i-9 )*(-T_new)/(Tide_M2*dt)
	end do

*     Sea level
	T_new= cos(               dph(i-9)*pi/180.)-
     *       cos(Tide_M2*(-dt) +dph(i-9)*pi/180.)
	dzext(i,jl)= dzext(i,jl)
     *            +DAMPL (i-9 )*(-T_new)/(Tide_M2*dt)
	end do

*     Denmark Strait

*     Velocity
	do j=42,47
	T_new= cos(               ph2(j-41)*pi/180.)-
     *       cos(Tide_M2*(-dt) +ph2(j-41)*pi/180.)
	do k=1,km2(7,j)
	ub(7,j,k)= ub(7,j,k)
     *+   AMPL2(j-41)*(-T_new)/(Tide_M2*dt)
      end do
*     Sea level
	T_new= cos(               dph2(j-41)*pi/180.)-
     *       cos(Tide_M2*(-dt) +dph2(j-41)*pi/180.)
	dzext(7,j)= dzext(7,j)
     *           +DAMPL2(j-41)*(-T_new)/(Tide_M2*dt)
	end do

*     Nares Strait

*     Velocity
	T_new= cos(               193.88*pi/180.)-
     *       cos(Tide_M2*(-dt) +193.88*pi/180.)
	do k=1,km2(11,28)
	ub(11,28,k)= ub(11,28,k) +4.00*(-T_new)/(Tide_M2*dt)
	end do
*     Sea level
	T_new= cos(               105.89*pi/180.)-
     *       cos(Tide_M2*(-dt) +105.89*pi/180.)
	dzext(11,28)= dzext(11,28) +37.39*(-T_new)/(Tide_M2*dt)


*     Alaska-Canadian
	T_new= cos(               215.31*pi/180.)-
     *       cos(Tide_M2*(-dt) +215.31*pi/180.)
	do k=1,km2(4,16)
	ub(4,16:17,k)= ub(4,16:17,k) +3.06*(-T_new)/(Tide_M2*dt)
	end do
*     Sea level
	T_new= cos(               275.39*pi/180.)-
     *       cos(Tide_M2*(-dt) +275.39*pi/180.)
	dzext(4,16)= dzext(4,16) +10.17*(-T_new)/(Tide_M2*dt)


*     Middle part of the Canadian 

*     Velocity
      do j=17,22
	T_new= cos(               PH3(j-17)*pi/180.)-
     *       cos(Tide_M2*(-dt) +PH3(j-17)*pi/180.)
	do k=1,km2(6,j)
	ub(6,j,k)= ub(6,j,k) +AMPL3(j-17)*(-T_new)/(Tide_M2*dt)
	end do
*     Sea level
	T_new= cos(               DPH3(j-17)*pi/180.)-
     *       cos(Tide_M2*(-dt) +DPH3(j-17)*pi/180.)
	dzext(6,j)= dzext(6,j) +DAMPL3(j-17)*(-T_new)/(Tide_M2*dt)
	end do

*     Bering Strait

*     Velocity
      do i=11,15
	T_new= cos(               PH4(i-10)*pi/180.)-
     *       cos(Tide_M2*(-dt) +PH4(i-10)*pi/180.)
	do k=1,km2(i,1)
	ub(i,1,k)= ub(i,1,k) +AMPL4(i-10)*(-T_new)/(Tide_M2*dt)
	end do
*     Sea level
	T_new= cos(               DPH4(i-10)*pi/180.)-
     *       cos(Tide_M2*(-dt) +DPH4(i-10)*pi/180.)
	dzext(i,1)= dzext(i,1) +DAMPL4(i-10)*(-T_new)/(Tide_M2*dt)
      end do

 
	end if   ! time=0

9999  continue


      DO j=1,jl
      DO i=1,il
      KK2=KM2(i,j)
      IF( KK2 .GT. 0) THEN
*     Pressure by snow on the ice should be compensated by level
	hi= 0.
      do mg=1,mgrad
	hi=hi+rosdry*hsnow(mg,i,j)/row ! snow mass
	hi=hi+roi   *hice (mg,i,j)/row ! ice mass
	end do
      dz(i,j) = -DBLE(hi)
c	dzm(i,j)= dz(i,j)
	END IF            ! kk2>0
	end do
	end do

	PME   =0.
	PME_DZ=0.

	Tsnow=-100.
	Tice =-100.
      

      Pice=2.75e8

c     ==============================================
c                    TIME LOOP
c     ==============================================

      DO iYEAR= iybeg, iybeg+nyear-1

	idglobal = 0

	DO imonth= 1,12

      ntide=0

c     Monthly means.
	uimonth= 0.
	vimonth= 0.
	dzmonth= 0.d0
	umonth = 0.
	vmonth = 0.
	wmonth = 0.
	Tmonth = 0.
	Smonth = 0.
	himonth= 0.
	aimonth= 0.
	hsmonth= 0.
	CH4month=0.
	O2month =0.
	biomonth=0.
      Flux_CH4_To_A_month= 0.


c     Ice mask for fast ice: ice stops in October-April.
*	do j=1,jl
*      do i=1,il
*	ice_mask(i,j)= ice_mask_p(i,j)
*      IF(imonth.GE.5 .AND. imonth.LE.9) then
*      if(nt3(i,j,1).EQ.1) ice_mask(i,j)=1
*      end if
*	end do
*	end do
	  
c     Precipitation according to the time. 
*     Initial data: Yang, D., 1999:  An improved precipitation
*     climatology for the Arctic Ocean, Geophys.
*     Res. Lett., 26(11), 1625-1628. OR:
*     Serreze, M.C., M.P. Clark and D.H. Bromwich, 2003: Monitoring
*     precipitation over the Arctic terrestrial drainage system:
*     Data requirements, shortcomings and applications of 
*     atmospheric reanalysis. Journal of Hydrometeorology, 4, 387-407.
c     Precipitation in mm/month --> in cm/sec.

c      open (41,file='pr_150703.dat',access='direct',
c     *      form='unformatted',recl=il*jl)
c	read(41,rec=imonth) Pr
c	close(41)
c      do j=1,jl
c      do i=1,il
c      Pr(i,j)= Pr(i,j)/(86400.*10.*FLOAT(MEC(imonth))) * 0.5 !!0.6666
c	end do
c	end do


c     Xie, P.P., and P.A. Arkin, 1996: Analyses of global monthly precipitation
c     using gauge observations, satellite estimates, and numerical model predictions.
c     Journal of Climate, 9, 840-858.
c     Xie, P.P. and P.A. Arkin, 1997: Global Precipitation: A 17-year monthly analysis
c     based on gauge observations, satellite estimates, and numerical model outputs. 
c     Bull. Amer. Met. Soc., 78, 2539-2558. 

      open (41,file='PRXA.dat',access='direct',
     *      form='unformatted',recl=il*jl)
	read(41,rec=imonth) Pr
	close(41)

      do j=1,jl
      do i=1,il
      Pr(i,j)= Pr(i,j)/864000. ! mm/day to cm/sec
	end do
	end do


c     River as rain. Comment it otherwise.
      call 
     &river_as_rain(imonth,River,il,jl,nt3,si,cg,il1,jl1,kl,r,hx,hy)

c     Cloud cover
      open (41,file='cloud.dat',access='direct',
     *      form='unformatted',recl=il*jl)
	read(41,rec=imonth) cloud
	close(41)


C     PHC 3.0 Monthly Mean Hydrology
      call TSInt(imonth)

c     Level at open boundary
      do i=11,23
	dzext(i,jl)= 5.*cos((imonth-10.)*pi/6.)  !+25. !!*(i-11.)/12. 
     *           + 1.0*dzext_mean(i,jl)                ! Norwegian Sea
	end do

      do j=43,46
	dzext(7,j) = 5.*cos((imonth-10.)*pi/6.)  -70. !!- 30.*(j-43.)/3.
     *          +  1.0*dzext_mean(7,j)                 ! Denmark Strait
	end do

	dzext(:,1) = 5.*cos((imonth-8. )*pi/6.)
     *           + 1.0*dzext_mean(:,1)         +50.    ! Bering Passage


*     --------------------- Steps - Days -----------------------------------

	DO iday=1,MEC(imonth)
	idglobal= idglobal +1

*     --------------------- Steps in the Day -------------------------------

	DO L=1,INT(24./dthr)

ccc	ntide= ntide +1

      GOTO 8888

*     ---------------- Tides -------------
      time_tide= MOD(REAL(time),(12.42*3600.)) 

*     Norwegian Sea

*     Velocity
      do i=11,22
	T_old= cos(Tide_M2*(time_tide   ) +ph(i-9)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +ph(i-9)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +ph(i-9)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +ph(i-9)*pi/180.)

	do k=1,km2(i,jl)
	ub(i,jl,k)= ub(i,jl,k)
     *+   AMPL (i-9 )*(T_old-T_new)/(Tide_M2*dt)
	end do
*     Sea level
	T_old= cos(Tide_M2*(time_tide   ) +dph(i-9)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +dph(i-9)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +dph(i-9)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +dph(i-9)*pi/180.)
	dzext(i,jl)= dzext(i,jl)
     *            +DAMPL (i-9 )*(T_old-T_new)/(Tide_M2*dt)
	end do

*     Denmark Strait
	do j=42,47
	T_old= cos(Tide_M2*(time_tide   ) +ph2(j-41)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +ph2(j-41)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +ph2(j-41)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +ph2(j-41)*pi/180.)
	do k=1,km2(7,j)
	ub(7,j,k)= ub(7,j,k)
     *+   AMPL2(j-41)*(T_old-T_new)/(Tide_M2*dt)
	end do
	T_old= cos(Tide_M2*(time_tide   ) +dph2(j-41)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +dph2(j-41)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +dph2(j-41)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +dph2(j-41)*pi/180.)
	dzext(7,j)= dzext(7,j)
     *+   DAMPL2(j-41)*(T_old-T_new)/(Tide_M2*dt)

	end do

*     Nares Strait
	T_old= cos(Tide_M2*(time_tide   ) +193.88*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +193.88*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +193.88*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +193.88*pi/180.)
	do k=1,km2(11,28)
	ub(11,28,k)= ub(11,28,k) +4.00*(T_old-T_new)/(Tide_M2*dt)
	end do

	T_old= cos(Tide_M2*(time_tide   ) +105.89*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +105.89*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +105.89*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +105.89*pi/180.)
	dzext(11,28)= dzext(11,28) +37.39*(T_old-T_new)/(Tide_M2*dt)

*     Alaska-Canadian
	T_old= cos(Tide_M2*(time_tide   ) +215.31*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +215.31*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +215.31*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +215.31*pi/180.)
	do k=1,km2(4,16)
	ub(4,16:17,k)= ub(4,16:17,k) +3.1*(T_old-T_new)/(Tide_M2*dt)
	end do

	T_old= cos(Tide_M2*(time_tide   ) +275.39*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +275.39*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +275.39*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +275.39*pi/180.)
	dzext(4,16)= dzext(4,16) +10.17*(T_old-T_new)/(Tide_M2*dt)


*     Middle part of the Canadian 
      do j=17,22
	T_old= cos(Tide_M2*(time_tide   ) +PH3(j-17)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +PH3(j-17)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +PH3(j-17)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +PH3(j-17)*pi/180.)
	do k=1,km2(6,j)
	ub(6,j,k)= ub(6,j,k) +AMPL3(j-17)*(T_old-T_new)/(Tide_M2*dt)
	end do
	T_old= cos(Tide_M2*(time_tide   ) +DPH3(j-17)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +DPH3(j-17)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +DPH3(j-17)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +DPH3(j-17)*pi/180.)
	dzext(6,j)= dzext(6,j) +DAMPL3(j-17)*(T_old-T_new)/(Tide_M2*dt)
	end do

*     Bering Strait
      do i=11,15
	T_old= cos(Tide_M2*(time_tide   ) +PH4(i-10)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +PH4(i-10)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +PH4(i-10)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +PH4(i-10)*pi/180.)
	do k=1,km2(i,1)
	ub(i,1,k)= ub(i,1,k) +AMPL4(i-10)*(T_old-T_new)/(Tide_M2*dt)
	end do
	T_old= cos(Tide_M2*(time_tide   ) +DPH4(i-10)*pi/180.)-
     *       cos(Tide_M2*(time_tide-dt) +DPH4(i-10)*pi/180.)
	T_new= cos(Tide_M2*(time_tide+dt) +DPH4(i-10)*pi/180.)-
     *       cos(Tide_M2*(time_tide   ) +DPH4(i-10)*pi/180.)
	dzext(i,1)= dzext(i,1) +DAMPL4(i-10)*(T_old-T_new)/(Tide_M2*dt)
      end do
*     ------------------------------------

8888  continue

      nst=nst +1



c     Density in situ

      do j=1,jl
      do i=1,il
      if( km2(i,j) .GT. 0) then
      do K=1,km2(i,j)
	ppp= 1.e-5* g*row*z(k)                 
      ro(i,j,k)= 1.e-3*REAL(sigma_t(t(i,j,k),s(i,j,k), ppp))
      end do
      end if
      end do
      end do

c     Vertical convective adjustment - unstable cells adjusted by pairs
**	call TS_convection


c     Potential density.

      do j=1,jl
      do i=1,il
      if( km2(i,j) .GT. 0) then
      do K=1,km2(i,j)
      ropot(i,j,k)= 1.e-3*REAL(sigma_t(t(i,j,k),s(i,j,k), 0.))
      end do
      end if
      end do
      end do


      call GMG_Visbeck_Coeffs (Agm) ! Variable GMG coefficient on the triangle


c     Baroclinic Pressure. 
c     m=1 - no Inverted Barometer, m=2 - Inverted Barometer. 

      call denbpr(2)

      call Forcing3(iyear,idglobal,L,dthr,TA,Q2m,WX,WY,Pa,il,jl,
     #              il1,jl1,si,co,hx,hy,r,om,serv,serv1) 

      call shortwave(SW,Pa,Q2m,idglobal,L,dthr,cloud,il,jl)

      call UB_SLO(ub,unept,vnept,nt3,km2,dzext,Pext,wx,wy,
     *            tobs,sobs,az,il,jl,il1,jl1,kl,klp,
     *            aice,u,v,uice,vice,CDgwd,mgrad,
     *            hx,hy,r,si,co,hz,z,g,om,row,roa)

c     ---------------------------------------------------------


c     Turbulence Closer Scheme Kantha and Clayson, 1994    


      Sclr2= Q2turb
	Sclr1= 0.
	SclrObs=Q2turb 
	Qcascad=0.
	call advect3D(Sclr,Sclr1,Sclr2,SclrObs,Qcascad,Agm)

*     Limitation Q2 > 0
	do j=1,jl
	do i=1,il
	kb=km2(i,j)
	if(kb.GT.0) then
	do k=1,kb
	Q2turb(i,j,k)= max( Sclr(i,j,k), 0.)
	end do
	end if
	end do
	end do

      call mixerKC


c	do j=1,jl
c	do i=1,il
c	do k=1,4
c	if(nt3(i,j,k) .NE. 0) then
c	az (i,j,k)= az (i,j,k) +1.e3
c	azt(i,j,k)= azt(i,j,k) +1.e3
c	azs(i,j,k)= azs(i,j,k) +1.e3
c	end if
c	end do
c	end do
c	end do

c     ------ Correction in river estuaries ---------------------

c     Ob and Yenisey, Puhr, Taymyra, Pyasina
c      i=34
c	do j=27,29
c	do k=1,min(4,km2(i,j)-1)    ! mixed upper 10m
c	AZS(i,j,k)= AZS(i,j,k)+ 1.e3
c	end do
c	end do

c     Lena, Khatanga, Olenek, Yana
c      i= 31
c	do j=11,14
c	do k=1,min(4,km2(i,j)-1)    ! mixed upper 10m
c	AZS(i,j,k)= AZS(i,j,k)+ 1.e3
c	end do
c	end do

c     Mackenzie
c      i=2
c	j=10
c	do k=1,min(4,km2(i,j)-1)    ! mixed upper 10m
c	AZS(i,j,k)= AZS(i,j,k)+ 1.e3
c	end do

c     Dvina and Mezen
c      i=34
c	j=43
c	do k=1,min(4,km2(i,j)-1)    ! mixed upper 10m
c	AZS(i,j,k)= AZS(i,j,k)+ 1.e3
c	end do

c     Pechora
c      i=35
c	j=36
c	do k=1,min(4,km2(i,j)-1)    ! mixed upper 10m
c	AZS(i,j,k)= AZS(i,j,k)+ 1.e3
c	end do

c     Kolyma+Indigirka
c      i=24
c	j=6
c	do k=1,min(4,km2(i,j)-1)    ! mixed upper 10m
c	AZS(i,j,k)= AZS(i,j,k)+ 1.e3
c	end do

*     ============== Sea Ice Dynamics ==========================

	uice2= uice
	vice2= vice

      call ice_dyn_evp(INT(dt/30.))  !!! 30 seconds internal time step

	uice2= uice
	vice2= vice
	call ice_inertion(1)
	uice1= uice
	vice1= vice
	call ice_inertion(2)

c     Ice/Snow Advection  

*     Divergence - used for sea ice drift.

      div_ice_tr= 0.0
	do j=1,jl
      do i=1,il
	if(km2(i,j).GT.0) then
      do mg=0,mgrad
      div_ice_tr(i,j)= div_ice_tr(i,j)+aice(mg,i,j)
      end do 
	end if
	end do
	end do

      Hice2 =Hice 
      Hsnow2=Hsnow
      Aice2 =Aice 
 
      call iceadvect(Aice,Aice1,Aice2,uice,vice,
     *RSice,ice_mask,nt3,si,KT,il,jl,il1,jl1,kl,0,mgrad,r,hx,hy,dt)

      call iceadvect(Hice,Hice1,Hice2,uice,vice,
     *RSice,ice_mask,nt3,si,KT,il,jl,il1,jl1,kl,1,mgrad,r,hx,hy,dt)

      call iceadvect(Hsnow,Hsnow1,Hsnow2,uice,vice,
     *RSice,ice_mask,nt3,si,KT,il,jl,il1,jl1,kl,1,mgrad,r,hx,hy,dt)


c     Small & negative values cut-off

	do mg=1,mgrad
	do j=1,jl
	do i=1,il
	if(km2(i,j).GT.0) then

      if(Hsnow(mg,i,j) .LT. Hsmin) then
	Hsnow(mg,i,j)= 0.
	Tsnow(mg,i,j)= 0.
	end if

      if(Aice(mg,i,j).LT.aimin .OR. Hice(mg,i,j).LT.Himin) then
	Aice (mg,i,j) =0.
	Hice (mg,i,j) =0.
	Tice (mg,i,j) =0.
	Hsnow(mg,i,j) =0.
	Tsnow(mg,i,j) =0.
	end if
	 
	end if  ! km2>0
	end do  ! i
	end do  ! j
	end do  ! mg

*     Divergence - used for sea ice drift.
	do j=1,jl
      do i=1,il
	if(km2(i,j).GT.0) then
      do mg=0,mgrad
      div_ice_tr(i,j)= div_ice_tr(i,j)-aice(mg,i,j)
      end do 
	div_ice_tr(i,j)= div_ice_tr(i,j)/dt
	end if
	end do
	end do


c     Ice categories redistribution
      Call Rebin

c     Ridging
      Call Redis

c     Ice categories redistribution
      Call Rebin


*     ============= End of Sea Ice Dynamics ===================



c     =========== Vertical T,S diffusion and Ice Thermodynamics
      call tsice

c      write(*,*) 'TS ice'
c      call varts

c     Ice categories redistribution
      Call Rebin

*     ================ T, S, transport ========================

c     Temperature and Salinity Advection/Diffusion
c      tm2= t
c      sm2= s
c	call advect3D(T,Tm1,Tm2,Tobs,Agm)
c	call advect3D(S,Sm1,Sm2,Sobs,Agm)

      call Cascading

      Sclr2= T
	Sclr1= 0.
	SclrObs=Tobs
	Qcascad = QTc 
	call advect3D(Sclr,Sclr1,Sclr2,SclrObs,Qcascad,Agm)
	T= Sclr


      Sclr2= S
	Sclr1= 0.
	SclrObs=Sobs 
	Qcascad=QSc
	call advect3D(Sclr,Sclr1,Sclr2,SclrObs,Qcascad,Agm)
	S= Sclr

*     Temperature limitation T > Tfr
c	do j=1,jl
c	do i=1,il
c	kb=km2(i,j)
c	if(kb.GT.0) then
c	do k=1,kb
c	T(i,j,k)= max( T(i,j,k), Tfr(S(i,j,k),0.))
c	end do
c	end if
c	end do
c	end do

c      call sadvect  
c      call tadvect  

c      if( MOD(nst,lprint) .EQ. 0) then
*      write(*,*) 'TS advect'
*      call varts
c      end if


c     Climate source effect, if any  --- 11 years of spin-up
      IF(NRestor .EQ. 1) then
      do j=1,jl
      do i=1,il
      kb=km2(i,j)
      if(kb.GT.0) then
      do k=1,kb
      Coef=dt/tscale(k)
cc	if(j.GE.47) Coef= 0.1
      Det= 1.+coef
cc      T(i,j,k)=(T(i,j,k)+coef*Tobs(i,j,k))/Det
      S(i,j,k)=(S(i,j,k)+coef*Sobs(i,j,k))/Det
      end do
      end if
      end do
      end do
	END IF

*     salinity limitation S > 0
	do j=1,jl
	do i=1,il
	kb=km2(i,j)
	if(kb.GT.0) then
	do k=1,kb
	S(i,j,k)= max( S(i,j,k), 0.)
	end do
	end if
	end do
	end do


*     ============== End of T, S, transport ====================


c     ================== Momentum transport ====================

      um2= u
      vm2= v
      call vhdif(1)
c     Vertical Water/Ice Momentum
***      call vvdif
      um1= u
      vm1= v
      call vhdif(2)


*     ------------------- NEPTUNE ------------------------------
                       call neptune
*     ----------------------------------------------------------

c      write(*,*) 'U advect'
c      call maxu


c     =================== Vertical Momentum ====================
      call vvdif

c      write(*,*) 'U Vdif'
c      call maxu



*     =============== End of momentum transport ================


c     ============ Velocity-Pressure Adjustment (No Level)
      um2= u
      vm2= v
      um1= u
      vm1= v
      call vadapt


c      write(*,*) 'U adapt'
c      call maxu

c     ============ Sea-Level Elevation Problem
      dzm= dz
      um1= u
      vm1= v
      call rside(ncond, rs)

      call dzcal4(dz,rs,diag,km2,adz,il,jl,il1,jl1,ilp,jlp,maxite,
     *            len,krylov,krp1,b,sl,r0,sk,ck,giv,w0,z0,v0)

*     =====  Ocean Velocity-Pressure Adjustment (Due to Level)
      call dzvel

c      write(*,*) 'U + grad Dz'
c      call maxu


*     ====== Vertical Velocity.
      call wcals

*     ++++++++++++++++++ Methane evolution +++++++++++++++++++++

**      call CH4_Atlantic_Boundary  ! Boundary conditions for methane 
	                            ! at Norwegian Sea
**      call O2_Boundary            ! Oxygen at open boundaries

cc      write(*,*) 'initial', CH4(3,16,15)

c     River inflow input
**      call CH4Rivers (imonth,CH4,nt3,km2,si,cg,il1,jl1,kl,r,hx,hy,hz,dt)

c     CH4 and O2 Advection/Diffusion


c      CH4m2= CH4
c	call advect3D(CH4,CH4m1,CH4m2,CH4obs,Agm)
c      O2m2= O2
c	call advect3D(O2,O2m1,O2m2,O2obs,Agm)

**      Sclr2= CH4
**	Sclr1= 0.
**	SclrObs=CH4Obs 
**	call advect3D(Sclr,Sclr1,Sclr2,SclrObs,Agm)
**	CH4= Sclr

**      Sclr2= O2
**	Sclr1= 0.
**	SclrObs=O2obs 
**	call advect3D(Sclr,Sclr1,Sclr2,SclrObs,Agm)
**	O2= Sclr


**      Sclr2= bio
**	Sclr1= 0.
**	SclrObs=1.0 
**	call advect3D(Sclr,Sclr1,Sclr2,SclrObs,Agm)
**	bio  = Sclr


cc      write(*,*) 'transport', CH4(3,16,15)

c     Vertical CH4 and O2 fluxes
*      call VertCH4
*     call VertO2
* 	call VertBIO

cc      write(*,*) 'vertical', CH4(3,16,15)

c     Methane oxidation
*      call CH4Oxidation

ccc	call bio_recovery  !!! Biology recovery

cc      write(*,*) 'oxidation', CH4(3,16,15)

*     methane limitation > 0
c	do j=1,jl
c	do i=1,il
c	kb=km2(i,j)
c	if(kb.GT.0) then
c	do k=1,kb
c	CH4(i,j,k)= max( CH4(i,j,k), 0.)
c	O2 (i,j,k)= max( O2 (i,j,k), 0.)
c	end do
c	end if
c	end do
c	end do


*     ++++++++++++++++++ End of Methane step +++++++++++++++++++


      time= time + DBLE(dt)


ccc      call tr(iYear,iMonth,u,v,km2,si,hz,il1,jl1,kl,0)


      if( MOD(nst,lprint) .EQ. 0) then
	write(*,*)'BG S model and data, 5m =',s(12,13,3),sobs(12,13,3)
      call maxu
      call dznorm(dz,dzm,nt3,ilp,jlp,il1,jl1,kl,cg,si)
      call dzerr(w,dz,dzm,PME,diag,rs,nt3,KT,cg,si,
     *              ilp,jlp,il1,jl1,il,jl,kl,dt)
      call CH4statistics
      call maxH
	write(*,*) 'Year',iYear,' month',imonth,' day', iday
      write(*,*)'-----------------------------------'
      end if

c     Monthly means.
	do j=0,jl1
      do i=0,il1

	uimonth(i,j)= uimonth(i,j) +uice(i,j)
	vimonth(i,j)= vimonth(i,j) +vice(i,j)
	dzmonth(i,j)= dzmonth(i,j) +dz  (i,j)


      Flux_CH4_To_A_month(i,j)= 
     *Flux_CH4_To_A_month(i,j) +Flux_CH4_To_A(i,j)


	do k=1,kl
	umonth(i,j,k)= umonth(i,j,k) +u(i,j,k)
	vmonth(i,j,k)= vmonth(i,j,k) +v(i,j,k)
	wmonth(i,j,k)= wmonth(i,j,k) +w(i,j,k)
	Tmonth(i,j,k)= Tmonth(i,j,k) +T(i,j,k)
	Smonth(i,j,k)= Smonth(i,j,k) +S(i,j,k)
	CH4month(i,j,k)= CH4month(i,j,k) +CH4(i,j,k)
	O2month(i,j,k) = O2month(i,j,k)  +O2(i,j,k)
	biomonth(i,j,k)= biomonth(i,j,k) +bio(i,j,k)
	end do            ! K loop.

	do mg=0,mgrad
	himonth(mg,i,j)= himonth(mg,i,j) +hice(mg,i,j)
	aimonth(mg,i,j)= aimonth(mg,i,j) +aice(mg,i,j)
	hsmonth(mg,i,j)= hsmonth(mg,i,j) +hsnow(mg,i,j)
	end do            ! Ice thickness gradations.

	end do
	end do

	END DO  ! steps
	END DO  ! days

	cc= 24.*FLOAT(MEC(imonth))/dthr

c     Monthly means.
	do j=0,jl1
      do i=0,il1

	uimonth(i,j)= uimonth(i,j)/cc
	vimonth(i,j)= vimonth(i,j)/cc
	dzmonth(i,j)= dzmonth(i,j)/DBLE(cc)
      Flux_CH4_To_A_month(i,j)= 
     *Flux_CH4_To_A_month(i,j) /cc
	do k=1,kl
	umonth(i,j,k)= umonth(i,j,k)/cc
	vmonth(i,j,k)= vmonth(i,j,k)/cc
	wmonth(i,j,k)= wmonth(i,j,k)/cc
	Tmonth(i,j,k)= Tmonth(i,j,k)/cc
	Smonth(i,j,k)= Smonth(i,j,k)/cc
	CH4month(i,j,k)= CH4month(i,j,k)/cc
	biomonth(i,j,k)= biomonth(i,j,k)/cc
	O2month(i,j,k)= O2month(i,j,k)*22.4e-6/cc  ! nmol\L -> mL/L
	end do            ! K loop.

	do mg=0,mgrad
	himonth(mg,i,j)= himonth(mg,i,j)/cc
	aimonth(mg,i,j)= aimonth(mg,i,j)/cc
	hsmonth(mg,i,j)= hsmonth(mg,i,j)/cc
	end do            ! Ice thikness gradations.

	end do
	end do

*     Sea level normalization.
c      sum=0.
c      area=0.
c      do j=1,jl
c      s0=si(j)
c      do i=1,il
c      if( nt3(i,j,1) .NE. 0) then
c      area= area + cg(abs(nt3(i,j,1)))*s0
c      sum = sum  + cg(abs(nt3(i,j,1)))*s0*REAL(dzmonth(i,j))
c      end if
c      end do
c      end do
c      dzmean= sum/area
c      do j=1,jl
c      do i=1,il
c      if(nt3(i,j,1) .NE. 0) then
c	dzmonth(i,j) = dzmonth(i,j) -DBLE(dzmean)
c	end if
c      end do
c      end do

*     Writing monthly mean data on disc.
      call Fdisc(umonth,vmonth,wmonth,tmonth,smonth,dzmonth,az,azt,
     *           uimonth,vimonth,himonth,aimonth,hsmonth,CH4month,
     *           Flux_CH4_To_A_month,O2month,biomonth,
     *           iYear,iMonth,il,jl,il1,jl1,ilp,jlp,kl,mgrad)
      call IArea(iYear,iMonth,Aimonth,nt3,il1,jl1,kl,
     *           mgrad,cg,si,R,hx,hy)
      call tr(iYear,iMonth,umonth,vmonth,km2,si,hz,il1,jl1,kl,1)
      call energ(iYear,iMonth,umonth,vmonth)

	END DO  ! months
	END DO  ! years
c     ================================================

      write(12,rec=1) u
      write(12,rec=2) v
      write(12,rec=3) w
      close(12)
      write(16,rec=1) dz
      close(16)
      write(14,rec=1) t
      write(14,rec=2) s
      close(14)
      write (28,rec=1) uice
      write (28,rec=2) vice
      close(28)

      write (26,rec=1) Hice
      write (26,rec=2) Hsnow
      close(26)
      write (27,rec=1) Aice
      close(27)
      write (29,rec=1) Tice
      write (29,rec=2) Tsnow
      close(29)

      open (32,file='time.dat')
      write(32,*) time
      close(32)

      write (44,rec=1) CH4
	close(44)
	close(45)

      write (46,rec=1) O2
	close(46)
	close(47)

      write (48,rec=1) bio
	close(48)


      stop
      end

      subroutine ktform(KT)
      real KT(6,13)
      do ntr=1,6
      do n  =1,13
      KT(ntr,n)=0.
      end do
      end do

      do ntr=1,6
      KT(ntr,1)=1.
      end do

      KT(1,2)=1.
      KT(2,2)=1.
      KT(3,2)=1.
      KT(6,2)=1.

      KT(1,3)=1.
      KT(2,3)=1.
      KT(4,3)=1.
      KT(5,3)=1.
      KT(6,3)=1.

      KT(3,4)=1.
      KT(4,4)=1.
      KT(5,4)=1.
      KT(6,4)=1.

      KT(1,5)=1.
      KT(2,5)=1.
      KT(3,5)=1.
      KT(4,5)=1.
      KT(5,5)=1.

      KT(1,6)=1.
      KT(2,6)=1.
      KT(6,6)=1.

      KT(4,7)=1.
      KT(5,7)=1.
      KT(6,7)=1.

      KT(3,8)=1.
      KT(4,8)=1.
      KT(5,8)=1.

      KT(1,9)=1.
      KT(2,9)=1.
      KT(3,9)=1.

      KT(1,10)=1.
      KT(2,10)=1.

      KT(3,11)=1.

      KT(4,12)=1.
      KT(5,12)=1.

      KT(6,13)=1.
      return
      end

      subroutine energ (iYear,iMonth,umonth,vmonth)
*     version 23.08.2012

*     Numbering begins with the 1948 year!

*     Watcon and NDP Fortrans Lfort=4, Compaq and UNIX - Lfort=1
	Parameter (Lfort=1)

      INCLUDE 'Slo2.fi'
      dimension umonth(0:il1,0:jl1,kl), vmonth(0:il1,0:jl1,kl) 
      do k10=1,klp
      ek(k10)=0.
      ss(k10)=0.
      end do

      kkk= 12*(iYear-1948) +iMonth

      rrhx= .5*r*r*hx*hy
      volume= 0.

      do k=1,kl
      a=0.
	ss(k)= 0.

      if( k.eq.1) then
        hzk1=0.
        hzk =hz(1)
      else
        hzk1=hz(k-1)
         if( k .eq. kl) then
         hzk=0.
         else
         hzk=hz(k)
         end if
      end if

      do j=1,jl
         s0=si(j)*rrhx
      do i=1,il

      kk2= km2(i,j)
      if( k .le. kk2) then
         n=abs(nt3(i,j,k))
         np=n
         if(k.lt.kk2) np=abs(nt3(i,j,k+1))
         geom= s0*(cg(n)*hzk1+cg(np)*hzk)/6.
         ss(k)=ss(k)+geom*(umonth(i,j,k)**2 +vmonth(i,j,k)**2)
         a= a +geom
      end if
      end do                  ! i loop
	end do                  ! j loop

	volume= volume +a

        ek(klp)=ek(klp)+ss(k)
        if( a.GT.0.) then
        ek(k)=ss(k)/a
        end if
      end do                  ! K loop

      ek(klp)=ek(klp)/volume

      kn= Lfort*klp
      open (unit=50,file='e.dat',status='old',access='direct',
     *      form='unformatted',recl=kn)
      write(50,rec=kkk) ek
      close(50)

      write(*,88) ek(1),ek(4),ek(7),ek(10),ek(16),ek(klp),iYear,iMonth
88    format(' e(1)=',f6.2,' e(4)=',f6.2,' e(7)=',f6.2,
     *' e(10)=',f6.2,' e(16)=',f6.2,' eInt=',f6.2, /' Time=',i4,
     *' years ',i2,' month')
      return
      end

      subroutine maxu
      include 'Slo2.fi'
      vmax= 0.
      do j=1,jl
      do i=1,il
      kmin=MIN(km2(i,j),km2(i+1,j),km2(i,j+1),km2(i+1,j+1))
        if( kmin .GT. 0) then
          do k=1,kmin
           u1=(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k)) /3.
           v1=(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k)) /3.
           u2=(u(i+1,j+1,k)+u(i+1,j,k)+u(i,j+1,k))/3.
           v2=(v(i+1,j+1,k)+v(i+1,j,k)+v(i,j+1,k))/3.
             q1= sqrt(u1**2 +v1**2)
             q2= sqrt(u2**2 +v2**2)
               if( MAX(q1,q2) .GT. vmax) then
                  vmax= MAX(q1,q2)
                  i0= i
                  j0= j
                  k0= k
               end if
           end do
        end if
      end do
      end do
      write(*,'(6H Vmax:,F7.2,4H at ,I2,2H, ,I2,2H, ,I2)')
     *      vmax,i0,j0,k0

      vmax= 0.
      do j=1,jl
      do i=1,il
      kmin=MIN(km2(i,j),km2(i+1,j),km2(i,j+1),km2(i+1,j+1))
	ai_mx=MAX(aice(0,i,j),aice(0,i+1,j),
     *          aice(0,i,j+1),aice(0,i+1,j+1))
        if(kmin.GT.0 .AND. ai_mx.LT.1.) then
           u1=(uice(i,j)+uice(i+1,j)+uice(i,j+1)) /3.
           v1=(vice(i,j)+vice(i+1,j)+vice(i,j+1)) /3.
           u2=(uice(i+1,j+1)+uice(i+1,j)+uice(i,j+1))/3.
           v2=(vice(i+1,j+1)+vice(i+1,j)+vice(i,j+1))/3.
             q1= sqrt(u1**2 +v1**2)
             q2= sqrt(u2**2 +v2**2)
               if( MAX(q1,q2) .GT. vmax) then
                  vmax= MAX(q1,q2)
                  i0= i
                  j0= j
               end if
        end if
      end do
      end do
      write(*,'(10H Vice max:,F7.2,4H at ,I2,2H, ,I2)')
     *      vmax,i0,j0

      return
      end

      subroutine maxH
      include 'Slo2.fi'
      Himax= -1000.
      Hsmax= -1000.
      Himean=0.
	HiPole=0.
	HsPole=0.
      Hsmean=0.
      Amean= 0.
	APole=0.

	do j=1,jl
      do i=1,il
      Hi=0.
      Hs=0.
      Ai=0.
      if( km2(i,j) .GT. 0) then

      do 1 m=1,mgrad
      Hi= Hi +Hice(m,i,j)
      Hs= Hs +Hsnow(m,i,j)
      Ai= Ai+Aice(m,i,j)
      coef=cg(abs(nt3(i,j,1)))/6.
		if(i.eq.17 .and. j.eq.24) then
		HiPole= HiPole +Hice(m,i,j)
		HsPole= HsPole +Hsnow(m,i,j)
		APole=  APole  +Aice(m,i,j)
		end if
      Himean= Himean +coef*Hice(m,i,j)
      Hsmean= Hsmean +coef*Hsnow(m,i,j)
1     Amean= Amean   +coef*Aice(m,i,j)
 

c     Mean Ice Thickness in the point.
      if (Ai .GT. 0.10) then
      Hi=Hi/Ai
      Hs=Hs/Ai
      else
      Hi=0.
      Hs=0.
      end if

      if(Hi .GT. Himax) then
      Himax= Hi
      i0imax= i
      j0imax= j
      end if

      if(Hs .GT. Hsmax) then
      Hsmax= Hs
      i0smax= i
      j0smax= j
      end if

      end if
      end do
      end do
      write(*,'(7H HImax:,F7.2,4H at ,I2,2H, ,I2)')
     *      himax,i0imax,j0imax
cc      write(*,'(8f10.1)') (Hice(MM,i0imax,j0imax), MM=1,mgrad)
cc      write(*,'(8f10.1)') (Aice(MM,i0imax,j0imax), MM=1,mgrad)
      write(*,'(7H HSmax:,F7.2,4H at ,I2,2H, ,I2)')
     *      hsmax,i0smax,j0smax

      if(Amean .GT. 1.e-3) then
      Himean= Himean/Amean
      Hsmean= Hsmean/Amean
      write(*,'(18H Hice, Hsnow mean:,2F7.2)') himean,hsmean
      end if
      if(APole .GT. 0.10) then
      HiPole= HiPole/APole
      HsPole= HsPole/APole
      write(*,*)'Hice, Hsnow at North Pole', hipole, hspole
      end if

      return
      end

      subroutine dzerr(w,dz,dzm,PME,diag,rs,nt3,KT,cg,si,
     *                          ilp,jlp,il1,jl1,il,jl,kl,dt)

*     version 11.05.2012


      dimension dz(-1:ilp,-1:jlp), diag(19,il,jl), rs(il,jl)
      dimension dzm(-1:ilp,-1:jlp), w(0:il1,0:jl1,kl),si(0:jl1)
	dimension PME(il,jl),nt3(0:il1,0:jl1,kl),cg(13),KT(6,13)
      real*8 dz, dzm, diag, rs, delta, dzold,dzmean,dzmmean
	real KT

	c3=1./3.

      delta= 0.d0
      do 4 j=1,jl
      do 4 i=1,il
      if(diag(10,i,j) .GT. 0.d0) then
      dzold=
     * (diag(1 ,i,j)*dz(i  ,j-2) + diag(2 ,i,j)*dz(i+1,j-2)+
     *  diag(3 ,i,j)*dz(i+2,j-2) + diag(4 ,i,j)*dz(i-1,j-1)+
     *  diag(5 ,i,j)*dz(i  ,j-1) + diag(6 ,i,j)*dz(i+1,j-1)+
     *  diag(7 ,i,j)*dz(i+2,j-1) + diag(8 ,i,j)*dz(i-2,j  )+
     *  diag(9 ,i,j)*dz(i-1,j  ) + diag(10,i,j)*dz(i  ,j  )+
     *  diag(11,i,j)*dz(i+1,j  ) + diag(12,i,j)*dz(i+2,j  )+
     *  diag(13,i,j)*dz(i-2,j+1) + diag(14,i,j)*dz(i-1,j+1)+
     *  diag(15,i,j)*dz(i  ,j+1) + diag(16,i,j)*dz(i+1,j+1)+
     *  diag(17,i,j)*dz(i-2,j+2) + diag(18,i,j)*dz(i-1,j+2)+
     *  diag(19,i,j)*dz(i  ,j+2) - rs(i,j) )
            IF( DABS(dzold) .GT. delta) THEN
            i0= i
            j0= j
            delta= DABS(dzold)
            END IF

      end if
4     continue

100   write(*,101) delta,i0,j0
101   format(' Level Misclosure=',D12.3,' i=',i3,', j=',i3)

      werrmax= -1.e+10
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
      if( diag(10,i,j) .GT. 0.d0) then
	n=abs(nt3(i,j,1))

	DzMean=( KT(1,n)*(2.d0*dz(i,j)+dz(i+1,j)+dz(i+1,j-1))
     +        +KT(2,n)*(2.d0*dz(i,j)+dz(i+1,j-1)+dz(i,j-1))
     +        +KT(3,n)*(2.d0*dz(i,j)+dz(i,j-1)+dz(i-1,j))
     +        +KT(4,n)*(2.d0*dz(i,j)+dz(i-1,j)+dz(i-1,j+1))
     +        +KT(5,n)*(2.d0*dz(i,j)+dz(i,j+1)+dz(i-1,j+1))
     +        +KT(6,n)*(2.d0*dz(i,j)+dz(i+1,j)+dz(i,j+1))
     *                                                 )/(4.*cg(n))

	DzmMean=( KT(1,n)*(2.d0*dzm(i,j)+dzm(i+1,j)+dzm(i+1,j-1))
     +         +KT(2,n)*(2.d0*dzm(i,j)+dzm(i+1,j-1)+dzm(i,j-1))
     +         +KT(3,n)*(2.d0*dzm(i,j)+dzm(i,j-1)+dzm(i-1,j))
     +         +KT(4,n)*(2.d0*dzm(i,j)+dzm(i-1,j)+dzm(i-1,j+1))
     +         +KT(5,n)*(2.d0*dzm(i,j)+dzm(i,j+1)+dzm(i-1,j+1))
     +         +KT(6,n)*(2.d0*dzm(i,j)+dzm(i+1,j)+dzm(i,j+1))
     *                                                 )/(4.*cg(n))

c	DzMean=( KT(1,n)*S1*(dz(i,j)+dz(i+1,j)+dz(i+1,j-1))
c     +        +KT(2,n)*S2*(dz(i,j)+dz(i+1,j-1)+dz(i,j-1))
c     +        +KT(3,n)*S3*(dz(i,j)+dz(i,j-1)+dz(i-1,j))
c     +        +KT(4,n)*S4*(dz(i,j)+dz(i-1,j)+dz(i-1,j+1))
c     +        +KT(5,n)*S5*(dz(i,j)+dz(i,j+1)+dz(i-1,j+1))
c     +        +KT(6,n)*S6*(dz(i,j)+dz(i+1,j)+dz(i,j+1))
c     *                                                 )/(3.*cg(n))

c	DzmMean=( KT(1,n)*S1*(dzm(i,j)+dzm(i+1,j)+dzm(i+1,j-1))
c     +         +KT(2,n)*S2*(dzm(i,j)+dzm(i+1,j-1)+dzm(i,j-1))
c     +         +KT(3,n)*S3*(dzm(i,j)+dzm(i,j-1)+dzm(i-1,j))
c     +         +KT(4,n)*S4*(dzm(i,j)+dzm(i-1,j)+dzm(i-1,j+1))
c     +         +KT(5,n)*S5*(dzm(i,j)+dzm(i,j+1)+dzm(i-1,j+1))
c     +         +KT(6,n)*S6*(dzm(i,j)+dzm(i+1,j)+dzm(i,j+1))
c     *                                                 )/(3.*cg(n))

      wmod= -REAL(dzmean-dzmmean)/dt
      werr= ABS(w(i,j,1)-wmod -PME(i,j)/dt)
      if (werr.GT.werrmax) then
      i0= i
      j0= j
      werrmax= werr
	wmodmax= wmod
      end if
      end if
      end do
      end do
      write(*,1) werrmax,wmodmax,w(i0,j0,1),i0,j0
	write(*,*) 'PME:', PME(i0,j0)/dt
1     format(' W Miscl: ',e9.3,' -dz/dt: ',e9.3,' w: ',e9.3,' i=',i2,
     #       ' j=',i2)
      return
      end

      subroutine dznorm(dz,dzm,nt3,ilp,jlp,il1,jl1,kl,cg,si)
      dimension dz(-1:ilp,-1:jlp),dzm(-1:ilp,-1:jlp),nt3(0:il1,0:jl1,kl)
      dimension cg(13), si(0:jl1)

      real*8 dz, dzm, area, sum, dzmean
      sum=0.d0
      area=0.d0
      do j=1,jl1-1
      s0=si(j)
      do i=1,il1-1
      if( nt3(i,j,1) .NE. 0) then
      area= area + DBLE(cg(abs(nt3(i,j,1)))*s0)
      sum = sum  + DBLE(cg(abs(nt3(i,j,1)))*s0)*dz (i,j)
      end if
      end do
      end do
      dzmean  = sum  /area
c      do i=1,il1-1
c      do j=1,jl1-1
c      if(nt3(i,j,1) .NE. 0) then
c     	dz (i,j) = dz (i,j) -dzmean
c     	dzm(i,j) = dzm(i,j) -dzmean
c      end if
c      end do
c      end do
      write(*,22) dzmean
22    format(1x,7hdzmean=,e9.3)
      return
      end

      subroutine IArea(iYear,iMonth,Aice,nt3,il1,jl1,kl,
     *                 mgrad,cg,si,R,hx,hy)
c     Total Ice Extent Calculation & Writing to Disc.
c     Numbering begins with the 1948 year.
c     Watcon and NDP Fortrans Lfort=1, Compaq and UNIX - Lfort=4
c     Version 04.04.04.
      Parameter (Lfort=1)
      dimension Aice(0:mgrad,0:il1,0:jl1),nt3(0:il1,0:jl1,kl),
     *          cg(13),si(0:jl1)

      Cmetric=hx*hy*R*R*1.e-10
      area= 0.
      ext = 0.

      do j=1,jl1-1
      s0=si(j)
      do i=1,il1-1
	mark= ABS(nt3(i,j,1))
      if( mark.GT.0) then

      asum=0.
      do m=1,mgrad
      asum=asum+Aice(m,i,j)
      end do

c     Only ice with compactness > 0.15 - Barry et. al., 1983.???
c     Data on compactness shows us values of 0.1
        if(asum .GT. 0.10) then 
	  ext= ext+s0*cg(mark)
	  area= area + s0*cg(mark)*asum
	  end if

      end if  ! mark>0
      end do  ! i
      end do  ! j

      area = 0.166667*Cmetric*area
      ext  = 0.166667*Cmetric*ext
      write(*,*)' Total Ice Area (km2):', area
      write(*,*)' Total Ice Extent(km2):', ext

*     Writing data on disc.

      kkk= 12*(iYear-1948) +iMonth

      open (unit=50,file='iarea.dat',status='old',access='direct',
     *      form='unformatted',recl=Lfort)
      write(50,rec=kkk) area
      close(50)
      open (unit=50,file='iext.dat',status='old',access='direct',
     *      form='unformatted',recl=Lfort)
      write(50,rec=kkk) ext
      close(50)
      return
      end

      subroutine Fdisc(u,v,w,t,s,dz,az,azt,uice,vice,hice,aice,hsnow,
     *                 CH4,Flux_CH4_To_A_month,O2,bio,
     *                 iYear,iMonth,il,jl,il1,jl1,ilp,jlp,kl,mgrad)

*     Monthly means or any other information if needed
*     version 17.02.2014 for big arrays methane + oxygen case

*     Watcon and NDP Fortrans Lfort=1, Compaq and UNIX - Lfort=4
      Parameter (Lfort=1)

*     The most important information writing.
*     File name format is: Fmonth.dat

*
*     Water parameters

      dimension u(0:il1,0:jl1,kl), v(0:il1,0:jl1,kl),
     *          t(0:il1,0:jl1,kl), s(0:il1,0:jl1,kl),
     *          CH4(0:il1,0:jl1,kl),Flux_CH4_To_A_month(0:il1,0:jl1),
     *          O2 (0:il1,0:jl1,kl), bio (0:il1,0:jl1,kl),
     *          dz(-1:ilp,-1:jlp), w(0:il1,0:jl1,kl),
     *          az(il,jl,kl),azt(il,jl,kl)

*     Ice/Snow parameters

      dimension uice(0:il1,0:jl1), vice(0:il1,0:jl1),
     *          HIce (0:mgrad,0:il1,0:jl1),Aice(0:mgrad,0:il1,0:jl1),
     *          HSnow(0:mgrad,0:il1,0:jl1)

      double precision dz
      character*72 name
      nnn=(il1+1)*(jl1+1)*kl*Lfort
      nnf=(il1+1)*(jl1+1)*Lfort
      nn2= Lfort*(il1+3)*(jl1+3)*2
      nn = Lfort*il*jl*kl
      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)


	nrec= 12*(iYear-1948) + 3*(iMonth-1) +1

      name='uomonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= nrec)   u
      write(31,rec= nrec+1) v
      write(31,rec= nrec+2) w
      close(31)

	nrec= 12*(iYear-1948) + 2*(iMonth-1) +1

      name='uimonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nuice)
      write(31,rec=nrec)   uice
      write(31,rec=nrec+1) vice
      close(31)

	nrec= 12*(iYear-1948) + 3*(iMonth-1) +1

      name='himonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nhice)
      write(31,rec= nrec)   Hice
      write(31,rec= nrec+1) Aice
      write(31,rec= nrec+2) Hsnow
      close(31)

	nrec= 12*(iYear-1948) + 2*(iMonth-1) +1

      name='tsmonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= nrec)   t
      write(31,rec= nrec+1) s
      close(31)


	nrec= 12*(iYear-1948) + iMonth

      name='ch4month.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= nrec) CH4
      close(31)

      name='o2month.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= nrec) O2
      close(31)

      name='biomonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= nrec) bio
      close(31)

      name='fluxch4month.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnf)
      write(31,rec= nrec) Flux_CH4_To_A_month
      close(31)

      name='dzmonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nn2)
      write(31,rec= nrec) dz
      close(31)

	nrec= 12*(iYear-1948) + 2*(iMonth-1) +1

      name='kzmonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nn)
      write(31,rec= nrec)   az
      write(31,rec= nrec+1) azt
      close(31)

      write(*,88) iYear, iMonth
88    format('<<<< Data are stored!!!', i4,i3,' >>>>')
      return
      end

      subroutine Neptune

*     Version 15.09.14 - Kazantsev et. al, 1998, R=5km, myu = R2.

	include 'Slo2.fi'
      common /diffus/ al,alt

	R_Neptune= 5.0e5  ! Neptune width scale (cm) - U~R2
	alpha= 1.0        ! alpha= hx* (gradH)/H -> hx / R_Neptune
	RN2 = R_Neptune**2
	tau= 1.e8/RN2       ! Al effective = 1.e8, mui = R**2, see Kazantsev. 
   
	do j=1,jl
      S0=Si(j)

	do i=1,il
      cor= -2.*om*Co(i)*S0

	if(km2(i,j).GT.0) then
	uneptF= -alpha*dt*Al*cor*unept(i,j)
	vneptF= -alpha*dt*Al*cor*vnept(i,j)

 	det= 1.0/(1.0 +dt*tau)
       
	do k=1,km2(i,j) 
     		u(i,j,k)= det*(u(i,j,k)+uneptF)
		v(i,j,k)= det*(v(i,j,k)+vneptF)
	end do
	
	endif   ! km2>0

	end do
	end do

	return
	end

      subroutine scalar_lamped(sclr1)

*     Lamped matrix for temperature and salinity transport scheme
*     Version 31.12.2013.


      INCLUDE 'Slo2.fi'
	dimension Sclr1(0:il1,0:jl1,kl)


      c3=1./3.
	c6=1./6.
	c4=0.25
	c36=1./36.
	c24=1./24.


	sclamp= 0.0
	Sm1   = 0.0

	do j=1,jl
      do i=1,il

	IF( nt3(i,j,1) .NE. 0) then
      do k=1,km2(i,j)
	Sclr1(i,j,k)= 1.0
	end do
	end if

	end do
	end do


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
	kb=km2(i,j)

	IF( kb .GT. 0) then


      do k=1,kb
	n=abs(nt3(i,j,k))

	if(k.LT.kb)then
	kp=k+1
	np=abs(nt3(i,j,kp))
	hzk=hz(k)

	sumkp=
     *c24*(KT(1,np)*S1*(2.*Sclr1(i,j,k)+Sclr1(i+1,j,k)+Sclr1(i+1,j-1,k))
     +    +KT(2,np)*S2*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i+1,j-1,k))
     +    +KT(3,np)*S3*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i-1,j  ,k))
     +    +KT(4,np)*S4*(2.*Sclr1(i,j,k)+Sclr1(i-1,j,k)+Sclr1(i-1,j+1,k))
     +    +KT(5,np)*S5*(2.*Sclr1(i,j,k)+Sclr1(i-1,j+1,k)+Sclr1(i,j+1,k))
     +    +KT(6,np)*S6*
     *                (2.*Sclr1(i,j,k)+Sclr1(i,j+1,k)+Sclr1(i+1,j  ,k)))


	sumkpp=
     *c24*(KT(1,np)*S1*
     *    (2.*Sclr1(i,j,kp)+Sclr1(i+1,j,kp)+Sclr1(i+1,j-1,kp))
     +    +KT(2,np)*S2*
     *    (2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i+1,j-1,kp))
     +    +KT(3,np)*S3*
     *    (2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i-1,j  ,kp))
     +    +KT(4,np)*S4*
     *    (2.*Sclr1(i,j,kp)+Sclr1(i-1,j,kp)+Sclr1(i-1,j+1,kp))
     +    +KT(5,np)*S5*
     *    (2.*Sclr1(i,j,kp)+Sclr1(i-1,j+1,kp)+Sclr1(i,j+1,kp))
     +    +KT(6,np)*S6*
     *    (2.*Sclr1(i,j,kp)+Sclr1(i,j+1,kp)+Sclr1(i+1,j  ,kp)))

	else
	hzk= 0.
	sumkp  =0.
	sumkpp =0.
	end if

	if(k.GT.1)then
	km=k-1
	hzk1=hz(km)

	sumk=
     *c24*(KT(1,n)*S1*(2.*Sclr1(i,j,k)+Sclr1(i+1,j,k)+Sclr1(i+1,j-1,k))
     +    +KT(2,n)*S2*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i+1,j-1,k))
     +    +KT(3,n)*S3*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i-1,j  ,k))
     +    +KT(4,n)*S4*(2.*Sclr1(i,j,k)+Sclr1(i-1,j,k)+Sclr1(i-1,j+1,k))
     +    +KT(5,n)*S5*(2.*Sclr1(i,j,k)+Sclr1(i-1,j+1,k)+Sclr1(i,j+1,k))
     +    +KT(6,n)*S6*(2.*Sclr1(i,j,k)+Sclr1(i,j+1,k)+Sclr1(i+1,j  ,k)))


	sumkm=
     *c24*(KT(1,n)*S1*
     *    (2.*Sclr1(i,j,km)+Sclr1(i+1,j,km)+Sclr1(i+1,j-1,km))
     +    +KT(2,n)*S2*
     *    (2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i+1,j-1,km))
     +    +KT(3,n)*S3*
     *    (2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i-1,j  ,km))
     +    +KT(4,n)*S4*
     *    (2.*Sclr1(i,j,km)+Sclr1(i-1,j,km)+Sclr1(i-1,j+1,km))
     +    +KT(5,n)*S5*
     *    (2.*Sclr1(i,j,km)+Sclr1(i-1,j+1,km)+Sclr1(i,j+1,km))
     +    +KT(6,n)*S6*
     *    (2.*Sclr1(i,j,km)+Sclr1(i,j+1,km)+Sclr1(i+1,j  ,km)))
      
	else
	hzk1=0.
	sumkm=0.
	sumk =0.
	end if

	sclamp(i,j,k)= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp))

c	if(i.eq.9.and.j.eq.41.and.k.eq.kb) then 
c      write(*,*) i,j,k, sclamp(i,j,k)
c	end if


	end do   ! k

      END IF   ! nt3 >0
	end do
	end do

	return
	end

      subroutine CH4statistics

*     14.10.2013
      include 'Slo2.fi'

        CH4max= -1000.
        do j=1,jl
        do i=1,il
           if( km2(i,j) .GT. 0) then
        do k=1,km2(i,j)
           if( CH4(i,j,k) .gt. CH4max) then
            CH4max= CH4(i,j,k)
            itmax= i
            jtmax= j
            ktmax= k
            endif
        end do
           endif
        end do
        end do

      write(*,*)'CH4= ',CH4max, 'i,j,k= ',itmax,jtmax,ktmax

      return
      end

