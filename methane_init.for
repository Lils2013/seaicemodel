      program initslo11


*     RFFI OFI-MProject Program with Sea Ice

*    ------------   Last modification date 07.06.15  --------------

*     Arctic Ocean model initialization, 1 x 1 degree.
*     PHC 3.0 initial temperature and salinity.
*     h.grd is wider, then model arrays, and array servgrd(51,51)
*     is essential. For a filtration arrays a,b are organized.

*     Methane is initialised according to Lein and Ivanov, 2009
*     Oxygen is 8 mL/L everywhere

*     Large arrays for cluster

*     Watcom and NDP Fortrans Lfort=4, Compaq and UNIX - Lfort=1
      Parameter(Lfort=1)

      parameter (il=35,jl=49,kl=40,ilp=il+2,jlp=jl+2,klp=kl+1)
      parameter (il1=il+1,jl1=jl+1)
      parameter (finord=23.,hxgr=1.,hygr=1.,re=.637e9, pi=3.1415926)
      parameter (alwest=-16.)

*     Number of sea ice thickness gradations.
      parameter (mgrad= 14)

*     General Ocean Arrays.

      dimension dz(-1:ilp,-1:jlp),t(0:il1,0:jl1,kl),dzext(il,jl),
     *                            s(0:il1,0:jl1,kl),
     *  si(0:jl1), co(0:jl1), z(klp), hz(kl), area(kl),
     *             u(0:il1,0:jl1,kl), v(0:il1,0:jl1,kl),
     *             az(il,jl,kl), ub(0:il1,0:jl1,kl), 
     *             unept(0:il1,0:jl1), vnept(0:il1,0:jl1)

*     Sea Ice Arrays.
      dimension uice(0:il1,0:jl1), Hice(0:mgrad,0:il1,0:jl1),
     *          Tice(0:mgrad,0:il1,0:jl1)

      dimension Flux_CH4_To_A_month(0:il1,0:jl1)

*     Model Geometry and Service Arrays.

      dimension servgrd(51,51), a(0:52,0:52), b(0:52,0:52)
	dimension serv(il,jl,kl)
      integer   nt3(0:il1,0:jl1,kl),km2(0:il1,0:jl1)
      double precision dz, time
cc    *,xmin,xmax,ymin,ymax,zmin,zmax

      character*72 udat,tsdat,dzdat,nt3dat,ubdat,km2dat,CH4dat,
     *             CH4obsdat,O2dat,O2obsdat,
     *             tsobserv, aidat,hisdat,uidat,tmpisdat,file,name
      data udat/'u.dat'/ ,tsdat/'t.dat'/,dzdat/'dz.dat'/,
     *     nt3dat/'nt.dat'/,ubdat/'ub.dat'/, km2dat/'km2.dat'/,
     *     tsobserv/'tsobs.dat'/,CH4dat/'CH4.dat'/, 
     *     CH4obsdat/'CH4obs.dat'/,
     *     O2dat/'O2.dat'/,O2obsdat/'O2obs.dat'/,
     *     aidat/'ai.dat'/, hisdat/'his.dat'/,uidat/'ui.dat'/,
     *     tmpisdat/'tmpis.dat'/


      om= .729e-4
      g= 980.
      r= re
      teta0=0.5*pi-finord*pi/180.
      al0  =alwest*pi/180.
      hx= hxgr*pi/180.
      hy= hygr*pi/180.
      do 30 j=1,jl
      teta=teta0+(j-1)*hy
	co(j)=cos(teta)
30    si(j)=sin(teta)

      do k=1,kl
      area(k)=0.
      end do

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

      do 33 k=1,klp
33    z(k)=100.*z(k)
      do 34 k=1,kl
34    hz(k)=(z(k+1)-z(k))

      nnn=(il+2)*(jl+2)*kl*Lfort
      nn=(il+2)*(jl+2)*Lfort
      nni=nn*mgrad
      nna=nn*(mgrad+1)
      nn0=il*jl*Lfort
      nnt=il*jl*kl*Lfort
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
      open (unit=44,file=CH4dat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=45,file=CH4obsdat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=46,file=O2dat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=47,file=O2obsdat,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=15,file=tsobserv,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=16,file=dzdat,status='old',access='direct',
     *      form='unformatted',recl=nn2)
      open (unit=17,file=km2dat,status='old',access='direct',
     *      form='unformatted',recl=nn)
      open (unit=33,file='neptune.dat',status='old',access='direct',
     *      form='unformatted',recl=nn)

c     Ice Data

      open (unit=26,file=hisdat,status='old',access='direct',
     *      form='unformatted',recl=nna)
      open (unit=27,file=aidat,status='old',access='direct',
     *      form='unformatted',recl=nna)
      open (unit=28,file=uidat,status='old',access='direct',
     *      form='unformatted',recl=nn)
      open (unit=29,file=tmpisdat,status='old',access='direct',
     *      form='unformatted',recl=nna)

      open(22,file='time.dat')

                time= 0.d0
                write(22,*) time
      write(*,*) 'Time is stored'

      call rd_grd('E:\methan\bottom\H5.grd',51,51,servgrd)

      call tfiltr(1,servgrd,a,b,51,51,52,52)
      call km2form(servgrd, km2, il1,jl1)
      call ntform (km2, nt3, il1, jl1, kl)


cccc	goto 7777
*     ---------------------------------------------
*             Marking of Liquid Boundary
*     ---------------------------------------------
      do k=1,kl
*     Canadian Arhipelago
      nt3(10,28,k)= -nt3(10,28,k)
	nt3(3,16,k) = -nt3(3,16,k)
	nt3(3,17,k) = -nt3(3,17,k)
	
****      nt3(1,12,k)= -nt3(1,12,k)

	do j=18,21
	nt3(6,j,k) = -nt3(6,j,k)
	end do

*     Bering Passage
      do i=11,15
      nt3(i,1,k)= -nt3(i,1,k)
      end do

*     Denmark Strait
      do j=43,46
      nt3(7,j,k)= -nt3(7,j,k)
      end do

*     Norwegian sea
      do i=11,23
      nt3(i,jl,k)= -nt3(i,jl,k)
      end do

      end do  ! K Loop

7777  continue

      Write(11,rec=1) nt3
      Write(17,rec=1) km2
      write(*,*)' Array K is stored'
      call print2(il1,jl1,km2,nt3(0,0,1))

*     Observed PHC 3.0 T, S.
*     December-January. 
*     Temperature in situ --> T potential.

	open (51,file='SPHC.dat',access='direct',form='unformatted',
     *      recl=il*jl*kl)
	read(51,rec=12) Serv
	close(51)
      do i=1,il
      do j=1,jl
      do k=1,kl
      s(i,j,k)= max(0.,serv(i,jl-j+1,k))      
      end do
      end do
      end do


	open (51,file='TPHC.dat',access='direct',form='unformatted',
     *      recl=il*jl*kl)
	read(51,rec=12) Serv
	close(51)

c     Polynomial coefficients.

      a1 =  8.654839e-6
      a2 = -1.416363e-6
      a3 = -7.382865e-9
      a4 = -8.382414e-6
      a5 =  2.839334e-8
      a6 =  1.778040e-8
      a7 =  1.711556e-10

      do i=1,il
      do j=1,jl
      do k=1,kl
	Tin = serv(i,jl-j+1,k)
	Tin2= Tin**2
	Tin3= Tin**3
	ppp = 1.e-5*g*1.020*z(k)
	sss = s(i,j,k)

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
      t(i,j,k)= Tin +Poly*ppp
     
      end do
      end do
      end do

*     Blanking and Freezing point check.

      do i=1,il
      do j=1,jl
      do k=1,kl
      if( km2(i,j) .GE. k) then

	if(S(i,j,k) .LT. 0. .OR. S(i,j,k) .GT. 100. ) 
     &write(*,*) 'Error in S', S(i,j,k),i,j,k

	if(T(i,j,k) .LT. -5. .OR. T(i,j,k) .GT. 100. ) 
     &write(*,*) 'Error in T', T(i,j,k),i,j,k


	sss= MAX(0.,S(i,j,k))
	ppp = 1.e-5*g*1.020*z(k)
      TFriz= Tfr(sss,ppp)
      if( t(i,j,k) .LT. tfriz) t(i,j,k)= tfriz
      
      if(sss .GT. 35.7) write(*,*) 'Large S', sss, i,j,k
      
      else
        t(i,j,k)= 0.
        s(i,j,k)= 0.
      END IF

      end do
      end do
      end do

*     Iitial conditions
      write(14,rec=1) t
      write(14,rec=2) s

      name='tsmonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= 1) t
      write(31,rec= 2) s
      close(31)


*     Methane

      T= 3.5
	S= 3.5

*     Boundary conditions in straits and passages
*     Russian book Lein and Ivanov, 2009

      S(:,jl,:)    = 2.4 ! North Atlantics
	S(:,1,: )    = 8.2 ! Bering Strait
      S(7,42:jl,:) = 2.4 ! Denmark Strait=North Atlantics

      write(44,rec=1) T
      write(45,rec=1) S
	close(44)
	close(45)

      name='ch4month.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= 1) t
      close(31)


*     Oxygen 8 mL/L = 8*1e6 /22.4L nmol/L

      T= 6.e6/22.4
	S= 6.e6/22.4

      write(46,rec=1) T
      write(47,rec=1) S
	close(46)
	close(47)

      name='o2month.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= 1) t
      close(31)


      Flux_CH4_To_A_month = 0.0
	 
      name='fluxch4month.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nn)
      write(31,rec= 1) Flux_CH4_To_A_month
      close(31)


*     Initial velocities.
      do i=0,il1
      do j=0,jl1
      do k=1,kl
      u(i,j,k)= 0.
      v(i,j,k)= 0.
      ub(i,j,k)= 0.
      end do
      end do
      end do

      write(12,rec=1) u
      write(12,rec=2) u
      write(12,rec=3) u

      name='uomonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      write(31,rec= 1) u
      write(31,rec= 2) u
      write(31,rec= 3) u
      close(31)

*     Initial turbulence corfficients

      name='kzmonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nnt)
      write(31,rec= 1) az
      write(31,rec= 2) az
      close(31)


*     Sea Ice Parameters. Initial velocities.

      do i=0,il1
      do j=0,jl1
      uice(i,j)= 0.
      end do
      end do
      write (28,rec=1) uice
      write (28,rec=2) uice


      name='uimonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nn)
      write(31,rec=1) uice
      write(31,rec=2) uice
      close(31)


*     Initial Ice Thickness.

      read(14,rec=1) t
      read(14,rec=2) s

      Hice=0.
      do m=0,mgrad
      do i=0,il1
      do j=0,jl1
      Hice(m,i,j) = 0.
	if( m.eq.9 .AND. km2(i,j) .GT. 0)then
	sss= max(0.,S(i,j,1))
      TFriz= Tfr(sss,0.)
      if( t(i,j,1) .LE. tfriz+1.0) Hice(m,i,j)= 300.
      end if
      end do
      end do
      end do

      write (26,rec=1) Hice

*     Initial sea level 
      
      dz(:,:) =0.d0
      do 3 i=0,il1
      do 3 j=0,jl1
      dz(i,j)= 0.d0
      if(hice(9,i,j) .GT. 0.) dz(i,j)= -300.d0
3     continue
      Write(16,rec=1) dz

      name='dzmonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nn2)
      write(31,rec= 1) dz
      close(31)
      

*     Initial Snow thickness.

      Hice= 0.
      do m=0,mgrad
      do i=0,il1
      do j=0,jl1
	if( m.eq.9 .AND. km2(i,j) .GT. 0)then
	sss= max(0.,S(i,j,1))
      TFriz= Tfr(sss,0.)
      if( t(i,j,1) .LE. tfriz+1.0) Hice(m,i,j)= 25.
      end if
      end do
      end do
      end do
      write (26,rec=2) Hice

*     Initial Compactness.

      Hice=0.
      do m=1,mgrad
      do i=0,il1
      do j=0,jl1
      Hice(m,i,j)= 0.
	if( m.eq.9 .AND. km2(i,j) .GT. 0)then
	sss= max(0.,S(i,j,1))
      TFriz= Tfr(sss,0.)
      if( t(i,j,1) .LE. tfriz+1.0) then
	Hice(m,i,j)= 0.99
	end if
	end if
      end do
      end do
      end do

	Hice(0,:,:) = 1.

	do j=1,jl
	do i=1,il
	do m=1,mgrad
	Hice(0,i,j)= Hice(0,i,j)-Hice(m,i,j)
	end do
	end do
	end do

      write (27,rec=1) Hice

*     Initial Snow/Ice Temperature.

      do m=0,mgrad
      do i=0,il1
      do j=0,jl1
      Tice(m,i,j) = -30.0
      end do
      end do
      end do
      write (29,rec=1) Tice
      write (29,rec=2) Tice


      hice= 0.

      name='himonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nna)
      write(31,rec= 1) hice
      write(31,rec= 2) hice
      write(31,rec= 3) hice
      close(31)


*     Model area geometry.
      x=hx*hy*r*r
      do 17 j=1,jl
      s0=si(j)
      do 17 i=1,il
      if( km2(i,j) .GT. 0) then
         do 18 k=1,km2(i,j)
         n= ABS(nt3(i,j,k))
         if( n.eq.0) cg=0.
         if( n.eq.1) cg=6.
         if( n.eq.3 .or. n.eq.5) cg=5.
         if( n.eq.2 .or. n.eq.4) cg=4.
         if( n.eq.10 .or. n.eq.12) cg=2.
         if( n.eq.11 .or. n.eq.13) cg=1.
         if( n.eq.6 .or. n.eq.7 .or. n.eq.8 .or. n.eq.9) cg=3.
         area(k)=area(k)+cg*x*s0/6.
 18      continue
      end if
 17   continue
      volume=0.
      do 23 k=2,kl
 23   volume=volume+area(k)*hz(k-1)

*
      do k=1,kl
      Write(*,*)' hz=', hz(k)
      end do
      Write(*,*)' volume=', volume
      do k=1,kl
      Write(*,*)' area=', area(k)
      end do

*     Number of points in the 2D model area: it's necessary for
*     GMRES, BiCGSTAB and DUMKA algorithms.
      len=0
      do i=1,il
      do j=1,jl
      if( km2(i,j) .GT. 0) then
      len= len+1
      end if
      end do
      end do
      write(*,*)' ************************************************'
      write(*,*)' LEN =', len
      write(*,*)' ************************************************'

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(26)
      close(27)
      close(28)
      close(29)

      stop
      end

      subroutine km2form( s, km2, il1,jl1)

*     Corrected version 11.03.2012

*     Arctic Ocean.
*     Note index shift in s(51,51)

      dimension s(51,51), km2(0:il1,0:jl1)

*     If ibottom=1, bottom gets shallower, otherwise ibottom=2.
      ibottom= 1

      il=il1-1
      jl=jl1-1
      do i=0,il1
      do j=0,jl1
      km2(i,j)= 0
      end do
      end do

      write(*,*)' iBottom =', ibottom
      kkl= ibottom +2
      write(*,*)' For the barotropization kkl=', kkl

      do 55 i=1,il
      ii= i+9

      do 55 j=1,jl
      IF( s(ii,jl-j+1) .GT. 0. .AND. s(ii,jl-j+1) .LT. 1.e+30) THEN

      if( ibottom .eq. 1) then
      IF(                            s(ii,jl-j+1).LT.10.)   km2(i,j)= 4
      IF(s(ii,jl-j+1).GE.10.   .AND. s(ii,jl-j+1).LT.15.)   km2(i,j)= 5
      IF(s(ii,jl-j+1).GE.15.   .AND. s(ii,jl-j+1).LT.20.)   km2(i,j)= 6
      IF(s(ii,jl-j+1).GE.20.   .AND. s(ii,jl-j+1).LT.27.5)  km2(i,j)= 7
      IF(s(ii,jl-j+1).GE.27.5  .AND. s(ii,jl-j+1).LT.37.5)  km2(i,j)= 8
      IF(s(ii,jl-j+1).GE.37.5  .AND. s(ii,jl-j+1).LT.50.)   km2(i,j)= 9
      IF(s(ii,jl-j+1).GE.50.   .AND. s(ii,jl-j+1).LT.65.)   km2(i,j)= 10
      IF(s(ii,jl-j+1).GE.65.   .AND. s(ii,jl-j+1).LT.82.5)  km2(i,j)= 11
      IF(s(ii,jl-j+1).GE.82.5  .AND. s(ii,jl-j+1).LT.100.)  km2(i,j)= 12
      IF(s(ii,jl-j+1).GE.100.  .AND. s(ii,jl-j+1).LT.125.)  km2(i,j)= 13
      IF(s(ii,jl-j+1).GE.125.  .AND. s(ii,jl-j+1).LT.150.)  km2(i,j)= 14
      IF(s(ii,jl-j+1).GE.150.  .AND. s(ii,jl-j+1).LT.175.)  km2(i,j)= 15
      IF(s(ii,jl-j+1).GE.175.  .AND. s(ii,jl-j+1).LT.200.)  km2(i,j)= 16
      IF(s(ii,jl-j+1).GE.200.  .AND. s(ii,jl-j+1).LT.225.)  km2(i,j)= 17
      IF(s(ii,jl-j+1).GE.225.  .AND. s(ii,jl-j+1).LT.250.)  km2(i,j)= 18
      IF(s(ii,jl-j+1).GE.250.  .AND. s(ii,jl-j+1).LT.275.)  km2(i,j)= 19
      IF(s(ii,jl-j+1).GE.275.  .AND. s(ii,jl-j+1).LT.300.)  km2(i,j)= 20
      IF(s(ii,jl-j+1).GE.300.  .AND. s(ii,jl-j+1).LT.325.)  km2(i,j)= 21
      IF(s(ii,jl-j+1).GE.325.  .AND. s(ii,jl-j+1).LT.350.)  km2(i,j)= 22
      IF(s(ii,jl-j+1).GE.350.  .AND. s(ii,jl-j+1).LT.375.)  km2(i,j)= 23
      IF(s(ii,jl-j+1).GE.375.  .AND. s(ii,jl-j+1).LT.410.)  km2(i,j)= 24
      IF(s(ii,jl-j+1).GE.410.  .AND. s(ii,jl-j+1).LT.450.)  km2(i,j)= 25
      IF(s(ii,jl-j+1).GE.450.  .AND. s(ii,jl-j+1).LT.500.)  km2(i,j)= 26
      IF(s(ii,jl-j+1).GE.500.  .AND. s(ii,jl-j+1).LT.600.)  km2(i,j)= 27
      IF(s(ii,jl-j+1).GE.600.  .AND. s(ii,jl-j+1).LT.700.)  km2(i,j)= 28
      IF(s(ii,jl-j+1).GE.700.  .AND. s(ii,jl-j+1).LT.850.)  km2(i,j)= 29
      IF(s(ii,jl-j+1).GE.850.  .AND. s(ii,jl-j+1).LT.1050.) km2(i,j)= 30
      IF(s(ii,jl-j+1).GE.1050. .AND. s(ii,jl-j+1).LT.1250.) km2(i,j)= 31
      IF(s(ii,jl-j+1).GE.1250. .AND. s(ii,jl-j+1).LT.1500.) km2(i,j)= 32
      IF(s(ii,jl-j+1).GE.1500. .AND. s(ii,jl-j+1).LT.1750.) km2(i,j)= 33
      IF(s(ii,jl-j+1).GE.1750. .AND. s(ii,jl-j+1).LT.2100.) km2(i,j)= 34
      IF(s(ii,jl-j+1).GE.2100. .AND. s(ii,jl-j+1).LT.2500.) km2(i,j)= 35
      IF(s(ii,jl-j+1).GE.2500. .AND. s(ii,jl-j+1).LT.3000.) km2(i,j)= 36
      IF(s(ii,jl-j+1).GE.3000. .AND. s(ii,jl-j+1).LT.3500.) km2(i,j)= 37
      IF(s(ii,jl-j+1).GE.3500. .AND. s(ii,jl-j+1).LT.4000.) km2(i,j)= 38
      IF(s(ii,jl-j+1).GE.4000. .AND. s(ii,jl-j+1).LT.4500.) km2(i,j)= 39
      IF(s(ii,jl-j+1).GE.4500.)                             km2(i,j)= 40

      else
      IF(                            s(ii,jl-j+1).LT.75.)  km2(i,j)= 4
      IF(s(ii,jl-j+1).GE.75.   .AND. s(ii,jl-j+1).LT.125.) km2(i,j)= 5
      IF(s(ii,jl-j+1).GE.125.  .AND. s(ii,jl-j+1).LT.175.) km2(i,j)= 6
      IF(s(ii,jl-j+1).GE.175.  .AND. s(ii,jl-j+1).LT.225.) km2(i,j)= 7
      IF(s(ii,jl-j+1).GE.225.  .AND. s(ii,jl-j+1).LT.275.) km2(i,j)= 8
      IF(s(ii,jl-j+1).GE.275.  .AND. s(ii,jl-j+1).LT.350.) km2(i,j)= 9
      IF(s(ii,jl-j+1).GE.350.  .AND. s(ii,jl-j+1).LT.450.) km2(i,j)= 10
      IF(s(ii,jl-j+1).GE.450.  .AND. s(ii,jl-j+1).LT.675.) km2(i,j)= 11
      IF(s(ii,jl-j+1).GE.675.  .AND. s(ii,jl-j+1).LT.875.) km2(i,j)= 12
      IF(s(ii,jl-j+1).GE.875.  .AND. s(ii,jl-j+1).LT.1500.)km2(i,j)= 13
      IF(s(ii,jl-j+1).GE.1500. .AND. s(ii,jl-j+1).LT.2500.)km2(i,j)= 14
      IF(s(ii,jl-j+1).GE.2500. .AND. s(ii,jl-j+1).LT.3500.)km2(i,j)= 15
      IF(s(ii,jl-j+1).GE.3500.                            )km2(i,j)= 16
      endif

      ELSE

      km2(i,j)= 0
      END IF
55    continue

*     Standart array modification

999   info= 1

      write(*,*)' There were modifications in the km2'

      km2(1 ,1 )     =MIN(km2(1 ,2   ),km2(2   ,1 ),km2(2,2))
      km2(1 ,jl)     =MIN(km2(1 ,jl-1),km2(2   ,jl),km2(2,jl-1))
      km2(il,1 )     =MIN(km2(il,2   ),km2(il-1,1 ),km2(il-1,2))
      km2(il,jl)     =MIN(km2(il,jl-1),km2(il-1,jl),km2(il-1,jl-1))


      do i=2,il-1
      if(km2(i-1,1).LT.km2(i,1).AND.km2(i+1,1).LT.km2(i,1))then
      j=1
      info= 0
      km2(i,1)= MAX(km2(i-1,1),km2(i+1,1))
      endif
      if(km2(i-1,jl).LT.km2(i,jl).AND.km2(i+1,jl).LT.km2(i,jl))then
      j= jl
      info= 0
      km2(i,jl)= MAX(km2(i-1,jl),km2(i+1,jl))
      endif
      if(km2(i,1).GT.km2(i,2))then
      j=1
      info= 0
      km2(i,1)= km2(i,2)
      endif
      if(km2(i,jl).GT.km2(i,jl-1))then
      j= jl
      info= 0
      km2(i,jl)= km2(i,jl-1)
      endif
      end do

      do j=2,jl-1
      if(km2(1,j-1).LT.km2(1,j).AND.km2(1,j+1).LT.km2(1,j))then
      i= 1
      info= 0
      km2(1,j)= MAX(km2(1,j-1),km2(1,j+1))
      endif
      if(km2(il,j-1).LT.km2(il,j).AND.km2(il,j+1).LT.km2(il,j))then
      i= il
      info= 0
      km2(il,j)= MAX(km2(il,j-1),km2(il,j+1))
      endif
      if(km2(1,j).GT.km2(2,j))then
      i= 1
      info= 0
      km2(1,j)= km2(2,j)
      endif
      if(km2(il,j).GT.km2(il-1,j))then
      i= il
      info= 0
      km2(il,j)= km2(il-1,j)
      endif
      end do

      do i=2,il-1
      do j=2,jl-1
      if(km2(i-1,j).LT.km2(i,j).AND.km2(i+1,j).LT.km2(i,j))then
      write(*,*) i, j
      info= 0
      km2(i,j)= MAX(km2(i-1,j),km2(i+1,j))
      endif
      if(km2(i,j-1).LT.km2(i,j).AND.km2(i,j+1).LT.km2(i,j))then
      info= 0
      km2(i,j)= MAX(km2(i,j-1),km2(i,j+1))
      endif
      end do
      end do

*     Singular points search
      write(*,*)' Singular points search...'

      do i=2,il-1
      do j=2,jl-1

      if(km2(i-1,j-1).LT.km2(i,j).AND.km2(i+1,j+1).LT.km2(i,j)
     &.and. km2(i+1,j-1).GE.km2(i,j) .and. km2(i-1,j+1).GE.km2(i,j)
     &                                                        )then
      write(*,*) i, j
      km2(i,j)= max(km2(i-1,j-1), km2(i+1,j+1))
      write(*,*) ' OK !!!'
      info= 0
      end if

      if(km2(i+1,j-1).LT.km2(i,j).AND.km2(i-1,j+1).LT.km2(i,j)
     &.and. km2(i-1,j-1).GE.km2(i,j) .and. km2(i+1,j+1).GE.km2(i,j)
     &                                                        )then
      write(*,*) i, j
      km2(i,j)= max(km2(i+1,j-1), km2(i-1,j+1))
      write(*,*) ' OK !!!'
      info= 0
      end if


      end do
      end do

	do i=1,il-1
	do j=1,jl-1
	if ((km2(i,j).lt.km2(i+1,j)).AND.(km2(i+1,j+1).lt.km2(i+1,j)).AND.
     *   (km2(i,j).lt.km2(i,j+1)).AND.
     *   (km2(i+1,j+1).lt.km2(i,j+1))) then
      write(*,*) 'singular 1', i, j
		km2(i+1,j)=km2(i,j)
	info=0
	end if

	if ((km2(i+1,j).lt.km2(i,j)).AND.(km2(i,j+1).lt.km2(i,j)).AND.
     *   (km2(i+1,j).lt.km2(i+1,j+1)).AND.
     *   (km2(i,j+1).lt.km2(i+1,j+1))) then
		km2(i,j)=km2(i+1,j)
      write(*,*) 'singular 2', i, j
	info=0
	end if

	enddo
	enddo


	if( info .EQ. 0) GOTO 999

      return
      end

      subroutine ntform(km2,nt3,il1,jl1,kl)
      dimension km2(0:il1,0:jl1), nt3(0:il1,0:jl1,kl)
*     Point type array.
*     Version 26.05.97 with wide arrays.
      do i=0,il1
      do j=0,jl1
      do k=1,kl
      nt3(i,j,k)=0
      end do
      end do
      end do

      do k=1,kl
      DO 11 i=1,il1-1
      DO 11 j=1,jl1-1
      IF(km2(i,j) .GE. k) THEN
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=1
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).LT.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=2
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).LT.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=3
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=4
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).LT.k)
     *       nt3(i,j,k)=5
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=6
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=7
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND.
     *   km2(i-1,j-1).GE.k                        )
     *       nt3(i,j,k)=8
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=9
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j-1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i+1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j-1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i+1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i-1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *   km2(i-1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i-1,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j+1).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j+1).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i-1,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *                           km2(I+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *                           km2(I+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      IF(km2(i-1,j+1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      IF(km2(i-1,j+1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      END IF
 11   CONTINUE
      end do

      return
      end

      subroutine print2(n,m,km2,nt)
      dimension  km2(0:n,0:m), nt(0:n,0:m)
c     Arrays are printed vertically.
1     format(1x,20i4)
2     format(1x,9i4)

      write(*,*)'   **** Km2 ****'
      write(*,1) ((km2(i,m-j), j=1 ,20), i=1,n-1)
      write(*,*)'  '
      write(*,1) ((km2(i,m-j), j=21,40), i=1,n-1)
      write(*,*)'  '
      write(*,2) ((km2(i,m-j), j=41,49), i=1,n-1)
      write(*,*)'  '

      write(*,*)'   **** Nt ****'
      write(*,1) ((nt(i,m-j), j=1 ,20), i=1,n-1)
      write(*,*)'  '
      write(*,1) ((nt(i,m-j), j=21,40), i=1,n-1)
      write(*,*)'  '
      write(*,2) ((nt(i,m-j), j=41,49), i=1,n-1)
      write(*,*)'  '

      return
      end
*
      subroutine print3(a,h,n,m,l)
      dimension a(0:n,0:m,l), h(0:n,0:m)
      do k=1,l,3
      write(*,*)
      write(*,*)' k=', k
         do i=0,n
         do j=0,m
         h(i,j)= a(i,j,k)
         end do
         end do
      write(*,1) ((h(i,m-j), j=1 ,20), i=1,n-1)
      write(*,*)'  '
      write(*,1) ((h(i,m-j), j=21,40), i=1,n-1)
      write(*,*)'  '
      write(*,2) ((h(i,m-j), j=41,49), i=1,n-1)
      write(*,*)'  '
      end do
1     format(1x,20f4.0)
2     format(1x,9f4.0)
      return
      end

      subroutine rd_grd (file,Nx,Ny,DATA)
      CHARACTER*(*) file
      dimension DATA(Nx,Ny)
      CHARACTER*4 Z
      integer ii,jj

      OPEN(9,FILE=file,status='old')
      READ(9,'(a4)') Z
      READ(9,*) ii,jj
      READ(9,*) Xmin, Xmax
      READ(9,*) Ymin, Ymax
      READ(9,*) Zmin, Zmax
      DO j=1,Ny
      READ(9,*) (data(i,j),i=1,Nx)
      ENDDO
      close(9)
      return
      end

      subroutine tfiltr(nfilter,serv,a,b,il,jl,ilp,jlp)
      dimension a(0:ilp,0:jlp), b(0:ilp,0:jlp), serv(il,jl)
*     Tukey COS filtration.

      if( nfilter .EQ. 0) goto 100
      do i=0,ilp
      do j=0,jlp
      a(i,j)= 0.
      b(i,j)= 0.
      end do
      end do
      do i=1,il
      do j=1,jl
      a(i,j)= serv(i,j)
      end do
      end do

*     Filtartion by Tukey filter

      do m=1,nfilter

      do i=1,il
      do j=1,jl
      if( a(i,j).GT.0. .AND. a(i,j) .LT. 1.e9) then

          c1=0.
          c2=0.
          c3=0.
          c4=0.
          d1=0.
          d2=0.
          d3=0.
          d4=0.

          if(a(i,j+1).GT.0..AND. a(i,j+1) .LT. 1.e9) c1=1.
          if(a(i+1,j).GT.0..AND. a(i+1,j) .LT. 1.e9) c2=1.
          if(a(i-1,j).GT.0..AND. a(i-1,j) .LT. 1.e9) c3=1.
          if(a(i,j-1).GT.0..AND. a(i,j-1) .LT. 1.e9) c4=1.

          if(a(i+1,j+1).GT.0..AND. a(i+1,j+1) .LT. 1.e9) d1=1.
          if(a(i+1,j-1).GT.0..AND. a(i+1,j-1) .LT. 1.e9) d2=1.
          if(a(i-1,j+1).GT.0..AND. a(i-1,j+1) .LT. 1.e9) d3=1.
          if(a(i-1,j-1).GT.0..AND. a(i-1,j-1) .LT. 1.e9) d4=1.

      csum=c1+c2+c3+c4
      dsum=d1+d2+d3+d4

      if(csum.gt.0. .and. dsum.gt.0.) then
      b(i,j)= 0.28*a(i,j)+0.52*(c1*a(i,j+1)+c2*a(i+1,j)+
     +                          c3*a(i-1,j)+c4*a(i,j-1))/csum +
     +                    0.2*(d1*a(i+1,j+1)+d2*a(i+1,j-1)+
     +                         d3*a(i-1,j+1)+d4*a(i-1,j-1))/dsum
      else
      b(i,j)= -100.
      end if

      end if
      end do
      end do

      do i=0,ilp
      do j=0,jlp
      a(i,j)= b(i,j)
      end do
      end do

      end do

      do i=1,il
      do j=1,jl
      if( a(i,j) .gt. 0..AND. a(i,j) .LT. 1.e9) then
      serv(i,j)= a(i,j)
      else
      serv(i,j)= -1.e10
      end if
      end do
      end do
100   return
      end

      function TFr(S,p)

*     Version 06.06.2012.
*     For global models it's recommended to use simplified formula 
*     to estimate the Tfr

	real S,p
	real*8 s8,p8,pn,pd

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
      sqs=SQRT(S8)


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