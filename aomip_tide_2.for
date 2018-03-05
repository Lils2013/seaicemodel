      program FEMAO1

      parameter (pi=3.1415926536)
      parameter (hxgr=1., hygr=1., dthr=24.0, finord=23., alwest=-16.)
      parameter (ncond=1)
      parameter (NRestor= 1)  ! T-S Restoring [0/1]
      parameter (lprint=24)

*     Watcom and NDP Fortrans Lfort=4, Compaq and UNIX - Lfort=1
      Parameter (Lfort=1)

      INCLUDE 'slo2.fi'
	dimension MEC(12)
      dimension rs(il,jl),adz(-1:ilp,-1:jlp),serv(il,jl),serv1(il,jl)

	data MEC/31,28,31,30,31,30,31,31,30,31,30,31/

      nnn=(il+2)*(jl+2)*kl*Lfort
      nn=(il+2)*(jl+2)*Lfort
      nni=nn*mgrad
      nna=nn*(mgrad+1)
      nn2= 2*Lfort*(il+4)*(jl+4)
      nd=38*il*jl*Lfort

                         dt=3600.*dthr
                         hx=hxgr*pi/180.
                         hy=hygr*pi/180.
                         om=.729e-4
                         r=.637e9
                         g=980.

      do 52 j=0,jl1
      teta=pi*.5 - finord*pi/180.+(j-1.)*hy
52    si(j)=sin(teta)
      do 521 i=0,il1
      teta= alwest*pi/180. +(i-1.)*hx
521   co(i)=cos(teta)



c     Global number of the day.
      day= 0.
c     Number of the year.


c      open(1,file='pa_methane1.dat')
c      open(2,file='ta_methane1.dat')
c      open(3,file='wind_methane1.dat')
c      open(11,file='pa_methane2.dat')
c      open(12,file='ta_methane2.dat')
c      open(13,file='wind_methane2.dat')
c      open(21,file='pa_methane3.dat')
c      open(22,file='ta_methane3.dat')
c      open(23,file='wind_methane3.dat')



c      open(10,file='ts_methane1.dat')
c      open(11,file='ts_methane2.dat')
c      open(12,file='ts_methane3.dat')

      open(10,file='cloud_methane1.dat')
      open(11,file='cloud_methane2.dat')
      open(12,file='cloud_methane3.dat')


c      open(2,file='ta_methane1.dat')
c      open(3,file='wind_methane1.dat')
      


c     ==============================================
c                    TIME LOOP
c     ==============================================

      DO iYEAR= 1948, 2011

	idglobal = 0

	DO imonth= 1,12

      write(*,*) iyear, imonth

c     Precipitation according to the time. 
*     Initial data: Yang, D., 1999:  An improved precipitation
*     climatology for the Arctic Ocean, Geophys.
*     Res. Lett., 26(11), 1625-1628. OR:
*     Serreze, M.C., M.P. Clark and D.H. Bromwich, 2003: Monitoring
*     precipitation over the Arctic terrestrial drainage system:
*     Data requirements, shortcomings and applications of 
*     atmospheric reanalysis. Journal of Hydrometeorology, 4, 387-407.
c     Precipitation in mm/month 

      open (41,file='pr_150703.dat',access='direct',
     *      form='unformatted',recl=il*jl)
	read(41,rec=imonth) Pr
	close(41)
***      write(1,*) imonth, Pr(6,22) *0.666

c     River as rain. Comment it otherwise.
      call 
     &river_as_rain(imonth,River,il,jl,nt3,si,cg,il1,jl1,kl,r,hx,hy)

c     Cloud cover
      open (41,file='cloud.dat',access='direct',
     *      form='unformatted',recl=il*jl)
	read(41,rec=imonth) cloud
	close(41)

      write(10,*) imonth, cloud(22,6)
      write(10,*) imonth, cloud(20,14)
      write(10,*) imonth, cloud(18,20)

C     PHC 3.0 Monthly Mean Hydrology
c      call TSInt(imonth)


c	DO iday=1,MEC(imonth)
c	idglobal= idglobal +1

*     Steps in the Day
c	L=1


c      nst=nst +1
c      call Forcing3(iyear,idglobal,L,dthr,TA,Q2m,WX,WY,Pa,il,jl,
c     #              il1,jl1,si,co,hx,hy,r,om,serv,serv1) 

c      write(1,*) iyear, idglobal, Pa(22,6)
c      write(2,*) iyear, idglobal, Ta(22,6) 
c      write(3,*) iyear, idglobal, Wx(22,6), Wy(22,6) 
c      write(11,*) iyear, idglobal, Pa(20,14)
c      write(12,*) iyear, idglobal, Ta(20,14) 
c      write(13,*) iyear, idglobal, Wx(20,14), Wy(20,14) 
c      write(21,*) iyear, idglobal, Pa(18,20)
c      write(22,*) iyear, idglobal, Ta(18,20) 
c      write(23,*) iyear, idglobal, Wx(18,20), Wy(18,20) 

cccc      call shortwave(SW,Pa,Q2m,idglobal,L,dthr,cloud,il,jl)


c	END DO  ! days


	END DO  ! months
	END DO  ! years
      close(21)
      close(22)
      close(23)
      close(13)
      close(12)
      close(11)
      close(3)
      close(2)
      close(1)


      stop
      end

