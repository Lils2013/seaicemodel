      Program MID

*     Mean Ice Drift velocity
*     Lfort=1 for Compaq and UNIX Fortrans, and
*     Lfort=4 for NDP, Watcom and MAC Fortrans.
      Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=16)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (FiN=25., aLW= -25., hxgr= 1., hygr= 1.)
      Parameter (pi=3.14159265, overfl= 1.701410e+38)
      Parameter (mgrad= 14)

      Integer*2 dummy
      Character*72 i2dat,idat,hdat

      Dimension km2(0:il1,0:jl1)
*     Ice/Snow parameters
      dimension uice(0:il1,0:jl1), vice(0:il1,0:jl1),
     *          Aice(0:mgrad,0:il1,0:jl1)

      FiS= FiN - hygr*(jlris-1.)
      ymin= FiS
      ymax= ymin +hygr*(jlris-1.)
      xmin= aLW
      xmax= xmin +hxgr*(ilris-1.)

      n2 = (il+2)*(jl+2)*Lfort
      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)
      idat= 'I000000'
      hdat= 'H000000'
      i2dat='km2.dat'

	open( unit=1, file='uimean.dat')

      Open ( Unit=49,File=i2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)

      Read(49,Rec=1) km2

      do iYear=1948,2002
	Write(idat(2:5),'(i4.4)') iYear
      Write(hdat(2:5),'(i4.4)') iYear

      
      do iMonth=1,12

	time = iyear +imonth/12.

	Write(idat(6:7),'(i2.2)') iMonth
      Write(hdat(6:7),'(i2.2)') iMonth

      Open ( Unit=31,File=hdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nhice)
      read(31,rec=2) Aice
      close(31)

      open (unit=31,file=idat,status='old',access='direct',
     *      form='unformatted',recl=nuice)
      read(31,rec=1) uice
      read(31,rec=2) vice
      close(31)

      uimean = 0.
      area   = 0.

      do i=1,il
	do j=1,jl
	IF( km2(i,j) .GT.0 .and. Aice(0,i,j).LT.0.1) THEN
	uimean= uimean + SQRT(uice(i,j)**2+vice(i,j)**2)
	area= area +1.
	END IF
	end do
	end do

	uimean= uimean/MAX(1.,area)

	write(1,*) time, uimean

	end do
	end do

      Stop
      End
