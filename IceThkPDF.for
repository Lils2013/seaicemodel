      Program Ice_Distrib

*     Ice Thickness PDF.

	Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=40, klp=kl+1, mgrad=14)
      Parameter (il1=il+1, jl1=jl+1)
      Character*72 i2dat,name,file,ntdat

      Dimension nt3(0:il1,0:jl1,kl),cg(13), 
     *          Aice (0:mgrad,0:il1,0:jl1),hmax(0:mgrad)

	Dimension PDF(mgrad) ! PD function for mgrad+1 thickness gradations
	Dimension area(mgrad)

      data cg/6.,4.,5.,4.,5.,4*3.,2.,1.,2.,1./

*     Ice thickness categories. USSR classification (up to 2m)
      hmax(0)=0.
      hmax(1)=10.
      hmax(2)=20.
      hmax(3)=30.
      hmax(4)=50.
      hmax(5)=70.
      hmax(6)=100.
      hmax(7)=150.
      hmax(8)=200.
      hmax(9)=300.
      hmax(10)=400.
      hmax(11)=500.
      hmax(12)=600.
	hmax(13)=1000.
      hmax(14)=90000.


      name='H000000'
      i2dat='km2.dat'
	ntdat='nt.dat'
	file='IceTH_PDE01.dat'
      n2 = (il+2)*(jl+2)*Lfort
      n3 = n2*kl
      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)

      open (unit=11,file=ntdat,status='old',access='direct',
     *      form='unformatted',recl=n3)
      read (11,rec=1) nt3
	close(11)


	do imonth=1,12

      write(file(10:11),'(i2.2)') imonth
      open(1,file=file)


      PDF(:)= 0.
      area(:)=0.

      do m=1,mgrad

      do iYear=1958,2011
      Time= iYear +iMonth/12.
      write(name(2:5),'(i4.4)') iYear
      write(name(6:7),'(i2.2)') iMonth
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nhice)
      read(31,rec=2) Aice
      close(31)

      do i=12,20    !!! Fram Strait
	  do j=34,38  !!! Fram Strait
c      do i=1, 20  !!!il
c	  do j=1,34 !!!jl  
	  if(ABS(nt3(i,j,1)) .GT. 0) then
        PDF(m) = PDF(m) +cg(abs(nt3(i,j,1)))*Aice(m,i,j)/6.
	  area(m)=area(m) +cg(abs(nt3(i,j,1)))            /6.
	  end if
	  end do
	  end do

	  end do ! year

	  write(1,*) 0.005*(hmax(m-1)+hmax(m)), 
     *  PDF(m)/MAX(1.e-9,area(m))

	end do   ! mgrad

	close(1)

	  end do ! month

	  end