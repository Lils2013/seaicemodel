      program framice

c     Version 05.05.2002.

      parameter(i1=15,i2=22,j0=33, finord=23.)
      parameter(Ldisc=15, mgrad=8)
      parameter(R=0.64, pi= 3.1415926, hx=pi/180.)
      Parameter (il=35, jl=49)
      Parameter (il1=il+1, jl1=jl+1)

	Parameter (Lfort=1)

*     Ice/Snow parameters
      dimension uice(0:il1,0:jl1), vice(0:il1,0:jl1),
     *          HIce (0:mgrad,0:il1,0:jl1),Aice(0:mgrad,0:il1,0:jl1),
     *          HSnow(0:mgrad,0:il1,0:jl1)

      character*72 idat, hdat
	idat= 'I000000'
      hdat= 'H000000'


	s0=SIN(pi*.5 - finord*pi/180.+(j0-1.)*hx)

      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)

      open(1, file='framimod.dat')
      open(2, file='framsnow.dat')

      write(*,*) 'Enter Starting Year-->'
      read(*,*) m1y
      write(*,*) 'Enter Ending Year-->'
      read(*,*) m2y

      do m=m1y,m2y
	do month=1,12

	time=REAL(m) +REAL(month-0.5)/12.

	Write(idat(2:5),'(i4.4)') m
      Write(hdat(2:5),'(i4.4)') m

	Write(idat(6:7),'(i2.2)') Month
      Write(hdat(6:7),'(i2.2)') Month


      open (unit=31,file=idat,status='old',access='direct',
     *      form='unformatted',recl=nuice)
      read(31,rec=1) uice
      read(31,rec=2) vice
      close(31)

      Open ( Unit=31,File=hdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nhice)
      read(31,rec=1) Hice
      read(31,rec=2) Aice
      read(31,rec=3) Hsnow
	close(31)

       
      Hfram=0.
      Afram=0.
	Flux= 0.
	SFram=0.
      do mg=1,mgrad
      do i=i1,i2
	if(i.eq.i1 .or. i.eq.i2) then
	c=0.5
	else
	c=1.
	end if
      Hfram= Hfram + C*Hice(mg,i,j0)
      Sfram= Sfram + C*Hsnow(mg,i,j0)
      Afram= Afram + C*Aice(mg,i,j0)
	Flux = Flux +  C*R*hx*s0*Hice(mg,i,j0)*vice(i,j0)
      end do
      end do

      if(Afram .GT. 0.001) then
      Hfram= 0.01*Hfram/Afram
	SFram= 0.01*SFram/Afram
      else
      Hfram=0.
	SFram=0.
      end if
      write(1,*) time, Hfram, Flux*3.*24.*36.*1.e-3
      write(2,*) time, SFram

      end do
	end do
      end
