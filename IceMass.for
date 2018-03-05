      Program Ice_Mass

*     Ice and Snow Volumes.

	Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=40, klp=kl+1, mgrad=14)
      Parameter (il1=il+1, jl1=jl+1)
      Character*72 i2dat,name,ntdat
	Real*4 IMass, SMass

      Dimension nt3(0:il1,0:jl1,kl),cg(13), 
     *          HIce (0:mgrad,0:il1,0:jl1),
     *          HSnow(0:mgrad,0:il1,0:jl1)
      data cg/6.,4.,5.,4.,5.,4*3.,2.,1.,2.,1./

      name='H000000'
      i2dat='km2.dat'
	ntdat='nt.dat'

      open(1,file='IceMass.dat')

      n2 = (il+2)*(jl+2)*Lfort
      n3 = n2*kl
      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)

      open (unit=11,file=ntdat,status='old',access='direct',
     *      form='unformatted',recl=n3)
      read (11,rec=1) nt3

      do iYear=1948,2011
	do imonth=1,12
      Time= iYear +iMonth/12.
      write(name(2:5),'(i4.4)') iYear
      write(name(6:7),'(i2.2)') iMonth
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nhice)
      read(31,rec=1) Hice
      read(31,rec=3) Hsnow
      close(31)

      IMass= 0.
	SMass=0.
      do i=1,il
	  do j=1,jl
	  if(ABS(nt3(i,j,1)) .GT. 0) then
	  do m=1,mgrad
      IMass =IMass +cg(abs(nt3(i,j,1)))*Hice (m,i,j)/6.
      SMass=SMass+cg(abs(nt3(i,j,1)))*Hsnow(m,i,j)/6.
	  end do
	  end if
	  end do
	  end do


	  write(1,*) time, 0.115*IMass, 0.115*SMass

	  end do ! month
	  end do ! year
	  end