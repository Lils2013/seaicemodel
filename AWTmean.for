      Program AWTemp
	  Parameter (Lfort=1)

	  Parameter (il=35, jl=49, kl=40, klp=kl+1)
      Parameter (il1=il+1, jl1=jl+1)
      Character*72 i2dat,tdat,ntdat

      Dimension T(0:il1,0:jl1,kl), nt3(0:il1,0:jl1,kl),
     *          km2(0:il1,0:jl1), cg(13), z(klp), hz(kl)
      data cg/6.,4.,5.,4.,5.,4*3.,2.,1.,2.,1./
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

      do 33 k=1,kl+1
33    z(k)=100.*z(k)
      do 34 k=1,kl
34    hz(k)=(z(k+1)-z(k))

      tdat= 't000000'
      i2dat='km2.dat'
	  ntdat='nt.dat'

      open(1,file='AWTmean.dat')

      n2 = (il+2)*(jl+2)*Lfort
      n3 = n2*kl
      Open ( Unit=49,File=i2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)
      Read(49,Rec=1) km2
      open (unit=11,file=ntdat,status='old',access='direct',
     *      form='unformatted',recl=n3)
      read (11,rec=1) nt3


      do iYear=1948,2011
      Write(tdat(2:5),'(i4.4)') iYear

      do iMonth=1,12
      Write(tdat(6:7),'(i2.2)') iMonth
      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=1) T

      Time= iYear +iMonth/12.
      Tmean= 0.
      volume=0.
      do i=1,il
	  do j=1,jl-15
	  if(km2(i,j) .GE. 31) then
	  do k=21,30
      Tmean=Tmean+cg(abs(nt3(i,j,k)))*(T(i,j,k)+T(i,j,k+1))*hz(k)/12.
      volume= volume +cg(abs(nt3(i,j,k)))*hz(k)/6.
	  end do
	  end if
	  end do
	  end do

	  Tmean = Tmean/(max(1.,volume))

	  write(1,*) time, Tmean

	  end do
	  end do
	  end