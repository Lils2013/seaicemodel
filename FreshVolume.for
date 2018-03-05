       Program Fresh_Volume

*     version 27.10.2009 

*     Time series of FWC [km3]
*     Lfort=1 for Compaq and UNIX Fortrans, and
*     Lfort=4 for NDP, Watcom and MAC Fortrans.
      Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=40)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (ilp= il+2, jlp= jl+2)
      Parameter (hxgr= 1., hygr= 1., finord=23.)
      Parameter (pi=3.14159265, r=.637e4)
      Character*72 i2dat,nt3dat,tdat,file,fileout

      Dimension s(0:il1,0:jl1,kl),km2(0:il1,0:jl1), 
     *          si(jl1),hz(kl),z(kl)
     
      hx=hxgr*pi/180.
      hy=hygr*pi/180.
	r2hy2=r*hy*r*hx

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

      do 4 k=1,kl
4     hz(k)=(z(k+1)-z(k))/1000.

      do 52 j=0,jl1
      teta=pi*.5 - finord*pi/180.+(j-1.)*hy
52    si(j)=sin(teta)

      i2dat='km2.dat'
      nt3dat='nt.dat'
      n2 = (il+2)*(jl+2)*Lfort
      n22= (il+4)*(jl+4)*2*Lfort
      n3 = n2*kl
      nm = il*jl*kl*Lfort
      Open ( Unit=49,File=i2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)
      Read(49,Rec=1) km2
      fileout= 'fwcontent.dat'
	open(unit=50, file=fileout)
      tdat= 't000000'

	write(*,*) 'Enter No of year start and stop-->'
	read(*,*) i0, i1

      DO iYear=i0,i1
	DO iMonth=1,12

      Write(tdat(2:5),'(i4.4)') iYear
      Write(tdat(6:7),'(i2.2)') iMonth

      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=2) s
      Close(32)

      FWC= 0.  

*     Sum over the volume
      do j=1,jl-1
      do i=1,il-1

      Do 33 k=1,kl-1
	
	depth=hz(k)

      IF(km2(i  ,j  ).GE. k+1 .and.
     *   km2(i+1,j  ).GE. k+1 .and. 
     *   km2(i  ,j+1).GE. k+1 .and. 
     *   km2(i+1,j+1).GE. k+1     ) THEN
      q= (s(i,j,k  )+2.*s(i+1,j,k  )+2.*s(i,j+1,k  )+s(i+1,j+1,k  ))/6.0
      qp=(s(i,j,k+1)+2.*s(i+1,j,k+1)+2.*s(i,j+1,k+1)+s(i+1,j+1,k+1))/6.0
	salt=0.5*(q+qp)
ccc      FWC= FWC+r2hy2*si(j)*depth*MAX(0.,(34.8-salt))/34.8
      FWC= FWC+r2hy2*si(j)*depth*(34.8-salt)/34.8
	END IF
33    CONTINUE ! K

      end do   ! i
	end do   ! j
*     End of integration

	time= REAL(iYear) +REAL(iMonth)/12. -1./24.
      write(50,*) time, FWC

	END DO   ! day
	END DO   ! year

      end
