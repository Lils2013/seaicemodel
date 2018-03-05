      program LengthScale
	Parameter (cK=0.4)    ! von Karman constant
	Parameter (cL0=1.e3)  ! macroscale for turbulence, set = 10m
      Parameter (cLmin=1.e2)
	dimension z(41), cL(40)
	    
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

      do 33 k=1,41
33    z(k)=100.*z(k)
	  
	do k=1,40
      cLs= max(cLmin, z(k))         ! distance to the surface
	cLb= max(cLmin, z(40)-z(k))   ! distance to the bottom
	cLd= cLs*cLb/(cLs+cLb)
	
	cL(k)= cL0*cK*cLd/( cK*cLd + cL0) ! Mellor and Yamada, 1982

	end do

	write(*,*) cL

	end
