      Subroutine MIXERKC

*     version 13.12.12.

      INCLUDE 'slo2.fi'
      dimension AM(kl),BM(kl),CM(kl),FM(kl),Rksi(kl)  !arrays for sweeping
	dimension cL(kl)

*     vertival mixing coefficients by Kantha & Clayson, 1994 -
*     version of Mellor and Yamada's level 2.5 model (Mellor and Yamada, 1982).
*     As decrribed in CLIO 3.0 manual.

*     References:
*     Mellor G.L. and T. Yamada (1982). Development of a turbulence closure 
*     model for geophysical fluid problems. Rev. Geophys. Spac. Phys. 20(4), 
*     851-875.
*     Kantha L.H. and C.A. Clayson (1994). An improved mixed layer model for 
*     geophysical applications. J. Geophys. Res. 99(C12), 25235-25266.
*     Galperin B., L.H. Kantha, S. Hassid and A. Rosati (1988). 
*     A quasi-equilibrium turbulent energy model for geophysical flows. 
*     J. Atmos. Sciences 45, 55-62.
*     Deleersnijder E. and P. Luyten (1994). On the practical advantages of
*     the quasiequilibrium version of the Mellor and Yamada level 2.5 
*     turbulence closure applied to marine modelling. Appl. Math. Modelling
*     18, 281-287

*     Surface winds:
*     Craig, P. D., and M. L. Banner, 1994: Modeling wave-enhanced turbulence
*     in the ocean surface layer. J. Phys. Oceanogr., 24, 2546–2559.

      parameter (cond=2)    ! B.c. 1- Dirichlet, 2 - Neumann.

      Parameter (Sq=0.2)    ! Mellor and Yamada, 1982

*     Original KC94
      Parameter
     *(A1=0.92,A2=0.74,B1=16.55,B2=10.1,C1=0.08,C2=0.7,C3=0.2)

*     Kantha, L. H., On an improved model for the turbulent pbl, J. Atmos. Sci., 
*     60, 2239–2246, 2003.
c      Parameter
c     *(A1=0.58,A2=0.62,B1=16.55,B2=11.6,C1=0.038,C2=0.7,C3=0.2)

	Parameter (C4=0.53)
	Parameter (cK=0.4)      ! von Karman constant
c	Parameter (cL0=2.e2)    ! macroscale for turbulence, set = 10m 
*	                        ! other proposals see in Mellor and Yamada, 1982
      Parameter (cLmin=1.e-1) ! Minimum length scale
	Parameter (B123=6.4943) ! B123= B1**(2/3) - for Dirichlet b.c.


      INCLUDE 'tparm.fi'

*     ---------------Constants  -------------
      GhL= 1.0/(A2*(B1+12.*A1+3.*B2*(1.-C3)))
      alpha = A2*(1.-6.*A1/B1)
	beta  = 3.*A2*(6.*A1 +B2*(1.-C3))
	gamma = A1*(1.-3.*C1-6.*A1/B1)
	delta = 9.*A1*(1.+2.*A1+A2-C2)
	epsil = 9.*A1*A2
	a_intern_waves = 1.0              ! 0.7 Mellor 1989, 08-1.4 Doronin 2000
*     --------------------------------------

      puny=1.e-12

      Rksi =0.

      do j=1,jl

	S0=Si(j)

	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then

*     Mean square variation of the ice thickness - to compute the length scale
	himean= 0.
	do mg=1,mgrad
	himean=himean+roi*hice(mg,i,j)+rosdry*hsnow(mg,i,j)
	enddo

      if(1.-aice(0,i,j) .GT. Aimin) then
      himean= himean/(1.-aice(0,i,j))
	else
	himean=0.
	end if

      sigmahice=0.
	do mg=1,mgrad
	sigmahice=sigmahice + 
     *Aice(mg,i,j)*(hice(mg,i,j)/max(Aimin,aice(mg,i,j)) -himean)**2
	enddo

	sigmahice= SQRT(sigmahice) +cLmin


*     ---------- Wind in [m/s] -----------
      wmod= SQRT( wx(i,j)**2 + wy(i,j)**2)
*     ------ Wind Stress Components -------
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      windx= cdrag*wx(i,j)
      windy= cdrag*wy(i,j)

      Aopen= MIN(Aice(0,i,j), 1.0)

*     Mean Drag coefficient
	drag=CDgwd(i,j)

	windx= Aopen*windx +drag*(1.-Aopen)*(uice(i,j)-u(i,j,1))
	windy= Aopen*windy +drag*(1.-Aopen)*(vice(i,j)-v(i,j,1))

	u_star= SQRT(sqrt(windx**2 + windy**2)/row)  ! Friction velocity

c      cL0=1.e3                                     ! CLIO 3.0 default
c	 cL0= 0.028*u_star/abs(2.*om*co(i)*S0)         ! McPhee book
	cL0= 0.7  *u_star/abs(2.*om*co(i)*S0)         ! Planetary boundary layer 

      do k=1,kb

*     Calculate master length scale

      cLs= max(sigmahice, z(k))          ! distance to the surface
	cLb= max(cLmin    , z(kb)-z(k))    ! distance to the bottom
	cLd= cLs*cLb/(cLs+cLb)
	
	cL(k)= cL0*cK*cLd/( cK*cLd + cL0) ! Mellor and Yamada, 1982

*     Limitation to avoid unphysical results

	IF(k.eq.1)then

	kp=2

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)

*      ro1 =Ro(i,j,kp)
*      ro2 =Ro(i,j,k )

      VB= g*(Ro1-Ro2)/hz(1)

	end if


	if (k.eq.kb) then

	km=k-1

	ppp= 1.e-5*g*row*(z(k)+z(km))/2.                   
      ro1= sigma_t(t(i,j,k ),s(i,j,k ), ppp)
      ro2= sigma_t(t(i,j,km),s(i,j,km), ppp)

*      ro1 =Ro(i,j,k )
*      ro2 =Ro(i,j,km)


      VB= g*(Ro1-Ro2)/hz(km)

	end if

	IF( k.GT.1 .and. k.LT.kb) then

	km=k-1
	kp=k+1

	ppp= 1.e-5*g*row*z(k)                   
      rokm12= 
     *sigma_t(0.5*(t(i,j,k)+t(i,j,km)),0.5*(s(i,j,k)+s(i,j,km)), ppp)

      rokp12= 
     *sigma_t(0.5*(t(i,j,k)+t(i,j,kp)),0.5*(s(i,j,k)+s(i,j,kp)), ppp)

	VB= 2.*g*(rokp12-rokm12)/(hz(k)+hz(km))

*      ro1 =Ro(i,j,kp)
*      ro2 =Ro(i,j,km)

	end if

	VB= VB/row

*     Limitation to avoid unphysical results
      IF(VB.GT.0.) then
      cL(k)=MIN(cL(k), c4*SQRT(q2turb(i,j,k)/max(puny, VB)))
	end if

*     Limitation to avoid low turbulence
      cL(k)=MAX(cL(k), cLmin)

      end do

*======================================================================

*                             Equation solution

*======================================================================

	do k=1,kb-1

	kp=k+1
	hzk =hz(k)
	n =abs(nt3(i,j,k ))
	np=abs(nt3(i,j,kp))

      IF(k.EQ.1) THEN
        hzk1 = 0.
      ELSE
	  hzk1 = hz(k-1)
      END IF

      CA=2.*DT/(HZK1+HZK)

*     Bouyoancy fluxes
      
	IF(k.eq.1)then

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)


*      ro1 =Ro(i,j,kp)
*      ro2 =Ro(i,j,k )

      VB= g*(Ro1-Ro2)/hz(1)/row

	umk=KT(1,n)*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))+
     *    KT(2,n)*(u(i,j,k)+u(i+1,j-1,k)+u(i,j-1,k))+
     *	KT(3,n)*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k))+
     *	KT(4,n)*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k))+
     *	KT(5,n)*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k))+
     *	KT(6,n)*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))

	umkp=KT(1,n)*(u(i,j,kp)+u(i+1,j,kp)+u(i+1,j-1,kp))+
     *     KT(2,n)*(u(i,j,kp)+u(i+1,j-1,kp)+u(i,j-1,kp))+
     *	 KT(3,n)*(u(i,j,kp)+u(i-1,j,kp)+u(i,j-1,kp))+
     *	 KT(4,n)*(u(i,j,kp)+u(i-1,j,kp)+u(i-1,j+1,kp))+
     *	 KT(5,n)*(u(i,j,kp)+u(i-1,j+1,kp)+u(i,j+1,kp))+
     *	 KT(6,n)*(u(i,j,kp)+u(i+1,j,kp)+u(i,j+1,kp))
	
	vmk=KT(1,n)*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))+
     *    KT(2,n)*(v(i,j,k)+v(i+1,j-1,k)+v(i,j-1,k))+
     *	KT(3,n)*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k))+
     *	KT(4,n)*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k))+
     *	KT(5,n)*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k))+
     *	KT(6,n)*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))

	vmkp=KT(1,n)*(v(i,j,kp)+v(i+1,j,kp)+v(i+1,j-1,kp))+
     *     KT(2,n)*(v(i,j,kp)+v(i+1,j-1,kp)+v(i,j-1,kp))+
     *	 KT(3,n)*(v(i,j,kp)+v(i-1,j,kp)+v(i,j-1,kp))+
     *	 KT(4,n)*(v(i,j,kp)+v(i-1,j,kp)+v(i-1,j+1,kp))+
     *	 KT(5,n)*(v(i,j,kp)+v(i-1,j+1,kp)+v(i,j+1,kp))+
     *	 KT(6,n)*(v(i,j,kp)+v(i+1,j,kp)+v(i,j+1,kp))

      uz= (umkp-umk)/hzk/cg(n)/3.
      vz= (vmkp-vmk)/hzk/cg(n)/3.

	else

	km=k-1

	ppp= 1.e-5*g*row*z(k)                   
      rokm12= 
     *sigma_t(0.5*(t(i,j,k)+t(i,j,km)),0.5*(s(i,j,k)+s(i,j,km)), ppp)

      rokp12= 
     *sigma_t(0.5*(t(i,j,k)+t(i,j,kp)),0.5*(s(i,j,k)+s(i,j,kp)), ppp)

	VB= 2.*g*(rokp12-rokm12)/(hz(k)+hz(km))/row


	umk=KT(1,n)*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))+
     *    KT(2,n)*(u(i,j,k)+u(i+1,j-1,k)+u(i,j-1,k))+
     *	KT(3,n)*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k))+
     *	KT(4,n)*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k))+
     *	KT(5,n)*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k))+
     *	KT(6,n)*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))
	umkm=KT(1,n)*(u(i,j,km)+u(i+1,j,km)+u(i+1,j-1,km))+
     *    KT(2,n)*(u(i,j,km)+u(i+1,j-1,km)+u(i,j-1,km))+
     *	KT(3,n)*(u(i,j,km)+u(i-1,j,km)+u(i,j-1,km))+
     *	KT(4,n)*(u(i,j,km)+u(i-1,j,km)+u(i-1,j+1,km))+
     *	KT(5,n)*(u(i,j,km)+u(i-1,j+1,km)+u(i,j+1,km))+
     *	KT(6,n)*(u(i,j,km)+u(i+1,j,km)+u(i,j+1,km))

	umkpp=KT(1,np)*(u(i,j,kp)+u(i+1,j,kp)+u(i+1,j-1,kp))+
     *     KT(2,np)*(u(i,j,kp)+u(i+1,j-1,kp)+u(i,j-1,kp))+
     *	 KT(3,np)*(u(i,j,kp)+u(i-1,j,kp)+u(i,j-1,kp))+
     *	 KT(4,np)*(u(i,j,kp)+u(i-1,j,kp)+u(i-1,j+1,kp))+
     *	 KT(5,np)*(u(i,j,kp)+u(i-1,j+1,kp)+u(i,j+1,kp))+
     *	 KT(6,np)*(u(i,j,kp)+u(i+1,j,kp)+u(i,j+1,kp))
	umkpk=KT(1,np)*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))+
     *     KT(2,np)*(u(i,j,k)+u(i+1,j-1,k)+u(i,j-1,k))+
     *	 KT(3,np)*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k))+
     *	 KT(4,np)*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k))+
     *	 KT(5,np)*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k))+
     *	 KT(6,np)*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))
	
	vmk=KT(1,n)*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))+
     *    KT(2,n)*(v(i,j,k)+v(i+1,j-1,k)+v(i,j-1,k))+
     *	KT(3,n)*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k))+
     *	KT(4,n)*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k))+
     *	KT(5,n)*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k))+
     *	KT(6,n)*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))
	vmkm=KT(1,n)*(v(i,j,km)+v(i+1,j,km)+v(i+1,j-1,km))+
     *    KT(2,n)*(v(i,j,km)+v(i+1,j-1,km)+v(i,j-1,km))+
     *	KT(3,n)*(v(i,j,km)+v(i-1,j,km)+v(i,j-1,km))+
     *	KT(4,n)*(v(i,j,km)+v(i-1,j,km)+v(i-1,j+1,km))+
     *	KT(5,n)*(v(i,j,km)+v(i-1,j+1,km)+v(i,j+1,km))+
     *	KT(6,n)*(v(i,j,km)+v(i+1,j,km)+v(i,j+1,km))

	vmkpp=KT(1,np)*(v(i,j,kp)+v(i+1,j,kp)+v(i+1,j-1,kp))+
     *     KT(2,np)*(v(i,j,kp)+v(i+1,j-1,kp)+v(i,j-1,kp))+
     *	 KT(3,np)*(v(i,j,kp)+v(i-1,j,kp)+v(i,j-1,kp))+
     *	 KT(4,np)*(v(i,j,kp)+v(i-1,j,kp)+v(i-1,j+1,kp))+
     *	 KT(5,np)*(v(i,j,kp)+v(i-1,j+1,kp)+v(i,j+1,kp))+
     *	 KT(6,np)*(v(i,j,kp)+v(i+1,j,kp)+v(i,j+1,kp))
	vmkpk=KT(1,np)*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))+
     *     KT(2,np)*(v(i,j,k)+v(i+1,j-1,k)+v(i,j-1,k))+
     *	 KT(3,np)*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k))+
     *	 KT(4,np)*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k))+
     *	 KT(5,np)*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k))+
     *	 KT(6,np)*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))

      uz= ((umkpp-umkpk)+(umk-umkm)) /(cg(np)*hzk+cg(n)*hzk1)/3.
      vz= ((vmkpp-vmkpk)+(vmk-vmkm)) /(cg(np)*hzk+cg(n)*hzk1)/3.

	end if

*     Vertical shift **2 in velocity points
      VS=uz*uz + vz*vz

*     coefficients and right hand side for vertical sweeping.

      IF( k.EQ.1) THEN

*     ---------- Wind in [m/s] -----------
      wmod= SQRT( wx(i,j)**2 + wy(i,j)**2)
*     ------ Wind Stress Components -------
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      windx= cdrag*wx(i,j)
      windy= cdrag*wy(i,j)

	u_star= SQRT(sqrt(windx**2 + windy**2)/row)  !!! Friction velocity
      Aopen= MIN(Aice(0,i,j), 1.0)

	Prod= 2.*dt*az (i,j,1)*VS
	Buoy= 2.*dt*azt(i,j,1)*VB *a_intern_waves
      Diss= 2.*dt*sqrt(q2turb(i,j,1))/(B1*cL(1))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,1)+q2turb(i,j,2)) 

      AM(1)=0.
      CM(1)=-0.5*CA*(cL(1)+cL(2))*Sq*sqrt(turb1)/hzk
      BM(1)=1.0 -CM(1)     + Pminus  
      FM(1)= q2turb(i,j,1) + Pplus
     +     + CA*Aopen*dt*100.*u_star**3  !! Surface waves

      ELSE ! k>1

	Prod= dt*(az (i,j,k-1)+az (i,j,k))*VS
	Buoy= dt*(azt(i,j,k-1)+azt(i,j,k))*VB *a_intern_waves
      Diss= 2.*dt*sqrt(q2turb(i,j,k))/(B1*cL(k))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,k)+q2turb(i,j,k+1)) 
	turb2 =0.5*(q2turb(i,j,k-1)+q2turb(i,j,k)) 

      AM(k)=-0.5*CA*(cL(k-1)+cL(k  ))*Sq*sqrt(turb2)/hzk1
      CM(k)=-0.5*CA*(cL(k  )+cL(k+1))*Sq*sqrt(turb1) /hzk
      BM(k)= 1.-CM(k)-AM(k) +Pminus
      FM(k)= q2turb(i,j,k)  +Pplus
      
	
	END IF ! K>1

      end do ! k

*     Bottom

      k=kb
	km=kb-1
      n=ABS(nt3(i,j,k))

	ppp= 1.e-5*g*row*(z(k)+z(km))/2.                   
      ro1= sigma_t(t(i,j,k ),s(i,j,k ), ppp)
      ro2= sigma_t(t(i,j,km),s(i,j,km), ppp)


*      ro1 =Ro(i,j,k)
*      ro2 =Ro(i,j,km)

      VB= g*(Ro1-Ro2)/hz(km)/row

	umk=KT(1,n)*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))+
     *    KT(2,n)*(u(i,j,k)+u(i+1,j-1,k)+u(i,j-1,k))+
     *	KT(3,n)*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k))+
     *	KT(4,n)*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k))+
     *	KT(5,n)*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k))+
     *	KT(6,n)*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))
	umkm=KT(1,n)*(u(i,j,km)+u(i+1,j,km)+u(i+1,j-1,km))+
     *    KT(2,n)*(u(i,j,km)+u(i+1,j-1,km)+u(i,j-1,km))+
     *	KT(3,n)*(u(i,j,km)+u(i-1,j,km)+u(i,j-1,km))+
     *	KT(4,n)*(u(i,j,km)+u(i-1,j,km)+u(i-1,j+1,km))+
     *	KT(5,n)*(u(i,j,km)+u(i-1,j+1,km)+u(i,j+1,km))+
     *	KT(6,n)*(u(i,j,km)+u(i+1,j,km)+u(i,j+1,km))

	vmk=KT(1,n)*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))+
     *    KT(2,n)*(v(i,j,k)+v(i+1,j-1,k)+v(i,j-1,k))+
     *	KT(3,n)*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k))+
     *	KT(4,n)*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k))+
     *	KT(5,n)*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k))+
     *	KT(6,n)*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))
	vmkm=KT(1,n)*(v(i,j,km)+v(i+1,j,km)+v(i+1,j-1,km))+
     *    KT(2,n)*(v(i,j,km)+v(i+1,j-1,km)+v(i,j-1,km))+
     *	KT(3,n)*(v(i,j,km)+v(i-1,j,km)+v(i,j-1,km))+
     *	KT(4,n)*(v(i,j,km)+v(i-1,j,km)+v(i-1,j+1,km))+
     *	KT(5,n)*(v(i,j,km)+v(i-1,j+1,km)+v(i,j+1,km))+
     *	KT(6,n)*(v(i,j,km)+v(i+1,j,km)+v(i,j+1,km))


      uz= (umk-umkm)/hz(km)/cg(n)/3.
      vz= (vmk-vmkm)/hz(km)/cg(n)/3.


*     Vertical shift **2 in velocity points
      VS=uz*uz + vz*vz

      CA=2.*DT/HZ(Kb-1)

	turb2 =0.5*(q2turb(i,j,kb-1)+q2turb(i,j,kb)) 

	Prod= 2.*dt*az (i,j,kb-1)*VS
	Buoy= 2.*dt*azt(i,j,kb-1)*VB *a_intern_waves
      Diss= 2.*dt*sqrt(q2turb(i,j,kb))/(B1*cL(kb))

	Pplus = Prod - Buoy
	Pminus= Diss

      AM(kb)=-0.5*CA*(cL(kb)+cL(kb-1))*Sq*sqrt(turb2)/hz(kb-1)
      CM(kb)= 0.
      BM(kb)= 1.-AM(kb)      +Pminus
      FM(kb)= q2turb(i,j,kb) +Pplus

      call FACTOR(KL,am,bm,cm,fm,rksi,1,kb)

      do k=1,kb
      q2turb(i,j,k)= MAX(0.,rksi(k))
      end do


      do k=1,kb-1

	kp=k+1

*     projection in the middle of the box in coefficients points
	turb =0.5*(q2turb(i,j,k)+q2turb(i,j,kp)) 

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)


*      ro1 =Ro(i,j,kp)
*      ro2 =Ro(i,j,k )

      VB= g*(Ro1-Ro2)/hz(k)/row


*     Following the arguments of Galperin et al. (1988), the stability functions
*     Su and Ss are only a function of Gh, which renders the model less prone to
*     oscillations than its original version, in which the stability functions 
*     also depend on the current shearing (e.g., Deleersnijder and Luyten,1994).

      cLz= 0.5*(cL(k)+cL(kp))

	Gh= -VB*(cLz**2)/max(puny,turb)   ! Richardson number

*     Limitation to avoid unphysical results 
      Gh= MIN(Gh,GhL)
	Gh= MAX(-0.28,Gh)

*     stability functions SFu, SFs, Kantha and Clayson (1994)

	SFs = alpha/(1.- beta*Gh)
      SFu = (gamma + delta*SFs*Gh)/(1.- epsil*Gh)

*     New coefficients

*     Here should be the layer beneath the ML
c      if(turb. GT. 2.e-2) then
      
      Az (i,j,k)= cLz*SFu*SQRT(turb) + AZbg
      AzT(i,j,k)= cLz*SFs*SQRT(turb) + AZTbg
      AzS(i,j,k)= cLz*SFs*SQRT(turb) + AZTbg


c	else ! Kantha and Clayson (1994), see also book 2000

*     Shear instability
*     Vertical shift
c      uz= u(i,j,k+1)-u(i,j,k)
c      vz= v(i,j,k+1)-v(i,j,k)
c      Ri= VB*hz(k)**2/MAX(puny,uz**2+vz**2)

c	if(Ri .GE. 0.) then    ! stable layer

c	if(Ri.LE.0.7) then      
c      add= 50.*(1.-(Ri/0.7)**2)**3
c      Az (i,j,k)= add + AZbg
c      AzT(i,j,k)= add + AZTbg
c      AzS(i,j,k)= add + AZTbg  
c	end if

c	else ! unstable layer
c      Az (i,j,k)= 50.+ AZbg
c      AzT(i,j,k)= 50.+ AZTbg
c      AzS(i,j,k)= 50.+ AZTbg  
c	end if                 ! stability criterium

c	end if  ! turb > 2e-2

*     Double diffusion

c	Tk12= 0.5*(T(i,j,k)+T(i,j,kp))
c	Sk12= 0.5*(S(i,j,k)+S(i,j,kp))
c	Pk12= 0.5e-5*g*row*(z(k)+z(kp)) !!!+10.*PA(i,j)


c      RoS=  sigma_t(Tk12,S(i,j,kp),Pk12)-sigma_t(Tk12,S(i,j,k),Pk12)
c      RoT=  sigma_t(T(i,j,kp),Sk12,Pk12)-sigma_t(T(i,j,k),Sk12,Pk12)

c	if(RoT .GE. 0.) then
c	R_Ro= RoS/max( puny,RoT)
c	else
c	R_Ro= RoS/min(-puny,RoT)
c	end if

c	IF(R_Ro .LE. 1.9 .AND. R_Ro .GE. 1.0) then
c	add = 50. * (1. - ((R_Ro-1.)/0.9)**2)**3
c      AzT(i,j,k)= AzT(i,j,k) + 0.7*add
c      AzS(i,j,k)= AzS(i,j,k) +     add
c	END IF

	end do ! k

	end if ! km2 >0
	end do ! i
	end do ! j

	return
	end

	subroutine TS_convection
      INCLUDE 'slo2.fi'
      INCLUDE 'tparm.fi'

*     version 28.10.12

*     Check instability
      do j=1,jl
	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then

	do iter=1,5

      do k=1,kb-1

	kp=k+1

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)

cccc      if( Ropot(i,j,k+1) .LT. Ropot(i,j,k)) then
      if( Ro1 .LT. Ro2 ) then

	write(*,*) 'Instability:',i,j,k, Ro1, Ro2

	if(k.eq.1) then
	Tmean = hz(k)*T(i,j,k) + (hz(k)+hz(k+1))*T(i,j,k+1)
	Tmean= Tmean /(2.*hz(k) +hz(k+1))
	end if

      if(k.GT.1 .and. k. LT. kb-1) then
	Tmean = (hz(k)+hz(k-1))*T(i,j,k) + (hz(k)+hz(k+1))*T(i,j,k+1)
	Tmean= Tmean /(hz(k-1) +2.*hz(k) +hz(k+1))
	end if

      if(k .eq. kb-1) then
	Tmean = (hz(k-1)+hz(k))*T(i,j,k) + hz(k)*T(i,j,k+1)
	Tmean= Tmean /(2.*hz(k) +hz(k-1))
	end if

	T(i,j,k  )= Tmean
	T(i,j,k+1)= Tmean


	if(k.eq.1) then
	Smean = hz(k)*S(i,j,k) + (hz(k)+hz(k+1))*S(i,j,k+1)
	Smean= Smean /(2.*hz(k) +hz(k+1))
	end if

      if(k.GT.1 .and. k. LT. kb-1) then
	Smean = (hz(k)+hz(k-1))*S(i,j,k) + (hz(k)+hz(k+1))*S(i,j,k+1)
	Smean= Smean /(hz(k-1) +2.*hz(k) +hz(k+1))
	end if

      if(k .eq. kb-1) then
	Smean = (hz(k-1)+hz(k))*S(i,j,k) + hz(k)*S(i,j,k+1)
	Smean= Smean /(2.*hz(k) +hz(k-1))
	end if

	S(i,j,k  )= Smean
	S(i,j,k+1)= Smean

	ppp= 1.e-5*g*row*z(k)                   
      ro(i,j,k)= sigma_t(t(i,j,k),s(i,j,k), ppp)
	ppp= 1.e-5*g*row*z(k+1)                   
      ro(i,j,k+1)= sigma_t(t(i,j,k+1),s(i,j,k+1), ppp)

      ropot(i,j,k)= sigma_t(t(i,j,k),s(i,j,k), 0.)
      ropot(i,j,k+1)= ropot(i,j,k)

	end if   ! Ropot is unstable

      end do   ! k

	end do   ! iter

	end if   ! kb>0

	end do   ! i
	end do   ! j


	return
	end

