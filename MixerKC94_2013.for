      Subroutine MIXERKC

*     version 16.04.13.

      INCLUDE 'slo2.fi'
      dimension AM(kl),BM(kl),CM(kl),FM(kl),Rksi(kl)  !arrays for sweeping
	dimension cL(kl)

cc	real*8 ro1,ro2 

*     vertival mixing coefficients by Kantha & Clayson, 1994 -
*     version of Mellor and Yamada's level 2.5 model (Mellor and Yamada, 1982).
*     As decrribed in CLIO 3.0 manual.

*     Also modifications for boundary conditions by ROMS

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
      GhL    = 1.0/(A2*(B1+12.*A1+3.*B2*(1.-C3)))
      alpha  = A2*(1.-6.*A1/B1)
	beta   = 3.*A2*(6.*A1 +B2*(1.-C3))
	gamma  = A1*(1.-3.*C1-6.*A1/B1)
	delta  = 9.*A1*(1.+2.*A1+A2-C2)
	epsil  = 9.*A1*A2
	a_intern_waves = 1.0 !!!0.7     ! 0.7 Mellor 1989, 0.8-1.4 Doronin 2000
	AZS_to_AZT     = 1.0
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

	u_star= SQRT(sqrt(windx**2 + windy**2)/row)  !!! Friction velocity

c	cL0= 0.028*u_star/abs(2.*om*co(i)*S0)         ! McPhee book
c	cL0= 0.7  *u_star/abs(2.*om*co(i)*S0)         ! Ekman layer depth
      cL0=1.e3                                      ! CLIO 3.0 default

      do k=1,kb

*     Calculate master length scale

      cLs= max(sigmahice, z(k))          ! distance to the surface
	cLb= max(cLmin    , z(kb)-z(k))    ! distance to the bottom
	cLd= cLs*cLb/(cLs+cLb)
	
	cL(k)= cL0*cK*cLd/( cK*cLd + cL0) ! Mellor and Yamada, 1982

*     Limitation to avoid unphysical results

	IF(k.eq.1)then

	kp=2


c	ppp= 5.e-6*g*row*(z(k)+z(kp))
c      ro1=sigma_t(t(i,j,kp),s(i,j,kp), ppp)
c      ro2=sigma_t(t(i,j,k ),s(i,j,k ), ppp)
c      VB= 1.e-3*g*REAL(Ro1-Ro2)/hz(1)/row

      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,k )
	call sound(T(i,j,k),S(i,j,k),PPP,Ro1,csound)
      VB= g*REAL(Ro1-Ro2)/hz(1)/row -(g/csound)**2

	end if

	if (k.eq.kb) then

	km=k-1


c	ppp= 5.e-6*g*row*(z(k)+z(km))
c      ro1=sigma_t(t(i,j,k ),s(i,j,k ), ppp)
c      ro2=sigma_t(t(i,j,km),s(i,j,km), ppp)
c      VB= 1.e-3*g*REAL(Ro1-Ro2)/hz(km)/row


      ro1 =Ro(i,j,k )
      ro2 =Ro(i,j,km)
	call sound(T(i,j,k),S(i,j,k),PPP,Ro1,csound)
      VB= g*REAL(Ro1-Ro2)/hz(km)/row -(g/csound)**2
	
	end if

	IF( k.GT.1 .and. k.LT.kb) then

	km=k-1
	kp=k+1

c	ppp= 1.e-5*g*row*z(k)
c      t12m= 0.5*(t(i,j,km)+t(i,j,k))
c	s12m= 0.5*(s(i,j,km)+s(i,j,k))
c      t12p= 0.5*(t(i,j,k)+t(i,j,kp))
c	s12p= 0.5*(s(i,j,k)+s(i,j,kp))
c      ro1=sigma_t(t12p,s12p, ppp)
c      ro2=sigma_t(t12m,s12m, ppp)
c	VB= 2.e-3*g*REAL(ro1-ro2)/(hz(k)+hz(km))/row

      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,km)
	call sound(T(i,j,k),S(i,j,k),PPP,Ro(i,j,k),csound)
	VB= g*REAL(ro1-ro2)/(hz(k)+hz(km))/row -(g/csound)**2

	end if


*     Limitation to avoid unphysical results
      IF(VB.GT.0.) then
      cL(k)=MIN(cL(k), C4*SQRT(q2turb(i,j,k)/max(puny, VB)))
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

      CA=2.*DT/(CG(n)*HZK1+CG(np)*HZK)

*     Bouyoancy fluxes
      
	IF(k.eq.1)then


c	ppp= 5.e-6*g*row*(z(k)+z(kp))
c      ro1=sigma_t(t(i,j,kp),s(i,j,kp), ppp)
c      ro2=sigma_t(t(i,j,k ),s(i,j,k ), ppp)
c      VB= 1.e-3*g*REAL(Ro1-Ro2)/hz(1)/row

      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,k )
      call sound(T(i,j,k),S(i,j,k),PPP,Ro1,csound)
      VB= g*REAL(Ro1-Ro2)/hzk/row -(g/csound)**2

	umk =u(i,j,kp)
	umkm=u(i,j,k)
	vmk =v(i,j,kp)
	vmkm=v(i,j,k)

      uz= (umk-umkm)/hzk
      vz= (vmk-vmkm)/hzk

	else

	km=k-1

c	ppp= 1.e-5*g*row*z(k)
c      t12m= 0.5*(t(i,j,km)+t(i,j,k))
c	s12m= 0.5*(s(i,j,km)+s(i,j,k))
c      t12p= 0.5*(t(i,j,k)+t(i,j,kp))
c	s12p= 0.5*(s(i,j,k)+s(i,j,kp))
c      ro1=sigma_t(t12p,s12p, ppp)
c      ro2=sigma_t(t12m,s12m, ppp)
c	VB= 2.e-3*g*REAL(ro1-ro2)/(hz(k)+hz(km))/row

      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,km)
	call sound(T(i,j,k),S(i,j,k),PPP,Ro(i,j,k),csound)
	VB= g*REAL(ro1-ro2)/(hz(k)+hz(km))/row -(g/csound)**2

	umk =u(i,j,kp)
	umkm=u(i,j,km)
	vmk =v(i,j,kp)
	vmkm=v(i,j,km)

      uz= (umk-umkm)/(hz(k)+hz(km))
      vz= (vmk-vmkm)/(hz(k)+hz(km))


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

*     Friction velocity open water
	u_star= SQRT(sqrt(windx**2 + windy**2)/row) 
      Aopen= MIN(Aice(0,i,j), 1.0)

	Prod= 2.*dt*az (i,j,1)*VS
	Buoy= 2.*dt*azt(i,j,1)*VB*a_intern_waves *AZS_to_AZT
      Diss= 2.*dt*sqrt(q2turb(i,j,1))/(B1*cL(1))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,1)+q2turb(i,j,2)) 

      AM(1)=0.
      CM(1)=-0.5*CA*CG(n)*(cL(1)+cL(2))*Sq*sqrt(turb1)/hzk
      BM(1)=1.0 -CM(1)     + Pminus  
      FM(1)= q2turb(i,j,1) + Pplus
     +     + CA*Aopen*dt*100.*u_star**3  ! Surface waves

      ELSE ! k>1

	Prod= dt*(az (i,j,k-1)+az (i,j,k))*VS
	Buoy= dt*(azt(i,j,k-1)+azt(i,j,k))*VB*a_intern_waves *AZS_to_AZT
      Diss= 2.*dt*sqrt(q2turb(i,j,k))/(B1*cL(k))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,k)+q2turb(i,j,k+1)) 
	turb2 =0.5*(q2turb(i,j,k-1)+q2turb(i,j,k)) 

      AM(k)=-0.5*CA*CG(n )*(cL(k-1)+cL(k  ))*Sq*sqrt(turb2)/hzk1
      CM(k)=-0.5*CA*CG(np)*(cL(k  )+cL(k+1))*Sq*sqrt(turb1) /hzk
      BM(k)= 1.-CM(k)-AM(k) +Pminus
      FM(k)= q2turb(i,j,k)  +Pplus
      
	
	END IF ! K>1

      end do ! k

*     Bottom

      k=kb
	km=kb-1
      n=ABS(nt3(i,j,k))

c	ppp= 5.e-6*g*row*(z(k)+z(km))
c      ro1=sigma_t(t(i,j,k ),s(i,j,k ), ppp)
c      ro2=sigma_t(t(i,j,km),s(i,j,km), ppp)
c      VB= 1.e-3*g*REAL(Ro1-Ro2)/hz(km)/row

      ro1 =Ro(i,j,k)
      ro2 =Ro(i,j,km)
	call sound(T(i,j,k),S(i,j,k),PPP,Ro1,csound)
      VB= g*REAL(Ro1-Ro2)/hz(km)/row - (g/csound)**2

	umk =u(i,j,k)
	umkm=u(i,j,km)
	vmk =v(i,j,k)
	vmkm=v(i,j,km)

      uz= (umk-umkm)/hz(km)
      vz= (vmk-vmkm)/hz(km)


*     Vertical shift **2 in velocity points
      VS=uz*uz + vz*vz

      CA=2.*DT/HZ(Kb-1)

	turb2 =0.5*(q2turb(i,j,kb-1)+q2turb(i,j,kb)) 

	Prod= 2.*dt*az (i,j,kb-1)*VS
	Buoy= 2.*dt*azt(i,j,kb-1)*VB*a_intern_waves *AZS_to_AZT
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
	hzk=hz(k)

	n =abs(nt3(i,j,kp))

*     projection in the middle of the box in coefficients points
	turb =0.5*(q2turb(i,j,k)+q2turb(i,j,kp)) 

c	ppp= 5.e-6*g*row*(z(k)+z(kp))
c      ro1=sigma_t(t(i,j,kp),s(i,j,kp), ppp)
c      ro2=sigma_t(t(i,j,k ),s(i,j,k ), ppp)
c      VB= 1.e-3*g*REAL(Ro1-Ro2)/hzk/row

      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,k )
	call sound(T(i,j,k),S(i,j,k),PPP,0.5*(ro1+ro2),csound)

      VB= g*REAL(Ro1-Ro2)/hzk/row -(g/csound)**2


*     Following the arguments of Galperin et al. (1988), the stability functions
*     Su and Ss are only a function of Gh, which renders the model less prone to
*     oscillations than its original version, in which the stability functions 
*     also depend on the current shearing (e.g., Deleersnijder and Luyten,1994).

      cLz= 0.5*(cL(k)+cL(kp))

	Gh= -VB*(cLz**2)/max(puny,turb)   ! Richardson number

*     Limitation to avoid unphysical results 
      Gh= MIN(Gh,GhL)
cccc	Gh= MAX(-0.0288,Gh)

*     stability functions SFu, SFs, Kantha and Clayson (1994)

	SFs = alpha/(1.- beta*Gh)
      SFu = (gamma + delta*SFs*Gh)/(1.- epsil*Gh)

*     ----------------  New coefficients  ---------------------

*     Here should be the layer beneath the ML
c      if(turb. GT. 2.e-2) then
      

      Az (i,j,k)= cLz*SFu*SQRT(turb) + AZbg
      AzT(i,j,k)= cLz*SFs*SQRT(turb) + AZTbg
      AzS(i,j,k)= AZS_to_AZT * (cLz*SFs*SQRT(turb) + AZTbg)

c	if(k.LE.7) then
c      Az (i,j,k)= 1.e5
c      AzT(i,j,k)= 1.e5
c      AzS(i,j,k)= 1.e5
c      end if



*     ---------------------------------------------------------

c	else ! Kantha and Clayson (1994), see also their book 2000

*     Shear instability
*     Vertical shift
c      uz= u(i,j,k+1)-u(i,j,k)
c      vz= v(i,j,k+1)-v(i,j,k)

c	umk=KT(1,n)*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))+
c     *    KT(2,n)*(u(i,j,k)+u(i+1,j-1,k)+u(i,j-1,k))+
c     *	KT(3,n)*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k))+
c     *	KT(4,n)*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k))+
c     *	KT(5,n)*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k))+
c     *	KT(6,n)*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))

c	umkp=KT(1,n)*(u(i,j,kp)+u(i+1,j,kp)+u(i+1,j-1,kp))+
c     *     KT(2,n)*(u(i,j,kp)+u(i+1,j-1,kp)+u(i,j-1,kp))+
c     *	 KT(3,n)*(u(i,j,kp)+u(i-1,j,kp)+u(i,j-1,kp))+
c     *	 KT(4,n)*(u(i,j,kp)+u(i-1,j,kp)+u(i-1,j+1,kp))+
c     *	 KT(5,n)*(u(i,j,kp)+u(i-1,j+1,kp)+u(i,j+1,kp))+
c     *	 KT(6,n)*(u(i,j,kp)+u(i+1,j,kp)+u(i,j+1,kp))
	
c	vmk=KT(1,n)*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))+
c     *    KT(2,n)*(v(i,j,k)+v(i+1,j-1,k)+v(i,j-1,k))+
c     *	KT(3,n)*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k))+
c     *	KT(4,n)*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k))+
c     *	KT(5,n)*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k))+
c     *	KT(6,n)*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))

c	vmkp=KT(1,n)*(v(i,j,kp)+v(i+1,j,kp)+v(i+1,j-1,kp))+
c     *     KT(2,n)*(v(i,j,kp)+v(i+1,j-1,kp)+v(i,j-1,kp))+
c     *	 KT(3,n)*(v(i,j,kp)+v(i-1,j,kp)+v(i,j-1,kp))+
c     *	 KT(4,n)*(v(i,j,kp)+v(i-1,j,kp)+v(i-1,j+1,kp))+
c     *	 KT(5,n)*(v(i,j,kp)+v(i-1,j+1,kp)+v(i,j+1,kp))+
c     *	 KT(6,n)*(v(i,j,kp)+v(i+1,j,kp)+v(i,j+1,kp))

c      uz= (umkp-umk)/hzk/cg(n)/3.
c      vz= (vmkp-vmk)/hzk/cg(n)/3.


c      Ri= VB*hzk**2/MAX(puny,uz**2+vz**2)

c	if(Ri .GE. 0.) then    ! stable layer

c	if(Ri .GE. 0. .and. Ri.LE.0.7) then      
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

*     Double diffusion - requires special treatment of bouyancy frequency!

c	Tk12= 0.5*(T(i,j,k)+T(i,j,kp))
c	Sk12= 0.5*(S(i,j,k)+S(i,j,kp))
c	Pk12= 0.5e-5*g*row*(z(k)+z(kp))


c      RoS=  sigma_t(Tk12,S(i,j,kp),Pk12)-sigma_t(Tk12,S(i,j,k),Pk12)
c      RoT=  sigma_t(T(i,j,kp),Sk12,Pk12)-sigma_t(T(i,j,k),Sk12,Pk12)

c	if(RoT .GE. 0.) then
c	R_Ro= RoS/max( puny,RoT)
c	else
c	R_Ro= RoS/min(-puny,RoT)
c	end if

c	IF(R_Ro .LE. 1.9 .AND. R_Ro .GE. 1.0) then

c	write(*,*) 'Double diffusion happened! i,j,k,R_ro,add', 
c     *                    i,j,k,R_Ro,add

c	add = 50. * (1. - ((R_Ro-1.)/0.9)**2)**3
c      AzT(i,j,k)= AzT(i,j,k) + 0.7*add
c      AzS(i,j,k)= AzS(i,j,k) +     add

c	END IF  ! R_Ro bounds

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

ccc	write(*,*) 'Instability:',i,j,k, Ro1, Ro2

	if(k.eq.1) then
	Tmean = hz(k)*T(i,j,k) + (hz(k)+hz(kp))*T(i,j,kp)
	Tmean= Tmean /(2.*hz(k) +hz(kp))
	end if

      if(k.GT.1 .and. k. LT. kb-1) then
	Tmean = (hz(k)+hz(k-1))*T(i,j,k) + (hz(k)+hz(kp))*T(i,j,kp)
	Tmean= Tmean /(hz(k-1) +2.*hz(k) +hz(kp))
	end if

      if(k .eq. kb-1) then
	Tmean = (hz(k-1)+hz(k))*T(i,j,k) + hz(k)*T(i,j,kp)
	Tmean= Tmean /(2.*hz(k) +hz(k-1))
	end if

	T(i,j,k )= Tmean
	T(i,j,kp)= Tmean


	if(k.eq.1) then
	Smean = hz(k)*S(i,j,k) + (hz(k)+hz(kp))*S(i,j,kp)
	Smean= Smean /(2.*hz(k) +hz(kp))
	end if

      if(k.GT.1 .and. k. LT. kb-1) then
	Smean = (hz(k)+hz(k-1))*S(i,j,k) + (hz(k)+hz(kp))*S(i,j,kp)
	Smean= Smean /(hz(k-1) +2.*hz(k) +hz(kp))
	end if

      if(k .eq. kb-1) then
	Smean = (hz(k-1)+hz(k))*S(i,j,k) + hz(k)*S(i,j,kp)
	Smean= Smean /(2.*hz(k) +hz(k-1))
	end if

	S(i,j,k )= Smean
	S(i,j,kp)= Smean

c	ppp= 1.e-5*g*row*z(k)                   
c      ro(i,j,k)= sigma_t(t(i,j,k),s(i,j,k), ppp)
c	ppp= 1.e-5*g*row*z(k+1)                   
c      ro(i,j,k+1)= sigma_t(t(i,j,k+1),s(i,j,k+1), ppp)

c      ropot(i,j,k)= sigma_t(t(i,j,k),s(i,j,k), 0.)
c      ropot(i,j,k+1)= ropot(i,j,k)

	end if   ! Ropot is unstable

      end do   ! k

	end do   ! iter

	end if   ! kb>0

	end do   ! i
	end do   ! j


	return
	end

      subroutine sound(T,S,P,Ro,c_sound)

!     sound speed as function of
!     salinity, potential temperature and pressure, as in
!     Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted JAOT

!     Reprogrammed for speed of sound and single precision T,S,P by Nikolay Iakovlev
!     14.12.2012
!
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C, ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!
!   rho              : in-situ density                    (kg m^-3)
!   rho_p            : partial derivative wrt p           (kg m^-3 dbar^-1)
!
!   check values     : sound(20,20,1000,...) gives rho_p =   4.317589133273301d-3
!

      real*8 
     *th,pth,s8,p8,rho,rho_p,th2,sqrts,aden,anum_p,aden_p,rec_aden

	th= dble(T)
	s8= dble(S)
	p8= dble(P)

      rho= 1.d3*(1.02d0 +DBLE(Ro))

      th2 = th*th
      sqrts = dsqrt(s8)

      aden =       1.0000000000000000d+00 + 	       
     *             th*( 7.2815210113327091d-03 + 	       
     *                  th*(-4.4787265461983921d-05 + 	       
     *                      th*( 3.3851002965802430d-07 + 	       
     *                          th*  1.3651202389758572d-10))) +       
     *             s8*( 1.7632126669040377d-03 - 	       
     *                 th*( 8.8066583251206474d-06 + 	       
     *                      th2*  1.8832689434804897d-10) +	       
     *                 sqrts*( 5.7463776745432097d-06 +
     *                         th2*  1.4716275472242334d-09))

      anum_p =     1.1798263740430364d-02 +
     *             th2*  9.8920219266399120d-08 +
     *             s8*    4.6996642771754730d-06

      aden_p =     6.7103246285651894d-06


      if(p8.ne.0.d0) then

      pth = p8*th

      aden = aden + p8*( 6.7103246285651894d-06 -
     *                   pth*(th2*  2.4461698007024582d-17 +
     *                         p8*  9.1534417604289062d-18))

      anum_p = anum_p -     p8*( 5.1724374150308704d-08 +
     *                          th2*  6.5842828015921320d-12)

      aden_p = aden_p -                                      
     *            pth*(th2*  4.8923396014049170d-17 +
     *                  p8*  2.7460325281286720d-17)
      end if


      rec_aden = 1.0d0/aden

      rho_p = (anum_p-aden_p*rho)*rec_aden


!
!      sound speed is 1.0d2/sqrt(rho_p) ! in m/sec
!


      c_sound= REAL(1.0d4/dsqrt(rho_p)) ! in cm/sec


      return
      end

