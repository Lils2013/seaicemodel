      Subroutine MIXERKC

*     version 07.09.12.

      INCLUDE 'slo2.fi'
      dimension AM(kl),BM(kl),CM(kl),FM(kl),Rksi(kl)  !arrays for sweeping
	dimension q2L(kl)

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
c      Parameter
c     *(A1=0.92,A2=0.74,B1=16.55,B2=10.1,C1=0.08,C2=0.7,C3=0.2)

*     Kantha, L. H., On an improved model for the turbulent pbl, J. Atmos. Sci., 
*     60, 2239–2246, 2003.
      Parameter
     *(A1=0.58,A2=0.62,B1=16.55,B2=11.6,C1=0.038,C2=0.7,C3=0.2)

	Parameter (C4=0.53)
	Parameter (cK=0.4)      ! von Karman constant
	Parameter (cL0=1.e2)    ! macroscale for turbulence, set = 1m (not 10)
*	                        ! other proposals see in Mellor and Yamada, 1982
      Parameter (cLmin=1.  )  ! Minimum lenhth scale - open water approx. 1cm
	Parameter (B123=6.4943) ! B123= B1**(2/3) - for Dirichlet b.c.


      INCLUDE 'tparm.fi'

*     ---------------Constants  -------------
      GhL= 1.0/(A2*(B1+12.*A1+3.*B2*(1.-C3)))
      alpha = A2*(1.-6.*A1/B1)
	beta  = 3.*A2*(6.*A1 +B2*(1.-C3))
	gamma = A1*(1.-3.*C1-6.*A1/B1)
	delta = 9.*A1*(1.+2.*A1+A2-C2)
	epsil = 9.*A1*A2
*     --------------------------------------

      puny=1.e-12

      Rksi(:)=0.


      do j=1,jl
	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then

*     Estimation of the length scale

      q2L(:)= q2turb(i,j,:)*cL(i,j,:)

      do k=1,kb
      cLs= max(cLmin, z(k))         ! distance to the surface
	cLb= max(cLmin, z(kb)-z(k))   ! distance to the bottom
	cLd= cLs*cLb/(cLs+cLb)
	
	cL(k)= cL0*cK*cLd/( cK*cLd + cL0) ! Mellor and Yamada, 1982

*     Limitation to avoid unphysical results

	IF(k.eq.1)then

	kp=2
      n= ABS(nt3(i,j,1))

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)

      VB= g*(Ro1-Ro2)/hz(1)

	else

	km=k-1

	if (k.eq.kb) then

      n= ABS(nt3(i,j,kb))

	ppp= 1.e-5*g*row*(z(k)+z(km))/2.                   
      ro1= sigma_t(t(i,j,k ),s(i,j,k ), ppp)
      ro2= sigma_t(t(i,j,km),s(i,j,km), ppp)

      VB= g*(Ro1-Ro2)/hz(km)

	else

	kp=k+1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

	ppp= 1.e-5*g*row*z(k)                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,km),s(i,j,km), ppp)

      VB= g*(Ro1-Ro2)/(hz(k)+hz(km))

	end if

	end if

	VB= VB/row

*     Limitation to avoid unphysical results
      IF(VB.GT.0.) then
      cL(k)=MIN(cL(k), c4*SQRT(q2turb(i,j,k)/max(puny, VB)))
	end if

*     Limitation to avoid low turbulence
      cL(k)=MAX(cL(k), cLmin)

      end do

*     Estimation of kinetic energy

	do k=1,kb-1

	kp=k+1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,kp))

	hzk =hz(k)

      IF(k.EQ.1) THEN
        hzk1 = 0.
      ELSE
	  hzk1 = hz(k-1)
      END IF

      CA=2.*DT/(CG(n)*HZK1+CG(np)*HZK)

*     Bouyoancy fluxes
      
	IF(k.eq.1)then

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)

      VB= g*(Ro1-Ro2)/hz(1)/row

      uz= (u(i,j,kp)-u(i,j,k))/hzk
      vz= (v(i,j,kp)-v(i,j,k))/hzk


	else

	km=k-1

	ppp= 1.e-5*g*row*z(k)                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,km),s(i,j,km), ppp)

      VB= g*(Ro1-Ro2)/(hzk+hzk1)/row

      uz= (u(i,j,kp)-u(i,j,km)) /(hzk+hzk1)
      vz= (v(i,j,kp)-v(i,j,km)) /(hzk+hzk1)

	end if

*     Vertical shift **2 in velocity points
      VS=uz*uz + vz*vz

      if(cond.eq.1)then

	if(k.eq.1) then

*     ---------- Wind in [m/s] -----------
      wmod= SQRT( wx(i,j)**2 + wy(i,j)**2)
*     ------ Wind Stress Components -------
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      windx= cdrag*wx(i,j)
      windy= cdrag*wy(i,j)

*     Mean Drag coefficient
	drag=CDgwd(i,j)
      Aopen= MIN(Aice(0,i,j), 1.0)
      Fx=(windx*Aopen+
     *                drag*(1.-Aopen)*(uice(i,j)-u(i,j,1)) )
      Fy=(windy*Aopen+
     *                drag*(1.-Aopen)*(vice(i,j)-v(i,j,1)) )

	u_star= sqrt(Fx**2 +Fy**2)/row  !!! Friction velocity

	
      q2turb(i,j,1)= B123* u_star**2  !!! ???? Check it!
	
*     coefficients and right hand side for vertical sweeping.

      else

	Prod= dt*(az (i,j,k-1)+az (i,j,k))*VS
	Buoy= dt*(azt(i,j,k-1)+azt(i,j,k))*VB
      Diss= 2.*dt*sqrt(q2turb(i,j,k))/(B1*cL(k))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,k)+q2turb(i,j,k+1)) 
	turb2 =0.5*(q2turb(i,j,k-1)+q2turb(i,j,k)) 


      AM(k)=-0.5*CA*CG(n) *(cL(k-1)+cL(k  ))*Sq*sqrt(turb2)/hzk1
      CM(k)=-0.5*CA*CG(np)*(cL(k  )+cL(k+1))*Sq*sqrt(turb1) /hzk
      BM(k)= 1.-CM(k)-AM(k) +Pminus
      FM(k)= q2turb(i,j,k)  +Pplus

	end if ! k= 1 
		
	else   ! b.c. at z=0

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
	Buoy= 2.*dt*azt(i,j,1)*VB
      Diss= 2.*dt*sqrt(q2turb(i,j,1))/(B1*cL(1))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,1)+q2turb(i,j,2)) 

      AM(1)=0.
      CM(1)=-0.5*CA*CG(n)*(cL(1)+cL(2))*Sq*sqrt(turb1)/hzk
      BM(1)=1.0 -CM(1)     + Pminus  
      FM(1)= q2turb(i,j,1) + Pplus
     +     + 2.*Aopen*dt*100.*(u_star)**3 /hzk  !! Surface waves

      ELSE ! k>1

	Prod= dt*(az (i,j,k-1)+az (i,j,k))*VS
	Buoy= dt*(azt(i,j,k-1)+azt(i,j,k))*VB
      Diss= 2.*dt*sqrt(q2turb(i,j,k))/(B1*cL(k))

	Pplus = Prod - Buoy
	Pminus= Diss

	turb1 =0.5*(q2turb(i,j,k)+q2turb(i,j,k+1)) 
	turb2 =0.5*(q2turb(i,j,k-1)+q2turb(i,j,k)) 

      AM(k)=-0.5*CA*CG(n) *(cL(k-1)+cL(k  ))*Sq*sqrt(turb2)/hzk1
      CM(k)=-0.5*CA*CG(np)*(cL(k  )+cL(k+1))*Sq*sqrt(turb1) /hzk
      BM(k)= 1.-CM(k)-AM(k) +Pminus
      FM(k)= q2turb(i,j,k)  +Pplus
      END IF

      end if ! b.c. at z=0

      end do ! k

*     Bottom

      k=kb
	km=kb-1
      n=ABS(nt3(i,j,kb))

	ppp= 1.e-5*g*row*(z(k)+z(km))/2.                   
      ro1= sigma_t(t(i,j,k ),s(i,j,k ), ppp)
      ro2= sigma_t(t(i,j,km),s(i,j,km), ppp)

      VB= g*(Ro1-Ro2)/hz(km)/row

      uz= (u(i,j,k)-u(i,j,km))/hz(km)
      vz= (v(i,j,k)-v(i,j,km))/hz(km)


*     Vertical shift **2 in velocity points
      VS=uz*uz + vz*vz

      CA=2.*DT/HZ(Kb-1)

	turb2 =0.5*(q2turb(i,j,kb-1)+q2turb(i,j,kb)) 

	Prod= 2.*dt*az (i,j,kb-1)*VS
	Buoy= 2.*dt*azt(i,j,kb-1)*VB
      Diss= 2.*dt*sqrt(q2turb(i,j,kb))/(B1*cL(kb))

	Pplus = Prod - Buoy
	Pminus= Diss

      AM(kb)=-0.5*CA*(cL(kb)+cL(kb-1))*Sq*sqrt(turb2)/hz(kb-1)
      CM(kb)= 0.
      BM(kb)= 1.-AM(kb)      +Pminus
      FM(kb)= q2turb(i,j,kb) +Pplus

      k1=1
	if(cond.eq.1) k1=2

      call FACTOR(kl,am,bm,cm,fm,rksi,k1,kb)
      do k=k1,kb
      q2turb(i,j,k)= MAX(0.,rksi(k))
      end do


      do k=1,kb-1

	kp=k+1
	n=abs(nt3(i,j,k))

*     projection in the middle of the box in coefficients points
	turb =0.5*(q2turb(i,j,k)+q2turb(i,j,k+1)) 

	ppp= 1.e-5*g*row*(z(k)+z(kp))/2.                   
      ro1= sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2= sigma_t(t(i,j,k ),s(i,j,k ), ppp)

      VB= g*(Ro1-Ro2)/hz(k)/row


*     Following the arguments of Galperin et al. (1988), the stability functions
*     Su and Ss are only a function of Gh, which renders the model less prone to
*     oscillations than its original version, in which the stability functions 
*     also depend on the current shearing (e.g., Deleersnijder and Luyten,1994).

      cLz= 0.5*(cL(k)+cL(k+1))

	Gh= -VB*(cLz**2)/max(puny,turb)

*     Limitation to avoid unphysical results 
      Gh= MIN(Gh,GhL)
	Gh= MAX(-0.28,Gh)

*     stability functions SFu, SFs, Kantha and Clayson (1994)

	SFs = alpha/(1.- beta*Gh)
      SFu = (gamma + delta*SFs*Gh)/(1.- epsil*Gh)

*     New coefficients

      if(turb. GT. 1.e-2) then
      
      Az (i,j,k)= cLz*SFu*SQRT(turb) + AZbg
      AzT(i,j,k)= cLz*SFs*SQRT(turb) + AZTbg
      AzS(i,j,k)= cLz*SFs*SQRT(turb) + AZTbg


	else ! Kantha and Clayson (1994)

*     Shear instability
*     Vertical shift
      uz= u(i,j,k+1)-u(i,j,k)
      vz= v(i,j,k+1)-v(i,j,k)
      Ri= VB*hz(k)**2/MAX(puny,uz**2+vz**2)

	if(Ri .GE. 0.) then ! stable layer

	if(Ri.LE.0.7) then      
      add= 50.*(1.-(Ri/0.7)**2)**3
      Az (i,j,k)= add + AZbg
      AzT(i,j,k)= add + AZTbg
      AzS(i,j,k)= add + AZTbg  
	end if

	else ! unstable layer
      Az (i,j,k)= 50.
      AzT(i,j,k)= 50.
      AzS(i,j,k)= 50.   
	end if

	end if  ! turb > 1e-2

	end do ! k

	end if ! km2 >0
	end do ! i
	end do ! j

	return
	end

	subroutine TS_convection
      INCLUDE 'slo2.fi'
      INCLUDE 'tparm.fi'

*     Check instability
      do j=1,jl
	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then

      nflag= 1

	do while (nflag. eq. 0)

      do k=1,kb-1

      if( Ropot(i,j,k+1) .LT. Ropot(i,j,k)) then

	write(*,*) 'Instability:',i,j,k, Ropot(i,j,k+1), Ropot(i,j,k)

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

      nflag= 1

	else

	nflag= 0

	end if   ! Ropot is unstable

      end do   ! k

	end do   ! nflag=0

	end if   ! kb>0

	end do   ! i
	end do   ! j


	return
	end

