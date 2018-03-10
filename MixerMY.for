      Subroutine MIXER

*     version 06.03.16.

      INCLUDE 'Slo2.fi'
      dimension AM(kl),BM(kl),CM(kl),FM(kl),Rksi(kl) ! for q2
      dimension AML(kl),BML(kl),CML(kl),FML(kl),RksiL(kl) ! for q2l
	dimension cL(kl),VB(kl),AzQ(kl)

*     Vertical mixing coefficients by Mellor and Yamada's 
*     level 2.5 model (Mellor and Yamada, 1982).
*     As described in ROMS manual:
*     Allen, J. S., P. A. Newberger, and J. Federiuk, 1995:
*     Upwelling circulation on the Oregon continental shelf. Part I: 
*     Response to Idealized Forcing, J. Phys. Oceanogr., 25, 1843-1866.

*     References:
*     Mellor G.L. and T. Yamada (1982). Development of a turbulence closure 
*     model for geophysical fluid problems. Rev. Geophys. Spac. Phys. 20(4), 
*     851-875.
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

*     Correction for E1,E3:
*     Burchard, H.: On the q2l Equation by Mellor and Yamada (1982),
*     J. Phys. Oceanogr., 31, 1377–1387, 2001.

*     Also:
*     K. O’Driscoll and V. Kamenkovich. The analysis of large-scale turbulence
*     characteristics in the Indonesian seas derived from a regional model based
*     on the Princeton Ocean Model. Ocean Sci., 8, 615–631, 2012. 
*     doi:10.5194/os-8-615-2012


*     --------------------------------------------------------------------------
*     NOTE: Spatial discretization should be indentical to discretizations for 
*           momentum and scalar approximations for averaged equations to be sure
*           that energy cascades from fluctuations to averages are correct.
*           Momentum production indentical to vertical viscosity, buyoancy 
*           similar to hydrostatic equations approximation.
*     --------------------------------------------------------------------------
     
*     Original
      Parameter
     *     (A1=0.92,A2=0.74,B1=16.6,B2=10.1,C1=0.08,C2=0.7,C3=0.2,
     *      E1=1.8,E2=1.33)
      Parameter (E3=E1)  ! E3=E1 or E3=5.0 Burchard,2001

*     Kantha, L. H., On an improved model for the turbulent pbl, J. Atmos. Sci., 
*     60, 2239–2246, 2003.
c      Parameter
c     *     (A1=0.58,A2=0.62,B1=16.55,B2=11.6,C1=0.038,C2=0.7,C3=0.2,
c     *      E1=1.8,E2=1.33)

	Parameter (C4=0.53)
	Parameter (cK=0.4)      ! von Karman constant
*	                        ! other proposals see in Mellor and Yamada, 1982
      Parameter (cLmin=1.0)   ! Minimum length scale, by default 10 cm
	Parameter (B123=6.4943) ! B123= B1**(2/3) -for ROMS, Allen et.al., 1995


      INCLUDE 'Tparm.fi'

*     ---------------Constants  -------------
      GhL    = 0.028
      alpha  = A2*(1.-6.*A1/B1)
	beta   = 3.*A2*B2 +18.*A1*A2
	gamma  = A1*(1.-3.*C1-6.*A1/B1)
	delta  = 9.*A1*(2.+A2)
	epsil  = 9.*A1*A2
	a_intern_waves = 1.0     ! 0.7 Mellor 1989, 0.8-1.4 Doronin 2000
	AZS_to_AZT     = 1.0
*     --------------------------------------

      puny=1.e-12
      
      Drag2= 1.0e-3  ! From Ibrayev (1.3), Clio 3.0, OPA 8.1 (1.0)
      
      Rksi =0.
      RksiL=0.

      do j=1,jl
	S0=Si(j)
	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then
      
*     Mean square variation of the ice thickness - to compute the length scale
	himean= 0.
	do mg=1,mgrad
	himean=himean+(roi*hice(mg,i,j)+rosdry*hsnow(mg,i,j))/row
	enddo
      if(1.-aice(0,i,j) .LT. Aimin) then
	himean=0.
	end if
      sigmahice= Aice(0,i,j) * himean**2
	do mg=1,mgrad
	sigmahice=sigmahice + 
     *Aice(mg,i,j)*(hice(mg,i,j)/max(Aimin,aice(mg,i,j)) -himean)**2
	enddo

	sigmahice= 2.*SQRT(sigmahice)
      

*     ---------- Wind in [m/s] -----------
      wmod= SQRT( wx(i,j)**2 + wy(i,j)**2)
*     ------ Wind Stress Components -------
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      windx= cdrag*wx(i,j)
      windy= cdrag*wy(i,j)

      Aopen= MIN(Aice(0,i,j), 1.0)

* ------------------------------------------------------
*     Master length scale using previous time step
*     And N

      do k=1,kb
      
      if(k.eq.1 .or. k.eq.kb) then
      cL(k)= 0.
      else
	cL(k)= q2L(i,j,k)/max(puny,q2turb(i,j,k))
	end if
* ------------------------------------------------------
	
	
* Limitation to avoid unphysical results
* if unstable density profile - then local depth or somth else...
* Blanke B. and P. Delecluse (1993). Variability of the tropical 
* Atlantic Ocean simulated by a general circulation model with
* two different mixed-layer physics. J. Phys. Oceanogr. 23, 1363-1388
	
	IF(k.LT.kb)then
	kp=k+1
		
	ppp= 5.e-6*g*row*(z(k)+z(kp)) 
	
      ro1=sigma_t(t(i,j,kp),s(i,j,kp), ppp)
      ro2=sigma_t(t(i,j,k ),s(i,j,k ), ppp)
      VB(k)= 1.e-3*g*REAL(Ro1-Ro2)/hz(k)/row
      
c      ro1=RoPot(i,j,kp)
c      ro2=RoPot(i,j,k )
c      VB(k)= g*(Ro1-Ro2)/hz(k)/row
	
      else
      VB(kb)= VB(kb-1)
	end if	
	
      cL(k)=MIN(cL(k), c4*SQRT(q2turb(i,j,k)/max(puny, VB(k))))

*     Limitation to avoid low turbulence
      cL(k)=MAX(cL(k), cLmin)
      
c      if(k.eq.1) then      
c      cL(1)=MAX(cL(1), sigmahice) ! Sea Ice makes turbulence
c      end if
      
      end do
      
      
      do k=1,kb-1   
      AzQ(k)= cK*Az(i,j,k)  
      end do
      
c      if(i.eq.17.and.j.eq.24) then
c      write(*,*) 'cL',(cL(k), k=1,kl)
c      write(*,*) 'VB',(VB(k), k=1,kl)
c      write(*,*) 'azt',(azt(i,j,k), k=1,kl)
c      write(*,*) 'q2',(q2turb(i,j,k), k=1,kl)
c      write(*,*) 'q2L',(q2L(i,j,k), k=1,kl)
c      endif

*======================================================================

*                             Equations solution

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

*     Vertical shift **2 in velocity points 
      
	IF(k.eq.1)then
	
      AZTVB= azt(i,j,1)*VB(1)
	
      uz2= az(i,j,k)*( (u(i,j,kp)-u(i,j,k ))/hzk )**2 
      vz2= az(i,j,k)*( (v(i,j,kp)-v(i,j,k ))/hzk )**2 
      if (i .eq. 11 .and. j .eq. 1) then
        !write(*,*) "uz2", uz2, vz2, az(i,j,k), u(i,j,kp), hzk
      end if
      VS= uz2 + vz2

	else

	km=k-1
	
	AZTVB= (cg(n )*azt(i,j,km)*hzk1*VB(km) 
     *       +cg(np)*azt(i,j,k )*hzk *VB(k ))
     *       /(cg(n)*hzk1+cg(np)*hzk)
	
      uz2= (cg(np)*az(i,j,k )*(u(i,j,kp)-u(i,j,k ))**2/hzk + 
     *      cg(n )*az(i,j,km)*(u(i,j,k )-u(i,j,km))**2/hzk1  )
     *                                 /(cg(np)*hzk+cg(n)*hzk1)
      vz2= (cg(np)*az(i,j,k )*(v(i,j,kp)-v(i,j,k ))**2/hzk + 
     *      cg(n )*az(i,j,km)*(v(i,j,k )-v(i,j,km))**2/hzk1  )
     *                                 /(cg(np)*hzk+cg(n)*hzk1)
      
      VS= uz2 + vz2
	end if


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
      
*     Mean Skin Drag coefficient under ice
	drag=CDw*(1.-Aopen)*delta_u(i,j)
      tau2= (windx*Aopen+  drag*(uice(i,j)-u(i,j,1)))**2 
     *     +(windy*Aopen+  drag*(vice(i,j)-v(i,j,1)))**2
      
	Prod= dt*VS
	Buoy= dt*AZTVB*a_intern_waves *AZS_to_AZT
      Diss= dt*sqrt(q2turb(i,j,1))/(B1*cL(1))

	Pplus = 2.*(Prod - Buoy)
	Pminus= 2.*Diss

      AM(1)=0.
      CM(1)=-CA*CG(n)*AzQ(1)/hzk
      BM(1)=1.0 -CM(1)     + Pminus
      FM(1)= q2turb(i,j,1) + Pplus  
     +      + 2.*dt*B123*tau2/hz(1)       ! Allen at.al., 1995
     +      + 2.*Aopen*dt*100.*u_star**3/hz(1)  ! Surface waves
     
	Pplus = cL(k)*(E1*Prod - E3*Buoy)
	Wall = 1.0 + E2*(cL(k)*
     * (1./max(cLmin,z(k))+1./max(cLmin,z(kb)-z(k)) ) )**2/(cK**2)
	Pminus= Wall*Diss
         
      AML(1)=0.
      CML(1)=-CA*CG(n)*AZQ(1)/hzk
      BML(1)= 1.-CML(1) +Pminus
      FML(1)= q2L(i,j,1)  +Pplus
     
     
      ELSE ! k>1

	Prod= dt*VS
	Buoy= dt*AZTVB*a_intern_waves *AZS_to_AZT
      Diss= dt*sqrt(q2turb(i,j,k))/(B1*cL(k))

	Pplus = 2.*(Prod - Buoy)
	Pminus= 2.*Diss

      AM(k)=-CA*CG(n )*AZQ(k  )/hzk1
      CM(k)=-CA*CG(np)*AZQ(k+1)/hzk
      BM(k)= 1.-CM(k)-AM(k) +Pminus
      FM(k)= q2turb(i,j,k)  +Pplus
      if (i .eq. 11 .and. j .eq. 1) then
        !write(*,*) "q2turb(i,j,k)", q2turb(i,j,k), Prod, Buoy
      end if
      
	Pplus = cL(k)*(E1*Prod - E3*Buoy)
	Wall = 1.0 + E2*(cL(k)*
     * (1./max(cLmin,z(k))+1./max(cLmin,z(kb)-z(k)) ) )**2/(cK**2)
	Pminus= Wall*Diss
      
      
      AML(k)=-CA*CG(n )*AZQ(k  )/hzk1
      CML(k)=-CA*CG(np)*AZQ(k+1)/hzk
      BML(k)= 1.-CML(k)-AML(k) +Pminus
      FML(k)= q2L(i,j,k)  +Pplus
      
      	
	END IF ! K>1

      end do ! k
            
*     Bottom

      k=kb
	km=kb-1
	hzk1=hz(km)
      n=ABS(nt3(i,j,k))
      
      AZTVB= azt(i,j,km)*VB(km)     

*     Vertical shift **2 in velocity points consistent with momentum equation

      uz2= az(i,j,km)*((u(i,j,k )-u(i,j,km))/hzk1)**2
      vz2= az(i,j,km)*((v(i,j,k )-v(i,j,km))/hzk1)**2
      VS= uz2 + vz2


      CA=2.*DT/HZ(Kb-1)
      
      BDrag = Drag2*sqrt(u(i,j,kb)**2 +v(i,j,kb)**2)
      tauxb = BDrag *u(i,j,kb) 
      tauyb = BDrag *v(i,j,kb) 

	Prod= dt*VS
	Buoy= dt*AZTVB*a_intern_waves *AZS_to_AZT
      Diss= dt*sqrt(q2turb(i,j,kb))/(B1*cL(kb))

	Pplus = 2.*(Prod - Buoy)
	Pminus= 2.*Diss
	
      AM(kb)=-CA*AZQ(kb-1)/hz(kb-1)
      CM(kb)= 0.
      BM(kb)= 1.0 -AM(kb)   +Pminus
      FM(kb)= q2turb(i,j,kb) +Pplus 
     * + 2.*dt*B123*(tauxb**2 +tauyb**2)/hz(kb-1) ! Allen, et.al.,1995
     
	Pplus = cL(kb)*(E1*Prod - E3*Buoy)
	Wall = 1.0 + E2*(cL(kb)*
     * (1./max(cLmin,z(kb))+1./cLmin ) )**2/(cK**2)
	Pminus= Wall*Diss
      
      AML(kb)=-CA*AZQ(kb-1)/hz(kb-1)
      CML(kb)=0.
      BML(kb)= 1.0 -AML(k) +Pminus
      FML(kb)= q2L(i,j,kb)  +Pplus
     
      call FACTOR(KL,am,bm,cm,fm,rksi,1,kb)

      
      do k=1,kb
      if (rksi(k) .ne. rksi(k)) then
        !write(*,*) rksi(k), "rksi(k)"
      end if
      q2turb(i,j,k)= MAX(1.e-4,rksi(k))
      end do
      
      if (i .eq. 11 .and. j .eq. 1) then
        !write(*,*) "matr", fm
      end if

      Rksi=0.
      call FACTOR(KL,amL,bmL,cmL,fmL,rksi,2,kb-1)

      q2L(i,j,1)=0.
      q2L(i,j,kb)=0.
      do k=2,kb-1
      q2L(i,j,k)= MAX(1.e-4,rksi(k))
      end do


*     --------------------------- Final calculations -----------------------------

      do k=1,kb-1

	kp=k+1
	hzk=hz(k)

	n =abs(nt3(i,j,kp))

*     projection in the middle of the box in coefficients points
	turb =0.5*(q2turb(i,j,k)+q2turb(i,j,kp)) 

*     Following the arguments of Galperin et al. (1988), the stability functions
*     Su and Ss are only a function of Gh, which renders the model less prone to
*     oscillations than its original version, in which the stability functions 
*     also depend on the current shearing (e.g., Deleersnijder and Luyten,1994).

*     Galperin B., L.H. Kantha, S. Hassid and A. Rosati (1988). 
*     A quasi-equilibrium turbulent energy model for geophysical flows. 
*     J. Atmos. Sciences 45, 55-62.

*     Deleersnijder E. and P. Luyten (1994). On the practical advantages of the 
*     quasiequilibrium version of the Mellor and Yamada level 2.5 turbulence closure
*     applied to marine modelling. Appl. Math. Modelling 18, 281-287.

      cLz= 0.5*( q2L(i,j,k )/max(puny,q2turb(i,j,k )) +
     *           q2L(i,j,kp)/max(puny,q2turb(i,j,kp))  )
     
      cLz=MIN(cLz, c4*SQRT( turb/max(puny, VB(k))) )

*     Limitation to avoid low turbulence
      cLz=MAX(cLz, cLmin)

c      if(k.eq.1) then      
c      cLz=MAX(cLz, sigmahice)
c      end if
	
      Gh= - VB(k)*(cLz**2)/max(puny,turb)   ! Richardson number
	
*     Limitation to avoid unphysical results 
      Gh= MIN(Gh,GhL)

*     stability functions SFu, SFs, Kantha and Clayson (1994)

	SFs = alpha/(1.- beta*Gh)
      SFu = (gamma + delta*SFs*Gh)/(1.- epsil*Gh)

*     ----------------  New coefficients  ---------------------

*     Here should be the layer beneath the ML

      Az (i,j,k)= cLz*SFu*SQRT(turb) + AZbg(k)  ! for all cases     
      AzT(i,j,k)= cLz*SFs*SQRT(turb) + AZTbg(k)
      AzS(i,j,k)= AZS_to_AZT * cLz*SFs*SQRT(turb) + AZTbg(k)
      if (AzT(i,j,k) .ne. AzT(i,j,k)) then
         !write(*,*) AzT(i,j,k), "AzT(i,j,k)", i, j, k
         !write(*,*) cLz, SFs, SQRT(turb), AZTbg(k), q2turb(i,j,k)
      end if
            

* --------------------------------------------------------------


*     Double diffusion 
*     - requires special treatment of bouyancy frequency!

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
      INCLUDE 'Slo2.fi'
      INCLUDE 'Tparm.fi'

*     version 19.01.2016

*     Convection under open water ot in deep ocean - large plumes, D=100m

	Rhx=R*hx
	Rhy=R*hy
	
*     Check instability
      do j=1,jl
      s0= si(j)
      R2hxhy=s0*Rhx*Rhy/6. 

	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then

      do k=1,kb-1

	kp=k+1

      deltaro= Ropot(i,j,kp)-Ropot(i,j,k)

      if( deltaro .LT. 0. ) then
      
**      np= abs(nt3(i,j,kp))
      np= abs(nt3(i,j,k))
      
c     Neutral bouyancy level
      
      do m= kp,kb-1
      if (Ropot(i,j,m  ).LE.Ropot(i,j,k) .and. 
     &    Ropot(i,j,m+1).GT.Ropot(i,j,k)) kn=m
      end do 
      if (Ropot(i,j,kb).LE.Ropot(i,j,k)) kn=kb
      
c      do m= k,kb-1
c 	   ppp= 1.e-5* g*row*z(k) + 1.e-4*PA(i,j)                
c         ro_new= 1.e-3*REAL(sigma_t(t(i,j,k),s(i,j,k), ppp))
c      if (Ro(i,j,m).LE.Ro_new .and. Ro(i,j,m+1).GT.Ro_new) kn= m+1
c      end do 
      
c      write(*,*) 'kn', k, kn 
     
c     Stokes equation v= 2/9 * r**2 *g * delta_ro/ miu ~ 1 cm/sec      
c     for r= 100m, u~r**2, max(u)=3-10cm/s according LES models
      u_stokes=  -deltaro*1.e5
      
c      write(*,*) 'U stokes ', u_stokes, i, j, k, deltaro
**      vconv = u_stokes *1.6e-3  ! 100 plumes D=200m per 50x50 sq.km
      vconv = u_stokes *1.e-2   ! 1 plume D=100m per 1x1 sq.km
**      vconv = u_stokes            ! 1 giant plume
      vconv = MIN(0.25*hz(k)/dt,vconv)       ! Stability criterium
      Uplume= S0*vconv*cg(np)/6.    ! Volume flux by plumes
   
      QTc(i,j,kn) = Qtc(i,j,kn) +Uplume*(T(i,j,k ) -T(i,j,kn))
	QSc(i,j,kn) = QSc(i,j,kn) +Uplume*(S(i,j,k ) -S(i,j,kn))
 
      do m= kn-1,k,-1 ! upwelling
      
      QTc(i,j,m ) = Qtc(i,j,m)  +Uplume*(T(i,j,m+1) -T(i,j,m ))
      QSc(i,j,m ) = QSc(i,j,m)  +Uplume*(S(i,j,m+1) -S(i,j,m ))
      
      end do      
      
c	write(*,*) 'Instability:',i,j,k, kn, Uplume, u_stokes, vconv
c	write(*,*) 'Profile S:', (S(i,j,m), m=1,kb)
      

	end if   ! Stability criterium deltaro <0

      end do   ! k


	end if   ! kb>0

	end do   ! i
	end do   ! j


	return
	end



      subroutine sound(T,S,P,c_sound)

!     sound speed as function of
!     salinity, potential temperature and pressure, as in
!     Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted JAOT

!     Reprogrammed for speed of sound and single precision T,S,P by Nikolay Iakovlev
!     17.04.2013
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

      real*8 anum, pp, t8,t2,pt,
     *th,pth,s8,p8,rho,rho_p,th2,sqrts,aden,anum_p,aden_p,rec_aden


      pp=DBLE(p) !!-10.1325d0 if atmospheric pressure is taken into account

      t8 = dble(t)
	s8 = dble(s)
      t2 = t8**2
	sqrts = sqrt(s8)

	th= dble(T)
	p8= dble(P)

*     Density in situ

      anum =      9.9984085444849347d+02 +    
     #       t8*( 7.3471625860981584d+00 +   
     #       t8*(-5.3211231792841769d-02 +   
     #       t8*  3.6492439109814549d-04)) + 
     #       s8*( 2.5880571023991390d+00 -   
     #       t8*  6.7168282786692355d-03 +   
     #       s8*  1.9203202055760151d-03) 

      aden =      1.0000000000000000d+00 +    
     #       t8*( 7.2815210113327091d-03 +    
     #       t8*(-4.4787265461983921d-05 +    
     #       t8*( 3.3851002965802430d-07 +    
     #       t8*  1.3651202389758572d-10))) + 
     #       s8*( 1.7632126669040377d-03 -    
     #       t8*( 8.8066583251206474d-06 +    
     #       t2*  1.8832689434804897d-10) +   
     #    sqrts*( 5.7463776745432097d-06 +    
     #       t2*  1.4716275472242334d-09))


      if(p.ne.0.0) then

      pt = pp*t8
                                    
      anum = anum +  pp*( 1.1798263740430364d-02 +   
     #               t2*  9.8920219266399117d-08 +    
     #               s8*  4.6996642771754730d-06 -    
     #               pp*( 2.5862187075154352d-08 +    
     #               t2*  3.2921414007960662d-12))    

      aden = aden +  pp*(     6.7103246285651894d-06 -   
     #               pt*(t2*  2.4461698007024582d-17 +   
     #               pp*      9.1534417604289062d-18))   

      end if


      rho = anum/aden

*     Speed of sound

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

