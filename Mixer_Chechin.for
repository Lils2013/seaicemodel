      Subroutine MIXER_Chechin

*     version 19.04.13.

      INCLUDE 'slo2.fi'
      INCLUDE 'tparm.fi'

      puny=1.e-12

      do j=1,jl
	do i=1,il

	kb=km2(i,j)
      if(kb .GT. 0) then

      cL0= 1.e3                            ! 10m CLIO 3.0 default


      do k=1,kb-1

	kp= k+1

	n= abs(nt3(i,j,k))

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

      uz= (umkp-umk)/hz(k)/cg(n)/3.
      vz= (vmkp-vmk)/hz(k)/cg(n)/3.

c	uZ= (u(i,j,kp)- u(i,j,k))/hz(k)
c	vZ= (v(i,j,kp)- v(i,j,k))/hz(k)

	uz = SQRT(uz**2 + vz**2)

	zk12= 0.5*(z(k)+z(kp))

	cL= 0.4*zk12*cL0/(cL0 +0.4*zk12)


	ppp= 5.e-6*g*row*(z(k)+z(kp))
      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,k)
	T12= 0.5*(T(i,j,k)+T(i,j,kp))
	S12= 0.5*(S(i,j,k)+S(i,j,kp))
	call sound(T12,S12,PPP,csound)
      VB= g*(Ro1-Ro2)/hz(k)/row        -(g/csound)**2

	Ri = VB/max(puny,(UZ**2))

cc	write(*,*) i,j,k, Ri


      IF( Ri . LE. 0.2) then                   

      Az (i,j,k)= (cL**2)*UZ*((1.-5. *Ri)**2)  
      AzT(i,j,k)= Az (i,j,k)

	IF( Ri . LT. 0.) then
      Az (i,j,k)= (cL**2)*UZ*SQRT(1.-16.*Ri)
      AzT(i,j,k)= Az (i,j,k)*(SQRT(SQRT(1.-16.*Ri)))
	end if

	else 

      Az (i,j,k)= 0.
      AzT(i,j,k)= 0.

	end if


cc	write(*,*) i,j,k, Az (i,j,k), AzT(i,j,k)

      Az (i,j,k)=  Az (i,j,k) + AZbg
      AzT(i,j,k)=  AzT(i,j,k) + AZTbg
      AzS(i,j,k)=  AzT(i,j,k) 

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

	ppp= 5.e-6*g*row*(z(k)+z(kp))
      ro1 =Ro(i,j,kp)
      ro2 =Ro(i,j,k)
	T12= 0.5*(T(i,j,k)+T(i,j,kp))
	S12= 0.5*(S(i,j,k)+S(i,j,kp))
	call sound(T12,S12,PPP,csound)
      VB= g*(Ro1-Ro2)/hz(k)/row -(g/csound)**2

      if( VB .LT. 0. ) then

	write(*,*) 'Instability:',iter,VB,i,j,k, Ro1, Ro2

	if(k.eq.1) then
	Tmean = T(i,j,k) +T(i,j,kp)
	Tmean= Tmean /2.
	end if

      if(k.GT.1 .and. k. LT. kb-1) then
	Tmean = T(i,j,k) +T(i,j,kp)
	Tmean= Tmean /2.
	end if

      if(k .eq. kb-1) then
	Tmean = T(i,j,k) + T(i,j,kp)
	Tmean= Tmean /2.
	end if

	T(i,j,k )= Tmean
	T(i,j,kp)= Tmean


	if(k.eq.1) then
	Smean = S(i,j,k) + S(i,j,kp)
	Smean= Smean /2.
	end if

      if(k.GT.1 .and. k. LT. kb-1) then
	Smean = S(i,j,k) + S(i,j,kp)
	Smean= Smean /2.
	end if

      if(k .eq. kb-1) then
	Smean = S(i,j,k) + S(i,j,kp)
	Smean= Smean /2.
	end if

	S(i,j,k )= Smean
	S(i,j,kp)= Smean

	ppp= 1.e-5* g*row*z(k)                 
      ro(i,j,k)= 1.e-3*REAL(sigma_t(t(i,j,k ),s(i,j,k ), ppp))
	ppp= 1.e-5* g*row*z(kp)                 
      ro(i,j,k)= 1.e-3*REAL(sigma_t(t(i,j,kp),s(i,j,kp), ppp))

	end if   ! Stability criterium

      end do   ! k

	end do   ! iter

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

      real*8 anum,
     *th,pth,s8,p8,rho,rho_p,th2,sqrts,aden,anum_p,aden_p,rec_aden


      pp=DBLE(p) !!!!-10.1325 ! if atmospheric pressure is taken into account

      t8 = dble(t)
	s8 = dble(s)
      t2 = t8**2
	sqrts = sqrt(s8)

	th= dble(T)
	s8= dble(S)
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

