      SUBROUTINE iceadvect(Scalar,Scalar1,Scalar2)

*******************************************************
*
*     Scalar 2D trasport and Horizontal diffusion
*     FEM Streamline Upwind Scheme Taylor-Galerkin
*     FCT by Loehner, et.al.
*     R. L¨ohner, K. Morgan, J. Peraire and M. Vahdati.
*     Finite element flux-corrected transport (FEM-FCT)
*     for the Euler and Navier-Stokes equations. 
*     Int. J. Numer. Meth. Fluids 7 (1987) 1093–1109.
*
*     Explicit time stepping - one step Taylor-Galerkin
*     Special treatment of open boundaries.
*     Version 12.08.2012
*
*******************************************************

      parameter (itermax=25, omega=1.0)  ! iteration parameters
      parameter (gamma_fct=1.0e-1)       ! Numerical diffusivity

      INCLUDE 'slo2.fi'

      Dimension Scalar (0:mgrad,0:il1,0:jl1),
     *          Scalar1(0:mgrad,0:il1,0:jl1),
     *          Scalar2(0:mgrad,0:il1,0:jl1) 

*     Auxiliary arrays for FCT calculations

	dimension Rplus(0:il1,0:jl1), Rminus(0:il1,0:jl1), alpha_e(6)
      dimension antidf(6,0:il1,0:jl1) !Antidif fluxes to node (i,j) 

      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      INCLUDE 'tparm.fi'

	Adiff= 0.0e0 ! Background diffusivity to suppress oscillations.
                   ! This oscillations occur primarility at solid bndrs.

*     Cutoffs for ice and snow.
      Aimin= 1.e-6 ! minimal scalar value


      c3 =1./3.
	c9 =c3*c3
	c6= 1./6.
	c12=0.5*c6
	c18= c6/3.
	c24=0.5*c12
      asr= hx/hy
	asr2=asr*asr
      HXX= 1./HX/R
      COEFX=HXX/18.
	HYY= 1./HY/R
      COEFY=HYY/18.
      A= 0.5*HXX*HXX


      RsICE  =0.
      Scalar =0.
      Scalar1=0.
	Rplus  =0.
	Rminus =0.
	antidf =0.


*---------------------------------------------------
*                 Right Hand Side
*---------------------------------------------------     

      DO J=1,JL
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
                     
      DO I=1,IL

	n= nt3(i,j,1)
      IF(n.GT.0) THEN

	do mg=1,mgrad

	call transp_ice(FXa,FYa,i,j,mg,n,uice,vice,Scalar2,KT,
     *                il1,jl1,mgrad)
      FXa =COEFX*FXa
      FYa =COEFY*FYa

      call DIFF_ICE(Scalar2, uice,vice,i,j,mg,n,DIFFa,il1,jl1,kl,
     &                mgrad,KT,Adiff,hx,hy,r,dt)
      Fa=-FXa-FYa+a*DIFFa
      RsICE(MG,i,j)=dt*Fa  

	end do  ! mgrad

      END IF
     
      end do
	end do



*----------------------------------------------------
*                  High Order solution
*                 Mass matrix invertion 
*----------------------------------------------------


      do mg=1,mgrad

      DO iter=1,itermax

	do j=1,jl
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
      r2s12= c6*S0
	
	do i=1,il
	n=nt3(i,j,1)
	IF( n .GT. 0) then

	sum =c24*
     *(KT(1,n)*S1*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i+1,j)+Scalar1(mg,i+1,j-1))
     ++KT(2,n)*S2*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i,j-1)+Scalar1(mg,i+1,j-1))
     ++KT(3,n)*S3*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i,j-1)+Scalar1(mg,i-1,j  ))
     ++KT(4,n)*S4*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i-1,j)+Scalar1(mg,i-1,j+1))
     ++KT(5,n)*S5*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i-1,j+1)+Scalar1(mg,i,j+1))
     ++KT(6,n)*S6*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i,j+1)+Scalar1(mg,i+1,j  )))


      Scalar(mg,i,j)= Scalar1(mg,i,j) 
     *            + omega*(RsICE(mg,i,j)-sum)/(r2s12*CG(n))

c	if(i.eq.9.and.j.eq.41.and.k.eq.kb) 
c     *write(*,*)iter,T(i,j,kb),sum,RsT(i,j,k)



      END IF   ! nt3 >0
	end do
	end do

	Scalar1(:,:,:)= Scalar(:,:,:)    ! We need only increment

	END DO   ! Mass matrix iteration

	end do   ! mg

*----------------------------------------------------------------
*                         Liquid points
*----------------------------------------------------------------

*     Open boundaries, n<0.
	call  icebc
     &(Scalar1,Scalar2,uice,vice,nt3,Si,ice_mask,
     &                                  mgrad,il1,jl1,kl,hx,hy,dt)
 

*----------------------------------------------------------------
*                 End of High Order solution
*----------------------------------------------------------------    


*----------------------------------------------------------------
*                 Low Order Scheme by Loehner
*----------------------------------------------------------------

      do mg=1,mgrad


	do j=1,jl
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
      r2s12= c6*S0
	
	do i=1,il
c	n=abs(nt3(i,j,1))
	n=nt3(i,j,1)
	IF( n .GT. 0) then


C     M*Scalar

	sum =c24*
     *(KT(1,n)*S1*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i+1,j)+Scalar2(mg,i+1,j-1))
     ++KT(2,n)*S2*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i,j-1)+Scalar2(mg,i+1,j-1))
     ++KT(3,n)*S3*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i,j-1)+Scalar2(mg,i-1,j  ))
     ++KT(4,n)*S4*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i-1,j)+Scalar2(mg,i-1,j+1))
     ++KT(5,n)*S5*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i-1,j+1)+Scalar2(mg,i,j+1))
     ++KT(6,n)*S6*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i,j+1)+Scalar2(mg,i+1,j  )) )

	Scalar(mg,i,j)= (1.-gamma_fct)*Scalar2(mg,i,j) 
     *+( RsICE(mg,i,j) +gamma_fct*sum )/(r2s12*CG(n))


      ELSE

	Scalar(mg,i,j)= Scalar2(mg,i,j)+ Scalar1(mg,i,j) 

      END IF   ! nt3 >0
	end do
	end do

	end do   ! mg

*----------------------------------------------------------------
*
*                          Flux Correction
*
*----------------------------------------------------------------


      do mg=1,mgrad


	antidf(:,:,:) = 0.0


	do j=1,jl
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
      r2s12= c6*S0
	
	do i=1,il

c	n=abs(nt3(i,j,1))
	n=nt3(i,j,1)
	IF(n .GT. 0) then

*     Antidiffusive fluxes - 6 in each node.
*     We need low-order solution and increment of the high-order one.

*     Element 1

      Sum1= c6*s0*KT(1,n)*Scalar2(mg,i,j) -
     -c24*KT(1,n)*S1*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i+1,j)+Scalar2(mg,i+1,j-1))

	Sum2= c6*s0*KT(1,n)*Scalar1(mg,i,j) -
     -c24*KT(1,n)*S1*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i+1,j)+Scalar1(mg,i+1,j-1))

      antidf(1,i,j)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*CG(n)) ! Lamped mass matrix inverted

*     Element 2

      Sum1= c6*s0*KT(2,n)*Scalar2(mg,i,j) -
     -c24*KT(2,n)*S2*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i,j-1)+Scalar2(mg,i+1,j-1))

	Sum2= c6*s0*KT(2,n)*Scalar1(mg,i,j) -
     -c24*KT(2,n)*S2*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i,j-1)+Scalar1(mg,i+1,j-1))

      antidf(2,i,j)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*CG(n)) ! Lamped mass matrix inverted

*     Element 3

      Sum1= c6*s0*KT(3,n)*Scalar2(mg,i,j) -
     -c24*KT(3,n)*S3*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i,j-1)+Scalar2(mg,i-1,j))

	Sum2= c6*s0*KT(3,n)*Scalar1(mg,i,j) -
     -c24*KT(3,n)*S3*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i,j-1)+Scalar1(mg,i-1,j))

      antidf(3,i,j)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*CG(n)) ! Lamped mass matrix inverted

*     Element 4

      Sum1= c6*s0*KT(4,n)*Scalar2(mg,i,j) -
     -c24*KT(4,n)*S4*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i-1,j)+Scalar2(mg,i-1,j+1))

	Sum2= c6*s0*KT(4,n)*Scalar1(mg,i,j) -
     -c24*KT(4,n)*S4*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i-1,j)+Scalar1(mg,i-1,j+1))

      antidf(4,i,j)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*CG(n)) ! Lamped mass matrix inverted

*     Element 5

      Sum1= c6*s0*KT(5,n)*Scalar2(mg,i,j) -
     -c24*KT(5,n)*S5*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i,j+1)+Scalar2(mg,i-1,j+1))

	Sum2= c6*s0*KT(5,n)*Scalar1(mg,i,j) -
     -c24*KT(5,n)*S5*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i,j+1)+Scalar1(mg,i-1,j+1))

      antidf(5,i,j)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*CG(n)) ! Lamped mass matrix inverted

*     Element 6

      Sum1= c6*s0*KT(6,n)*Scalar2(mg,i,j) -
     -c24*KT(6,n)*S6*
     *(2.*Scalar2(mg,i,j)+Scalar2(mg,i+1,j)+Scalar2(mg,i,j+1))

	Sum2= c6*s0*KT(6,n)*Scalar1(mg,i,j) -
     -c24*KT(6,n)*S6*
     *(2.*Scalar1(mg,i,j)+Scalar1(mg,i+1,j)+Scalar1(mg,i,j+1))

      antidf(6,i,j)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*CG(n)) ! Lamped mass matrix inverted

      END IF   ! nt3 > 0
	end do
	end do


*     Sum of negative and positive fluxes in each node

	do j=1,jl
	do i=1,il

c	n=abs(nt3(i,j,1))
      n=nt3(i,j,1)
	IF( n .GT. 0) then

	Pplus =0.
	Pminus=0.

      do nelem=1,6
	if(antidf(nelem,i,j) .GE. 0.) then 
	Pplus= Pplus + antidf(nelem,i,j)
	else
	Pminus=Pminus+ antidf(nelem,i,j)
	end if
	end do

*     Admitted max/min on the stencil

      Tmax= -1.e19
	Tmin=  1.e19

      if(KT(1,n) .GT. 0.) then	
	Tmax= 
     &max(Tmax, Scalar(mg,i,j),Scalar(mg,i+1,j),Scalar(mg,i+1,j-1))
	Tmin= 
     &min(Tmin, Scalar(mg,i,j),Scalar(mg,i+1,j),Scalar(mg,i+1,j-1))
*     Additional limitation like Zalesak
	Tmax= 
     &max(Tmax, Scalar2(mg,i,j),Scalar2(mg,i+1,j),Scalar2(mg,i+1,j-1))
	Tmin= 
     &min(Tmin, Scalar2(mg,i,j),Scalar2(mg,i+1,j),Scalar2(mg,i+1,j-1))
	end if

      if(KT(2,n) .GT. 0.) then	
	Tmax= 
     &max(Tmax, Scalar(mg,i,j),Scalar(mg,i,j-1),Scalar(mg,i+1,j-1))
	Tmin= 
     &min(Tmin, Scalar(mg,i,j),Scalar(mg,i,j-1),Scalar(mg,i+1,j-1))
*     Additional limitation like Zalesak
	Tmax= 
     &max(Tmax, Scalar2(mg,i,j),Scalar2(mg,i,j-1),Scalar2(mg,i+1,j-1))
	Tmin= 
     &min(Tmin, Scalar2(mg,i,j),Scalar2(mg,i,j-1),Scalar2(mg,i+1,j-1))
	end if

      if(KT(3,n) .GT. 0.) then	
	Tmax= 
     &max(Tmax, Scalar(mg,i,j),Scalar(mg,i-1,j),Scalar(mg,i,j-1))
	Tmin= 
     &min(Tmin, Scalar(mg,i,j),Scalar(mg,i-1,j),Scalar(mg,i,j-1))
*     Additional limitation like Zalesak
	Tmax= 
     &max(Tmax, Scalar2(mg,i,j),Scalar2(mg,i-1,j),Scalar2(mg,i,j-1))
	Tmin= 
     &min(Tmin, Scalar2(mg,i,j),Scalar2(mg,i-1,j),Scalar2(mg,i,j-1))
	end if

      if(KT(4,n) .GT. 0.) then	
	Tmax= 
     &max(Tmax, Scalar(mg,i,j),Scalar(mg,i-1,j),Scalar(mg,i-1,j+1))
	Tmin= 
     &min(Tmin, Scalar(mg,i,j),Scalar(mg,i-1,j),Scalar(mg,i-1,j+1))
*     Additional limitation like Zalesak
	Tmax= 
     &max(Tmax, Scalar2(mg,i,j),Scalar2(mg,i-1,j),Scalar2(mg,i-1,j+1))
	Tmin= 
     &min(Tmin, Scalar2(mg,i,j),Scalar2(mg,i-1,j),Scalar2(mg,i-1,j+1))
	end if

      if(KT(5,n) .GT. 0.) then	
	Tmax= 
     &max(Tmax, Scalar(mg,i,j),Scalar(mg,i,j+1),Scalar(mg,i-1,j+1))
	Tmin= 
     &min(Tmin, Scalar(mg,i,j),Scalar(mg,i,j+1),Scalar(mg,i-1,j+1))
*     Additional limitation like Zalesak
	Tmax= 
     &max(Tmax, Scalar2(mg,i,j),Scalar2(mg,i,j+1),Scalar2(mg,i-1,j+1))
	Tmin= 
     &min(Tmin, Scalar2(mg,i,j),Scalar2(mg,i,j+1),Scalar2(mg,i-1,j+1))
	end if

      if(KT(6,n) .GT. 0.) then	
	Tmax= 
     &max(Tmax, Scalar(mg,i,j),Scalar(mg,i+1,j),Scalar(mg,i,j+1))
	Tmin= 
     &min(Tmin, Scalar(mg,i,j),Scalar(mg,i+1,j),Scalar(mg,i,j+1))
*     Additional limitation like Zalesak
	Tmax= 
     &max(Tmax, Scalar2(mg,i,j),Scalar2(mg,i+1,j),Scalar2(mg,i,j+1))
	Tmin= 
     &min(Tmin, Scalar2(mg,i,j),Scalar2(mg,i+1,j),Scalar2(mg,i,j+1))
	end if
	

*     Admitted increment in the node

      Qplus = Tmax -Scalar(mg,i,j)
	Qminus= Tmin -Scalar(mg,i,j)

*     Antidiffusive flux input in the node 

      if(Pplus. GT. 1.e-19) then
      Rplus(i,j) = min(1., Qplus/Pplus)
	else
	Rplus(i,j) = 0.
	end if

      if(Pminus. LT. -1.e-19) then
      Rminus(i,j) = min(1., Qminus/Pminus)
	else
	Rminus(i,j) = 0.
	end if


      END IF   ! nt3  > 0
	end do
	end do


*---------------------------------------------------------------     
*                         Solution update
*---------------------------------------------------------------

	do j=1,jl
	do i=1,il
	n=nt3(i,j,1)
	IF( n .GT. 0) then

      alpha_e(:)=1.0

*     Weights alpha

*     Element #1

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(1,i,j) .GT. 0.) then
      a1= min(a1,Rplus (i,j))
	else
      a1= min(a1,Rminus(i,j))
      end if
      if(antidf(3,i+1,j) .GT. 0.) then
      a2= min(a2,Rplus (i+1,j))
	else
      a2= min(a2,Rminus(i+1,j))
      end if
      if(antidf(5,i+1,j-1) .GT. 0.) then
      a3= min(a3,Rplus (i+1,j-1))
	else
      a3= min(a3,Rminus(i+1,j-1))
      end if

	alpha_e(1)= min(a1,a2,a3)


*     Element #2

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(2,i,j) .GT. 0.) then
      a1= min(a1,Rplus (i,j))
	else
      a1= min(a1,Rminus(i,j))
      end if
      if(antidf(4,i+1,j-1) .GT. 0.) then
      a2= min(a2,Rplus (i+1,j-1))
	else
      a2= min(a2,Rminus(i+1,j-1))
      end if
      if(antidf(6,i,j-1) .GT. 0.) then
      a3= min(a3,Rplus (i,j-1))
	else
      a3= min(a3,Rminus(i,j-1))
      end if

	alpha_e(2)= min(a1,a2,a3)


*     Element #3

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(3,i,j) .GT. 0.) then
      a1= min(a1,Rplus (i,j))
	else
      a1= min(a1,Rminus(i,j))
      end if
      if(antidf(5,i,j-1) .GT. 0.) then
      a2= min(a2,Rplus (i,j-1))
	else
      a2= min(a2,Rminus(i,j-1))
      end if
      if(antidf(1,i-1,j) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j))
	else
      a3= min(a3,Rminus(i-1,j))
      end if

	alpha_e(3)= min(a1,a2,a3)


*     Element #4

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(4,i,j) .GT. 0.) then
      a1= min(a1,Rplus (i,j))
	else
      a1= min(a1,Rminus(i,j))
      end if
      if(antidf(2,i-1,j+1) .GT. 0.) then
      a2= min(a2,Rplus (i-1,j+1))
	else
      a2= min(a2,Rminus(i-1,j+1))
      end if
      if(antidf(6,i-1,j) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j))
	else
      a3= min(a3,Rminus(i-1,j))
      end if

	alpha_e(4)= min(a1,a2,a3)


*     Element #5

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(5,i,j) .GT. 0.) then
      a1= min(a1,Rplus (i,j))
	else
      a1= min(a1,Rminus(i,j))
      end if
      if(antidf(3,i,j+1) .GT. 0.) then
      a2= min(a2,Rplus (i,j+1))
	else
      a2= min(a2,Rminus(i,j+1))
      end if
      if(antidf(1,i-1,j+1) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j+1))
	else
      a3= min(a3,Rminus(i-1,j+1))
      end if

	alpha_e(5)= min(a1,a2,a3)


*     Element #6

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(6,i,j) .GT. 0.) then
      a1= min(a1,Rplus (i,j))
	else
      a1= min(a1,Rminus(i,j))
      end if
      if(antidf(2,i,j+1) .GT. 0.) then
      a2= min(a2,Rplus (i,j+1))
	else
      a2= min(a2,Rminus(i,j+1))
      end if
      if(antidf(4,i+1,j) .GT. 0.) then
      a3= min(a3,Rplus (i+1,j))
	else
      a3= min(a3,Rminus(i+1,j))
      end if

	alpha_e(6)= min(a1,a2,a3)


*----------------------------------------------------------------
      do nelem=1,6

      Scalar(mg,i,j)=Scalar(mg,i,j)+alpha_e(nelem)*antidf(nelem,i,j)

	end do
*----------------------------------------------------------------


      END IF   ! nt3 > 0
	end do
	end do
	end do   ! mg

2000  continue


      RETURN
      END

      subroutine transp_ice(FXTK,FYTK,i,j,mg,n,u,v,t,KT,il,jl,mgrad)
      dimension u(0:il,0:jl),v(0:il,0:jl),t(0:mgrad,0:il,0:jl),KT(6,13)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
*     Version 16.07.2012.

      FXTK= 0.
      FYTK= 0.

*---------------     Triangle # 1.  ---------------------------------
      FXTK=FXTK+(U(I,J)+U(I+1,J)+U(I+1,J-1))*
     *          (T(mg,i  ,J  )+T(mg,i+1,J  )+T(mg,i+1,J-1))
     *  *KT(1,n)
*---------------     Triangle # 2.  -----------------------------------
      FYTK=FYTK-S2*(V(I,J)+V(I,J-1)+V(I+1,J-1))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J-1)+T(mg,i+1,J-1))
     *  *KT(2,n)
*---------------     Triangle # 3.  -------------------------------
      FXTK=FXTK-(U(I,J)+U(I-1,J)+U(I,J-1))*
     *          (T(mg,i  ,J  )+T(mg,i-1,J  )+T(mg,i  ,J-1))
     *  *KT(3,n)
      FYTK=FYTK-S3*(V(I,J)+V(I,J-1)+V(I-1,J))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J-1)+T(mg,i-1,J  ))
     *  *KT(3,n)
*---------------     Triangle # 4.  ---------------------------------
      FXTK=FXTK-(U(I,J)+U(I-1,J)+U(I-1,J+1))*
     *          (T(mg,i  ,J  )+T(mg,i-1,J  )+T(mg,i-1,J+1))
     *  *KT(4,n)
*---------------     Triangle # 5.  ----------------------------------
      FYTK=FYTK+S5*(V(I,J)+V(I,J+1)+V(I-1,J+1))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J+1)+T(mg,i-1,J+1))
     *  *KT(5,n)
*---------------     Triangle # 6.  -----------------------------------
      FXTK=FXTK+(U(I,J)+U(I+1,J)+U(I,J+1))*
     *          (T(mg,i  ,J   )+T(mg,i+1,J  )+T(mg,i  ,J+1))
     *  *KT(6,n)
      FYTK=FYTK+S6*(V(I,J)+V(I,J+1)+V(I+1,J))*
     *             (T(mg,i  ,J  )+T(mg,i  ,J+1)+T(mg,i+1,J))
     *  *KT(6,n)

      return
      end

	subroutine icebc
     &(aice,aice2,u,v,nt3,Si,ice_mask,mgrad,il1,jl1,kl,hx,hy,dt)

*     Version 14.08.2012 for FCT scheme
*     Aice is the time increment - NOTE IT!

	dimension aice(0:mgrad,0:il1,0:jl1),aice2(0:mgrad,0:il1,0:jl1),
     &          nt3(0:il1,0:jl1,kl), Si(0:jl1),
     &          ice_mask(0:il1,0:jl1),
     &          u(0:il1,0:jl1),v(0:il1,0:jl1)

	do j=1,jl1-1
	S0=Si(j)
      do i=1,il1-1

	IF(ice_mask(i,j).EQ.1) then
	ni=nt3(i,j,1)

	do mg=1,mgrad
	if(ni.EQ.-6) then

      delta_t=aice(mg,i+1,j)
	delta_x=(aice(mg,i+1,j) + aice2(mg,i+1,j)
     *        -aice(mg,i+2,j)-aice2(mg,i+2,j))/(hx*S0)
	sss= delta_t*(aice2(mg,i+1,j+1)-aice2(mg,i+1,j-1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i+1,j)-aice2(mg,i+1,j-1))/hy
	else
	delta_y=(aice2(mg,i+1,j+1)-aice2(mg,i+1,j))/hy
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)            !-u(i+1,j)/(hx*S0)
	Ry=-delta_t*delta_y*D/hy                 !+v(i+1,j)/hy
	D=1./(1.+Rx)

*	if(i.eq.7.and.j.eq.44)
c	write(*,*)i,j,mg,Rx,Ry,aice2(mg,i,j),aice(mg,i+1,j)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i+1,j)+aice2(mg,i+1,j))
     $                             +Ry*(aice2(mg,i,j-1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i+1,j)+aice2(mg,i+1,j))
     $                             -Ry*(aice2(mg,i,j+1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
	end if  ! n=-6

	if(ni.EQ.-7) then

      delta_t=aice(mg,i,j+1)
	delta_x=(aice(mg,i,j+1)+aice2(mg,i,j+1)
     *        -aice(mg,i,j+2)-aice2(mg,i,j+2))/hy
	sss= delta_t*(aice2(mg,i+1,j+1)-aice2(mg,i-1,j+1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i,j+1)-aice2(mg,i-1,j+1))/(hx*S0)
	else
	delta_y=(aice2(mg,i+1,j+1)-aice2(mg,i,j+1))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy                  !-v(i,j+1)/hy
	Ry=-delta_t*delta_y*D/(hx*S0)             !-u(i,j+1)/(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j+1)+aice2(mg,i,j+1))
     $                             +Ry*(aice2(mg,i-1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j+1)+aice2(mg,i,j+1))
     $                             -Ry*(aice2(mg,i+1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if

	end if  ! n=-7

	if(ni.EQ.-8) then

      delta_t=aice(mg,i-1,j)
	delta_x=(aice(mg,i-1,j)+ aice2(mg,i-1,j)
     *        -aice(mg,i-2,j)- aice2(mg,i-2,j))/(hx*S0)
	sss= delta_t*(aice2(mg,i-1,j+1)-aice2(mg,i-1,j-1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i-1,j)-aice2(mg,i-1,j-1))/hy
	else
	delta_y=(aice2(mg,i-1,j+1)-aice2(mg,i-1,j))/hy
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)             !+u(i-1,j)/(S0*hx)
	Ry=-delta_t*delta_y*D/hy                  !-v(i-1,j)/hy
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i-1,j)+aice2(mg,i-1,j))
     $                             +Ry*(aice2(mg,i,j-1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i-1,j)+aice2(mg,i-1,j))
     $                             -Ry*(aice2(mg,i,j+1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
	end if  ! n=-8

	if(ni.EQ.-9) then
      delta_t=aice(mg,i,j-1)
	delta_x=(aice(mg,i,j-1)+aice2(mg,i,j-1)
     *        -aice(mg,i,j-2)-aice2(mg,i,j-2))/hy
	sss= delta_t*(aice2(mg,i+1,j-1)-aice2(mg,i-1,j-1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i,j-1)-aice2(mg,i-1,j-1))/(hx*S0)
	else
	delta_y=(aice2(mg,i+1,j-1)-aice2(mg,i,j-1))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy                   !+v(i,j-1)/hy
	Ry=-delta_t*delta_y*D/(hx*S0)              !+u(i,j-1)/(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j-1)+aice2(mg,i,j-1))
     &                             +Ry*(aice2(mg,i-1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j-1)+aice2(mg,i,j-1))
     &                             -Ry*(aice2(mg,i+1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
	end if  ! n=-9

	if(ni.EQ.-10) then

	if(nt3(i+1,j,1).GT.0.AND.nt3(i,j-1,1).LT.0)then !transport to W

      delta_t=aice(mg,i+1,j)
	delta_x=(aice(mg,i+1,j)+aice2(mg,i+1,j)
     *        -aice(mg,i+2,j)-aice2(mg,i+2,j))/(hx*S0)
	sss= delta_t*(aice2(mg,i+1,j)-aice2(mg,i+1,j-1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i+1,j)-aice2(mg,i+1,j-1))/hy
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)             !-u(i+1,j)/(hx*S0)
	Ry=-delta_t*delta_y*D/hy                  !+v(i+1,j)/hy    
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i+1,j)+aice2(mg,i+1,j))
     $                             +Ry*(aice2(mg,i,j-1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i+1,j)+aice2(mg,i+1,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if

      end if  ! To W

	if(nt3(i+1,j,1).LT.0.AND.nt3(i,j-1,1).GT.0)then !transport to S
      delta_t=aice(mg,i,j-1)
	delta_x=(aice(mg,i,j-1)+aice2(mg,i,j-1)
     *        -aice(mg,i,j-2)-aice2(mg,i,j-2))/hy
	sss= delta_t*(aice2(mg,i+1,j-1)-aice2(mg,i,j-1))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(aice2(mg,i+1,j-1)-aice2(mg,i,j-1))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy                  !+v(i,j-1)/hy
	Ry=-delta_t*delta_y*D/(hx*S0)             !+u(i,j-1)/(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j-1)+aice2(mg,i,j-1)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j-1)+aice2(mg,i,j-1))
     &                             -Ry*(aice2(mg,i+1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
      end if  ! To S

	end if  ! n=-10

	if(ni.EQ.-11) then

	if(nt3(i-1,j,1).GT.0.AND.nt3(i,j-1,1).LT.0)then !transport to E
      delta_t=aice(mg,i-1,j)
	delta_x=(aice(mg,i-1,j)+aice2(mg,i-1,j)
     *        -aice(mg,i-2,j)-aice2(mg,i-2,j))/(hx*S0)
	sss= delta_t*(aice2(mg,i-1,j)-aice2(mg,i-1,j-1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i-1,j)-aice2(mg,i-1,j-1))/hy
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)             !+u(i-1,j)/(S0*hx)
	Ry=-delta_t*delta_y*D/hy                  !-v(i-1,j)/hy
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i-1,j)+aice2(mg,i-1,j))
     $                             +Ry*(aice2(mg,i,j-1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i-1,j)+aice2(mg,i-1,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
      
	end if ! To E

	if(nt3(i-1,j,1).LT.0.AND.nt3(i,j-1,1).GT.0)then !transport to S
      delta_t=aice(mg,i,j-1)
	delta_x=(aice(mg,i,j-1)+aice2(mg,i,j-1)
     *        -aice(mg,i,j-2)-aice2(mg,i,j-2))/hy
	sss= delta_t*(aice2(mg,i,j-1)-aice2(mg,i-1,j-1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i,j-1)-aice2(mg,i-1,j-1))/(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy                 !+v(i,j-1)/hy
	Ry=-delta_t*delta_y*D/(hx*S0)            !+u(i,j-1)/(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j-1)+aice2(mg,i,j-1))
     &                             +Ry*(aice2(mg,i-1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j-1)+aice2(mg,i,j-1)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
      
	end if  ! To S

	end if  ! n=-11

	if(ni.EQ.-12) then

	if(nt3(i-1,j,1).GT.0.AND.nt3(i,j+1,1).LT.0)then !transport to E
      delta_t=aice(mg,i-1,j)
	delta_x=(aice(mg,i-1,j)+aice2(mg,i-1,j)
     *        -aice(mg,i-2,j)-aice2(mg,i-2,j))/(hx*S0)
	sss= delta_t*(aice2(mg,i-1,j+1)-aice2(mg,i-1,j))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(aice2(mg,i-1,j+1)-aice2(mg,i-1,j))/hy
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)            !+u(i-1,j)/(S0*hx)
	Ry=-delta_t*delta_y*D/hy                 !-v(i-1,j)/hy 
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i-1,j)+aice2(mg,i-1,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i-1,j)+aice2(mg,i-1,j))
     $                             -Ry*(aice2(mg,i,j+1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
      
	end if   ! To E

	if(nt3(i-1,j,1).LT.0.AND.nt3(i,j+1,1).GT.0)then !transport to N
      delta_t=aice(mg,i,j+1)
	delta_x=(aice(mg,i,j+1)+aice2(mg,i,j+1)
     *        -aice(mg,i,j+2)-aice2(mg,i,j+2))/hy
	sss= delta_t*(aice2(mg,i,j+1)-aice2(mg,i-1,j+1))

	if(sss.GT.0.)then
	delta_y=(aice2(mg,i,j+1)-aice2(mg,i-1,j+1))/(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy                  !-v(i,j+1)/hy
	Ry=-delta_t*delta_y*D/(hx*S0)             !-u(i,j+1)/(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j+1)+aice2(mg,i,j+1))
     $                             +Ry*(aice2(mg,i-1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j+1)+aice2(mg,i,j+1)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
      
	end if  ! To N

	end if  ! n=-12

	if(ni.EQ.-13) then

	if(nt3(i+1,j,1).GT.0.AND.nt3(i,j+1,1).LT.0)then !transport to W

      delta_t=aice(mg,i+1,j)
	delta_x=(aice(mg,i+1,j)+aice2(mg,i+1,j)
     *        -aice(mg,i+2,j)-aice2(mg,i+2,j))/(hx*S0)

	sss= delta_t*(aice2(mg,i+1,j+1)-aice2(mg,i+1,j))

	if(sss.GT.0.)then
	delta_y= 0.
	else
	delta_y=(aice2(mg,i+1,j+1)-aice2(mg,i+1,j))/hy
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)             !-u(i+1,j)/(hx*S0)
	Ry=-delta_t*delta_y*D/hy                  !+v(i+1,j)/hy     
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i+1,j)+aice2(mg,i+1,j)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i+1,j)+aice2(mg,i+1,j))
     $                             -Ry*(aice2(mg,i,j+1)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if
	else
      aice(mg,i,j)=0.
	end if

      end if

	if(nt3(i+1,j,1).LT.0.AND.nt3(i,j+1,1).GT.0)then !transport to N
      delta_t=aice(mg,i,j+1)
	delta_x=(aice(mg,i,j+1)+aice2(mg,i,j+1)
     *        -aice(mg,i,j+2)-aice2(mg,i,j+2))/hy
	sss= delta_t*(aice2(mg,i+1,j+1)-aice2(mg,i,j+1))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(aice2(mg,i+1,j+1)-aice2(mg,i,j+1))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy                 !-v(i,j+1)/hy
	Ry=-delta_t*delta_y*D/(hx*S0)            !-u(i,j+1)/(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j+1)+aice2(mg,i,j+1)))
     $             -aice2(mg,i,j)
	else
	aice(mg,i,j)=D*(aice2(mg,i,j)+Rx*(aice(mg,i,j+1)+aice2(mg,i,j+1))
     $                             -Ry*(aice2(mg,i+1,j)-aice2(mg,i,j)))
     $             -aice2(mg,i,j)
	end if

	else
      aice(mg,i,j)=0.
	end if
      
	end if  ! To N

	end if  ! n=-13

      end do  ! mgrad
      end if  ! mask=1
	end do  ! j
	end do  ! i

	return
	end


      SUBROUTINE DIFF_ICE(T,u,v,i,j,mg,n,DIFT,il,jl,kl,mgrad,
     &                      KT,A,hx,hy,R,dt)
*     Version 16.07.2012.
*     Time approximation by Taylor-Glerkin.
*     A is the "physical" or rather "background" diffusivity.
*     Boundary conditions are taken into account in the Main Program.

      dimension u(0:il,0:jl), v(0:il,0:jl)
      dimension T(0:mgrad,0:il,0:jl), KT(6,13)
	dimension a11(6),a12(6),a22(6)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      real KT

      c3=1./3.

      DIFT=0.

      C= 0.5*dt

      u1=(u(i,j)+u(i+1,j)+u(i+1,j-1))*c3
      v1=(v(i,j)+v(i+1,j)+v(i+1,j-1))*c3
      a11(1)= C*u1*u1 +a
      a12(1)= C*u1*v1 
      a22(1)= C*v1*v1 +a


c     -------------------------------

      u2=(u(i,j)+u(i,j-1)+u(i+1,j-1))*c3
      v2=(v(i,j)+v(i,j-1)+v(i+1,j-1))*c3
      a11(2)= C*u2*u2 +a
      a12(2)= C*u2*v2  
      a22(2)= C*v2*v2 +a

c     -------------------------------

      u3=(u(i,j)+u(i-1,j)+u(i,j-1))*c3
      v3=(v(i,j)+v(i-1,j)+v(i,j-1))*c3
      a11(3)= C*u3*u3 +a
      a12(3)= C*u3*v3  
      a22(3)= C*v3*v3 +a

c     -------------------------------

      u4=(u(i,j)+u(i-1,j)+u(i-1,j+1))*c3
      v4=(v(i,j)+v(i-1,j)+v(i-1,j+1))*c3
      a11(4)= C*u4*u4 +a
      a12(4)= C*u4*v4   
      a22(4)= C*v4*v4 +a

c     -------------------------------

      u5=(u(i,j)+u(i-1,j+1)+u(i,j+1))*c3
      v5=(v(i,j)+v(i-1,j+1)+v(i,j+1))*c3
      a11(5)= C*u5*u5 +a
      a12(5)= C*u5*v5  
      a22(5)= C*v5*v5 +a

c     -------------------------------

      u6=(u(i,j)+u(i+1,j)+u(i,j+1))*c3
      v6=(v(i,j)+v(i+1,j)+v(i,j+1))*c3
      a11(6)= C*u6*u6 +a
      a12(6)= C*u6*v6  
      a22(6)= C*v6*v6 +a

c     -------------------------------

c     Triangle # 1.
      DIFT= DIFT+KT(1,n)*(a11(1)*(T(mg,I+1,J)-T(mg,I,J))/S0
     &               +asr*a12(1)*(T(mg,i+1,j)-T(mg,i+1,j-1)))
c     Triangle # 2.
      DIFT= DIFT-KT(2,n)*(a22(2)*S2*(T(mg,I,J)-T(mg,I,J-1))*ASR2
     &               +asr*a12(2)*(T(mg,i+1,j-1)-T(mg,i,j-1)))
c     Triangle # 3.
      DIFT= DIFT+KT(3,n)*(a11(3)*(T(mg,I-1,J)-T(mg,I,J))/S0
     &               +asr*a12(3)*(T(mg,i,j-1)-T(mg,i,j))
     &               -asr*a12(3)*(T(mg,i,j)-T(mg,i-1,j))
     &                -S3*a22(3)*(T(mg,I,J)-T(mg,I,J-1))*ASR2)
c     Triangle # 4.
      DIFT= DIFT+KT(4,n)*(a11(4)*(T(mg,I-1,J)-T(mg,I,J))/S0
     &               +asr*a12(4)*(T(mg,i-1,j)-T(mg,i-1,j+1)))
c     Triangle # 5.
      DIFT= DIFT-KT(5,n)*(a22(5)*S5*(T(mg,I,J)-T(mg,I,J+1))*ASR2
     &               +asr*a12(5)*(T(mg,i-1,j+1)-T(mg,i,j+1)))
c     Triangle # 6.
      DIFT= DIFT+KT(6,n)*(a11(6)*(T(mg,I+1,J)-T(mg,I,J))/S0
     &               +asr*a12(6)*(T(mg,i,j+1)-T(mg,i,j))   
     &               -asr*a12(6)*(T(mg,i,j)-T(mg,i+1,j))
     &                -S6*a22(6)*(T(mg,I,J)-T(mg,I,J+1))*ASR2)

      RETURN
      END
