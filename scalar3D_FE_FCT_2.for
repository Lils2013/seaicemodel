      SUBROUTINE advect3D(Sclr,Sclr1,Sclr2,SclrObs)

*******************************************************
*
*     Sclr 3D trasport and Horizontal diffusion
*     FEM Streamline Upwind Scheme Taylor-Galerkin
*     FCT by Loehner, et.al.
*     R. L�ohner, K. Morgan, J. Peraire and M. Vahdati.
*     Finite element flux-corrected transport (FEM-FCT)
*     for the Euler and Navier-Stokes equations. 
*     Int. J. Numer. Meth. Fluids 7 (1987) 1093�1109.
*
*     Explicit time stepping - one step Taylor-Galerkin
*     Special treatment of open boundaries.
*     Version 27.08.2012
*
*******************************************************

      parameter (itermax=10, omega=1.5)    ! Iteration parameters
      parameter (gamma_fct=1.25e-1)        ! Numerical scalar diffusivity
*     The last parameter may be variable, see 
*     G.E. Georghiou, R. Morrow and A.C. Metaxas, An improved finite-element
*     fluxcorrected transport algorithm. J. Comput. Phys. 148 (1999) 605�620.
*     The general ide ia that gamma_fct=0.5(1-c)c, 
*     c = Courant number on the cluster, gamma_fct<=1/8=0.125


      INCLUDE 'slo2.fi'

	dimension Sclr (0:il1,0:jl1,kl),Sclr1(0:il1,0:jl1,kl),
     *          Sclr2(0:il1,0:jl1,kl),SclrObs(0:il1,0:jl1,kl)

*     Auxiliary arrays for FCT calculations

	dimension Rplus(0:il1,0:jl1,kl), Rminus(0:il1,0:jl1,kl), 
     *          alpha_e(12)
      dimension antidf(12,0:il1,0:jl1,kl) !Antidif fluxes to node (i,j,k) 

      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      common /diffus/ al,alt
	common /consts/ c3,c6,c9,c12,c18,Rhx,R2hx2
      INCLUDE 'tparm.fi'

      Tdamp_In = 2.*dt          ! Time to restore at inflow
      Tdamp_Out=180.*24.*3600.  ! and outflow points.

	c6= 1./6.
	c3= 1./3.
	c9 =c3*c3
	c2=0.5
	c4=0.25
	c12=0.5*c6
	c18=c3*c6
	c24=0.5*c12
      asr=hx/hy
	asr2=asr*asr
	R2=R*R
	Rhx=R*hx
	Rhy=R*hy
	R2hx2=c18*Rhx*Rhy 

      COEFX=c18/Rhx
      COEFY=c18/Rhy


      RsT(:,:,:)    =0.
      Sclr(:,:,:)   =0.
      Sclr1(:,:,:)  =0.
	Rplus(:,:,:)  =0.
	Rminus(:,:,:) =0.
	antidf(:,:,:,:) =0.


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
      A=.5/(S0*Rhx*Rhy)
      r2s12= c12*S0
                     
      DO I=1,IL
      KB = KM2(I,J)
      IF(nt3(i,j,1).GT.0) THEN

*-------------------------------------------------
*         Surface and Nonbottom points
*-------------------------------------------------

      DO K=1, KB -1
      n= nt3(i,j,k)
      np=nt3(i,j,k+1)
	hzk = hz(k)

      IF( k .EQ. 1) THEN
        hzk1= 0.
      ELSE
        hzk1= hz(k-1)
      END IF
	 
      CA= dt

      call DIFFUW_TG(Sclr2,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,il,jl,dt,KT,ALT,hx,hy,r)

      if( k .GT. 1) then 
      call transp(FXTK,FYTK,i,j,k,n,u,v,Sclr2,KT,il1,jl1,kl)
      else
      FXTK =0.
      FYTK =0.
      end if
      call transp(FXTK1,FYTK1,i,j,k,np,u,v,Sclr2,KT,il1,jl1,kl)

      FXTK =COEFX*FXTK
      FYTK =COEFY*FYTK
      FXTK1=COEFX*FXTK1
      FYTK1=COEFY*FYTK1
      DIFT =A*DIFT
      DIFTP=A*DIFTP
      FXT=.5*(FXTK*HZK1 +FXTK1*HZK)
      FYT=.5*(FYTK*HZK1 +FYTK1*HZK)

      IF( k .EQ. 1) THEN
*     -------------- Surface -----------------
      FZT=-r2s12*cg(n)*((Sclr2(I,J,2)+Sclr2(I,J,1))*W(I,J,2)
     *                 -2.*Sclr2(I,J,1)*(W(i,j,1)-PME(i,j)/dt) 
     *                                         )
      ELSE
*     --------------  Deep water -------------
      FZT=-r2s12*(cg(np)*(Sclr2(I,J,K+1)+Sclr2(I,J,K))*W(I,J,K+1)
     -           -cg( n)*(Sclr2(I,J,K  )+Sclr2(I,J,K-1))*W(I,J,K) 
     * )
      END IF

      FT=-FXT-FYT+FZT  +DIFT+DIFTP

      RsT(i,j,k)= CA*FT

*     First guess for mas matrix inversion
	Sclr1(i,j,k)= RsT(i,j,k)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))


c	if(i.eq.19.and.j.eq.5.and.k.eq.1) write(*,*)'init',Sclr1(i,j,k)

      end do

*     --------------- Bottom ------------------
      N = NP
      CA= dt   
      call DIFFUW_TG(Sclr2,u,v,w,hz,i,j,kb,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,il,jl,dt,KT,ALT,hx,hy,r)
      DIFT=A*DIFT
      call transp(FXTK ,FYTK ,i,j,kb  ,n,u,v,Sclr2,KT,il1,jl1,kl)

      FXTK=COEFX*FXTK
      FYTK=COEFY*FYTK
      FXT=.5*FXTK*HZ(KB-1)
      FYT=.5*FYTK*HZ(KB-1)

      FZT=r2s12*cg(n)*(Sclr2(I,J,KB-1)+Sclr2(I,J,KB))*W(I,J,KB)

      FT=-FXT-FYT+FZT  +DIFT

      RsT(I,J,KB)= CA*FT

*     First guess for mas matrix inversion
	Sclr1(i,j,kb)= RsT(i,j,kb)/( r2s12*CG(n)*HZ(KB-1))

cc	if(i.eq.9.and.j.eq.41.and.k.eq.kb)write(*,*)'init:',Sclr1(i,j,kb)
      END IF
     
      end do
	end do

*----------------------------------------------------------------
*                         Liquid points
*----------------------------------------------------------------

      call TSbc_inc(Tdamp_In,Tdamp_Out,
     *dt,Sclr1,Sclr2,Sclrobs,u,v,w,nt3,Si,il1,jl1,kl,hx,hy,hz,R)

*----------------------------------------------------------------
*                      Mass matrix invertion 
*----------------------------------------------------------------

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
      r2s12= c12*S0
	
	do i=1,il

	kb=km2(i,j)

	IF( nt3(i,j,1) .GT. 0) then

      do k=1,kb
	n=nt3(i,j,k)

	sumk=c24*
     *(KT(1,n)*S1*(2.*Sclr1(i,j,k)+Sclr1(i+1,j,k)+Sclr1(i+1,j-1,k))
     ++KT(2,n)*S2*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i+1,j-1,k))
     ++KT(3,n)*S3*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i-1,j  ,k))
     ++KT(4,n)*S4*(2.*Sclr1(i,j,k)+Sclr1(i-1,j,k)+Sclr1(i-1,j+1,k))
     ++KT(5,n)*S5*(2.*Sclr1(i,j,k)+Sclr1(i-1,j+1,k)+Sclr1(i,j+1,k))
     ++KT(6,n)*S6*(2.*Sclr1(i,j,k)+Sclr1(i,j+1,k)+Sclr1(i+1,j  ,k)))

	if(k.LT.kb)then
	kp=k+1
	np=nt3(i,j,kp)
	hzk=hz(k)

	sumkp=c24*
     *(KT(1,np)*S1*(2.*Sclr1(i,j,k)+Sclr1(i+1,j,k)+Sclr1(i+1,j-1,k))
     ++KT(2,np)*S2*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i+1,j-1,k))
     ++KT(3,np)*S3*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i-1,j  ,k))
     ++KT(4,np)*S4*(2.*Sclr1(i,j,k)+Sclr1(i-1,j,k)+Sclr1(i-1,j+1,k))
     ++KT(5,np)*S5*(2.*Sclr1(i,j,k)+Sclr1(i-1,j+1,k)+Sclr1(i,j+1,k))
     ++KT(6,np)*S6*(2.*Sclr1(i,j,k)+Sclr1(i,j+1,k)+Sclr1(i+1,j  ,k)))
	sumkpp=c24*
     *(KT(1,np)*S1*(2.*Sclr1(i,j,kp)+Sclr1(i+1,j,kp)+Sclr1(i+1,j-1,kp))
     ++KT(2,np)*S2*(2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i+1,j-1,kp))
     ++KT(3,np)*S3*(2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i-1,j  ,kp))
     ++KT(4,np)*S4*(2.*Sclr1(i,j,kp)+Sclr1(i-1,j,kp)+Sclr1(i-1,j+1,kp))
     ++KT(5,np)*S5*(2.*Sclr1(i,j,kp)+Sclr1(i-1,j+1,kp)+Sclr1(i,j+1,kp))
     ++KT(6,np)*
     *         S6*(2.*Sclr1(i,j,kp)+Sclr1(i,j+1,kp)+Sclr1(i+1,j  ,kp)))

	else
	hzk= 0.
	np=n
	sumkp  =0.
	sumkpp =0.
	end if

	if(k.GT.1)then
	km=k-1
	hzk1=hz(km)
	sumkm=c24*
     *(KT(1,n)*S1*(2.*Sclr1(i,j,km)+Sclr1(i+1,j,km)+Sclr1(i+1,j-1,km))
     ++KT(2,n)*S2*(2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i+1,j-1,km))
     ++KT(3,n)*S3*(2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i-1,j  ,km))
     ++KT(4,n)*S4*(2.*Sclr1(i,j,km)+Sclr1(i-1,j,km)+Sclr1(i-1,j+1,km))
     ++KT(5,n)*S5*(2.*Sclr1(i,j,km)+Sclr1(i-1,j+1,km)+Sclr1(i,j+1,km))
     ++KT(6,n)*S6*(2.*Sclr1(i,j,km)+Sclr1(i,j+1,km)+Sclr1(i+1,j  ,km)))
 
	else
	hzk1=0.
	sumkm=0.
	end if

cc	sum= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp))
	sum= 0.5*(HZK1*sumk +HZK*sumkp)

      Sclr(i,j,k)= Sclr1(i,j,k) 
     *        + omega*(RsT(i,j,k)-sum)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

c	if(i.eq.12.and.j.eq.9.and.k.eq.23) 
c     *write(*,*)iter,Sclr(i,j,k),
c     *(RsT(i,j,k)-sum)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

	end do   ! k
  
      END IF   ! nt3 >0
	end do
	end do

	do k=1,kl
      do j=1,jl
	do i=1,il
	if(nt3(i,j,k) .GT.0) then
	Sclr1(i,j,k)= Sclr(i,j,k)    ! We need only increment
	end if
	end do
	end do
	end do

	END DO   ! Mass matrix iteration

	do k=1,kl
      do j=1,jl
	do i=1,il
	if(nt3(i,j,k) .LT. 0) then
	Sclr(i,j,k)= Sclr1(i,j,k) +Sclr2(i,j,k)    ! Liquid boundaries
	end if
	end do
	end do
	end do


*----------------------------------------------------------------
*                 End of High Order solution
*----------------------------------------------------------------   
 
*----------------------------------------------------------------
*                 Low Order Scheme by Loehner
*----------------------------------------------------------------

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
      r2s12= c12*S0
	
	do i=1,il
	kb=km2(i,j)
	IF( nt3(i,j,1) .GT. 0) then

      do k=1,kb
	n=nt3(i,j,k)

	sumk=c24*
     * (KT(1,n)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
     + +KT(6,n)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k)))

	if(k.LT.kb)then
	kp=k+1
	np=nt3(i,j,kp)
	hzk=hz(k)

	sumkp=c24*
     * (KT(1,np)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
     + +KT(2,np)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
     + +KT(3,np)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
     + +KT(4,np)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
     + +KT(5,np)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
     + +KT(6,np)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k)))
	sumkpp=c24*
     *(KT(1,np)*S1*(2.*Sclr2(i,j,kp)+Sclr2(i+1,j,kp)+Sclr2(i+1,j-1,kp))
     ++KT(2,np)*S2*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i+1,j-1,kp))
     ++KT(3,np)*S3*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i-1,j  ,kp))
     ++KT(4,np)*S4*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j,kp)+Sclr2(i-1,j+1,kp))
     ++KT(5,np)*S5*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j+1,kp)+Sclr2(i,j+1,kp))
     ++KT(6,np)*
     *         S6*(2.*Sclr2(i,j,kp)+Sclr2(i,j+1,kp)+Sclr2(i+1,j  ,kp)))

	else
	hzk= 0.
	np=n
	sumkp  =0.
	sumkpp =0.
	end if

	if(k.GT.1)then
	km=k-1
	hzk1=hz(km)
	sumkm=c24*
     *(KT(1,n)*S1*(2.*Sclr2(i,j,km)+Sclr2(i+1,j,km)+Sclr2(i+1,j-1,km))
     ++KT(2,n)*S2*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i+1,j-1,km))
     ++KT(3,n)*S3*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i-1,j  ,km))
     ++KT(4,n)*S4*(2.*Sclr2(i,j,km)+Sclr2(i-1,j,km)+Sclr2(i-1,j+1,km))
     ++KT(5,n)*S5*(2.*Sclr2(i,j,km)+Sclr2(i-1,j+1,km)+Sclr2(i,j+1,km))
     ++KT(6,n)*S6*(2.*Sclr2(i,j,km)+Sclr2(i,j+1,km)+Sclr2(i+1,j  ,km)))
 
	else
	hzk1=0.
	sumkm=0.
	end if

c	sum= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp)) ! M*T
	sum= 0.5*(HZK1*sumk +HZK*sumkp) ! M*T lumped in vertical

	Sclr(i,j,k)= (1.-gamma_fct)*Sclr2(i,j,k) 
     *+(RsT(i,j,k) +gamma_fct*sum )/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

c	if(i.eq.3.and.j.eq.11.and.k.eq.19) write(*,*) Sclr(i,j,k), 
c     *Sclr2(i,j,k), RsT(i,j,k), 
c     *sum/( r2s12*(CG(n)*HZK1+CG(np)*HZK)), 
c     *cg(n),cg(np),hzk1,hzk,n,np,kb,sumk,sumkm,sumkp,sumkpp

	end do   ! k

      END IF   ! nt3 >0
	end do   ! i
	end do   ! j

c      call TSbc(Tdamp_In,Tdamp_Out,
c     *  dt,Sclr,Sclr2,SclrObs,u,v,w,nt3,Si,il1,jl1,kl,hx,hy,hz,R)

***      GOTO 2000

*----------------------------------------------------------------
*
*                          Flux Correction
*
*----------------------------------------------------------------

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
      r2s12= c12*S0
	
	do i=1,il
c	IF( nt3(i,j,1) .GT. 0) then
	IF( nt3(i,j,1) .NE. 0) then

	kb=km2(i,j)

      do k=1,kb
	n=abs(nt3(i,j,k))

	hzk =hz(k)

*     Antidiffusive fluxes - 12 in each node.
*     We need low-order solution and increment of the high-order one.


      if(k.GT.1) then

	km=k-1
	hzk1=hz(km)

*     Element 1

	sumk=c24*
     *KT(1,n)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
	sumkm=c24*
     *KT(1,n)*S1*(2.*Sclr2(i,j,km)+Sclr2(i+1,j,km)+Sclr2(i+1,j-1,km))

cc      Sum1= r2s12*KT(1,n)*Sclr2(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum1= r2s12*KT(1,n)*Sclr2(i,j,k)*HZK1 - 0.5*HZK1*sumk 

	sumk=c24*
     *KT(1,n)*S1*(2.*Sclr1(i,j,k)+Sclr1(i+1,j,k)+Sclr1(i+1,j-1,k))
	sumkm=c24*
     *KT(1,n)*S1*(2.*Sclr1(i,j,km)+Sclr1(i+1,j,km)+Sclr1(i+1,j-1,km))

cc	Sum2= r2s12*KT(1,n)*Sclr1(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm)
	Sum2= r2s12*KT(1,n)*Sclr1(i,j,k)*HZK1 - 0.5*HZK1*sumk

      antidf(1,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted


*     Element 2

	sumk=c24*
     *KT(2,n)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
	sumkm=c24*
     *KT(2,n)*S2*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i+1,j-1,km))

cc      Sum1= r2s12*KT(2,n)*Sclr2(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum1= r2s12*KT(2,n)*Sclr2(i,j,k)*HZK1 - 0.5*HZK1*sumk 

	sumk=c24*
     *KT(2,n)*S2*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i+1,j-1,k))
	sumkm=c24*
     *KT(2,n)*S2*(2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i+1,j-1,km))

cc	Sum2= r2s12*KT(2,n)*Sclr1(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm)
	Sum2= r2s12*KT(2,n)*Sclr1(i,j,k)*HZK1 - 0.5*HZK1*sumk

      antidf(2,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted


*     Element 3

	sumk=c24*
     *KT(3,n)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
	sumkm=c24*
     *KT(3,n)*S3*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i-1,j  ,km))

cc      Sum1= r2s12*KT(3,n)*Sclr2(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum1= r2s12*KT(3,n)*Sclr2(i,j,k)*HZK1 - 0.5*HZK1*sumk 

	sumk=c24*
     *KT(3,n)*S3*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i-1,j  ,k))
	sumkm=c24*
     *KT(3,n)*S3*(2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i-1,j  ,km))

cc      Sum2= r2s12*KT(3,n)*Sclr1(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum2= r2s12*KT(3,n)*Sclr1(i,j,k)*HZK1 - 0.5*HZK1*sumk 

      antidf(3,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

*     Element 4

	sumk=c24*
     *KT(4,n)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
	sumkm=c24*
     *KT(4,n)*S4*(2.*Sclr2(i,j,km)+Sclr2(i-1,j,km)+Sclr2(i-1,j+1,km))

cc      Sum1= r2s12*KT(4,n)*Sclr2(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum1= r2s12*KT(4,n)*Sclr2(i,j,k)*HZK1 - 0.5*HZK1*sumk 

	sumk=c24*
     *KT(4,n)*S4*(2.*Sclr1(i,j,k)+Sclr1(i-1,j,k)+Sclr1(i-1,j+1,k))
	sumkm=c24*
     *KT(4,n)*S4*(2.*Sclr1(i,j,km)+Sclr1(i-1,j,km)+Sclr1(i-1,j+1,km))

cc      Sum2= r2s12*KT(4,n)*Sclr1(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum2= r2s12*KT(4,n)*Sclr1(i,j,k)*HZK1 - 0.5*HZK1*sumk 

      antidf(4,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

*     Element 5

	sumk=c24*
     *KT(5,n)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
	sumkm=c24*
     *KT(5,n)*S5*(2.*Sclr2(i,j,km)+Sclr2(i-1,j+1,km)+Sclr2(i,j+1,km))

cc      Sum1= r2s12*KT(5,n)*Sclr2(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum1= r2s12*KT(5,n)*Sclr2(i,j,k)*HZK1 - 0.5*HZK1*sumk 

	sumk=c24*
     *KT(5,n)*S5*(2.*Sclr1(i,j,k)+Sclr1(i-1,j+1,k)+Sclr1(i,j+1,k))
	sumkm=c24*
     *KT(5,n)*S5*(2.*Sclr1(i,j,km)+Sclr1(i-1,j+1,km)+Sclr1(i,j+1,km))

cc      Sum2= r2s12*KT(5,n)*Sclr1(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum2= r2s12*KT(5,n)*Sclr1(i,j,k)*HZK1 - 0.5*HZK1*sumk 

      antidf(5,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

*     Element 6

	sumk=c24*
     *KT(6,n)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k))
	sumkm=c24*
     *KT(6,n)*S6*(2.*Sclr2(i,j,km)+Sclr2(i,j+1,km)+Sclr2(i+1,j  ,km))

cc      Sum1= r2s12*KT(6,n)*Sclr2(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum1= r2s12*KT(6,n)*Sclr2(i,j,k)*HZK1 - 0.5*HZK1*sumk 

	sumk=c24*
     *KT(6,n)*S6*(2.*Sclr1(i,j,k)+Sclr1(i,j+1,k)+Sclr1(i+1,j  ,k))
	sumkm=c24*
     *KT(6,n)*S6*(2.*Sclr1(i,j,km)+Sclr1(i,j+1,km)+Sclr1(i+1,j  ,km))

cc      Sum2= r2s12*KT(6,n)*Sclr1(i,j,k)*HZK1 - c6*HZK1*(2.*sumk+sumkm) 
      Sum2= r2s12*KT(6,n)*Sclr1(i,j,k)*HZK1 - 0.5*HZK1*sumk 

      antidf(6,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

      end if

	if(k.LT.kb) then
*     Lower part

	kp=k+1
	np=abs(nt3(i,j,kp))

*     Element 7

	sumk=c24*
     *KT(1,np)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
	sumkp=c24*
     *KT(1,np)*S1*(2.*Sclr2(i,j,kp)+Sclr2(i+1,j,kp)+Sclr2(i+1,j-1,kp))

cc      Sum1= r2s12*KT(1,np)*Sclr2(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum1= r2s12*KT(1,np)*Sclr2(i,j,k)*HZK - 0.5*HZK*sumk 

	sumk=c24*
     *KT(1,np)*S1*(2.*Sclr1(i,j,k)+Sclr1(i+1,j,k)+Sclr1(i+1,j-1,k))
	sumkp=c24*
     *KT(1,np)*S1*(2.*Sclr1(i,j,kp)+Sclr1(i+1,j,kp)+Sclr1(i+1,j-1,kp))

cc	Sum2= r2s12*KT(1,np)*Sclr1(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp)
	Sum2= r2s12*KT(1,np)*Sclr1(i,j,k)*HZK - 0.5*HZK*sumk

      antidf(7,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted


*     Element 8

	sumk=c24*
     *KT(2,np)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
	sumkp=c24*
     *KT(2,np)*S2*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i+1,j-1,kp))

cc      Sum1= r2s12*KT(2,np)*Sclr2(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum1= r2s12*KT(2,np)*Sclr2(i,j,k)*HZK - 0.5*HZK*sumk 

	sumk=c24*
     *KT(2,np)*S2*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i+1,j-1,k))
	sumkp=c24*
     *KT(2,np)*S2*(2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i+1,j-1,kp))

cc	Sum2= r2s12*KT(2,np)*Sclr1(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp)
	Sum2= r2s12*KT(2,np)*Sclr1(i,j,k)*HZK - 0.5*HZK*sumk

      antidf(8,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted


*     Element 9

	sumk=c24*
     *KT(3,np)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
	sumkp=c24*
     *KT(3,np)*S3*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i-1,j  ,kp))

cc      Sum1= r2s12*KT(3,np)*Sclr2(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum1= r2s12*KT(3,np)*Sclr2(i,j,k)*HZK - 0.5*HZK*sumk 

	sumk=c24*
     *KT(3,np)*S3*(2.*Sclr1(i,j,k)+Sclr1(i,j-1,k)+Sclr1(i-1,j  ,k))
	sumkp=c24*
     *KT(3,np)*S3*(2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i-1,j  ,kp))

cc      Sum2= r2s12*KT(3,np)*Sclr1(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum2= r2s12*KT(3,np)*Sclr1(i,j,k)*HZK - 0.5*HZK*sumk 

      antidf(9,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

*     Element 10

	sumk=c24*
     *KT(4,np)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
	sumkp=c24*
     *KT(4,np)*S4*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j,kp)+Sclr2(i-1,j+1,kp))

cc      Sum1= r2s12*KT(4,np)*Sclr2(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum1= r2s12*KT(4,np)*Sclr2(i,j,k)*HZK - 0.5*HZK*sumk 

	sumk=c24*
     *KT(4,np)*S4*(2.*Sclr1(i,j,k)+Sclr1(i-1,j,k)+Sclr1(i-1,j+1,k))
	sumkp=c24*
     *KT(4,np)*S4*(2.*Sclr1(i,j,kp)+Sclr1(i-1,j,kp)+Sclr1(i-1,j+1,kp))

cc      Sum2= r2s12*KT(4,np)*Sclr1(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum2= r2s12*KT(4,np)*Sclr1(i,j,k)*HZK - 0.5*HZK*sumk 

      antidf(10,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

*     Element 11

	sumk=c24*
     *KT(5,np)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
	sumkp=c24*
     *KT(5,np)*S5*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j+1,kp)+Sclr2(i,j+1,kp))

cc      Sum1= r2s12*KT(5,np)*Sclr2(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum1= r2s12*KT(5,np)*Sclr2(i,j,k)*HZK - 0.5*HZK*sumk 

	sumk=c24*
     *KT(5,np)*S5*(2.*Sclr1(i,j,k)+Sclr1(i-1,j+1,k)+Sclr1(i,j+1,k))
	sumkp=c24*
     *KT(5,np)*S5*(2.*Sclr1(i,j,kp)+Sclr1(i-1,j+1,kp)+Sclr1(i,j+1,kp))

cc      Sum2= r2s12*KT(5,np)*Sclr1(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum2= r2s12*KT(5,np)*Sclr1(i,j,k)*HZK - 0.5*HZK*sumk 

      antidf(11,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

*     Element 12

	sumk=c24*
     *KT(6,np)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k))
	sumkp=c24*
     *KT(6,np)*S6*(2.*Sclr2(i,j,kp)+Sclr2(i,j+1,kp)+Sclr2(i+1,j  ,kp))

cc      Sum1= r2s12*KT(6,np)*Sclr2(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum1= r2s12*KT(6,np)*Sclr2(i,j,k)*HZK - 0.5*HZK*sumk 

	sumk=c24*
     *KT(6,np)*S6*(2.*Sclr1(i,j,k)+Sclr1(i,j+1,k)+Sclr1(i+1,j  ,k))
	sumkp=c24*
     *KT(6,np)*S6*(2.*Sclr1(i,j,kp)+Sclr1(i,j+1,kp)+Sclr1(i+1,j  ,kp))

cc      Sum2= r2s12*KT(6,np)*Sclr1(i,j,k)*HZK - c6*HZK*(2.*sumk+sumkp) 
      Sum2= r2s12*KT(6,np)*Sclr1(i,j,k)*HZK - 0.5*HZK*sumk 

      antidf(12,i,j,k)= (gamma_fct*Sum1+Sum2) 
     *         / (r2s12*(CG(n)*HZK1+CG(np)*HZK)) ! Lamped mass matrix inverted

	end if


	end do   ! k

      END IF   ! nt3 ne 0
	end do
	end do


*     Sum of negative and positive fluxes in each node

	do j=1,jl
	do i=1,il
c	IF( nt3(i,j,1) .GT. 0) then
	IF( nt3(i,j,1) .NE. 0) then

	kb=km2(i,j)

      do k=1,kb

      n=abs(nt3(i,j,k))

	Pplus =0.
	Pminus=0.

      do nelem=1,12
	if(antidf(nelem,i,j,k) .GE. 0.) then 
	Pplus= Pplus + antidf(nelem,i,j,k)
	else
	Pminus=Pminus+ antidf(nelem,i,j,k)
	end if
	end do

*     Admitted max/min on the stencil

      Tmax= -1.e19
	Tmin=  1.e19

      if(KT(1,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i+1,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i+1,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i+1,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i+1,j-1,k ))
	end if

      if(KT(2,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i,j-1,k ),Sclr(i+1,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i,j-1,k ),Sclr(i+1,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i,j-1,k ),Sclr2(i+1,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i,j-1,k ),Sclr2(i+1,j-1,k ))
	end if

      if(KT(3,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i,j-1,k ))
	end if

      if(KT(4,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i-1,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i-1,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i-1,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i-1,j+1,k ))
	end if

      if(KT(5,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i,j+1,k ),Sclr(i-1,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i,j+1,k ),Sclr(i-1,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i,j+1,k ),Sclr2(i-1,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i,j+1,k ),Sclr2(i-1,j+1,k ))
	end if

      if(KT(6,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i,j+1,k ))
	end if
	

	      
*     Admitted increment in the node

      Qplus = Tmax -Sclr(i,j,k)
	Qminus= Tmin -Sclr(i,j,k)

*     Antidiffusive flux input in the node 

      if(Pplus. GT. 1.e-19) then
      Rplus(i,j,k) = min(1., Qplus/Pplus)
	else
	Rplus(i,j,k) = 0.
	end if

      if(Pminus. LT. -1.e-19) then
      Rminus(i,j,k) = min(1., Qminus/Pminus)
	else
	Rminus(i,j,k) = 0.
	end if

	end do   ! k

      END IF   ! nt3 > 0
	end do
	end do

*---------------------------------------------------------------     
*                         Solution update
*---------------------------------------------------------------

	do j=1,jl
	do i=1,il
	IF( nt3(i,j,1) .GT. 0) then
c	IF( nt3(i,j,1) .NE. 0) then

	kb=km2(i,j)

      do k=1,kb
	n=abs(nt3(i,j,k))

      alpha_e(:)=1.0

*     Weights alpha


*     Element #1

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(1,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(3,i+1,j,k) .GT. 0.) then
      a2= Rplus (i+1,j,k)
	else
      a2= Rminus(i+1,j,k)
      end if
      if(antidf(5,i+1,j-1,k) .GT. 0.) then
      a3= Rplus (i+1,j-1,k)
	else
      a3= Rminus(i+1,j-1,k)
      end if

      if(k.GT.1) then

      km=k-1

      if(antidf(1,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(3,i+1,j,km) .GT. 0.) then
      a2= min(a2,Rplus (i+1,j,km))
	else
      a2= min(a2,Rminus(i+1,j,km))
      end if
      if(antidf(5,i+1,j-1,km) .GT. 0.) then
      a3= min(a3,Rplus (i+1,j-1,km))
	else
      a3= min(a3,Rminus(i+1,j-1,km))
      end if

	end if ! k>1

	alpha_e(1)= min(a1,a2,a3)


*     Element #2

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(2,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(4,i+1,j-1,k) .GT. 0.) then
      a2= Rplus (i+1,j-1,k)
	else
      a2= Rminus(i+1,j-1,k)
      end if
      if(antidf(6,i,j-1,k) .GT. 0.) then
      a3= Rplus (i,j-1,k)
	else
      a3= Rminus(i,j-1,k)
      end if

      if(k.GT.1) then

      km= k-1

      if(antidf(2,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(4,i+1,j-1,km) .GT. 0.) then
      a2= min(a2,Rplus (i+1,j-1,km))
	else
      a2= min(a2,Rminus(i+1,j-1,km))
      end if
      if(antidf(6,i,j-1,k) .GT. 0.) then
      a3= min(a3,Rplus (i,j-1,km))
	else
      a3= min(a3,Rminus(i,j-1,km))
      end if
	end if ! k>1

	alpha_e(2)= min(a1,a2,a3)


*     Element #3

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(3,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(5,i,j-1,k) .GT. 0.) then
      a2= Rplus (i,j-1,k)
	else
      a2= Rminus(i,j-1,k)
      end if
      if(antidf(1,i-1,j,k) .GT. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      if(k.GT.1) then

      km= k-1

      if(antidf(3,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(5,i,j-1,km) .GT. 0.) then
      a2= min(a2,Rplus (i,j-1,km))
	else
      a2= min(a2,Rminus(i,j-1,km))
      end if
      if(antidf(1,i-1,j,k) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if
	end if ! k>1

	alpha_e(3)= min(a1,a2,a3)


*     Element #4

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(4,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(2,i-1,j+1,k) .GT. 0.) then
      a2= Rplus (i-1,j+1,k)
	else
      a2= Rminus(i-1,j+1,k)
      end if
      if(antidf(6,i-1,j,k) .GT. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      if(k.GT.1) then

      km= k-1

      if(antidf(4,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(2,i-1,j+1,km) .GT. 0.) then
      a2= min(a2,Rplus (i-1,j+1,km))
	else
      a2= min(a2,Rminus(i-1,j+1,km))
      end if
      if(antidf(6,i-1,j,k) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if
	end if ! k>1

	alpha_e(4)= min(a1,a2,a3)


*     Element #5

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(5,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(3,i,j+1,k) .GT. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(1,i-1,j+1,k) .GT. 0.) then
      a3= Rplus (i-1,j+1,k)
	else
      a3= Rminus(i-1,j+1,k)
      end if

      if(k.GT.1) then

      km= k-1

      if(antidf(5,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(3,i,j+1,km) .GT. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(1,i-1,j+1,k) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j+1,km))
	else
      a3= min(a3,Rminus(i-1,j+1,km))
      end if
	end if ! k>1

	alpha_e(5)= min(a1,a2,a3)


*     Element #6

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(6,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(2,i,j+1,k) .GT. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(4,i+1,j,k) .GT. 0.) then
      a3= Rplus (i+1,j,k)
	else
      a3= Rminus(i+1,j,k)
      end if

      if(k.GT.1) then

      km= k-1

      if(antidf(6,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(2,i,j+1,km) .GT. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(4,i+1,j,k) .GT. 0.) then
      a3= min(a3,Rplus (i+1,j,km))
	else
      a3= min(a3,Rminus(i+1,j,km))
      end if
	end if ! k>1

	alpha_e(6)= min(a1,a2,a3)


*     Element #7

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(7,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(9,i+1,j,k) .GT. 0.) then
      a2= Rplus (i+1,j,k)
	else
      a2= Rminus(i+1,j,k)
      end if
      if(antidf(11,i+1,j-1,k) .GT. 0.) then
      a3= Rplus (i+1,j-1,k)
	else
      a3= Rminus(i+1,j-1,k)
      end if

      if(k.LT.kb) then

      km= k+1

      if(antidf(7,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(9,i+1,j,km) .GT. 0.) then
      a2= min(a2,Rplus (i+1,j,km))
	else
      a2= min(a2,Rminus(i+1,j,km))
      end if
      if(antidf(11,i+1,j-1,km) .GT. 0.) then
      a3= min(a3,Rplus (i+1,j-1,km))
	else
      a3= min(a3,Rminus(i+1,j-1,km))
      end if

	end if ! k>1

	alpha_e(7)= min(a1,a2,a3)


*     Element #8

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(8,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(10,i+1,j-1,k) .GT. 0.) then
      a2= Rplus (i+1,j-1,k)
	else
      a2= Rminus(i+1,j-1,k)
      end if
      if(antidf(12,i,j-1,k) .GT. 0.) then
      a3= Rplus (i,j-1,k)
	else
      a3= Rminus(i,j-1,k)
      end if

      if(k.LT.kb) then

      km= k+1

      if(antidf(8,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(10,i+1,j-1,km) .GT. 0.) then
      a2= min(a2,Rplus (i+1,j-1,km))
	else
      a2= min(a2,Rminus(i+1,j-1,km))
      end if
      if(antidf(12,i,j-1,k) .GT. 0.) then
      a3= min(a3,Rplus (i,j-1,km))
	else
      a3= min(a3,Rminus(i,j-1,km))
      end if
	end if ! k>1

	alpha_e(8)= min(a1,a2,a3)


*     Element #9

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(9,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(11,i,j-1,k) .GT. 0.) then
      a2= Rplus (i,j-1,k)
	else
      a2= Rminus(i,j-1,k)
      end if
      if(antidf(7,i-1,j,k) .GT. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      if(k.LT.kb) then

      km= k+1

      if(antidf(9,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(11,i,j-1,km) .GT. 0.) then
      a2= min(a2,Rplus (i,j-1,km))
	else
      a2= min(a2,Rminus(i,j-1,km))
      end if
      if(antidf(7,i-1,j,k) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if
	end if ! k>1

	alpha_e(9)= min(a1,a2,a3)


*     Element #10

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(10,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(8,i-1,j+1,k) .GT. 0.) then
      a2= Rplus (i-1,j+1,k)
	else
      a2= Rminus(i-1,j+1,k)
      end if
      if(antidf(12,i-1,j,k) .GT. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      if(k.LT.kb) then

      km= k+1

      if(antidf(10,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(8,i-1,j+1,km) .GT. 0.) then
      a2= min(a2,Rplus (i-1,j+1,km))
	else
      a2= min(a2,Rminus(i-1,j+1,km))
      end if
      if(antidf(12,i-1,j,k) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if
	end if ! k>1

	alpha_e(10)= min(a1,a2,a3)


*     Element #11

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(11,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(9,i,j+1,k) .GT. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(7,i-1,j+1,k) .GT. 0.) then
      a3= Rplus (i-1,j+1,k)
	else
      a3= Rminus(i-1,j+1,k)
      end if

      if(k.LT.kb) then

      km= k+1

      if(antidf(11,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(9,i,j+1,km) .GT. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(7,i-1,j+1,k) .GT. 0.) then
      a3= min(a3,Rplus (i-1,j+1,km))
	else
      a3= min(a3,Rminus(i-1,j+1,km))
      end if
	end if ! k>1

	alpha_e(11)= min(a1,a2,a3)


*     Element #12

	a1=1.0
	a2=1.0
	a3=1.0

      if(antidf(12,i,j,k) .GT. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(8,i,j+1,k) .GT. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(10,i+1,j,k) .GT. 0.) then
      a3= Rplus (i+1,j,k)
	else
      a3= Rminus(i+1,j,k)
      end if

      if(k.LT.kb) then

      km= k+1

      if(antidf(12,i,j,km) .GT. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(8,i,j+1,km) .GT. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(10,i+1,j,k) .GT. 0.) then
      a3= min(a3,Rplus (i+1,j,km))
	else
      a3= min(a3,Rminus(i+1,j,km))
      end if
	end if ! k>1

	alpha_e(12)= min(a1,a2,a3)


*----------------------------------------------------------------
      do nelem=1,12

      Sclr(i,j,k)= Sclr(i,j,k) + alpha_e(nelem)*antidf(nelem,i,j,k)

	end do
*----------------------------------------------------------------

	end do   ! k

c     This part if only hi-res solution at the boundary is not oscillating

c	ELSE

c	if(nt3(i,j,1) .LT. 0) then
c	Sclr(i,j,:)= Sclr1(i,j,:) + Sclr2(i,j,:)
c	end if   ! nt3 < 0

      END IF   ! nt3 > 0
	end do
	end do

2000  RETURN
      END

      subroutine transp(FXTK,FYTK,i,j,k,n,u,v,t,KT,il,jl,kl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl),
     *          t(0:il,0:jl,kl), KT(6,13)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

*     Version 24.08.2012.

      FXTK= 0.
      FYTK= 0.
*---------------     Triangle # 1.  ---------------------------------
      FXTK=FXTK+(U(I,J,K)+U(I+1,J,K)+U(I+1,J-1,K))*
     *          (T(I  ,J  ,K)+T(I+1,J  ,K)+T(I+1,J-1,K))
     * *KT(1,n)
*---------------     Triangle # 2.  -----------------------------------
      FYTK=FYTK-S2*(V(I,J,K)+V(I,J-1,K)+V(I+1,J-1,K))*
     *             (T(I  ,J  ,K)+T(I  ,J-1,K)+T(I+1,J-1,K))
     *  *KT(2,n)
*---------------     Triangle # 3.  -------------------------------
      FXTK=FXTK-(U(I,J,K)+U(I-1,J,K)+U(I,J-1,K))*
     *          (T(I  ,J  ,K)+T(I-1,J  ,K)+T(I  ,J-1,K))
     *  *KT(3,n)
      FYTK=FYTK-S3*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K))*
     *             (T(I  ,J  ,K)+T(I  ,J-1,K)+T(I-1,J  ,K))
     *  *KT(3,n)
*---------------     Triangle # 4.  ---------------------------------
      FXTK=FXTK-(U(I,J,K)+U(I-1,J,K)+U(I-1,J+1,K))*
     *          (T(I  ,J  ,K)+T(I-1,J  ,K)+T(I-1,J+1,K))
     *  *KT(4,n)
*---------------     Triangle # 5.  ----------------------------------
      FYTK=FYTK+S5*(V(I,J,K)+V(I,J+1,K)+V(I-1,J+1,K))*
     *             (T(I  ,J  ,K)+T(I  ,J+1,K)+T(I-1,J+1,K))
     *  *KT(5,n)
*---------------     Triangle # 6.  -----------------------------------
      FXTK=FXTK+(U(I,J,K)+U(I+1,J,K)+U(I,J+1,K))*
     *          (T(I  ,J  ,K)+T(I+1,J  ,K)+T(I  ,J+1,K))
     *  *KT(6,n)
      FYTK=FYTK+S6*(V(I,J,K)+V(I,J+1,K)+V(I+1,J,K))*
     *             (T(I  ,J  ,K)+T(I  ,J+1,K)+T(I+1,J  ,K))
     *  *KT(6,n)
      return
      end

	subroutine TSbc(Tdamp_In,Tdamp_Out,
     *           dt,T,Tm2,Tobs,u,v,w, nt3,Si,il,jl,kl,hx,hy,hz,R)

*     Version 14/11/2006 with vertical advection

      dimension T(0:il,0:jl,kl), Tm2(0:il,0:jl,kl), hz(kl),
     &          Tobs(0:il,0:jl,kl), nt3(0:il,0:jl,kl),Si(0:jl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)

	do j=1,jl-1
	S0=Si(j)
      do i=1,il-1
	do k=1,kl

	np=nt3(i,j,k)

      FZ = 0.

c      if(np .NE. 0) then
c	if( k.LT. kl) then
c      FZ= FZ + min(0., w(i,j,k+1))*(Tm2(i,j,k+1)-Tm2(i,j,k))/hz(k)
c	end if
c	if (k.GT.1) then
c      FZ= FZ + max(0., w(i,j,k  ))*(Tm2(i,j,k)-Tm2(i,j,k-1))/hz(k-1)
c	end if

c	end if  ! np .ne. 0

	if(np.EQ.-6) then

      delta_t=(T(i+1,j,k)-Tm2(i+1,j,k))
	delta_x=(T(i+1,j,k)-T(i+2,j,k))/(hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))/(hy)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))/(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/(hy)

c	umean= (2.*u(i,j,k)+2.*u(i+1,j,k)+u(i+1,j-1,k)+u(i,j+1,k))/6.
c	vmean= (2.*v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)
c     &+v(i+1,j-1,k)+v(i,j+1,k))/6.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-6

	if(np.EQ.-9) then

      delta_t=(T(i,j-1,k)-Tm2(i,j-1,k))
	delta_x=(T(i,j-1,k)-T(i,j-2,k))/(hy)
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))/(hx*S0)
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx +dt*v(i,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-9

	if(np.EQ.-8) then

      delta_t=(T(i-1,j,k)-Tm2(i-1,j,k))
	delta_x=(T(i-1,j,k)-T(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k)) /(hy)
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/ (hy) 

c	Rx= Rx +dt*u(i-1,j,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-8

	if(np.EQ.-7) then

      delta_t=(T(i,j+1,k)-Tm2(i,j+1,k))
	delta_x=(T(i,j+1,k)-T(i,j+2,k)) /(hy)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(hx*S0)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx -dt*v(i,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-7

	if(np.EQ.-10) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to W

      delta_t=(T(i+1,j,k)-Tm2(i+1,j,k))
	delta_x=(T(i+1,j,k)-T(i+2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k)) /(hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /(hy)      

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/3.
c	vmean= (v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/3.

c      Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Rx +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j-1,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j-1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=(T(i,j-1,k)-Tm2(i,j-1,k))
	delta_x=(T(i,j-1,k)-T(i,j-2,k)) /(hy) 
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx +dt*v(i+1,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i+1,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To S

	end if  ! n=-10


	if(np.EQ.-11) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to E
      delta_t=(T(i-1,j,k)-Tm2(i-1,j,k))
	delta_x=(T(i-1,j,k)-T(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))/(hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/(hy)

c	Rx= Rx +dt*u(i-1,j-1,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j-1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=(T(i,j-1,k)-Tm2(i,j-1,k))
	delta_x=(T(i,j-1,k)-T(i,j-2,k)) /(hy) 
	sss= delta_t*(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k)) /(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx +dt*v(i-1,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i-1,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j-1,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
    
	end if  ! To S

	end if  ! n=-11

	if(np.EQ.-12) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to E
      delta_t=(T(i-1,j,k)-Tm2(i-1,j,k))
	delta_x=(T(i-1,j,k)-T(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/(hy)

c	Rx= Rx +dt*u(i-1,j+1,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j+1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i-1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=(T(i,j+1,k)-Tm2(i,j+1,k))
	delta_x=(T(i,j+1,k)-T(i,j+2,k)) /(hy) 
	sss= delta_t*(Tm2(i,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hy)
	Ry=-delta_t*delta_y*D /(hx*S0)

c	Rx= Rx -dt*v(i-1,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i-1,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To N

	end if  ! n=-12

	if(np.EQ.-13) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to W

      delta_t=(T(i+1,j,k)-Tm2(i+1,j,k))
	delta_x=(T(i+1,j,k)-T(i+2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k)) /(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /(hy) 

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))/3.
c	vmean= (v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))/3.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j+1,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j+1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i+1,j,k)
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=(T(i,j+1,k)-Tm2(i,j+1,k))
	delta_x=(T(i,j+1,k)-T(i,j+2,k)) /(hy) 
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i,j+1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hy)
	Ry=-delta_t*delta_y*D /(hx*S0)

c	Rx= Rx -dt*v(i+1,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i+1,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*T(i,j+1,k)
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
	end if

	else
      T(i,j,k)=Tm2(i,j,k) +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if

      end if  ! To N
	end if  ! n=-13

	end do
	end do
	end do

	Return
	End

	subroutine TSbc_inc(Tdamp_In,Tdamp_Out,
     *           dt,T,Tm2,Tobs,u,v,w, nt3,Si,il,jl,kl,hx,hy,hz,R)

*     T is the time increment of the scalar 

*     Version 21/08/2012 with vertical advection

      dimension T(0:il,0:jl,kl), Tm2(0:il,0:jl,kl), hz(kl),
     &          Tobs(0:il,0:jl,kl), nt3(0:il,0:jl,kl),Si(0:jl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)

	do j=1,jl-1
	S0=Si(j)
      do i=1,il-1
	do k=1,kl

	np=nt3(i,j,k)

      if(np .NE. 0) then

      FZ = 0.
c	if( k.LT. kl) then
c      FZ= FZ + min(0., w(i,j,k+1))*(Tm2(i,j,k+1)-Tm2(i,j,k))/hz(k)
c	end if
c	if (k.GT.1) then
c      FZ= FZ + max(0., w(i,j,k  ))*(Tm2(i,j,k)-Tm2(i,j,k-1))/hz(k-1)
c	end if

	end if  ! np .ne. 0

	if(np.EQ.-6) then

      delta_t=T(i+1,j,k)
	delta_x=(T(i+1,j,k)+Tm2(i+1,j,k)
     *        -T(i+2,j,k)-Tm2(i+2,j,k))/(R*hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))/(R*hy)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))/(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/(R*hy)

c	umean= (2.*u(i,j,k)+2.*u(i+1,j,k)+u(i+1,j-1,k)+u(i,j+1,k))/6.
c	vmean= (2.*v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)
c     &+v(i+1,j-1,k)+v(i,j+1,k))/6.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j,k)/(R*hy)
	Rx= Rx -dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i+1,j,k)+Tm2(i+1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i+1,j,k)+Tm2(i+1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-6

	if(np.EQ.-9) then

      delta_t= T(i,j-1,k)
	delta_x=(T(i,j-1,k)+Tm2(i,j-1,k)
     *        -T(i,j-2,k)-Tm2(i,j-2,k))/(R*hy)
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))/(R*hx*S0)
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

c	Rx= Rx +dt*v(i,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i,j-1,k)/(R*hx*S0)
	Rx= Rx +dt*v(i,j,k)/(R*hy)
	Ry= Ry +dt*u(i,j,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j-1,k)+Tm2(i,j-1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j-1,k)+Tm2(i,j-1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-9

	if(np.EQ.-8) then

      delta_t=T(i-1,j,k)
	delta_x=(T(i-1,j,k)+Tm2(i-1,j,k)
     *        -T(i-2,j,k)-Tm2(i-2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k)) /(R*hy)
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/ (R*hy) 

c	Rx= Rx +dt*u(i-1,j,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j,k)/(R*hy)
	Rx= Rx +dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i-1,j,k)+Tm2(i-1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i-1,j,k)+Tm2(i-1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-8

	if(np.EQ.-7) then

      delta_t=T(i,j+1,k)
	delta_x=(T(i,j+1,k)+Tm2(i,j+1,k)
     *        -T(i,j+2,k)-Tm2(i,j+2,k)) /(R*hy)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(R*hx*S0)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

c	Rx= Rx -dt*v(i,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i,j+1,k)/(R*hx*S0)
	Rx= Rx -dt*v(i,j,k)/(R*hy)
	Ry= Ry -dt*u(i,j,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j+1,k)+Tm2(i,j+1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j+1,k)+Tm2(i,j+1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-7

	if(np.EQ.-10) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to W

      delta_t=T(i+1,j,k)
	delta_x=(T(i+1,j,k)+Tm2(i+1,j,k)
     *        -T(i+2,j,k)-Tm2(i+2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k)) /(R*hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hx*S0)
	Ry=-delta_t*delta_y*D /(R*hy)      

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/3.
c	vmean= (v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/3.

c      Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Rx +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j-1,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j-1,k)/(R*hy)
	Rx= Rx -dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i+1,j,k)+Tm2(i+1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i+1,j,k)+Tm2(i+1,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=T(i,j-1,k)
	delta_x=(T(i,j-1,k)+Tm2(i,j-1,k)
     *        -T(i,j-2,k)-Tm2(i,j-2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

c	Rx= Rx +dt*v(i+1,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i+1,j-1,k)/(R*hx*S0)
	Rx= Rx +dt*v(i,j,k)/(R*hy)
	Ry= Ry +dt*u(i,j,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j-1,k)+Tm2(i,j-1,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j-1,k)+Tm2(i,j-1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To S

	end if  ! n=-10


	if(np.EQ.-11) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to E
      delta_t=T(i-1,j,k)
	delta_x=(T(i-1,j,k)+Tm2(i-1,j,k)
     *        -T(i-2,j,k)-Tm2(i-2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))/(R*hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/(R*hy)

c	Rx= Rx +dt*u(i-1,j-1,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j-1,k)/(R*hy)
	Rx= Rx +dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i-1,j,k)+Tm2(i-1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i-1,j,k)+Tm2(i-1,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=T(i,j-1,k)
	delta_x=(T(i,j-1,k)+Tm2(i,j-1,k)
     *        -T(i,j-2,k)-Tm2(i,j-2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k)) /(R*hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

c	Rx= Rx +dt*v(i-1,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i-1,j-1,k)/(R*hx*S0)
	Rx= Rx +dt*v(i,j,k)/(R*hy)
	Ry= Ry +dt*u(i,j,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j-1,k)+Tm2(i,j-1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j-1,k)+Tm2(i,j-1,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
    
	end if  ! To S

	end if  ! n=-11

	if(np.EQ.-12) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to E
      delta_t=T(i-1,j,k)
	delta_x=(T(i-1,j,k)+Tm2(i-1,j,k)
     *        -T(i-2,j,k)-Tm2(i-2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hx*S0)
	Ry=-delta_t*delta_y*D/(R*hy)

c	Rx= Rx +dt*u(i-1,j+1,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j+1,k)/(R*hy)
	Rx= Rx +dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i-1,j,k)+Tm2(i-1,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i-1,j,k)+Tm2(i-1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=T(i,j+1,k)
	delta_x=(T(i,j+1,k)+Tm2(i,j+1,k)
     *        -T(i,j+2,k)-Tm2(i,j+2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(R*hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hy)
	Ry=-delta_t*delta_y*D /(R*hx*S0)

c	Rx= Rx -dt*v(i-1,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i-1,j+1,k)/(R*hx*S0)
	Rx= Rx -dt*v(i,j,k)/(R*hy)
	Ry= Ry -dt*u(i,j,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j+1,k)+Tm2(i,j+1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j+1,k)+Tm2(i,j+1,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To N

	end if  ! n=-12

	if(np.EQ.-13) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to W

      delta_t=T(i+1,j,k)
	delta_x=(T(i+1,j,k)+Tm2(i+1,j,k)
     *        -T(i+2,j,k)-Tm2(i+2,j,k)) /(R*hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k)) /(R*hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hx*S0)
	Ry=-delta_t*delta_y*D /(R*hy) 

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))/3.
c	vmean= (v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))/3.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j+1,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j+1,k)/(R*hy)
	Rx= Rx -dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i+1,j,k)+Tm2(i+1,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i+1,j,k)+Tm2(i+1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=T(i,j+1,k)
	delta_x=(T(i,j+1,k)+Tm2(i,j+1,k)
     *        -T(i,j+2,k)-Tm2(i,j+2,k)) /(R*hy) 
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i,j+1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(R*hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(R*hy)
	Ry=-delta_t*delta_y*D /(R*hx*S0)

c	Rx= Rx -dt*v(i+1,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i+1,j+1,k)/(R*hx*S0)
	Rx= Rx -dt*v(i,j,k)/(R*hy)
	Ry= Ry -dt*u(i,j,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j+1,k)+Tm2(i,j+1,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(T(i,j+1,k)+Tm2(i,j+1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if

      end if  ! To N
	end if  ! n=-13

	end do
	end do
	end do

	Return
	End

	subroutine TSbc_inc_iter(i,j,k,Tdamp_In,Tdamp_Out,
     *           dt,T,Tm1,Tm2,Tobs,u,v,w, nt3,Si,il,jl,kl,hx,hy,hz,R)

*     T is the time increment of the scalar 

*     Version 21/08/2012 with vertical advection

      dimension T(0:il,0:jl,kl),Tm1(0:il,0:jl,kl),Tm2(0:il,0:jl,kl),
     &          hz(kl),
     &          Tobs(0:il,0:jl,kl), nt3(0:il,0:jl,kl),Si(0:jl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)

	S0=Si(j)

	np=nt3(i,j,k)

      if(np .NE. 0) then

      FZ = 0.
c	if( k.LT. kl) then
c      FZ= FZ + min(0., w(i,j,k+1))*(Tm2(i,j,k+1)-Tm2(i,j,k))/hz(k)
c	end if
c	if (k.GT.1) then
c      FZ= FZ + max(0., w(i,j,k  ))*(Tm2(i,j,k)-Tm2(i,j,k-1))/hz(k-1)
c	end if

	end if  ! np .ne. 0

	if(np.EQ.-6) then

      delta_t=Tm1(i+1,j,k)
	delta_x=(Tm1(i+1,j,k)+Tm2(i+1,j,k)
     *        -Tm1(i+2,j,k)-Tm2(i+2,j,k))/(hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))/(hy)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))/(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/(hy)

c	umean= (2.*u(i,j,k)+2.*u(i+1,j,k)+u(i+1,j-1,k)+u(i,j+1,k))/6.
c	vmean= (2.*v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)
c     &+v(i+1,j-1,k)+v(i,j+1,k))/6.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

	Rx= Rx -dt*u(i+1,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i+1,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i+1,j,k)+Tm2(i+1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i+1,j,k)+Tm2(i+1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-6

	if(np.EQ.-9) then

      delta_t= Tm1(i,j-1,k)
	delta_x=(Tm1(i,j-1,k)+Tm2(i,j-1,k)
     *        -Tm1(i,j-2,k)-Tm2(i,j-2,k))/(hy)
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))/(hx*S0)
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

	Rx= Rx +dt*v(i,j-1,k)/(R*hy)
	Ry= Ry +dt*u(i,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j-1,k)+Tm2(i,j-1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j-1,k)+Tm2(i,j-1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-9

	if(np.EQ.-8) then

      delta_t=Tm1(i-1,j,k)
	delta_x=(Tm1(i-1,j,k)+Tm2(i-1,j,k)
     *        -Tm1(i-2,j,k)-Tm2(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k)) /(hy)
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/ (hy) 

	Rx= Rx +dt*u(i-1,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i-1,j,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i-1,j,k)+Tm2(i-1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i-1,j,k)+Tm2(i-1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-8

	if(np.EQ.-7) then

      delta_t=Tm1(i,j+1,k)
	delta_x=(Tm1(i,j+1,k)+Tm2(i,j+1,k)
     *        -Tm1(i,j+2,k)-Tm2(i,j+2,k)) /(hy)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(hx*S0)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

	Rx= Rx -dt*v(i,j+1,k)/(R*hy)
	Ry= Ry -dt*u(i,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j+1,k)+Tm2(i,j+1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j+1,k)+Tm2(i,j+1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
	end if  ! n=-7

	if(np.EQ.-10) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to W

      delta_t=Tm1(i+1,j,k)
	delta_x=(Tm1(i+1,j,k)+Tm2(i+1,j,k)
     *        -Tm1(i+2,j,k)-Tm2(i+2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k)) /(hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /(hy)      

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/3.
c	vmean= (v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/3.

c      Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Rx +dt*vmean/(R*hy)

	Rx= Rx -dt*u(i+1,j-1,k)/(R*hx*S0)
	Ry= Ry +dt*v(i+1,j-1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i+1,j,k)+Tm2(i+1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i+1,j,k)+Tm2(i+1,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=Tm1(i,j-1,k)
	delta_x=(Tm1(i,j-1,k)+Tm2(i,j-1,k)
     *        -Tm1(i,j-2,k)-Tm2(i,j-2,k)) /(hy) 
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(R*hy)
	Ry=-delta_t*delta_y*D/(R*hx*S0)

	Rx= Rx +dt*v(i+1,j-1,k)/(R*hy)
	Ry= Ry +dt*u(i+1,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j-1,k)+Tm2(i,j-1,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j-1,k)+Tm2(i,j-1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To S

	end if  ! n=-10


	if(np.EQ.-11) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j-1,k).LT.0)then !transport to E
      delta_t=Tm1(i-1,j,k)
	delta_x=(Tm1(i-1,j,k)+Tm2(i-1,j,k)
     *        -Tm1(i-2,j,k)-Tm2(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))/(hy)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/(hy)

	Rx= Rx +dt*u(i-1,j-1,k)/(R*hx*S0)
	Ry= Ry -dt*v(i-1,j-1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i-1,j,k)+Tm2(i-1,j,k))
     &                      +Ry*(Tm2(i,j-1,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i-1,j,k)+Tm2(i-1,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j-1,k).GT.0)then !transport to S

      delta_t=Tm1(i,j-1,k)
	delta_x=(Tm1(i,j-1,k)+Tm2(i,j-1,k)
     *        -Tm1(i,j-2,k)-Tm2(i,j-2,k)) /(hy) 
	sss= delta_t*(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k)) /(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hy)
	Ry=-delta_t*delta_y*D/(hx*S0)

	Rx= Rx +dt*v(i-1,j-1,k)/(R*hy)
	Ry= Ry +dt*u(i-1,j-1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j-1,k)+Tm2(i,j-1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j-1,k)+Tm2(i,j-1,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
    
	end if  ! To S

	end if  ! n=-11

	if(np.EQ.-12) then

	if(nt3(i-1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to E
      delta_t=Tm1(i-1,j,k)
	delta_x=(Tm1(i-1,j,k)+Tm2(i-1,j,k)
     *        -Tm1(i-2,j,k)-Tm2(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/(hy)

	Rx= Rx +dt*u(i-1,j+1,k)/(R*hx*S0)
	Ry= Ry -dt*v(i-1,j+1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)
	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i-1,j,k)+Tm2(i-1,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i-1,j,k)+Tm2(i-1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To E

	if(nt3(i-1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=Tm1(i,j+1,k)
	delta_x=(Tm1(i,j+1,k)+Tm2(i,j+1,k)
     *        -Tm1(i,j+2,k)-Tm2(i,j+2,k)) /(hy) 
	sss= delta_t*(Tm2(i,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hy)
	Ry=-delta_t*delta_y*D /(hx*S0)

	Rx= Rx -dt*v(i-1,j+1,k)/(R*hy)
	Ry= Ry -dt*u(i-1,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j+1,k)+Tm2(i,j+1,k))
     &                      +Ry*(Tm2(i-1,j,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j+1,k)+Tm2(i,j+1,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
      
	end if  ! To N

	end if  ! n=-12

	if(np.EQ.-13) then

	if(nt3(i+1,j,k).GT.0.AND.nt3(i,j+1,k).LT.0)then !transport to W

      delta_t=Tm1(i+1,j,k)
	delta_x=(Tm1(i+1,j,k)+Tm2(i+1,j,k)
     *        -Tm1(i+2,j,k)-Tm2(i+2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k)) /(hy)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /(hy) 

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))/3.
c	vmean= (v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))/3.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

	Rx= Rx -dt*u(i+1,j+1,k)/(R*hx*S0)
	Ry= Ry +dt*v(i+1,j+1,k)/(R*hy)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i+1,j,k)+Tm2(i+1,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i+1,j,k)+Tm2(i+1,j,k))
     &                      -Ry*(Tm2(i,j+1,k)-Tm2(i,j,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if
     
	end if  ! To W

	if(nt3(i+1,j,k).LT.0.AND.nt3(i,j+1,k).GT.0)then !transport to N
      delta_t=Tm1(i,j+1,k)
	delta_x=(Tm1(i,j+1,k)+Tm2(i,j+1,k)
     *        -Tm1(i,j+2,k)-Tm2(i,j+2,k)) /(hy) 
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i,j+1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hy)
	Ry=-delta_t*delta_y*D /(hx*S0)

	Rx= Rx -dt*v(i+1,j+1,k)/(R*hy)
	Ry= Ry -dt*u(i+1,j+1,k)/(R*hx*S0)

      if(Rx .GE. 0.)then
	D=1./(1.+Rx)

	if(Ry.GE.0.)then
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j+1,k)+Tm2(i,j+1,k))
     &           -dt*FZ +dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	else
	T(i,j,k)=D*(Tm2(i,j,k)+Rx*(Tm1(i,j+1,k)+Tm2(i,j+1,k))
     &                      -Ry*(Tm2(i+1,j,k)-Tm2(i,j,k))
     &            -dt*FZ+dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_Out)
     &         -Tm2(i,j,k)
	end if

	else
      T(i,j,k)= dt*(Tobs(i,j,k)-Tm2(i,j,k))/Tdamp_In
	end if

      end if  ! To N
	end if  ! n=-13


cccc	if(nt3(i,j,k) .LT. 0) T(i,j,k)= Tobs(i,j,k)-Tm2(i,j,k) !!!!!!!!!!!!!!

	Return
	End

      SUBROUTINE DIFFUW_TG(T,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &                     il,jl,kl,ilm,jlm,dt,KT,A,hx,hy,R)

*     Version 16.11.2010.

*     Galerkin artificial diffusion approximation by Taylor-Galerkin
*     A is the "physical" or rather "background" diffusivity.
*     Boundary conditions are taken into account in the Main Program.

      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)
      dimension T(0:il,0:jl,kl), KT(6,13), 
     *          hz(kl)
	real a11(6),a12(6),a22(6)
	real a13(6),a23(6),a33(6)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
	common /consts/ c3,c6,c9,c12,c18,Rhx,R2hx2
      real KT

	TAU   = 0.5*dt

	Rhy=Rhx/asr
	km=k-1
	kp=k+1	
      DIFT=0.
      DIFTP=0.
      Q1K=T(i,j,k)+T(i+1,j,k)+T(i+1,j-1,k)
      Q2K=T(i,j,k)+T(i+1,j-1,k)+T(i,j-1,k)
      Q3K=T(i,j,k)+T(i-1,j,k)+T(i,j-1,k)
      Q4K=T(i,j,k)+T(i-1,j,k)+T(i-1,j+1,k)
      Q5K=T(i,j,k)+T(i,j+1,k)+T(i-1,j+1,k)
      Q6K=T(i,j,k)+T(i+1,j,k)+T(i,j+1,k)

      IF( K .GT. 1) THEN
C     UPPER HALF

      Q1Km=T(i,j,km)+T(i+1,j,km)+T(i+1,j-1,km)
      Q2Km=T(i,j,km)+T(i+1,j-1,km)+T(i,j-1,km)
      Q3Km=T(i,j,km)+T(i-1,j,km)+T(i,j-1,km)
      Q4Km=T(i,j,km)+T(i-1,j,km)+T(i-1,j+1,km)
      Q5Km=T(i,j,km)+T(i,j+1,km)+T(i-1,j+1,km)
      Q6Km=T(i,j,km)+T(i+1,j,km)+T(i,j+1,km)

c     Triangle # 1.

      u1=c6*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)+
     +    u(i,j,km)+u(i+1,j,km)+u(i+1,j-1,km))
      v1=c6*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k)+
     +    v(i,j,km)+v(i+1,j,km)+v(i+1,j-1,km))
      w1=c3*(w(i,j,k)+w(i+1,j,k)+w(i+1,j-1,k))

      a11(1)= TAU*u1*u1 +A
      a12(1)= TAU*u1*v1 
      a22(1)= TAU*v1*v1 +A
      a13(1)= TAU*u1*w1
      a23(1)= TAU*w1*v1 
      a33(1)= TAU*w1*w1

c     -------------------------------

c     Triangle # 2.

      u2=c6*(u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k)+
     +    u(i,j,km)+u(i,j-1,km)+u(i+1,j-1,km))
      v2=c6*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k)+
     +    v(i,j,km)+v(i,j-1,km)+v(i+1,j-1,km))
      w2=c3*(w(i,j,k)+w(i,j-1,k)+w(i+1,j-1,k))

      a11(2)= TAU*u2*u2 +A
      a12(2)= TAU*u2*v2 
      a22(2)= TAU*v2*v2 +A
      a13(2)= TAU*u2*w2
      a23(2)= TAU*w2*v2 
      a33(2)= TAU*w2*w2

c     -------------------------------

c     Triangle # 3.

      u3=c6*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k)+
     +    u(i,j,km)+u(i-1,j,km)+u(i,j-1,km))
      v3=c6*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k)+
     +    v(i,j,km)+v(i-1,j,km)+v(i,j-1,km))
      w3=c3*(w(i,j,k)+w(i-1,j,k)+w(i,j-1,k))

      a11(3)= TAU*u3*u3 +A
      a12(3)= TAU*u3*v3 
      a22(3)= TAU*v3*v3 +A
      a13(3)= TAU*u3*w3
      a23(3)= TAU*w3*v3 
      a33(3)= TAU*w3*w3

c     -------------------------------

c     Triangle # 4.

      u4=c6*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k)+
     +    u(i,j,km)+u(i-1,j,km)+u(i-1,j+1,km))
      v4=c6*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k)+
     +    v(i,j,km)+v(i-1,j,km)+v(i-1,j+1,km))
      w4=c3*(w(i,j,k)+w(i-1,j,k)+w(i-1,j+1,k))

      a11(4)= TAU*u4*u4 +A
      a12(4)= TAU*u4*v4 
      a22(4)= TAU*v4*v4 +A
      a13(4)= TAU*u4*w4
      a23(4)= TAU*w4*v4 
      a33(4)= TAU*w4*w4

c     -------------------------------

c     Triangle # 5.

      u5=c6*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k)+
     +    u(i,j,km)+u(i-1,j+1,km)+u(i,j+1,km))
      v5=c6*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)+
     +    v(i,j,km)+v(i-1,j+1,km)+v(i,j+1,km))
      w5=c3*(w(i,j,k)+w(i-1,j+1,k)+w(i,j+1,k))

      a11(5)= TAU*u5*u5 +A
      a12(5)= TAU*u5*v5 
      a22(5)= TAU*v5*v5 +A
      a13(5)= TAU*u5*w5
      a23(5)= TAU*w5*v5 
      a33(5)= TAU*w5*w5

c     -------------------------------

c     Triangle # 6.

      u6=c6*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k)+
     +    u(i,j,km)+u(i+1,j,km)+u(i,j+1,km))
      v6=c6*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k)+
     +    v(i,j,km)+v(i+1,j,km)+v(i,j+1,km))
      w6=c3*(w(i,j,k)+w(i+1,j,k)+w(i,j+1,k))

      a11(6)= TAU*u6*u6 +A
      a12(6)= TAU*u6*v6 
      a22(6)= TAU*v6*v6 +A
      a13(6)= TAU*u6*w6
      a23(6)= TAU*w6*v6 
      a33(6)= TAU*w6*w6

c     -------------------------------

c     Triangle # 1.

c     Horizontal part of the operator
      DIFT= DIFT +KT(1,n)*0.5*(
     &     a11(1)*(T(I+1,J,K)-T(I  ,J,  K))/S0
     &+asr*a12(1)*(T(i+1,j,k)-T(i+1,j-1,k))   )*hz(km)  

c     dT/dz*dfi/dlambda
      DIFT= DIFT +KT(1,n)*Rhx*c6*a13(1)*(Q1K-Q1Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(1,n)*Rhx*c12*a13(1)*(T(i+1,j,k )-T(i,j,k)+
     &                                   T(i+1,j,km)-T(i,j,km))

c     dfi/dz*dT/dteta
      DIFT=DIFT-KT(1,n)*Rhy*c12*S1*a23(1)*(T(i+1,j,k )-T(i+1,j-1,k)+
     &                                     T(i+1,j,km)-T(i+1,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(1,n)*R2hx2*S1*a33(1)*(Q1K-Q1Km)/hz(km)

c     Triangle # 2.

c     Horizontal part of the operator
      DIFT= DIFT -KT(2,n)*0.5*(
     &asr2*a22(2)*S2*(T(I,J,K)-T(I,J-1,K))
     &+asr*a12(2)*   (T(i+1,j-1,k)-T(i,j-1,k)))*hz(km)


c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(2,n)*Rhx*c12*a13(2)*(T(i+1,j-1,k )-T(i,j-1,k)+
     &                                    T(i+1,j-1,km)-T(i,j-1,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT -KT(2,n)*Rhy*S2*c6*a23(2)*asr*(Q2K-Q2Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(2,n)*Rhy*c12*S2*a23(2)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,km)-T(i,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(2,n)*R2hx2*S2*a33(2)*(Q2K-Q2Km)/hz(km)

c     Triangle # 3.

c     Horizontal part of the operator
      DIFT= DIFT +KT(3,n)*0.5*(
     &         a11(3)*(T(I-1,J,K)-T(I,J,K))/S0
     &    +asr*a12(3)*(T(i,j-1,k)-T(i  ,j,k))
     &    -asr*a12(3)*(T(i,j  ,k)-T(i-1,j,k))
     &-S3*asr2*a22(3)*(T(I,J  ,K)-T(I,J-1,K)))*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT -KT(3,n)*Rhx*c6*a13(3)*(Q3K-Q3Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(3,n)*Rhx*c12*a13(3)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,km)-T(i-1,j,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT -KT(3,n)*Rhy*S3*c6*a23(3)*asr*(Q3K-Q3Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(3,n)*Rhy*c12*S3*a23(3)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,km)-T(i,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(3,n)*R2hx2*S3*a33(3)*(Q3K-Q3Km)/hz(km)
c      write(*,*) 'point 4'

c     Triangle # 4.

c     Horizontal part of the operator
      DIFT= DIFT +KT(4,n)*0.5*(
     &     a11(4)*(T(I-1,J,K)-T(I,J,K))/S0
     &+asr*a12(4)*(T(i-1,j,k)-T(i-1,j+1,k)))*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT -KT(4,n)*Rhx*c6*a13(4)*(Q4K-Q4Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(4,n)*Rhx*c12*a13(4)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,km)-T(i-1,j,km))

c     dfi/dz*dT/dteta
      DIFT=DIFT-KT(4,n)*Rhy*c12*S4*a23(4)*(T(i-1,j+1,k )-T(i-1,j,k)+
     &                                     T(i-1,j+1,km)-T(i-1,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(4,n)*R2hx2*S4*a33(4)*(Q4K-Q4Km)/hz(km)

c     Triangle # 5.

c     Horizontal part of the operator
      DIFT= DIFT -KT(5,n)*0.5*(
     &asr2*S5*a22(5)*(T(I,J,K)-T(I,J+1,K))
     &   +asr*a12(5)*(T(i-1,j+1,k)-T(i,j+1,k)))*hz(km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(5,n)*Rhx*c12*a13(5)*(T(i,j+1,k )-T(i-1,j+1,k)+
     &                                    T(i,j+1,km)-T(i-1,j+1,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT +KT(5,n)*Rhy*S5*c6*a23(5)*asr*(Q5K-Q5Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(5,n)*Rhy*c12*S5*a23(5)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,km)-T(i,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(5,n)*R2hx2*S5*a33(5)*(Q5K-Q5Km)/hz(km)

c     Triangle # 6.

c     Horizontal part of the operator
      DIFT= DIFT +KT(6,n)*0.5*(
     &         a11(6)*(T(I+1,J,K)-T(I,J,K))/S0
     &    +asr*a12(6)*(T(i,j+1,k)-T(i,j,k))   
     &    -asr*a12(6)*(T(i,j,k)-T(i+1,j,k))
     *-asr2*S6*a22(6)*(T(I,J,K)-T(I,J+1,K)))*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT +KT(6,n)*Rhx*c6*a13(6)*(Q6K-Q6Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(6,n)*Rhx*c12*a13(6)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,km)-T(i,j,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT +KT(6,n)*Rhy*S6*c6*a23(6)*asr*(Q6K-Q6Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(6,n)*Rhy*c12*S6*a23(6)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,km)-T(i,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(6,n)*R2hx2*S6*a33(6)*(Q6K-Q6Km)/hz(km)

      END IF

      IF( K .LT. KB) THEN 
C     LOWER HALF

      Q1Kp=T(i,j,kp)+T(i+1,j,kp)+T(i+1,j-1,kp)
      Q2Kp=T(i,j,kp)+T(i+1,j-1,kp)+T(i,j-1,kp)
      Q3Kp=T(i,j,kp)+T(i-1,j,kp)+T(i,j-1,kp)
      Q4Kp=T(i,j,kp)+T(i-1,j,kp)+T(i-1,j+1,kp)
      Q5Kp=T(i,j,kp)+T(i,j+1,kp)+T(i-1,j+1,kp)
      Q6Kp=T(i,j,kp)+T(i+1,j,kp)+T(i,j+1,kp)


c     Triangle # 1.

      u1=c6*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)+
     +    u(i,j,kp)+u(i+1,j,kp)+u(i+1,j-1,kp))
      v1=c6*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k)+
     +    v(i,j,kp)+v(i+1,j,kp)+v(i+1,j-1,kp))
      w1=c3*(w(i,j,kp)+w(i+1,j,kp)+w(i+1,j-1,kp))

      a11(1)= TAU*u1*u1 +A
      a12(1)= TAU*u1*v1 
      a22(1)= TAU*v1*v1 +A
      a13(1)= TAU*u1*w1
      a23(1)= TAU*w1*v1 
      a33(1)= TAU*w1*w1

c     -------------------------------

c     Triangle # 2.

      u2=c6*(u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k)+
     +    u(i,j,kp)+u(i,j-1,kp)+u(i+1,j-1,kp))
      v2=c6*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k)+
     +    v(i,j,kp)+v(i,j-1,kp)+v(i+1,j-1,kp))
      w2=c3*(w(i,j,kp)+w(i,j-1,kp)+w(i+1,j-1,kp))

      a11(2)= TAU*u2*u2 +A
      a12(2)= TAU*u2*v2 
      a22(2)= TAU*v2*v2 +A
      a13(2)= TAU*u2*w2
      a23(2)= TAU*w2*v2 
      a33(2)= TAU*w2*w2

c     -------------------------------

c     Triangle # 3.

      u3=c6*(u(i,j,k)+u(i-1,j,k)+u(i,j-1,k)+
     +    u(i,j,kp)+u(i-1,j,kp)+u(i,j-1,kp))
      v3=c6*(v(i,j,k)+v(i-1,j,k)+v(i,j-1,k)+
     +    v(i,j,kp)+v(i-1,j,kp)+v(i,j-1,kp))
      w3=c3*(w(i,j,kp)+w(i-1,j,kp)+w(i,j-1,kp))

      a11(3)= TAU*u3*u3 +A
      a12(3)= TAU*u3*v3 
      a22(3)= TAU*v3*v3 +A
      a13(3)= TAU*u3*w3
      a23(3)= TAU*w3*v3 
      a33(3)= TAU*w3*w3

c     -------------------------------

c     Triangle # 4.

      u4=c6*(u(i,j,k)+u(i-1,j,k)+u(i-1,j+1,k)+
     +    u(i,j,kp)+u(i-1,j,kp)+u(i-1,j+1,kp))
      v4=c6*(v(i,j,k)+v(i-1,j,k)+v(i-1,j+1,k)+
     +    v(i,j,kp)+v(i-1,j,kp)+v(i-1,j+1,kp))
      w4=c3*(w(i,j,kp)+w(i-1,j,kp)+w(i-1,j+1,kp))

      a11(4)= TAU*u4*u4 +A
      a12(4)= TAU*u4*v4 
      a22(4)= TAU*v4*v4 +A
      a13(4)= TAU*u4*w4
      a23(4)= TAU*w4*v4 
      a33(4)= TAU*w4*w4

c     -------------------------------

c     Triangle # 5.

      u5=c6*(u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k)+
     +    u(i,j,kp)+u(i-1,j+1,kp)+u(i,j+1,kp))
      v5=c6*(v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)+
     +    v(i,j,kp)+v(i-1,j+1,kp)+v(i,j+1,kp))
      w5=c3*(w(i,j,kp)+w(i-1,j+1,kp)+w(i,j+1,kp))

      a11(5)= TAU*u5*u5 +A
      a12(5)= TAU*u5*v5 
      a22(5)= TAU*v5*v5 +A
      a13(5)= TAU*u5*w5
      a23(5)= TAU*w5*v5 
      a33(5)= TAU*w5*w5

c     -------------------------------

c     Triangle # 6.

      u6=c6*(u(i,j,k)+u(i+1,j,k)+u(i,j+1,k)+
     +    u(i,j,kp)+u(i+1,j,kp)+u(i,j+1,kp))
      v6=c6*(v(i,j,k)+v(i+1,j,k)+v(i,j+1,k)+
     +    v(i,j,kp)+v(i+1,j,kp)+v(i,j+1,kp))
      w6=c3*(w(i,j,kp)+w(i+1,j,kp)+w(i,j+1,kp))

      a11(6)= TAU*u6*u6 +A
      a12(6)= TAU*u6*v6 
      a22(6)= TAU*v6*v6 +A
      a13(6)= TAU*u6*w6
      a23(6)= TAU*w6*v6 
      a33(6)= TAU*w6*w6

c     -------------------------------

c     Triangle # 1.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(1,np)*(a11(1)*(T(I+1,J,K)-T(I,J,K))/S0
     &            +asr*a12(1)*(T(i+1,j,k)-T(i+1,j-1,k)))*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP +KT(1,np)*Rhx*c6*a13(1)*(Q1Kp-Q1K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(1,np)*Rhx*c12*a13(1)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,kp)-T(i,j,kp))

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(1,np)*Rhy*c12*S1*a23(1)*(T(i+1,j,k )-T(i+1,j-1,k)+
     &                                       T(i+1,j,kp)-T(i+1,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(1,np)*R2hx2*S1*a33(1)*(Q1K-Q1Kp)/hz(k)

c     Triangle # 2.

c     Horizontal part of the operator
      DIFTP= DIFTP-0.5*KT(2,np)*(a22(2)*S2*(T(I,J,K)-T(I,J-1,K))*ASR2
     &           +        asr*a12(2)*(T(i+1,j-1,k)-T(i,j-1,k)))*hz(k)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(2,np)*Rhx*c12*a13(2)*
     &                                   (T(i+1,j-1,k )-T(i,j-1,k)+
     &                                    T(i+1,j-1,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP -KT(2,np)*Rhy*S2*c6*a23(2)*asr*(Q2Kp-Q2K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(2,np)*Rhy*c12*S2*a23(2)*
     &                                      (T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(2,np)*R2hx2*S2*a33(2)*(Q2K-Q2Kp)/hz(k)

c     Triangle # 3.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(3,np)*(a11(3)*(T(I-1,J,K)-T(I,J,K))/S0
     &            +asr*a12(3)*(T(i,j-1,k)-T(i,j,k))
     &            -asr*a12(3)*(T(i,j,k)-T(i-1,j,k))
     &             -S3*a22(3)*(T(I,J,K)-T(I,J-1,K))*ASR2)*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP -KT(3,np)*Rhx*c6*a13(3)*(Q3Kp-Q3K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(3,np)*Rhx*c12*a13(3)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,kp)-T(i-1,j,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP -KT(3,np)*Rhy*S3*c6*a23(3)*asr*(Q3Kp-Q3K)

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(3,np)*Rhy*c12*S3*a23(3)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(3,np)*R2hx2*S3*a33(3)*(Q3K-Q3Kp)/hz(k)

c     Triangle # 4.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(4,np)*(a11(4)*(T(I-1,J,K)-T(I,J,K))/S0
     &            +asr*a12(4)*(T(i-1,j,k)-T(i-1,j+1,k)))*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP -KT(4,np)*Rhx*c6*a13(4)*(Q4Kp-Q4K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(4,np)*Rhx*c12*a13(4)*(T(i,j,k )-T(i-1,j,k)+
     &                                    T(i,j,kp)-T(i-1,j,kp))

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(4,np)*Rhy*c12*S4*a23(4)*
     &                                  (T(i-1,j+1,k )-T(i-1,j,k)+
     &                                   T(i-1,j+1,kp)-T(i-1,j,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(4,np)*R2hx2*S4*a33(4)*(Q4K-Q4Kp)/hz(k)


c     Triangle # 5.

c     Horizontal part of the operator
      DIFTP= DIFTP-0.5*KT(5,np)*(a22(5)*S5*(T(I,J,K)-T(I,J+1,K))*ASR2
     &           +        asr*a12(5)*(T(i-1,j+1,k)-T(i,j+1,k)))*hz(k)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(5,np)*Rhx*c12*a13(5)*
     &                                   (T(i,j+1,k )-T(i-1,j+1,k)+
     &                                    T(i,j+1,kp)-T(i-1,j+1,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP +KT(5,np)*Rhy*S5*c6*a23(5)*asr*(Q5Kp-Q5K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(5,np)*Rhy*c12*S5*a23(5)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,kp)-T(i,j,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(5,np)*R2hx2*S5*a33(5)*(Q5K-Q5Kp)/hz(k)


c     Triangle # 6.

c     Horizontal part of the operator
      DIFTP= DIFTP+0.5*KT(6,np)*(a11(6)*(T(I+1,J,K)-T(I,J,K))/S0
     &               +a12(6)*asr*(T(i,j+1,k)-T(i,j,k))   
     &               -a12(6)*asr*(T(i,j,k)-T(i+1,j,k))
     *               -S6*a22(6)*(T(I,J,K)-T(I,J+1,K))*ASR2)*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP +KT(6,np)*Rhx*c6*a13(6)*(Q6Kp-Q6K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(6,np)*Rhx*c12*a13(6)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,kp)-T(i,j,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP +KT(6,np)*Rhy*S6*c6*a23(6)*asr*(Q6Kp-Q6K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(6,np)*Rhy*c12*S6*a23(6)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,kp)-T(i,j,kp))

c     dT/dz*dfi/dz 
      DIFTP= DIFTP -KT(6,np)*R2hx2*S6*a33(6)*(Q6K-Q6Kp)/hz(k)

      END IF

      RETURN
      END
