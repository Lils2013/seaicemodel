      SUBROUTINE advect3D(Sclr,Sclr1,Sclr2,SclrObs,Qcascad,Agm)

*******************************************************
*
*     Sclr 3D trasport and Horizontal diffusion
*     FEM Streamline Upwind Scheme Taylor-Galerkin
*     FCT by Loehner, et.al.
*     R. L¨hner, K. Morgan, J. Peraire and M. Vahdati.
*     Finite element flux-corrected transport (FEM-FCT)
*     for the Euler and Navier-Stokes equations. 
*     Int. J. Numer. Meth. Fluids 7 (1987) 1093–1109.
*
*     Explicit time stepping - one step Taylor-Galerkin
*     Special treatment of open boundaries.
*     Version 13.10.2015
*
*******************************************************

      parameter (itermax=3, omega=1.0)   ! Iteration parameters
      parameter (gamma_fct=5.e-1)        ! Numerical scalar diffusivity
*     The last parameter may be variable, see 
*     G.E. Georghiou, R. Morrow and A.C. Metaxas, An improved finite-element
*     fluxcorrected transport algorithm. J. Comput. Phys. 148 (1999) 605–620.
*     The general idea is that gamma_fct=0.5(1-c)c, 
*     c = Courant number on the cluster, gamma_fct<=1/8=0.125

      INCLUDE 'Slo2.fi'

	dimension Sclr (0:il1,0:jl1,kl),Sclr1(0:il1,0:jl1,kl),
     *          Sclr2(0:il1,0:jl1,kl),SclrObs(0:il1,0:jl1,kl)  
      dimension Qcascad(il,jl,kl)

*     Auxiliary arrays for FCT calculations

	dimension Rplus(il,jl,kl), Rminus(il,jl,kl), 
     *          alpha_e(12)
      dimension antidf(12,il,jl,kl) !Antidif fluxes to node (i,j,k) 

	real Agm(2,0:il1,0:jl1)  ! Variable GMG coefficient on the triangle

      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      common /diffus/ al,alt
	common /consts/ c3,c6,c12,c18,Rhx,Rhy,R2hxhy
      INCLUDE 'Tparm.fi'

      Tdamp_In = dt             ! Time to restore at inflow
      Tdamp_Out= 30.*24.*3600.  ! and outflow points.

	c6= 1./6.
	c3= 1./3.
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
	R2hxhy=c3*Rhx*Rhy 
	
      A= 0.5/(Rhx*Rhy)


      COEFX=c18/Rhx
      COEFY=c18/Rhy


      RsT    =0.
      Sclr   =0.
      Sclr1  =0.
	Rplus  =0.
	Rminus =0.
	antidf =0.


*---------------------------------------------------
*         Right Hand Side with upwinding
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
      KB = KM2(I,J)
      IF(nt3(i,j,1).GT.0) THEN

*-------------------------------------------------
*         Surface and Nonbottom points
*-------------------------------------------------

      DO K=1, KB -1
      
	hzk = hz(k)
	km= k-1
	kp= k+1
	
      n= nt3(i,j,k)
      np=nt3(i,j,kp)

      IF( k .EQ. 1) THEN
        hzk1= 0.
      ELSE
        hzk1= hz(km)
      END IF
	 
      call DIFFUW_TG(Sclr2,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,dt,KT,ALT)

      if( k .GT. 1) then 
      call transp(FXTK ,FYTK ,i,j,k ,n,u,v,Sclr2,KT,il1,jl1,kl)
      call transp(FXTKm,FYTKm,i,j,km,n,u,v,Sclr2,KT,il1,jl1,kl)
      else
      FXTK =0.
      FYTK =0.
      FXTKm =0.
      FYTKm =0.
      end if
      call transp(FXTK1 ,FYTK1 ,i,j,k ,np,u,v,Sclr2,KT,il1,jl1,kl)
      call transp(FXTK1p,FYTK1p,i,j,kp,np,u,v,Sclr2,KT,il1,jl1,kl)

      FXTK =COEFX*FXTK
      FYTK =COEFY*FYTK
      FXTKm =COEFX*FXTKm
      FYTKm =COEFY*FYTKm
      FXTK1=COEFX*FXTK1
      FYTK1=COEFY*FYTK1
      FXTK1p=COEFX*FXTK1p
      FYTK1p=COEFY*FYTK1p
      FXT=c6*((2.*FXTK+FXTKm)*HZK1 +(2.*FXTK1+FXTK1p)*HZK)
      FYT=c6*((2.*FYTK+FYTKm)*HZK1 +(2.*FYTK1+FYTK1p)*HZK)

      IF( k .EQ. 1) THEN
*     -------------- Surface -----------------

      Sum1= c24*
     * (KT(1,n)*S1*(2.*Sclr2(i,j,1)+Sclr2(i+1,j,1)+Sclr2(i+1,j-1,1))
     + +KT(2,n)*S2*(2.*Sclr2(i,j,1)+Sclr2(i,j-1,1)+Sclr2(i+1,j-1,1))
     + +KT(3,n)*S3*(2.*Sclr2(i,j,1)+Sclr2(i,j-1,1)+Sclr2(i-1,j  ,1))
     + +KT(4,n)*S4*(2.*Sclr2(i,j,1)+Sclr2(i-1,j,1)+Sclr2(i-1,j+1,1))
     + +KT(5,n)*S5*(2.*Sclr2(i,j,1)+Sclr2(i-1,j+1,1)+Sclr2(i,j+1,1))
     + +KT(6,n)*S6*(2.*Sclr2(i,j,1)+Sclr2(i,j+1,1)+Sclr2(i+1,j  ,1)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*Sclr2(i,j,2)+Sclr2(i+1,j,2)+Sclr2(i+1,j-1,2))
     + +KT(2,n)*S2*(2.*Sclr2(i,j,2)+Sclr2(i,j-1,2)+Sclr2(i+1,j-1,2))
     + +KT(3,n)*S3*(2.*Sclr2(i,j,2)+Sclr2(i,j-1,2)+Sclr2(i-1,j  ,2))
     + +KT(4,n)*S4*(2.*Sclr2(i,j,2)+Sclr2(i-1,j,2)+Sclr2(i-1,j+1,2))
     + +KT(5,n)*S5*(2.*Sclr2(i,j,2)+Sclr2(i-1,j+1,2)+Sclr2(i,j+1,2))
     + +KT(6,n)*S6*(2.*Sclr2(i,j,2)+Sclr2(i,j+1,2)+Sclr2(i+1,j  ,2)))

      FZT=-0.5*((Sum1+Sum2)*W(i,j,2)-2.*Sum1*(W(i,j,1)-PME(i,j)/dt))
**      FZT=-0.5*((Sum1+Sum2)*W(i,j,2)-2.*Sum1*W(i,j,1))

      ELSE
*     --------------  Deep water -------------

      Sum1= c24*
     * (KT(1,n)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
     + +KT(6,n)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k)))

      Sum1m= c24*
     * (KT(1,n)*S1*(2.*Sclr2(i,j,km)+Sclr2(i+1,j,km)+Sclr2(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*Sclr2(i,j,km)+Sclr2(i-1,j,km)+Sclr2(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*Sclr2(i,j,km)+Sclr2(i-1,j+1,km)+Sclr2(i,j+1,km))
     + +KT(6,n)*S6*(2.*Sclr2(i,j,km)+Sclr2(i,j+1,km)+Sclr2(i+1,j  ,km))
     *                                                              )

      Sum2= c24*
     * (KT(1,np)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
     + +KT(2,np)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
     + +KT(3,np)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
     + +KT(4,np)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
     + +KT(5,np)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
     + +KT(6,np)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k)))

      Sum2p= c24*
     *(KT(1,np)*S1*(2.*Sclr2(i,j,kp)+Sclr2(i+1,j,kp)+Sclr2(i+1,j-1,kp))
     ++KT(2,np)*S2*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i+1,j-1,kp))
     ++KT(3,np)*S3*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i-1,j  ,kp))
     ++KT(4,np)*S4*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j,kp)+Sclr2(i-1,j+1,kp))
     ++KT(5,np)*S5*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j+1,kp)+Sclr2(i,j+1,kp))
     ++KT(6,np)*S6*(2.*Sclr2(i,j,kp)+Sclr2(i,j+1,kp)+Sclr2(i+1,j  ,kp))
     *                                                               )

      FZT=-0.5*( (Sum2+Sum2p)*W(i,j,Kp) -(Sum1+Sum1m)*W(i,j,K) )
      END IF

      FT=-FXT-FYT+FZT  +A*(DIFT+DIFTP)

*     ------------------------------------------------------------------
*               Gent & McWilliams transport, if necessary
*     ------------------------------------------------------------------

      call SkewFlux(Sclr2,Ropot,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &              il1,jl1,kl,klp,KT,z,Agm)
      DIFT = A*DIFT
      DIFTP= A*DIFTP

	FT = FT +DIFT +DIFTP

*     ------------------------------------------------------------------

      RsT(i,j,k)= dt*FT +dt*Qcascad(i,j,k)/(Rhx*Rhy)

*     First guess for mas matrix inversion
cc	Sclr1(i,j,k)= RsT(i,j,k)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))
	Sclr1(i,j,k)= 0. !!!RsT(i,j,k)/sclamp(i,j,k)


c	if(i.eq.19.and.j.eq.5.and.k.eq.1) write(*,*)'init',Sclr1(i,j,k)

      end do

*     --------------- Bottom ------------------
      N = NP
      call DIFFUW_TG(Sclr2,u,v,w,hz,i,j,kb,N,NP,DIFT,DIFTP,KB,
     &               il1,jl1,kl,dt,KT,ALT)
      call transp(FXTK  ,FYTK  ,i,j,kb  ,n,u,v,Sclr2,KT,il1,jl1,kl)
      call transp(FXTKm ,FYTKm ,i,j,kb-1,n,u,v,Sclr2,KT,il1,jl1,kl)

      FXTK=COEFX*FXTK
      FYTK=COEFY*FYTK
      FXTKm=COEFX*FXTKm
      FYTKm=COEFY*FYTKm
      FXT=c6*(2.*FXTK+FXTKm)*HZ(KB-1)
      FYT=c6*(2.*FYTK+FYTKm)*HZ(KB-1)

	k=kb
	km=kb-1

      Sum1= c24*
     * (KT(1,n)*S1*(2.*Sclr2(i,j,k)+Sclr2(i+1,j,k)+Sclr2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*Sclr2(i,j,k)+Sclr2(i,j-1,k)+Sclr2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*Sclr2(i,j,k)+Sclr2(i-1,j,k)+Sclr2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*Sclr2(i,j,k)+Sclr2(i-1,j+1,k)+Sclr2(i,j+1,k))
     + +KT(6,n)*S6*(2.*Sclr2(i,j,k)+Sclr2(i,j+1,k)+Sclr2(i+1,j  ,k)))

      Sum2= c24*
     *(KT(1,n)*S1*(2.*Sclr2(i,j,km)+Sclr2(i+1,j,km)+Sclr2(i+1,j-1,km))
     ++KT(2,n)*S2*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i+1,j-1,km))
     ++KT(3,n)*S3*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i-1,j  ,km))
     ++KT(4,n)*S4*(2.*Sclr2(i,j,km)+Sclr2(i-1,j,km)+Sclr2(i-1,j+1,km))
     ++KT(5,n)*S5*(2.*Sclr2(i,j,km)+Sclr2(i-1,j+1,km)+Sclr2(i,j+1,km))
     ++KT(6,n)*S6*(2.*Sclr2(i,j,km)+Sclr2(i,j+1,km)+Sclr2(i+1,j  ,km)))

      FZT=0.5*(Sum1+Sum2)*W(i,j,KB)

      FT=-FXT-FYT+FZT  +A*DIFT

*     ------------------------------------------------------------------
*               Gent & McWilliams transport, if necessary.
*          Some believe that it should't be at the boundaries
*     ------------------------------------------------------------------

      call SkewFlux(Sclr2,Ropot,hz,i,j,kb,N,NP,DIFT,DIFTP,KB,
     &              il1,jl1,kl,klp,KT,z,Agm)
      DIFT = A*DIFT
	FT = FT +DIFT
*     ------------------------------------------------------------------

      RsT(I,J,KB)= dt*FT +dt*Qcascad(i,j,kb)/(Rhx*Rhy)

*     First guess for mass matrix inversion
cc	Sclr1(i,j,kb)= RsT(i,j,kb)/( r2s12*CG(n)*HZ(KB-1))
	Sclr1(i,j,kb)= 0. !!!!RsT(i,j,kb)/sclamp(i,j,kb)

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

	sum= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp))

      Sclr(i,j,k)= Sclr1(i,j,k) + omega*(RsT(i,j,k)-sum)/sclamp(i,j,k)

cccc     *        + omega*(RsT(i,j,k)-sum)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

c	if(i.eq.12.and.j.eq.9.and.k.eq.23) 
c     *write(*,*)iter,Sclr(i,j,k),
cc     *(RsT(i,j,k)-sum)/( r2s12*(CG(n)*HZK1+CG(np)*HZK))
c     *(RsT(i,j,k)-sum)/sclamp(i,j,k)

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

	END DO   ! iter, Mass matrix iteration

	do k=1,kl
      do j=1,jl
	do i=1,il
	if(nt3(i,j,k) .LT. 0) then  ! Liquid boundaries
****	if(nt3(i,j,k) .NE. 0) then  ! All points
	Sclr(i,j,k)= Sclr1(i,j,k) +Sclr2(i,j,k)    
	end if
	end do
	end do
	end do

c	RSmean = 0.
	
c	do i=1,il
c	do j=1,jl
c	kb=km2(i,j)
c	if(kb .GT. 0 .and. nt3(i,j,1).gt.0) then
c	do k=1,kb
c	RSmean = RSmean + RsT(i,j,k)
c	end do
c	end if
c	end do
c	end do
	
c      write(*,*)'RSmean', Sclr(17,24,1), RSmean
	


***      GOTO 2000 !!!!!!!!!!!!!

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

	sum= c6*(HZK1*(2.*sumk+sumkm) +HZK*(2.*sumkp+sumkpp)) ! M*T

	Sclr(i,j,k)= (1.-gamma_fct)*Sclr2(i,j,k) 
     *            +(RsT(i,j,k) +gamma_fct*sum )/sclamp(i,j,k)

ccc     *+(RsT(i,j,k) +gamma_fct*sum )/( r2s12*(CG(n)*HZK1+CG(np)*HZK))

c	if(i.eq.3.and.j.eq.11.and.k.eq.19) write(*,*) Sclr(i,j,k), 
c     *Sclr2(i,j,k), RsT(i,j,k), 
c     *sum/( r2s12*(CG(n)*HZK1+CG(np)*HZK)), 
c     *cg(n),cg(np),hzk1,hzk,n,np,kb,sumk,sumkm,sumkp,sumkpp

	end do   ! k

      END IF   ! nt3 >0
	end do   ! i
	end do   ! j


****      call TSbc(Tdamp_In,Tdamp_Out,
****     *  dt,Sclr,Sclr2,SclrObs,u,v,w,nt3,Si,il1,jl1,kl,hx,hy,hz,R)

****      GOTO 2000

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
	
	do i=1,il

	kb=km2(i,j)
	IF( nt3(i,j,1) .GT. 0) then

      do k=1,kb
	n=nt3(i,j,k)

	if(k.GT.1)then
	km=k-1
	hzk1=hz(km)
	else
	hzk1=0.
	km=k
	end if

	if(k.LT.kb)then
	kp=k+1
	np=nt3(i,j,kp)
	hzk=hz(k)
	else
	kp=k
	hzk= 0.
	np=n
	end if

	QQQ=1./sclamp(i,j,k) ! Lamped mass matrix inverted

*     Antidiffusive fluxes - 12 in each node.
*     We need low-order solution and increment of the high-order one.


      if(k.GT.1) then

*     Element 1

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i+1,j,k )+Sclr2(i+1,j-1,k ))
	sumkm=c24*(2.*Sclr2(i,j,km)+Sclr2(i+1,j,km)+Sclr2(i+1,j-1,km))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkm) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i+1,j,k )+Sclr1(i+1,j-1,k ))
	sumkm=c24*(2.*Sclr1(i,j,km)+Sclr1(i+1,j,km)+Sclr1(i+1,j-1,km))

	Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkm)

      antidf(1,i,j,k)= KT(1,n)*HZK1*S1*(gamma_fct*Sum1+Sum2) *QQQ 


*     Element 2

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i,j-1,k )+Sclr2(i+1,j-1,k ))
	sumkm=c24*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i+1,j-1,km))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkm) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i,j-1,k )+Sclr1(i+1,j-1,k ))
	sumkm=c24*(2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i+1,j-1,km))

	Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkm)

      antidf(2,i,j,k)= KT(2,n)*S2*HZK1*(gamma_fct*Sum1+Sum2) *QQQ 


*     Element 3

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i,j-1,k )+Sclr2(i-1,j  ,k ))
	sumkm=c24*(2.*Sclr2(i,j,km)+Sclr2(i,j-1,km)+Sclr2(i-1,j  ,km))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkm) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i,j-1,k )+Sclr1(i-1,j  ,k ))
	sumkm=c24*(2.*Sclr1(i,j,km)+Sclr1(i,j-1,km)+Sclr1(i-1,j  ,km))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkm) 

      antidf(3,i,j,k)= KT(3,n)*S3*HZK1*(gamma_fct*Sum1+Sum2) *QQQ 

*     Element 4

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i-1,j,k )+Sclr2(i-1,j+1,k ))
	sumkm=c24*(2.*Sclr2(i,j,km)+Sclr2(i-1,j,km)+Sclr2(i-1,j+1,km))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkm) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i-1,j,k )+Sclr1(i-1,j+1,k ))
	sumkm=c24*(2.*Sclr1(i,j,km)+Sclr1(i-1,j,km)+Sclr1(i-1,j+1,km))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkm) 

      antidf(4,i,j,k)= KT(4,n)*S4*HZK1*(gamma_fct*Sum1+Sum2) *QQQ 

*     Element 5

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i-1,j+1,k )+Sclr2(i,j+1,k ))
	sumkm=c24*(2.*Sclr2(i,j,km)+Sclr2(i-1,j+1,km)+Sclr2(i,j+1,km))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkm) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i-1,j+1,k )+Sclr1(i,j+1,k ))
	sumkm=c24*(2.*Sclr1(i,j,km)+Sclr1(i-1,j+1,km)+Sclr1(i,j+1,km))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkm) 

      antidf(5,i,j,k)= KT(5,n)*S5*HZK1*(gamma_fct*Sum1+Sum2) *QQQ 

*     Element 6

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i,j+1,k )+Sclr2(i+1,j  ,k ))
	sumkm=c24*(2.*Sclr2(i,j,km)+Sclr2(i,j+1,km)+Sclr2(i+1,j  ,km))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkm) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i,j+1,k )+Sclr1(i+1,j  ,k ))
	sumkm=c24*(2.*Sclr1(i,j,km)+Sclr1(i,j+1,km)+Sclr1(i+1,j  ,km))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkm) 

      antidf(6,i,j,k)= KT(6,n)*S6*HZK1*(gamma_fct*Sum1+Sum2) *QQQ 

      end if

	if(k.LT.kb) then
*     Lower part

*     Element 7

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i+1,j,k )+Sclr2(i+1,j-1,k ))
	sumkp=c24*(2.*Sclr2(i,j,kp)+Sclr2(i+1,j,kp)+Sclr2(i+1,j-1,kp))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkp) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i+1,j,k )+Sclr1(i+1,j-1,k ))
	sumkp=c24*(2.*Sclr1(i,j,kp)+Sclr1(i+1,j,kp)+Sclr1(i+1,j-1,kp))

	Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkp)

      antidf(7,i,j,k)= KT(1,np)*S1*HZK*(gamma_fct*Sum1+Sum2) *QQQ 


*     Element 8

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i,j-1,k )+Sclr2(i+1,j-1,k ))
	sumkp=c24*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i+1,j-1,kp))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkp) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i,j-1,k )+Sclr1(i+1,j-1,k ))
	sumkp=c24*(2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i+1,j-1,kp))

	Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkp)

      antidf(8,i,j,k)= KT(2,np)*S2*HZK*(gamma_fct*Sum1+Sum2) *QQQ 


*     Element 9

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i,j-1,k )+Sclr2(i-1,j  ,k ))
	sumkp=c24*(2.*Sclr2(i,j,kp)+Sclr2(i,j-1,kp)+Sclr2(i-1,j  ,kp))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkp) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i,j-1,k )+Sclr1(i-1,j  ,k ))
	sumkp=c24*(2.*Sclr1(i,j,kp)+Sclr1(i,j-1,kp)+Sclr1(i-1,j  ,kp))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkp) 

      antidf(9,i,j,k)= KT(3,np)*S3*HZK*(gamma_fct*Sum1+Sum2) *QQQ 

*     Element 10

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i-1,j,k )+Sclr2(i-1,j+1,k ))
	sumkp=c24*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j,kp)+Sclr2(i-1,j+1,kp))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkp) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i-1,j,k )+Sclr1(i-1,j+1,k ))
	sumkp=c24*(2.*Sclr1(i,j,kp)+Sclr1(i-1,j,kp)+Sclr1(i-1,j+1,kp))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkp) 

      antidf(10,i,j,k)= KT(4,np)*S4*HZK*(gamma_fct*Sum1+Sum2) *QQQ 

*     Element 11

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i-1,j+1,k )+Sclr2(i,j+1,k ))
	sumkp=c24*(2.*Sclr2(i,j,kp)+Sclr2(i-1,j+1,kp)+Sclr2(i,j+1,kp))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkp) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i-1,j+1,k )+Sclr1(i,j+1,k ))
	sumkp=c24*(2.*Sclr1(i,j,kp)+Sclr1(i-1,j+1,kp)+Sclr1(i,j+1,kp))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkp) 

      antidf(11,i,j,k)= KT(5,np)*S5*HZK*(gamma_fct*Sum1+Sum2) *QQQ 

*     Element 12

	sumk =c24*(2.*Sclr2(i,j,k )+Sclr2(i,j+1,k )+Sclr2(i+1,j  ,k ))
	sumkp=c24*(2.*Sclr2(i,j,kp)+Sclr2(i,j+1,kp)+Sclr2(i+1,j  ,kp))

      Sum1= c12*Sclr2(i,j,k) - c6*(2.*sumk+sumkp) 

	sumk =c24*(2.*Sclr1(i,j,k )+Sclr1(i,j+1,k )+Sclr1(i+1,j  ,k ))
	sumkp=c24*(2.*Sclr1(i,j,kp)+Sclr1(i,j+1,kp)+Sclr1(i+1,j  ,kp))

      Sum2= c12*Sclr1(i,j,k) - c6*(2.*sumk+sumkp) 

      antidf(12,i,j,k)= KT(6,np)*S6*HZK*(gamma_fct*Sum1+Sum2) *QQQ 

	end if



c	IF(i.eq.25.and.j.eq.30.and.k.eq.1) then
c	write(*,*) antidf(:,i,j,k)
c	write(*,*) Sclr2(i,j,k), Sclr(i,j,k), Sclr2(i,j,k)+Sclr1(i,j,k)
c      end if



	end do   ! k

      END IF   ! nt3 ne 0
	end do
	end do

*     Sum of negative and positive fluxes in each node

	do j=1,jl
	do i=1,il
	kb=km2(i,j)

	IF( kb .GT. 0) then

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

      if(k.GT.1) then ! upper half

	km= k-1

      if(KT(1,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,km),Sclr(i+1,j,km),Sclr(i+1,j-1,km))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i+1,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,km),Sclr(i+1,j,km),Sclr(i+1,j-1,km))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i+1,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,km),Sclr2(i+1,j,km),Sclr2(i+1,j-1,km))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i+1,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,km),Sclr2(i+1,j,km),Sclr2(i+1,j-1,km))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i+1,j-1,k ))
	end if

      if(KT(2,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,km),Sclr(i,j-1,km),Sclr(i+1,j-1,km))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i,j-1,k ),Sclr(i+1,j-1,k ))
      Tmin= min(Tmin, Sclr(i,j,km),Sclr(i,j-1,km),Sclr(i+1,j-1,km))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i,j-1,k ),Sclr(i+1,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,km),Sclr2(i,j-1,km),Sclr2(i+1,j-1,km))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i,j-1,k ),Sclr2(i+1,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,km),Sclr2(i,j-1,km),Sclr2(i+1,j-1,km))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i,j-1,k ),Sclr2(i+1,j-1,k ))
	end if

      if(KT(3,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,km),Sclr(i-1,j,km),Sclr(i,j-1,km))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,km),Sclr(i-1,j,km),Sclr(i,j-1,km))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,km),Sclr2(i-1,j,km),Sclr2(i,j-1,km))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,km),Sclr2(i-1,j,km),Sclr2(i,j-1,km))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i,j-1,k ))
	end if

      if(KT(4,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,km),Sclr(i-1,j,km),Sclr(i-1,j+1,km))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i-1,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,km),Sclr(i-1,j,km),Sclr(i-1,j+1,km))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i-1,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,km),Sclr2(i-1,j,km),Sclr2(i-1,j+1,km))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i-1,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,km),Sclr2(i-1,j,km),Sclr2(i-1,j+1,km))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i-1,j+1,k ))
	end if

      if(KT(5,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,km),Sclr(i,j+1,km),Sclr(i-1,j+1,km))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i,j+1,k ),Sclr(i-1,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,km),Sclr(i,j+1,km),Sclr(i-1,j+1,km))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i,j+1,k ),Sclr(i-1,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,km),Sclr2(i,j+1,km),Sclr2(i-1,j+1,km))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i,j+1,k ),Sclr2(i-1,j+1,k ))
c      Tmin= min(Tmin, Sclr2(i,j,km),Sclr2(i,j+1,km),Sclr2(i-1,j+1,km))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i,j+1,k ),Sclr2(i-1,j+1,k ))
	end if

      if(KT(6,n) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,km),Sclr(i+1,j,km),Sclr(i,j+1,km))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,km),Sclr(i+1,j,km),Sclr(i,j+1,km))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,km),Sclr2(i+1,j,km),Sclr2(i,j+1,km))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,km),Sclr2(i+1,j,km),Sclr2(i,j+1,km))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i,j+1,k ))
	end if
	
	Tmax= max(Tmax, Sclr2(i,j,km))
      Tmax= max(Tmax, Sclr2(i,j,k ))
	Tmin= min(Tmin, Sclr2(i,j,km))
      Tmin= min(Tmin, Sclr2(i,j,k ))  	
	
	end if
	
	if(k.LT.kb) then ! lower half

      kp= k+1
      np=abs(nt3(i,j,kp))

      if(KT(1,np) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,kp),Sclr(i+1,j,kp),Sclr(i+1,j-1,kp))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i+1,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,kp),Sclr(i+1,j,kp),Sclr(i+1,j-1,kp))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i+1,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,kp),Sclr2(i+1,j,kp),Sclr2(i+1,j-1,kp))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i+1,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,kp),Sclr2(i+1,j,kp),Sclr2(i+1,j-1,kp))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i+1,j-1,k ))
	end if

      if(KT(2,np) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,kp),Sclr(i,j-1,kp),Sclr(i+1,j-1,kp))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i,j-1,k ),Sclr(i+1,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,kp),Sclr(i,j-1,kp),Sclr(i+1,j-1,kp))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i,j-1,k ),Sclr(i+1,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,kp),Sclr2(i,j-1,kp),Sclr2(i+1,j-1,kp))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i,j-1,k ),Sclr2(i+1,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,kp),Sclr2(i,j-1,kp),Sclr2(i+1,j-1,kp))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i,j-1,k ),Sclr2(i+1,j-1,k ))
	end if

      if(KT(3,np) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,kp),Sclr(i-1,j,kp),Sclr(i,j-1,kp))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i,j-1,k ))
	Tmin= min(Tmin, Sclr(i,j,kp),Sclr(i-1,j,kp),Sclr(i,j-1,kp))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i,j-1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,kp),Sclr2(i-1,j,kp),Sclr2(i,j-1,kp))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i,j-1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,kp),Sclr2(i-1,j,kp),Sclr2(i,j-1,kp))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i,j-1,k ))
	end if

      if(KT(4,np) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,kp),Sclr(i-1,j,kp),Sclr(i-1,j+1,kp))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i-1,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,kp),Sclr(i-1,j,kp),Sclr(i-1,j+1,kp))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i-1,j,k ),Sclr(i-1,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,kp),Sclr2(i-1,j,kp),Sclr2(i-1,j+1,kp))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i-1,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,kp),Sclr2(i-1,j,kp),Sclr2(i-1,j+1,kp))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i-1,j,k ),Sclr2(i-1,j+1,k ))
	end if

      if(KT(5,np) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,kp),Sclr(i,j+1,kp),Sclr(i-1,j+1,kp))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i,j+1,k ),Sclr(i-1,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,kp),Sclr(i,j+1,kp),Sclr(i-1,j+1,kp))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i,j+1,k ),Sclr(i-1,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,kp),Sclr2(i,j+1,kp),Sclr2(i-1,j+1,kp))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i,j+1,k ),Sclr2(i-1,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,kp),Sclr2(i,j+1,kp),Sclr2(i-1,j+1,kp))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i,j+1,k ),Sclr2(i-1,j+1,k ))
	end if

      if(KT(6,np) .GT. 0.) then	
	Tmax= max(Tmax, Sclr(i,j,kp),Sclr(i+1,j,kp),Sclr(i,j+1,kp))
	Tmax= max(Tmax, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i,j+1,k ))
	Tmin= min(Tmin, Sclr(i,j,kp),Sclr(i+1,j,kp),Sclr(i,j+1,kp))
	Tmin= min(Tmin, Sclr(i,j,k ),Sclr(i+1,j,k ),Sclr(i,j+1,k ))
*     Additional limitation like Zalesak
c	Tmax= max(Tmax, Sclr2(i,j,kp),Sclr2(i+1,j,kp),Sclr2(i,j+1,kp))
c	Tmax= max(Tmax, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i,j+1,k ))
c	Tmin= min(Tmin, Sclr2(i,j,kp),Sclr2(i+1,j,kp),Sclr2(i,j+1,kp))
c	Tmin= min(Tmin, Sclr2(i,j,k ),Sclr2(i+1,j,k ),Sclr2(i,j+1,k ))
	end if
	
	Tmax= max(Tmax, Sclr2(i,j,kp))
      Tmax= max(Tmax, Sclr2(i,j,k ))
	Tmin= min(Tmin, Sclr2(i,j,kp))
      Tmin= min(Tmin, Sclr2(i,j,k )) 
      
           
	end if
	      
*     Admitted increment in the node

      Qplus = Tmax -Sclr(i,j,k)
	Qminus= Tmin -Sclr(i,j,k)

*     Antidiffusive flux input in the node 

      if(Pplus. GT. 0.) then
      Rplus(i,j,k) = min(1., Qplus/max(1.e-9,Pplus))
	else
	Rplus(i,j,k) = 0.
	end if

      if(Pminus. LT. 0.) then
      Rminus(i,j,k) = min(1., Qminus/min(-1.e-9,Pminus))
	else
	Rminus(i,j,k) = 0.
	end if

c	IF(i.eq.25.and.j.eq.30.and.k.eq.1) then
c	write(*,*) antidf(:,i,j,k)
c	write(*,*) Rplus(i,j,k), Qplus,  Pplus, Tmax, Sclr(i,j,k)
c	write(*,*) Rminus(i,j,k),Qminus, Pminus,Tmin, Sclr(i,j,k)
c      end if



	end do   ! k

      END IF   ! kb > 0
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


	n=nt3(i,j,k)

	if(k.LT.kb) then
	np=nt3(i,j,k+1)
	else
      np=n
	end if

      alpha_e =0.0

*     Weights alpha


*     Element #1

      if( KT(1,n) .GT. 0. .and. k.GT.1 ) then

      if(antidf(1,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(3,i+1,j,k) .GE. 0.) then
      a2= Rplus (i+1,j,k)
	else
      a2= Rminus(i+1,j,k)
      end if
      if(antidf(5,i+1,j-1,k) .GE. 0.) then
      a3= Rplus (i+1,j-1,k)
	else
      a3= Rminus(i+1,j-1,k)
      end if

      km=k-1

      if(antidf(7,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(9,i+1,j,km) .GE. 0.) then
      a2= min(a2,Rplus (i+1,j,km))
	else
      a2= min(a2,Rminus(i+1,j,km))
      end if
      if(antidf(11,i+1,j-1,km) .GE. 0.) then
      a3= min(a3,Rplus (i+1,j-1,km))
	else
      a3= min(a3,Rminus(i+1,j-1,km))
      end if

	alpha_e(1)= min(a1,a2,a3)

	end if ! #1



*     Element #2


      if( KT(2,n) .GT. 0. .and. k.GT.1 ) then

      if(antidf(2,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(4,i+1,j-1,k) .GE. 0.) then
      a2= Rplus (i+1,j-1,k)
	else
      a2= Rminus(i+1,j-1,k)
      end if
      if(antidf(6,i,j-1,k) .GE. 0.) then
      a3= Rplus (i,j-1,k)
	else
      a3= Rminus(i,j-1,k)
      end if

      km= k-1

      if(antidf(8,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(10,i+1,j-1,km) .GE. 0.) then
      a2= min(a2,Rplus (i+1,j-1,km))
	else
      a2= min(a2,Rminus(i+1,j-1,km))
      end if
      if(antidf(12,i,j-1,km) .GE. 0.) then
      a3= min(a3,Rplus (i,j-1,km))
	else
      a3= min(a3,Rminus(i,j-1,km))
      end if

	alpha_e(2)= min(a1,a2,a3)

	end if ! #2



*     Element #3

      if( KT(3,n) .GT. 0. .and. k.GT.1 ) then

      if(antidf(3,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(5,i,j-1,k) .GE. 0.) then
      a2= Rplus (i,j-1,k)
	else
      a2= Rminus(i,j-1,k)
      end if
      if(antidf(1,i-1,j,k) .GE. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      km= k-1

      if(antidf(9,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(11,i,j-1,km) .GE. 0.) then
      a2= min(a2,Rplus (i,j-1,km))
	else
      a2= min(a2,Rminus(i,j-1,km))
      end if
      if(antidf(7,i-1,j,km) .GE. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if

	alpha_e(3)= min(a1,a2,a3)

	end if ! #3



*     Element #4


      if( KT(4,n) .GT. 0. .and. k.GT.1 ) then

      if(antidf(4,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(2,i-1,j+1,k) .GE. 0.) then
      a2= Rplus (i-1,j+1,k)
	else
      a2= Rminus(i-1,j+1,k)
      end if
      if(antidf(6,i-1,j,k) .GE. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      km= k-1

      if(antidf(10,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(8,i-1,j+1,km) .GE. 0.) then
      a2= min(a2,Rplus (i-1,j+1,km))
	else
      a2= min(a2,Rminus(i-1,j+1,km))
      end if
      if(antidf(12,i-1,j,km) .GE. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if

	alpha_e(4)= min(a1,a2,a3)

	end if ! #4


*     Element #5

      if( KT(5,n) .GT. 0. .and. k.GT.1 ) then

      if(antidf(5,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(3,i,j+1,k) .GE. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(1,i-1,j+1,k) .GE. 0.) then
      a3= Rplus (i-1,j+1,k)
	else
      a3= Rminus(i-1,j+1,k)
      end if

      km= k-1

      if(antidf(11,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(9,i,j+1,km) .GE. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(7,i-1,j+1,km) .GE. 0.) then
      a3= min(a3,Rplus (i-1,j+1,km))
	else
      a3= min(a3,Rminus(i-1,j+1,km))
      end if

	alpha_e(5)= min(a1,a2,a3)

	end if ! #5


*     Element #6
 
      if( KT(6,n) .GT. 0. .and. k.GT.1 ) then

      if(antidf(6,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(2,i,j+1,k) .GE. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(4,i+1,j,k) .GE. 0.) then
      a3= Rplus (i+1,j,k)
	else
      a3= Rminus(i+1,j,k)
      end if

      km= k-1

      if(antidf(12,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(8,i,j+1,km) .GE. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(10,i+1,j,km) .GE. 0.) then
      a3= min(a3,Rplus (i+1,j,km))
	else
      a3= min(a3,Rminus(i+1,j,km))
      end if

       alpha_e(6)= min(a1,a2,a3)

	end if ! #6


*     Element #7

      if( KT(1,np) .GT. 0. .and. k.LT.kb) then

      if(antidf(7,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(9,i+1,j,k) .GE. 0.) then
      a2= Rplus (i+1,j,k)
	else
      a2= Rminus(i+1,j,k)
      end if
      if(antidf(11,i+1,j-1,k) .GE. 0.) then
      a3= Rplus (i+1,j-1,k)
	else
      a3= Rminus(i+1,j-1,k)
      end if

      km= k+1

      if(antidf(1,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(3,i+1,j,km) .GE. 0.) then
      a2= min(a2,Rplus (i+1,j,km))
	else
      a2= min(a2,Rminus(i+1,j,km))
      end if
      if(antidf(5,i+1,j-1,km) .GE. 0.) then
      a3= min(a3,Rplus (i+1,j-1,km))
	else
      a3= min(a3,Rminus(i+1,j-1,km))
      end if

	alpha_e(7)= min(a1,a2,a3)

	end if ! #7


*     Element #8

      if( KT(2,np) .GT. 0. .and. k.LT.kb) then

      if(antidf(8,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(10,i+1,j-1,k) .GE. 0.) then
      a2= Rplus (i+1,j-1,k)
	else
      a2= Rminus(i+1,j-1,k)
      end if
      if(antidf(12,i,j-1,k) .GE. 0.) then
      a3= Rplus (i,j-1,k)
	else
      a3= Rminus(i,j-1,k)
      end if

      km= k+1

      if(antidf(2,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(4,i+1,j-1,km) .GE. 0.) then
      a2= min(a2,Rplus (i+1,j-1,km))
	else
      a2= min(a2,Rminus(i+1,j-1,km))
      end if
      if(antidf(6,i,j-1,km) .GE. 0.) then
      a3= min(a3,Rplus (i,j-1,km))
	else
      a3= min(a3,Rminus(i,j-1,km))
      end if

	alpha_e(8)= min(a1,a2,a3)

	end if ! #8


*     Element #9

      if( KT(3,np) .GT. 0. .and. k.LT.kb) then

      if(antidf(9,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(11,i,j-1,k) .GE. 0.) then
      a2= Rplus (i,j-1,k)
	else
      a2= Rminus(i,j-1,k)
      end if
      if(antidf(7,i-1,j,k) .GE. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      km= k+1

      if(antidf(3,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(5,i,j-1,km) .GE. 0.) then
      a2= min(a2,Rplus (i,j-1,km))
	else
      a2= min(a2,Rminus(i,j-1,km))
      end if
      if(antidf(1,i-1,j,km) .GE. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if

	alpha_e(9)= min(a1,a2,a3)

	end if ! #9


*     Element #10

      if( KT(4,np) .GT. 0. .and. k.LT.kb) then

      if(antidf(10,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(8,i-1,j+1,k) .GE. 0.) then
      a2= Rplus (i-1,j+1,k)
	else
      a2= Rminus(i-1,j+1,k)
      end if
      if(antidf(12,i-1,j,k) .GE. 0.) then
      a3= Rplus (i-1,j,k)
	else
      a3= Rminus(i-1,j,k)
      end if

      km= k+1

      if(antidf(4,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(2,i-1,j+1,km) .GE. 0.) then
      a2= min(a2,Rplus (i-1,j+1,km))
	else
      a2= min(a2,Rminus(i-1,j+1,km))
      end if
      if(antidf(6,i-1,j,km) .GE. 0.) then
      a3= min(a3,Rplus (i-1,j,km))
	else
      a3= min(a3,Rminus(i-1,j,km))
      end if

	alpha_e(10)= min(a1,a2,a3)

	end if ! #10


*     Element #11

      if( KT(5,np) .GT. 0. .and. k.LT.kb) then

      if(antidf(11,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(9,i,j+1,k) .GE. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(7,i-1,j+1,k) .GE. 0.) then
      a3= Rplus (i-1,j+1,k)
	else
      a3= Rminus(i-1,j+1,k)
      end if

      km= k+1

      if(antidf(5,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(3,i,j+1,km) .GE. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(1,i-1,j+1,km) .GE. 0.) then
      a3= min(a3,Rplus (i-1,j+1,km))
	else
      a3= min(a3,Rminus(i-1,j+1,km))
      end if

	alpha_e(11)= min(a1,a2,a3)

	end if ! #11


*     Element #12

      if( KT(6,np) .GT. 0. .and. k.LT.kb) then

      if(antidf(12,i,j,k) .GE. 0.) then
      a1= Rplus (i,j,k)
	else
      a1= Rminus(i,j,k)
      end if
      if(antidf(8,i,j+1,k) .GE. 0.) then
      a2= Rplus (i,j+1,k)
	else
      a2= Rminus(i,j+1,k)
      end if
      if(antidf(10,i+1,j,k) .GE. 0.) then
      a3= Rplus (i+1,j,k)
	else
      a3= Rminus(i+1,j,k)
      end if

      km= k+1

      if(antidf(6,i,j,km) .GE. 0.) then
      a1= min(a1,Rplus (i,j,km))
	else
      a1= min(a1,Rminus(i,j,km))
      end if
      if(antidf(2,i,j+1,km) .GE. 0.) then
      a2= min(a2,Rplus (i,j+1,km))
	else
      a2= min(a2,Rminus(i,j+1,km))
      end if
      if(antidf(4,i+1,j,km) .GE. 0.) then
      a3= min(a3,Rplus (i+1,j,km))
	else
      a3= min(a3,Rminus(i+1,j,km))
      end if

	alpha_e(12)= min(a1,a2,a3)

	end if ! #12


c	IF(i.eq.24.and.j.eq.29.and.k.eq.1) then
c	write(*,*) alpha_e(:)
c	end if

*----------------------------------------------------------------
      ADD = 0.
      do nelem=1,12
      ADD= ADD + alpha_e(nelem)*antidf(nelem,i,j,k)
	end do
*----------------------------------------------------------------

      Sclr(i,j,k)= Sclr(i,j,k) +ADD

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

*     Version 29.06.2015.

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



	subroutine TSbc_inc(Tdamp_In,Tdamp_Out,
     *           dt,T,Tm2,Tobs,u,v,w, nt3,Si,il,jl,kl,hx,hy,hz,R)

*     T is the time increment of the scalar 

*     Version 14/08/2015

      dimension T(0:il,0:jl,kl), Tm2(0:il,0:jl,kl), hz(kl),
     &          Tobs(0:il,0:jl,kl), nt3(0:il,0:jl,kl),Si(0:jl)
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)

	puny= 1.e-19

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
     *        -T(i+2,j,k)-Tm2(i+2,j,k))/(hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))/hy
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))/hy
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/hy

c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i,j+1,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i,j+1,k))/4.

c	Rx= -dt*umean/(R*hx*S0)
c	Ry= +dt*vmean/(R*hy)

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)
	
	Rx= Rx -dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i,j,k)/(R*hy)

c	Rx= -dt*u(i,j,k)/(R*hx*S0)
c	Ry= +dt*v(i,j,k)/(R*hy)

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


*     Upwind scheme
c      un= (2.*v(i,j,k) +v(i-1,j,k) +v(i+1,j,k))*0.25
c	T(i,j,k)=dt*((un+abs(un))*(Tm2 (i,j,k)-Tm2(i,j-1,k)) -
c     +             (un-abs(un))*(Tobs(i,j,k)-Tm2(i,j  ,k)) )/(R*hy)




      delta_t= T(i,j-1,k)
	delta_x=(T(i,j-1,k)+Tm2(i,j-1,k)
     *        -T(i,j-2,k)-Tm2(i,j-2,k))/hy
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))/(hx*S0)
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(hx*S0)
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx +dt*v(i,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i,j-1,k)/(R*hx*S0)
	Rx= Rx +dt*v(i,j,k)/(R*hy)
	Ry= Ry +dt*u(i,j,k)/(R*hx*S0)


c	Rx= +dt*v(i,j,k)/(R*hy)
c	Ry= +dt*u(i,j,k)/(R*hx*S0)


c	umean= (2.*u(i,j,k)+u(i-1,j,k)+u(i+1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i-1,j,k)+v(i+1,j,k))/4.

c	Ry= Rx +dt*umean/(R*hx*S0)
c	Rx= Ry +dt*vmean/(R*hy)

c	Ry= +dt*umean/(R*hx*S0)
c	Rx= +dt*vmean/(R*hy)


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
     *        -T(i-2,j,k)-Tm2(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k)) /hy
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /hy
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/hy 

c	Rx= Rx +dt*u(i-1,j,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j,k)/(R*hy)
	Rx= Rx +dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i,j,k)/(R*hy)

c	Rx= +dt*u(i,j,k)/(R*hx*S0)
c	Ry= -dt*v(i,j,k)/(R*hy)

c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i,j+1,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i,j+1,k))/4.

c	Rx= Rx +dt*umean/(R*hx*S0)
c	Ry= Ry -dt*vmean/(R*hy)

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
     *        -T(i,j+2,k)-Tm2(i,j+2,k)) /hy
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(hx*S0)
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(hx*S0)
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx -dt*v(i,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i,j+1,k)/(R*hx*S0)
	Rx= Rx -dt*v(i,j,k)/(R*hy)
	Ry= Ry -dt*u(i,j,k)/(R*hx*S0)

c	Rx= -dt*v(i,j,k)/(R*hy)
c	Ry= -dt*u(i,j,k)/(R*hx*S0)

c	umean= (2.*u(i,j,k)+u(i-1,j,k)+u(i+1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i-1,j,k)+v(i+1,j,k))/4.

c	Ry= -dt*umean/(R*hx*S0)
c	Rx= -dt*vmean/(R*hy)


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
     *        -T(i+2,j,k)-Tm2(i+2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i+1,j,k)-Tm2(i+1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i+1,j,k)-Tm2(i+1,j-1,k)) /hy
	else
	delta_y=0.
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /hy      

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/3.
c	vmean= (v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/3.

c      Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Rx +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j-1,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j-1,k)/(R*hy)
	Rx= Rx -dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i,j,k)/(R*hy)

c	Rx= -dt*u(i,j,k)/(R*hx*S0)
c	Ry= +dt*v(i,j,k)/(R*hy)


c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i+1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i+1,j,k))/4.

c	Ry=  dt*vmean/(R*hy)
c	Rx= -dt*umean/(R*hx*S0)

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
     *        -T(i,j-2,k)-Tm2(i,j-2,k)) /hy 
	sss= delta_t*(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j-1,k)-Tm2(i,j-1,k))/(hx*S0)
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx +dt*v(i+1,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i+1,j-1,k)/(R*hx*S0)
	Rx= Rx +dt*v(i,j,k)/(R*hy)
	Ry= Ry +dt*u(i,j,k)/(R*hx*S0)


c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i+1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i+1,j,k))/4.

c	Rx= +dt*vmean/(R*hy)
c	Ry= +dt*umean/(R*hx*S0)

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
     *        -T(i-2,j,k)-Tm2(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i-1,j,k)-Tm2(i-1,j-1,k))/hy
	else
	delta_y=0.
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/hy

c	Rx= Rx +dt*u(i-1,j-1,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j-1,k)/(R*hy)
	Rx= Rx +dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i,j,k)/(R*hy)

c	Rx= +dt*u(i,j,k)/(R*hx*S0)
c	Ry= -dt*v(i,j,k)/(R*hy)

c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i-1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i-1,j,k))/4.

c	Ry= -dt*vmean/(R*hy)
c	Rx= +dt*umean/(R*hx*S0)


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
     *        -T(i,j-2,k)-Tm2(i,j-2,k)) /hy 
	sss= delta_t*(Tm2(i,j-1,k)-Tm2(i-1,j-1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j-1,k)-Tm2(i-1,j-1,k)) /(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/hy
	Ry=-delta_t*delta_y*D/(hx*S0)

c	Rx= Rx +dt*v(i-1,j-1,k)/(R*hy)
c	Ry= Ry +dt*u(i-1,j-1,k)/(R*hx*S0)
	Rx= Rx +dt*v(i,j,k)/(R*hy)
	Ry= Ry +dt*u(i,j,k)/(R*hx*S0)

c	Rx= +dt*v(i,j,k)/(R*hy)
c	Ry= +dt*u(i,j,k)/(R*hx*S0)

c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i-1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i-1,j,k))/4.

c	Rx=  dt*vmean/(R*hy)
c	Ry= +dt*umean/(R*hx*S0)


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
     *        -T(i-2,j,k)-Tm2(i-2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i-1,j+1,k)-Tm2(i-1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i-1,j+1,k)-Tm2(i-1,j,k)) /hy
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D/(hx*S0)
	Ry=-delta_t*delta_y*D/hy

c	Rx= Rx +dt*u(i-1,j+1,k)/(R*hx*S0)
c	Ry= Ry -dt*v(i-1,j+1,k)/(R*hy)
	Rx= Rx +dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry -dt*v(i,j,k)/(R*hy)

c	Rx= +dt*u(i,j,k)/(R*hx*S0)
c	Ry= -dt*v(i,j,k)/(R*hy)

c	umean= (2.*u(i,j,k)+u(i,j+1,k)+u(i-1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j+1,k)+v(i-1,j,k))/4.

c	Ry= -dt*vmean/(R*hy)
c	Rx= +dt*umean/(R*hx*S0)


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
     *        -T(i,j+2,k)-Tm2(i,j+2,k)) /hy 
	sss= delta_t*(Tm2(i,j+1,k)-Tm2(i-1,j+1,k))

	if(sss.GT.0.)then
	delta_y=(Tm2(i,j+1,k)-Tm2(i-1,j+1,k)) /(hx*S0)
	else
	delta_y=0.
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /hy
	Ry=-delta_t*delta_y*D /(hx*S0)

c	Rx= Rx -dt*v(i-1,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i-1,j+1,k)/(R*hx*S0)
	Rx= Rx -dt*v(i,j,k)/(R*hy)
	Ry= Ry -dt*u(i,j,k)/(R*hx*S0)
	
c	Rx= -dt*v(i,j,k)/(R*hy)
c	Ry= -dt*u(i,j,k)/(R*hx*S0)

c	umean= (2.*u(i,j,k)+u(i,j+1,k)+u(i-1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j+1,k)+v(i-1,j,k))/4.

c	Rx= -dt*vmean/(R*hy)
c	Ry= -dt*umean/(R*hx*S0)


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
     *        -T(i+2,j,k)-Tm2(i+2,j,k)) /(hx*S0)
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i+1,j,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i+1,j,k)) /hy
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /hy 

c	umean= (u(i,j,k)+u(i+1,j,k)+u(i,j+1,k))/3.
c	vmean= (v(i,j,k)+v(i+1,j,k)+v(i,j+1,k))/3.

c	Rx= Rx -dt*umean/(R*hx*S0)
c	Ry= Ry +dt*vmean/(R*hy)

c	Rx= Rx -dt*u(i+1,j+1,k)/(R*hx*S0)
c	Ry= Ry +dt*v(i+1,j+1,k)/(R*hy)
	Rx= Rx -dt*u(i,j,k)/(R*hx*S0)
	Ry= Ry +dt*v(i,j,k)/(R*hy)

c	Rx= -dt*u(i,j,k)/(R*hx*S0)
c	Ry= +dt*v(i,j,k)/(R*hy)

c	umean= (2.*u(i,j,k)+u(i,j+1,k)+u(i+1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j+1,k)+v(i+1,j,k))/4.

c	Ry= +dt*vmean/(R*hy)
c	Rx= -dt*umean/(R*hx*S0)


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
     *        -T(i,j+2,k)-Tm2(i,j+2,k)) /hy 
	sss= delta_t*(Tm2(i+1,j+1,k)-Tm2(i,j+1,k))

	if(sss.GT.0.)then
	delta_y=0.
	else
	delta_y=(Tm2(i+1,j+1,k)-Tm2(i,j+1,k)) /(hx*S0)
	end if

      D= 1./MAX(puny,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /hy
	Ry=-delta_t*delta_y*D /(hx*S0)

c	Rx= Rx -dt*v(i+1,j+1,k)/(R*hy)
c	Ry= Ry -dt*u(i+1,j+1,k)/(R*hx*S0)
	Rx= Rx -dt*v(i,j,k)/(R*hy)
	Ry= Ry -dt*u(i,j,k)/(R*hx*S0)

c	Rx= -dt*v(i,j,k)/(R*hy)
c	Ry= -dt*u(i,j,k)/(R*hx*S0)

c	umean= (2.*u(i,j,k)+u(i,j-1,k)+u(i-1,j,k))/4.
c	vmean= (2.*v(i,j,k)+v(i,j-1,k)+v(i-1,j,k))/4.

c	Rx= -dt*vmean/(R*hy)
c	Ry= -dt*umean/(R*hx*S0)


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

c      if(i.eq.32.and.j.eq.32.and.k.eq.13) then
c	write(*,*) 'Transport',np,delta_t,delta_x,delta_y,Rx,Ry,D, 
c     &  T(i,j,k), Tm2(i,j,k), Tobs(i,j,k), nt3(i+1,j,k),nt3(i,j-1,k)
c	end if


	end do
	end do
	end do

	Return
	End


      SUBROUTINE DIFFUW_TG(T,u,v,w,hz,i,j,k,N,NP,DIFT,DIFTP,KB,
     &                     il,jl,kl,dt,KT,A)

*     Version 27.07.2015.

*     Galerkin artificial diffusion approximation by Taylor-Galerkin
*     A is the "physical" or rather "background" diffusivity.
*     Boundary conditions are taken into account in the Main Program.

      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), w(0:il,0:jl,kl)
      dimension T(0:il,0:jl,kl), KT(6,13), 
     *          hz(kl)
	real a11(6),a12(6),a21(6),a22(6)
	real a13(6),a23(6),a33(6)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
	common /consts/ c3,c6,c12,c18,Rhx,Rhy,R2hxhy
      real KT

	TAU  = 0.5*dt

	km=k-1
	kp=k+1	
      DIFT=0.
      DIFTP=0.
      Q1K=c3*(T(i,j,k)+T(i+1,j,k)+T(i+1,j-1,k))
      Q2K=c3*(T(i,j,k)+T(i+1,j-1,k)+T(i,j-1,k))
      Q3K=c3*(T(i,j,k)+T(i-1,j,k)+T(i,j-1,k))
      Q4K=c3*(T(i,j,k)+T(i-1,j,k)+T(i-1,j+1,k))
      Q5K=c3*(T(i,j,k)+T(i,j+1,k)+T(i-1,j+1,k))
      Q6K=c3*(T(i,j,k)+T(i+1,j,k)+T(i,j+1,k))

      W1K=0.25*(2.*T(i,j,k)+T(i+1,j,k)+T(i+1,j-1,k))
      W2K=0.25*(2.*T(i,j,k)+T(i+1,j-1,k)+T(i,j-1,k))
      W3K=0.25*(2.*T(i,j,k)+T(i-1,j,k)+T(i,j-1,k))
      W4K=0.25*(2.*T(i,j,k)+T(i-1,j,k)+T(i-1,j+1,k))
      W5K=0.25*(2.*T(i,j,k)+T(i,j+1,k)+T(i-1,j+1,k))
      W6K=0.25*(2.*T(i,j,k)+T(i+1,j,k)+T(i,j+1,k))

      IF( K .GT. 1) THEN
C     UPPER HALF

      Q1Km=c3*(T(i,j,km)+T(i+1,j,km)+T(i+1,j-1,km))
      Q2Km=c3*(T(i,j,km)+T(i+1,j-1,km)+T(i,j-1,km))
      Q3Km=c3*(T(i,j,km)+T(i-1,j,km)+T(i,j-1,km))
      Q4Km=c3*(T(i,j,km)+T(i-1,j,km)+T(i-1,j+1,km))
      Q5Km=c3*(T(i,j,km)+T(i,j+1,km)+T(i-1,j+1,km))
      Q6Km=c3*(T(i,j,km)+T(i+1,j,km)+T(i,j+1,km))

      W1Km=0.25*(2.*T(i,j,km)+T(i+1,j,km)+T(i+1,j-1,km))
      W2Km=0.25*(2.*T(i,j,km)+T(i+1,j-1,km)+T(i,j-1,km))
      W3Km=0.25*(2.*T(i,j,km)+T(i-1,j,km)+T(i,j-1,km))
      W4Km=0.25*(2.*T(i,j,km)+T(i-1,j,km)+T(i-1,j+1,km))
      W5Km=0.25*(2.*T(i,j,km)+T(i,j+1,km)+T(i-1,j+1,km))
      W6Km=0.25*(2.*T(i,j,km)+T(i+1,j,km)+T(i,j+1,km))

c     Triangle # 1.

      u1=c6*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)+
     +    u(i,j,km)+u(i+1,j,km)+u(i+1,j-1,km))
      v1=c6*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k)+
     +    v(i,j,km)+v(i+1,j,km)+v(i+1,j-1,km))
      w1=c3*(w(i,j,k)+w(i+1,j,k)+w(i+1,j-1,k))

      a11(1)= TAU*u1*u1 +A
      a12(1)= TAU*u1*v1
      a21(1)= TAU*u1*v1
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
      a21(2)= TAU*u2*v2
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
      a21(3)= TAU*u3*v3
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
      a21(4)= TAU*u4*v4
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
      a21(5)= TAU*u5*v5
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
      a21(6)= TAU*u6*v6
      a22(6)= TAU*v6*v6 +A
      a13(6)= TAU*u6*w6
      a23(6)= TAU*w6*v6 
      a33(6)= TAU*w6*w6

c     -------------------------------

c     Triangle # 1.

c     Horizontal part of the operator
      DIFT= DIFT +KT(1,n)*c6*(
     &a11(1)*(2.*(T(I+1,J,K)-T(I,J,K))+T(I+1,J,Km)-T(I,J,Km))/asr/S1
     &+a12(1)*(2.*(T(i+1,j,k)-T(i+1,j-1,k))+T(i+1,j,km)-T(i+1,j-1,km))
     &                        )*hz(km)  

c     dT/dz*dfi/dlambda
      DIFT= DIFT +KT(1,n)*Rhy*0.5*a13(1)*(Q1K-Q1Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(1,n)*Rhy*c6*a13(1)*(T(i+1,j,k )-T(i,j,k)+
     &                                   T(i+1,j,km)-T(i,j,km))

c     dfi/dz*dT/dteta
      DIFT=DIFT-KT(1,n)*Rhx*c6*S1*a23(1)*(T(i+1,j,k )-T(i+1,j-1,k)+
     &                                    T(i+1,j,km)-T(i+1,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(1,n)*R2hxhy*S1*a33(1)*(W1K-W1Km)/hz(km)

c     Triangle # 2.

c     Horizontal part of the operator
      DIFT= DIFT -KT(2,n)*c6*(
     &a22(2)*S2*(2.*(T(I,J,K)-T(I,J-1,K))+T(I,J,Km)-T(I,J-1,Km))*asr
     &+a21(2)*(2.*(T(i+1,j-1,k)-T(i,j-1,k))+T(i+1,j-1,km)-T(i,j-1,km))
     &                        )*hz(km)


c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(2,n)*Rhy*c6*a13(2)*(T(i+1,j-1,k )-T(i,j-1,k)+
     &                                   T(i+1,j-1,km)-T(i,j-1,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT -KT(2,n)*Rhx*S2*0.5*a23(2)*(Q2K-Q2Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(2,n)*Rhx*c6*S2*a23(2)*(T(i,j,k )-T(i,j-1,k)+
     &                                      T(i,j,km)-T(i,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(2,n)*R2hxhy*S2*a33(2)*(W2K-W2Km)/hz(km)

c     Triangle # 3.

c     Horizontal part of the operator
      DIFT= DIFT +KT(3,n)*c6*(
     & a11(3)*(2.*(T(I-1,J,K)-T(I,J,K)) +T(I-1,J,Km)-T(I,J,Km) )/asr/S3
     &+a12(3)*(2.*(T(i,j-1,k)-T(i  ,j,k)) +T(i,j-1,km)-T(i  ,j,km))
     &-a21(3)*(2.*(T(i,j  ,k)-T(i-1,j,k)) +T(i,j  ,km)-T(i-1,j,km))
     &-S3*asr*a22(3)*(2.*(T(I,J,K)-T(I,J-1,K)) +T(I,J,Km)-T(I,J-1,Km))
     &                        )*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT -KT(3,n)*Rhx*0.5*a13(3)*(Q3K-Q3Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(3,n)*Rhx*c6*a13(3)*(T(i,j,k )-T(i-1,j,k)+
     &                                   T(i,j,km)-T(i-1,j,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT -KT(3,n)*Rhy*S3*0.5*a23(3)*(Q3K-Q3Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(3,n)*Rhy*c6*S3*a23(3)*(T(i,j,k )-T(i,j-1,k)+
     &                                      T(i,j,km)-T(i,j-1,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(3,n)*R2hxhy*S3*a33(3)*(W3K-W3Km)/hz(km)

c     Triangle # 4.

c     Horizontal part of the operator
      DIFT= DIFT +KT(4,n)*c6*(
     & a11(4)*(2.*(T(I-1,J,K)-T(I,J,K))+T(I-1,J,Km)-T(I,J,Km))/asr/S4
     &+a12(4)*(2.*(T(i-1,j,k)-T(i-1,j+1,k))+T(i-1,j,km)-T(i-1,j+1,km))
     &                        )*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT -KT(4,n)*Rhy*0.5*a13(4)*(Q4K-Q4Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(4,n)*Rhy*c6*a13(4)*(T(i,j,k )-T(i-1,j,k)+
     &                                   T(i,j,km)-T(i-1,j,km))

c     dfi/dz*dT/dteta
      DIFT=DIFT-KT(4,n)*Rhx*c6*S4*a23(4)*(T(i-1,j+1,k )-T(i-1,j,k)+
     &                                    T(i-1,j+1,km)-T(i-1,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(4,n)*R2hxhy*S4*a33(4)*(W4K-W4Km)/hz(km)

c     Triangle # 5.

c     Horizontal part of the operator
      DIFT= DIFT -KT(5,n)*c6*(
     &asr*S5*a22(5)*(2.*(T(I,J,K)-T(I,J+1,K))+T(I,J,Km)-T(I,J+1,Km))
     &+a21(5)*(2.*(T(i-1,j+1,k)-T(i,j+1,k))+T(i-1,j+1,km)-T(i,j+1,km))
     &                        )*hz(km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(5,n)*Rhy*c6*a13(5)*(T(i,j+1,k )-T(i-1,j+1,k)+
     &                                   T(i,j+1,km)-T(i-1,j+1,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT +KT(5,n)*Rhx*S5*0.5*a23(5)*(Q5K-Q5Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(5,n)*Rhx*c6*S5*a23(5)*(T(i,j+1,k )-T(i,j,k)+
     &                                      T(i,j+1,km)-T(i,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(5,n)*R2hxhy*S5*a33(5)*(W5K-W5Km)/hz(km)

c     Triangle # 6.

c     Horizontal part of the operator
      DIFT= DIFT +KT(6,n)*c6*(
     & a11(6)*(2.*(T(I+1,J,K)-T(I,J,K)) +T(I+1,J,Km)-T(I,J,Km))/asr/S6
     &+a12(6)*(2.*(T(i,j+1,k)-T(i,j,k)) +T(i,j+1,km)-T(i,j,km))   
     &-a21(6)*(2.*(T(i,j,k)-T(i+1,j,k)) +T(i,j,km)-T(i+1,j,km))
     *-asr*S6*a22(6)*(2.*(T(I,J,K)-T(I,J+1,K)) +T(I,J,Km)-T(I,J+1,Km))
     &                        )*hz(km)

c     dT/dz*dfi/dlambda
      DIFT= DIFT +KT(6,n)*Rhy*0.5*a13(6)*(Q6K-Q6Km)

c     dT/dlambda*dfi/dz
      DIFT= DIFT -KT(6,n)*Rhy*c6*a13(6)*(T(i+1,j,k )-T(i,j,k)+
     &                                   T(i+1,j,km)-T(i,j,km))

c     dT/dz*dfi/dteta
      DIFT= DIFT +KT(6,n)*Rhx*S6*0.5*a23(6)*(Q6K-Q6Km)

c     dfi/dz*dT/dteta
      DIFT= DIFT -KT(6,n)*Rhx*c6*S6*a23(6)*(T(i,j+1,k )-T(i,j,k)+
     &                                      T(i,j+1,km)-T(i,j,km))

c     dT/dz*dfi/dz 
      DIFT= DIFT -KT(6,n)*R2hxhy*S6*a33(6)*(W6K-W6Km)/hz(km)

      END IF

      IF( K .LT. KB) THEN 
C     LOWER HALF

      Q1Kp=c3*(T(i,j,kp)+T(i+1,j,kp)+T(i+1,j-1,kp))
      Q2Kp=c3*(T(i,j,kp)+T(i+1,j-1,kp)+T(i,j-1,kp))
      Q3Kp=c3*(T(i,j,kp)+T(i-1,j,kp)+T(i,j-1,kp))
      Q4Kp=c3*(T(i,j,kp)+T(i-1,j,kp)+T(i-1,j+1,kp))
      Q5Kp=c3*(T(i,j,kp)+T(i,j+1,kp)+T(i-1,j+1,kp))
      Q6Kp=c3*(T(i,j,kp)+T(i+1,j,kp)+T(i,j+1,kp))

      W1Kp=0.25*(2.*T(i,j,kp)+T(i+1,j,kp)+T(i+1,j-1,kp))
      W2Kp=0.25*(2.*T(i,j,kp)+T(i+1,j-1,kp)+T(i,j-1,kp))
      W3Kp=0.25*(2.*T(i,j,kp)+T(i-1,j,kp)+T(i,j-1,kp))
      W4Kp=0.25*(2.*T(i,j,kp)+T(i-1,j,kp)+T(i-1,j+1,kp))
      W5Kp=0.25*(2.*T(i,j,kp)+T(i,j+1,kp)+T(i-1,j+1,kp))
      W6Kp=0.25*(2.*T(i,j,kp)+T(i+1,j,kp)+T(i,j+1,kp))

c     Triangle # 1.

      u1=c6*(u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)+
     +    u(i,j,kp)+u(i+1,j,kp)+u(i+1,j-1,kp))
      v1=c6*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k)+
     +    v(i,j,kp)+v(i+1,j,kp)+v(i+1,j-1,kp))
      w1=c3*(w(i,j,kp)+w(i+1,j,kp)+w(i+1,j-1,kp))

      a11(1)= TAU*u1*u1 +A
      a12(1)= TAU*u1*v1
      a21(1)= TAU*u1*v1
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
      a21(2)= TAU*u2*v2
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
      a21(3)= TAU*u3*v3
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
      a21(4)= TAU*u4*v4
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
      a21(5)= TAU*u5*v5
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
      a21(6)= TAU*u6*v6
      a22(6)= TAU*v6*v6 +A
      a13(6)= TAU*u6*w6
      a23(6)= TAU*w6*v6 
      a33(6)= TAU*w6*w6

c     -------------------------------

c     Triangle # 1.

c     Horizontal part of the operator
      DIFTP= DIFTP+c6*KT(1,np)*(
     & a11(1)*(2.*(T(I+1,J,K)-T(I,J,K)) +T(I+1,J,Kp)-T(I,J,Kp))/asr/S1
     &+a12(1)*(2.*(T(i+1,j,k)-T(i+1,j-1,k)) +T(i+1,j,kp)-T(i+1,j-1,kp))
     &                          )*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP +KT(1,np)*Rhy*0.5*a13(1)*(Q1Kp-Q1K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(1,np)*Rhy*c6*a13(1)*(T(i+1,j,k )-T(i,j,k)+
     &                                     T(i+1,j,kp)-T(i,j,kp))

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(1,np)*Rhx*c6*S1*a23(1)*(T(i+1,j,k )-T(i+1,j-1,k)+
     &                                       T(i+1,j,kp)-T(i+1,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(1,np)*R2hxhy*S1*a33(1)*(W1K-W1Kp)/hz(k)

c     Triangle # 2.

c     Horizontal part of the operator
      DIFTP= DIFTP-c6*KT(2,np)*(
     & a22(2)*S2*(2.*(T(I,J,K)-T(I,J-1,K))+T(I,J,Kp)-T(I,J-1,Kp))*ASR
     &+a21(2)*(2.*(T(i+1,j-1,k)-T(i,j-1,k))+T(i+1,j-1,kp)-T(i,j-1,kp))
     &                          )*hz(k)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(2,np)*Rhy*c6*a13(2)*
     &                                   (T(i+1,j-1,k )-T(i,j-1,k)+
     &                                    T(i+1,j-1,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP -KT(2,np)*Rhx*S2*0.5*a23(2)*(Q2Kp-Q2K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(2,np)*Rhx*c6*S2*a23(2)*
     &                                      (T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(2,np)*R2hxhy*S2*a33(2)*(W2K-W2Kp)/hz(k)

c     Triangle # 3.

c     Horizontal part of the operator
      DIFTP= DIFTP+c6*KT(3,np)*(
     & a11(3)*(2.*(T(I-1,J,K)-T(I,J,K))+T(I-1,J,Kp)-T(I,J,Kp))/asr/S3
     &+a12(3)*(2.*(T(i,j-1,k)-T(i,j,k))+T(i,j-1,kp)-T(i,j,kp))
     &-a21(3)*(2.*(T(i,j,k)-T(i-1,j,k))+T(i,j,kp)-T(i-1,j,kp))
     &-S3*a22(3)*(2.*(T(I,J,K)-T(I,J-1,K))+T(I,J,Kp)-T(I,J-1,Kp))*ASR
     & )*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP -KT(3,np)*Rhy*0.5*a13(3)*(Q3Kp-Q3K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(3,np)*Rhy*c6*a13(3)*(T(i,j,k )-T(i-1,j,k)+
     &                                      T(i,j,kp)-T(i-1,j,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP -KT(3,np)*Rhx*S3*0.5*a23(3)*(Q3Kp-Q3K)

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(3,np)*Rhx*c6*S3*a23(3)*(T(i,j,k )-T(i,j-1,k)+
     &                                       T(i,j,kp)-T(i,j-1,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(3,np)*R2hxhy*S3*a33(3)*(W3K-W3Kp)/hz(k)

c     Triangle # 4.

c     Horizontal part of the operator
      DIFTP= DIFTP+c6*KT(4,np)*(
     & a11(4)*(2.*(T(I-1,J,K)-T(I,J,K)) +T(I-1,J,Kp)-T(I,J,Kp))/asr/S4
     &+a12(4)*(2.*(T(i-1,j,k)-T(i-1,j+1,k))+T(i-1,j,kp)-T(i-1,j+1,kp))
     &                          )*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP -KT(4,np)*Rhy*0.5*a13(4)*(Q4Kp-Q4K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(4,np)*Rhy*c6*a13(4)*(T(i,j,k )-T(i-1,j,k)+
     &                                      T(i,j,kp)-T(i-1,j,kp))

c     dfi/dz*dT/dteta
      DIFTP=DIFTP+KT(4,np)*Rhx*c6*S4*a23(4)*
     &                                  (T(i-1,j+1,k )-T(i-1,j,k)+
     &                                   T(i-1,j+1,kp)-T(i-1,j,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(4,np)*R2hxhy*S4*a33(4)*(W4K-W4Kp)/hz(k)


c     Triangle # 5.

c     Horizontal part of the operator
      DIFTP= DIFTP-c6*KT(5,np)*(
     & a22(5)*S5*(2.*(T(I,J,K)-T(I,J+1,K)) +T(I,J,Kp)-T(I,J+1,Kp))*ASR
     &+a21(5)*(2.*(T(i-1,j+1,k)-T(i,j+1,k))+T(i-1,j+1,kp)-T(i,j+1,kp))
     &                          )*hz(k)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(5,np)*Rhy*c6*a13(5)*
     &                                   (T(i,j+1,k )-T(i-1,j+1,k)+
     &                                    T(i,j+1,kp)-T(i-1,j+1,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP +KT(5,np)*Rhx*S5*0.5*a23(5)*(Q5Kp-Q5K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(5,np)*Rhx*c6*S5*a23(5)*(T(i,j+1,k )-T(i,j,k)+
     &                                         T(i,j+1,kp)-T(i,j,kp))

c     dT/dz*dfi/dz 
      DIFTP=DIFTP-KT(5,np)*R2hxhy*S5*a33(5)*(W5K-W5Kp)/hz(k)


c     Triangle # 6.

c     Horizontal part of the operator
      DIFTP= DIFTP+c6*KT(6,np)*(
     & a11(6)*(2.*(T(I+1,J,K)-T(I,J,K)) +T(I+1,J,Kp)-T(I,J,Kp))/asr/S6
     &+a12(6)*(2.*(T(i,j+1,k)-T(i,j,k)) +T(i,j+1,kp)-T(i,j,kp))   
     &-a21(6)*(2.*(T(i,j,k)-T(i+1,j,k)) +T(i,j,kp)-T(i+1,j,kp))
     *-S6*a22(6)*(2.*(T(I,J,K)-T(I,J+1,K)) +T(I,J,Kp)-T(I,J+1,Kp))*ASR
     &                          )*hz(k)

c     dT/dz*dfi/dlambda
      DIFTP= DIFTP +KT(6,np)*Rhy*0.5*a13(6)*(Q6Kp-Q6K)

c     dT/dlambda*dfi/dz
      DIFTP= DIFTP +KT(6,np)*Rhy*c6*a13(6)*(T(i+1,j,k )-T(i,j,k)+
     &                                    T(i+1,j,kp)-T(i,j,kp))

c     dT/dz*dfi/dteta
      DIFTP= DIFTP +KT(6,np)*Rhx*S6*0.5*a23(6)*(Q6Kp-Q6K)

c     dfi/dz*dT/dteta
      DIFTP= DIFTP +KT(6,np)*Rhx*c6*S6*a23(6)*(T(i,j+1,k )-T(i,j,k)+
     &                                       T(i,j+1,kp)-T(i,j,kp))

c     dT/dz*dfi/dz 
      DIFTP= DIFTP -KT(6,np)*R2hxhy*S6*a33(6)*(W6K-W6Kp)/hz(k)

      END IF

      RETURN
      END

