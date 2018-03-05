      SUBROUTINE VVDIF
******************************************************
*     Vertical Diffusion of Ocean Momentum.
*     Free Upper Surface.
*     Nonlinear Metric Advection Approximation.
*     Version for standing alone ice and ocean models.
*     Explicit Scheme for ice-ocean and ocean-bottom drags.
*     Gravity Wave Drag for Ice-Ocean at the bottom of ML
*
*     Version 03.10.2015.
*
******************************************************
      INCLUDE 'Slo2.fi'
      
c     Dimensions for three-diagonal solver
      dimension AM(kl),BM(kl),CM(kl),FMx(kl),FMy(kl),Rksi(kl)
                  
      INCLUDE 'Tparm.fi'

*     ------- Additional Numerical viscosity ----------------
      Anumer= 0.

c     -----Parameters for quadratic bottom drag  ------------
**      Drag2= 2.6e-3  ! Kowalik&Polyakov
      Drag2= 1.0e-3  ! From Ibrayev (1.3), Clio 3.0, OPA 8.1 (1.0)
	gamma= 0.
c	gamma=10.
	cosg=cosd(gamma)
	sing=sind(gamma)
c     -----Parameters for linear bootom drag ----------------
*     Tfric - Damping time in the bottom layer.
*     Tfric=4 days, Numer.simulat.of topographic Rossby waves
*     along the East Greenland front. W.Maslowski. JGR, v.101,
*     No C4, pp.8775-8787, April 1996.
      Tfric= 4.0 *24.*3600.
c     -------------------------------------------------------

	c3=1./3.

      DO 1 J=1,JL
      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S2=c3*(2.*SM+S0)
      S3=S1
      S4=c3*(2.*S0+SP)
      S5=c3*(2.*SP+S0)
      S6=S4
      DO 1 I=1,IL

      KB = KM2(I,J)
      IF(KB.GT.0) THEN

*     ---------- Wind in [m/s] ------------
      wmod= SQRT( wx(i,j)**2 + wy(i,j)**2)
*     ------ Wind Stress Components -------
      cdrag= roa*Cdn(wmod,1)*wmod*1.e4/row
      windx= cdrag*wx(i,j)
      windy= cdrag*wy(i,j)
      Aopen= MIN(Aice(0,i,j), 1.0)
      
*     -------- Mixed layer thickness ------

*     Mixed layer thickness

      ro0=ropot(i,j,1) 

      kmix= kb
      
      do k=1,kb-1
      
c      ro1=ropot(i,j,k)   
c      ro2=ropot(i,j,k+1)
c      hzk= hz(k)   

c      VB= g*(ro2-ro1)/hzk/row
c      VB= SQRT(MAX(0.,VB))
      
c      if(vb .GT. 1.e-2) then
c      kmix= k
c      exit
c      end if


      ro2=ropot(i,j,k)
      kmix= k
       
      if( abs(ro2-ro0) .GT. 1.25e-4) exit ! Levitus, 1982

      end do	! k 
   
*     ---------------------------------------
*     Surface and interior of the ocean
*     ---------------------------------------

      DO 2 K=1, KB -1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))
            
	hzk  = hz(k)

      IF(k.EQ.1) THEN
        hzk1 = 0.
      ELSE
	  hzk1 = hz(k-1)
      END IF

      CA=2.*DT/(CG(n)*HZK1+CG(np)*HZK)

      IF( k.EQ.1) THEN
      
*     --------- viscous friction ------

*     Mean Drag coefficient under ice
	drag= Cdw*delta_u(i,j)*(1.-Aopen)

*     Levitating ice
      FMx(1)=u(i,j,1)+CA*cg(n)*(windx*Aopen+
     *                drag* uice(i,j) )
      FMy(1)=v(i,j,1)+CA*cg(n)*(windy*Aopen+
     *                drag* vice(i,j) )

      AM(1)=0.
      CM(1)=-CA*CG(n)*(Az(i,j,1) +Anumer)/hzk
      BM(1)=1.0-CM(1) +CA*CG(n)*drag 
            
      if(k.eq.kmix) then
      BM(1) = BM(1)  +CA*CG(n)*CDgwd(i,j)*(1.-Aopen)  
      FMx(1)= FMx(1)+
     +        CA*cg(n)*CDgwd(i,j)*(1.-Aopen)*uice(i,j)
      FMy(1)= FMy(1)+
     +        CA*cg(n)*CDgwd(i,j)*(1.-Aopen)*vice(i,j)
      end if

      ELSE ! k>1

      AM(k)=-CA*CG(n) *(az(i,j,k-1) +Anumer)/hzk1
      CM(k)=-CA*CG(np)*(az(i,j,k)   +Anumer)/hzk
      BM(k)= 1.-CM(k)-AM(k)
      FMx(k)=u(i,j,k)
      FMy(k)=v(i,j,k)
      
      if(k.eq.kmix) then
      BM(k) = BM(k) +CA*(CG(n)+CG(np))*CDgwd(i,j)*(1.-Aopen)   
      FMx(k)= FMx(k)+
     +  CA*(cg(n)+cg(np))*CDgwd(i,j)*(1.-Aopen)*uice(i,j)
      FMy(k)= FMy(k)+
     +  CA*(cg(n)+cg(np))*CDgwd(i,j)*(1.-Aopen)*vice(i,j)
      end if
      
      END IF

2     CONTINUE

*     ---------- Bottom ---------
      CA=2.*DT/hz(kb-1)


c     ---- Quadratic Bottom Drag ------
cc      BDrag=  Drag2*sqrt(u(i,j,kb)**2 +v(i,j,kb)**2)
      BDrag=  Drag2*sqrt(um2(i,j,kb)**2 +vm2(i,j,kb)**2)
cc	BDrag= MIN(BDRag, 0.25*hz(kb-1)/DT)

      AM(kb)=-CA*(az(i,j,kb-1)+Anumer)/hz(kb-1)
      BM(kb)= 1.-AM(kb)  +CA*BDrag 
      CM(kb)= 0.
      FMx(kb)= u(i,j,kb) !!!-CA*BDrag*(u(i,j,kb)*cosg+v(i,j,kb)*sing)
      FMy(kb)= v(i,j,kb) !!!-CA*BDrag*(v(i,j,kb)*cosg-u(i,j,kb)*sing)
      
      if(kb.eq.kmix) then
      BM(kb) = BM(kb) +CA*CDgwd(i,j)*(1.-Aopen)
      FMx(kb)=FMx(kb) +CA*CDgwd(i,j)*(1.-Aopen)*uice(i,j) 
      FMy(kb)=FMy(kb) +CA*CDgwd(i,j)*(1.-Aopen)*vice(i,j)
      end if
      

      call FACTOR(kl,am,bm,cm,fmx,rksi,1,kb)
      do k=1,kb
      u(i,j,k)= rksi(k)
      end do

      call FACTOR(kl,am,bm,cm,fmy,rksi,1,kb)
      do k=1,kb
      v(i,j,k)= rksi(k)
      end do

      END IF   ! kb>0
1     CONTINUE

      RETURN
      END

      subroutine vHdif(M)
c-----------------------------------------------------------
c     Horizontal diffusion and 3D transport of the momentum.
C     Free Upper Surface.
c     Modified Euler Scheme. M - number of step (1,2).
c     First step with 0.5*dt - Swansea time scheme
c     Metric terms in spherical coordinates.
c                Version 23.12.2015.
c-----------------------------------------------------------
      INCLUDE 'Slo2.fi'
      
*     Service arrays for Bi-Harmonic friction
      Grad2U(0:il1,0:jl1,kl), Grad2V(0:il1,0:jl1,kl)
      
CCC	double precision dzmean, dzmmean
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      common /diffus/ al,alt

	c6=1./6.
	c12= 0.5*c6
	c24= 0.5*c12

      HXX= 1./HX
      COEFX=HXX/18./R
      HYY= 1./HY
      COEFY=HYY/18./R
      R2=R*R
      asr=hx/hy
	asr2=asr*asr
      A=.5*AL*HXX*HXX/R2  ! For Standard Galerkin Scheme.
      
	Transp_Metric_U= 0.
	Transp_Metric_V= 0.
      
      S_Metric_Diff_U= 0. 
      S_Metric_Diff_V= 0. 
      
*     Friction 


                           DO 1 J=1,JL
	S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=(2.*S0+SM)/3.
      S2=(2.*SM+S0)/3.
      S3=S1
      S4=(2.*S0+SP)/3.
      S5=(2.*SP+S0)/3.
      S6=S4
      r2s6 = c6*s0

	COS_A= 0.5*(SP-SM)/hy
	CTG2 = (COS_A/S0)**2

                           DO 1 I=1,IL

      KB = KM2(I,J)
      IF(KB.GT.0) THEN
*     ---------------------------------------
*        Surface and interior of the ocean
*     ---------------------------------------
      n= ABS(nt3(i,j,1))

      DO K=1,KB-1
 	kp=k+1
	km=k-1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,kp))

	hzk = hz(k)

      IF( k .EQ. 1) THEN
	hzk1= 0.
      ELSE
      hzk1= hz(k-1)
      END IF
	 
      CA= REAL(m)*dt/( r2s6*(CG(n)*HZK1+CG(np)*HZK))
***      CA= 2.*dt/( r2s6*(CG(n)*HZK1+CG(np)*HZK))

      IF( M.EQ. 2) then
c     Standart Galerkin Scheme
      CALL DIFF(um2,I,J,K,N,NP,DIFU,DIFUP,KB,il1,jl1,kl,KT)
      CALL DIFF(vm2,I,J,K,N,NP,DIFV,DIFVP,KB,il1,jl1,kl,KT)
*     First metric
      CALL F_Metric(um2,vm2,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
     *                                   il1,jl1,kl,KT,COS_A,hx)
*     Second Metric
      S_Metric_Diff_U= AL*(1.-CTG2)*Um2(i,j,k)/R2
      S_Metric_Diff_V= AL*(1.-CTG2)*Vm2(i,j,k)/R2
      
      DIFU = A*DIFU
      DIFUP= A*DIFUP
      DIFV = A*DIFV
      DIFVP= A*DIFVP
            
      FU= .5*(DIFU*HZK1+DIFUP*HZ(K))  
      FV= .5*(DIFV*HZK1+DIFVP*HZ(K))  
      
      Grad2U(I,J,K)= CA*FU 
     *         + 0.5*REAL(m)*dt*(Transp_Metric_U  +S_Metric_Diff_U)
      Grad2V(I,J,K)= CA*FV 
     *         + 0.5*REAL(m)*dt*(Transp_Metric_V  +S_Metric_Diff_V)

      end do ! K
      
c     ----------------- Bottom ---------------
      N = NP
	k=kb
	km=kb-1

      CA=REAL(m)*DT/(R2S6*CG(n)*HZ(KB-1))


      IF(M.EQ.2) then
c     Standart Galerkin Scheme
      CALL DIFF(um2,I,J,K,N,NP,DIFU,DIFUP,KB,il1,jl1,kl,KT)
      CALL DIFF(vm2,I,J,K,N,NP,DIFV,DIFVP,KB,il1,jl1,kl,KT)
*     First metric
      CALL F_Metric(um2,vm2,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
     *                                    il1,jl1,kl,KT,COS_A,hx)
*     Second Metric
      S_Metric_Diff_U= AL*(1.-CTG2)*Um2(i,j,k)/R2
      S_Metric_Diff_V= AL*(1.-CTG2)*Vm2(i,j,k)/R2

      DIFU = A*DIFU
      DIFV = A*DIFV

      FU= .5*DIFU*HZ(Km)  
      FV= .5*DIFV*HZ(Km)  

      U(I,J,K)=Um2(I,J,K)+CA*FU 
     #          + 0.5*REAL(m)*dt*(Transp_Metric_U  +S_Metric_Diff_U)
      V(I,J,K)=Vm2(I,J,K)+CA*FV 
     #          + 0.5*REAL(m)*dt*(Transp_Metric_V  +S_Metric_Diff_V)
      
      
      end if
      end do ! i
      end do ! j

      
*    ------------------------ Time Steps -----------------      

                           DO 1 J=1,JL
	S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=(2.*S0+SM)/3.
      S2=(2.*SM+S0)/3.
      S3=S1
      S4=(2.*S0+SP)/3.
      S5=(2.*SP+S0)/3.
      S6=S4
      r2s6 = c6*s0

	COS_A= 0.5*(SP-SM)/hy
	CTG2 = (COS_A/S0)**2

                           DO 1 I=1,IL

      KB = KM2(I,J)
      IF(KB.GT.0) THEN
*     ---------------------------------------
*        Surface and interior of the ocean
*     ---------------------------------------
      n= ABS(nt3(i,j,1))

      DO 2 K=1,KB-1
 	kp=k+1
	km=k-1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,kp))

	hzk = hz(k)

      IF( k .EQ. 1) THEN
	hzk1= 0.
      ELSE
      hzk1= hz(k-1)
      END IF
	 
      CA= REAL(m)*dt/( r2s6*(CG(n)*HZK1+CG(np)*HZK))
***      CA= 2.*dt/( r2s6*(CG(n)*HZK1+CG(np)*HZK))

      FXUK =0.
      FXUKm =0.
      FXUK1=0.
      FXUK1p=0.
      FYUK =0.
      FYUKm =0.
      FYUK1=0.
      FYUK1p=0.
      FXVK =0.
      FXVKm =0.
      FXVK1=0.
      FXVK1p=0.
      FYVK =0.
      FYVKm =0.
      FYVK1=0.
      FYVK1p=0.

      IF( M.EQ. 2) then
c     Standart Galerkin Scheme
      CALL DIFF(um2,I,J,K,N,NP,DIFU,DIFUP,KB,il1,jl1,kl,KT)
      CALL DIFF(vm2,I,J,K,N,NP,DIFV,DIFVP,KB,il1,jl1,kl,KT)
*     First metric
      CALL F_Metric(um2,vm2,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
     *                                   il1,jl1,kl,KT,COS_A,hx)
*     Second Metric
      S_Metric_Diff_U= AL*(1.-CTG2)*Um2(i,j,k)/R2
      S_Metric_Diff_V= AL*(1.-CTG2)*Vm2(i,j,k)/R2

	ELSE
	
c     Standart Galerkin Scheme
c      CALL DIFF(um2,I,J,K,N,NP,DIFU,DIFUP,KB,il1,jl1,kl,KT)
c      CALL DIFF(vm2,I,J,K,N,NP,DIFV,DIFVP,KB,il1,jl1,kl,KT)
*     First metric
c      CALL F_Metric(um2,vm2,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
c     *                                   il1,jl1,kl,KT,COS_A,hx)
*     Second Metric
c      S_Metric_Diff_U= AL*(1.-CTG2)*Um2(i,j,k)/R2
c      S_Metric_Diff_V= AL*(1.-CTG2)*Vm2(i,j,k)/R2
	

      DIFU=0.
	DIFV=0.
	DIFUP=0.
	DIFVP=0.
	S_Metric_Diff_U=0.
	S_Metric_Diff_V=0.

	END IF

      if( k .GT. 1) then
      if(M.EQ.2) then
      call transp(FXUK,FYUK,i,j,k,n,um2,vm2,um1,KT,il1,jl1,kl)
      call transp(FXVK,FYVK,i,j,k,n,um2,vm2,vm1,KT,il1,jl1,kl)
      call transp(FXUKm,FYUKm,i,j,km,n,um2,vm2,um1,KT,il1,jl1,kl)
      call transp(FXVKm,FYVKm,i,j,km,n,um2,vm2,vm1,KT,il1,jl1,kl)
      else   ! M=2
      call transp(FXUK,FYUK,i,j,k,n,um2,vm2,um2,KT,il1,jl1,kl)  
      call transp(FXVK,FYVK,i,j,k,n,um2,vm2,vm2,KT,il1,jl1,kl)     
      call transp(FXUKm,FYUKm,i,j,km,n,um2,vm2,um2,KT,il1,jl1,kl)  
      call transp(FXVKm,FYVKm,i,j,km,n,um2,vm2,vm2,KT,il1,jl1,kl)     
      end if ! M=1
      end if  ! K>1

      if(M.EQ.2) then
      call transp(FXUK1,FYUK1,i,j,k,np,um2,vm2,um1,KT,il1,jl1,kl)
      call transp(FXVK1,FYVK1,i,j,k,np,um2,vm2,vm1,KT,il1,jl1,kl)
      call transp(FXUK1p,FYUK1p,i,j,kp,np,um2,vm2,um1,KT,il1,jl1,kl)
      call transp(FXVK1p,FYVK1p,i,j,kp,np,um2,vm2,vm1,KT,il1,jl1,kl)
      else   ! M=2
      call transp(FXUK1,FYUK1,i,j,k,np,um2,vm2,um2,KT,il1,jl1,kl)
      call transp(FXVK1,FYVK1,i,j,k,np,um2,vm2,vm2,KT,il1,jl1,kl)
      call transp(FXUK1p,FYUK1p,i,j,kp,np,um2,vm2,um2,KT,il1,jl1,kl)
      call transp(FXVK1p,FYVK1p,i,j,kp,np,um2,vm2,vm2,KT,il1,jl1,kl)
      end if ! M=1

      FXUK = COEFX*FXUK
      FYUK = COEFY*FYUK
      FXUK1= COEFX*FXUK1
      FYUK1= COEFY*FYUK1
      FXUKm = COEFX*FXUKm
      FYUKm = COEFY*FYUKm
      FXUK1p= COEFX*FXUK1p
      FYUK1p= COEFY*FYUK1p

      FXVK = COEFX*FXVK
      FYVK = COEFY*FYVK
      FXVK1= COEFX*FXVK1
      FYVK1= COEFY*FYVK1
      FXVKm = COEFX*FXVKm
      FYVKm = COEFY*FYVKm
      FXVK1p= COEFX*FXVK1p
      FYVK1p= COEFY*FYVK1p

      DIFU = A*DIFU
      DIFUP= A*DIFUP
      DIFV = A*DIFV
      DIFVP= A*DIFVP

      FXU=c6*((2.*FXUK+FXUKm)*HZK1 +(2.*FXUK1+FXUK1p)*HZK)
      FYU=c6*((2.*FYUK+FYUKm)*HZK1 +(2.*FYUK1+FYUK1p)*HZK)
      FXV=c6*((2.*FXVK+FXVKm)*HZK1 +(2.*FXVK1+FXVK1p)*HZK)
      FYV=c6*((2.*FYVK+FYVKm)*HZK1 +(2.*FYVK1+FYVK1p)*HZK)

      FZU= 0.
      FZV= 0.

      IF( k .EQ. 1) THEN
*     -------------- Surface -----------------
 
 
      IF( M .EQ. 2) THEN

      Sum1= c24*
     * (KT(1,n)*S1*(2.*um1(i,j,1)+um1(i+1,j,1)+um1(i+1,j-1,1))
     + +KT(2,n)*S2*(2.*um1(i,j,1)+um1(i,j-1,1)+um1(i+1,j-1,1))
     + +KT(3,n)*S3*(2.*um1(i,j,1)+um1(i,j-1,1)+um1(i-1,j  ,1))
     + +KT(4,n)*S4*(2.*um1(i,j,1)+um1(i-1,j,1)+um1(i-1,j+1,1))
     + +KT(5,n)*S5*(2.*um1(i,j,1)+um1(i-1,j+1,1)+um1(i,j+1,1))
     + +KT(6,n)*S6*(2.*um1(i,j,1)+um1(i,j+1,1)+um1(i+1,j  ,1)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*um1(i,j,2)+um1(i+1,j,2)+um1(i+1,j-1,2))
     + +KT(2,n)*S2*(2.*um1(i,j,2)+um1(i,j-1,2)+um1(i+1,j-1,2))
     + +KT(3,n)*S3*(2.*um1(i,j,2)+um1(i,j-1,2)+um1(i-1,j  ,2))
     + +KT(4,n)*S4*(2.*um1(i,j,2)+um1(i-1,j,2)+um1(i-1,j+1,2))
     + +KT(5,n)*S5*(2.*um1(i,j,2)+um1(i-1,j+1,2)+um1(i,j+1,2))
     + +KT(6,n)*S6*(2.*um1(i,j,2)+um1(i,j+1,2)+um1(i+1,j  ,2)))

      FZu=-0.5*((Sum1+Sum2)*W(i,j,2)-2.*Sum1*(W(i,j,1)-PME(i,j)/dt))

      Sum1= c24*
     * (KT(1,n)*S1*(2.*vm1(i,j,1)+vm1(i+1,j,1)+vm1(i+1,j-1,1))
     + +KT(2,n)*S2*(2.*vm1(i,j,1)+vm1(i,j-1,1)+vm1(i+1,j-1,1))
     + +KT(3,n)*S3*(2.*vm1(i,j,1)+vm1(i,j-1,1)+vm1(i-1,j  ,1))
     + +KT(4,n)*S4*(2.*vm1(i,j,1)+vm1(i-1,j,1)+vm1(i-1,j+1,1))
     + +KT(5,n)*S5*(2.*vm1(i,j,1)+vm1(i-1,j+1,1)+vm1(i,j+1,1))
     + +KT(6,n)*S6*(2.*vm1(i,j,1)+vm1(i,j+1,1)+vm1(i+1,j  ,1)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*vm1(i,j,2)+vm1(i+1,j,2)+vm1(i+1,j-1,2))
     + +KT(2,n)*S2*(2.*vm1(i,j,2)+vm1(i,j-1,2)+vm1(i+1,j-1,2))
     + +KT(3,n)*S3*(2.*vm1(i,j,2)+vm1(i,j-1,2)+vm1(i-1,j  ,2))
     + +KT(4,n)*S4*(2.*vm1(i,j,2)+vm1(i-1,j,2)+vm1(i-1,j+1,2))
     + +KT(5,n)*S5*(2.*vm1(i,j,2)+vm1(i-1,j+1,2)+vm1(i,j+1,2))
     + +KT(6,n)*S6*(2.*vm1(i,j,2)+vm1(i,j+1,2)+vm1(i+1,j  ,2)))

      FZv=-0.5*((Sum1+Sum2)*W(i,j,2)-2.*Sum1*(W(i,j,1)-PME(i,j)/dt))

*     Metric in Transport
	Transp_Metric_U= -Um2(i,j,1)*Vm1(i,j,1)*COS_A/(R*S0)
	Transp_Metric_V= +Um2(i,j,1)*Um1(i,j,1)*COS_A/(R*S0)

      ELSE

      Sum1= c24*
     * (KT(1,n)*S1*(2.*um2(i,j,1)+um2(i+1,j,1)+um2(i+1,j-1,1))
     + +KT(2,n)*S2*(2.*um2(i,j,1)+um2(i,j-1,1)+um2(i+1,j-1,1))
     + +KT(3,n)*S3*(2.*um2(i,j,1)+um2(i,j-1,1)+um2(i-1,j  ,1))
     + +KT(4,n)*S4*(2.*um2(i,j,1)+um2(i-1,j,1)+um2(i-1,j+1,1))
     + +KT(5,n)*S5*(2.*um2(i,j,1)+um2(i-1,j+1,1)+um2(i,j+1,1))
     + +KT(6,n)*S6*(2.*um2(i,j,1)+um2(i,j+1,1)+um2(i+1,j  ,1)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*um2(i,j,2)+um2(i+1,j,2)+um2(i+1,j-1,2))
     + +KT(2,n)*S2*(2.*um2(i,j,2)+um2(i,j-1,2)+um2(i+1,j-1,2))
     + +KT(3,n)*S3*(2.*um2(i,j,2)+um2(i,j-1,2)+um2(i-1,j  ,2))
     + +KT(4,n)*S4*(2.*um2(i,j,2)+um2(i-1,j,2)+um2(i-1,j+1,2))
     + +KT(5,n)*S5*(2.*um2(i,j,2)+um2(i-1,j+1,2)+um2(i,j+1,2))
     + +KT(6,n)*S6*(2.*um2(i,j,2)+um2(i,j+1,2)+um2(i+1,j  ,2)))

      FZu=-0.5*((Sum1+Sum2)*W(i,j,2)-2.*Sum1*(W(i,j,1)-PME(i,j)/dt))

      Sum1= c24*
     * (KT(1,n)*S1*(2.*vm2(i,j,1)+vm2(i+1,j,1)+vm2(i+1,j-1,1))
     + +KT(2,n)*S2*(2.*vm2(i,j,1)+vm2(i,j-1,1)+vm2(i+1,j-1,1))
     + +KT(3,n)*S3*(2.*vm2(i,j,1)+vm2(i,j-1,1)+vm2(i-1,j  ,1))
     + +KT(4,n)*S4*(2.*vm2(i,j,1)+vm2(i-1,j,1)+vm2(i-1,j+1,1))
     + +KT(5,n)*S5*(2.*vm2(i,j,1)+vm2(i-1,j+1,1)+vm2(i,j+1,1))
     + +KT(6,n)*S6*(2.*vm2(i,j,1)+vm2(i,j+1,1)+vm2(i+1,j  ,1)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*vm2(i,j,2)+vm2(i+1,j,2)+vm2(i+1,j-1,2))
     + +KT(2,n)*S2*(2.*vm2(i,j,2)+vm2(i,j-1,2)+vm2(i+1,j-1,2))
     + +KT(3,n)*S3*(2.*vm2(i,j,2)+vm2(i,j-1,2)+vm2(i-1,j  ,2))
     + +KT(4,n)*S4*(2.*vm2(i,j,2)+vm2(i-1,j,2)+vm2(i-1,j+1,2))
     + +KT(5,n)*S5*(2.*vm2(i,j,2)+vm2(i-1,j+1,2)+vm2(i,j+1,2))
     + +KT(6,n)*S6*(2.*vm2(i,j,2)+vm2(i,j+1,2)+vm2(i+1,j  ,2)))

      FZv=-0.5*((Sum1+Sum2)*W(i,j,2)-2.*Sum1*(W(i,j,1)-PME(i,j)/dt))

*     Metric in Transport
	Transp_Metric_U= -Um2(i,j,1)*Vm2(i,j,1)*COS_A/(R*S0)
	Transp_Metric_V= +Um2(i,j,1)*Um2(i,j,1)*COS_A/(R*S0)
      END IF

      ELSE


*     ------------------ Interior ----------------

      IF( M .EQ. 2) THEN


      Sum1= c24*
     * (KT(1,n)*S1*(2.*um1(i,j,k)+um1(i+1,j,k)+um1(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*um1(i,j,k)+um1(i,j-1,k)+um1(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*um1(i,j,k)+um1(i,j-1,k)+um1(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*um1(i,j,k)+um1(i-1,j,k)+um1(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*um1(i,j,k)+um1(i-1,j+1,k)+um1(i,j+1,k))
     + +KT(6,n)*S6*(2.*um1(i,j,k)+um1(i,j+1,k)+um1(i+1,j  ,k)))

      Sum1m= c24*
     * (KT(1,n)*S1*(2.*um1(i,j,km)+um1(i+1,j,km)+um1(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*um1(i,j,km)+um1(i,j-1,km)+um1(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*um1(i,j,km)+um1(i,j-1,km)+um1(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*um1(i,j,km)+um1(i-1,j,km)+um1(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*um1(i,j,km)+um1(i-1,j+1,km)+um1(i,j+1,km))
     + +KT(6,n)*S6*(2.*um1(i,j,km)+um1(i,j+1,km)+um1(i+1,j  ,km)))

      Sum2= c24*
     * (KT(1,np)*S1*(2.*um1(i,j,k)+um1(i+1,j,k)+um1(i+1,j-1,k))
     + +KT(2,np)*S2*(2.*um1(i,j,k)+um1(i,j-1,k)+um1(i+1,j-1,k))
     + +KT(3,np)*S3*(2.*um1(i,j,k)+um1(i,j-1,k)+um1(i-1,j  ,k))
     + +KT(4,np)*S4*(2.*um1(i,j,k)+um1(i-1,j,k)+um1(i-1,j+1,k))
     + +KT(5,np)*S5*(2.*um1(i,j,k)+um1(i-1,j+1,k)+um1(i,j+1,k))
     + +KT(6,np)*S6*(2.*um1(i,j,k)+um1(i,j+1,k)+um1(i+1,j  ,k)))

      Sum2p= c24*
     * (KT(1,np)*S1*(2.*um1(i,j,kp)+um1(i+1,j,kp)+um1(i+1,j-1,kp))
     + +KT(2,np)*S2*(2.*um1(i,j,kp)+um1(i,j-1,kp)+um1(i+1,j-1,kp))
     + +KT(3,np)*S3*(2.*um1(i,j,kp)+um1(i,j-1,kp)+um1(i-1,j  ,kp))
     + +KT(4,np)*S4*(2.*um1(i,j,kp)+um1(i-1,j,kp)+um1(i-1,j+1,kp))
     + +KT(5,np)*S5*(2.*um1(i,j,kp)+um1(i-1,j+1,kp)+um1(i,j+1,kp))
     + +KT(6,np)*S6*(2.*um1(i,j,kp)+um1(i,j+1,kp)+um1(i+1,j  ,kp)))

      FZU=-0.5*( (Sum2+Sum2p)*W(i,j,Kp) -(Sum1+Sum1m)*W(i,j,K) )


      Sum1= c24*
     * (KT(1,n)*S1*(2.*vm1(i,j,k)+vm1(i+1,j,k)+vm1(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*vm1(i,j,k)+vm1(i,j-1,k)+vm1(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*vm1(i,j,k)+vm1(i,j-1,k)+vm1(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*vm1(i,j,k)+vm1(i-1,j,k)+vm1(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*vm1(i,j,k)+vm1(i-1,j+1,k)+vm1(i,j+1,k))
     + +KT(6,n)*S6*(2.*vm1(i,j,k)+vm1(i,j+1,k)+vm1(i+1,j  ,k)))

      Sum1m= c24*
     * (KT(1,n)*S1*(2.*vm1(i,j,km)+vm1(i+1,j,km)+vm1(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*vm1(i,j,km)+vm1(i,j-1,km)+vm1(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*vm1(i,j,km)+vm1(i,j-1,km)+vm1(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*vm1(i,j,km)+vm1(i-1,j,km)+vm1(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*vm1(i,j,km)+vm1(i-1,j+1,km)+vm1(i,j+1,km))
     + +KT(6,n)*S6*(2.*vm1(i,j,km)+vm1(i,j+1,km)+vm1(i+1,j  ,km)))

      Sum2= c24*
     * (KT(1,np)*S1*(2.*vm1(i,j,k)+vm1(i+1,j,k)+vm1(i+1,j-1,k))
     + +KT(2,np)*S2*(2.*vm1(i,j,k)+vm1(i,j-1,k)+vm1(i+1,j-1,k))
     + +KT(3,np)*S3*(2.*vm1(i,j,k)+vm1(i,j-1,k)+vm1(i-1,j  ,k))
     + +KT(4,np)*S4*(2.*vm1(i,j,k)+vm1(i-1,j,k)+vm1(i-1,j+1,k))
     + +KT(5,np)*S5*(2.*vm1(i,j,k)+vm1(i-1,j+1,k)+vm1(i,j+1,k))
     + +KT(6,np)*S6*(2.*vm1(i,j,k)+vm1(i,j+1,k)+vm1(i+1,j  ,k)))

      Sum2p= c24*
     * (KT(1,np)*S1*(2.*vm1(i,j,kp)+vm1(i+1,j,kp)+vm1(i+1,j-1,kp))
     + +KT(2,np)*S2*(2.*vm1(i,j,kp)+vm1(i,j-1,kp)+vm1(i+1,j-1,kp))
     + +KT(3,np)*S3*(2.*vm1(i,j,kp)+vm1(i,j-1,kp)+vm1(i-1,j  ,kp))
     + +KT(4,np)*S4*(2.*vm1(i,j,kp)+vm1(i-1,j,kp)+vm1(i-1,j+1,kp))
     + +KT(5,np)*S5*(2.*vm1(i,j,kp)+vm1(i-1,j+1,kp)+vm1(i,j+1,kp))
     + +KT(6,np)*S6*(2.*vm1(i,j,kp)+vm1(i,j+1,kp)+vm1(i+1,j  ,kp)))

      FZV=-0.5*( (Sum2+Sum2p)*W(i,j,Kp) -(Sum1+Sum1m)*W(i,j,K) )

*     Metric in Transport
	Transp_Metric_U= -Um2(i,j,k)*Vm1(i,j,k)*COS_A/(R*S0)
	Transp_Metric_V= +Um2(i,j,k)*Um1(i,j,k)*COS_A/(R*S0)


      ELSE

      Sum1= c24*
     * (KT(1,n)*S1*(2.*um2(i,j,k)+um2(i+1,j,k)+um2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*um2(i,j,k)+um2(i,j-1,k)+um2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*um2(i,j,k)+um2(i,j-1,k)+um2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*um2(i,j,k)+um2(i-1,j,k)+um2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*um2(i,j,k)+um2(i-1,j+1,k)+um2(i,j+1,k))
     + +KT(6,n)*S6*(2.*um2(i,j,k)+um2(i,j+1,k)+um2(i+1,j  ,k)))

      Sum1m= c24*
     * (KT(1,n)*S1*(2.*um2(i,j,km)+um2(i+1,j,km)+um2(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*um2(i,j,km)+um2(i,j-1,km)+um2(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*um2(i,j,km)+um2(i,j-1,km)+um2(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*um2(i,j,km)+um2(i-1,j,km)+um2(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*um2(i,j,km)+um2(i-1,j+1,km)+um2(i,j+1,km))
     + +KT(6,n)*S6*(2.*um2(i,j,km)+um2(i,j+1,km)+um2(i+1,j  ,km)))

      Sum2= c24*
     * (KT(1,np)*S1*(2.*um2(i,j,k)+um2(i+1,j,k)+um2(i+1,j-1,k))
     + +KT(2,np)*S2*(2.*um2(i,j,k)+um2(i,j-1,k)+um2(i+1,j-1,k))
     + +KT(3,np)*S3*(2.*um2(i,j,k)+um2(i,j-1,k)+um2(i-1,j  ,k))
     + +KT(4,np)*S4*(2.*um2(i,j,k)+um2(i-1,j,k)+um2(i-1,j+1,k))
     + +KT(5,np)*S5*(2.*um2(i,j,k)+um2(i-1,j+1,k)+um2(i,j+1,k))
     + +KT(6,np)*S6*(2.*um2(i,j,k)+um2(i,j+1,k)+um2(i+1,j  ,k)))

      Sum2p= c24*
     * (KT(1,np)*S1*(2.*um2(i,j,kp)+um2(i+1,j,kp)+um2(i+1,j-1,kp))
     + +KT(2,np)*S2*(2.*um2(i,j,kp)+um2(i,j-1,kp)+um2(i+1,j-1,kp))
     + +KT(3,np)*S3*(2.*um2(i,j,kp)+um2(i,j-1,kp)+um2(i-1,j  ,kp))
     + +KT(4,np)*S4*(2.*um2(i,j,kp)+um2(i-1,j,kp)+um2(i-1,j+1,kp))
     + +KT(5,np)*S5*(2.*um2(i,j,kp)+um2(i-1,j+1,kp)+um2(i,j+1,kp))
     + +KT(6,np)*S6*(2.*um2(i,j,kp)+um2(i,j+1,kp)+um2(i+1,j  ,kp)))

      FZU=-0.5*( (Sum2+Sum2p)*W(i,j,Kp) -(Sum1+Sum1m)*W(i,j,K) )


      Sum1= c24*
     * (KT(1,n)*S1*(2.*vm2(i,j,k)+vm2(i+1,j,k)+vm2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*vm2(i,j,k)+vm2(i,j-1,k)+vm2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*vm2(i,j,k)+vm2(i,j-1,k)+vm2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*vm2(i,j,k)+vm2(i-1,j,k)+vm2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*vm2(i,j,k)+vm2(i-1,j+1,k)+vm2(i,j+1,k))
     + +KT(6,n)*S6*(2.*vm2(i,j,k)+vm2(i,j+1,k)+vm2(i+1,j  ,k)))

      Sum1m= c24*
     * (KT(1,n)*S1*(2.*vm2(i,j,km)+vm2(i+1,j,km)+vm2(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*vm2(i,j,km)+vm2(i,j-1,km)+vm2(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*vm2(i,j,km)+vm2(i,j-1,km)+vm2(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*vm2(i,j,km)+vm2(i-1,j,km)+vm2(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*vm2(i,j,km)+vm2(i-1,j+1,km)+vm2(i,j+1,km))
     + +KT(6,n)*S6*(2.*vm2(i,j,km)+vm2(i,j+1,km)+vm2(i+1,j  ,km)))

      Sum2= c24*
     * (KT(1,np)*S1*(2.*vm2(i,j,k)+vm2(i+1,j,k)+vm2(i+1,j-1,k))
     + +KT(2,np)*S2*(2.*vm2(i,j,k)+vm2(i,j-1,k)+vm2(i+1,j-1,k))
     + +KT(3,np)*S3*(2.*vm2(i,j,k)+vm2(i,j-1,k)+vm2(i-1,j  ,k))
     + +KT(4,np)*S4*(2.*vm2(i,j,k)+vm2(i-1,j,k)+vm2(i-1,j+1,k))
     + +KT(5,np)*S5*(2.*vm2(i,j,k)+vm2(i-1,j+1,k)+vm2(i,j+1,k))
     + +KT(6,np)*S6*(2.*vm2(i,j,k)+vm2(i,j+1,k)+vm2(i+1,j  ,k)))

      Sum2p= c24*
     * (KT(1,np)*S1*(2.*vm2(i,j,kp)+vm2(i+1,j,kp)+vm2(i+1,j-1,kp))
     + +KT(2,np)*S2*(2.*vm2(i,j,kp)+vm2(i,j-1,kp)+vm2(i+1,j-1,kp))
     + +KT(3,np)*S3*(2.*vm2(i,j,kp)+vm2(i,j-1,kp)+vm2(i-1,j  ,kp))
     + +KT(4,np)*S4*(2.*vm2(i,j,kp)+vm2(i-1,j,kp)+vm2(i-1,j+1,kp))
     + +KT(5,np)*S5*(2.*vm2(i,j,kp)+vm2(i-1,j+1,kp)+vm2(i,j+1,kp))
     + +KT(6,np)*S6*(2.*vm2(i,j,kp)+vm2(i,j+1,kp)+vm2(i+1,j  ,kp)))

      FZV=-0.5*( (Sum2+Sum2p)*W(i,j,Kp) -(Sum1+Sum1m)*W(i,j,K) )

*     Metric in Transport
	Transp_Metric_U= -Um2(i,j,k)*Vm2(i,j,k)*COS_A/(R*S0)
	Transp_Metric_V= +Um2(i,j,k)*Um2(i,j,k)*COS_A/(R*S0)
      END IF
      END IF

      FU= .5*(DIFU*HZK1+DIFUP*HZ(K))  +FZU-FXU-FYU
      FV= .5*(DIFV*HZK1+DIFVP*HZ(K))  +FZV-FXV-FYV
      U(I,J,K)=Um2(I,J,K)+CA*FU 
     *         + 0.5*REAL(m)*dt*(Transp_Metric_U  +S_Metric_Diff_U)
      V(I,J,K)=Vm2(I,J,K)+CA*FV 
     *         + 0.5*REAL(m)*dt*(Transp_Metric_V  +S_Metric_Diff_V)

2     CONTINUE

c     ----------------- Bottom ---------------
      N = NP
	k=kb
	km=kb-1

      CA=REAL(m)*DT/(R2S6*CG(n)*HZ(KB-1))

      FXUK= 0.
      FYUK= 0.
      FXVK= 0.
      FYVK= 0.
      FXUKm= 0.
      FYUKm= 0.
      FXVKm= 0.
      FYVKm= 0.


      IF(M.EQ.2) then
c     Standart Galerkin Scheme
      CALL DIFF(um2,I,J,K,N,NP,DIFU,DIFUP,KB,il1,jl1,kl,KT)
      CALL DIFF(vm2,I,J,K,N,NP,DIFV,DIFVP,KB,il1,jl1,kl,KT)
*     First metric
      CALL F_Metric(um2,vm2,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
     *                                    il1,jl1,kl,KT,COS_A,hx)
*     Second Metric
      S_Metric_Diff_U= AL*(1.-CTG2)*Um2(i,j,k)/R2
      S_Metric_Diff_V= AL*(1.-CTG2)*Vm2(i,j,k)/R2

      ELSE

c     Standart Galerkin Scheme
c      CALL DIFF(um2,I,J,K,N,NP,DIFU,DIFUP,KB,il1,jl1,kl,KT)
c      CALL DIFF(vm2,I,J,K,N,NP,DIFV,DIFVP,KB,il1,jl1,kl,KT)
*     First metric
c      CALL F_Metric(um2,vm2,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
c     *                                    il1,jl1,kl,KT,COS_A,hx)
*     Second Metric
c      S_Metric_Diff_U= AL*(1.-CTG2)*Um2(i,j,k)/R2
c      S_Metric_Diff_V= AL*(1.-CTG2)*Vm2(i,j,k)/R2

      DIFU=0.
	DIFV=0.
	DIFUP=0.
	DIFVP=0.
	S_Metric_Diff_U=0.
	S_Metric_Diff_V=0.

	END IF

      IF( M .EQ. 1) THEN
      call transp(FXUK,FYUK,i,j,k,n,um2,vm2,um2,KT,il1,jl1,kl)
      call transp(FXVK,FYVK,i,j,k,n,um2,vm2,vm2,KT,il1,jl1,kl)
      call transp(FXUKm,FYUKm,i,j,km,n,um2,vm2,um2,KT,il1,jl1,kl)
      call transp(FXVKm,FYVKm,i,j,km,n,um2,vm2,vm2,KT,il1,jl1,kl)

      Sum1= c24*
     * (KT(1,n)*S1*(2.*um2(i,j,k)+um2(i+1,j,k)+um2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*um2(i,j,k)+um2(i,j-1,k)+um2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*um2(i,j,k)+um2(i,j-1,k)+um2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*um2(i,j,k)+um2(i-1,j,k)+um2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*um2(i,j,k)+um2(i-1,j+1,k)+um2(i,j+1,k))
     + +KT(6,n)*S6*(2.*um2(i,j,k)+um2(i,j+1,k)+um2(i+1,j  ,k)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*um2(i,j,km)+um2(i+1,j,km)+um2(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*um2(i,j,km)+um2(i,j-1,km)+um2(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*um2(i,j,km)+um2(i,j-1,km)+um2(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*um2(i,j,km)+um2(i-1,j,km)+um2(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*um2(i,j,km)+um2(i-1,j+1,km)+um2(i,j+1,km))
     + +KT(6,n)*S6*(2.*um2(i,j,km)+um2(i,j+1,km)+um2(i+1,j  ,km)))

      FZU=0.5*(Sum1+Sum2)*W(i,j,K)

      Sum1= c24*
     * (KT(1,n)*S1*(2.*vm2(i,j,k)+vm2(i+1,j,k)+vm2(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*vm2(i,j,k)+vm2(i,j-1,k)+vm2(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*vm2(i,j,k)+vm2(i,j-1,k)+vm2(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*vm2(i,j,k)+vm2(i-1,j,k)+vm2(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*vm2(i,j,k)+vm2(i-1,j+1,k)+vm2(i,j+1,k))
     + +KT(6,n)*S6*(2.*vm2(i,j,k)+vm2(i,j+1,k)+vm2(i+1,j  ,k)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*vm2(i,j,km)+vm2(i+1,j,km)+vm2(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*vm2(i,j,km)+vm2(i,j-1,km)+vm2(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*vm2(i,j,km)+vm2(i,j-1,km)+vm2(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*vm2(i,j,km)+vm2(i-1,j,km)+vm2(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*vm2(i,j,km)+vm2(i-1,j+1,km)+vm2(i,j+1,km))
     + +KT(6,n)*S6*(2.*vm2(i,j,km)+vm2(i,j+1,km)+vm2(i+1,j  ,km)))

      FZV=0.5*(Sum1+Sum2)*W(i,j,K)


*     Metric in Transport
	Transp_Metric_U= -Um2(i,j,k)*Vm2(i,j,k)*COS_A/(R*S0)
	Transp_Metric_V= +Um2(i,j,k)*Um2(i,j,k)*COS_A/(R*S0)

      ELSE

      call transp(FXUK,FYUK,i,j,k,n,um2,vm2,um1,KT,il1,jl1,kl)
      call transp(FXVK,FYVK,i,j,k,n,um2,vm2,vm1,KT,il1,jl1,kl)
      call transp(FXUKm,FYUKm,i,j,km,n,um2,vm2,um1,KT,il1,jl1,kl)
      call transp(FXVKm,FYVKm,i,j,km,n,um2,vm2,vm1,KT,il1,jl1,kl)

      Sum1= c24*
     * (KT(1,n)*S1*(2.*um1(i,j,k)+um1(i+1,j,k)+um1(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*um1(i,j,k)+um1(i,j-1,k)+um1(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*um1(i,j,k)+um1(i,j-1,k)+um1(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*um1(i,j,k)+um1(i-1,j,k)+um1(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*um1(i,j,k)+um1(i-1,j+1,k)+um1(i,j+1,k))
     + +KT(6,n)*S6*(2.*um1(i,j,k)+um1(i,j+1,k)+um1(i+1,j  ,k)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*um1(i,j,km)+um1(i+1,j,km)+um1(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*um1(i,j,km)+um1(i,j-1,km)+um1(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*um1(i,j,km)+um1(i,j-1,km)+um1(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*um1(i,j,km)+um1(i-1,j,km)+um1(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*um1(i,j,km)+um1(i-1,j+1,km)+um1(i,j+1,km))
     + +KT(6,n)*S6*(2.*um1(i,j,km)+um1(i,j+1,km)+um1(i+1,j  ,km)))

      FZU=0.5*(Sum1+Sum2)*W(i,j,K)

      Sum1= c24*
     * (KT(1,n)*S1*(2.*vm1(i,j,k)+vm1(i+1,j,k)+vm1(i+1,j-1,k))
     + +KT(2,n)*S2*(2.*vm1(i,j,k)+vm1(i,j-1,k)+vm1(i+1,j-1,k))
     + +KT(3,n)*S3*(2.*vm1(i,j,k)+vm1(i,j-1,k)+vm1(i-1,j  ,k))
     + +KT(4,n)*S4*(2.*vm1(i,j,k)+vm1(i-1,j,k)+vm1(i-1,j+1,k))
     + +KT(5,n)*S5*(2.*vm1(i,j,k)+vm1(i-1,j+1,k)+vm1(i,j+1,k))
     + +KT(6,n)*S6*(2.*vm1(i,j,k)+vm1(i,j+1,k)+vm1(i+1,j  ,k)))

      Sum2= c24*
     * (KT(1,n)*S1*(2.*vm1(i,j,km)+vm1(i+1,j,km)+vm1(i+1,j-1,km))
     + +KT(2,n)*S2*(2.*vm1(i,j,km)+vm1(i,j-1,km)+vm1(i+1,j-1,km))
     + +KT(3,n)*S3*(2.*vm1(i,j,km)+vm1(i,j-1,km)+vm1(i-1,j  ,km))
     + +KT(4,n)*S4*(2.*vm1(i,j,km)+vm1(i-1,j,km)+vm1(i-1,j+1,km))
     + +KT(5,n)*S5*(2.*vm1(i,j,km)+vm1(i-1,j+1,km)+vm1(i,j+1,km))
     + +KT(6,n)*S6*(2.*vm1(i,j,km)+vm1(i,j+1,km)+vm1(i+1,j  ,km)))

      FZV=0.5*(Sum1+Sum2)*W(i,j,K)

*     Metric in Transport
	Transp_Metric_U= -Um2(i,j,k)*Vm1(i,j,k)*COS_A/(R*S0)
	Transp_Metric_V= +Um2(i,j,k)*Um1(i,j,k)*COS_A/(R*S0)

      END IF

      DIFU = A*DIFU
      DIFV = A*DIFV

      FXUK=COEFX*FXUK
      FYUK=COEFY*FYUK
      FXVK=COEFX*FXVK
      FYVK=COEFY*FYVK
      FXUKm=COEFX*FXUKm
      FYUKm=COEFY*FYUKm
      FXVKm=COEFX*FXVKm
      FYVKm=COEFY*FYVKm

      FXU=c6*(2.*FXUK+FXUKm)*HZ(Km)
      FYU=c6*(2.*FYUK+FYUKm)*HZ(Km)
      FXV=c6*(2.*FXVK+FXVKm)*HZ(Km)
      FYV=c6*(2.*FYVK+FYVKm)*HZ(Km)

      FU= .5*DIFU*HZ(Km)  +FZU-FXU-FYU
      FV= .5*DIFV*HZ(Km)  +FZV-FXV-FYV

      U(I,J,K)=Um2(I,J,K)+CA*FU 
     #          + 0.5*REAL(m)*dt*(Transp_Metric_U  +S_Metric_Diff_U)
      V(I,J,K)=Vm2(I,J,K)+CA*FV 
     #          + 0.5*REAL(m)*dt*(Transp_Metric_V  +S_Metric_Diff_V)
      END IF
1     CONTINUE

      RETURN
      END

      SUBROUTINE DIFF(T,i,j,k,N,NP,DIFT,DIFTP,KB,il,jl,kl,KT)

*     Standart Galerkin diffusion approximation for
*     constant coefficients.
*     Boundary conditions - Climate restoring term -
*     is taken into account in the Main Program.
*     Version 30.06.2015.

      dimension T(0:il,0:jl,kl), KT(6,13)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      real KT
      DIFT=0.
      DIFTP=0.
      IF( K .GT. 1) THEN
C     UPPER HALF
      DIFT= DIFT+KT(1,n)*(T(I+1,J,K)-T(I,J,K))/S1
     *       -KT(2,n)*S2*(T(I,J,K)-T(I,J-1,K))*ASR2
      DIFT= DIFT+KT(3,n)*((T(I-1,J,K)-T(I,J,K))/S3
     *               -S3*(T(I,J,K)-T(I,J-1,K))*ASR2)
      DIFT= DIFT+KT(4,n)*(T(I-1,J,K)-T(I,J,K))/S4
     *       -KT(5,n)*S5*(T(I,J,K)-T(I,J+1,K))*ASR2
      DIFT= DIFT+KT(6,n)*((T(I+1,J,K)-T(I,J,K))/S6
     *               -S6*(T(I,J,K)-T(I,J+1,K))*ASR2)
      END IF
      IF( K .LT. KB) THEN
C     LOWER HALF
      DIFTP= DIFTP+KT(1,np)*(T(I+1,J,K)-T(I,J,K))/S1
     *         -KT(2,np)*S2*(T(I,J,K)-T(I,J-1,K))*ASR2
      DIFTP= DIFTP+KT(3,np)*((T(I-1,J,K)-T(I,J,K))/S3
     *               -S3*(T(I,J,K)-T(I,J-1,K))*ASR2)
      DIFTP= DIFTP+KT(4,np)*(T(I-1,J,K)-T(I,J,K))/S4
     *         -KT(5,np)*S5*(T(I,J,K)-T(I,J+1,K))*ASR2
      DIFTP= DIFTP+KT(6,np)*((T(I+1,J,K)-T(I,J,K))/S6
     *               -S6*(T(I,J,K)-T(I,J+1,K))*ASR2)
      END IF
      RETURN
      END

      SUBROUTINE F_Metric(u,v,i,j,k,N,NP,DIFU,DIFUP,DIFV,DIFVP,KB,
     *                                       il,jl,kl,KT,COS_A,hx)
*     First metric terms for momentum viscosity.
*     Version 30.06.2015.
      dimension u(0:il,0:jl,kl), v(0:il,0:jl,kl), KT(6,13)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      real KT

      IF( K .GT. 1) THEN
C     UPPER HALF
      DIFU=DIFU-2.*COS_A*hx*(KT(1,n)*(v(i+1,j  ,k)-v(i,j  ,k))/S1
     *                      +KT(2,n)*(v(i+1,j-1,k)-v(i,j-1,k))/S2
     *                      +KT(3,n)*(v(i  ,j  ,k)-v(i-1,j,k))/S3
     *                      +KT(4,n)*(v(i,  j,  k)-v(i-1,j,k))/S4
     *                      +KT(5,n)*(v(i,j+1,k)-v(i-1,j+1,k))/S5
     *                      +KT(6,n)*(v(i+1,j,k)-v(i  ,j  ,k))/S6 )

      DIFV=DIFV+2.*COS_A*hx*(KT(1,n)*(u(i+1,j  ,k)-u(i,j  ,k))/S1
     *                      +KT(2,n)*(u(i+1,j-1,k)-u(i,j-1,k))/S2
     *                      +KT(3,n)*(u(i  ,j  ,k)-u(i-1,j,k))/S3
     *                      +KT(4,n)*(u(i  ,j  ,k)-u(i-1,j,k))/S4
     *                      +KT(5,n)*(u(i,j+1,k)-u(i-1,j+1,k))/S5
     *                      +KT(6,n)*(u(i+1,j,k)-u(i  ,j  ,k))/S6 )
      END IF
      IF( K .LT. KB) THEN
C     LOWER HALF
      DIFUP=DIFUP-2.*COS_A*hx*(KT(1,np)*(v(i+1,j,k)-v(i,j,k))/S1
     *                      +KT(2,np)*(v(i+1,j-1,k)-v(i,j-1,k))/S2
     *                      +KT(3,np)*(v(i  ,j  ,k)-v(i-1,j,k))/S3
     *                      +KT(4,np)*(v(i  ,j  ,k)-v(i-1,j,k))/S4
     *                      +KT(5,np)*(v(i,j+1,k)-v(i-1,j+1,k))/S5
     *                      +KT(6,np)*(v(i+1,j,k)-v(i  ,j  ,k))/S6 )

      DIFVP=DIFVP+2.*COS_A*hx*(KT(1,np)*(u(i+1,j,k)-u(i,j,k))/S1
     *                      +KT(2,np)*(u(i+1,j-1,k)-u(i,j-1,k))/S2
     *                      +KT(3,np)*(u(i  ,j  ,k)-u(i-1,j,k))/S3
     *                      +KT(4,np)*(u(i  ,j  ,k)-u(i-1,j,k))/S4
     *                      +KT(5,np)*(u(i,j+1,k)-u(i-1,j+1,k))/S5
     *                      +KT(6,np)*(u(i+1,j,k)-u(i  ,j  ,k))/S6 )
      END IF


      RETURN
	END

      subroutine Vadapt
c-------------------------------------------------------
c     Velocity component caused by baroclinic pressure
c     and Coriolis force and explicit level gradients.
c     Crank-Nikolson time scheme.
c     Version 21.05.2015.
c-------------------------------------------------------

      INCLUDE 'Slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

	c6= 1./6.

      HXX= 1./HX
      HYY= 1./HY
      Ax= G*Hxx
      Ay= G*Hyy
      asr=hx/hy

      do j=1,jl
      s0=si(j)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=(2.*S0+SM)/3.
      S2=(2.*SM+S0)/3.
      S3=S1
      S4=(2.*S0+SP)/3.
      S5=(2.*SP+S0)/3.
      S6=S4

      Rs0= r*s0

      do i=1,il
      cor=-2.*co(i)*s0*Om
      CB= 0.5*cor*DT
      CDc=1./(1.+cb**2)
      kb=km2(i,j)

      IF( kb .GT. 0) then


      do k=1,kb-1

      IF( k .EQ. 1) THEN
      hzk1= 0.
      ELSE
      hzk1= hz(k-1)
      END IF

      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA= 12.*DT/(Rs0*(cg(n)*hzk1+cg(np)*hz(k)))
      CD= CA/(1.+CB**2)

*     Coriolis
	fu= um1(i,j,k)+cb*vm1(i,j,k)
	fv= vm1(i,j,k)-cb*um1(i,j,k)
      u(i,j,k)= CDc*(fu+cb*fv)
      v(i,j,k)= CDc*(fv-cb*fu)

      call gradcal(i,j,k,N ,derx ,dery ,PBclin,il1,jl1,kl,KT)
      call gradcal(i,j,k,NP,derxp,deryp,PBclin,il1,jl1,kl,KT)

      if(k.GT.1) then
      call gradcal(i,j,k-1,N ,derxm ,derym ,PBclin,il1,jl1,kl,KT)
	else
	derxm=0.
	derym=0.
	end if

      call gradcal(i,j,k+1,NP,derxpp,derypp,PBclin,il1,jl1,kl,KT)

      fu= c6*Hxx*( hzk1*(2.*derx+derxm) +hz(k)*(2.*derxp+derxpp))
      fv= c6*Hyy*( hzk1*(2.*dery+derym) +hz(k)*(2.*deryp+derypp))
      
      
c      fu= 0.5*Hxx*( hzk1*derx +hz(k)*derxp)
c      fv= 0.5*Hyy*( hzk1*dery +hz(k)*deryp)

*     Explicit barotropic pressure
      call graddz(i,j,N ,derx ,dery ,dz,ilp,jlp,KT)
      call graddz(i,j,NP,derxp,deryp,dz,ilp,jlp,KT)
      fu= fu+ 0.25*Ax*( hzk1*derx +hz(k)*derxp)
      fv= fv+ 0.25*Ay*( hzk1*dery +hz(k)*deryp)

      U(I,J,K)=U(I,J,K)+CD*(FU+CB*FV)
      V(I,J,K)=V(I,J,K)+CD*(FV-CB*FU)

      end do

      CA= 12.*DT/(Rs0*(cg(np)*hz(kb-1)))
      CD= CA/(1.+CB**2)

*     Coriolis
	fu= um1(i,j,kb)+cb*vm1(i,j,kb)
	fv= vm1(i,j,kb)-cb*um1(i,j,kb)
      u(i,j,kb)= CDc*(fu+cb*fv)
      v(i,j,kb)= CDc*(fv-cb*fu)

      call gradcal(i,j,kb  ,NP,derx ,dery ,PBclin,il1,jl1,kl,KT)
      call gradcal(i,j,kb-1,NP,derxm,derym,PBclin,il1,jl1,kl,KT)

      fu= c6*Hxx*hz(kb-1)*(2.*derx+derxm)
      fv= c6*Hyy*hz(kb-1)*(2.*dery+derym)
      
c      fu= 0.5*Hxx* hz(kb-1)*derx
c      fv= 0.5*Hyy* hz(kb-1)*dery
      

*     Explicit barotropic pressure
      call graddz(i,j,NP ,derx ,dery ,dz,ilp,jlp,KT)
      fu= fu+ 0.25*Ax*hz(kb-1)*derx 
      fv= fv+ 0.25*Ay*hz(kb-1)*dery 

      U(I,J,Kb)=U(I,J,Kb)+CD*(FU+CB*FV)
      V(I,J,Kb)=V(I,J,Kb)+CD*(FV-CB*FU)

      end if

      end do
      end do

      return
      end

      subroutine gradcal(i,j,k,N,derx,dery,P,il1,jl1,kl,KT)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      dimension P(0:il1,0:jl1,kl), KT(6,13)
      real KT
      derx=-(   KT(1,n)*(P(i+1,j  ,k)-P(i  ,j  ,k))+
     +          KT(6,n)*(P(i+1,j  ,k)-P(i  ,j  ,k))+
     +          KT(3,n)*(P(i  ,j  ,k)-P(i-1,j  ,k))+
     +          KT(4,n)*(P(i  ,j  ,k)-P(i-1,j  ,k))+
     +          KT(2,n)*(P(i+1,j-1,k)-P(i  ,j-1,k))+
     +          KT(5,n)*(P(i  ,j+1,k)-P(i-1,j+1,k))) /6.
      dery=-(s1*KT(1,n)*(P(i+1,j  ,k)-P(i+1,j-1,k))+
     +       s2*KT(2,n)*(P(i  ,j  ,k)-P(i  ,j-1,k))+
     +       s3*KT(3,n)*(P(i  ,j  ,k)-P(i  ,j-1,k))+
     +       s4*KT(4,n)*(P(i-1,j+1,k)-P(i-1,j  ,k))+
     +       s5*KT(5,n)*(P(i  ,j+1,k)-P(i  ,j  ,k))+
     +       s6*KT(6,n)*(P(i  ,j+1,k)-P(i  ,j  ,k))) /6.
      return
      end

      SUBROUTINE DENBPR(m)

*     Version 03.06.2015.
*
*     Water Pressure (without sea level).
*     Atmosphere pressure is taken into account if M=2.
*     If M=1 - ordinary assumtion Pbclin(z=0)=0.
*     Pressure in [Pa] -> 10.*Pa(i,j) in [dyn/cm2].
*     Snow Pressure Added.
*     Result is in the PBclin.


      INCLUDE 'Slo2.fi'
	Include 'Tparm.fi'

      DO j=1,jl
      DO i=1,il
      KK2=KM2(i,j)
      IF( KK2 .GT. 0) THEN


	IF(m.EQ.2) then ! Inverted Barometer Effect
	 PBclin(i,j,1)=10.*PA(i,j)
	ELSE
	 PBclin(i,j,1)=0.
	END IF ! m

*     Pressure by snow and ice.
      hs= 0.
      hi= 0.
      do mg=1,mgrad
	hs=hs+rosdry*hsnow(mg,i,j) ! snow mass
	hi=hi+roi   *hice (mg,i,j) ! ice mass
	end do


      PBclin(i,j,1)= PBclin(i,j,1) +g*(hs+hi) !-g*PME_dz(i,j)

*     Integration by Z
      DO k=2,KK2
      PBclin(i,j,k)=PBclin(i,j,k-1) +
     &              0.5*g*hz(k-1)*(ro(i,j,k-1)+ro(i,j,k))/row
      end do

	end if ! km2>0
      end do
      end do

      return
      end

      function sigma_t(t,s,p)

	real*8 pp,t8,s8,t2,sqrts,anum,aden,pt,rho_from_theta

*     version 14.06.2013

C     SEA WATER DENSITY DEVIATION FROM 1.02 [GR/CM**3]
C     AS FUNCTION OF T[grad C](potential),S[PPT],P[MPa]
C     By D. Brydon et al.:"A new approximation of the equation of state for
C                   seawater, suitable for numerical ocean models."
C     In: J.Geophys.Res.,v.104,No.C1, p.1537-1540, 1999.

C     VARIANT 2: -2<T<40;0<S<42;0<P<100.

*      P= P -0.101325 !! if atmospheric pressure is taken into account


*      P2= P*P
*	T2= T*T
*      sigma_t =-2.0092E-02  +5.07043E-04*P-5.43283E-07*P2
*     # + ( 5.10768E-05-3.69119E-06*P+6.54837E-09*P2)*T
*     # + ( 8.05999E-04-9.34012E-07*P+1.38777E-09*P2)*S
*     # + (-7.40849E-06+5.33243E-08*P-1.01563E-10*P2)*T2
*     # + (-3.01036E-06+1.75145E-08*P-2.34892E-11*P2)*T*S
*     # + ( 3.32267E-08-3.25887E-10*P+4.98612E-13*P2)*T2*T
*     # + ( 3.21931E-08-1.65849E-10*P+2.17612E-13*P2)*T2*S

******* ANOTHER POSSIBILITY

!     in-situ density from potential temperature, as in 
!     Jackett, D. R., McDougall, T. J., Feistel, R., Wright, 
!     D. G., and Griffies, S. M.: 
!     Algorithms for density, potential temperature, conservative
!     temperature, and freezing temperature of seawater, 
!     Journal of Atmospheric and Oceanic Technology, 23, 1709–1728, 2006.
!
!    s                : salinity                           (psu)
!    t                : potential temperature              (deg C, ITS-90)
!    p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!
!    rho_from_theta   : in-situ density                    (kg m^-3)
!
!    check value      : rho_from_theta(20,20,1000) = 1017.728868019642

      pp=DBLE(p) !!-10.1325d0 !!if atm. pressure is taken into account

      t8 = dble(t)
	s8 = dble(s)
      t2 = t8**2
	sqrts = sqrt(s8)

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


      rho_from_theta = anum/aden

	sigma_t= REAL(rho_from_theta -1.02d3)

      return
      end

      real function CDN(W,media)

c     Version 13.06.13.
c     Air-water/ice drag coefficient. [W]=m/s.
c     Media=1 - water and media=2 - ice.
      if(media.EQ.1) then

c      Cdn= (1.1+0.04*W)*1.e-3 ! AOMIP Drag. Scheme 1.
 
*Large, W. G., and S. Pond, 1981: Open ocean momentum flux measurements
*in moderate to strong winds. J. Phys. Oceanogr., 11, 324-336.

	Cdn= 2.7e-3 /max(0.5, W) +1.42e-4 + 7.64e-5 *W

c     My old parameterization based on Bowden
c      IF( w.LE.5.) Cdn=1.11e-3
c      IF( w.GT.5. .AND. w.LE.10.)
c     *  Cdn=1.11e-3 +(w-5.)*0.2*(1.45-1.11)*1.e-3
c      IF( w.GT.10. .AND. w.LE.15.)
c     *  Cdn=1.45e-3 +(w-10.)*0.2*(1.77-1.45)*1.e-3
c      IF( w.GT.15. .AND. w.LE.20.)
c     *  Cdn=1.77e-3 +(w-15.)*0.2*(2.07-1.77)*1.e-3
c      IF( w.GT.20. .AND. w.LE.25.)
c     *  Cdn=2.07e-3 +(w-20.)*0.2*(2.36-2.07)*1.e-3
c      IF( w.GT.25. .AND. w.LE.30.)
c     *  Cdn=2.36e-3 +(w-25.)*0.2*(2.65-2.36)*1.e-3
c      IF( w.GT.30. .AND. w.LE.35.)
c     *  Cdn=2.65e-3 +(w-30.)*0.2*(2.95-2.65)*1.e-3
c      IF( w.GT.35. .AND. w.LE.40.)
c     *  Cdn=2.95e-3 +(w-35.)*0.2*(3.25-2.95)*1.e-3
c      IF( w.GT.40.) Cdn=3.25e-3
 
	end if

	if(media.EQ.2) then

cc      Cdn= 7.00 e-3           ! Canadian archipelago estimates
                              ! Steiner, 1999
cc	Cdn= 2.75 e-3           ! SIMIP Results for ice.
cc      Cdn= 3.00 e-3           ! Steiner, 2001.
cc	Cdn= 1.925e-3           ! AWI model - 0.35*CDw
cc	Cdn= 1.605e-3           ! CAM5
cc	Cdn= 1.895e-3           ! ECHAM5
      Cdn= 2.00 e-3         ! Overland, 1985 (see Zhang 2003)
**	Cdn= 3.00 e-3           ! Overland, 1985 (see Leepparanta)

cc      Cdn= (1.1+0.04*W)*1.e-3 ! AOMIP Drag.

	end if

	if(media.EQ.3) then
      Cdn=(1.-EXP(-2.e-3*W**2))/MAX(1.,W**2) ! Andrey's Scheme 2.
	end if

      return
      end
