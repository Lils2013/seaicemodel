      SUBROUTINE Bio_Vert_Gravity (wg,Sclr,QScalar)
      
* vertival turbulent diffusion and vertival gravity deposition of the bio-scalar
* N. Iakovlev for I. Chernov's version of BFM      

c     Version 13.03.2015 Black Friday!

      INCLUDE 'Slo2.fi'
      REAL AM(KL), BM(KL), CM(KL), FM(KL), RKSI(KL)
      Dimension Sclr(0:il1,0:jl1,kl), Qscalar(il,jl), Wg(il,jl,kl)

	c3=1./3.
	c12=0.25*c3
	
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

      KB = KM2(I,J)
	n=ABS(nt3(i,j,1))
	
*     -----------------------------------------------

      IF( KB .GT. 0) THEN

*     -------------- Ocean Surface ------------------

      CA=2.*DT/hz(1)
      
*     ---------- Vertical turbulent diffusion -------

      AM(1)= 0.
      CM(1)= -CA*AZT(i,j,1)/hz(1)
      BM(1)= 1. -CM(1) 
      FM(1)= Sclr(i,j,1)+CA*QScalar(i,j)
      
*     ---------- Vertical gravity deposition --------
     
      CM(1)= CM(1) +0.5*CA*S0*Wg(i,j,1)
      BM(1)= BM(1) +0.5*CA*S0*Wg(i,j,1)
      
*     ------------------- Deep water ----------------
      DO K=2, KB -1
      n= ABS(nt3(i,j,k))
      np=ABS(nt3(i,j,k+1))

      CA=2.*DT/(CG(n)*HZ(K-1)+CG(np)*HZ(K))

*     ---------- Vertical turbulent diffusion -------
      AM(K)= -CA*CG( n)*AZT(i,j,k-1)/HZ(K-1)
      CM(K)= -CA*CG(np)*AZT(i,j,k  )/HZ(k)
      BM(K)= 1. -AM(K) -CM(K)
      FM(K)= Sclr(i,j,k) 
      
*     ---------- Vertical gravity deposition --------

      AM(K)= AM(K) -0.5*S0*CA*cg(n )*Wg(i,j,k-1)
      CM(K)= CM(K) +0.5*S0*CA*cg(np)*Wg(i,j,k)
      BM(K)= BM(K) +0.5*S0*CA*(cg(np)*Wg(i,j,k)-cg(n)*Wg(i,j,k-1))
      
      END DO
      
*     ----------------- Ocean bottom ----------------

      CA=2.*DT
      Q=1./HZ(KB-1)

      AM(KB)= -CA*Q*Q*AZT(i,j,kb-1)
      BM(KB)= 1. -AM(KB)
      CM(KB)= 0.
      FM(KB)= Sclr(i,j,kb)
      
*     ---------- Vertical gravity deposition --------

      AM(KB)= AM(KB) -0.5*CA*Q*S0*Wg(i,j,kb-1)
      BM(KB)= BM(KB) +0.5*CA*Q*S0*Wg(i,j,kb-1)
      

*     3 Diagonal matrix factorization.
*     Vertical temperature distribution.
      CALL FACTOR(KL,AM,BM,CM,FM,RKSI,1,KB)
      
      Do k=1,kb
      Sclr(i,j,k)= max(0.,RKSI(K))
      End do

      END IF  ! kb>0
      END DO  ! i
      END DO  ! j


      RETURN
      END