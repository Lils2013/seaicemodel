      subroutine ice_inertion(M)
c-----------------------------------------------------------
c     Horizontal diffusion and 2D transport of the ice momentum.
c     Modified Euler Scheme (Swansea). M - number of step (1,2).
c     Metric terms in spherical coordinates.
c                Version 14.02.2016.
c-----------------------------------------------------------
      INCLUDE 'Slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

	A_ice= 5.e7  ! Numerical Ice Viscosity to suppress oscillations

	c6=1./6.

      HXX= 1./HX
      COEFX=R*HXX/18.
      HYY= 1./HY
      COEFY=R*HYY/18.
      R2=R*R
      asr=hx/hy
	asr2=asr*asr
      A=.5*A_ice*HXX*HXX  ! For Standart Galerkin Scheme.

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
      r2s06= c6*r2*s0

	COS_A= 0.5*(SP-SM)/hy
	CTG2 = (COS_A/S0)**2

                           DO 1 I=1,IL

      ni=nt3(i,j,1)

ccc      if(ice_mask(i,j).EQ.1) then 
ccc      if(ice_mask(i,j).EQ.1 .and. n.LE.1 .and. n.NE.0) then 
      if(ice_mask(i,j).EQ.1 .and. ni.EQ.1) then ! Internal points

      n= ABS(ni)
      CA= 0.5*REAL(M)*DT/(R2S06*CG(n))

      FXUK =0.
      FYUK =0.
      FXVK =0.
      FYVK =0.



      if(M.EQ.1) then

	DIFU=0.
	DIFV=0.

      S_Metric_Diff_U= 0.
      S_Metric_Diff_V= 0.

      call 
     *transp_momentum_ice(FXUK,FYUK,i,j,n,uice2,vice2,uice2,KT,il1,jl1)
      call 
     *transp_momentum_ice(FXVK,FYVK,i,j,n,uice2,vice2,vice2,KT,il1,jl1)

*     Metric in Transport
	Transp_Metric_U= -Uice2(i,j)*Vice2(i,j)*COS_A/(R*S0)
	Transp_Metric_V= +Uice2(i,j)*Uice2(i,j)*COS_A/(R*S0)

      else   ! M=2

      CALL DIFF_IceMomentum(uice2,I,J,N,DIFU,il1,jl1,KT)
      CALL DIFF_IceMomentum(vice2,I,J,N,DIFV,il1,jl1,KT)

*     First metric
      CALL F_Metric_Ice(uice2,vice2,i,j,N,DIFU,DIFV,il1,jl1,KT,COS_A,hx)

*     Second Metric
      S_Metric_Diff_U= A_ice*(1.-CTG2)*Uice2(i,j)/R2
      S_Metric_Diff_V= A_ice*(1.-CTG2)*Vice2(i,j)/R2

      call 
     *transp_momentum_ice(FXUK,FYUK,i,j,n,uice2,vice2,uice1,KT,il1,jl1)
      call 
     *transp_momentum_ice(FXVK,FYVK,i,j,n,uice2,vice2,vice1,KT,il1,jl1)

*     Metric in Transport
	Transp_Metric_U= -Uice2(i,j)*Vice1(i,j)*COS_A/(R*S0)
	Transp_Metric_V= +Uice2(i,j)*Uice1(i,j)*COS_A/(R*S0)

      end if ! M=1,2

      FXU = COEFX*FXUK
      FYU = COEFY*FYUK
      FXV = COEFX*FXVK
      FYV = COEFY*FYVK

      DIFU = A*DIFU
      DIFV = A*DIFV

      FU=-FXU-FYU+DIFU
      FV=-FXV-FYV+DIFV

      Uice(I,J)=Uice2(I,J)+CA*FU + 
     *          0.5*real(m)*dt*(Transp_Metric_U+S_Metric_Diff_U)
      Vice(I,J)=Vice2(I,J)+CA*FV + 
     *          0.5*real(m)*dt*(Transp_Metric_V+S_Metric_Diff_V)

	end if  !  Ice motion area ice_mask=1

1     CONTINUE

      RETURN
      END


      subroutine transp_momentum_ice(FXTK,FYTK,i,j,n,u,v,t,KT,il,jl)
      dimension u(0:il,0:jl),v(0:il,0:jl),t(0:il,0:jl),KT(6,13)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
*     Version 14.02.2016.

*---------------     Triangle # 1.  ---------------------------------
      FXTK=FXTK+ 0.75*
     *((2.*U(I,J)+U(I+1,J)+U(I+1,J-1))*T(i  ,J  )
     *   +(U(I,J)+2.*U(I+1,J)+U(I+1,J-1))*T(i+1,J  )
     *   +(U(I,J)+U(I+1,J)+2.*U(I+1,J-1))*T(i+1,J-1)  )
     *  *KT(1,n)
*---------------     Triangle # 2.  -----------------------------------
      FYTK=FYTK-0.75*S2*
     *((2.*V(I,J)+V(I,J-1)+V(I+1,J-1))*T(i  ,J  )
     *   +(V(I,J)+2.*V(I,J-1)+V(I+1,J-1))*T(i  ,J-1)
     *   +(V(I,J)+V(I,J-1)+2.*V(I+1,J-1))*T(i+1,J-1)  )
     *  *KT(2,n)
*---------------     Triangle # 3.  -------------------------------
      FXTK=FXTK-0.75*
     *((2.*U(I,J)+U(I-1,J)+U(I,J-1))*T(i  ,J  )
     *   +(U(I,J)+2.*U(I-1,J)+U(I,J-1))*T(i-1,J)
     *   +(U(I,J)+U(I-1,J)+2.*U(I,J-1))*T(i  ,J-1)  )
     *  *KT(3,n)
     
      FYTK=FYTK-0.75*S3*
     *((2.*V(I,J)+V(I,J-1)+V(I-1,J))*T(i  ,J  )
     *   +(V(I,J)+2.*V(I,J-1)+V(I-1,J))*T(i  ,J-1)
     *   +(V(I,J)+V(I,J-1)+2.*V(I-1,J))*T(i-1,J  )  )
     *  *KT(3,n)
*---------------     Triangle # 4.  ---------------------------------
      FXTK=FXTK-0.75*
     *((2.*U(I,J)+U(I-1,J)+U(I-1,J+1))*T(i  ,J  )
     *   +(U(I,J)+2.*U(I-1,J)+U(I-1,J+1))*T(i-1,J  )
     *   +(U(I,J)+U(I-1,J)+2.*U(I-1,J+1))*T(i-1,J+1)  )
     *  *KT(4,n)
*---------------     Triangle # 5.  ----------------------------------
      FYTK=FYTK+0.75*S5*
     *((2.*V(I,J)+V(I,J+1)+V(I-1,J+1))*T(i  ,J  )
     *   +(V(I,J)+2.*V(I,J+1)+V(I-1,J+1))*T(i  ,J+1)
     *   +(V(I,J)+V(I,J+1)+2.*V(I-1,J+1))*T(i-1,J+1)  )
     *  *KT(5,n)
*---------------     Triangle # 6.  -----------------------------------
      FXTK=FXTK+0.75*
     *( (2.*U(I,J)+U(I+1,J)+U(I,J+1))*T(i  ,J   )
     *  +(U(I,J)+2.*U(I+1,J)+U(I,J+1))*T(i+1,J  )
     *  +(U(I,J)+U(I+1,J)+2.*U(I,J+1))*T(i  ,J+1)  )
     *  *KT(6,n)
     
      FYTK=FYTK+0.75*S6*
     *((2.*V(I,J)+V(I,J+1)+V(I+1,J))*T(i  ,J  )
     * +(V(I,J)+2.*V(I,J+1)+V(I+1,J))*T(i  ,J+1)
     * +(V(I,J)+V(I,J+1)+2.*V(I+1,J))*T(i+1,J)   )
     *  *KT(6,n)

      return
      end

      SUBROUTINE F_Metric_Ice(u,v,i,j,N,DIFU,DIFV,il,jl,KT,COS_A,hx)
*     First metric terms for ice momentum viscosity.
*     Version 21.01.2007.
      dimension u(0:il,0:jl), v(0:il,0:jl), KT(6,13)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      real KT

      DIFU=DIFU+2.*COS_A*hx*(KT(1,n)*(v(i+1,j)-v(i,j))/S1
     *                      +KT(2,n)*(v(i+1,j-1)-v(i,j-1))/S2
     *                      +KT(3,n)*(v(i,j)-v(i-1,j))/S3
     *                      +KT(4,n)*(v(i,j)-v(i-1,j))/S4
     *                      +KT(5,n)*(v(i,j+1)-v(i-1,j+1))/S5
     *                      +KT(6,n)*(v(i+1,j)-v(i,j))/S6)

      DIFV=DIFV-2.*COS_A*hx*(KT(1,n)*(u(i+1,j)-u(i,j))/S1
     *                      +KT(2,n)*(u(i+1,j-1)-u(i,j-1))/S2
     *                      +KT(3,n)*(u(i,j)-u(i-1,j))/S3
     *                      +KT(4,n)*(u(i,j)-u(i-1,j))/S4
     *                      +KT(5,n)*(u(i,j+1)-u(i-1,j+1))/S5
     *                      +KT(6,n)*(u(i+1,j)-u(i,j))/S6)
      RETURN
	END

      SUBROUTINE DIFF_IceMomentum(T,i,j,N,DIFT,il,jl,KT)

*     Standart Galerkin diffusion approximation for
*     constant coefficients.
*     Version 02.02.2016.

      dimension T(0:il,0:jl), KT(6,13)
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      real KT
      DIFT=0.
      DIFT= DIFT+KT(1,n)*(T(I+1,J)-T(I,J))/S1
     *       -S2*KT(2,n)*(T(I,J)-T(I,J-1))*ASR2
      DIFT= DIFT+KT(3,n)*((T(I-1,J)-T(I,J))/S3
     *               -S3*(T(I,J)-T(I,J-1))*ASR2)
      DIFT= DIFT+KT(4,n)*(T(I-1,J)-T(I,J))/S4
     *       -S5*KT(5,n)*(T(I,J)-T(I,J+1))*ASR2
      DIFT= DIFT+KT(6,n)*((T(I+1,J)-T(I,J))/S6
     *               -S6*(T(I,J)-T(I,J+1))*ASR2)
      RETURN
      END
