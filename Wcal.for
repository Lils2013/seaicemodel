      SUBROUTINE WCALS

*     W is piecewise constant
*     Version 24/06/15 for Flather b.c. at open boundary

      INCLUDE 'Slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
	Include 'Tparm.fi'

      AsR= Hx/Hy
	asr2=asr*asr
	c3=1./3.
	HxR=Hx*R

	c3=1./3.

      SM=0.
      SP=0.
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
	HxRS0=HxR*S0

      DO 1 I=1,IL
      L3= KM2(I,J)
      IF( L3 .GT. 0) THEN
     
*     Pressure by snow and ice   
	hi= 0.
      do mg=1,mgrad
	hi=hi+rosdry*hsnow(mg,i,j)/row ! normalized snow mass
	hi=hi+roi   *hice (mg,i,j)/row ! normalized ice mass
	end do
	
*     Atmospheric pressure	
*	hi= hi +10.*PA(i,j)/(g*row)
		
*     Near the bottom
      nk= nt3(i,j,L3)
	nb=1
      IF(nk.LT.0) nb=-1

      Call Diverg(abs(nk),nb,KT,i,j,L3,divk,u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
      Call Diverg(abs(nk),nb,KT,i,j,L3-1,divkm,u,v,ub,dz,dzext,hi,
     *            km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
	
	DivK=c3*(2.*divK+divkm)

      nk= abs(nk)
	geomk=
     * (KT(1,nk)*S1+KT(2,nk)*S2+KT(3,nk)*S3
     + +KT(4,nk)*S4+KT(5,nk)*S5+KT(6,nk)*S6)


      W(i,j,L3)= -0.5*Hz(L3-1)*DivK/(geomk*HxR)
      
*     Next Horizons
      DO K=L3,3,-1

      nkm1= ABS(nt3(i,j,k-1))
      nk  = ABS(nt3(i,j,k  ))
      
	geomk=
     * (KT(1,nk)*S1+KT(2,nk)*S2+KT(3,nk)*S3
     + +KT(4,nk)*S4+KT(5,nk)*S5+KT(6,nk)*S6)

	geomkm1=
     * (KT(1,nkm1)*S1+KT(2,nkm1)*S2+KT(3,nkm1)*S3
     + +KT(4,nkm1)*S4+KT(5,nkm1)*S5+KT(6,nkm1)*S6)
      

      call Diverg(nkm1,nb,KT,i,j,k-1,divk, u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
      call Diverg(nkm1,nb,KT,i,j,k-2,divkm, u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
	DivK=c3*(2.*divK+divkm)


      call Diverg(nk  ,nb,KT,i,j,k-1,divkp,u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
      call Diverg(nk  ,nb,KT,i,j,k,divkpp,u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
	DivKP=c3*(2.*divKP+divkpp)

      W(i,j,K-1)= geomk*W(i,j,K)/geomkm1
     *            -0.5*(Hz(K-2)*DivK+Hz(K-1)*DivKP)/
     /                                          (geomkm1*HxR)
     
      END DO  ! K 

*     Ocean Surface.
 
      nk=ABS(nt3(i,j,1))
      
	geomk=
     * (KT(1,nk)*S1+KT(2,nk)*S2+KT(3,nk)*S3
     + +KT(4,nk)*S4+KT(5,nk)*S5+KT(6,nk)*S6)
      

      call Diverg(nk,nb,KT,i,j,1,divkp, u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
      call Diverg(nk,nb,KT,i,j,2,divkpp, u,v,ub,dz,dzext,hi,km2,z,
     *            ilp,jlp,il1,jl1,il,jl,kl,klp,g)
	DivKP=c3*(2.*divKP+divkpp)

      W(i,j,1)= W(i,j,2)-0.5*Hz(1)*DivKP/(geomk*HxR)
	 
      END IF
1     CONTINUE
      RETURN
      END

      subroutine Diverg(m,nb,KT,i,j,k,divk,u,v,ub,dz,dzext,hi,km2,z,
     *                  ilp,jlp,il1,jl1,il,jl,kl,klp,g)

*     Version 20/8/10 for Flather b.c. at open boundary modified for ice

      double precision dz(-1:ilp,-1:jlp)
      dimension u(0:il1,0:jl1,kl),v(0:il1,0:jl1,kl)
      dimension ub(0:il1,0:jl1,kl)
      dimension KT(6,13)
 	dimension dzext(il,jl), km2(0:il1,0:jl1), z(klp)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

	coef= 1.0e0

      DivK=   - kt(1,m)*(U(I,J,K)+U(I+1,J,K)+U(I+1,J-1,K))
     +     + S2*kt(2,m)*(V(I,J,K)+V(I,J-1,K)+V(I+1,J-1,K))*AsR
     +        + kt(3,m)*(U(I,J,K)+U(I-1,J,K)+U(I  ,J-1,K))
     +     + S3*kt(3,m)*(V(I,J,K)+V(I-1,J,K)+V(I  ,J-1,K))*AsR
     +        + kt(4,m)*(U(I,J,K)+U(I-1,J,K)+U(I-1,J+1,K))
     -     - S5*kt(5,m)*(V(I,J,K)+V(I-1,J+1,K)+V(I,J+1,K))*AsR
     -        - kt(6,m)*(U(I,J,K)+U(I+1,J,K)+U(I  ,J+1,K))
     -     - S6*kt(6,m)*(V(I,J,K)+V(I+1,J,K)+V(I  ,J+1,K))*AsR

c     Liquid boundaries
      IF( nb .EQ. -1) THEN
      add= -coef*(dzext(i,j)-REAL(dz(i,j))-hi ) * SQRT(g/z(km2(i,j)))
cc      add= -coef*(-REAL(dz(i,j))-hi ) * SQRT(g/z(km2(i,j)))

      IF(m.EQ.6) DivK= DivK-6.*(Ub(i,j,K)+add)
      IF(m.EQ.7) DivK= DivK-6.*S0*AsR*(Ub(i,j,K)+add)
      IF(m.EQ.8) DivK= DivK-6.*(Ub(i,j,K)+add)
      IF(m.EQ.9) DivK= DivK-6.*S0*AsR*(Ub(i,j,K)+add)
      IF(m.GT.9) DivK= DivK-3.*(1.+S0*AsR)*(Ub(i,j,K)+add)
      END IF

      return
      end
