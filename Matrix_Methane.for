      Program DiagSLO_Methane

*     Version 10.02.2014 for Flather b.c. at open boundary
*     Matrix elements calculation for the sea level equation.
*     All of the 19 diagonals are calculated.
*     Hmean - Depth scale (in cm).
*     Version for Polar Cape with shifted coordinates.
*     Note Coriolis parameter!

      double precision pi, finord, alwest, hxgr,hygr, dthr, Hmean
      Parameter (il=35, jl=49, kl=40, ilp=il+2, jlp=jl+2, klp=kl+1)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (hxgr=1.d0, hygr=1.d0, finord=23.d0)
      Parameter (Pi=3.1415926d0, Ncond=1, alwest=-16.d0)
      Parameter (Hmean= 1.d5)

*     Rec Length in different Fortrans: 1 for Compaq and UNIX,
*     4 for Watcom and NDP.
	Parameter(LFort=1)

      real KT(6,13)
      Dimension  dz(-1:ilp,-1:jlp),km2(0:il1,0:jl1),
     *          nt3(0:il1,0:jl1,kl), hz(kl), z(klp), co(-1:ilp),
     *          si(-1:jlp), nt1(-1:ilp,-1:jlp,kl)
      Dimension u(7,kl),v(7,kl)
      Dimension km21(-1:ilp,-1:jlp)
      Dimension Diag(19,il,jl), cg(13)
      double precision diag,dz,u,v,si,co,dt,hx,hy,om,r,g,asr,z,hz,a
      double precision rhx,rrhx,rhy,rrhy,ax,ay,s0,sm,sp,s1,s2,s3,s4,
     *  s5,s6,f,dzmean,teta,schem, divk,divkp
      double precision sum, cg
      common /cg/ cg
      character*72 i2ddim, diagdat, nt3dat
      data i2ddim/'km2.dat'/ diagdat/'diag.dat'/ nt3dat/'nt.dat'/
      data cg/6.d0,4.d0,5.d0,4.d0,5.d0,4*3.d0,2.d0,1.d0,2.d0,1.d0/

      write(*,*)' Enter time step in Hours (Double Prec.)--->'
      read(*,*) dthr


cc      dthr= 1.d0
                         dt=3600.d0*dthr
                         hx=hxgr*pi/180.d0
                         hy=hygr*pi/180.d0
                         om=.729d-4
                         r=.637d9
                         g=980.d0
        AsR= hxgr/hygr

        z(1)= 0.d0
        z(2)= 2.5d0
        z(3)= 5.d0
        z(4)= 7.5d0
        z(5)= 10.d0
        z(6)= 15.d0
        z(7)= 20.d0
        z(8)= 27.5d0
        z(9)= 37.5d0
        z(10)=50.d0
        z(11)=65.d0
        z(12)=82.5d0
        z(13)=100.d0
        z(14)=125.d0
        z(15)=150.d0
        z(16)=175.d0
        z(17)=200.d0
        z(18)=225.d0
        z(19)=250.d0
        z(20)=275.d0
        z(21)=300.d0
        z(22)=325.d0
        z(23)=350.d0
        z(24)=375.d0
        z(25)=410.d0
        z(26)=450.d0
        z(27)=500.d0
        z(28)=600.d0
        z(29)=700.d0
        z(30)=850.d0
        z(31)=1050.d0
        z(32)=1250.d0
        z(33)=1500.d0
        z(34)=1750.d0
        z(35)=2100.d0
        z(36)=2500.d0
        z(37)=3000.d0
        z(38)=3500.d0
        z(39)=4000.d0
        z(40)=4500.d0
        z(41)= 9999999.d0



      do 33 k=1,kl+1
33    z(k)=100.d0*z(k)
      do 34 k=1,kl
34    hz(k)=(z(k+1)-z(k))
      do 52 j=-1,jlp
      teta=pi*.5d0 - finord*pi/180.d0+(j-1.d0)*hy
      a= DSIN(teta)
52    si(j)= a
      do 521 i=-1,ilp
      teta= alwest*pi/180.d0+(i-1.d0)*hx
      a= DCOS(teta)
521   co(i)= a

      n3= (il+2)*(jl+2)*kl*LFort
      nnn=(il+2)*(jl+2)*LFort
      nd= 38*il*jl*LFort
      open (unit=21,file=nt3dat,status='old',access='direct',
     *      form='unformatted',recl=n3)
      open (unit=49,file=i2ddim,status='old',access='direct',
     *      form='unformatted',recl=nnn)
      open (unit=50,file=diagdat,status='old',access='direct',
     *      form='unformatted',recl=nd)
      read (21,rec=1) nt3
      read (49,rec=1) km2
      close (21)
      close (49)


      Schem= 2.d0
      info= 1

      do 9 m=1,19
      do 9 i=1,il
      do 9 j=1,jl
9     diag(m,i,j)= 0.d0

      do 30 i=-1,ilp
      do 30 j=-1,jlp
      dz(i,j)= 0.d0
      km21(i,j)=0
      do 30 k=1,kl
30    nt1(i,j,k)= 0

      do 31 i=1,il
      do 31 j=1,jl
      km21(i,j)=km2(i,j)
      do 31 k=1,kl
31    nt1(i,j,k)= nt3(i,j,k)

      RHx = R*Hx
      RRHx= 1.d0/(RHx)
      Ax= G/Hx
      RHy = R*Hy
      RRHy= 1.d0/(RHy)
      Ay= G/Hy

      call ktform(KT)

      do 1 j=1,jL

         S0=Si(j)
         SM=Si(j-1)
         SP=Si(j+1)
                  S1=(2.d0*S0+SM)/3.d0
                  S2=(2.d0*SM+S0)/3.d0
                  S3=S1
                  S4=(2.d0*S0+SP)/3.d0
                  S5=(2.d0*SP+S0)/3.d0
                  S6=S4
      do 1 i=1,iL

      IF( abs(nt1(i,j,1)) .GT. 0) THEN

      DO 100 mdiag=1,19

      IF( mdiag .EQ. 1) dz(i  ,j-2)=1.d0
      IF( mdiag .EQ. 2) dz(i+1,j-2)=1.d0
      IF( mdiag .EQ. 3) dz(i+2,j-2)=1.d0
      IF( mdiag .EQ. 4) dz(i-1,j-1)=1.d0
      IF( mdiag .EQ. 5) dz(i  ,j-1)=1.d0
      IF( mdiag .EQ. 6) dz(i+1,j-1)=1.d0
      IF( mdiag .EQ. 7) dz(i+2,j-1)=1.d0
      IF( mdiag .EQ. 8) dz(i-2,j  )=1.d0
      IF( mdiag .EQ. 9) dz(i-1,j  )=1.d0
      IF( mdiag .EQ. 10) dz(i  ,j  )=1.d0
      IF( mdiag .EQ. 11) dz(i+1,j  )=1.d0
      IF( mdiag .EQ. 12) dz(i+2,j  )=1.d0
      IF( mdiag .EQ. 13) dz(i-2,j+1)=1.d0
      IF( mdiag .EQ. 14) dz(i-1,j+1)=1.d0
      IF( mdiag .EQ. 15) dz(i  ,j+1)=1.d0
      IF( mdiag .EQ. 16) dz(i+1,j+1)=1.d0
      IF( mdiag .EQ. 17) dz(i-2,j+2)=1.d0
      IF( mdiag .EQ. 18) dz(i-1,j+2)=1.d0
      IF( mdiag .EQ. 19) dz(i  ,j+2)=1.d0

*     Barotropic velocities.

           do 3 m=1,7
           do 3 k=1,kl
           u(m,k)= 0.d0
3          v(m,k)= 0.d0
      if( nt1(i,j-1,1) .NE. 0)
     * call barotrop(1,i,j-1,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      if( nt1(i+1,j-1,1) .NE. 0)
     * call barotrop(2,i+1,j-1,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      if( nt1(i-1,j,1) .NE. 0)
     * call barotrop(3,i-1,j,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      if( nt1(i,j,1) .NE. 0)
     * call barotrop(4,i,j,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      if( nt1(i+1,j,1) .NE. 0)
     * call barotrop(5,i+1,j,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      if( nt1(i-1,j+1,1) .NE. 0)
     * call barotrop(6,i-1,j+1,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      if( nt1(i,j+1,1) .NE. 0)
     * call barotrop(7,i,j+1,schem,nt1,km21,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)

*     Substitution of barotropic velocities
*           to continuity equation.

      nb= 1
      IF(nt1(i,j,1) .LT. 0) nb=-1
      N=abs(nt1(i,j,1))

      IF( Ncond .EQ. 1) THEN

	DzMean=  KT(1,n)*S1*(2.d0*dz(i,j)+dz(i+1,j)+dz(i+1,j-1))
     +        +KT(2,n)*S2*(2.d0*dz(i,j)+dz(i+1,j-1)+dz(i,j-1))
     +        +KT(3,n)*S3*(2.d0*dz(i,j)+dz(i,j-1)+dz(i-1,j))
     +        +KT(4,n)*S4*(2.d0*dz(i,j)+dz(i-1,j)+dz(i-1,j+1))
     +        +KT(5,n)*S5*(2.d0*dz(i,j)+dz(i,j+1)+dz(i-1,j+1))
     +        +KT(6,n)*S6*(2.d0*dz(i,j)+dz(i+1,j)+dz(i,j+1))

      F= DzMean/(24.d0*DT)

      ELSE
      F= 0.d0
      END IF

      kb=km21(i,j)
      DO 55 K=1,kb-1


      n= abs(nt1(i,j,k+1))
      call diverg
     *           (n,nb,i,j,k  ,divk ,u,v,dz,km21,z,ilp,jlp,kl,klp,g,
     *                                           s0,s2,s3,s5,s6,asr)

      call diverg
     *           (n,nb,i,j,k+1,divkp,u,v,dz,km21,z,ilp,jlp,kl,klp,g,
     *                                           s0,s2,s3,s5,s6,asr)

55    F= F - RRHx*hz(k)*(divk+divkp)/12.d0

      diag(mdiag,i,j)= F/Hmean

*     Level restoring

      IF( mdiag .EQ. 1) dz(i  ,j-2)=0.d0
      IF( mdiag .EQ. 2) dz(i+1,j-2)=0.d0
      IF( mdiag .EQ. 3) dz(i+2,j-2)=0.d0
      IF( mdiag .EQ. 4) dz(i-1,j-1)=0.d0
      IF( mdiag .EQ. 5) dz(i  ,j-1)=0.d0
      IF( mdiag .EQ. 6) dz(i+1,j-1)=0.d0
      IF( mdiag .EQ. 7) dz(i+2,j-1)=0.d0
      IF( mdiag .EQ. 8) dz(i-2,j  )=0.d0
      IF( mdiag .EQ. 9) dz(i-1,j  )=0.d0
      IF( mdiag .EQ. 10) dz(i  ,j  )=0.d0
      IF( mdiag .EQ. 11) dz(i+1,j  )=0.d0
      IF( mdiag .EQ. 12) dz(i+2,j  )=0.d0
      IF( mdiag .EQ. 13) dz(i-2,j+1)=0.d0
      IF( mdiag .EQ. 14) dz(i-1,j+1)=0.d0
      IF( mdiag .EQ. 15) dz(i  ,j+1)=0.d0
      IF( mdiag .EQ. 16) dz(i+1,j+1)=0.d0
      IF( mdiag .EQ. 17) dz(i-2,j+2)=0.d0
      IF( mdiag .EQ. 18) dz(i-1,j+2)=0.d0
      IF( mdiag .EQ. 19) dz(i  ,j+2)=0.d0

100   CONTINUE

c     Diagonal domination check.
      sum= DABS(diag(1,i,j))+DABS(diag(2,i,j))+DABS(diag(3,i,j))+
     +     DABS(diag(4,i,j))+DABS(diag(5,i,j))+DABS(diag(6,i,j))+
     +     DABS(diag(7,i,j))+DABS(diag(8,i,j))+DABS(diag(9,i,j))+
     +     DABS(diag(11,i,j))+DABS(diag(12,i,j))+DABS(diag(13,i,j))+
     +     DABS(diag(14,i,j))+DABS(diag(15,i,j))+DABS(diag(16,i,j))+
     +     DABS(diag(17,i,j))+DABS(diag(18,i,j))+DABS(diag(19,i,j))
      if( sum .GT. DABS(diag(10,i,j)) ) info= 0

      END IF
1     Continue

      if( info .eq. 1) then
      write(*,*)' Matrix with diagonal domination'
      else
      write(*,*)' Matrix with no diagonal domination'
      end if

      write(50,rec=1) diag
      close(50)
7777  STOP
      END

      subroutine barotrop(m,i,j,schem,nt,km2,u,v,dz,hz,ilp,jlp,
     *            kl,si,co,om,R,dt,hx,hy,g)
      dimension dz(-1:ilp,-1:jlp)
      common /cg/ cg(13)
      dimension nt(-1:ilp,-1:jlp,kl),hz(kl),
     *          km2(-1:ilp,-1:jlp),si(-1:jlp), co(-1:ilp),
     *          u(7,kl),v(7,kl)
      double precision fu,fv,hzk1,derx,dery,derxp,deryp,cg,
     *          s0,s1,s2,s3,s4,s5,s6, schem,hz,u,v,dz,si,co,r,dt,
     *          hx,hy, sm,sp, r2,rhx,rhy,rrhx,rrhy,hxx,hyy,
     *          ca,cb,cd,g,om,cor

      cor= -2.d0*om*co(i)*si(j)

      r2= r*r
      rhx=r*hx
      rrhx=1.d0/rhx
      rhy=r*hy
      rrhy=1.d0/rhy
      hxx=1.d0/hx
      hyy=1.d0/hy

      KB = KM2(I,J)

         S0=Si(j)
         SM=Si(j-1)
         SP=Si(j+1)
                  S1=(2.d0*S0+SM)/3.d0
                  S2=(2.d0*SM+S0)/3.d0
                  S3=S1
                  S4=(2.d0*S0+SP)/3.d0
                  S5=(2.d0*SP+S0)/3.d0
                  S6=S4

      DO 2 K=1,KB-1
      n=ABS(nt(i,j,k))
      np=ABS(nt(i,j,k+1))

      IF( k .EQ. 1) THEN
        hzk1= 0.d0
      ELSE
        hzk1= hz(k-1)
      END IF

      CA= 24.d0*DT/(schem*R*S0*(CG(n)*HZK1+CG(np)*HZ(K)))
      CB= cor*dt/schem
      CD= CA/(1.d0 +CB**2)

      call gradcal(i,j,n ,derx ,dery ,dz,ilp,jlp,s1,s2,s3,s4,s5,s6)
      call gradcal(i,j,np,derxp,deryp,dz,ilp,jlp,s1,s2,s3,s4,s5,s6)

      fu= .25d0*g*Hxx*( hzk1*derx +hz(k)*derxp)
      fv= .25d0*g*Hyy*( hzk1*dery +hz(k)*deryp)

      U(m,K)= CD*(FU +CB*FV)
      V(m,K)= CD*(FV -CB*FU)

2     continue

      CA= 24.d0*DT/(schem*R*CG(np)*S0)
      CB= dt*cor/schem
      CD= CA/(1.d0 +CB**2)

      fu= .25d0*g*Hxx*derxp
      fv= .25d0*g*Hyy*deryp

      U(m,KB)= CD*(FU +CB*FV)
      V(m,KB)= CD*(FV -CB*FU)

      return
      end

      subroutine gradcal(i,j,N,derx,dery,dz,ilp,jlp,s1,s2,s3,s4,s5,s6)
      dimension dz(-1:ilp,-1:jlp)
      double precision dz, derx,dery,s1,s2,s3,s4,s5,s6

      if (N.eq.1) then
      derx=-(2.d0*(dz(i+1,j  )-dz(i-1,j  ))+
     +          dz(i+1,j-1)-   dz(i  ,j-1)+
     +          dz(i  ,j+1)-   dz(i-1,j+1)) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s3*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))+
     +       s6*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      else
      if (N.eq.2) then
      derx=-(2.d0*(dz(i+1,j  )-dz(i  ,j  ))+
     +       dz(i  ,j  )-dz(i-1,j  )+
     +       dz(i+1,j-1)-dz(i  ,j-1)) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s3*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s6*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.3) then
      derx=-(2.d0*(dz(i+1,j  )-dz(i  ,j  ))+
     +           dz(i  ,j+1)-dz(i-1,j+1)+
     +           dz(i+1,j-1)-dz(i  ,j-1)+
     +           dz(i  ,j  )-dz(i-1,j  )) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))+
     +       s6*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.4) then
      derx=-(2.d0*(dz(i  ,j  )-dz(i-1,j  ))+
     +           dz(i+1,j  )-dz(i  ,j  )+
     +           dz(i  ,j+1)-dz(i-1,j+1)) /6.d0
      dery=-(s3*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))+
     +       s6*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.5) then
      derx=-(dz(i+1,j  )-dz(i  ,j  )+
     +       2.d0*(dz(i  ,j  )-dz(i-1,j  ))+
     +       dz(i+1,j-1)-dz(i  ,j-1)+
     +       dz(i  ,j+1)-dz(i-1,j+1)) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s3*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.6) then
      derx=-(2.d0*(dz(i+1,j  )-dz(i  ,j  ))+
     +           dz(i+1,j-1)-dz(i  ,j-1)) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s6*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.7) then
      derx=-(dz(i+1,j  )-dz(i  ,j  )+
     +       dz(i  ,j  )-dz(i-1,j  )+
     +       dz(i  ,j+1)-dz(i-1,j+1)) /6.d0
      dery=-(s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))+
     +       s6*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.8) then
      derx=-(2.d0*(dz(i  ,j  )-dz(i-1,j  ))+
     +           dz(i  ,j+1)-dz(i-1,j+1)) /6.d0
      dery=-(s3*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.9) then
      derx=-(dz(i+1,j  )-dz(i  ,j  )+
     +       dz(i  ,j  )-dz(i-1,j  )+
     +       dz(i+1,j-1)-dz(i  ,j-1)) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))+
     +       s3*(dz(i  ,j  )-dz(i  ,j-1))) /6.d0
      end if
      if (N.eq.10) then
      derx=-(dz(i+1,j  )-dz(i  ,j  )+
     +       dz(i+1,j-1)-dz(i  ,j-1)) /6.d0
      dery=-(s1*(dz(i+1,j  )-dz(i+1,j-1))+
     +       s2*(dz(i  ,j  )-dz(i  ,j-1))) /6.d0
      end if
      if (N.eq.11) then
      derx=-(dz(i  ,j  )-dz(i-1,j  )) /6.d0
      dery=-(s3*(dz(i  ,j  )-dz(i  ,j-1))) /6.d0
      end if
      if (N.eq.12) then
      derx=-(dz(i  ,j  )-dz(i-1,j  )+
     +       dz(i  ,j+1)-dz(i-1,j+1)) /6.d0
      dery=-(s4*(dz(i-1,j+1)-dz(i-1,j  ))+
     +       s5*(dz(i  ,j+1)-dz(i  ,j  ))) /6.d0
      end if
      if (N.eq.13) then
      derx=-(dz(i+1,j  )-dz(i  ,j  )) /6.d0
      dery=-s6*(dz(i  ,j+1)-dz(i  ,j  )) /6.d0
      end if
      end if
      return
      end

      subroutine Diverg
     *(m,nb,i,j,k,divk,u,v,dz,km21,z,ilp,jlp,kl,klp,g,
     *s0,s2,s3,s5,s6,asr)

*     Version 20/10/09 for Flather b.c. at open boundary

      dimension u(7,kl),v(7,kl)
	dimension km21(-1:ilp,-1:jlp), z(klp)
      dimension dz(-1:ilp,-1:jlp)
      double precision dz,u,v,divk,s0,s2,s3,s5,s6,asr,add,z,g,coef


      coef= 1.0d0

**	DZ_Trelax= 6.4e-8
***	DZ_Trelax= SQRT(g*z(km21(i,j))) * 1.e-7 !! Horizontal resolution


      IF( m .EQ. 1) THEN
         DivK= -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
     +         +U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
     +         +U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
     -         -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 2) THEN
         DivK= -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
     +         +U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
     -         -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 3) THEN
         DivK= -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
     +         +U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
     -         -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 4) THEN
         DivK=  U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
     +         +U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
     -         -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 5) THEN
         DivK= -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
     +         +U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
     +         +U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 6) THEN
         DivK=  -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
     -          -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 7) THEN
         DivK=   U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
     -          -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 8) THEN
         DivK=   U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
     +         + U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
      END IF
      IF( m .EQ. 9) THEN
         DivK=  -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
     +          +U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
      END IF
      IF( m .EQ.10) THEN
         DivK=  -U(4,K)-U(5,K)-U(2,K) +S2*(V(4,K)+V(1,K)+V(2,K))*AsR
      END IF
      IF( m .EQ.11) THEN
         DivK=   U(4,K)+U(3,K)+U(1,K) +S3*(V(4,K)+V(3,K)+V(1,K))*AsR
      END IF
      IF( m .EQ.12) THEN
         DivK=   U(4,K)+U(3,K)+U(6,K) -S5*(V(4,K)+V(6,K)+V(7,K))*AsR
      END IF
      IF( m .EQ.13) THEN
         DivK=  -U(4,K)-U(5,K)-U(7,K) -S6*(V(4,K)+V(5,K)+V(7,K))*AsR
      END IF

*     Open boundary radiation conditin by Flather,1976.

      IF( nb .EQ. -1) THEN

***      add= -dz(i,j)* DBLE(Dz_Trelax)
      add= -coef*dz(i,j)* SQRT(g/z(km21(i,j)))

      IF(m.GE.10)
     *           DivK=DivK +3.d0*(1.d0+S0*AsR)*add
      IF(m.EQ.6) DivK=DivK +6.d0*add
      IF(m.EQ.7) DivK=DivK +S0*6.d0*AsR*add
      IF(m.EQ.8) DivK=DivK +6.d0*add
      IF(m.EQ.9) DivK=DivK +S0*6.d0*AsR*add
      
	END IF


      return
      end

      subroutine ktform(KT)
      real KT(6,13)
      do ntr=1,6
      do n  =1,13
      KT(ntr,n)=0.
      end do
      end do

      do ntr=1,6
      KT(ntr,1)=1.
      end do

      KT(1,2)=1.
      KT(2,2)=1.
      KT(3,2)=1.
      KT(6,2)=1.

      KT(1,3)=1.
      KT(2,3)=1.
      KT(4,3)=1.
      KT(5,3)=1.
      KT(6,3)=1.

      KT(3,4)=1.
      KT(4,4)=1.
      KT(5,4)=1.
      KT(6,4)=1.

      KT(1,5)=1.
      KT(2,5)=1.
      KT(3,5)=1.
      KT(4,5)=1.
      KT(5,5)=1.

      KT(1,6)=1.
      KT(2,6)=1.
      KT(6,6)=1.

      KT(4,7)=1.
      KT(5,7)=1.
      KT(6,7)=1.

      KT(3,8)=1.
      KT(4,8)=1.
      KT(5,8)=1.

      KT(1,9)=1.
      KT(2,9)=1.
      KT(3,9)=1.

      KT(1,10)=1.
      KT(2,10)=1.

      KT(3,11)=1.

      KT(4,12)=1.
      KT(5,12)=1.

      KT(6,13)=1.
      return
      end




