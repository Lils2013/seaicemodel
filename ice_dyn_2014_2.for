      subroutine ice_dyn_evp(nstep)

*     Version 03.10.2015.

*     Gravity wave Drag for ice-ocean interaction
	include 'Slo2.fi'
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
	dimension serv(0:il1,0:jl1)
	integer abs_n
	include 'Tparm.fi'

c     -----Parameters for ice-ocean drag  ---------
	d_ice= 100.e2   !! mean icefloe diameter, 100-300m
      C100= row/d_ice

*     Geometry and cut-offs
      Ax= 1.0/Hx/R
      Ay= 1.0/Hy/R
      asr=hx/hy
      asr2=asr**2
      di= 10.0   !! Minimal dynamical ice mass(thickness in cm).

*     Time steps according to CICE recommendations 1:40:120.(2h)
*     Parameters for EVP integration
      Tdamp=dt/3.
	dte  =dt/REAL(nstep)

	serv = 0.

c     ----- Surface pressure and initial data for drift -----

	do j=1,jl
	do i=1,il

      kb= km2(i,j)

	if (kb.GT.0) then
	
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
      
      ro1=ropot(i,j,kmix)   
      ro2=ropot(i,j,kmix+1)
      hzk= hz(kmix)   

      VB= g*(ro2-ro1)/hzk/row
      VB= SQRT(MAX(0.,VB))      
      
c      if(i.eq.17.and.j.eq.24) then
c	write(*,*) kmix, vb
c	end if	
	
c     Sea ice mean draft (!), not thickness
	himean= 0.0
	do mg=1,mgrad
	himean=himean+roi*hice (mg,i,j)+rosdry*hsnow(mg,i,j)
	enddo
	himean = himean/row

c     Sea Ice draft variance
      hi_variance = 0.0
	do m=1,mgrad
	IF(aice(m,i,j) .GT. aimin) then
	total_mass= roi*hice(m,i,j)+rosdry*hsnow(m,i,j)
	hi= total_mass/MAX(aimin,aice(m,i,j))/row


*      For ridges parameterization!!!! BE CAREFUL!
cccc       if(hi.ge.200.) hi= 3.*hi
*      -------------------------------------------


	hi_variance= hi_variance +  Aice(m,i,j)* (hi-himean)**2
	end if   ! aice(m) > 0
	end do   ! mgrad


c     ---------- Gravity wave Drag formal coefficient --------------
      
      IF( kmix .GT. 7) then ! Thick upper mixed layer - no wave drag
      
      tauw =0.
      
      else
      
      param_nyu= 3.14*delta_u(i,j)/(max(1.e-9,VB) *d_ice)
      
      if(param_nyu .LE. 0.2)then  ! K0/Kc   
      tauw= 0.375* 3.14*C100 * hi_variance * VB
      else
      tauw= 0.075* row *hi_variance * VB**2 /max(1.e-9, delta_u(i,j))
      end if
      
      end if   ! VB >0
      

c     ---------------------------------------------------------

	CDgwd(i,j) = 0. !!!tauw

c     ---------------------------------------------------------


c      if(VB30 .GE. 1.e-9) then
c      if(i.eq.17.and.j.eq.24) then
c	write(*,*) i,j,Himean,hi_variance,param_nyu,
c     *VB, delta_u, tauw, row*Cdw*delta_u(i,j)
c	end if

      end if ! kb > 0
	end do
	end do


*====================================================================

* -------------------------------------------------------------------

*     Simplified ice stregth

c	do j=1,jl
c	do i=1,il
c      kb= km2(i,j)

c	if (kb.GT.0) then

*     Sea ice mean thickness
c	himean= 0.0
c	do mg=1,mgrad
c	himean=himean+hice(mg,i,j)
c	enddo
c	himean = himean/max(1.e-12,(1.-Aice(0,i,j)))
	
c	Pice(i,j)= 2.75e5*himean*EXP(-20.*aice(0,i,j)) ! Hibler 1979
c	Pice(i,j)= 1.00e5*himean*EXP(-20.*aice(0,i,j)) ! LIM3,dyn/cm2
c	Pice(i,j)= 2.00e5*himean*EXP(-20.*aice(0,i,j)) ! HadCEM,dyn/cm2
      
c      end if ! kb > 0
c	end do
c	end do

* -------------------------------------------------------------------

*====================================================================

*============     Total surface pressure on ice     =================
	do j=1,jl
	do i=1,il
      kb= km2(i,j)
	if (kb.GT.0) then

c     Sea ice-snow mean draft
	himean= 0.0
	do mg=1,mgrad
	himean=himean+roi*hice (mg,i,j)+rosdry*hsnow(mg,i,j)
**	himean=himean + rosdry*hsnow(mg,i,j) ! only snow
	enddo
	
      serv(i,j)=10.*PA(i,j) + g*(row*REAL(dz(i,j)) +himean)
cc      serv(i,j)= g*(row*REAL(dz(i,j)) +himean)

      end if ! kb > 0
	end do
	end do
      
*====================================================================      

*     Time integration with small time steps.

      DO miter=1,nstep

      DO 41 j=1,jL
      S0=Si(j)
      SM=Si(j-1)
      SP=Si(j+1)
      S1=(2.*S0+SM)/3.
      S2=(2.*SM+S0)/3.
      S3=S1
      S4=(2.*S0+SP)/3.
      S5=(2.*SP+S0)/3.
      S6=S4
	  	  
	DO 41 i=1,iL

      n=nt3(i,j,1)
	abs_n=ABS(n)
*     Rheology is only in the interior of the ocean. 
*     Some auxillary functions are everywhere.
      if(abs_n.GT.0) call EVP(i,j,miter,nstep,fu_evp,fv_evp,dte,Tdamp)


ccc      if(ice_mask(i,j).EQ.1) then 
ccc      if(ice_mask(i,j).EQ.1 .and. n.LE.1 .and. n.NE.0) then 
      if(ice_mask(i,j).EQ.1 .and. n.EQ.1) then
*     Points of Ice Motion. When no ice - formally effective thickness di.

      cor= -2.*om*co(i)*S0
      C_geometry=6.0/(s0*cg(abs_n))

      total_mass=0.0
      hi= 0.0
	do m=1,mgrad
	hi= hi +hice(m,i,j)
	total_mass= total_mass +roi*hice(m,i,j)+rosdry*hsnow(m,i,j)
	end do

*     Effective Thickness for Air and Ocean Drags.
	eff_thick=total_mass/MAX(aimin,1.0-aice(0,i,j))

*     Mean ice thickness
	hi= hi/MAX(aimin,1.0-aice(0,i,j))

*     Wind Stress Components Over Ice - Note Cdn second parameter
c      wmod=SQRT((wx(i,j)-0.01*uice2(i,j))**2+
c     +          (wy(i,j)-0.01*vice2(i,j))**2)
      wmod=SQRT(wx(i,j)**2+wy(i,j)**2)

*     Air drag

*     AOMIP - ice=water
c	Cdn_ice=Cdn(wmod,1)

*     MYCOM-ROMS
c	Cdn_ice=2.2e-3*0.5*(1.-cos(2*pi*min( 0.01*hi+.1, 0.5)))

*     Variable drag - dependence on ice thickness and wind 
	IF(hi.GE.300.)then 
	Cdn_ice=Cdn(wmod,2)
      else
	Cdn_ice=(hi*Cdn(wmod,2) +(300.-hi)*Cdn(wmod,1))/300.
	endif

	cdrag=roa*1.e4*Cdn_ice*wmod
      windxi= cdrag*wx(i,j)
      windyi= cdrag*wy(i,j)
c      windxi= cdrag*(wx(i,j)-0.01*uice2(i,j))
c      windyi= cdrag*(wy(i,j)-0.01*vice2(i,j))

c     ----------------Wave radiation pressure -------------------------
      Q=1.+1.6e-4*(row/roa)*SQRT(Aice(0,i,j)*(1.-Aice(0,i,j)))
      windxi= Q*windxi
      windyi= Q*windyi
c     -----------------------------------------------------------------

      ro0=ropot(i,j,1) 

      kb= km2(i,j)
      kmix= kb
      
      do k=1,kb-1
      ro2=ropot(i,j,k)
       
      if( abs(ro2-ro0) .GT. 1.25e-4) then ! Levitus, 1982
      kmix= k
      exit
      end if

      end do	! k
      
	CC=(row*Cdw*delta_u(i,j)+CDgwd(i,j))*dte/MAX(di,eff_thick)
      a11= 1.0+CC
      a22= a11
	a12= -0.5*cor*dte 
	a21= -a12
	det= 1.0/(a11**2+a12**2)

*     Water Friction.
      Fu= uice2(i,j)-a12*vice2(i,j)+
     + (row*Cdw*delta_u(i,j)*u(i,j,1)+CDgwd(i,j)*u(i,j,kmix))
     * *dte/MAX(di,eff_thick)
	Fv= vice2(i,j)+a12*uice2(i,j)+
     + (row*Cdw*delta_u(i,j)*v(i,j,1)+CDgwd(i,j)*v(i,j,kmix))
     * *dte/MAX(di,eff_thick)
     

*     Pressure at Sea Surface
      call gradPA(i,j,abs_n,derx,dery,serv,il1,jl1,KT)
      Fu= Fu+dte*Ax*C_geometry*derx/MAX(di,eff_thick)
      Fv= Fv+dte*Ay*C_geometry*dery/MAX(di,eff_thick)

*     Wind
      Fu= Fu +dte*windxi/MAX(di,eff_thick)
      Fv= Fv +dte*windyi/MAX(di,eff_thick)

*     Rheology
         Fu= Fu +dte*fu_evp/MAX(di,total_mass)/R/S0
         Fv= Fv +dte*fv_evp/MAX(di,total_mass)/R/S0

      uice(i,j)= det*(a22*Fu-a12*Fv)
      vice(i,j)= det*(a11*Fv-a21*Fu)

      else

      uice(i,j)=0.
      vice(i,j)=0.

      end if

      
cccc      write(*,*) i,j,uice(i,j),vice(i,j)

41    CONTINUE

c      DDD=0.
c      do j=1,jl
c	do i=1,il
c	delta_u= SQRT((uice2(i,j)-u(i,j,1))**2+(vice2(i,j)-v(i,j,1))**2)
c      delta_u= delta_u +1.e-3
c	if(CDgwd(i,j)/delta_u .GT. DDD) then
c	DDD=CDgwd(i,j)/delta_u
c      ro1=ro(i,j,1)
c      ro2=ro(i,j,2)
c      VB= 2.*g*(ro2-ro1)/hz(1)/(2.+ro1+ro2)
c      VB= SQRT(MAX(0.,VB))
c	i0=i
c	j0=j
c	end if
c	end do
c	end do

ccc      write(*,*) DDD,i0,j0,S(i0,j0,1),S(i0,j0,2),T(i0,j0,1),T(i0,j0,2)

*     Open boundary - Radiation
      call uicebc(uice,uice2,nt3,ice_mask,Si,il1,jl1,kl,hx,hy)
      call uicebc(vice,vice2,nt3,ice_mask,Si,il1,jl1,kl,hx,hy)

c     Reordering the time steps and regularization
      do i=1,il
      do j=1,jl
      uice2(i,j)= uice(i,j)
      vice2(i,j)= vice(i,j)
      end do
      end do

      END DO             ! miter, small time steps

      RETURN
      END

	subroutine uicebc(u,u2,nt3,ice_mask,Si,il1,jl1,kl,hx,hy)
*     Radiation boundary condition.
*     Version 14.07.06.

	dimension u(0:il1,0:jl1),u2(0:il1,0:jl1),
     &          nt3(0:il1,0:jl1,kl), Si(0:jl1),
     &          ice_mask(0:il1,0:jl1)

	do j=1,jl1-1
	S0=Si(j)
      do i=1,il1-1

	if(ice_mask(i,j) .EQ. 1) then

	ni=nt3(i,j,1)

	if(ni.EQ.-6) then

      delta_t=(u(i+1,j)-u2(i+1,j))
	delta_x=(u(i+1,j)-u(i+2,j)) /(hx*S0)
	sss= delta_t*(u2(i+1,j+1)-u2(i+1,j-1))

	if(sss.GT.0.)then
	delta_y=(u2(i+1,j)-u2(i+1,j-1)) /hy
	else
	delta_y=(u2(i+1,j+1)-u2(i+1,j)) /hy
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /hy
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	u(i,j)=D*(u2(i,j)+Rx*u(i+1,j)+Ry*(u2(i,j-1)-u2(i,j)))
	else
	u(i,j)=D*(u2(i,j)+Rx*u(i+1,j)-Ry*(u2(i,j+1)-u2(i,j)))
	end if

	else
      u(i,j)=   u(i+1,j)
	end if
	end if  ! n=-6

	if(ni.EQ.-7) then

      delta_t=(u(i,j+1)-u2(i,j+1))
	delta_x=(u(i,j+1)-u(i,j+2)) /hy
	sss= delta_t*(u2(i+1,j+1)-u2(i-1,j+1))

	if(sss.GT.0.)then
	delta_y=(u2(i,j+1)-u2(i-1,j+1)) /(hx*S0)
	else
	delta_y=(u2(i+1,j+1)-u2(i,j+1)) /(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /hy
	Ry=-delta_t*delta_y*D /(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	u(i,j)=D*(u2(i,j)+Rx*u(i,j+1)+Ry*(u2(i-1,j)-u2(i,j)))
	else
	u(i,j)=D*(u2(i,j)+Rx*u(i,j+1)-Ry*(u2(i+1,j)-u2(i,j)))
	end if

	else
      u(i,j)=  u(i,j+1)
	end if

	end if  ! n=-7

	if(ni.EQ.-8) then

      delta_t=(u(i-1,j)-u2(i-1,j))
	delta_x=(u(i-1,j)-u(i-2,j)) /(hx*S0)
	sss= delta_t*(u2(i-1,j+1)-u2(i-1,j-1))

	if(sss.GT.0.)then
	delta_y=(u2(i-1,j)-u2(i-1,j-1)) /hy
	else
	delta_y=(u2(i-1,j+1)-u2(i-1,j)) /hy
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /(hx*S0)
	Ry=-delta_t*delta_y*D /hy
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	u(i,j)=D*(u2(i,j)+Rx*u(i-1,j)+Ry*(u2(i,j-1)-u2(i,j)))
	else
	u(i,j)=D*(u2(i,j)+Rx*u(i-1,j)-Ry*(u2(i,j+1)-u2(i,j)))
	end if

	else
      u(i,j)= u(i-1,j)
	end if
	end if  ! n=-8

	if(ni.EQ.-9) then
      delta_t=(u(i,j-1)-u2(i,j-1))
	delta_x=(u(i,j-1)-u(i,j-2)) /hy
	sss= delta_t*(u2(i+1,j-1)-u2(i-1,j-1))

	if(sss.GT.0.)then
	delta_y=(u2(i,j-1)-u2(i-1,j-1)) /(hx*S0)
	else
	delta_y=(u2(i+1,j-1)-u2(i,j-1)) /(hx*S0)
	end if

      D= 1./MAX(1.e-12,delta_x**2 + delta_y**2)
	Rx=-delta_t*delta_x*D /hy
	Ry=-delta_t*delta_y*D /(hx*S0)
	D=1./(1.+Rx)

      if(Rx .GE. 0.)then

	if(Ry.GE.0.)then
	u(i,j)=D*(u2(i,j)+Rx*u(i,j-1)+Ry*(u2(i-1,j)-u2(i,j)))
	else
	u(i,j)=D*(u2(i,j)+Rx*u(i,j-1)-Ry*(u2(i+1,j)-u2(i,j)))
	end if

	else
      u(i,j)= u(i,j-1)
	end if
	end if  ! n=-9

      end if  ! mask=1

	end do  ! j
	end do  ! i

	return
	end


      subroutine gradPA(i,j,N,derx,dery,serv,il1,jl1,KT)
      dimension KT(6,13), serv(0:il1,0:jl1)
      real KT
      common /si/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2

      derx=-(   KT(1,n)*(serv(i+1,j  )-serv(i  ,j  ))+
     +          KT(2,n)*(serv(i+1,j-1)-serv(i  ,j-1))+
     +          KT(3,n)*(serv(i  ,j  )-serv(i-1,j  ))+
     +          KT(4,n)*(serv(i  ,j  )-serv(i-1,j  ))+
     +          KT(6,n)*(serv(i+1,j  )-serv(i  ,j  ))+
     +          KT(5,n)*(serv(i  ,j+1)-serv(i-1,j+1)))/6.

      dery=-(s1*KT(1,n)*(serv(i+1,j  )-serv(i+1,j-1))+
     +       s2*KT(2,n)*(serv(i  ,j  )-serv(i  ,j-1))+
     +       s3*KT(3,n)*(serv(i  ,j  )-serv(i  ,j-1))+
     +       s4*KT(4,n)*(serv(i-1,j+1)-serv(i-1,j  ))+
     +       s5*KT(5,n)*(serv(i  ,j+1)-serv(i  ,j  ))+
     +       s6*KT(6,n)*(serv(i  ,j+1)-serv(i  ,j  )))/6.
      return
      end


      subroutine EVP(i,j,itime,nstep,fu,fv,dte,Tdamp)
*     Force caused by rheology stresses.
*     Version 4.4 23.01.2013.
*     EVP sea ice rheology based on CICE v.3.0 time scheme.
*     Integration over ALL AREA.
*     New pressure is the constant piecewise function as well.
*     NOTE the special approximation of the COS function:
*     cos = (sin(j+1)-sin(j))/hy
*     All stresses and coefficients are constants at triagles.
*     Regularization either by Harder or by Hunke.
*     nreg={1,2} - {Hunke, Harder}
	include 'Slo2.fi'
	dimension e11(6),e12(6),e22(6),s11(6),s12(6),s22(6)
	include 'Tparm.fi'

	n=abs(nt3(i,j,1))

*     Method of regularization
      nreg= 1

      S0=Si(j)
      SM= Si(j-1)
      SP= Si(j+1)
      S1=(2.*S0+SM)/3.
      S2=(2.*SM+S0)/3.
      S3=S1
      S4=(2.*S0+SP)/3.
      S5=(2.*SP+S0)/3.
      S6=S4

	det1=1.0+0.5*dte/Tdamp
	det2=1.0+0.5*dte*extr2/Tdamp

	Fu=0.
	Fv=0.
*     Damping constant 615 kg/m2 ~ 61 cm.
**	C_Damp= 61.5*Tdamp*R*R*hx*hy*S0/dte/dte 
	C_Damp= 2800.*Tdamp*R*R*hx*hy*S0/dte/dte ! 6m, extr=2
**	C_Damp= 5600.*Tdamp*R*R*hx*hy*S0/dte/dte ! 12m, extr=2


c     All integrals are over the area, covered by ice.
c     Triangle assumed be covered by ice, if all the
c     corners are covered by sufficiantly compacted ice.

c 1.

      if(KT(1,n) .EQ. 1) then

      e11(1)=( (uice2(i+1,j  )-uice2(i  ,j  ))/hx
     #+(vice2(i,j)+vice2(i+1,j)+vice2(i+1,j-1))*(S0-SM)/(3.*hy)
     #         )/(R*S1)
	e22(1)=(vice2(i+1,j)-vice2(i+1,j-1))/(R*hy)
	e12(1)=0.5*(
     #      (uice2(i+1,j)-uice2(i+1,j-1))/(R*hy) +
     #      ( (vice2(i+1,j  )-vice2(i  ,j  ))/hx 
     #-(uice2(i,j)+uice2(i+1,j)+uice2(i+1,j-1))*(S0-SM)/(3.*hy)
     #            )/(R*S1))

      delta=(e11(1)+e22(1))**2+((e11(1)-e22(1))**2+4.*e12(1)**2)/extr2	
      delta=SQRT(delta)
	Pmean= (Pice(i,j)+Pice(i+1,j)+Pice(i+1,j-1))/3.

*     Regularization
      if(nreg .EQ. 1) then
      delta=MAX(dmin,delta)                  ! E. Hunke
      Pmean = MIN(Pmean,C_Damp*delta)        ! E. Hunke
	else
      Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
      delta=delta+dmin                       ! M. Harder.
      end if
	d1=delta

c     Stress tensor components.

      sigma1(1,i,j)=(sigma1(1,i,j)+
     +0.5*dte*Pmean*((e11(1)+e22(1))/delta-1.0)/Tdamp)/det1

      sigma2(1,i,j)=(sigma2(1,i,j)+
     +0.5*dte*Pmean*(e11(1)-e22(1))/delta/Tdamp )/det2

	sigma12(1,i,j)=(sigma12(1,i,j)+
     +0.5*dte*Pmean*e12(1)/delta/Tdamp )/det2  

c      sigma2(1,i,j)=(sigma2(1,i,j)+
c     +0.5*dte*Pmean*(e11(1)-e22(1))/delta/Tdamp/extr2 )/det1

c	sigma12(1,i,j)=(sigma12(1,i,j)+
c     +0.5*dte*Pmean*e12(1)/delta/Tdamp/extr2 )/det1  


	s11(1)= 0.5*(sigma1(1,i,j)+sigma2(1,i,j))
      s12(1)= sigma12(1,i,j)
      s22(1)= 0.5*(sigma1(1,i,j)-sigma2(1,i,j))

	else

	d1=0.

	s11(1)= 0.
      s12(1)= 0.
      s22(1)= 0.

	end if


c 2.

      if(KT(2,n) .EQ. 1) then

      e11(2)=( (uice2(i+1,j-1)-uice2(i  ,j-1))/hx
     #+(vice2(i,j)+vice2(i+1,j-1)+vice2(i,j-1))*(S0-SM)/(3.*hy)
     #        )/(R*S2)

	e22(2)=(vice2(i,j)-vice2(i,j-1))/(R*hy)
	e12(2)=0.5*(
     #      (uice2(i,j)-uice2(i,j-1))/(R*hy) +
     #      ( (vice2(i+1,j-1)-vice2(i,j-1))/hx
     #-(uice2(i,j)+uice2(i+1,j-1)+uice2(i,j-1))*(S0-SM)/(3.*hy)
     #            )/(R*S2))

      delta=(e11(2)+e22(2))**2+((e11(2)-e22(2))**2+4.*e12(2)**2)/extr2	
      delta=SQRT(delta)
	Pmean= (Pice(i,j)+Pice(i,j-1)+Pice(i+1,j-1))/3.

*     Regularization
      if(nreg .EQ. 1) then
      delta=MAX(dmin,delta)                  ! E. Hunke
      Pmean = MIN(Pmean,C_Damp*delta)        ! E. Hunke
	else
      Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
      delta=delta+dmin                       ! M. Harder.
      end if
	d2=delta

c     Stress tensor components.

      sigma1(2,i,j)=(sigma1(2,i,j)+
     +0.5*dte*Pmean*((e11(2)+e22(2))/delta-1.0)/Tdamp)/det1

      sigma2(2,i,j)=(sigma2(2,i,j)+
     +0.5*dte*Pmean*(e11(2)-e22(2))/delta/Tdamp )/det2

	sigma12(2,i,j)=(sigma12(2,i,j)+
     +0.5*dte*Pmean*e12(2)/delta/Tdamp )/det2  

c      sigma2(2,i,j)=(sigma2(2,i,j)+
c     +0.5*dte*Pmean*(e11(2)-e22(2))/delta/Tdamp/extr2 )/det1

c	sigma12(2,i,j)=(sigma12(2,i,j)+
c     +0.5*dte*Pmean*e12(2)/delta/Tdamp/extr2 )/det1  


	s11(2)= 0.5*(sigma1(2,i,j)+sigma2(2,i,j))
      s12(2)= sigma12(2,i,j)
      s22(2)= 0.5*(sigma1(2,i,j)-sigma2(2,i,j))

	else

	d2=0.

	s11(2)= 0.
      s12(2)= 0.
      s22(2)= 0.

	end if


c 3.

      if(KT(3,n) .EQ. 1) then

      e11(3)=( (uice2(i  ,j  )-uice2(i-1,j  ))/hx 
     #+(vice2(i,j)+vice2(i-1,j)+vice2(i,j-1))*(S0-SM)/(3.*hy)
     #        )/(R*S3)
	e22(3)=(vice2(i,j)-vice2(i,j-1))/(R*hy)
	e12(3)=0.5*(
     #      (uice2(i,j)-uice2(i,j-1))/(R*hy) +
     #      ( (vice2(i,j  )-vice2(i-1,j  ))/hx
     #-(uice2(i,j)+uice2(i-1,j)+uice2(i,j-1))*(S0-SM)/(3.*hy)
     #            )/(R*S3))

      delta=(e11(3)+e22(3))**2 +
     *     ((e11(3)-e22(3))**2+4.*e12(3)**2)/extr2	
      delta=SQRT(delta)
	Pmean= (Pice(i,j)+Pice(i-1,j)+Pice(i,j-1))/3.

*     Regularization
      if(nreg .EQ. 1) then
      delta=MAX(dmin,delta)                  ! E. Hunke
      Pmean = MIN(Pmean,C_Damp*delta)        ! E. Hunke
	else
      Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
      delta=delta+dmin                       ! M. Harder.
      end if
	d3=delta

c     Stress tensor components.

      sigma1(3,i,j)=(sigma1(3,i,j)+
     +0.5*dte*Pmean*((e11(3)+e22(3))/delta-1.0)/Tdamp)/det1

      sigma2(3,i,j)=(sigma2(3,i,j)+
     +0.5*dte*Pmean*(e11(3)-e22(3))/delta/Tdamp )/det2

	sigma12(3,i,j)=(sigma12(3,i,j)+
     +0.5*dte*Pmean*e12(3)/delta/Tdamp )/det2  

c      sigma2(3,i,j)=(sigma2(3,i,j)+
c     +0.5*dte*Pmean*(e11(3)-e22(3))/delta/Tdamp/extr2 )/det1

c	sigma12(3,i,j)=(sigma12(3,i,j)+
c     +0.5*dte*Pmean*e12(3)/delta/Tdamp/extr2 )/det1  

	s11(3)= 0.5*(sigma1(3,i,j)+sigma2(3,i,j))
      s12(3)= sigma12(3,i,j)
      s22(3)= 0.5*(sigma1(3,i,j)-sigma2(3,i,j))

	else

	d3=0.

	s11(3)= 0.
      s12(3)= 0.
      s22(3)= 0.

	end if


c 4.

      if(KT(4,n) .EQ. 1) then

      e11(4)=( (uice2(i  ,j  )-uice2(i-1,j  ))/hx
     #+(vice2(i,j)+vice2(i-1,j)+vice2(i-1,j+1))*(SP-S0)/(3.*hy)
     #        )/(R*S4)
	e22(4)=(vice2(i-1,j+1)-vice2(i-1,j))/(R*hy)
	e12(4)=0.5*(
     #      (uice2(i-1,j+1)-uice2(i-1,j))/(R*hy) +
     #      ( (vice2(i,j  )-vice2(i-1,j  ))/hx
     #-(uice2(i,j)+uice2(i-1,j)+uice2(i-1,j+1))*(SP-S0)/(3.*hy)
     #            )/(R*S4))

      delta=(e11(4)+e22(4))**2 +
     *     ((e11(4)-e22(4))**2+4.*e12(4)**2)/extr2	
      delta=SQRT(delta)
	Pmean= (Pice(i,j)+Pice(i-1,j)+Pice(i-1,j+1))/3.

*     Regularization
      if(nreg .EQ. 1) then
      delta=MAX(dmin,delta)                  ! E. Hunke
      Pmean = MIN(Pmean,C_Damp*delta)        ! E. Hunke
	else
      Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
      delta=delta+dmin                       ! M. Harder.
      end if
	d4=delta

c     Stress tensor components.

      sigma1(4,i,j)=(sigma1(4,i,j)+
     +0.5*dte*Pmean*((e11(4)+e22(4))/delta-1.0)/Tdamp)/det1

      sigma2(4,i,j)=(sigma2(4,i,j)+
     +0.5*dte*Pmean*(e11(4)-e22(4))/delta/Tdamp )/det2

	sigma12(4,i,j)=(sigma12(4,i,j)+
     +0.5*dte*Pmean*e12(4)/delta/Tdamp )/det2  

c      sigma2(4,i,j)=(sigma2(4,i,j)+
c     +0.5*dte*Pmean*(e11(4)-e22(4))/delta/Tdamp/extr2 )/det1

c	sigma12(4,i,j)=(sigma12(4,i,j)+
c     +0.5*dte*Pmean*e12(4)/delta/Tdamp/extr2 )/det1  

	s11(4)= 0.5*(sigma1(4,i,j)+sigma2(4,i,j))
      s12(4)= sigma12(4,i,j)
      s22(4)= 0.5*(sigma1(4,i,j)-sigma2(4,i,j))

	else

	d4=0.

	s11(4)= 0.
      s12(4)= 0.
      s22(4)= 0.

	end if


c 5.

      if(KT(5,n) .EQ. 1) then

      e11(5)=( (uice2(i,j+1)-uice2(i-1,j+1))/hx
     #+(vice2(i,j)+vice2(i-1,j+1)+vice2(i,j+1))*(SP-S0)/(3.*hy)
     #        )/(R*S5)
	e22(5)=(vice2(i,j+1)-vice2(i,j))/(R*hy)
	e12(5)=0.5*(
     #      (uice2(i,j+1)-uice2(i,j))/(R*hy) +
     #      ( (vice2(i,j+1)-vice2(i-1,j+1))/hx
     #-(uice2(i,j)+uice2(i,j+1)+uice2(i-1,j+1))*(SP-S0)/(3.*hy)
     #            )/(R*S5))

      delta=(e11(5)+e22(5))**2 +
     *     ((e11(5)-e22(5))**2+4.*e12(5)**2)/extr2	
      delta=SQRT(delta)
	Pmean= (Pice(i,j)+Pice(i,j+1)+Pice(i-1,j+1))/3.

*     Regularization
      if(nreg .EQ. 1) then
      delta=MAX(dmin,delta)                  ! E. Hunke
      Pmean = MIN(Pmean,C_Damp*delta)        ! E. Hunke
	else
      Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
      delta=delta+dmin                       ! M. Harder.
      end if
	d5=delta

c     Stress tensor components.

      sigma1(5,i,j)=(sigma1(5,i,j)+
     +0.5*dte*Pmean*((e11(5)+e22(5))/delta-1.0)/Tdamp)/det1

      sigma2(5,i,j)=(sigma2(5,i,j)+
     +0.5*dte*Pmean*(e11(5)-e22(5))/delta/Tdamp )/det2

	sigma12(5,i,j)=(sigma12(5,i,j)+
     +0.5*dte*Pmean*e12(5)/delta/Tdamp )/det2  

c      sigma2(5,i,j)=(sigma2(5,i,j)+
c     +0.5*dte*Pmean*(e11(5)-e22(5))/delta/Tdamp/extr2 )/det1

c	sigma12(5,i,j)=(sigma12(5,i,j)+
c     +0.5*dte*Pmean*e12(5)/delta/Tdamp/extr2 )/det1  

	s11(5)= 0.5*(sigma1(5,i,j)+sigma2(5,i,j))
      s12(5)= sigma12(5,i,j)
      s22(5)= 0.5*(sigma1(5,i,j)-sigma2(5,i,j))

	else

	d5=0.

	s11(5)= 0.
      s12(5)= 0.
      s22(5)= 0.

	end if


c 6.  

      if(KT(6,n) .EQ. 1) then
      e11(6)=( (uice2(i+1,j  )-uice2(i  ,j  ))/hx
     #+(vice2(i,j)+vice2(i+1,j)+vice2(i,j+1))*(SP-S0)/(3.*hy)
     #        )/(R*S6)
	e22(6)=(vice2(i,j+1)-vice2(i,j))/(R*hy)
	e12(6)=0.5*(
     #      (uice2(i,j+1)-uice2(i,j))/(R*hy) +
     #      ( (vice2(i+1,j  )-vice2(i  ,j  ))/hx
     #-(uice2(i,j)+uice2(i+1,j)+uice2(i,j+1))*(SP-S0)/(3.*hy)
     #            )/(R*S6))

      delta=(e11(6)+e22(6))**2 +
     *     ((e11(6)-e22(6))**2+4.*e12(6)**2)/extr2	
      delta=SQRT(delta)
	Pmean= (Pice(i,j)+Pice(i+1,j)+Pice(i,j+1))/3.

*     Regularization
      if(nreg .EQ. 1) then
      delta=MAX(dmin,delta)                  ! E. Hunke
      Pmean = MIN(Pmean,C_Damp*delta)        ! E. Hunke
	else
      Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
      delta=delta+dmin                       ! M. Harder.
      end if
	d6=delta

c     Stress tensor components.

      sigma1(6,i,j)=(sigma1(6,i,j)+
     +0.5*dte*Pmean*((e11(6)+e22(6))/delta-1.0)/Tdamp)/det1

      sigma2(6,i,j)=(sigma2(6,i,j)+
     +0.5*dte*Pmean*(e11(6)-e22(6))/delta/Tdamp )/det2

	sigma12(6,i,j)=(sigma12(6,i,j)+
     +0.5*dte*Pmean*e12(6)/delta/Tdamp )/det2  

c      sigma2(6,i,j)=(sigma2(6,i,j)+
c     +0.5*dte*Pmean*(e11(6)-e22(6))/delta/Tdamp/extr2 )/det1

c	sigma12(6,i,j)=(sigma12(6,i,j)+
c     +0.5*dte*Pmean*e12(6)/delta/Tdamp/extr2 )/det1  

	s11(6)= 0.5*(sigma1(6,i,j)+sigma2(6,i,j))
      s12(6)= sigma12(6,i,j)
      s22(6)= 0.5*(sigma1(6,i,j)-sigma2(6,i,j))

	else

	d6=0.

	s11(6)= 0.
      s12(6)= 0.
      s22(6)= 0.

	end if

*-----------------------------------------------------------------

      if(n.eq.1) then

      Fu=  0.5*((s11(1)+s11(6))-(s11(3)+s11(4)))/hx
     #   + 0.5*((s12(5)*S5+s12(6)*S6)-(s12(2)*S2+s12(3)*S3))/hy
     #   + ( (s12(1)+s12(2)+s12(3))*(S0-SM)+
     #       (s12(4)+s12(5)+s12(6))*(SP-S0) )/(6.*hy)

      Fv=  0.5*((s12(1)+s12(6))-(s12(3)+s12(4)))/hx
     #   + 0.5*((s22(5)*S5+s22(6)*S6)-(s22(2)*S2+s22(3)*S3))/hy
     #   - ( (s11(1)+s11(2)+s11(3))*(S0-SM)+
     #       (s11(4)+s11(5)+s11(6))*(SP-S0) )/(6.*hy)

	else

	Fu=0.
	Fv=0.

	end if

*------Delta and divergence in nodes for ice ridging calculations.
      IF(itime.EQ.nstep) THEN
	delta_ice(i,j)= (d1+d2+d3+d4+d5+d6) /cg(n)
	div_ice(i,j)=(e11(1)+e11(2)+e11(3)+e11(4)+e11(5)+e11(6)+
     +              e22(1)+e22(2)+e22(3)+e22(4)+e22(5)+e22(6)) /cg(n)

cc	if(i.eq.17.and.j.eq.24)write(*,*)'F=',Fu,Fv
	END IF

      return
      end
