      Program SWNoCloud
      parameter (il=35, jl=49,dthr=1.0)
      Parameter (ilris=51,jlris=51)
       dimension serv2(ilris,jlris)
	dimension SWmean(il,jl)

	character*80 file
      real*8 zmin, zmax


	SWmean(:,:)= 0.

	write(*,*) 'Enter day'
	read(*,*)   m1day

      do i=1,il
	do j=1,jl

	do mhour=1,24
	call shortwave(SW,m1day,mhour,dthr,i,j)
      SWmean(i,j)= SWmean(i,j) +SW
	end do

      SWmean(i,j) = SWmean(i,j)/24.
	end do
	end do

	serv2(:,:)= 1.701410e+38

      zmin=+99999999.
      zmax=-99999999.

	Do i=10,44
      Do j=1,jl
      serv2(i,j)= SWmean(i-9,j)
      zmin= min (zmin,serv2(i,j))
      zmax= max (zmax,serv2(i,j))
      end do
      end do

        file='SW000.grd'
        write(file(3:5),'(i3.3)') m1day
	open(1,file=file)
        write(1,'(a4)') 'DSAA'
        write(1,*) ilris, jlris
        write(1,*) -25.0, 25.0
        write(1,*) -25.0, 25.0
        write(1,*) zmin, zmax
        do j=1,jlris
	write(1,*) (serv2(i,j), i=1,ilris)
        end do
	close(1)
	end		


      subroutine shortwave(SW,m1day,mhour,dthr,i,j)

c     Simple version with no Clouds and Vapor.

c     Short wave radiation calculation. Cloud is the cloud cover.
c     Zillman, 1972, Parkinson & Washington, 1979.

	parameter (S=1353.*1.e3)  ! Solar constant 
	parameter (hx=1.,hy=1.)   ! Grid size in degrees.
	parameter(pi=3.14159265)
	parameter(Alpha=0., Beta=-90., Gamma=0.)
	dimension R_Forward(3,3)
	parameter(FiSmodel=-25., aLoWmodel=-16.)


c     Rotation matrix into geographical coordinates
      cosa= COS(Alpha*pi/180.)
	sina= SIN(Alpha*pi/180.)
      cosb= COS(Beta*pi/180.)
	sinb= SIN(Beta*pi/180.)
      cosg= COS(Gamma*pi/180.)
	sing= SIN(Gamma*pi/180.)
 
      R_Forward(1,1)=  cosg*cosb*cosa - sing*sina
      R_Forward(1,2)= -sing*cosb*cosa - cosg*sina
      R_Forward(1,3)=  sinb*cosa
      R_Forward(2,1)=  cosg*cosb*sina + sing*cosa
      R_Forward(2,2)= -sing*cosb*sina + cosg*cosa
      R_Forward(2,3)=  sinb*sina
      R_Forward(3,1)= -cosg*sinb
      R_Forward(3,2)=  sing*sinb
      R_Forward(3,3)=  cosb



	Xin= aLoWmodel+(i-1.)*hx
	Yin= FiSmodel +(j-1.)*hy
	call rotation(Xin,Yin,Xout,Fi,R_Forward)


      write(*,*) i,j,xout,fi

c     Vapor pressure [Pa]
      evapor= 0.

c     Declination angle
      delta= 23.44*cos((172.-REAL(m1day))*pi/180.)
	  sd=sind(delta)
	  cd=cosd(delta)

c     Solar time at the location
	stime=dthr*(REAL(mhour)-0.5) ! GMT
	stime=stime +Xout/15.

	CZ= MAX(0.,sind(fi)*sd+cosd(fi)*cd*cos((12.-stime)*pi/12.))
	Q0=S*(CZ**2)/((CZ+2.7)*evapor*1.e-5 +1.085*CZ+0.10)

      SW=Q0

	return
	end

	subroutine Rotation(Xin,Yin,Xout,Yout,R)
	dimension R(3,3)
	parameter (Pi=3.1415926)

c     A coordinate transformation from coordinates
c     with rotated grid, Euler angles (X_IN_ANGLE,Y_IN_ANGLE,0),
c     to geographical coordinates (long, teta).
c     R - Euler matrix of rotation.
c     (Xin, Yin) - initial model coordinates,
c     (Xout,Yout) - geographical coordinates.

c     avoid troble of an exactly zero angle by adding the offset.      
c     and converting to radians.
      eps= 1.e-5
	  Xin= Pi*(Xin +eps)/180.
	  Yin= Pi*(Yin +eps)/180.

c     Spherical coordinates to cartesian coordinates
      xx= COS(Yin)*COS(Xin)
	yy= COS(Yin)*SIN(Xin)
	zz= SIN(Yin)

c     New cartesian coordinates
      Xnew= R(1,1) * xx
     &     +R(1,2) * yy
     &     +R(1,3) * zz

      Ynew= R(2,1) * xx
     &     +R(2,2) * yy
     &     +R(2,3) * zz

      Znew= R(3,1) * xx
     &     +R(3,2) * yy
     &     +R(3,3) * zz

c     New angles teta_new, costn, phi_new
      teta_new = ASIN(Znew)
	costn= MAX(eps,SQRT(1.0 - Znew**2))

	IF      ((Xnew.GT.0.) .AND. (Ynew.GT.0.)) THEN

	     IF(Xnew.LT.Ynew) THEN
		 phi_new= ACOS(Xnew/costn)
		 ELSE
		 phi_new= ASIN(Ynew/costn)
		 ENDIF

      ELSEIF ((Xnew.LT.0.) .AND. (Ynew.GT.0.)) THEN
      
	      IF(ABS(Xnew).LT.Ynew) THEN
		  phi_new= pi-ACOS(ABS(Xnew)/costn)
		  ELSE
		  phi_new= pi-ASIN(    Ynew /costn)
		  ENDIF

      ELSEIF ((Xnew.LT.0.) .AND. (Ynew.LT.0.)) THEN

	      IF(ABS(Xnew).LT.ABS(Ynew)) THEN
		  phi_new= -pi+ACOS(ABS(Xnew)/costn)
		  ELSE
		  phi_new= -pi+ASIN(ABS(Ynew)/costn)
		  ENDIF

      ELSEIF ((Xnew.GT.0.) .AND. (Ynew.LT.0.)) THEN

	      IF( Xnew.LT.ABS(Ynew)) THEN
		  phi_new= -ACOS(ABS(Xnew)/costn)
		  ELSE
		  phi_new= -ASIN(ABS(Ynew)/costn)
		  ENDIF

	ENDIF

c     New spherical coordinates.

      Xout= Phi_new  *180./Pi
	Yout= teta_new *180./Pi
ccc        IF( Xout.LT. 0.) Xout= Xout +360.

c     avoid troble of an exactly zero angle by subtracting the offset.      
	Xout= Xout -eps
	Yout= Yout -eps

	END
