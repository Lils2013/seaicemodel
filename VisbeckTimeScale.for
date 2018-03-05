      Program VisbeckTimeScale

      Parameter (kbeg=1, kfin=16)
      Parameter (Smax  = 5.e-3)   ! Maximum slope

	Parameter (il=35, jl=49, kl=16,klp=kl+1)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (ilp= il+2, jlp= jl+2)
      Parameter (ilris=51,jlris=51)
      Parameter (FiN=25., aLW= -25., hxgr= 1., hygr= 1.)
      Parameter (pi=3.14159265, overfl= 1.701410e+38)

      INTEGER*4 lenr
      external lenr
      Integer*2 wr_grd
      external wr_grd

      Integer*2 dummy
	
	real Agm(il,jl)   ! Variable coefficient at the triangle

      Dimension km2(0:il1,0:jl1), serv(ilris,jlris), Ro(il,jl,kl),
     *          T(0:il1,0:jl1,kl),S(0:il1,0:jl1,kl),serv2(il,jl,kl),
     *          z(klp),hz(kl),si(0:jl1)
      character*72 file

      hx=hxgr*pi/180.
      hy=hygr*pi/180.
      r=.637e9
      g=980.


*     Vertical grid
        z(1)= 0.
        z(2)= 10.
        z(3)= 25.
        z(4)= 50.
        z(5)= 100.
        z(6)= 150.
        z(7)= 200.
        z(8)= 250.
        z(9)= 300.
        z(10)=400.
        z(11)=500.
        z(12)=750.
        z(13)=1000.
        z(14)=2000.
        z(15)=3000.
        z(16)=4000.
        z(17)= 9999999.
*     -----------------
      do 33 k=1,klp
33    z(k)=100.*z(k)
      do 34 k=1,kl
34    hz(k)=(z(k+1)-z(k))

      do 52 j=0,jl1
      teta=pi*.5 - finord*pi/180.+(j-1.)*hy
52    si(j)=sin(teta)

      FiS= FiN - hygr*(jlris-1.)

	Rhx=R*hx
	Rhy=R*Hy

      open (unit=17,file='km2.dat',status='old',access='direct',
     *      form='unformatted',recl=(il+2)*(jl+2))
      read (17,rec=1) km2


	DO iMonth=1,12

	open (51,file='TPHC.dat',access='direct',form='unformatted',
     *      recl=il*jl*kl)
	read(51,rec=imonth) Serv2
	close(51)

      do i=1,il
      do j=1,jl
      do k=1,kl
      t(i,j,k)= serv2(i,jl-j+1,k)
      end do
      end do
      end do
*     GIN Sea
	do k=1,kl
      t(10,jl,k)= t(11,jl,k)
      end do

	open (51,file='SPHC.dat',access='direct',form='unformatted',
     *      recl=il*jl*kl)
	read(51,rec=imonth) Serv2
	close(51)

      do i=1,il
      do j=1,jl
      do k=1,kl
      s(i,j,k)= serv2(i,jl-j+1,k)
      end do
      end do
      end do
*     GIN Sea
	do k=1,kl
      s(10,jl,k)= s(11,jl,k)
      end do

	do i=1,il
	do j=1,jl
	if(km2(i,j).GT.0) then
	do k=1,km2(i,j)
	Ro(i,j,k) = sigma_t(t(i,j,k),s(i,j,k),0.) ! Potential density
	end do
	end if
	end do
	end do


*     Eddy diffusion coefficient calculation

      do j=1,jl-1

      S0=SI(J)
      SM=SI(J-1)
      SP=SI(J+1)
      S1=c3*(2.*S0+SM)
      S3=c3*(2.*SP+S0)

      do i=1,il-1

	Agm(i,j)= 0.
	Depth     = 0.

	K4= min(km2(i,j),km2(i+1,j),km2(i,j+1),km2(i+1,j+1))

      if(K4 .GT. kbeg+1) then

	do k= kbeg, min(kfin,K4-1)

	kp= k+1

c     Triangle # 1.

      Ro6K =Ro(i,j,k) +Ro(i+1,j,k) +Ro(i,j+1,k)
      Ro6Kp=Ro(i,j,kp)+Ro(i+1,j,kp)+Ro(i,j+1,kp)

      agm13= -1.5*hz(k)*(Ro(i+1,j,k )-Ro(i,j,k )+
     &                 Ro(i+1,j,kp)-Ro(i,j,kp))/
     &       (max(1.e-19,(Ro6kp-Ro6k))*Rhx*S1)
	agm23= -1.5*hz(k)*(Ro(i,j+1,k )-Ro(i,j,k )+
     &                 Ro(i,j+1,kp)-Ro(i,j,kp))/
     &       (max(1.e-19,(Ro6kp-Ro6k))*Rhy)
	Smodul1= MIN(Smax,sqrt(agm13**2+agm23**2))


	VB1= SQRT(max(0., g*(Ro6KP-Ro6K)/3.))

      Agm(i,j)= Agm(i,j) + Smodul1*VB1*SQRT(hz(k))
	Depth     = Depth +hz(k)
	end do  !  k
	else
	Agm(i,j)= 0.
	Depth= 1.
	end if

	Agm(i,j)= 1.e6*Agm(i,j)/Depth

	end do
	end do


*     Plot the Time scale
      ymin= FiS
      ymax= ymin +hygr*(jlris-1.)
      xmin= aLW
      xmax= xmin +hxgr*(ilris-1.)

      zmin= overfl
      zmax=-overfl

      file= 'TVisbeck00.grd'
      Write(file(9:10),'(i2.2)') iMonth

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

      Do i=10,44
      Do j=1,jl
      km=km2(i-9,jl-j+1)
      IF( km .GT. kbeg .and. Agm(i-9,jl-j+1).GT.0.) THEN
      Q= Agm(i-9,jl-j+1)
      serv(i,j)= Q
      zmin= min (zmin,Q)
      zmax= max (zmax,Q)
      ELSE
      serv (i,j)= overfl
      END IF
      end do
      end do
	nx=ilris
	ny=jlris
      dummy=
     *wr_grd(file,nx,ny,xmin,xmax,ymin,ymax,zmin,zmax,serv,nx)

      END DO  ! Month

      end


      function sigma_t(t,s,p)

C     SEA WATER DENSITY DEVIATION FROM 1.02 [GR/CM**3]
C     AS FUNCTION OF T[°C](potential),S[PPT],P[MPa]
C     By Bryden et al.:"A new approximation of the equation of state for
C                   seawater, suitable for numerical ocean models."
C     In: J.Geophys.Res.,v.104,No.C1, p.1537-1540, 1999.

C     VARIANT 2: -2<T<40;0<S<42;0<P<100.
      sigma_t = !!!-2.0092E-02 ! +(5.07043E-04*P-5.43283E-07*P*P)
     # + ( 5.10768E-05-3.69119E-06*P+6.54837E-09*P*P)*T
     # + ( 8.05999E-04-9.34012E-07*P+1.38777E-09*P*P)*S
     # + (-7.40849E-06+5.33243E-08*P-1.01563E-10*P*P)*T*T
     # + (-3.01036E-06+1.75145E-08*P-2.34892E-11*P*P)*T*S
     # + ( 3.32267E-08-3.25887E-10*P+4.98612E-13*P*P)*T*T*T
     # + ( 3.21931E-08-1.65849E-10*P+2.17612E-13*P*P)*T*T*S

      return
      end


      Function wr_grd(file,Nx,Ny,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,
     *                DATA,IDATA)
      INTEGER*2 wr_grd
      CHARACTER*(*) file
      real DATA (IDATA,Ny)

      if (Nx .gt. IDATA) then
         print *,'WR_GRD: Nx > IDATA. Nx=',Nx,', IDATA=',IDATA
         WR_GRD = -2
         goto 9999
      endif

         OPEN(9,FILE=file,FORM='FORMATTED',err=888)
         WRITE(9,'(a4)',err=888) 'DSAA'
         WRITE(9,*,err=888) Nx,Ny
         WRITE(9,*,err=888) Xmin, Xmax
         WRITE(9,*,err=888) Ymin, Ymax
         WRITE(9,*,err=888) Zmin, Zmax
         DO I=1,Ny
            WRITE(9,*,err=888) (data(j,I),j=1,Nx)
         ENDDO

      WR_GRD = 0
      goto 9999

  888 continue
      print *,'WR_GRD: Ошибка пpи записи во внешний файл ',
     +      file(1:lenr(file)), ' !'
      WR_GRD = -2

 9999 close (9)
      return
      end

      integer*4 function lenr(str)
      character*(*) str

      integer*2 i

      lenr = 0
      do 2 i=len(str),1,-1
         if (str(i:i) .ne. ' ')  then
            lenr = i
            goto 3
         end if
    2 continue

    3 continue
      return
      end

