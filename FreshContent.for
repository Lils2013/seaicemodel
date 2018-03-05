       Program Fresh_Content
*     Lfort=1 for Compaq and UNIX Fortrans, and
*     Lfort=4 for NDP, Watcom and MAC Fortrans.
      Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=40)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (ilp= il+2, jlp= jl+2)
      Parameter (ilris=51,jlris=51,ilrism=ilris-1,jlrism=jlris-1)
      Parameter (FiN=25., aLW= -25., hxgr= 1., hygr= 1., finord=23.)
      Parameter (pi=3.14159265, overfl= 1.701410e+38)

      INTEGER*4 lenr
      external lenr
      Integer*2 wr_grd
      external wr_grd

      Integer*2 dummy
      Character*72 udat,i2dat,nt3dat,tdat,zdat,file,mdat,idat,hdat

      Dimension u(0:il1,0:jl1,kl),v(0:il1,0:jl1,kl),dz(-1:ilp,-1:jlp),
     *          km2(0:il1,0:jl1), coefvd(il,jl,kl),nt3(0:il1,0:jl1,kl),
     *serv(ilris,jlris,kl),serv_i(ilris,jlris), si(jl1), hz(kl),z(kl),
     *cg(13)
      data cg/6.,4.,5.,4.,5.,4*3.,2.,1.,2.,1./
                        
	r=.637e9
      hy=hygr*pi/180.
	r2hy2=(r*hy)**2 /1.e15

*     Vertical grid
        z(1)= 0.
        z(2)= 2.5
        z(3)= 5.
        z(4)= 7.5
        z(5)= 10.
        z(6)= 15.
        z(7)= 20.
        z(8)= 27.5
        z(9)= 37.5
        z(10)=50.
        z(11)=65.
        z(12)=82.5
        z(13)=100.
        z(14)=125.
        z(15)=150.
        z(16)=175.
        z(17)=200.
        z(18)=225.
        z(19)=250.
        z(20)=275.
        z(21)=300.
        z(22)=325.
        z(23)=350.
        z(24)=375.
        z(25)=410.
        z(26)=450.
        z(27)=500.
        z(28)=600.
        z(29)=700.
        z(30)=850.
        z(31)=1050.
        z(32)=1250.
        z(33)=1500.
        z(34)=1750.
        z(35)=2100.
        z(36)=2500.
        z(37)=3000.
        z(38)=3500.
        z(39)=4000.
        z(40)=4500.
        z(41)= 9999999.

      do 4 k=1,kl
4     hz(k)=(z(k+1)-z(k))

      do 52 j=0,jl1
      teta=pi*.5 - finord*pi/180.+(j-1.)*hy
52    si(j)=sin(teta)


      FiS= FiN - hygr*(jlris-1.)

      write(*,*) ' Enter Year and Month ->'
      read(*,*) iYear, iMonth

      tdat= 't000000'
      i2dat='km2.dat'
      nt3dat='nt.dat'

      Write(tdat(2:5),'(i4.4)') iYear
      Write(tdat(6:7),'(i2.2)') iMonth

      ymin= FiS
      ymax= ymin +hygr*(jlris-1.)
      xmin= aLW
      xmax= xmin +hxgr*(ilris-1.)

      n2 = (il+2)*(jl+2)*Lfort
      n3 = n2*kl
      nm = il*jl*kl*Lfort

      Open ( Unit=49,File=i2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)
      Read(49,Rec=1) km2

      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=2) u
      Close(32)

	open (unit=11,file=nt3dat,status='old',access='direct',
     *      form='unformatted',recl=n3)
      read (11,rec=1) nt3
      close(11)

      file= 'fwc00000000.grd'
      Write(file(4:7),'(i4.4)') iYear
      Write(file(8:9),'(i2.2)') iMonth

      Do 33 k=1,kl-1

	depth=hz(k)

      write(file(10:11),'(i2.2)') k
      zmin= overfl
      zmax=-overfl
      tmean= 0.
      do i=1,ilris
      do j=1,jlris
      serv(i,j,k)=overfl
      end do
      end do

           Do 44 i=10,44
           Do 44 j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k+1) THEN
	     n= abs(nt3(i-9,jl-j+1,k+1))
           q=  u(i-9,jl-j+1,k)
           qp= u(i-9,jl-j+1,k+1)
	     s=(q+qp)/2.
           serv(i,j,k)= depth*MAX(0.,(34.8-s))/34.8
           zmin= min (zmin,serv(i,j,k))
           zmax= max (zmax,serv(i,j,k))
           ELSE
           serv(i,j,k)= overfl
           END IF
44         CONTINUE
      do i=1,ilris
      do j=1,jlris
      serv_i(i,j)= serv(i,j,k)
      end do
      end do

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,
     *serv_i,ilris)

33    CONTINUE

      do i=10,44
	do j=1,jl
      serv_i(i,j)=serv(i,j,1)
      do k=2,kl-1
	if(serv(i,j,k).LT.overfl) serv_i(i,j)=serv_i(i,j) + serv(i,j,k)
	end do
	end do
	end do
      write(file(10:11),'(i2.2)') kl
      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,(0.),(40.),
     *serv_i,ilris)

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
