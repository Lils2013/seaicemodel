      Program TOPO_Ridge

*     version 12.02.2014

*     Large arrays for cluster

*     Lfort=1 for Compaq and UNIX Fortrans, and
*      Lfort=4 for NDP, Watcom and MAC Fortrans.
      Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=40)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (ilp= il+2, jlp= jl+2)
      Parameter (ilris=51,jlris=51,ilrism=ilris-1,jlrism=jlris-1)
      Parameter (FiN=25., aLW= -25., hxgr= 1., hygr= 1.)
      Parameter (pi=3.14159265, overfl= 1.701410e+38)
      Parameter (mgrad= 14)

      INTEGER*4 lenr
      external lenr
      Integer*2 wr_grd
      external wr_grd

      Integer*2 dummy
      Character*72 file,name

      Dimension km2(0:il1,0:jl1), serv(ilris,jlris)

*     Ice/Snow parameters
      dimension HIce (0:mgrad,0:il1,0:jl1),Aice(0:mgrad,0:il1,0:jl1)

      FiS= FiN - hygr*(jlris-1.)

      write(*,*) ' Enter Year and Month ->'
      read(*,*) iYear, iMonth


      ymin= FiS
      ymax= ymin +hygr*(jlris-1.)
      xmin= aLW
      xmax= xmin +hxgr*(ilris-1.)

      nnn=(il1+1)*(jl1+1)*kl*Lfort
      nnf=(il1+1)*(jl1+1)*Lfort
      nn2= Lfort*(il1+3)*(jl1+3)*2
      nn = Lfort*il*jl*kl
      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)



      Open ( Unit=49,File='km2.dat',Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nnf)
      Read(49,Rec=1) km2
	close(49)


      write(*,*) 'Write Mean Ice/Snow Thickness & Compactness? -->'
      read(*,*) ansis

      IF( ansis .EQ. 1) THEN

	nrec= 12*(iYear-1948) + 3*(iMonth-1) +1

      name='himonth.dat'
      open (unit=31,file=name,status='old',access='direct',
     *      form='unformatted',recl=nhice)
      read(31,rec= nrec)   Hice
      read(31,rec= nrec+1) Aice
      close(31)

      do i=1,il1-1
      do j=1,jl1-1
      
      Himean=0.
      Ai=0.
      
      if( km2(i,j) .GT. 0) then

      do m=0,mgrad
      Himean= Himean +Hice(m,i,j)
ccc      Ai= Ai+Aice(m,i,j)
      end do

c     Mean Thickness in the point.
      if (Aice(0,i,j) .GT. 1.e-6) then
      Himean =Himean    !!!!/Ai
      else
      Himean=0.
      end if
      
      sigma2=0.
*     Ice thickness Variance 
      do m=0,mgrad
      hi= Hice(m,i,j)/max(1.e-6, Aice(m,i,j))
      sigma2= sigma2 +Aice(m,i,j)*(Hi-Himean)**2
      end do
      
      if(i.eq.17 .and. j.eq.24) write(*,*) himean, sqrt(sigma2)
  
  
      ARidge=0.
      Aice(0,i,j)=0.
      
      do m=1,mgrad
      
      hi= Hice(m,i,j)/max(1.e-6, Aice(m,i,j))

      if(hi .GT. himean + 30. .and. Himean .GT. 100.)then
cc      if(hi .GT. himean + 0.3*sqrt(sigma2) .and. Himean .GT. 100.)then
cc      if( hi .GT. 500. )then

cc      ARidge= ARidge +Aice(m,i,j)
      
      Aice(0,i,j)= Aice(0,i,j) +1.e6*Aice(m,i,j)/(hi**2)      
      end if
      
      end do
      
      Aice(0,i,j) =SQRT(Aice(0,i,j))
cc      Aice(0,i,j) = ARidge
    
c      if(Himean .GT. 1.e-6) 
cc     * Aice(0,i,j)= 1.e3*sqrt(ARidge)/300.
cc     * Aice(0,i,j)= 1.e3*sqrt(ARidge)/(himean + 2.*SQRT(sigma2))
  
      else
      Aice(0,i,j)=0.
      end if

      end do
      end do

      zmin= overfl
      zmax=-overfl

      file= 'AiR000000.grd'
      Write(file(4:7),'(i4.4)') iYear
      Write(file(8:9),'(i2.2)') iMonth

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

      Do i=10,44
      Do j=1,jl
      km=km2(i-9,jl-j+1)
      IF( km .GT. 0) THEN
      Q= Aice(0,i-9,jl-j+1)
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

      END IF

      Stop
      End

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
