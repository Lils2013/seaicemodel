      Program TOPO_M

*     version 15.11.2013

*     Программа для выдачи результатов работы coupled
*     модели в формате, пригодном для SURFER-a.
*     Программа работает с выдачей подпрограммы FDisc.

*     Веpсия 2.0 для Северного Ледовитого Океана.
*     Вертикальная скорость рисуется в центре отрезка по вертикали.
*     Скорость осредняется по каждому треугольнику.
*     При рисовании уровня есть две возможности осреднения:
*     в центр квадрата и фильтром по шаблону.


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
      Character*72 udat,i2dat,tdat,zdat,file,mdat,ndat,idat,hdat,
     *             CH4dat,FCH4dat, O2dat

      Dimension u(0:il1,0:jl1,kl),v(0:il1,0:jl1,kl),dz(-1:ilp,-1:jlp),
     *          km2(0:il1,0:jl1), coefvd(il,jl,kl),
     *          serv(ilris,jlris), Fch4(0:il1,0:jl1)

*     Ice/Snow parameters
      dimension uice(0:il1,0:jl1), vice(0:il1,0:jl1),
     *          HIce (0:mgrad,0:il1,0:jl1),Aice(0:mgrad,0:il1,0:jl1),
     *          HSnow(0:mgrad,0:il1,0:jl1)

      Dimension servm(ilrism,jlrism)
      Double precision dz

      FiS= FiN - hygr*(jlris-1.)

      write(*,*) ' Enter Year and Month ->'
      read(*,*) iYear, iMonth

*     Формирование имени файла

      udat= 'u000000'
      tdat= 't000000'
      zdat= 'z000000'
      mdat= 'M000000'
      ndat= 'N000000'
      idat= 'I000000'
      hdat= 'H000000'
      i2dat='km2.dat'
      CH4dat='CH4000000'
      O2dat='O2000000'
      FCH4dat='FluxCH4000000'

      Write(udat(2:5),'(i4.4)') iYear
      Write(tdat(2:5),'(i4.4)') iYear
      Write(zdat(2:5),'(i4.4)') iYear
      Write(mdat(2:5),'(i4.4)') iYear
      Write(ndat(2:5),'(i4.4)') iYear
	Write(idat(2:5),'(i4.4)') iYear
      Write(hdat(2:5),'(i4.4)') iYear
      Write(CH4dat(4:7),'(i4.4)') iYear
      Write(O2dat(3:6),'(i4.4)') iYear

      Write(udat(6:7),'(i2.2)') iMonth
      Write(tdat(6:7),'(i2.2)') iMonth
      Write(zdat(6:7),'(i2.2)') iMonth
      Write(mdat(6:7),'(i2.2)') iMonth
      Write(ndat(6:7),'(i2.2)') iMonth
	Write(idat(6:7),'(i2.2)') iMonth
      Write(hdat(6:7),'(i2.2)') iMonth
      Write(CH4dat(8:9),'(i2.2)') iMonth
      Write(O2dat(7:8),'(i2.2)') iMonth

      write(FCH4dat(8:11),'(i4.4)') iYear
      write(FCH4dat(12:13),'(i2.2)') iMonth


      ymin= FiS
      ymax= ymin +hygr*(jlris-1.)
      xmin= aLW
      xmax= xmin +hxgr*(ilris-1.)

      n2 = (il+2)*(jl+2)*Lfort
      n22= (il+4)*(jl+4)*2*Lfort
      n3 = n2*kl
      nm = il*jl*kl*Lfort
      nuice= Lfort*(il1+1)*(jl1+1)
      nhice= nuice*(mgrad+1)

      Open ( Unit=49,File=i2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)

      Read(49,Rec=1) km2

      write(*,*) ' Number of SURFER version: DOS/WINDOWS [1/2] -->'
      read(*,*) nver

      write(*,*) ' Write Velocity ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN

      Open ( Unit=31,File=udat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(31,Rec=1) u
      Read(31,Rec=2) v
      Close(31)

      file='u00000000.dat'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth
      umin= overfl
      umax=-overfl
      umean= 0.
      m= 0
      iris= il
      jris= jl

      Do 1 k=1,kl
      Write(file(8:9),'(i2.2)') k
      Open(12,FILE=file)

      Do 2 j=1,jL
       y1= FiS +hygr*(j-0.33333)
       y2= FiS +hygr*(j-0.66666)
       ymean= FiS+hygr*(j-0.5)
         Do 2 i=1,iL
         kmean= MIN(km2(i,jl-j+1),km2(i+1,jl-j+1),km2(i,jl-j),
     *              km2(i+1,jl-j))

         IF( kmean .GE. k) THEN

         um1= (u(i,jl-j  ,k)+u(i+1,jl-j,k)+u(i  ,jl-j+1,k))/3.
         vm1= (v(i,jl-j  ,k)+v(i+1,jl-j,k)+v(i  ,jl-j+1,k))/3.
         um2= (u(i,jl-j+1,k)+u(i+1,jl-j,k)+u(i+1,jl-j+1,k))/3.
         vm2= (v(i,jl-j+1,k)+v(i+1,jl-j,k)+v(i+1,jl-j+1,k))/3.

         umod1= SQRT(um1**2 +vm1**2)
         umod2= SQRT(um2**2 +vm2**2)

         umean= umean + 0.5*(umod1+umod2)
         m= m+1

         if(umod1.GT.umax) then
         umax= umod1
         i0=i
         j0=j
         k0=k
         end if

         if(umod2.GT.umax) then
         umax= umod2
         i0=i
         j0=j
         k0=k
         end if

*************************************************************
*        Ограничение размера стрелочки при рисовании.
               IF( umod .GT. 50.) umod= 50.
*************************************************************

               IF( umod1 .GT. 1.e-4) THEN
               x1= aLW + hxgr*(i-0.66666)   +hxgr*9.
         uvalue1= 180./pi*ATAN2(-vm1, um1)

         Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x1,y1,umod1,uvalue1
               END IF

               IF( umod2 .GT. 1.e-4) THEN
               x2= aLW + hxgr*(i-0.33333)   +hxgr*9.
         uvalue2= 180./pi*ATAN2(-vm2, um2)
      if(nver.EQ.2) then
         Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x2,y2,umod2,uvalue2
      end if
               END IF

         END IF
2     Continue

      umean= umean/(m)

      write(*,*)' umean=', umean

      x= 20.
      y= 22.
      umod= 6.
      uvalue= 0.
      Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x,y,umod,uvalue

      x= 20.
      y= 19.
      umod= 2.
      uvalue= 0.
      Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x,y,umod,uvalue
      if(nver.EQ.2) then
      x= 20.
      y= 19.
      umod= 2.
      uvalue= 0.
      Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x,y,umod,uvalue
      end if
      Close(12)

1     Continue

      write(*,*) 'umax= ', umax, i0,j0,k0

      END IF


      write(*,*) 'Write Mean Ice/Snow Thickness & Compactness? -->'
      read(*,*) ansis

      IF( ansis .EQ. 1) THEN
      Open ( Unit=31,File=hdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nhice)
      read(31,rec=1) Hice
      read(31,rec=2) Aice
      read(31,rec=3) Hsnow
      close(31)

      do i=1,il1-1
      do j=1,jl1-1
      Hi=0.
      Hs=0.
      Ai=0.
      if( km2(i,j) .GT. 0) then

      do m=1,mgrad
      Hi= Hi +Hice(m,i,j)
      Hs= Hs +Hsnow(m,i,j)
      Ai= Ai+Aice(m,i,j)
      end do

c     Mean Thickness in the point.
      if (Ai .GT. 1.e-1) then
      Hice(0,i,j) =Hi/Ai
      Hsnow(0,i,j)=Hs/Ai
      Aice(0,i,j)= Ai
      else
      Hice(0,i,j)=0.
      Hsnow(0,i,j)=0.
      Aice(0,i,j)=0.
      end if

      else
      Hice(0,i,j)=0.
      Hsnow(0,i,j)=0.
      Aice(0,i,j)=0.
      end if

      end do
      end do

      zmin= overfl
      zmax=-overfl

      file= 'Hi000000.grd'
      Write(file(3:6),'(i4.4)') iYear
      Write(file(7:8),'(i2.2)') iMonth

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

      Do i=10,44
      Do j=1,jl
      km=km2(i-9,jl-j+1)
      IF( km .GT. 0) THEN
      Q= Hice(0,i-9,jl-j+1)
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


      file= 'Hs000000.grd'
      Write(file(3:6),'(i4.4)') iYear
      Write(file(7:8),'(i2.2)') iMonth

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

      Do i=10,44
      Do j=1,jl
      km=km2(i-9,jl-j+1)
      IF( km .GT. 0) THEN
      Q= Hsnow(0,i-9,jl-j+1)
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

      zmin= overfl
      zmax=-overfl

      file= 'Ai000000.grd'
      Write(file(3:6),'(i4.4)') iYear
      Write(file(7:8),'(i2.2)') iMonth

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

      write(*,*) ' Write Ice Velocity ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN

*     Assumed that in Ai(0,i,j) is the array to cut off open water!!!

      open (unit=31,file=idat,status='old',access='direct',
     *      form='unformatted',recl=nuice)
      read(31,rec=1) uice
      read(31,rec=2) vice
      close(31)

      file='ui000000.dat'
      Write(file(3:6),'(i4.4)') iYear
      Write(file(7:8),'(i2.2)') iMonth
      umin= overfl
      umax=-overfl
      umean= 0.
      m= 0
      iris= il
      jris= jl

      Open(12,FILE=file)

      Do 22 j=1,jL
       y1= FiS +hygr*(j-1.)
         Do 22 i=1,iL
         kmean= km2(i,jl-j+1)

         IF( kmean.GT.0 .and. Aice(0,i,jL-j+1).GT. 0.1) THEN

         um1= uice(i,jl-j+1)
         vm1= vice(i,jl-j+1)

         umod= SQRT(um1**2 +vm1**2)

         umean= umean +umod
         m= m+1

         if(umod.GT.umax) then
         umax= umod
         i0=i
         j0=j
         end if

*************************************************************
*        Ограничение размера стрелочки при рисовании.
               IF( umod .GT. 50.) umod= 50.
*************************************************************

               IF( umod .GT. 1.e-2) THEN
               x1= aLW + hxgr*(i-1.)   +hxgr*9.
         uvalue= 180./pi*ATAN2(-vm1, um1)

         Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x1,y1,umod,uvalue
               END IF


         END IF
22    continue
      umean= umean/(m)

      write(*,*)' umean=', umean

      x= 20.
      y= 22.
      umod= 20.
      uvalue= 0.
      Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x,y,umod,uvalue

      x= 20.
      y= 19.
      umod= 10.
      uvalue= 0.
      Write(12,'(f8.4,1x,f8.4,e15.7,e17.7)') x,y,umod,uvalue
      Close(12)

      write(*,*) 'uicemax= ', umax, i0,j0

      END IF

      write(*,*) 'Write Temperature T ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN
      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=1) u
      Close(32)

      file= 't00000000.grd'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth

      Do 3 k=1,kl

      write(file(8:9),'(i2.2)') k

      zmin= overfl
      zmax=-overfl
      tmean= 0.

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do 4 i=10,44
           Do 4 j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= u(i-9,jl-j+1,k) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
4          CONTINUE

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

3     CONTINUE
      END IF

      write(*,*) 'Write Salinity S ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN
      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=2) u
      Close(32)

      file= 's00000000.grd'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth

      Do 33 k=1,kl

      write(file(8:9),'(i2.2)') k
      zmin= overfl
      zmax=-overfl
      tmean= 0.
      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do 44 i=10,44
           Do 44 j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= u(i-9,jl-j+1,k) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
44         CONTINUE

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

33    CONTINUE
      END IF

      write(*,*) 'Write methane CH4 ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN
      Open ( Unit=32,File=CH4dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=1) u
      Close(32)

      file= 'ch400000000.grd'
      Write(file(4:7),'(i4.4)') iYear
      Write(file(8:9),'(i2.2)') iMonth

      Do k=1,kl

      write(file(10:11),'(i2.2)') k

      zmin= overfl
      zmax=-overfl
      tmean= 0.

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do i=10,44
           Do j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= u(i-9,jl-j+1,k) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
           end do
           end do

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

      end do
      END IF

      write(*,*) 'Write Oxygen ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN
      Open ( Unit=32,File=O2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=1) u
      Close(32)

      file= 'O200000000.grd'
      Write(file(3:6),'(i4.4)') iYear
      Write(file(7:8),'(i2.2)') iMonth

      Do k=1,kl

      write(file(9:10),'(i2.2)') k

      zmin= overfl
      zmax=-overfl
      tmean= 0.

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do i=10,44
           Do j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= u(i-9,jl-j+1,k) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
           end do
           end do

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

      end do
      END IF


      write(*,*) 'Write methane CH4 flux to air ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN
      Open ( Unit=32,File=FCH4dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)
      Read(32,Rec=1) Fch4
      Close(32)

      file= 'fch4000000.grd'
      Write(file(5:8),'(i4.4)') iYear
      Write(file(9:10),'(i2.2)') iMonth

      zmin= overfl
      zmax=-overfl
      tmean= 0.

      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do i=10,44
           Do j=1,jl
           IF( km2(i-9,jl-j+1) .GE. 1) THEN
           q= Fch4(i-9,jl-j+1) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
           end do
           end do

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

      END IF


      write(*,*) 'Write Vertical Coefficient M,N? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN


      Open ( Unit=32,File=mdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nm)
      Read(32,Rec=1) coefvd
      Close(32)

      file= 'm00000000.grd'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth

      Do k=1,kl

      write(file(8:9),'(i2.2)') k
      zmin= overfl
      zmax=-overfl
      tmean= 0.
      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do  i=10,44
           Do  j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= coefvd(i-9,jl-j+1,k) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
           end do
	     end do

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

      end do

      Open ( Unit=32,File=ndat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nm)
      Read(32,Rec=1) coefvd
      Close(32)

      file= 'n00000000.grd'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth

      Do k=1,kl

      write(file(8:9),'(i2.2)') k
      zmin= overfl
      zmax=-overfl
      tmean= 0.
      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do  i=10,44
           Do  j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= coefvd(i-9,jl-j+1,k) +tmean
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
           end do
	     end do

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

      end do

      END IF


      write(*,*) 'Write Vertical velocity W ? -->'
      read(*,*) answ

      IF( answ .EQ. 1) THEN
      Open ( Unit=32,File=udat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=3) u
      Close(32)

      file= 'w00000000.grd'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth

      Do 93 k=1,kl

      write(file(8:9),'(i2.2)') k

      zmin= overfl
      zmax=-overfl
      tmean= 0.
      do i=1,ilris
      do j=1,jlris
      serv(i,j)= overfl
      end do
      end do

           Do 94 i=10,44
           Do 94 j=1,jl
           IF( km2(i-9,jl-j+1) .GE. k) THEN
           q= 1000.*u(i-9,jl-j+1,k)
           serv(i,j)= q
           zmin= min (zmin,q)
           zmax= max (zmax,q)
           ELSE
           serv(i,j)= overfl
           END IF
94         CONTINUE

      dummy=
     *wr_grd(file,ilris,jlris,xmin,xmax,ymin,ymax,zmin,zmax,serv,ilris)

93    CONTINUE
      END IF

*     Это надо делать только для версии SURFER под DOS.
*     Marking of the Upwelling Zone

      if(nver.EQ.1) then

      file='w00000000.dat'
      Write(file(2:5),'(i4.4)') iYear
      Write(file(6:7),'(i2.2)') iMonth
      Do 81 k=1,kl
      Write(file(8:9),'(i2.2)') k
      Open(12,FILE=file)

      Do 82 j=1,jl
      y= FiS +hygr*(j-1.)
      Do 82 i=1,il

      q= 1000.*u(i,jl-j+1,k)

         IF( km2(i,jl-j+1) .GE. k .AND. q .LT. 0.0) THEN
         x= aLW + hxgr*(i-1.)  +hxgr*9.
         Write(12,'(f8.4,1x,f8.4,1x,a2)') x,y,'59'
         END IF

82    Continue
      Close(12)

81    Continue
      end if

      write(*,*) 'Write Sea Level Dz ? -->'
      read(*,*) ansdz

      IF( ansdz .EQ. 1) THEN
      Open ( Unit=32,File=zdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n22)
      Read(32,Rec=1) dz
      Close(32)

      zmin= overfl
      zmax=-overfl

      file= 'dz000000.grd'
      Write(file(3:6),'(i4.4)') iYear
      Write(file(7:8),'(i2.2)') iMonth

      do i=1,ilris
      do j=1,jlris
        serv(i,j)= overfl
      end do
      end do

      write(*,*)'You have two possible ways of interpolation:'
      write(*,*)'interpolate in the center of the square [1],'
      write(*,*)'or use more sophisticated methods [2].'
      write(*,*)'Enter No of filtration method --->'
      read(*,*) ans

             m=1

      IF( ans .EQ. 1) THEN
      do i=1,ilrism
      do j=1,jlrism
        servm(i,j)= overfl
      end do
      end do
             m=2
      ymin= ymin+.5*hygr
      ymax= ymax -.5*hygr
      xmin= xmin +.5*hxgr
      xmax= xmax -.5*hxgr
      END IF

      Do 14 i=10,44-m+1
      Do 14 j=1,jl-m+1

*     Определение целочисленной глубины
      IF(ans.EQ.1) THEN
*     При интерполяции в центр квадрата
      km=min(km2(i-9,jl-j+1),km2(i-8,jl-j+1),
     #       km2(i-9,jl-j  ),km2(i-8,jl-j  ))
      ELSE
*     При интерполяции в расчетные точки
      km=km2(i-9,jl-j+1)
      END IF

      IF( km .GT. 0) THEN
        IF(ans.EQ.1) THEN

*       Интерполяция по четырем точкам в центр квадрата

        Q= 0.25*(dz(i-8,jl-j+1)+ dz(i-9,jl-j  )+
     #           dz(i-8,jl-j  )+ dz(i-9,jl-j+1))
        servm(i,j)= Q
        ELSE

*       Интерполяция в расчетные точки "естественным" фильтром

        k1=min(km2(i-9 ,jl-j+1),km2(i-8,jl-j+1),
     #         km2(i-9 ,jl-j  ),km2(i-8,jl-j  ))
        k2=min(km2(i-9 ,jl-j+2),km2(i-8,jl-j+2),
     #         km2(i-9 ,jl-j+1),km2(i-8,jl-j+1))
        k3=min(km2(i-10,jl-j+2),km2(i-9,jl-j+2),
     #         km2(i-10,jl-j+1),km2(i-9,jl-j+1))
        k4=min(km2(i-10,jl-j+1),km2(i-9,jl-j+1),
     #         km2(i-10,jl-j  ),km2(i-9,jl-j  ))

*           Интерполяция в центры квадратов
                         dz1=0.
                         dz2=0.
                         dz3=0.
                         dz4=0.
                         a1=0.
                         a2=0.
                         a3=0.
                         a4=0.

        if( k1.GT.0) then
        dz1=0.25*(dz(i-8,jl-j+1)+ dz(i-9,jl-j  )+
     #            dz(i-8,jl-j  )+ dz(i-9,jl-j+1))
        a1=1.
        end if
        if( k2.GT.0) then
        dz2=0.25*(dz(i-8,jl-j+2)+ dz(i-9,jl-j+1)+
     #            dz(i-8,jl-j+1)+ dz(i-9,jl-j+2))
        a2=1.
        end if
        if( k3.GT.0) then
        dz3=0.25*(dz(i-10,jl-j+2)+ dz(i-9,jl-j+1)+
     #            dz(i-10,jl-j+1)+ dz(i-9,jl-j+2))
        a3=1.
        end if
        if( k4.GT.0) then
        dz4=0.25*(dz(i-10,jl-j+1)+ dz(i-9,jl-j  )+
     #            dz(i-10,jl-j  )+ dz(i-9,jl-j+1))
        a4=1.
        end if

*               Суммирование по квадратам

        Q= (dz1+dz2+dz3+dz4)/(a1+a2+a3+a4)
        serv(i,j)= Q

        END IF
                   zmin= min (zmin,Q)
                   zmax= max (zmax,Q)
      ELSE
         IF( ans .EQ. 1) THEN
         servm(i,j)= overfl
         ELSE
         serv (i,j)= overfl
         END IF
      END IF
14    CONTINUE

      IF(ans.EQ.1) THEN
      nx=ilrism
      ny=jlrism
      dummy=
     *wr_grd(file,nx,ny,xmin,xmax,ymin,ymax,zmin,zmax,servm,nx)

      ELSE

      nx=ilris
      ny=jlris
      dummy=
     *wr_grd(file,nx,ny,xmin,xmax,ymin,ymax,zmin,zmax,serv,nx)
      END IF

      END IF
      Close(49)
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
