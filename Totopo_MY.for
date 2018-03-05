      Program TOPO_MY

*     ┬√ўшёыхэшх ьэюуюыхЄэшї ёЁхфэхьхё ўэ√ї яюыхщ.
*     ┬хЁёш  схч ъю¤ЇшЎхэЄют тхЁЄшъры№эющ ЄєЁсєыхэЄэюёЄш 
*     ш схч тхЁЄшъры№эющ ёъюЁюёЄш
     
      Parameter (Lfort=1)

	Parameter (il=35, jl=49, kl=16)
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
      Character*72 udat,i2dat,tdat,zdat,file,mdat,idat,hdat

      Dimension u(0:il1,0:jl1,kl),v(0:il1,0:jl1,kl),dz(-1:ilp,-1:jlp),
     *          km2(0:il1,0:jl1), coefvd(il,jl,kl),
     *          serv(ilris,jlris),
     *          umy(12,0:il1,0:jl1,kl),vmy(12,0:il1,0:jl1,kl)

*     Ice/Snow parameters
      dimension uice(0:il1,0:jl1), vice(0:il1,0:jl1),
     *          HIce (0:mgrad,0:il1,0:jl1),Aice(0:mgrad,0:il1,0:jl1),
     *          HSnow(0:mgrad,0:il1,0:jl1)
      dimension uiceMY(12,0:il1,0:jl1), viceMY(12,0:il1,0:jl1),
     *          HIceMY (12,0:mgrad,0:il1,0:jl1),
     *          AiceMY(12,0:mgrad,0:il1,0:jl1),
     *          HSnowMY(12,0:mgrad,0:il1,0:jl1), dzMY(12,-1:ilp,-1:jlp)


      Dimension servm(ilrism,jlrism)
      Double precision dz

      FiS= FiN - hygr*(jlris-1.)

      udat= 'u000000'
      tdat= 't000000'
      zdat= 'z000000'
      mdat= 'M000000'
      idat= 'I000000'
      hdat= 'H000000'
      i2dat='km2.dat'

      Write(udat(2:5),'(i4.4)') iYear
      Write(tdat(2:5),'(i4.4)') iYear
      Write(zdat(2:5),'(i4.4)') iYear
      Write(mdat(2:5),'(i4.4)') iYear
	Write(idat(2:5),'(i4.4)') iYear
      Write(hdat(2:5),'(i4.4)') iYear

      Write(udat(6:7),'(i2.2)') iMonth
      Write(tdat(6:7),'(i2.2)') iMonth
      Write(zdat(6:7),'(i2.2)') iMonth
      Write(mdat(6:7),'(i2.2)') iMonth
	Write(idat(6:7),'(i2.2)') iMonth
      Write(hdat(6:7),'(i2.2)') iMonth

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

      umy(:,:,:,:) = 0.
      vmy(:,:,:,:) = 0.

	do iMonth=1,12
	do iYear= 1958, 2010

      Write(udat(2:5),'(i4.4)') iYear
      Write(udat(6:7),'(i2.2)') iMonth

      Open ( Unit=31,File=udat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(31,Rec=1) u
      Read(31,Rec=2) v
      Close(31)
      umy(iMonth,:,:,:) = umy(iMonth,:,:,:) +u(:,:,:)/REAL(2010-1958+1)
      vmy(iMonth,:,:,:) = vmy(iMonth,:,:,:) +v(:,:,:)/REAL(2010-1958+1)

	End do   ! iYear
	End do   ! iMonth

	Do iMonth=1,12

      file='umy0000.dat'
      Write(file(4:5),'(i2.2)') iMonth
      umin= overfl
      umax=-overfl
      umean= 0.
      m= 0
      iris= il
      jris= jl

      Do k=1,kl
      Write(file(6:7),'(i2.2)') k
      Open(12,FILE=file)

      Do j=1,jL
       y1= FiS +hygr*(j-0.33333)
       y2= FiS +hygr*(j-0.66666)
       ymean= FiS+hygr*(j-0.5)
         Do i=1,iL
         kmean= MIN(km2(i,jl-j+1),km2(i+1,jl-j+1),km2(i,jl-j),
     *              km2(i+1,jl-j))

         IF( kmean .GE. k) THEN

         um1= (umy(imonth,i,jl-j  ,k)+umy(imonth,i+1,jl-j,k)+
     *         umy(imonth,i  ,jl-j+1,k))/3.
         vm1= (vmy(imonth,i,jl-j  ,k)+vmy(imonth,i+1,jl-j,k)+
     *         vmy(imonth,i  ,jl-j+1,k))/3.
         um2= (umy(imonth,i,jl-j+1,k)+umy(imonth,i+1,jl-j,k)+
     *         umy(imonth,i+1,jl-j+1,k))/3.
         vm2= (vmy(imonth,i,jl-j+1,k)+vmy(imonth,i+1,jl-j,k)+
     *         vmy(imonth,i+1,jl-j+1,k))/3.

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
               IF( umod .GT. 100.) umod= 100.
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
      End do    ! i,j
      End do    ! Loop



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

	End do    ! K
	End do    ! iMonth



      write(*,*) 'umax= ', umax, i0,j0,k0

      END IF


      write(*,*)'Write Mean Ice/Snow Thickness, Compactness, Drift? -->'
      read(*,*) ansis

      IF( ansis .EQ. 1) THEN

	HiceMY(:,:,:,:)=  0.
	AiceMY(:,:,:,:)=  0.
	HsnowMY(:,:,:,:)= 0.
	UiceMY(:,:,:)=  0.
	ViceMY(:,:,:)=  0.

	do iMonth=1,12
	do iYear= 1958, 2010
      Write(hdat(2:5),'(i4.4)') iYear
      Write(hdat(6:7),'(i2.2)') iMonth
	Write(idat(2:5),'(i4.4)') iYear
	Write(idat(6:7),'(i2.2)') iMonth

      Open ( Unit=31,File=hdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=nhice)
      read(31,rec=1) Hice
      read(31,rec=2) Aice
      read(31,rec=3) Hsnow
      close(31)

      open (unit=31,file=idat,status='old',access='direct',
     *      form='unformatted',recl=nuice)
      read(31,rec=1) uice
      read(31,rec=2) vice
      close(31)

	HiceMY(imonth,:,:,:)=  HiceMY(imonth,:,:,:)  + Hice(:,:,:)/
     *REAL(2010-1958+1)
	AiceMY(imonth,:,:,:)=  AiceMY(imonth,:,:,:)  + Aice(:,:,:)/
     *REAL(2010-1958+1)
	HsnowMY(imonth,:,:,:)= HsnowMY(imonth,:,:,:) + Hsnow(:,:,:)/
     *REAL(2010-1958+1)

	UiceMY(imonth,:,:)=  UiceMY(imonth,:,:)  + Uice(:,:)/
     *REAL(2010-1958+1)
	ViceMY(imonth,:,:)=  ViceMY(imonth,:,:)  + Vice(:,:)/
     *REAL(2010-1958+1)

	end do
	end do

	Do iMonth=1,12

      do i=1,il1-1
      do j=1,jl1-1
      Hi=0.
      Hs=0.
      Ai=0.
      if( km2(i,j) .GT. 0) then

      do m=1,mgrad
      Hi= Hi +HiceMY(imonth,m,i,j)
      Hs= Hs +HsnowMY(imonth,m,i,j)
      Ai= Ai+AiceMY(imonth,m,i,j)
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

      file='Hicemy00.grd'
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


      file='Hsnwmy00.grd'
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

      file='Aicemy00.grd'
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

      file='Uicemy00.dat'
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

         um1= uiceMY(imonth,i,jl-j+1)
         vm1= viceMY(imonth,i,jl-j+1)

         umod= SQRT(um1**2 +vm1**2)

         umean= umean +umod
         m= m+1

         if(umod.GT.umax) then
         umax= umod
         i0=i
         j0=j
         end if

*************************************************************
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

	end do ! Month

      END IF

      write(*,*) 'Write Temperature T ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN

      umy(:,:,:,:) = 0.

	do iMonth=1,12
	do iYear= 1958, 2010

      Write(tdat(2:5),'(i4.4)') iYear
      Write(tdat(6:7),'(i2.2)') iMonth

      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=1) u
      Close(32)

      umy(iMonth,:,:,:) = umy(iMonth,:,:,:) +u(:,:,:)/REAL(2010-1958+1)

	end do
	end do

	Do iMonth=1,12

      file= 'tMY0000.grd'
      Write(file(4:5),'(i2.2)') iMonth

      Do k=1,kl

      write(file(6:7),'(i2.2)') k

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
           q= umy(imonth,i-9,jl-j+1,k) 
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
	end do
      END IF

      write(*,*) 'Write Salinity S ? -->'
      read(*,*) ans

      IF( ans .EQ. 1) THEN

      umy(:,:,:,:) = 0.

	do iMonth=1,12
	do iYear= 1958, 2010

      Write(tdat(2:5),'(i4.4)') iYear
      Write(tdat(6:7),'(i2.2)') iMonth

      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=2) u
      Close(32)

      umy(iMonth,:,:,:) = umy(iMonth,:,:,:) +u(:,:,:)/REAL(2010-1958+1)

	end do
	end do

	Do iMonth=1,12

      file= 'sMY0000.grd'
      Write(file(4:5),'(i2.2)') iMonth

      Do k=1,kl

      write(file(6:7),'(i2.2)') k

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
           q= umy(imonth,i-9,jl-j+1,k) 
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
	end do
      END IF



      write(*,*) 'Write Sea Level Dz ? -->'
      read(*,*) ansdz

      IF( ansdz .EQ. 1) THEN

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

      dzmy(:,:,:) = 0.

	do iMonth=1,12
	do iYear= 1958, 2010

      Write(zdat(2:5),'(i4.4)') iYear
      Write(zdat(6:7),'(i2.2)') iMonth


      Open ( Unit=32,File=zdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n22)
      Read(32,Rec=1) dz
      Close(32)

      dzmy(iMonth,:,:) = dzmy(iMonth,:,:) +REAL(dz(:,:))/
     *REAL(2010-1958+1)

	end do
	end do

	do iMonth=1,12

      do i=1,il
	do j=1,jl
	
	Hi=0.
	Hs=0.
	do mg=1,mgrad
cccc      Hi= Hi +HiceMY(imonth,mg,i,j)
      Hs= Hs +HsnowMY(imonth,mg,i,j)
      end do

	dzmy(imonth,i,j) = dzmy(imonth,i,j) + 0.9*Hi+0.3*Hs

	end do
	end do

*     Mean level

      dzmean= 0.
	npp= 0.
      do i=1,il
	do j=1,jl
	IF( km2(i,j) .GT. 0) then
	dzmean= dzmean +dzmy(imonth,i,j)
	npp= npp +1
	END IF
	end do
	end do
	dzmean= dzmean/REAL(npp)

      do i=1,il
	do j=1,jl
	IF( km2(i,j) .GT. 0) then
	dzmy(imonth,i,j)= dzmy(imonth,i,j) - dzmean
	END IF
	end do
	end do

	end do ! month

	do imonth=1,12

      zmin= overfl
      zmax=-overfl
      file= 'dzmy00.grd'
      Write(file(5:6),'(i2.2)') iMonth

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

        Q= 0.25*(dzmy(imonth,i-8,jl-j+1)+ dzmy(imonth,i-9,jl-j  )+
     #           dzmy(imonth,i-8,jl-j  )+ dzmy(imonth,i-9,jl-j+1))
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
        dz1=0.25*(dzmy(imonth,i-8,jl-j+1)+ dzmy(imonth,i-9,jl-j  )+
     #            dzmy(imonth,i-8,jl-j  )+ dzmy(imonth,i-9,jl-j+1))
        a1=1.
        end if
        if( k2.GT.0) then
        dz2=0.25*(dzmy(imonth,i-8,jl-j+2)+ dzmy(imonth,i-9,jl-j+1)+
     #            dzmy(imonth,i-8,jl-j+1)+ dzmy(imonth,i-9,jl-j+2))
        a2=1.
        end if
        if( k3.GT.0) then
        dz3=0.25*(dzmy(imonth,i-10,jl-j+2)+ dzmy(imonth,i-9,jl-j+1)+
     #            dzmy(imonth,i-10,jl-j+1)+ dzmy(imonth,i-9,jl-j+2))
        a3=1.
        end if
        if( k4.GT.0) then
        dz4=0.25*(dzmy(imonth,i-10,jl-j+1)+ dzmy(imonth,i-9,jl-j  )+
     #            dzmy(imonth,i-10,jl-j  )+ dzmy(imonth,i-9,jl-j+1))
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

      END DO
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
