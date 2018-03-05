      Program mainFCT

*     Version 16.07.2012

	parameter (dthr=1.0)
      parameter (pi=3.1415926)
      parameter (hxgr=1., hygr=1., finord=15., alwest=-15.)

      INCLUDE 'Slo2.fi'
	dimension MEC(12), aiday(il,jl), serv(il,jl)
      common /sin/ s0,sm,sp,s1,s2,s3,s4,s5,s6,asr,asr2
      data cg/6.,4.,5.,4.,5.,4*3.,2.,1.,2.,1./
	data MEC/31,28,31,30,31,30,31,31,30,31,30,31/


      dt=3600.*dthr 
	hx=hxgr*pi/180.
      hy=hygr*pi/180.
      r=.637e9

	km2(:,:)= 0
      do j=1,jl
	do i=1,il
	km2(i,j)= 4
	end do
	end do


	aice =0.
c	aice1=0.
	aice2=0.

      do j=1,jl
	do i=1,il
	uice(i,j) =10. !!!*(31.-i)/REAL(il)
	vice(i,j) =0.
	end do
	end do

	aice(1,3:7,14:18)= 1.
	aice2(1,3:7,14:18)= 1.

	call ntform (km2, nt3, il1, jl1, kl)
	do j=1,jl
	nt3(il,j,1)= -nt3(il,j,1)
	end do

      call ktform(KT)

      do 52 j=0,jl1
      teta=pi*.5 - finord*pi/180.+(j-1.)*hy
52    si(j)=sin(teta)


c     Ice mask for fast ice
	do j=1,jl
      do i=1,il
	ice_mask(i,j)= 1
cc      if(nt3(i,j,1).EQ.1) ice_mask(i,j)=1
	end do
	end do

      iDayglobal=0

	DO imonth= 1,12
	DO iday=1,MEC(imonth)

	idayglobal= idayglobal +1
c     Day means.
	aiday(:,:)= 0.

*     Steps in the Day
	DO L=1,INT(24./dthr)

      Aice2 = Aice
      call iceadvect(Aice,Aice1,Aice2)

      do j=1,jl
	do i=1,il
	aiday(i,j)= aiday(i,j) +aice(1,i,j)
	end do
	end do

	END DO ! day


c     Day means.
	aiday(:,:)= aiday(:,:)/(24./dthr)

      if( MOD(iDayGlobal,10) .EQ. 0) then
	write(*,*)' Day=', iDayGlobal
      call Fdisc(aiday,serv,iDayGlobal,il,jl)
	end if

	END DO ! month
	END DO ! year

	stop 
	end

      subroutine ntform(km2,nt3,il1,jl1,kl)
      dimension km2(0:il1,0:jl1), nt3(0:il1,0:jl1,kl)
*     Point type array.
*     Version 26.05.97 with wide arrays.
      do i=0,il1
      do j=0,jl1
      do k=1,kl
      nt3(i,j,k)=0
      end do
      end do
      end do

      do k=1,kl
      DO 11 i=1,il1-1
      DO 11 j=1,jl1-1
      IF(km2(i,j) .GE. k) THEN
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=1
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).LT.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=2
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).LT.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=3
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=4
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND. km2(i+1,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k .AND. km2(i+1,j+1).LT.k)
     *       nt3(i,j,k)=5
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=6
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=7
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j+1).GE.k .AND.
     *   km2(i-1,j-1).GE.k                        )
     *       nt3(i,j,k)=8
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=9
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j-1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i+1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j-1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i+1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *                           km2(i+1,j-1).GE.k)
     *       nt3(i,j,k)=10
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i-1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *   km2(i  ,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *   km2(i-1,j+1).LT.k .AND. km2(i  ,j-1).GE.k .AND.
     *   km2(i-1,j-1).GE.k)
     *       nt3(i,j,k)=11
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j  ).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i-1,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j+1).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).GE.k .AND. km2(i+1,j+1).LT.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i-1,j-1).LT.k .AND.
     *   km2(i-1,j+1).GE.k)
     *       nt3(i,j,k)=12
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *                           km2(I+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      IF(km2(i-1,j  ).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *                           km2(I+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      IF(km2(i-1,j+1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i  ,j-1).LT.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      IF(km2(i-1,j+1).LT.k .AND. km2(i+1,j  ).GE.k .AND.
     *   km2(i  ,j+1).GE.k .AND. km2(i+1,j-1).LT.k .AND.
     *                           km2(i+1,j+1).GE.k)
     *       nt3(i,j,k)=13
      END IF
 11   CONTINUE
      end do

      return
      end

      subroutine rwr2(a,b,il1,jl1,mg)
      dimension a(0:mg,0:il1,0:jl1), b(0:mg,0:il1,0:jl1)
      do j=0,jl1
      do i=0,il1
      do m=0,mg
      b(m,i,j)=a(m,i,j)
      end do
      end do
      end do
      return
      end

      subroutine Fdisc(aice,serv,iDay,il,jl)
*     Watcon and NDP Fortrans Lfort=1, Compaq and UNIX - Lfort=4
      Parameter (Lfort=1)
      Integer*2 dummy
*     The most important information writing.
*     File name is:  F xxxx 
*                      day 
*
      dimension Aice(il,jl), serv(il,jl)
      character*72 name
      nhice= Lfort*jl*jl
      name='a0000.grd'
      write(name(2:5),'(i4.4)') iday

      do j=1,jl
	do i=1,il
      Q= Aice(i,jl-j+1)
      serv(i,j)= Q
      zmin= min (zmin,Q)
      zmax= max (zmax,Q)
      end do
      end do
       dummy=
     *wr_grd(name,il,jl,-15.,15.,-15.,15.,zmin,zmax,serv,il)
      write(*,88) iday
88    format('<<<< Data are stored!!!', i4  ' >>>>')
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

