      program fast_ice
	  parameter (il=35, jl=49, il1=il+1, jl1=jl+1, kl=16)
	  dimension ice_mask(0:il1,0:jl1)
	  dimension nt3(0:il1,0:jl1,kl), servgrd(51,51), h(il,jl)
      nnn=(il+2)*(jl+2)*kl
      nn=(il+2)*(jl+2)
      open (unit=11,file='nt.dat',status='old',access='direct',
     *      form='unformatted',recl=nnn)
      read(11,rec=1) nt3
	  close(11)

      call rd_grd('D:\metan\bottom\h.grd',51,51,servgrd)
      do i=1,il
      ii= i+9
	  do j=1,jl
	  h(i,j)=servgrd(ii,jl-j+1)  ! depth in [m]
	  end do
	  end do

	  do i=1,il
	  do j=1,jl
	  ice_mask(i,j)=0
cc	  IF( h(i,j).GE.25. .AND. nt3(i,j,1).LE.1) ice_mask(i,j)=1
	  IF( nt3(i,j,1).LE.1) ice_mask(i,j)=1
	  IF( nt3(i,j,1).EQ.0) ice_mask(i,j)=0
	  end do
	  end do

      open (unit=12,file='icemask.dat',status='old',access='direct',
     *      form='unformatted',recl=nn)
	  write(12, rec=1) ice_mask

	  close(12)

	write(*,111) ((ice_mask(i,j), i=1,35), j=1,49)
111   format (i2)

	  end
 
      subroutine rd_grd (file,Nx,Ny,DATA)
      CHARACTER*(*) file
      dimension DATA(Nx,Ny)
      CHARACTER*4 Z
      integer ii,jj

      OPEN(9,FILE=file,status='old')
      READ(9,'(a4)') Z
      READ(9,*) ii,jj
      READ(9,*) Xmin, Xmax
      READ(9,*) Ymin, Ymax
      READ(9,*) Zmin, Zmax
      DO j=1,Ny
      READ(9,*) (data(i,j),i=1,Nx)
      ENDDO
      close(9)
      return
      end

      
