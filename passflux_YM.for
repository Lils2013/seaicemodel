      program passdraw
	Parameter (Lfort=1)
      parameter (lenerg=15*12, dthr= 2.*(lenerg), dt=dthr*3600.)
      dimension y(7), ym(7)

	n=55*12

      nnn=Lfort*7
      open(31,file='passage.dat',form='unformatted',
     *        access='direct',recl=nnn)
      open(32,file='passgraf.dat')
      open(33,file='passgraf_YM.dat')

c      open(33,file='time.dat' )
c      read(33,*) time
c      close(33)

c      n= INT(time/dt)
      write(*,*) 'n=', n

	do ky=1,55

	ym(1)= (ky)
      do j=2,7
	ym(j)=0.
	end do

      do 1 month=1,12
	i= 12*(ky-1)+month
      read(31,rec= i) y

      do j=2,7
	ym(j)= ym(j) +y(j)/12.
	end do

cc      y(1)= y(1)/(24.*3600.*365.)
      write(32,100) y
1     continue

      write(33,100) ym
	end do

100   format(7f9.5)

      close(31)
      close(32)
      close(33)
      stop
      end
