      program passdraw
	Parameter (Lfort=1)
      parameter (lenerg=24*30, dthr= 1.*(lenerg), dt=dthr*3600.)
      dimension y(10)

      nnn=10*Lfort

      open(31,file='passage.dat',form='unformatted',
     *        access='direct',recl=nnn)
      open(32,file='passgraf.dat')

c      open(33,file='time.dat' )
c      read(33,*) time
c      close(33)

      n=  66*12 !!!+1   + INT(time/dt)
      write(*,*) 'n=', n

      do 1 i=1,n
      read(31,rec= i) y
ccc      y(1)= y(1)/(24.*3600.*365.)
      write(32,100) y
1     continue

100   format(10f9.3)

      close(31)
      close(32)
      stop
      end
