      program FramStrait
      Parameter (il=35, jl=49, kl=16)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (Lfort=1)

      Character*72 udat,i2dat,tdat,file

      Dimension v(0:il1,0:jl1,kl), km2(0:il1,0:jl1), z(kl)

        z(1)= 0.
        z(2)= -10.
        z(3)= -25.
        z(4)= -50.
        z(5)= -100.
        z(6)= -150.
        z(7)= -200.
        z(8)= -250.
        z(9)= -300.
        z(10)=-400.
        z(11)=-500.
        z(12)=-750.
        z(13)=-1000.
        z(14)=-2000.
        z(15)=-3000.
        z(16)=-4000.

      write(*,*) ' Enter year and month ->'
      read(*,*) ny, nmon

      udat= 'u000000'
      tdat= 't000000'
      i2dat='km2.dat'

      Write(udat(2:5),'(i4.4)') ny
      Write(tdat(2:5),'(i4.4)') ny
      Write(udat(6:7),'(i2.2)') nmon
      Write(tdat(6:7),'(i2.2)') nmon

      n2 = (il+2)*(jl+2)*Lfort
      n3 = n2*kl
      nm = il*jl*kl*Lfort

      Open ( Unit=49,File=i2dat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n2)

      Read(49,Rec=1) km2

*     Velocity

      Open ( Unit=31,File=udat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(31,Rec=2) v
      Close(31)

      file='vf000000.dat'
      Write(file(3:6),'(i4.4)') ny
      Write(file(7:8),'(i2.2)') nmon
      Open(12,FILE=file)
      
	Do 1 i=13,19
	x=110.*(i-13)
      Do 1 k=1,kl
      IF( km2(i,35).GE. k) THEN
      vm1= -v(i,35,k)
	write(12,*) x, z(k), vm1
	END IF
1     Continue
      close(12)

*     Temperature

      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=1) v
      Close(32)

      file='tf000000.dat'
      Write(file(3:6),'(i4.4)') ny
      Write(file(7:8),'(i2.2)') nmon
      Open(12,FILE=file)
      
	Do 2 i=13,19
	x=110.*(i-13)
      Do 2 k=1,kl
      IF( km2(i,35).GE. k) THEN
	write(12,*) x, z(k), v(i,35,k)
	END IF
2     Continue
      close(12)

*     Salinity

      Open ( Unit=32,File=tdat,Status='OLD',Access='DIRECT',
     *       Form='unformatted',Recl=n3)
      Read(32,Rec=2) v
      Close(32)

      file='sf000000.dat'
      Write(file(3:6),'(i4.4)') ny
      Write(file(7:8),'(i2.2)') nmon
      Open(12,FILE=file)
      
	Do 3 i=13,19
	x=110.*(i-13)
      Do 3 k=1,kl
      IF( km2(i,35).GE. k) THEN
	write(12,*) x, z(k), v(i,35,k)
	END IF
3     Continue
      close(12)

      stop
	end