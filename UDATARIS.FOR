      Program UdatRis

*     Программа для выдачи результатов расчетов скорости
*     на широте 75 для сравнения с данными наблюдений. 
*     Программа работает с выдачей подпрограммы FDisc.

*     После работы программы надо поработать программой GRID!

      Parameter (il=35, jl=49, kl=16)
      Parameter (il1=il+1, jl1=jl+1)
      Parameter (pi=3.14159265, overfl= 1.701410e+38)

      Character*72 udat,file

      Dimension u(0:il1,0:jl1,kl),v(0:il1,0:jl1,kl),km2(3),z(14)
      DATA km2/11,13,14/

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

      write(*,*) ' Enter year, month and day ->'
      read(*,*) nyear, month, nday

*     Формирование имени файла

      udat= 'u0000000'

      Write(udat(2:4),'(i3.3)') nyear
      Write(udat(5:6),'(i2.2)') month
      Write(udat(7:8),'(i2.2)') nday

      n2 = 4.*(il1+1)*(jl1+1)
      n3 = n2*kl

      Open ( Unit=31,File=udat,Status='OLD',Access='DIRECT',
     *       Form='Unformatted',Recl=n3)
      Read(31,Rec=1) u
      Read(31,Rec=2) v
      Close(31)

      file='us00.dat'
      Write(file(3:4),'(i2.2)') month
      Open(12,FILE=file)

      cosa=cos(pi*30./180.)
      sina=sin(pi*30./180.)

      Do 1 i=1,3
      x= 1000.*(i-1.)
      Do 1 k=1,14

      IF( km2(i) .GE. k) THEN
      vmod= v(i+12,39,k)*cosa - u(i+12,39,k)*sina
      Write(12,'(f8.0,1x,f8.0,e15.7,e17.7)') x,z(k),vmod
ccc      ELSE
ccc      Write(12,'(f8.0,1x,f8.0,e15.7,e17.7)') x,z(k),overfl
      END IF

1     Continue

      Stop
      End
