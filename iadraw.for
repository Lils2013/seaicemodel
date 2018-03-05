      program iadraw
c     Watcon and NDP Fortrans Lfort=1, Compaq and UNIX - Lfort=4
      Parameter (Lfort=1)
      parameter (iYear=2011,iMonth=12)
      open(31,file='iarea.dat',form='unformatted',
     *        access='direct',recl=Lfort)
      open(41,file='iext.dat',form='unformatted',
     *        access='direct',recl=Lfort)
      open(32,file='iasurf.dat')

      kkk= 12*(iYear-1948) +iMonth
      x = 1948.
      ya= 0.
      ye= 0.
      write(32,*) x, ya, ye

      do k=1,kkk
      read(31,rec= k) ya
      read(41,rec= k) ye
      write(32,*) 1948+(k)/12., ya, ye
      end do

      close(31)
      close(41)
      close(32)
      stop
      end
