      subroutine tr(iYear,iMonth,u,v,km2,si,hz,il1,jl1,kl,ndisk)

*     Version for Methane Project

*     Watcom and NDP Fortrans Lfort=4, Compaq and UNIX - Lfort=1
      Parameter(Lfort=1)
      dimension u(0:il1,0:jl1,kl),v(0:il1,0:jl1,kl)
     *         ,km2(0:il1,0:jl1),hz(kl),si(0:jl1)

      pi=3.1415926
      hx=.64e-3*pi/180.

*     Calculation of transports through some passages

      kk= 12*(iYear-1948) +iMonth

*     1. Fram Strait

      FS= 0.
      j=35
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=13,18
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      FS= FS+ (c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +      +  c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      FS= FS*hx

*     2. Karskiye Vorota

      VK= 0.
      j=33
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=33,34
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      VK= VK+ (c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +      +  c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      VK= VK*hx

*     3. ZFI-Novaya Zemlya

      ZFINZ= 0.
      j=30
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=26,28
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      ZFINZ= ZFINZ+ (c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +            +  c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      ZFINZ= ZFINZ*hx

*     4. Norway-Spitsbergen

      SNA= 0.
      j=36
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=22,24
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SNa= SNa + (c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +         +  c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      i=25
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SNa= SNa + (c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do

c      SNA= SNA*hx

      SNB= 0.
      i=25
      do j=37,42
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SNb= SNb +((u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/6.
     +         + (u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      j=36
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SNb= SNb +((u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/6.)*dz1
      end do

c      SNB= SNB*hx

      SN= (SNA -SNB)*hx

*     5. Spitsbergen-ZFI

      SZFIa= 0.
      j=34
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=22,23
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZFIa= SZFIa +(c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +             + c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      i=24
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZFIa= SZFIa +(c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do

c      SZFIa= SZFIa*hx

      SZFIb= 0.
      i=24
      do j=31,33
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZFIb= SZFIb +((u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/6.
     +             + (u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      j=34
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZFIb= SZFIb +((u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k))/6.)*dz1
      end do

c      SZFIB= SZFIB*hx

      SZFI= (SZFIA+SZFIB)*hx

*     6. ZFI-Severnaya Zemlya

      SZZFIa= 0.
      i=25
      do j=25,28
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZZFIa= SZZFIa +((u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/6.
     +               + (u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      j=24
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZZFIa= SZZFIa +((u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/6.)*dz1
      end do

c      SZZFIa= SZZFIa*hx

      SZZFIb= 0.
      j=24
      c2=(2.*si(j)+si(j-1))/3.
      i=25
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SZZFIb= SZZFIb +(c2*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.)*dz1
      end do

c      SZZFIb= SZZFIb*hx
      SZZFI = (SZZFIa+SZZFIb)*hx

*     7. Norwegian Sea

      SNorw= 0.
      j=49
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=11,22
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      SNorw= SNorw+ (c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +      +  c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      SNorw= SNOrw*hx

*     8. Denmark Strait

      DenS= 0.
      i=7
      do j=44,46
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      DEnS= DenS +((u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))/6.
     +               + (u(i,j,k)+u(i,j-1,k)+u(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      DEnS= DenS*hx

*     9. Bering Passage

      Ber= 0.
      j=2
      c1=(2.*si(j)+si(j-1))/3.
      c2=(2.*si(j-1)+si(j))/3.
      do i=11,15
      klast=MIN(km2(i,j-1),km2(i+1,j-1),km2(i,j),km2(i+1,j))
	IF(klast .GT. 0) then
      do k=1,klast
      if(k.eq.1) dz1= 0.5*hz(1)
      if(k.gt.1.and.k.lt.klast) dz1= 0.5*(hz(k)+hz(k-1))
      if(k.eq.klast) dz1= 0.5*hz(k-1)

      Ber= Ber+ (c1*(v(i,j,k)+v(i+1,j,k)+v(i+1,j-1,k))/6.
     +      +  c2*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k))/6.)*dz1
      end do
	end if
      end do

      Ber= Ber*hx

      write(*,999) FS,VK,ZFINZ,-SN,SZFI,SZZFI,SNorw,DenS,Ber
999   format(' Transports(Sv):', 9f7.3)



c     Some balances: inflow
c     ---------------------
c     Kara Sea
      B1= SZZFI-VK-ZFINZ
c     Barents Sea
      B2= -SN+VK+SZFI+ZFINZ
c     GIN Sea
      B3= 0.9+0.7 +(SN+FS)
ccc      B3= SN+FS
      write(*,998) B1,B2,B3
998   format(' Kara Sea:',e10.4,' Bar.Sea:',e10.4,' GIN Sea:',e10.4)
cc      write(*,*)'  '


      if(ndisk.eq.1) then
      timeS=REAL(iYear-1948) +iMonth/12.
      open(10,file='passage.dat',
     *     form='unformatted',access='direct',
     *     recl=10*Lfort)
      write(10,rec=kk) timeS,FS,VK,ZFINZ,-SN,SZFI,SZZFI,SNorw,DenS,Ber
      close(10)
	end if

      return
      end
