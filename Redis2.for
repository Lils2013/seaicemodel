      Subroutine Rebin
*     Version 1.7 19.10.10.
*     Redistribution to ice thickness bins.
      include 'Slo2.fi'
      dimension qi(mgrad),qs(mgrad)
      integer ni,ns
	include 'Tparm.fi'

      ni=nlice
      ns=ni+nlsno
      do 1 j=1,jl
      do 1 i=1,il

      IF( km2(i,j) .GT. 0) then

c	if(i.eq.8.and.j.eq.43) write(*,*) 'A', Hice(:,i,j)

c     Open water.
      aice(0,i,j)= 1.0
      do mg=1,mgrad
      aice(0,i,j)= aice(0,i,j)-aice(mg,i,j)
      end do
      IF( (1.-Aice(0,i,j)).GT. Aimin) then

c     Initial Heat content of ice and snow
      TFreez= TFr(S(i,j,1),0.)
      do mg=1,mgrad
	if(hice(mg,i,j) .GT. Himin) then
      qi(mg)= (tice(mg,i,j)+Tfreez)*hice(mg,i,j)
      qs(mg)= (tice(mg,i,j)+tsnow(mg,i,j))*hsnow(mg,i,j)
	else
	qi(mg)=0.
	qs(mg)=0.
	end if
      end do

c     First gradation water is eliminated by convergion
c     in a case of "negative" open water concentration.
      if( aice(0,i,j) .LT. 0.) then
      aice(1,i,j)= aice(0,i,j)+aice(1,i,j)
      aice(0,i,j)= 0.
      end if

c     Concentration deficit is passed to next thicker category.
      do mg=1,(mgrad-1)
      if( aice(mg,i,j).LT.0.) then
      aice(mg+1,i,j)=aice(mg+1,i,j)+aice(mg,i,j)
      aice(mg,i,j)= 0.
      hice(mg+1,i,j)=hice(mg+1,i,j)+hice(mg,i,j)
      hice(mg,i,j)= 0.
      hsnow(mg+1,i,j)=hsnow(mg+1,i,j)+hsnow(mg,i,j)
      hsnow(mg,i,j)= 0.
      qi(mg+1)=qi(mg+1)+qi(mg)
      qi(mg)= 0.
      qs(mg+1)=qs(mg+1)+qs(mg)
      qs(mg)= 0.
      end if
      end do

c     Reestablishing thickness distribution.

c     To next thicker category.
      do mg=1,(mgrad-1)
      if(hice(mg,i,j).GT.hmax(mg)*aice(mg,i,j)) then
      
      if (i .eq. 14 .and. j .eq. 11) then
        write(*,*) "hice(mg+,i,j)",hice(mg,i,j),mg,i,j,k
        write(*,*) "aice(mg+,i,j)", aice(mg,i,j)
      end if
      aice(mg+1,i,j)=aice(mg+1,i,j)+aice(mg,i,j)
      aice(mg,i,j)= 0.
      hice(mg+1,i,j)=hice(mg+1,i,j)+hice(mg,i,j)
      hice(mg,i,j)= 0.
      hsnow(mg+1,i,j)=hsnow(mg+1,i,j)+hsnow(mg,i,j)
      hsnow(mg,i,j)= 0.
      qi(mg+1)=qi(mg+1)+qi(mg)
      qi(mg)= 0.
      qs(mg+1)=qs(mg+1)+qs(mg)
      qs(mg)= 0.
      do k=1,ni
          TiceFE(mg+1,i,j,k)=TiceFE(mg,i,j,k)
          if (i .eq. 14 .and. j .eq. 11) then
        write(*,*) "TiceFE(mg+1,i,j,k)",TiceFE(mg+1,i,j,k),mg+1,i,j,k
      end if
      enddo
      do k=ni+1,ns
          TsnowFE(mg+1,i,j,k-ni)=TsnowFE(mg,i,j,k-ni)
          if (i .eq. 14 .and. j .eq. 11) then
       write(*,*) "TsnowFE()",TsnowFE(mg+1,i,j,k-ni),mg+1,i,j,k
      end if
      enddo
      end if
      end do

c     To next thinner category.
      do mg=mgrad,2,-1
      if( hice(mg,i,j) .LT. hmax(mg-1)*aice(mg,i,j)) then
      
      if (i .eq. 14 .and. j .eq. 11) then
        write(*,*) "hice(mg-,i,j)",hice(mg,i,j),mg,i,j,k
      end if
      aice(mg-1,i,j)=aice(mg-1,i,j)+aice(mg,i,j)
      aice(mg,i,j)= 0.
      hice(mg-1,i,j)=hice(mg-1,i,j)+hice(mg,i,j)
      hice(mg,i,j)= 0.
      hsnow(mg-1,i,j)=hsnow(mg-1,i,j)+hsnow(mg,i,j)
      hsnow(mg,i,j)= 0.
      qi(mg-1)=qi(mg-1)+qi(mg)   ! Heat conservation law
      qi(mg)= 0.
      qs(mg-1)=qs(mg-1)+qs(mg)
      qs(mg)= 0.
      do k=1,ni
          TiceFE(mg-1,i,j,k)=TiceFE(mg,i,j,k)
          if (i .eq. 14 .and. j .eq. 11) then
        write(*,*) "TiceFE(mg-1,i,j,k)",TiceFE(mg-1,i,j,k),mg-1,i,j,k
      end if
      enddo
      do k=ni+1,ns
          TsnowFE(mg-1,i,j,k-ni)=TsnowFE(mg,i,j,k-ni)
          if (i .eq. 14 .and. j .eq. 11) then
        write(*,*) "TsnowFE()",TsnowFE(mg-1,i,j,k-ni),mg-1,i,j,k
      end if
      enddo
      end if
      end do

c     New ice and snow temperatures.
      do mg=1,mgrad
      if( hsnow(mg,i,j) .GT. Hsmin) then
      tsnow(mg,i,j)= qs(mg)/hsnow(mg,i,j)- tice(mg,i,j)
      else
      tsnow(mg,i,j)= 0.
      end if
      if( hice(mg,i,j) .GT. Himin) then
      tice(mg,i,j)= qi(mg)/hice(mg,i,j) - TFreez
      else
      tice(mg,i,j)= 0.
      end if
      end do

c     Small floes cut off.
      do mg=1,mgrad
      if(hice(mg,i,j).LT.Himin .OR. aice(mg,i,j).LT.aimin) then
      aice(mg,i,j)=0.
      hice(mg,i,j)=0.
      Tice(mg,i,j)=0.
	hsnow(mg,i,j)=0.
	Tsnow(mg,i,j)=0.
	end if
      if(hsnow(mg,i,j).LT.Hsmin) then
	hsnow(mg,i,j)=0.
	Tsnow(mg,i,j)=0.
      end if
      end do

c     Open water.
      aice(0,i,j)= 1.0
      do mg=1,mgrad
      aice(0,i,j)= aice(0,i,j)-aice(mg,i,j)
      end do
      if( aice(0,i,j) .LT. 0.) aice(0,i,j)=0.

      END IF   ! Sufficient Ice Compactness in the Node.
	
      end if ! Km2>0
1     continue 
      return
      end

