      subroutine Redis
*     Version 15.01.13 based on the CICE3.0-CICE4.1 models.

!---!-------------------------------------------------------------------
!---! Computes changes in the thickness distribution due to divergence
!---!  and shear.  
!---!
!---! From Flato and Hibler (1995), the rate of energy dissipation is
!---!     M = P * [Cs*0.5*(Delta-eps1) - min(eps1,0)]
!---! where P is the ice strength
!---!       Cs accounts for energy lost in shear without ridging 
!---!       eps1 = divergence rate
!---!       eps2 = shear rate
!---!       Delta = sqrt (eps1^2 + eps2^2/e^2)   
!---!       e = aspect ratio of elliptical yield curve
!---! See Flato and Hibler (1995) and Stern et al. (1995) for details.
!---!
!---! Note: The expression in brackets is the NET area of ice and water
!---! removed; i.e. the open water area removed, plus the area of ice
!---! ridged, minus the area of new ridges.  The NET area removed, armvnet, 
!---! is less than the TOTAL area removed, armvtot: 
!---!              armvnet = armvtot * aksumn, where aksumn <=1.
!---!------------------------------------------------------------------

      parameter
     &(Gstar=0.15       ! used to compute participation function
     &,Hstar=25.e2      ! determines mean thickness of ridged ice (cm): 25-100m.
     &,Cf = 17.         ! ratio of ridging work to PE change in ridging
     &,Cs = 0.25        ! fraction of shear energy contrbtng to ridging: 0.-0.5
     &,fsnowrdg=0.5     ! snow fraction that survives in ridging: 0-0.5
     &,Gstari=1./Gstar  ! 1/Gstar
     &,amaxraft= 1.e2 ) ! Maximum thickness of rafted ice 1 m by default.

      include 'slo2.fi'
      dimension qi(mgrad),qs(mgrad)

      real     
     & armvtot         ! area removed due to ridging and closing
     &,armvnet         ! net area removed 
     &,ardg(mgrad)     ! area of ice ridged
     &,anew(mgrad)     ! area of new ridge
     &,afrac           ! fraction of category area ridged
     &,virdg(mgrad)    ! ice volume of ridging ice
     &,vsrdg(mgrad)    ! snow volume of ridging ice
     &,eirdg(mgrad)    ! ice energy of ridging ice
     &,esrdg(mgrad)    ! snow energy of ridging ice
     &,msnow_ocn       ! mass of snow added to ocean (kg m-2)
     &,farea           ! fraction of new ridge area going to n2
     &,fvol            ! fraction of new ridge volume going to n2
     &,hL, hR          ! left and right limits of integration
     &,dhr             ! hrmax - hrmin
     &,dhr2            ! hrmax^2 - hrmin^2
     &,asurp(mgrad)    ! excess ridging area when ardg > ain
     &,Gsum(-1:mgrad)  ! Function g(h)
                           
      real     
     &athornn (0:mgrad) ! participation function;fraction of ridging/
                        ! closing associated w/ category n
      real     
     & hrminn(mgrad)  ! minimum ridge thickness
     &,hrmaxn(mgrad)  ! maximum ridge thickness
     &,krdgn(mgrad)   ! mean ridge thickness/thickness of ridging ice

      real
     &   aksumn            ! ratio of area removed to area ridged
     &,  ardgtot           ! total ice area ridged
     &,  anewtot           ! total area of new ridges
     &,  vrdgtot           ! total ice volume ridged

	include 'tparm.fi'

	Cp=0.5*g*roi*(row-roi)/row
	hz1 = HZ(1)

   
*       Variables needed for both ridging and strength

      do j=1,jl
      do i=1,il
      IF( km2(i,j) .GT. 0) then

      if (aice(0,i,j) .lt. Gstar) then  

c-----------------------------------------------------------------
c Compute the participation function athornn(n); this is analogous to
c a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
c
c               area lost from category n due to ridging/closing
c  athornn(n) = ------------------------------------------------
c                  total area lost due to ridging/closing
c
c Assume b(h) = (2/Gstar) * (1 - G(h)/Gstar). 
c The expressions for athornn are found by integrating b(h)g(h) between
c the category boundaries.
c-----------------------------------------------------------------
        Gsum(-1)= 0.
        do n = 0, mgrad
          Gsum(n) = Gsum(n-1) + aice(n,i,j)
        enddo
        if (Gsum(mgrad) .gt. 1.) then  ! renormalize Gsum
          do n = 0, mgrad
            Gsum(n) = Gsum(n) / Gsum(mgrad)
          enddo
        endif

        do n = 0, mgrad
          if (Gsum(n) .lt. Gstar) then
            athornn(n) = Gstari * (Gsum(n)-Gsum(n-1)) * 
     &              (2. - (Gsum(n-1)+Gsum(n))*Gstari)
          elseif (Gsum(n-1) .lt. Gstar) then
            athornn(n) = Gstari * (Gstar-Gsum(n-1)) * 
     &              (2. - (Gsum(n-1)+Gstar)*Gstari)
          else
            athornn(n) = 0.
          endif
        enddo
        
c      -----------------------------------------------------------------
c      Compute max and min ridged ice thickness for each ridging category.
c      Assume ice is uniformly distributed between hrmin and hrmax.
c      aksumn = net/total = weighted sum of net area loss terms:
c      total area removed = area of ice that ridges
c      net area removed = total - area of new ridges
c      -----------------------------------------------------------------
        aksumn = athornn(0)

        dPE = 0.   ! change in potential energy due to ridging 

        do n = 1, mgrad
          if (athornn(n) .GT. 0. .and. aice(n,i,j) .GT. aimin) then
            hi = Hice(n,i,j)/aice(n,i,j)
            hrminn(n)= min(2.*hi, hi+amaxraft)  ! amaxraft - max rafted ice
*                                        ! Hibler (1980) + modification CICE4.1
            hrmaxn(n)= 2.*sqrt(Hstar*hi) ! Hibler (1980)
            if (hrmaxn(n) .lt. hrminn(n)) hrmaxn(n) = hrminn(n) + himin
c                       This can happen for large hin and small Hstar
            krdgn(n) = 0.5 * (hrmaxn(n) + hrminn(n)) / hi
            aksumn = aksumn + athornn(n) * (1. - 1./krdgn(n))

            dPE = dPE - athornn(n)*(hi**2) ! PE loss from ridging ice
            dPE = dPE +(athornn(n)/(3.*krdgn(n))) * 
     &           (hrmaxn(n)**3 - hrminn(n)**3) / (hrmaxn(n)-hrminn(n))
c                                           PE gain from new ridge

          else                ! neglect ridging for this category
            hrminn(n) = 0.
            hrmaxn(n) = 0.
            krdgn(n)  = 1.
          endif
        enddo         

*==================== Ice Strength ===================================

      dPE = Cp * dPE / aksumn  ! Cp = (g/2) * (rhow-rhoi) * (rhoi/rhow)

      Pice(i,j) = Cf * dPE     ! ice strength P, kg/s^2
c                                Cf accounts for frictional dissipation
c                                See Flato and Hibler (1995)

c	if(i.eq.17 .and. j.eq.24) write(*,*) Pice(i,j)

*==================== Ridging ========================================

c     Initial Heat content of ice and snow
      TFreez= TFr(S(i,j,1),0.)
	hsi   = 0.
      do mg=1,mgrad
	if(hice(mg,i,j) .GT. Himin) then
c	hsi= hsi +
c     #    Aice(mg,i,j)*(roi*hice(mg,i,j)+rosdry*hsnow(mg,i,j))
	qi(mg)= (tice(mg,i,j)+Tfreez)*hice(mg,i,j)
      qs(mg)= (tice(mg,i,j)+tsnow(mg,i,j))*hsnow(mg,i,j)
	else
	qi(mg)=0.
	qs(mg)=0.
	end if
      end do

c	hsi= MIN(hsi/row,800.)
         
c       net ice/water area removed, from Hibler (1980)
c       armvnet = 0.5*(Delta_ice(i,j)-div_ice_tr(i,j)) * dt
c       net ice/water area removed, from Flato and Hibler (1995)
        armvnet =(  0.5*Cs*(Delta_ice(i,j)-abs(div_ice(i,j))) - 
     &              min(div_ice_tr(i,j),0.)   ) * dt

      armvnet=max(armvnet,0.)

c     Make sure enough ice will be removed that the new fractional 
c     area of ice satisfies aice <= 1
	Atotal=0.
	do m=0,mgrad !!!!
cccc	do m=1,mgrad
	Atotal= Atotal+aice(m,i,j)
	end do
        if (atotal.gt.1.0 .and. armvnet.lt.atotal-1.0) then
          armvnet = atotal - 1.0
        endif

        ! total ice/water area removed
        armvtot = armvnet / aksumn


c      if(i.eq.17.and.j.eq.24) then
c	write(*,*) armvtot, aksumn
c	write(*,*) delta_ice(i,j), div_ice(i,j), div_ice_tr(i,j)
c	end if

c       area to be ridged in each category
        do n = 1, mgrad
          asurp(n) = 0.0       !  Initialize
        enddo

c       Initialize diagnostics
        ardgtot = 0.           ! area of ridging ice
        anewtot = 0.           ! area of new ridge
        vrdgtot = 0.           ! ridged ice volume

c     -------------------------------------------------------------
c      Compute the area, volume, and energy of ice ridging in each
c      category, along with the area of the resulting ridge.
c     -------------------------------------------------------------

        do n1 = 1, mgrad
          if (athornn(n1) .GT. 0.) then
            ardg(n1) = armvtot * athornn(n1) + asurp(n1) 
                                          ! area of ridging ice
            if (ardg(n1) .gt. aice(n1,i,j)) then
              if (n1 .lt. mgrad) then
                asurp(n1+1) = (ardg(n1)-aice(n1,i,j))*(1.-1./krdgn(n1))
     &             / (1.-1./krdgn(n1+1)) ! increase ridging in next cat
              endif
              ardg(n1) = aice(n1,i,j)
            endif
            anew(n1) = ardg(n1) / krdgn(n1) !area of new ridge
            afrac = ardg(n1) / aice(n1,i,j) !fraction of cat n1 ridging 
            virdg(n1) = hice(n1,i,j) * afrac
            vsrdg(n1) = hsnow(n1,i,j) * afrac
            eirdg(n1) = qi(n1) * afrac
            esrdg(n1) = qs(n1) * afrac
            ardgtot = ardgtot + ardg(n1)
            anewtot = anewtot + anew(n1)
            vrdgtot = vrdgtot + virdg(n1)
          else                ! athornn(n1) = 0
            ardg(n1) = 0.
          endif
        enddo

c     ---------------------------------------------------------------
c      For each ridging category n1:
c      (1) Remove area, volume and energy from n1.
c      (2) Compute the fraction of the area and volume going to each
c          category n2.  Assume Hibler (1980) redistribution function.
c      (3) Add area, volume, and energy to n2.
c     ---------------------------------------------------------------

        do n1 = 1, mgrad
          if (ardg(n1) .gt. aimin) then   

c     Remove area, volume, and energy from ridging category.
            aice(n1,i,j) = aice(n1,i,j) - ardg(n1)
            hice(n1,i,j) = hice(n1,i,j) - virdg(n1)
            hsnow(n1,i,j)= hsnow(n1,i,j)- vsrdg(n1)
            qi(n1) = qi(n1) - eirdg(n1)
            qs(n1) = qs(n1) - esrdg(n1)

c     Place part of the snow lost by ridging into the ocean.
            msnow_ocn = rosdry*vsrdg(n1)*(1.-fsnowrdg)
c     Fresh water flux change
            PME(i,j)=PME(i,j) + msnow_ocn/row
c     Heat to melt the snow --> to chill water
            Qdelta= -Qsnow*msnow_ocn 
            T(i,j,1)= T(i,j,1) +2.*Qdelta/(hz1*row*Cpw)
            
c     Compute the fraction of ridged ice area and volume going to 
c     each thickness category, and transfer accordingly.
            dhr = hrmaxn(n1) - hrminn(n1)
            dhr2 = hrmaxn(n1)*hrmaxn(n1) - hrminn(n1)*hrminn(n1)
               
            do n2 = 1, mgrad
                  
              if (hrminn(n1).ge.hmax(n2) .or. 
     $           hrmaxn(n1).le.hmax(n2-1)) then
                hL = 0.
                hR = 0.
              else
                hL = max(hrminn(n1),hmax(n2-1))
                hR = min(hrmaxn(n1),hmax(n2))
              endif
               
              if (hR-hL .gt. aimin) then
                farea = (hR-hL) / dhr 
c     fraction of ridged ice in cat n1 that moves to cat n2
              aice(n2,i,j) = aice(n2,i,j) + farea*anew(n1)  

              fvol = (hR*hR - hL*hL) / dhr2 ! associated ice volume
              hice(n2,i,j) =hice(n2,i,j) +fvol*virdg(n1)
              hsnow(n2,i,j)=hsnow(n2,i,j)+fvol*vsrdg(n1)*fsnowrdg 
              qi(n2) = qi(n2) + fvol*eirdg(n1)
              qs(n2) = qs(n2) + fvol*esrdg(n1)*fsnowrdg
              endif         ! hR-hL > aimin

            enddo            ! loop over n2 (new ridges)            
          endif              ! ardg(n1) > aimin
        enddo                ! loop over n1 (ridging categories)

      endif                  ! ai0 < Gstar 


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
	hsnow(mg,i,j)=0.
      end if
      end do

c     Open water.
      aice(0,i,j)= 1.0
      do mg=1,mgrad
      aice(0,i,j)= aice(0,i,j)-aice(mg,i,j)
      end do
      if( aice(0,i,j) .LT. 0.) aice(0,i,j)=0.

      end if
      enddo
      enddo

      return
      end 
