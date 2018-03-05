      subroutine Redis

*     Version 28.01.16 based on the CICE4.1 model.
*     Exponential ice thickness distribution function.

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

!     !-----------------------------------------------------------------
!     ! b(h) = exp(-G(h)/astar)
!     ! apartic(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)]. 
!     ! The expression for apartic is found by integrating b(h)g(h)
!     ! between the category boundaries.
!     !-----------------------------------------------------------------


!     !----------------------------------------------------------------- 
!     ! The ridge ITD is a negative exponential: 
!     ! 
!     !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin 
!     ! 
!     ! where hrmin is the minimum thickness of ridging ice and 
!     ! hrexp is the e-folding thickness.
!     ! 
!     ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
!     ! That is, the minimum ridge thickness results from rafting,
!     !  unless the ice is thicker than maxraft.
!     !
!     ! Also, assume that hrexp = mu_rdg*sqrt(hi).
!     ! The parameter mu_rdg is tuned to give e-folding scales mostly
!     !  in the range 2-4 m as observed by upward-looking sonar.
!     !
!     ! Values of mu_rdg in the right column give ice strengths
!     !  roughly equal to values of Hstar in the left column
!     !  (within ~10 kN/m for typical ITDs):
!     !
!     !   Hstar     mu_rdg
!     !
!     !     25        3.0  m**0.5   30.cm**0.5
!     !     50        4.0
!     !     75        5.0
!     !    100        6.0
!     !----------------------------------------------------------------- 

      real mu_rdg
      parameter 
     &(Cf = 17.    ! ratio of ridging work to PE change in ridging 17.0
     &,Cs = 0.25   ! fraction of shear energy contrbtng to ridging: 0.-0.5
     &,fsnowrdg=0.5    ! snow fraction that survives in ridging: 0-0.5
     &,amaxraft= 1.e2  ! Maximum thickness of rafted ice 1 m by default.
     &,astar  = 0.05   ! e-folding scale for G(h) participation, ~0.05 
     &,mu_rdg= 30.     ! e-folding scale of ridged ice, ~30 (cm^0.5)
     &,maxiter= 2    ) ! number of iterations

      include 'Slo2.fi'

      real     
     & armvtot         ! area removed due to ridging and closing
     &,armvnet         ! net area removed 
     &,ardg(mgrad)     ! area of ice ridged
     &,anew(mgrad)     ! area of new ridge
     &,afrac           ! fraction of category area ridged
     &,virdg(mgrad)    ! ice volume of ridging ice
     &,vsrdg(mgrad)    ! snow volume of ridging ice
     &,msnow_ocn       ! mass of snow added to ocean (kg m-2)
     &,farea           ! fraction of new ridge area going to n2
     &,fvol            ! fraction of new ridge volume going to n2
     &,asurp(mgrad)    ! excess ridging area when ardg > ain
     &,Gsum(-1:mgrad)  ! Function g(h)
                           
      real     
     &athornn (0:mgrad) ! participation function;fraction of ridging/
                        ! closing associated w/ category n
      real     
     & hrminn(mgrad)  ! minimum ridge thickness
     &,hrexp (mgrad)  ! e-folding thickness
     &,krdgn (mgrad)  ! mean ridge thickness/thickness of ridging ice

      real
     &   aksumn            ! ratio of area removed to area ridged
     &,  ardgtot           ! total ice area ridged
     &,  anewtot           ! total area of new ridges
     &,  vrdgtot           ! total ice volume ridged


	include 'Tparm.fi'

	
	astari = 1./astar
      xtmp = 1.0/(1.0 - exp(-astari))

	Cp=0.5*g*roi*(row-roi)/row
	hz1 = HZ(1)
	
	STWR(:,:) =0.

   
*     Variables needed for both ridging and strength

      do j=1,jl
      do i=1,il
      IF( km2(i,j) .GT. 0) then


	LOOP_ITER: do iter=1,maxiter

	hrminn =0.
      krdgn  =1.
      hrexp  =0.


      Gsum(-1)= 0.
      do n = 0, mgrad
          Gsum(n) = Gsum(n-1) + aice(n,i,j)
      enddo
      if (Gsum(mgrad) .gt. 1.) then  ! renormalize Gsum
          do n = 0, mgrad
            Gsum(n) = Gsum(n) / Gsum(mgrad)
          enddo
      endif


      !-----------------------------------------------------------------
      ! b(h) = exp(-G(h)/astar)
      ! athornn(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)]. 
      ! The expression for athornn(n) is found by integrating b(h)g(h)
      ! between the category boundaries.
      !-----------------------------------------------------------------

      ! precompute exponential terms using Gsum as work array


         do n = -1, mgrad
         Gsum(n) = exp(-Gsum(n)*astari) * xtmp
         enddo                  

         do n = 0, mgrad
         athornn(n) = Gsum(n-1) - Gsum(n)
         enddo                  


        aksumn = athornn(0)

        dPE = 0.   ! change in potential energy due to ridging 

        do n = 1, mgrad
          if (athornn(n) .GT. 0. and. aice(n,i,j).GT.aimin) then
            hi = Hice(n,i,j)/aice(n,i,j)
            hi = max(hi,himin)
            hrminn(n)= min(2.*hi, hi+amaxraft)  ! amaxraft - max rafted
            hrexp(n) = mu_rdg * sqrt(hi)

            krdgn(n) = (hrminn(n) + hrexp(n)) / hi
            aksumn = aksumn + athornn(n) * (1. - 1./krdgn(n))

            h2rdg =      hrminn(n)*hrminn(n)
     &              + 2.*hrminn(n)*hrexp(n)
     &              + 2.*hrexp(n)*hrexp(n)
            dh2rdg = -hi*hi + h2rdg/krdgn(n)
            dPE = dPE + athornn(n) * dh2rdg

	     endif
        enddo         

*==================== Ice Strength ===================================

      dPE = Cp * dPE / aksumn  ! Cp = (g/2) * (rhow-rhoi) * (rhoi/rhow)

      Pice(i,j) = Cf * dPE     ! ice strength P, kg/s^2
c                                Cf accounts for frictional dissipation
c                                See Flato and Hibler (1995)

***	if(i.eq.17 .and. j.eq.24) write(*,*) aksumn, Pice(i,j)

*==================== Ridging ========================================

         
c     net ice/water area removed, from Hibler (1980)
c     armvnet = 0.5*(Delta_ice(i,j)-div_ice_tr(i,j)) * dt

c     net ice/water area removed, from Flato and Hibler (1995)

c      armvnet =(  0.5*Cs*(Delta_ice(i,j)-abs(div_ice(i,j))) - 
c     &              min(div_ice(i,j),0.)   ) * dt

      armvnet =(  0.5*Cs*(Delta_ice(i,j)-abs(div_ice(i,j))) - 
     &              min(div_ice_tr(i,j),0.)   ) * dt

      armvnet=max(armvnet,0.)

	if(div_ice_tr(i,j).LT.0.)
     &                  armvnet= max(armvnet,-div_ice_tr(i,j)*dt)

c     Make sure enough ice will be removed that the new fractional 
c     area of ice satisfies aice <= 1

	Atotal=0.
	do m=0,mgrad
	Atotal= Atotal+aice(m,i,j)
	end do
        if (atotal.gt.1.0 .and. armvnet.lt.atotal-1.0) then
          armvnet = atotal - 1.0
        endif

*     Total ice/water area removed
                 armvtot = armvnet / aksumn


c      if(i.eq.17.and.j.eq.24) then
c	write(*,*) armvnet, armvtot, aksumn
c	write(*,*) delta_ice(i,j), div_ice(i,j), div_ice_tr(i,j)
c	end if

c     area to be ridged in each category
        do n = 1, mgrad
          asurp(n) = 0.0       !  Initialize
        enddo

c     Initialize diagnostics
        ardgtot = 0.           ! area of ridging ice
        anewtot = 0.           ! area of new ridge
        vrdgtot = 0.           ! ridged ice volume

c     -------------------------------------------------------------
c      Compute the area, volume, and energy of ice ridging in each
c      category, along with the area of the resulting ridge.
c     -------------------------------------------------------------

        do n1 = 1, mgrad
          if (athornn(n1).GT.0. .and. aice(n1,i,j).GT.aimin) then
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
c          category n2.  
c      (3) Add area, volume, and energy to n2.
c     ---------------------------------------------------------------

        do n1 = 1, mgrad
          if (ardg(n1) .GT. 0.) then

c     Remove area, volume, and energy from ridging category.
            aice(n1,i,j) = aice(n1,i,j) - ardg(n1)
            hice(n1,i,j) = hice(n1,i,j) - virdg(n1)
            hsnow(n1,i,j)= hsnow(n1,i,j)- vsrdg(n1)

           hi1  = hrminn(n1)
           hexp = hrexp (n1)

c     Place part of the snow lost by ridging into the ocean.
            msnow_ocn = rosdry*vsrdg(n1)*(1.-fsnowrdg)
c     Fresh water flux change
            STWR(i,j)= STWR(i,j)+ msnow_ocn/row
            
c     Compute the fraction of ridged ice area and volume going to 
c     each thickness category, and transfer accordingly.
               
           do n2 = 1, mgrad

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category n2.
      !-----------------------------------------------------------------

      if (n2 .LT. mgrad) then

           if (hi1 .GE. hmax(n2)) then
                        farea = 0.
                        fvol  = 0.
           else
                        hL = max (hi1, hmax(n2-1))
                        hR = hmax(n2)

cccc	write(*,*)'hL,hR,n2',hL,hR,n2,hi1,hexp

                        expL = exp(-(hL-hi1)/hexp)
                        expR = exp(-(hR-hi1)/hexp)
           farea = expL - expR
           fvol  = ((hL + hexp)*expL- (hR + hexp)*expR) / (hi1 + hexp)
           endif

      else             ! n2 = mgrad

                     hL = max (hi1, hmax(n2-1))
                     expL = exp(-(hL-hi1)/hexp)
           farea = expL
           fvol  = (hL + hexp)*expL / (hi1 + hexp)

      endif            ! n2 < mgrad
                  
c     fraction of ridged ice in cat n1 that moves to cat n2

              aice (n2,i,j) = aice(n2,i,j) + farea*anew(n1)  
              hice (n2,i,j) = hice(n2,i,j) + fvol*virdg(n1)
              hsnow(n2,i,j)= hsnow(n2,i,j) + fvol*vsrdg(n1)*fsnowrdg 
          
		  enddo            ! loop over n2 (new ridges)  
		  
		  endif            ! cells with ridging for n1
        enddo              ! loop over n1 (ridging categories)



c     Small floes cut off.
      do mg=1,mgrad
      if(hice(mg,i,j).LT.Himin .OR. aice(mg,i,j).LT.aimin) then
      aice(mg,i,j)=0.
      hice(mg,i,j)=0.
      Tice(mg,i,j)=0.
	hsnow(mg,i,j)=0.
	Tsnow(mg,i,j)=0.
      end if
      end do


c     Check the total area
      ATotal= 0.
      do mg=1,mgrad
      ATotal= ATotal +aice(mg,i,j)
      !write(*,*) ATotal, "ATOTAL"
      end do

      if( ATotal .LE. 1.0) 	exit LOOP_ITER

	if(iter .EQ. maxiter .and. ATotal .GT. 1.0+aimin) then 
      write(*,*) 'max iterations exceeded', i,j, ATotal
      end if
     

      end do LOOP_ITER   ! iterations

c     Open water.
      aice(0,i,j)= 1.0
      do mg=1,mgrad
      aice(0,i,j)= aice(0,i,j)-aice(mg,i,j)
      end do
      if( aice(0,i,j) .LT. 0.) aice(0,i,j)=0.

      end if   ! km2 > 0
      enddo    ! j
      enddo    ! i

      return
      end 
