      subroutine cascading
      
*     Version 25.01.2016      
      
* Computation of source terms due to bottom slope cascading
* Campin J.M. and H. Goosse (1999). A parameterization of density 
* driven downsloping flow for coarse resolution model in z-coordinate. 
* Tellus 51A,412-430      

      include 'Slo2.fi'

      row = 1.020     ! Water density, g/cm3

      do j=1,jl
      S0= si(j)
      do i=1,il
      
      kb= km2(i,j)
      
      if(kb .GT. 0) then
      
c      if(nt3(i,j,1).GT.0 .and. nt3(i+1,j,1).GT.0 .and.
c     *   nt3(i-1,j,1).GT.0  .and. nt3(i,j-1,1).GT.0 .and.
c     *   nt3(i,j+1,1).GT.0 ) then 
            
      if(km2(i+1,j) .gt. kb) then
         if( ro(i,j,kb) .GT. ro(i+1,j,kb)) then
         
c        find the neutral bouyancy level
         kn= km2(i+1,j)
         do k= kb+1, km2(i+1,j)
	   ppp= 1.e-5* g*row*z(k)               
         ro_new= 1.e-3*REAL(sigma_t(t(i,j,kb),s(i,j,kb), ppp))
         if(ro_new.LT.ro(i+1,j,k) .and. ro_new.GE.ro(i+1,j,k-1)) kn=k
         end do
      
c        compute the slope velocity
         deltaro = (ro(i,j,kb)-ro(i+1,j,kb))/row
         uslope= 0.5e4*g*deltaro*(z(km2(i+1,j))-z(kb))/(hx*S0)
         
c        Volume transport                  
         Vcascade = uslope * hz(kb-1)*hy
         
c        Compute sources
         QTc(i,j,kb)= QTc(i,j,kb)+ Vcascade*(T(i+1,j,kb)-T(i,j,kb))
         QSc(i,j,kb)= QSc(i,j,kb)+ Vcascade*(S(i+1,j,kb)-S(i,j,kb))
                         
         QTc(i+1,j,kn)= QTc(i+1,j,kn)+ Vcascade*(T(i,j,kb)-T(i+1,j,kn))
         QSc(i+1,j,kn)= QSc(i+1,j,kn)+ Vcascade*(S(i,j,kb)-S(i+1,j,kn))
         
         do k=kn-1,kb,-1       
         QTc(i+1,j,k)= QTc(i+1,j,k)+ Vcascade*(T(i+1,j,k+1)-T(i+1,j,k))
         QSc(i+1,j,k)= QSc(i+1,j,k)+ Vcascade*(S(i+1,j,k+1)-S(i+1,j,k))
         end do
         
        
c         Qmean= QTc(i,i,kb)
c         do k=kb,kn
c         Qmean = Qmean +QTc(i+1,j,k)
c         end do
c         write(*,*) 'C1:', i,j,kb,kn, Qmean

         end if   ! density i+1,j
      end if   ! slope i+1,j    
      
      
      if(km2(i-1,j) .gt. kb) then
         if( ro(i,j,kb) .GT. ro(i-1,j,kb)) then
         
c        find the neutral bouyancy level
         kn= km2(i-1,j)
         do k= kb+1, km2(i-1,j)
	   ppp= 1.e-5* g*row*z(k)             
         ro_new= 1.e-3*REAL(sigma_t(t(i,j,kb),s(i,j,kb), ppp))
         if(ro_new.LT.ro(i-1,j,k) .and. ro_new.GE.ro(i-1,j,k-1)) kn= k
         end do
      
c        compute the slope velocity
         deltaro = (ro(i,j,kb)-ro(i-1,j,kb))/row
         uslope= 0.5e4*g*deltaro*(z(km2(i-1,j))-z(kb))/(hx*S0)
         
c        Volume transport                  
         Vcascade = uslope * hz(kb-1)*hy
         
c        Compute sources
         QTc(i,j,kb)= QTc(i,j,kb)+ Vcascade*(T(i-1,j,kb)-T(i,j,kb))
         QSc(i,j,kb)= QSc(i,j,kb)+ Vcascade*(S(i-1,j,kb)-S(i,j,kb))
                
         QTc(i-1,j,kn)= QTc(i-1,j,kn)+ Vcascade*(T(i,j,kb)-T(i-1,j,kn))
         QSc(i-1,j,kn)= QSc(i-1,j,kn)+ Vcascade*(S(i,j,kb)-S(i-1,j,kn))
         
         do k=kn-1,kb,-1       
         QTc(i-1,j,k)= QTc(i-1,j,k)+ Vcascade*(T(i-1,j,k+1)-T(i-1,j,k))
         QSc(i-1,j,k)= QSc(i-1,j,k)+ Vcascade*(S(i-1,j,k+1)-S(i-1,j,k))
         end do
         
ccc        write(*,*) 'cascading', i,j,deltaro, uslope, Vcascade, kb,kn

         end if   ! density i-1,j
      end if   ! slope i-1,j
      
      if(km2(i,j+1) .gt. kb) then
         if( ro(i,j,kb) .GT. ro(i,j+1,kb)) then
         
c        find the neutral bouyancy level
         kn= km2(i,j+1)
         do k= kb+1, km2(i,j+1)
	   ppp= 1.e-5* g*row*z(k)               
         ro_new= 1.e-3*REAL(sigma_t(t(i,j,kb),s(i,j,kb), ppp))
         if(ro_new.LT.ro(i,j+1,k) .and. ro_new.GE.ro(i,j+1,k-1)) kn= k
         end do
      
c        compute the slope velocity
         deltaro = (ro(i,j,kb)-ro(i,j+1,kb))/row
         uslope= 0.5e4*g*deltaro*(z(km2(i,j+1))-z(kb))/hy
         
c        Volume transport                  
         Vcascade = uslope * hz(kb-1)*hx*S0
         
c        Compute sources
         QTc(i,j,kb)= QTc(i,j,kb)+ Vcascade*(T(i,j+1,kb)-T(i,j,kb))
         QSc(i,j,kb)= QSc(i,j,kb)+ Vcascade*(S(i,j+1,kb)-S(i,j,kb))
                
         QTc(i,j+1,kn)= QTc(i,j+1,kn)+ Vcascade*(T(i,j,kb)-T(i,j+1,kn))
         QSc(i,j+1,kn)= QSc(i,j+1,kn)+ Vcascade*(S(i,j,kb)-S(i,j+1,kn))
         
         do k=kn-1,kb,-1       
         QTc(i,j+1,k)= QTc(i,j+1,k)+ Vcascade*(T(i,j+1,k+1)-T(i,j+1,k))
         QSc(i,j+1,k)= QSc(i,j+1,k)+ Vcascade*(S(i,j+1,k+1)-S(i,j+1,k))
         end do

         end if   ! density i,j+1
      end if   ! slope i,j+1
            
      if(km2(i,j-1) .gt. kb) then
         if( ro(i,j,kb) .GT. ro(i,j-1,kb)) then
         
c        find the neutral bouyancy level
         kn= km2(i,j-1)
         do k= kb+1, km2(i,j-1)
	   ppp= 1.e-5* g*row*z(k)            
         ro_new= 1.e-3*REAL(sigma_t(t(i,j,kb),s(i,j,kb), ppp))
         if(ro_new.LT.ro(i,j-1,k) .and. ro_new.GE.ro(i,j-1,k-1)) kn= k
         end do
      
c        compute the slope velocity
         deltaro = (ro(i,j,kb)-ro(i,j-1,kb))/row
         uslope= 0.5e4*g*deltaro*(z(km2(i,j-1))-z(kb))/hy
         
c        Volume transport                  
         Vcascade = uslope * hz(kb-1)*hx*S0
         
c        Compute sources
         QTc(i,j,kb)= QTc(i,j,kb)+ Vcascade*(T(i,j-1,kb)-T(i,j,kb))
         QSc(i,j,kb)= QSc(i,j,kb)+ Vcascade*(S(i,j-1,kb)-S(i,j,kb))
                
         QTc(i,j-1,kn)= QTc(i,j-1,kn)+ Vcascade*(T(i,j,kb)-T(i,j-1,kn))
         QSc(i,j-1,kn)= QSc(i,j-1,kn)+ Vcascade*(S(i,j,kb)-S(i,j-1,kn))
         
         do k=kn-1,kb,-1       
         QTc(i,j-1,k)= QTc(i,j-1,k)+ Vcascade*(T(i,j-1,k+1)-T(i,j-1,k))
         QSc(i,j-1,k)= QSc(i,j-1,k)+ Vcascade*(S(i,j-1,k+1)-S(i,j-1,k))
         end do

         end if   ! density i,j-1
      end if   ! slope i,j-1    
         
      end if   ! model area
      
      end do
      end do
      
      return
      end
      