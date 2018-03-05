      Program TF
	  write(*,*) 'Enter S and P'
	  read(*,*) S, p
	  
	  sqs=SQRT(S)
      Pn= 2.518052e-3 +
     +    ( -5.854586e-2 +2.298e-3*sqs -3.008634e-4*S
     +      +1.184586e-11*p                          )*S +
     +    (-7.002353e-4 +8.414961e-9*p)*p

	Pd= 1.0 +1.363248e-6 *S*S*sqs +
     +       (-3.849327e-5 +9.168654e-10*p)*p

      TFr= Pn/Pd ! Pure water

	TFr=TFr - 2.518052e-3 + 0.014286e-3*S ! Air saturated water

	write(*,*) 'Tf=', Tfr

	stop
	end
