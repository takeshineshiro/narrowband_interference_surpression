function  [p] = ss_pe (snr_db,lc,a,w0)
  snr=10^(snr_db/10);
  sigma  =1 ;
  Eb     = 2*sigma^2*snr;
  E_chip = Eb/lc;
  N      = 10000;
  num_of_err = 0 ;
  for i = 1:N
     temp  = rand ;
       if (temp<0.5)
           data  =  -1 ;
       else
           data  = 1 ;
       end
      
       for j=  1:lc
           temp = rand ;
           if (temp<0.5)
               pn_seq(j) = -1 ;
           else
               pn_seq(j) =  1 ;
               
           end
       end
       
       trans_sig  = sqrt(E_chip)*repeat_data.*pn_seq ;
       noise  =  sigma*randn(1,lc);
       n= (i-1)*lc+1 :i*lc ;
       interference = a*sin(w0*n);
       
       rec_sig  = trans_sig +niose +interference ;
       
       temp = rec_sig.*pn_seq ;
       
       decision_var  =  sum(temp);
        
       if (decision_var <0)
            decision   =  -1 ;
       else
           decision    =  1  ;
       end
       
       if (decision ~=data)
           num_of_err  = num_of_err  +1 ;
           ends
       
       
      
      



       end

       
       p=num_of_err/N ;
       