clc;
clear all ;

tp = load('home/wong/nbi_interference/coe/whs.txt'); 

tpm = tp';

sign_ss =  [];

for i = 1:length(tpm)
    if(tpm(i)>260)
         sign_ss(i)  =  -(1024-tpm(i));
    else
        sign_ss(i)   =  tpm(i);
    end
    
end

  fs = 62e6 ;
  
  figure (1);
  
  fft_s   =  fft(sign_ss);
  abs_s   =  abs(fft_s).^2/length(fft_s);
   db_s   =  10*log10(abs_s) ;
 length_s =  [0:length(fft_s)-1]*fs/length(fft_s);
  ans_s   =  abs_s(1:length(length_s));
  
  
  plot(length_s,abs_s) ;
  
  xlabel ('Hz');
  ylabel ('power_spectrum');
  title  ('power_spectrum');
  
  
  
  figure(2);
  
  plot(length_s,db_s);
  xlabel('Hz');
  ylabel('power_spectrum(dB)');
  title ('power_spectrum(dB)');
  
  
  
  
  
    

