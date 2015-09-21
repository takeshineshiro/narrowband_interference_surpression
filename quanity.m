%%%%%%% this for  the quanity  of  the nbi_interference %%%%%%%%%%%%%
clc;
clear  all ;
fb  = 50  ;                             %symbol_rate 
T   = 1*60;                             %time_second
num = T*fb;                             % num
chip_rate = 10.23e6 ;                   %chip_rate 
fs  = 61.38e6 ;                         %fs 
fc  = 15.48e6 ;                         %fc 
fft_N = 1024 ;                          %fft
over_sample = fs/chip_rate ;            %oversample
roolfactor  = 0.75;                     %factor
delay  = 3;                             %delay 
j= sqrt(-1);                            %hehe
snr_db = -5 ;                           %snr
isr  = [40:70];                         %isr
step = 5;                               %isr_step
inter_mod  =1 ;                         %interference type   1:cos_nbi  2:bpsk_nbi   3:sw_nbi
inter_num  =1 ;                         %interference  num  <4 
inter_fc   = fc;                        %interference  primitive fc
delay_data = delay*over_sample+1;       %cutoff
sigma  =1 ;                             %noise energy

snr   =  10^(snr_db/10);
eb    =  2*snr^2*sigma ;

e_chip = eb/1023;
e_over = e_chip/over_sample;  



%%%%%base_signal%%%%%%%%%%%%%%
 bit_gen  =  randn(1,num);
 bit_gen  =  bit_gen  >0 ;                                                 % source_bit_nz
 bit_gen  =  2*bit_gen -1;
 
 coeff_0  =  [1,0,1,0,0,0,0,0,0,1];                                        % m_0 coefficient  G1_poly= 1+x^3+x^10 gps_c/a
 coeff_1  =  [1,1,1,0,0,1,0,1,1,1];                                        % m_1 coefficient  G2_poly= 1+x^2+x^3+x^6+x^8+x^9+x^10  gps_c/a
 
 pn_code  = prn_code(coeff_0,coeff_1);  
 
 pn_fix   = pn_code (3,:);
 pn_fix   = 2*pn_fix-1 ;                                                   % pn_code_nz
 
 
 ss=  [];
 
 
 for i= 1:length(bit_gen)
      temp(1:1023)            =  bit_gen(i);
      ss((i-1)*1023+1:i*1023) =  temp.*pn_fix;                             % signal spread
         
 end
 
 for i =  1:length(ss)
     ss_os ((i-1)*over_sample+1:i*over_sample) = ss(i);                    % over_sample
         
 end
 
 % h_mod      =  modem.pskmod('M',2);
 % mod_signal = modulate(h_mod,ss) ;
 % %modulate
 
 for i= 1:length(ss)
     
     if(ss(i)>0)
         mod_signal(i)     =  ss(i)+sqrt(-1)*0  ;
     else
         mod_signal(i)     = ss(i)+sqrt(-1)*0   ;
     
     end
     
 end
 
 
 
 hn_sqrt    =  rcosine (1,over_sample,'sqrt',roolfactor,delay);
 i_pulse    =  conv(ss_os,hn_sqrt);
 i_pulse    =  i_pulse(delay_data:end-delay_data+1);
 pulse_base =  i_pulse ;
 
 n          = [0:length(pulse_base)-1];
 
 fc_ss      = pulse_base.*exp(j*2*pi*fc*n/fs);                             %fc 
 
 ss_awgn    = fc_ss+sigma*rand(1,length(fc_ss));                           % add  guass  white noise
 
 
 
 figure(2);
  
   fft_awgn      =  fft(ss_awgn ,fft_N);
   
   awgn_spec     =  abs(fft_awgn).^2/fft_N;                                % ds_spectrum
   
   length_awgn   =  [0:length(fft_awgn)-1]*fs/length(fft_awgn);
   
   awgn_spec     =  awgn_spec(1:length(length_awgn));
   
   awgn_spec     =   10*log10(awgn_spec);
   
   
   plot(length_awgn,awgn_spec);
   
   xlabel('Hz');
   ylabel('power spectrum(dB)');
   title ('power spectrum(dB)'));
   
   
   
   
   
   
   
   
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 

