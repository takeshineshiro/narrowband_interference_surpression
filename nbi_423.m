%%%%%%simulation for nbi %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%GPS接收窄带干扰时域LMS滤波%%%%%%
  %%%author:wong  email:takeshineshiro@126.com
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  clc;
  clear  all ;
  
  symbol_rate = 1000;
  pn_rate     = 1.024e6 ;
  fc          = 2.048e6 ;
  fs          = 4.096e6 ;
  
 over_sample = fs/symbol_rate ;
 N           = 10 ;
 T           = N/symbol_rate ;
 fj          = fc-pn_rate/2 ;
 fft_N       = 4096 ;
 lamda       = 1  ;
 M           = 1 ;
 theta       = 4 ;
 
 j  = sqrt(-1);
 
 K = 100;                   % filter tap num
 wn=zeros(1,K);             % filter tap 
 error =  [];               % filter  error
 mu  = 0.0001;              % lms
 
 
     snr_db  =  0 ;                                                      % snr    here fix noise
 
      sigma   = 1 ;
      snr     = 10^(snr_db/10);
      eb      = 2*snr^2*sigma;
      e_chip  = eb/1023 ;
      e_over  = e_chip/over_sample;
      
      
      isr_db  = 40 ;
      isr     = 10^(isr_db/10);                                            % single  interference
      e_i     = isr*eb;
      a_coe   = sqrt(e_i);
      a_over  = sqrt(e_i/1023/over_sample);
      
      
        f0     = 1e6 ;
       bw      =  1e6   ;                                                  % lfm     interference
       delta_f =  bw/(2*T);
      
       fj_0    =  0.5e6 ;
       fj_1    =  1e6 ;
       aj_0    = sqrt(e_i/1023/over_sample/2);                             %two   interference
       aj_1    =  sqrt(e_i/1023/over_sample/2);
       
       
       fj_c0   = 0.5e6 ;
       fj_c1   =  1e6  ;
       fj_c2   =  1.3e6;
       
       aj_c0   = sqrt(e_i/1023/over_sample/3);
       aj_c1   = sqrt(e_i/1023/over_sample/3);                             % three  interference
       aj_c2   = sqrt(e_i/1023/over_sample/3);
       
       
           
       psk_rate = 0.256e6 ;                                                % psk  interference
       psk_over = fs/psk_rate; 
       fc_psk   = 1e6 ;
       psk_num  = N*1023*(over_sample/psk_over);
       
       psk_gen  = randn(1,psk_num);
       psk_gen  = psk_gen >0 ;
       psk_gen  = 2*psk_gen-1;
       psk_sam  = [];
       
        for i = 1:length(psk_gen)
         psk_sam((i-1)*psk_over+1:i*psk_over)  = psk_gen(i);  
           
       end
       psk_n    = [0:length(psk_sam)-1];
       psk_fc   = a_over*psk_sam.*exp(j*2*pi*fc_psk*psk_n/fs);
       
       
       
        %%%%%%%%%%%%%%%%%%%trans%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bit_gen   =  randn(1,N);
        bit_gen   =  bit_gen >0 ;
        bit_gen   =  2*bit_gen-1;                                          % source_bit_nz
        
        coeff_0   = [1,0,1,0,0,0,0,0,0,1];                                 % m_0 coefficient G1_poly=1+x^3+x^10 gps_c/a
        coeff_1   = [1,1,1,0,0,1,0,1,1,1];                                 % m_1 coefficient G2_poly=1+x^2+x^3+x^6+x^8+x^9+x^10 gps_c/a
        
        pn_code   = prn_code(coeff_0,coeff_1);
        pn_fix    = pn_code(3,:);
        pn_fix    = 2*pn_fix-1;                                            % pn_code_nz
        
        
        %spread_base = ds_mod(bit_gen ,pn_fix);
        
        %h_mod       = modem.pskmod('M',2);
        %mod_signal  = modulate(h_mod,spread_base);
        % ss         = real(mod_signal);
        
        
        ss = [];
        
         for i  =1:length(bit_gen)
             temp(1:1023)            = bit_gen(i);
             ss((i-1)*1023+1:i*1023) = temp.*pn_fix ;                     % signal spread
            
            
        end
        
         for i= 1:length(ss)
            
            ss_os((i-1)*over_sample+1:i*over_sample)  = sqrt(a_over)*ss(i);  % oversample
            
         end
         
         
         hh =  fir1(100,0.30);
         pulse_base  = filter(hh,1,ss_os);                                 % base impulse
         
         
                                                 % rcosine  coeffcient
      % rr     = rcosine(1,over_sample,'sqrt',roolfactor,delay);      
      % I_sent = rcosflt(real(mod_signal),1,over_sample,'filter',rr);
      % I_sent = I_sent(over_sample*delay+1:end-over_sample*delay)';
      
      
      % ss_os       = I_sent ;
      % pulse_base  = ss_os  ;
      
      % hn_sqrt    = rcosine(1,over_sample,'sqrt',0.75,delay);
      % i_pulse    = conv(ss_os,hn_sqrt);
      % i_pulse    = i_pulse(delay_data+1:end-delay_data);
      
      % ss_os      = i_pulse;
      
           n       =    [0:length(ss_os)-1];
      
      fc_ss   =  pulse_base.*exp(j*2*pi*fc*n/fs);                          % fc
      
      ss_awgn = fc_ss+sigma*randn(1,length(fc_ss));                        % add guass white noise
      
      interference  = a_over*cos(2*pi*fj*n/fs) ;                           % single interference
      
      interference_1 = a_over*cos(2*pi*f0*n/fs+2*pi*delta_f*(n/fs).^2);    % lfm  interference
      
      interference_2 = aj_0*cos(2*pi*fj_0*n/fs)+aj_1*cos(2*pi*fj_1*n/fs);  % two  interference
      
      interference_3 = aj_c0*cos(2*pi*fj_c0*n/fs)+aj_c1*cos(2*pi*fj_c1*n/fs)+aj_c2*cos(2*pi*fj_c2*n/fs);  % three interference
      
      interference_4 = psk_fc ;                                            % bpsk  interference
      
      
      ss_total      = ss_awgn + interference_3  ;                          % add   interference
      
      
        figure(1);
      
      fft_awgn     =  fft(ss_awgn ,fft_N);
      
      awgn_spec    = abs(fft_awgn).^2/fft_N;                               % ds   spectrum
      
      length_awgn  = [0:length(fft_awgn)-1]*fs/length(fft_awgn);
      
      awgn_spec    = awgn_spec(1:length(length_awgn));
      
      awgn_spec    = 10*log10(awgn_spec);
      
      
      plot(length_awgn,awgn_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      legend('未加干扰时信号功率普');
      title('未加干扰时信号功率普');
      
      
      figure(2);
      
      ss_fft    =  fft(ss_total,fft_N);
      ss_spec   = abs(ss_fft).^2/fft_N;
      length_ss = [0:length(ss_fft)-1]*fs/length(ss_fft) ;
      
      ss_spec    = ss_spec(1:length(length_ss));
      ss_spec    = 10*log10(ss_spec);
      
      plot(length_ss,ss_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      legend('加干扰时信号功率普');
      title ('加干扰时信号功率普');
      
      
      %%%%%%%%%%%%supression%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for  i =K:length(ss_total)
            
          s_in     = real(ss_total((i-K)+1:i));
          y_out(i) =  real(s_in)*wn';
          error(i) = real(ss_awgn(i))- y_out(i);               %lms
          wn       = wn+mu*error(i)*s_in;
               
          
      end
      
      
      
      mse = error(K:end).^2;                                    %mse
      
      
      
      figure(3);
      
      t_s  = [0:length(ss_total)-1]/fs ;
      plot(t_s,real(ss_total));
      
      xlabel('时间');
      legend('接收信号时域波形');
      title ('接收信号时域波形');
      
      
       figure(4);
      
      t_y  = [0:length(y_out)-1]/fs ;
      plot(t_y,y_out);
      
      xlabel('时间');
      legend('接收信号LMS滤波输出时域波形');
      title ('接收信号LMS滤波输出时域波形');
      
      
       figure(5);
      
      t_e  = [0:length(mse)-1]/fs ;
      plot(t_e,mse);
      
      xlabel('时间');
      legend('接收信号LMS滤波输出均方误差');
      title ('接收信号LMS滤波输出均方误差');
      
      
      hn  =  wn;
      
      lms_s  =  filter(hn,1,ss_total);      % filter
      
      figure(6);
      
      hn_fft    = fft(hn,fft_N);
      hn_spec   = abs(hn_fft).^2/fft_N; 
      length_hn = [0:length(hn_fft)-1]*fs/length(hn_fft);   % supression spectrum
      hn_spec   = hn_spec(1:length(length_hn));
      hn_spec   = 10*log10(hn_spec);
      
      plot(length_hn,hn_spec);
      xlabel('Hz');
      ylabel('power spectrum/dB');
      legend('LMS滤波器功率普');
      title ('LMS滤波器功率普');
      
      
      
      figure(7);
      
      lms_fft    =  fft(lms_s,fft_N);
      lms_spec   = abs(lms_fft).^2/fft_N;
      length_lms = [0:length(lms_fft)-1]*fs/length(lms_fft);  % supression spectrum
      lms_spec   = lms_spec(1:length(length_lms));
      lms_spec   = 10*log10(lms_spec);
      
      plot(length_lms,lms_spec);
      
      xlabel('Hz');
      ylabel('power spectrum/dB');
      legend('LMS滤波后信号功率普');
      title ('LMS滤波后信号功率普');
      
      %%%%%%%%%%%%%rev%%%%%%%%%%%%%%%%%%%%
        s_rev     = lms_s.*exp(-j*2*pi*fc*n/fs);
        
        rev_base  =  filter(hh,1,s_rev);               %ddc
        rev_base  =  real(rev_base);
        
        rev_des   = [];
        
        len_base  = length(rev_base)/over_sample ;
        
        
        for  i  =  1:len_base
            
            rev_des(i) = sum(rev_base((i-1)*over_sample+1:i*over_sample))/over_sample;
                    
            
        end
        
        
        % rev_des =  rev_base(1:over_sample:end);
        
        
        for i = 1:N
            
            de_spread(i)  = sum(rev_des((i-1)*1023+1:i*1023).*pn_fix);
            
            
        end
        
        
        bit_decode  =  [];
        
        
        
        
        
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
         
       
      
 