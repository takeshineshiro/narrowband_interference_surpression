%%%%%%%%%%%%%%simulation for  nbi%%%%%%%%%%%%%%%%%%%%%%
  %%%%GPS接收窄带干扰抑制新方法-张兰时频结合%%%%%%
  %%%%author:wong   email :takeshineshiro@126.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clc;
 clear all ;
 
 symbol_rate = 1000;                                                       % 1k   bps
 pn_rate     = 1.024e6;
 fc          = 2.048e6;
 fs          = 4.096e6;
 over_sample = fs/symbol_rate ;
 N           = 500 ;
 T           = N/symbol_rate ;
 fj          = fc-pn_rate/2 ;
 fft_N       = 4096 ;
 lamda       = 1  ;
 M           = 1 ;
 theta       = 4 ;
 
  j           = sqrt(-1);
 
      snr_db  =  0 ;                                                       % snr    here fix noise
 
      sigma   = 1 ;
      snr     = 10^(snr_db/10);
      eb      = 2*snr^2*sigma;
      e_chip  = eb/1023 ;
      e_over  = e_chip/over_sample;
      
      isr_db  = 25 ;
      isr     = 10^(isr_db/10);                                            % single  interference
      e_i     = isr*eb;
      a_coe   = sqrt(e_i);
      a_over  = sqrt(e_i/1023/over_sample);
      
       f0      =  1e6 ;
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
        
         % h_mod   = modem.pskmod('M',2);
         % bit_psk = modulate(h_mod,bit_gen);
         
         
        
        coeff_0   = [1,0,1,0,0,0,0,0,0,1];                                 % m_0 coefficient G1_poly=1+x^3+x^10 gps_c/a
        coeff_1   = [1,1,1,0,0,1,0,1,1,1];                                 % m_1 coefficient G2_poly=1+x^2+x^3+x^6+x^8+x^9+x^10 gps_c/a
        
        pn_code   = prn_code(coeff_0,coeff_1);
        pn_fix    = pn_code(3,:);
        pn_fix    = 2*pn_fix-1;                                            % pn_code_nz
        
             ss = [];
        
        for i  =1:length(bit_gen)
             temp(1:1023)            = bit_gen(i);
             ss((i-1)*1023+1:i*1023) = temp.*pn_fix ;                     % signal spread
            
            
        end
        
         for i= 1:length(ss)
            
            ss_os((i-1)*over_sample+1:i*over_sample)  = sqrt(a_over)*ss(i);  % oversample
            
         end
         
          n          =  [0:length(ss_os)-1];
          
       fc_ss         =  pulse_base.*exp(j*2*pi*fc*n/fs);                    % fc
      
      ss_awgn        = fc_ss+sigma*randn(1,length(fc_ss));                  % add guass white noise
      
      interference   = a_over*cos(2*pi*fj*n/fs) ;                           % single interference
      
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
      
      %%%%%%%%%%rev%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % rev_ss  = ss_total.*exp(-j*2*pi*fc*n/fs);            % ddc
     
     rev_ss  = ss_total ;
     
      Max    = ceil(length(rev_ss)/fft_N);
      
      rev_ex = [rev_ss, zeros(1,Max*fft_N-length(rev_ss))];
      
      
      %%%%%%%%%%%%%%%%%%%%%supresssion%%%%%%%%%%%%%%%%%
        black_coe  = blackman(fft_N,'periodic');                           % black   window
        black_coe  = hamming(fft_N,'periodic');                            % hamming window
        black_coe  = black_coe';
        
        win_ss   = [];
        black_ex = [];
        
        
           for i= 1:Max
           
            black_ex  = [black_ex black_coe];                              % black value seq
            
        end
        
       black_delay  = [zeros(1,fft_N/2) black_ex(1:end-fft_N/2)];          % black value delay 
       
       
       win_ss       =  rev_ex.*black_ex ;                                  % add  black window
       
       
       fft_ss  = [];
       
       for   i  =  1:Max
            fft_temp    = fft(win_ss((i-1)*fft_N+1),fft_N);
            abs_temp    = abs(fft_temp);                                   % fft sequence
            fft_ss      = [fft_ss abs_temp];    
       end
       
        %%%%%%%follow  this paper  to define threshold%%%%%%%%%%
         % rev_awgn  =  ss_awgn.*exp(-j*2*pi*fc*n/fs);    % ddc
          
         rev_awgn   = ss_awgn ;
         
         awgn_ex    =  [rev_awgn,zeros(1,Max*fft_N-length(rev_awgn))];
         
         th_min     = compute_mean(awgn_ex,black_ex,Max,fft_N);            % threshold  without interference
         
         
         fft_pre   =  [];
         
         for i  = 2:Max
             
             ampti_temp    =  fft_ss((i-M-1)*fft_N+1:(i-1)*fft_N);         % auto threshold
             th(i-1)       =  th_min +theta*sum(ampti_temp)/(M*fft_N);     % th = th_min +theta*sum()/M*nfft
                       
             
         end
         
         flag  =  [];
         
         
         for i  =  fft_N+1 :length(fft_ss)
                  th_index   =  floor((i-1)/fft_N);
                if(fft_ss(i) > th(th_index))                               % detection
                     flag  = [flag  i];
                  
                end
             
         end
         
    flag_temp   = rem (flag,fft_N);
    
    unique_index   = unique(flag_temp);                                    % unique
    
    
   inter_fre       =  unique_index*fs/fft_N ;                              % detect frequency
   
   
   %%%%%%%%%%%%%%%corresponding  define  the filter  elments %%%%%%%%%%
      num  = 1;
      
      flag_lower(num)   =  unique_index(1);
      flag_upper(num)   =  unique_index(1);
      
      f_cen  =  [];
      bw     =  [];
      
      
      for   i  = 2:length(unique_index)
           
           if(unique_index(i)-unique_index(i-1) >20)
                 num  =   num   +1 ;
                 flag_lower(num)    =  unique_index(i) ;                   % detect num
                 flag_upper(num-1)  =  unique_index(i-1) ;
           elseif(i==length(unique_index))
                 flag_upper(num)   = unique_index(i);
           else
                  num  = num;
           end
          
          
          
      end
    
  
      bw_inter   =  (flag_upper-flag_lower)*fs/fft_N/2 ;                   % interference bw
      fc_inter   =   flag_lower*fs/fft_N+bw_inter ;                        % interference fc
   
      flag_low_fs  =  flag_lower/fft_N;
      flag_upper_fs =  flag_upper/fft_N;
      
      
   %%%%%%%%%%%%%%%%%%%%%%%define  the  fir/iir  filter%%%%%%%
      wn= [];
      
      for i =  1:num/2
          
          wn_temp  = [flag_low_fs(i) flag_upper_fs(i)];
          wn       = [wn wn_temp];
          
      end
      
      
      wn=  wn*2;
      
      
      hh =  fir1(5000,wn,'stop');
      
      figure(3);
      
      
      hh_fft   =  fft(hh,fft_N);
      hh_spec  =  abs(hh_fft).^2/fft_N ;
      length_hh=  [0:length(hh_fft)-1]*fs/length(hh_fft);                  % supression  spectrum
      
      hh_spec  =  hh_spec(1:length(length_hh));
      hh_spec  = 10*log10(hh_spec);
      
      
      plot(length_hh,hh_spec);
      
      xlabel('Hz');
      ylabel('power spectrum/db');
      legend('FIR滤波器功率普');
      title('FIR滤波器功率普');
      
      
      ss_filter  =   filter (hh,1,ss_total);
      
      
      figure(4);
      
      ss_filter_fft   = fft(ss_filter,fft_N);
      
      ss_filter_spec  = abs(ss_filter_fft).^2/fft_N;
      
      len_ss_filter   = [0:length(ss_filter_fft)-1]*fs/length(ss_filter_fft);
      
      ss_filter_spec  = ss_filter_spec(1:length(len_ss_filter));
      
      ss_filter_spec  = 10*log10(ss_filter_spec);
      
      
      plot(len_ss_filter,ss_filter_spec);
      
       xlabel('Hz');
      ylabel('power spectrum/db');
      legend('抗窄带处理FIR滤波后信号功率普');
      title('抗窄带处理FIR滤波后信号功率普');
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
       
      