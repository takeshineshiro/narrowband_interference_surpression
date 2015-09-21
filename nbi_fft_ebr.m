%%%%%%%%%%simulation for nbi%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%GPS接收窄带干扰抑制频域滤波%%%%%%
  %%%author:wong   email:takeshineshiro@126.com%%%%
 %%%%%%%%%%%%%%%%%it  works %%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
 snr_h       = [-5:5];
 M           = 1 ;
 theta       = 4 ;
 
 roolfactor  = 0.50;
 delay       =3;
 delay_data  = delay*over_sample+1;
 j           = sqrt(-1);
 bit_error   = zeros(1,length(snr_h));
 pb          = zeros(1,length(snr_h));
 loop_num    = 1 ;
 
 for   snr_index = 1:length(snr_h)
   
     for  loop_index = 1:loop_num
         
      snr_db  = snr_h(snr_index);                                          % snr  here  fix noise
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
      
      
        f0     =  0.5e6 ;
        bw     =  1e6   ;                                                  % lfm     interference
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
        
        
        for i = 1:length(ss)
              if(ss(i)>0)
                   mod_signal(i)  = ss(i)+sqrt(-1)*0;
              else
                  mod_signal(i)   = ss(i)+sqrt(-1)*0;
            
            
              end
        end
        
       hn_sqrt   = rcosine(1,over_sample,'sqrt',0.70,delay);
       i_pulse   = conv(ss_os,hn_sqrt);
       i_pulse   = i_pulse(delay_data:end-delay_data+1);
       pulse_base= i_pulse;
       
       
       
       % hh= fir1(100,0.50);
       
       % pulse_base = filter(hh,1,ss_os);   % base impulse
       
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
      
      %%%%%%%%%%%%%rev%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       rev_ss  = ss_total  ;
       Max     = ceil(length(rev_ss)/fft_N);
       rev_ex  = [rev_ss,zeros(1,Max*fft_N-length(rev_ss))];               % rev  extend
       rev_len = length(rev_ex);
       
       delay_ex = [zeros(1,fft_N/2)  rev_ss(1:rev_len-fft_N/2)];           % delay fft_N/2
       
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
       
       delay_win    = delay_ex.*black_delay ;                              % add black window delay
       
       
       figure(3);
       plot(real(ss_total));
       xlabel('时间');
       ylabel('幅度');
       legend('原接收信号未加窗波形');
       title ('原接收信号未加窗波形');
       
       figure(4);
       plot(real(win_ss));
       xlabel('时间');
       ylabel('幅度');
       legend('原接收信号加窗后波形');
       title ('原接收信号加窗后波形');
       
       figure(5);
       plot(real(delay_win+win_ss));
       xlabel('时间');
       ylabel('幅度');
       legend('延迟N／2加窗后两信号相加波形');
       title('延迟N／2加窗后两信号相加波形');
       
       
       fft_ss     =  [];
       
       fft_buffer =  [];
       
       for i = 1:Max
           fft_temp   =  fft(win_ss((i-1)*fft_N+1:i*fft_N),fft_N);
           fft_buffer =  [fft_buffer fft_temp];
           abs_temp   =  abs(fft_temp);                                    % fft_sequence
           fft_ss     =  [fft_ss  abs_temp];
                  
       end
       
       
       fft_delay_ss     = [];
       fft_delay_buffer = [];
       
       
       
       for  i= 1:Max
            fft_delay        =  fft(delay_win((i-1)*fft_N+1:i*fft_N),fft_N);
            fft_delay_buffer =  [fft_delay_buffer  fft_delay];             % fft_sequence delay
            abs_delay        =  abs(fft_delay);
            fft_delay_ss     =  [fft_delay_ss abs_delay];
                  
       end
       
       
    
       %%%%%%%%%%%%follow  fre  domain   detection %%%%%%%%%%%%%%%
            %%% 1. assign 0   2. limit to threshold%%
                  %% th = 5*mu%%%%%
                  
            for  i  =  1:Max
                   mu_x(i)   =  sum (fft_ss((i-1)*fft_N+1:i*fft_N).^2)/fft_N;
                   th_x(i)   =  4*mu_x(i);                                    % threshold  cal
                                
            end
         
          for  i  = 1:Max
                  mu_y(i)   =  sum (fft_delay_ss((i-1)*fft_N+1:i*fft_N).^2)/fft_N ;                                 
                  th_y(i)   =  4*mu_y(i);                                   % threshold delay cal
              
          end
       
       % for i =1: length(fft_ss)
       %  index_x  =  floor ((i-1)/fft_N)  +1 ;
       %  if(fft_ss(i)^2 > th_x(index_x))
       %      fft_ss(i)     = 0 ;                     %  define   0              
       %      fft_ss(i)     = th_x(index_x);          % define th
       %   
       %      fft_buffer(i) = 0 ;                     % fft define to 0          
       %      fft_buffer(i) =(sqrt(th_x(index_x))/fft_ss(i))*fft_buffer(i);
       %   end                                       % fft amplitude  define  to  threshold                                 
       
       % end 
       
       
       
       for  i   =  1:length(fft_ss)
            index_x  =  floor((i-1)/fft_N)+1;
            if(fft_ss(i)^2  <= th_x(index_x))
                 fft_buffer(i)  = fft_buffer(i);
            elseif (fft_ss(i)^2  <= 8*th_x(index_x))
                  fft_buffer(i)  = (1/8)*fft_buffer(i);
           elseif(fft_ss(i)^2  <= 64*th_x(index_x))
                  fft_buffer(i)  = (1/64)*fft_buffer(i);                   % another paper define
           elseif(fft_ss(i)^2  <= 512*th_x(index_x))
                  fft_buffer(i)  = (1/512)*fft_buffer(i);
           else
                    
                  fft_buffer(i)  = (1/4096)*fft_buffer(i);
                      
                end
       end
                
     %  for   i  =  1: length(fft_delay_ss)    
     %      index_y  =  floor((i-1)/fft_N)  +1 ;
     %   if(fft_delay_ss(i)^2  > th_y(index_y))
     %         fft_delay_ss(i)     =  0;              % define  0
     %         fft_delay_ss(i)     =  th_y(index_y);  % define  th
     %         fft_delay_buffer(i) =  0;              %  fft define to 0
     %         fft_delay_buffer(i) =                  % fft amplitude
     %                                              define to threshold
     %         (sqrt(th_x(index_y))/fft_delay_ss(i))*fft_delay_buffer(i);
     
     %    end
     %    end
     
     
     for  i  = 1:length(fft_delay_ss)
            index_y  =   floor((i-1)/fft_N) +1 ;                           % another paper define  delay
            
             if (fft_delay_ss(i)^2   <= th_y(index_y))
                   fft_delay_buffer(i)  = fft_delay_buffer(i);
             elseif (fft_delay_ss(i)^2   <= 8*th_y(index_y))
                   fft_delay_buffer(i)  = (1/8)*fft_delay_buffer(i);
             elseif (fft_delay_ss(i)^2   <= 64*th_y(index_y))
                   fft_delay_buffer(i)  = (1/64)*fft_delay_buffer(i);
             elseif (fft_delay_ss(i)^2   <= 512*th_y(index_y))
                   fft_delay_buffer(i)  = (1/512)*fft_delay_buffer(i);
             else
                  fft_delay_buffer(i)  = (1/4096)*fft_delay_buffer(i);     
              
             end       
         
            
     end
       
       
    ifft_buffer   =  [];
    
    for  i  = 1:Max
        
         ifft_ss     =  ifft(fft_buffer((i-1)*fft_N+1:i*fft_N),fft_N);
         ifft_buffer =  [ifft_buffer  ifft_ss];                            % ifft
        
        
    end
       
     ifft_delay_buffer  =  [];
     
     for i  =  1:Max
         ifft_delay_ss     =  ifft(fft_delay_buffer((i-1)*fft_N+1:i*fft_N),fft_N);
         ifft_delay_buffer =  [ifft_delay_buffer  ifft_delay_ss];          % delay  ifft
                
     end
       
    ifft_delay_buffer  = ifft_delay_buffer(fft_N/2+1:end);                 % cut_delay_head
    
    append_cut_num  = length(ifft_buffer)-length(ifft_delay_buffer);
    ifft_delay_qq   =  [ifft_delay_buffer,zeros(1,append_cut_num)];
    
    
    ifft_merge  =  [];
    
    %for i = 1:Max
    %  ifft_cutoff        =  ifft_buffer((i-1)*fft_N+(fft_N/4)+1:(i-1)*fft_N+(fft_N/4)+(fft_N/2));
    %  ifft_delay_cutoff  =  ifft_delay_buffer((i-1)*fft_N+(fft_N/4)+1:(i-1)*fft_N+(fft_N/4)+(fft_N/2));
    %  ifft_merge         =  [ifft_merge  ifft_delay_cutoff  ifft_cutoff;
    % end
    
    
    ifft_merge_1  = [];
    
    % for  i  =  1:Max
    %       ifft_cutoff_1       =  ifft_buffer((i-1)*fft_N+1:(i-1)*fft_N+(fft_N/2));
    %       ifft_delay_cutoff_1 =  ifft_delay_buffer((i-1)*fft_N+1:(i-1)*fft_N+(fft_N/2));
    %       ifft_merge_1        =  [ifft_merge_1 ifft_cutoff_1  ifft_delay_cutoff_1  ];   
    %   end
    
    
    
    ifft_buffer_delay  =  [zeros(1,fft_N/2) ifft_buffer(1,length(ifft_buffer)-fft_N/2)];  % ifft delay
    
    %  ifft_total  =  ifft_buffer_delay+ifft_delay_buffer     %add
    
    ifft_total  =  ifft_buffer + ifft_delay_qq ;              % just add
    
    %   ifft_total   =   ifft_buffer ;
    
    figure(6);
    
    plot(real(ifft_total));
    xlabel('时间');
    ylabel('幅度');
    legend('经过抗窄带处理两信号相加后波形');
    title('经过抗窄带处理两信号相加后波形');
    
    
    figure(7);
    
    ss_ifft         =  fft(ifft_total,fft_N);
    ss_supre_spec   =  abs(ss_ifft).^2/fft_N;
    length_ss_supre = [0:length(ss_ifft)-1]*fs/length(ss_ifft); 
    
    ss_supre_spec   = ss_supre_spec(1:length(length_ss_supre));       % spresssion  spectrum
    ss_supre_spec   =  10*log10(ss_supre_spec);
    
    plot(length_ss_supre,ss_supre_spec);
    xlabel ('Hz');
    ylabel('power spectrum/db');
    legend('干扰抑制后信号功率普');
    title ('干扰抑制后信号功率普');
    
    
%%%%%%%%%%decode%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev_supres_s  =  ifft_total ;
rev_supres_n  = [0:length(rev_supres_s)-1];
rev_supres_m  = rev_supres_s.*exp(-j*2*pi*fc*rev_supres_n/fs);
%s_filter     = filter(hh,1,rev_supres_m);
i_demo_s      = real(rev_supres_m);
q_demo_s      = imag(rev_supres_m);
i_match       = conv(i_demo_s,hn_sqrt);
q_match       = conv(q_demo_s,hn_sqrt);
i_match       = i_match(delay_data:end-delay_data);
q_match       = q_match(delay_data:end-delay_data);
base_match    = i_match+j*q_match;

figure(8);
filter_ifft     = fft(base_match,fft_N);
filter_ddc_spec = abs(filter_ifft)/fft_N;
length_filter_s = [0:length(filter_ifft)-1]*fs/length(filter_ifft);
filter_ddc_spec = filter_ddc_spec(1:length(length_filter_s));
filter_ddc_spec =  10*log10(filter_ddc_spec);
plot(length_filter_s,filter_ddc_spec);
xlabel('Hz');
ylabel('power spectrum /db');
legend('DDC信号频谱');
title ('DDC信号频谱');




%%%%%%%%%decide%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hn_sqrt      = rcosine (1,over_sample,'sqrt',0.75,delay);  %equal method
% i_pulse_rev  = conv(rev_supres_m,hn_sqrt);
% i_pulse_rev  = i_pulse_rev(delay_data+1:end-delay_data);


s_desample  =  base_match(1:over_sample:end);

decode_buffer  =  [];

for  i =  1:N
    temp = sum(s_desample((i-1)*1023+1:i*1023).*pn_fix);
    decode_buffer  =  [decode_buffer temp];
    
    
end


for  i  = 1:length(decode_buffer)
      if(decode_buffer(i)>0)
           bit_out(i)  = 1 ;
      else
            bit_out(i)  = -1 ;
      end
     
end

bit_error(snr_index) = bit_error(snr_index)+sum(bit_out~=bit_gen);

     end
     
     
     % ideal error probability
     
     pb(snr_index)= 0.5*erfc(sqrt(10.^(snr_h(snr_index)/10)));
     
 end
 
 ber=  bit_error/N/loop_num;
 
 figure(9);
 
  semilogy(snr_h,pb,'b');
  
  hold  on ;
  grid  on ;
  
  semilogy(snr_h,ber,'r');
  
  title('simulation');
  legend('theoritical case','real simulation');
  xlabel('eb/n0(dB)');
  
  
  
  
  % for  i =1:length(bit_out)
  %    if(bit_out(i)-bit_gen(i))
  %       bit_error  = bit_error+1 ;
  %  end
  % end
  
  
  i_match_ze    =  s_desample<0;
  
  decode_bit    =  ds_demod(i_match_ze,pn_fix);
  
  
  
  bit_buf  = [];
  
  
  for   i  = 1:N
        s_dess   =  s_desample((i-1)*1023+1:i*1023).*pn_fix;
        s_sum    =  sum(s_dess);
        
        bit_buf  = [bit_buf  s_sum]  ;
      
      
  end
  
    real_bit   =  real(bit_buf);
    
    bit_decision = [];
    
    for  i  =1:length(real_bit)
          if(real_bit(i)>0)
               bit_decsion(i)   = 1 ;
          else
               bit_decsion(i)   =- 1 ;
          end
        
        
        
    end
    
    
    
    
    %%%follow  this  paper  to  define  threshold %%%%%%%%%
      % rev_awgn  =  ss_awgn.*exp(-j*2*pi*fc*n/fs);           %ddc
      
      rev_awgn   =  ss_awgn  ;
      
      awgn_ex    =  [rev_awgn, zeros(1,Max*fft_N-length(rev_awgn))];
      
      th_min     = compute_mean(awgn_ex,black_ex,Max,fft_N);       %  threshold  without interference
      
      fft_pre    = [];
      
      for  i = 2:Max
          
          ampti_temp  =   fft_ss((i-M-1)*fft_N+1:(i-1)*fft_N);          %  auto_threshold
          th(i-1)     =    th_min +theta*sum(ampti_temp)/(M*fft_N);     % th = th_min + theta*sum()/M*Nfft
          
          
      end
      
      flag = [];
      
      for i  =  fft_N+1  :length(fft_ss)
             th_index  =  floor((i-1)/fft_N);
             
             if(fft_ss(i)> th(th_index))                    % detection
                 flag =  [flag  i];
             end
             
      end
      
      
      
      flag_temp     =  rem(flag ,fft_N);
      
      unique_index  =  unique(flag_temp);                   % unique
      
      
      interf_fre    =  unique_index*fs/fft_N ;              % detect frequeny
      
      
      %%%%%%%inter  deal with %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
      
      %%%%%%%%%corresponding  define  the filter%%%%%%%%%%%
         num  = 1;
         flag_lower(num)   =  unique_index(1);
         flag_upper(num)   =  unique_index(1);
         
         f_cen    =  [];
         bw       =  [];                              % detect  num
         
         for i  = 2:length(unique_index)
               
              if (unique_index(i)-unique_index(i-1)>20)
                    num  =  num  +1 ;
                    flag_lower(num)    = unique_index(i); 
                    flag_upper(num-1)  = unique_index(i-1);
              elseif(i==length(unique_index))
                    flag_upper(num)    =  unique_index(i);
              else
                     num   = num;
              end
                         
             
         end
         
         
         bw_inter   =  (flag_upper-flag_lower)*fs/fft_N/2  ;       % interference  bw
         fc_inter   =  flag_lower*fs/fft_N +bw_inter;              % interference fc
         
         
         %%%%%%%%%%%define   the  filter %%%%%%%
          hn = [];
          
      
      
      
      
      
      
      
      
      
      
    
  
  
  
  
  
  
     






























    
    
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    
       
       
       
       
       
       
       
       
        
        
        
        
        
        
        
       
       
       
       
       
       
       
       
       
       
       
       
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
       
       
       
       
       
       
       
              
              
              
        
        
        
        
        
        
        
        
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
         
         
         
         
         
     end
     
     
     
     
     
     
 end
 
 
 
 
 
 
 