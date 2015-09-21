function [th_min]  = compute_mean(awgn_ex,black_ex,Max,fft_N)
 % here  we  define   the mean  threshold  without interference
 
 win_awgn = awgn_ex.*black_ex;          % add   black_window
 
 fft_awgn  =  [];
 
 
 
 
 for   i  =  1:Max
      fft_temp_0   =  fft(win_awgn((i-1)*fft_N+1:i*fft_N),fft_N);
      abs_temp_0   =  abs(fft_temp_0);
      
      fft_awgn     =  [fft_awgn  abs_temp_0];
      
     
     
 end
 
 
 length_awgn   =  length(fft_awgn);
 
 th_min   =  sum(fft_awgn)/length_awgn;
 
 
 