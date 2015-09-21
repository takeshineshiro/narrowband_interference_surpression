clc;
clear  all  ;

xlLoadChipScopeData('wp.prn');
signal = dout;
signal = signal';
fs = 48.96e6 ;

figure(1);

fft_s  =  fft(signal);

abs_s  =  abs(fft_s).^2/length(fft_s);

db_s   =  10log10(abs_s);

length_s  = [0:length(fft_s)-1]*fs/length(fft_s);

abs_s   =  abs_s(1:length(length_s));

plot(length_s,abs_s);
xlabel('Hz');
ylabel('power spectrum');
title ('power');











