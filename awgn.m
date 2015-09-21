function  [i_wgn, q_wgn]  = awgn(i_pulse,q_pulse,sigma)
  % designed by  wong    takeshineshiro@126.com
  len  = length(i_pulse);
  
  i_wgn  = i_pulse+sigma*randn(1,len);
  q_wgn  = q_pulse+sigma*randn(1,len);
  
  