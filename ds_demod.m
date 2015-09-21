function  [decode_bit]  = ds_demod(rev_data, pn_fix)
  %decode the  spread  code
  
  decode_bit = [];
  decode_len = length(rev_data)/length(pn_fix);
  
  for i  = 1:decode_len
       
      decode_bit(i)   = mod(sum(rev_data((i-1)*1023+1:1023*i)-1.*pn_fix),2);
      
  end