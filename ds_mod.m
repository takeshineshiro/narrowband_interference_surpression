function  [spread_base]   = ds_mod (bit_gen, pn_fix)
    % design  by  wong   takeshineshiro@126.com
    % this module  for base_signal  to spread
    % follow this algorithm : y(k)= pn_fix(k)xor x(i);
    % mod(sum(pn_fix(k)+x(i)),2) ==pn_fix(k) xor x(i);
    
    spread_base =  [];
    
    % for i  = 1:length(bit_gen)
    %  base_buffer_i    =  mod(pn_fix+bit_gen(i),2);
    %  spread_base      =  [spread_base ,base_buffer_i];
    % end
    
    
    length_pn   =  length(pn_fix);
    
    bit_buffer   = [];
    xor_buffer   =  [];
    
    for  i  =1:length(bit_gen)
        bit_buffer (1:length(pn_fix))  = bit_gen(i);
        xor_buffer     =  bit_buffer.*pn_fix ;
        spread_base    = [spread_base xor_buffer ];
        
        
        
    end
    
    
    
  
