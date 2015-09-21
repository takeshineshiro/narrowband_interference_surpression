function  [pn_data]  = prn_code(coeff_0,coeff_1)
   % design  by  wong   takeshineshiro@126.com
   % generate  the prn code
   
  state_0  = [1,0,1,1,0,1,1,1,0,1];
  state_1  = [1,1,1,1,1,1,1,1,1,1];
  
  N=2^length(coeff_0)-1;
  
  m0_dout = zeros(1,N);
  m1_dout = zeros(1,N);
  
  for i= 1:N
       m0_dout(i)   =  state_0(10);                                        %m0  generate
       new_s0       = mod(sum(coeff_0.*state_0),2);
       state_0(2:10)= state_0(1:9);
       state_0(1)   = new_s0;
         
  end
  
  for i= 1:N                     
      m1_dout(i)   = state_1(10);                                          %m1 generate
      new_s1       = mod(sum(coeff_1.*state_1),2);
      state_1(2:10)= state_1(1:9);
      state_1(1)   = new_s1;
        
  end
  
  for  i  =  1:N
      pn_data(i,:)   = xor(m0_dout,m1_dout);                                % xor
      
      % pn_data(i,:) = 2*pn_data(i,:)-1;    % change  to nz_code
      
      % m1_dout      =  shift(m1_dout,1,0); % shift
      
      m1_dout        =  circshift(m1_dout,[0 1]);
      
      
      
      
  end
  
  