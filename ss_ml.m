function [seq]  = ss_ml(connection)
   
   m= length(connection) ;
   L = 2^m-1;
   reg = [zeros(1,m-10) 1];
   seq(1) = reg(m);
   
   for i = 2:L
       new_reg(1) =  connection(1)*seq(i-1) ;
       for j= 2:m
           new_reg(j) = reg(j-1)+connection(j)*seq(i-1) ;
       end
        reg   =  new_reg  ;
        seq(i) = reg(m);
   end