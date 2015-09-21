function [outregi] = shift(inregi,shiftr,shiftu)
  %designed by wong  takeshineshiro@126.com
  %inregi : vector  or matrix
  %shiftr : the amount of shift  to the  right 
  %shiftu : the amount of shift  to the  top 
  %outregi:  register output
  
  [h,v]   = size(inregi);
  outregi =  inregi;
  shiftr  =  rem(shiftr,v);
  shiftu  =  rem(shiftu,h);
  
  if (shiftr >0) 
       outregi(:,1:shiftr)   = inregi(:,v-shiftr+1:v);
       outregi(:,1+shiftr:v) = inregi(:,1:v-shiftr); 
       
  elseif shiftr <0 
      outregi(:,1:v+shiftr)  = inregi(:,1-shiftr:v);
      outregi(:,v+shiftr+1:v)= inregi(:,1:-shiftr);
  end
  
  
  inregi = outregi;
  
  if shiftu  >0 
       outregi(1:h-shiftu,:)   = inregi(1+shiftu:h,:);
       outregi(h-shiftu+1:h,:) = inregi(1:shiftu,:);
       
  elseif shiftu  <0 
       outregi(1:-shiftu,:)   =  inregi(h+shiftu+1:h,:);
       outregi(1-shiftu:h,:)  =  inregi(1:h+shiftu,:);
  end
  
  
       
  
  
  
  
  
  
  