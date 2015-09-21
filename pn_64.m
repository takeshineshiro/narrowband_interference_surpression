echo  on ;

connnection1 = [1,0,1,0,0];
connection2  = [1,1,1,0,1];
seq1         = ss_ml(connnection1);
seq2         = ss_ml(connection2);

L= 2^length(connnection1)-1;

for  shift_amount =0:L-1
    temp = [seq2(shift_amount+1:L), seq2(1:shift_amount)];
    gold_seq(shift_amount+1,:)= (seq1+temp)-floor((seq1+temp)/2);
    
    echo  off ;
      
end

echo  on ;

max_cross_corr =  0;


for i=  1:L-1
    
    for j= i+1:L
        
        c1 = 2*gold_seq(i,:)-1;
        c2 = 2*gold_seq(j,:)-1;
    for m= 0:L-1
        shift_c2  =  [c2(m+1:L), c2(1:m)];
        corr      =  abs(sum(c1.*shift_c2));
        
        if(corr>max_cross_corr)
             max_cross_corr = corr;
        end
        
        echo  off ;
        
        
    end
        
        
    
    
    
end