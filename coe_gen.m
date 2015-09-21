hh= hamming(1024,'periodic');
hn= hh';
hn=hn*1024*16;
hn=fix(hn);

hm = [];


fid = fopen('/homw/wong/hamming.coe','w');


if(fid~=-1)
        fprintf(fid,'%s\n','memory_initialization_radix=16;');
        fprintf(fid,'%s\n','memory_initialization_vetor=');
        
        
        for i  = 1:length(hn)
              
            if(i==length(hn))
                fprintf(fid,'%4x;\n',hn(i));
            else
               fprintf(fid,'%4x,\n',hn(i)); 
                
            end
            
        end
        
        
    
    
end

fclose(fid);
