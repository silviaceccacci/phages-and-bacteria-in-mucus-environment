function [fullName] = printMTR(fileName,h)

    fullName = [fileName '.mtr'];

    stream = fopen(fullName,'w');

    if(size(h,2)==1)
        fprintf(stream,'%d %d\n',size(h,1),size(h,2));
        fprintf(stream,'%f\n', h');
    else
        fprintf(stream,'%d %d\n',size(h,1),size(h,2));
        fprintf(stream,'%f %f %f\n', h');        
    end    
    
    fclose(stream);


end