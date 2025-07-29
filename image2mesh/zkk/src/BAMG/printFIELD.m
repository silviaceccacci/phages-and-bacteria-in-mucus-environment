function [fullName] = printFIELD(fileName,scalar)

    fullName = [fileName '.bb'];

    stream = fopen(fullName,'w');
    
    nbsol = 1;
    nbv = length(scalar);
        
    fprintf(stream,['2 ' int2str(nbsol) ' ' int2str(nbv) ' 2\n']);
    fprintf(stream,'%f \n', scalar);
    
    fclose(stream);


end