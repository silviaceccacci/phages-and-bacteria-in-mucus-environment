function [fullName] = printSOL(fileName,u)

% 2 nbsol nbv 2
% ((Uij , ∀i ∈ {1, ..., nbsol}) , ∀j ∈ {1, ..., nbv})

    if(strcmp(fileName(end-2:end),'.bb'))
        fullName = fileName;
    else
        fullName = [fileName '.bb'];
    end
    
    stream = fopen(fullName,'w');
    
    nbsol = 1;
    nbv   = length(u);

    fprintf(stream,'2 %d %d 2\n',nbsol,nbv);
    fprintf(stream,'%f\n', u');
    
    fclose(stream);
    
end