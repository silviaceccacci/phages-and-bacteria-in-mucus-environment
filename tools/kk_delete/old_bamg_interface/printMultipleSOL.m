function [] = printMultipleSOL(fullName,u)

% 2 nbsol nbv 2
% ((Uij , ∀i ∈ {1, ..., nbsol}) , ∀j ∈ {1, ..., nbv})

    %fullName = [fileName '.bb'];

    stream = fopen(fullName,'w');
    
    nbsol = size(u,2);
    nbv   = size(u,1);

    fprintf(stream,'2 %d %d 2\n',nbsol,nbv);
    
    string_sol = '';
    for i=1:nbsol
        string_sol = [string_sol ' %f'];
    end
    string_sol = [string_sol '\n'];
    fprintf(stream,string_sol, u');
        
    fclose(stream);
    
end