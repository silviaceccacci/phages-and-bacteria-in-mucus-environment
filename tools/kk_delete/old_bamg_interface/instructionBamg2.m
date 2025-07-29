
%creates the instructions to Bamg

function [adaptativeStep]=instructionBamg2(bamg_exe,mesh_input,mesh_output,...
    metric_input,solution_input, pressure_input,combination_input,...
    data,adapt_variable,solution_output)
    %parameters BAMG
    switch data.iso
        case '1'
            text_iso = ' -iso ';
        case '0'
            text_iso = ' ';
        otherwise
            error('Wrong value for iso');
    end
    
    switch data.hmin
        case '0'
            text_hmin = ' ';
        otherwise
            text_hmin = [' -hmin ' data.hmin];
    end
    
    switch data.hmax
        case '0'
            text_hmax = ' ';
        otherwise
            text_hmax = [' -hmax ' data.hmax];
    end
    
    switch data.anisomax
        case '0'
            text_anisomax = ' ';
        otherwise
            text_anisomax = [' -anisomax ' data.anisomax];
    end
    
    %metrics BAMG
    if adapt_variable.comb == 1
        comanda_comb = [' -Mbb ' combination_input];
        comanda_d = ' ';
        comanda_u = ' ';
        comanda_p = ' ';
        
    elseif adapt_variable.comb == 0
        comanda_comb = ' ';
        switch adapt_variable.d
            case true
                comanda_d = [' -Mbb ' metric_input];
            case false
                comanda_d = ' ';
            otherwise
                error('d');
        end
    
        switch adapt_variable.u
            case true
                comanda_u = [' -Mbb ' solution_input];
            case false
                comanda_u = ' ';
            otherwise
                error('u');
        end
    
        switch adapt_variable.p
            case true
                comanda_p = [' -Mbb ' pressure_input];
            case false
                comanda_p = ' ';
            otherwise
                error('p');
        end
        
    else
        error('comb');
    end
    
    adaptativeStep = [bamg_exe ' -b ' mesh_input ...
                        ' -coef ' data.coef ' -ratio ' data.ratio ' -errg ' ...
                        data.errg ' -err ' data.err text_iso text_hmin ...
                        text_hmax text_anisomax ' -nbv ' data.nbv ...
                        comanda_comb comanda_d comanda_u comanda_p ...
                        ' -o '  mesh_output ...
                        ' -rbb ' solution_input ' -wbb ' solution_output ];

end