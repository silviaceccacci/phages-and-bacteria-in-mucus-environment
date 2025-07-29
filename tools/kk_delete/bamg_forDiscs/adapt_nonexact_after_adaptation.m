function [Xh,Th,u,p,ux,uy] = adapt_nonexact_after_adaptation(X,T,data,adapt_variable,...
    metric_input,solution_input, pressure_input,combination_input,...
    ux_input,uy_input)

    [bamg_exe]=get_bamg_exe();

    fileName_mesh_input = './temps/temp_mesh_input';
    mesh.X = X';
    mesh.T = T;
    mesh_input = printMSH(fileName_mesh_input,mesh);
    
    fileName_mesh_output = './temps/temp_mesh_output';
    mesh_output = './temps/temp_mesh_output.msh';

    fn_sol_output   = './temps/temp_solution_output';
    solution_output = './temps/temp_solution_output.bb';
    fn_press_output = './temps/temp_pressure_output';
    pressure_output = './temps/temp_pressure_output.bb';

    adaptativeStep = instructionBamg2(bamg_exe,mesh_input,mesh_output,...
        metric_input,solution_input,pressure_input,combination_input,...
        data,adapt_variable,solution_output);

    adaptativeStep = [adaptativeStep ' > ./temps/out_bamg_adapt.txt'];%PETA NO SE XQ
    disp(adaptativeStep)
    eval(adaptativeStep) %aquest funcio canvia directament arxius
    
    %interpolate addinitioal variables
    if(isempty(ux_input))
        adaptativeStep = [ bamg_exe  ... 
            ' -b ' mesh_input        ...
            ' -r ' mesh_output       ...
            ' -rbb ' pressure_input  ...
            ' -wbb ' pressure_output ...
            ];
        adaptativeStep = [ adaptativeStep ' > ./temps/out_bamg_interp_press.txt'];
        eval(adaptativeStep)
        ux = [];
        uy = [];
    else
        fn_ux_out = './temps/temp_ux_out';
        ux_out = './temps/temp_ux_out.bb';    
        fn_uy_out = './temps/temp_uy_out';
        uy_out = './temps/temp_uy_out.bb';
        adaptativeStep = [ bamg_exe  ... 
            ' -b ' mesh_input        ...
            ' -r ' mesh_output       ...
            ' -rbb ' pressure_input  ...
            ' -wbb ' pressure_output ...
            ];
        adaptativeStep = [ adaptativeStep ' > ./temps/out_bamg_interp_press.txt'];
        eval(adaptativeStep)
        adaptativeStep = [ bamg_exe  ... 
            ' -b ' mesh_input        ...
            ' -r ' mesh_output       ...
            ' -rbb ' ux_input  ...
            ' -wbb ' ux_out ...
            ];
        adaptativeStep = [ adaptativeStep ' > ./temps/out_bamg_interp_ux.txt'];
        eval(adaptativeStep)
        adaptativeStep = [ bamg_exe  ... 
            ' -b ' mesh_input        ...
            ' -r ' mesh_output       ...
            ' -rbb ' uy_input  ...
            ' -wbb ' uy_out ...
            ];
        adaptativeStep = [ adaptativeStep ' > ./temps/out_bamg_interp_uy.txt'];
        eval(adaptativeStep)
        ux = readBAMGbb(fn_ux_out);
        uy = readBAMGbb(fn_uy_out);
        delete(ux_out);
        delete(uy_out);
    end
    
    mesh = readBAMGmsh(fileName_mesh_output);
    
    Xh = mesh.X';
    Th = mesh.T;
    u  = readBAMGbb(fn_sol_output);
    p  = readBAMGbb(fn_press_output);
    
    delete(mesh_input);
    delete(mesh_output);
    delete(solution_output);
    delete(pressure_output);
    
end


%     if adapt_variable.comb == 1
%         f_file = printSOL(fileName,U);
%         msh = './temps/z_temp.msh';
%     else
%         f_file = printSOL(fileName,U);
%         msh = './temps/z_temp.msh';
%     end
%     fileNameP = './temps/z_temp_sol_p.bb';
%     
%     %create the instruction to Bamg
%     if adapt_variable.comb == 1
%         adaptativeStep = instructionBamg3(bamg_exe,msh,data,f_file);
%     else
%         if adapt_variable.d == true && adapt_variable.u == true && adapt_variable.p == true
%             adaptativeStep = instructionBamg4(bamg_exe,msh,data,f_file,fileNamesol,fileNameP);
% 
%         elseif adapt_variable.d == true && adapt_variable.u == true && adapt_variable.p == false
%             adaptativeStep = instructionBamg2(bamg_exe,msh,data,f_file,fileNamesol);
% 
%         elseif adapt_variable.d == true && adapt_variable.u == false && adapt_variable.p == true
%             adaptativeStep = instructionBamg2(bamg_exe,msh,data,f_file,fileNameP);
% 
%         elseif adapt_variable.d == true && adapt_variable.u == false && adapt_variable.p == false
%             adaptativeStep = instructionBamg3(bamg_exe,msh,data,f_file,fileNamesol);
% 
%         elseif adapt_variable.d == false && adapt_variable.u == true && adapt_variable.p == true
%             adaptativeStep = instructionBamg2(bamg_exe,msh,data,fileNamesol,fileNameP);
% 
%         elseif adapt_variable.d == false && adapt_variable.u == true && adapt_variable.p == false
%             adaptativeStep = instructionBamg3(bamg_exe,msh,data,fileNamesol);
% 
%         elseif adapt_variable.d == false && adapt_variable.u == false && adapt_variable.p == true
%             adaptativeStep = instructionBamg3(bamg_exe,msh,data,fileNameP);
% 
%         else
%             error('wrong combination');
%         end
%     end
    



%     disp('aqui hi ha error segons quina combinacio de variables')
%     disp('no se exactament que es, REVISAR comandes/fitxers/etc')
%     fprintf('\n\n\n\n\n\n\n\n')
%     disp('aquest dona error... quan intento p i u nomes, pero llegeix la malla antiga')
%     disp('com que es diuen igual, no et dones compte si escriu molt per pantalla')
%     disp('es necessari fer mes segur aixo, escrivint el resultat amb un nom diferent')
%     disp('el mateix passa si el resultat no es llegeix auqi sino despres duna estona en una altra funcio')
%     error('aaaa')
    
